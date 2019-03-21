from Bio import Seq
import gzip
import sys
from Bio.Phylo.PAML import baseml
import multiprocessing
from multiprocessing import Pool
from ete3 import PhyloTree
import subprocess
import os
from MAFSeq import MAFSeq
from MAFPair import MAFPair
import numpy
import copy
from Bio import SeqIO

fsa_path = "/Genomics/kocherlab/berubin/local/src/fsa_long/bin"
raxml_path = "/Genomics/kocherlab/berubin/local/src/standard-RAxML"
paml_path = "/Genomics/kocherlab/berubin/local/src/paml4.9e/bin"
trimal_path = "/Genomics/kocherlab/berubin/local/src/trimal/source"

HIT_OVERLAP = 0
INTRON_SIZE = 500
FILL_CHAR = "N"

def read_mafs_overlord(base_dir, num_threads, outspecies_list, maf_dir, inspecies, gff_dir, target_scafs):
    scaf_list = []
    reader = open(target_scafs, 'rU')
    for line in reader:
        scaf_list.append(line.strip())
    pool = multiprocessing.Pool(processes = num_threads)
    work_list = []
    for outspecies in outspecies_list:
        maf_file = "%s/%s.%s.sing.maf.gz" % (maf_dir, inspecies, outspecies)
        work_list.append([maf_file, inspecies, outspecies, gff_dir, scaf_list])
    maf_dic_list = pool.map_async(read_mafs_worker, work_list).get(99999999)
    pairs_dic = {}
    for x in range(len(outspecies_list)):
        pairs_dic[outspecies_list[x]] = maf_dic_list[x]
    nearest_dic = maf_dic_list[0]
    return pairs_dic

def read_mafs_worker(param_list):
    maf_file = param_list[0]
    inspecies = param_list[1]
    outspecies = param_list[2]
    gff_dir = param_list[3]
    scaf_list = param_list[4]
    maf_list = read_maf(maf_file, inspecies, outspecies)
    maf_dic = mafs_by_scaf(maf_list, inspecies)
    maf_dic = only_syntenic(maf_dic, outspecies)
    processed_maf_dic = {}
    for scaf, maf_list in maf_dic.items():
        if scaf not in scaf_list:
            continue
        maf_list = sort_hsps(maf_list)
        merged_maf_list = merge_overlord(maf_list, inspecies, outspecies)
        merged_maf_list = remove_contained(maf_list, inspecies)
        ingff = "%s/%s.gff" % (gff_dir, inspecies)
        outgff = "%s/%s.gff" % (gff_dir, outspecies)
        maf_list = mask_proteins(merged_maf_list, inspecies, ingff)
        maf_list = mask_proteins(merged_maf_list, outspecies, outgff)
        processed_maf_dic[scaf] = maf_list
    return processed_maf_dic

def merge_HSPs(my_hsp, temp_hsp_list, target_species, out_species):
    i = 0
    included = False
    border_size = INTRON_SIZE
    query_overlap = HIT_OVERLAP
    add_hsp = False
    while i < len(temp_hsp_list):
        my_hsp_target = my_hsp.get_maf(target_species)
        my_hsp_out = my_hsp.get_maf(out_species)
        cur_hsp = temp_hsp_list[i]
        cur_hsp_target = cur_hsp.get_maf(target_species)
        cur_hsp_out = cur_hsp.get_maf(out_species)
        if my_hsp_target.scaf != cur_hsp_target.scaf:
            i += 1
            continue
        if my_hsp_target.strand != cur_hsp_target.strand:
            i += 1
            continue
        if my_hsp_out.scaf != cur_hsp_out.scaf:
            i += 1
            continue
        if my_hsp_out.strand != cur_hsp_out.strand:
            i += 1
            continue

        if my_hsp_out.strand == "+":
            if my_hsp_target.start <= cur_hsp_target.start and my_hsp_target.end >= cur_hsp_target.start - border_size:
                if my_hsp_target.end >= cur_hsp_target.end:
                    new_hsp = my_hsp
                    my_hsp = new_hsp
                    del temp_hsp_list[i]
                    i = -1
                else:
                    if my_hsp_out.end <= cur_hsp_out.start + query_overlap and my_hsp_out.end + border_size > cur_hsp_out.start:
                        blank_len = cur_hsp_target.start - (my_hsp_target.start + my_hsp_target.orig_len)
                        out_blank_len = cur_hsp_out.start - (my_hsp_out.start + my_hsp_out.orig_len)
                        target_gaps = 0
                        out_gaps = 0
                        if blank_len > out_blank_len:
                            out_gaps = blank_len - out_blank_len
                        elif out_blank_len > blank_len:
                            target_gaps = out_blank_len - blank_len
                        target_seq = my_hsp_target.seq + FILL_CHAR*blank_len + "-"*target_gaps + cur_hsp_target.seq
                        out_seq = my_hsp_out.seq + FILL_CHAR*out_blank_len + "-"*out_gaps + cur_hsp_out.seq 
                        target_maf = MAFSeq(my_hsp_target.species, my_hsp_target.scaf, my_hsp_target.start, cur_hsp_target.end - my_hsp_target.start, target_seq, my_hsp_target.strand, my_hsp_target.scaf_len, cur_hsp_target.end)
                        out_maf = MAFSeq(my_hsp_out.species, my_hsp_out.scaf, my_hsp_out.start, cur_hsp_out.end - my_hsp_out.start, out_seq, my_hsp_out.strand, my_hsp_out.scaf_len, cur_hsp_out.end)
                        new_hsp = MAFPair(target_maf.species, target_maf, out_maf.species, out_maf)
                        
                        my_hsp = new_hsp
                        add_hsp = True
                        insert_pos = i
                        del temp_hsp_list[i]
                        i = -1
            elif my_hsp_target.end >= cur_hsp_target.end and my_hsp_target.start <= cur_hsp_target.end + border_size:

                if cur_hsp_out.end <= my_hsp_out.start + query_overlap and cur_hsp_out.end + border_size > my_hsp_out.start:
                    blank_len = my_hsp_target.start - (cur_hsp_target.start + cur_hsp_target.orig_len)
                    out_blank_len = my_hsp_out.start - (cur_hsp_out.start + cur_hsp_out.orig_len)
                    target_gaps = 0
                    out_gaps = 0
                    if blank_len > out_blank_len:
                        out_gaps = blank_len - out_blank_len
                    elif out_blank_len > blank_len:
                        target_gaps = out_blank_len - blank_len
                    target_seq = cur_hsp_target.seq + FILL_CHAR*blank_len + "-"*target_gaps + my_hsp_target.seq
                    out_seq = cur_hsp_out.seq + FILL_CHAR*out_blank_len + "-"*out_gaps + my_hsp_out.seq 
                    target_maf = MAFSeq(my_hsp_target.species, my_hsp_target.scaf, cur_hsp_target.start, my_hsp_target.end - cur_hsp_target.start, target_seq, my_hsp_target.strand, my_hsp_target.scaf_len, my_hsp_target.end)
                    out_maf = MAFSeq(my_hsp_out.species, my_hsp_out.scaf, cur_hsp_out.start, my_hsp_out.end - cur_hsp_out.start, out_seq, my_hsp_out.strand, my_hsp_out.scaf_len, my_hsp_out.end)
                    new_hsp = MAFPair(target_maf.species, target_maf, out_maf.species, out_maf)

                    my_hsp = new_hsp
                    add_hsp = True
                    insert_pos = i
                    del temp_hsp_list[i]
                    i = -1
        elif my_hsp_out.strand == "-":

            if my_hsp_target.start <= cur_hsp_target.start and my_hsp_target.end >= cur_hsp_target.start - border_size:

                if my_hsp_target.end >= cur_hsp_target.end:

                    new_hsp = my_hsp
                    my_hsp = new_hsp
                    del temp_hsp_list[i]
                    i = -1
                else:
                    if my_hsp_out.start >= cur_hsp_out.end - query_overlap and cur_hsp_out.end + border_size > my_hsp_out.start:

                        blank_len = cur_hsp_target.start - (my_hsp_target.start + my_hsp_target.orig_len)
                        out_blank_len = my_hsp_out.start - (cur_hsp_out.start + cur_hsp_out.orig_len)
                        target_gaps = 0
                        out_gaps = 0
                        if blank_len > out_blank_len:
                            out_gaps = blank_len - out_blank_len
                        elif out_blank_len > blank_len:
                            target_gaps = out_blank_len - blank_len
                        target_seq = my_hsp_target.seq + FILL_CHAR*blank_len + "-"*target_gaps + cur_hsp_target.seq
                        out_seq = my_hsp_out.seq + FILL_CHAR*out_blank_len + "-"*out_gaps + cur_hsp_out.seq 
                        target_maf = MAFSeq(my_hsp_target.species, my_hsp_target.scaf, my_hsp_target.start, cur_hsp_target.end - my_hsp_target.start, target_seq, my_hsp_target.strand, my_hsp_target.scaf_len, cur_hsp_target.end)

                        out_maf = MAFSeq(my_hsp_out.species, my_hsp_out.scaf, cur_hsp_out.start, my_hsp_out.end - cur_hsp_out.start, out_seq, cur_hsp_out.strand, my_hsp_out.scaf_len, my_hsp_out.end)
                        new_hsp = MAFPair(target_maf.species, target_maf, out_maf.species, out_maf)
                        my_hsp = new_hsp
                        add_hsp = True
                        insert_post = i
                        del temp_hsp_list[i]
                        i = -1
            elif my_hsp_target.end >= cur_hsp_target.end and my_hsp_target.start <= cur_hsp_target.end + border_size:
                if cur_hsp_out.start >= my_hsp_out.end - query_overlap and my_hsp_out.end + border_size > cur_hsp_out.end:
                    blank_len = my_hsp_target.start - (cur_hsp_target.start + cur_hsp_target.orig_len)
                    out_blank_len = cur_hsp_out.start - (my_hsp_out.start + my_hsp_out.orig_len)
                    target_gaps = 0
                    out_gaps = 0
                    if blank_len > out_blank_len:
                        out_gaps = blank_len - out_blank_len
                    elif out_blank_len > blank_len:
                        target_gaps = out_blank_len - blank_len
                    target_seq = cur_hsp_target.seq + FILL_CHAR*blank_len + "-"*target_gaps + my_hsp_target.seq
                    out_seq = cur_hsp_out.seq + FILL_CHAR*out_blank_len + "-"*out_gaps + my_hsp_out.seq 
                    target_maf = MAFSeq(my_hsp_target.species, my_hsp_target.scaf, cur_hsp_target.start, my_hsp_target.end - cur_hsp_target.start, target_seq, my_hsp_target.strand, my_hsp_target.scaf_len, my_hsp_target.end)
                    out_maf = MAFSeq(my_hsp_out.species, my_hsp_out.scaf, my_hsp_out.start, cur_hsp_out.end - my_hsp_out.start, out_seq, cur_hsp_out.strand, my_hsp_out.scaf_len, cur_hsp_out.end)
                    new_hsp = MAFPair(target_maf.species, target_maf, out_maf.species, out_maf)
                    my_hsp = new_hsp
                    add_hsp = True
                    insert_post = i
                    del temp_hsp_list[i]
                    i = -1                
        i += 1
    temp_hsp_list.append(my_hsp)
    return temp_hsp_list


def overlap(mytuple, tuplelist):
    included = False
    i = 0
    merged = []
    while i < len(tuplelist):
        cur_tuple = tuplelist[i]
        if mytuple[0] <= cur_tuple[0] and mytuple[1] >= cur_tuple[0]-1000:
            if mytuple[1] >= cur_tuple[1]:
                
                new_tuple = mytuple
                mytuple = new_tuple
                del tuplelist[i]
                i = -1
            else:
                new_tuple = (mytuple[0], cur_tuple[1])
                mytuple = new_tuple
                del tuplelist[i]
                i = -1
        elif mytuple[1] >= cur_tuple[1] and mytuple[0] <= cur_tuple[1] + 1000:
            if mytuple[0] >= cur_tuple[0]:
                new_tuple = (cur_tuple[0], mytuple[1])
                mytuple = new_tuple
                del tuplelist[i]
                i = -1   
        elif mytuple[0] >= cur_tuple[0] and mytuple[1] <= cur_tuple[1]:
            included = True
            break
        i += 1
    if not included:
        tuplelist.append(mytuple)
    return tuplelist

def distance(seq1, seq2):
    diff_count = 0.0
    base_count = 0.0
    for x in range(len(seq1) - 1):
        if seq1[x] == "-" or seq2[x] == "-":
            continue
        if seq1[x] == "N" or seq2[x] == "N":
            continue
        base_count += 1
        if seq1[x].upper() != seq2[x].upper():
            diff_count += 1
    return diff_count, base_count


def sliding_window(maf_list, window_size, window_step, inspecies, outspecies, output_dir, cur_scaf):
    target_scaf_dists = {}
    out_scaf_dists = {}
    outfile = open("%s/%s_vs_%s_%s_%s_%s_distances.txt" % (output_dir, inspecies, outspecies, window_size, window_step, cur_scaf), 'w')    
    outfile.write("inspecies\tinscaf\tincoord\toutspecies\toutscaf\toutcoord\tdxy\n")
    for mafpair in maf_list:
        print mafpair
        distance_tuples = mafpair.sliding_distance(window_size, window_step, inspecies, outspecies)
    for mafpair in maf_list:
        print "printspot"
        print mafpair
        inmaf = mafpair.get_maf(inspecies)
        outmaf = mafpair.get_maf(outspecies)
        for dist in mafpair.distances:
            print dist
            outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (inspecies, inmaf.scaf, dist[0], outspecies, outmaf.scaf, dist[1], dist[2]))
        outfile.flush()
    outfile.close()

def sliding_window_transfer_frame(maf_dic, frame_maf_dic, window_size, window_step, inspecies, outspecies, framespecies, output_dir):
    target_scaf_dists = {}
    out_scaf_dists = {}
    outfile = open("%s/%s_vs_%s_%s_%s_distances_frame_%s.txt" % (output_dir, inspecies, outspecies, window_size, window_step, framespecies), 'w')    
    outfile.write("inspecies\tinscaf\tincoord\toutspecies\toutscaf\toutcoord\tdxy\tframescaf\tframecoord\n")
    for scaf, maf_list in maf_dic.items():
        for mafpair in maf_list:
            distance_tuples = mafpair.sliding_distance(window_size, window_step, inspecies, outspecies)
        for mafpair in maf_list:
            mafpair.transfer_coords(frame_maf_dic[scaf], inspecies, framespecies)
        for mafpair in maf_list:
            inmaf = mafpair.get_maf(inspecies)
            outmaf = mafpair.get_maf(outspecies)
            for dist in mafpair.framed_distances:
                
                outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (inspecies, inmaf.scaf, dist[0], outspecies, outmaf.scaf, dist[1], dist[2], dist[3], dist[4]))
    outfile.close()


def write_distances(maf_list, window_size, window_step, inspecies, outspecies, output_dir):
    target_scaf_dists, out_scaf_dists = sliding_window(maf_list, window_size, window_step, inspecies, outspecies)
    outfile = open("%s/%s_vs_%s_%s_%s_distances.txt" (output_dir, inspecies, outspecies, window_size, window_step), 'w')
    

def read_genome(genome_file):
    reader = SeqIO.parse(genome_file, format = 'fasta')
    seq_dic = {}
    for rec in reader:
        seq_dic[rec.id.split(":")[1]] = str(rec.seq)
    return seq_dic

def read_maf(inmaf, inspecies, outspecies):
    reader = gzip.open(inmaf, 'rb')
    maf_blocks = []
    seq_count = 0
    for line in reader:
        if line.startswith("#"):
            continue
        if line.startswith("a"):
            if seq_count > 0:
                if cur_score > 1000:
                    if seq1.species == inspecies and seq2.species == outspecies:
                        maf_blocks.append(MAFPair(inspecies, seq1, outspecies, seq2))
                    elif seq1.species == outspecies and seq2.species == inspecies:
                        maf_blocks.append(MAFPair(inspecies, seq2, outspecies, seq1))
            cur_line = line.split()
            cur_score = int(cur_line[1].split("score=")[1])
            cur_mismap = cur_line[2]
            seq_count = 0
        if line.startswith("s"):
            cur_line = line.split()
            cur_name = cur_line[1].split(":")
            if seq_count == 0:
                cur_species = inspecies
            else:
                cur_species = outspecies
            cur_scaf = cur_line[1]
            cur_start = int(cur_line[2])
            cur_len = int(cur_line[3])
            cur_strand = cur_line[4]
            cur_scaf_len = int(cur_line[5])
            cur_seq = cur_line[6].upper()

            if cur_strand == "-":
                cur_end = cur_scaf_len - cur_start
                cur_start = cur_end - cur_len
            else:
                cur_end = cur_start + cur_len

            if seq_count == 0:
                seq1 = MAFSeq(cur_species, cur_scaf, cur_start, cur_len, cur_seq, cur_strand, cur_scaf_len, cur_end)
            else:
                seq2 = MAFSeq(cur_species, cur_scaf, cur_start, cur_len, cur_seq, cur_strand, cur_scaf_len, cur_end)
            seq_count += 1
    return maf_blocks

def restrict_region(in_maf_list, inspecies):
    #hardcoded restriction of region to include in the pipeline. For testing.
    new_maf_list = []
    for in_maf in in_maf_list:
        target = in_maf.get_maf(inspecies)
        if target.start > 998000 and target.start < 1000000:
            new_maf_list.append(in_maf)
    return new_maf_list

def mafs_by_scaf(maf_blocks, inspecies):
    scaf_dic = {}
    for maf_pair in maf_blocks:
        cur_scaf = maf_pair.get_maf(inspecies).scaf 
        if cur_scaf not in scaf_dic.keys():
            scaf_dic[cur_scaf] = []
        scaf_dic[cur_scaf].append(maf_pair)
    return scaf_dic

def only_syntenic(maf_dic, outspecies):
    new_maf_dic = {}
    for scaf, maf_list in maf_dic.items():
        syntenics = syntenic_scaf(maf_list, outspecies)
        new_maf_dic[scaf] = []
        for mafpair in maf_list:
            if mafpair.get_maf(outspecies).scaf in syntenics:
                new_maf_dic[scaf].append(mafpair)
    
    return new_maf_dic

def syntenic_scaf(maf_list, outspecies):
    other_scaf_dic = {}
    for maf_pair in maf_list:
        other_scaf = maf_pair.get_maf(outspecies).scaf
        if other_scaf not in other_scaf_dic.keys():
            other_scaf_dic[other_scaf] = 0
        other_scaf_dic[other_scaf] += maf_pair.get_maf(outspecies).orig_len
    most_scaf = ""
    most_scaf_bases = 0
    for scaf, count in other_scaf_dic.items():
        if count > most_scaf_bases:
            most_scaf_bases = count
            most_scaf = scaf
    other_scaf_props = {}
    for scaf, count in other_scaf_dic.items():
        other_scaf_props[scaf] = 1.0 * count / most_scaf_bases
    syntenics = []
    for scaf, prop in other_scaf_props.items():
        if prop > 0.05:
            syntenics.append(scaf)
    return syntenics

def remove_shorts(hsp_list, target_species):
    new_hsp_list = []
    for hsp in hsp_list:
        mafseq = hsp.get_maf(target_species)
        if mafseq.orig_len > 1000:
            new_hsp_list.append(hsp)
    return new_hsp_list

def sort_hsps(hsp_list):
    new_hsp_list = []
    hsp_count = len(hsp_list)
    while True:
        if len(hsp_list) == 0:
            break
        smallest_hsp = hsp_list[0]
        for hsp in hsp_list:
            if hsp.less_than(smallest_hsp):
                smallest_hsp = hsp
        new_hsp_list.append(smallest_hsp)
        hsp_list.remove(smallest_hsp)
        if len(new_hsp_list) == hsp_count:
            break
    return new_hsp_list

def remove_contained(hsp_list, target_species):
    new_hsp_list = []
    for hsp in hsp_list:
        cur_hsp = hsp.get_maf(target_species)
        keep_hsp = True
        for contained_hsp in hsp_list:
            cur_contained_hsp = contained_hsp.get_maf(target_species)
            if cur_contained_hsp.start < cur_hsp.start and cur_contained_hsp.end > cur_hsp.end:
                keep_hsp = False
                break
        if keep_hsp:
            new_hsp_list.append(hsp)
    return new_hsp_list

def merge_overlord(maf_list, inspecies, outspecies):
    
    x = 0
    temp_maf_list = copy.deepcopy(maf_list)
    cur_maf_len = len(maf_list)
    while x < cur_maf_len:
        temp_maf_list = merge_HSPs(temp_maf_list[x], temp_maf_list, inspecies, outspecies)
        if len(temp_maf_list) < cur_maf_len:
            x = 0
            cur_maf_len = len(temp_maf_list)
        else:
            x += 1

    temp_maf_list = remove_shorts(temp_maf_list, inspecies)
    temp_maf_list = sort_hsps(temp_maf_list)
    return temp_maf_list

def mask_proteins(maf_list, species, gff_file):
    masked_maf_list = []
    reader = open(gff_file, 'rU')
    cds_dic = {}
    for line in reader:
        if line.startswith("#"):
            continue
        if line.startswith(">"):
            break
        cur_line = line.split()
        cur_type = cur_line[2]
        if cur_type == "CDS":
            cur_scaf = cur_line[0].replace("_scaf_", "_scaff_")
            if cur_scaf not in cds_dic.keys():                
                cds_dic[cur_scaf] = []
            cds_dic[cur_scaf].append((int(cur_line[3]), int(cur_line[4])))
    for mafpair in maf_list:
        target_maf = mafpair.get_maf(species)
        seq_len = len(target_maf.seq)
        if target_maf.scaf in cds_dic.keys():
            for cds in cds_dic[target_maf.scaf]:
                if target_maf.strand == "+" or target_maf.strand == "-":
                    if cds[1] > target_maf.start and cds[0] < target_maf.end:
                        if cds[0] > target_maf.start:
                            start_coord = liftover(cds[0] - target_maf.start, target_maf.seq, cds[0] - target_maf.start)
                        elif cds[0] <= target_maf.start:
                            start_coord = 1
                        if cds[1] < target_maf.end:
                            end_coord = liftover(cds[1] - target_maf.start, target_maf.seq, start_coord)
                        elif cds[1] >= target_maf.end:
                            end_coord = seq_len - 1
                        insert_seq = target_maf.seq[start_coord-1:end_coord].upper()
                        insert_seq = insert_seq.replace("A", "a").replace("C", "c").replace("T", "t").replace("G","g")
                        target_maf.seq = target_maf.seq[:start_coord-1] + insert_seq + target_maf.seq[end_coord:]
                        if len(target_maf.seq) != seq_len:
                            print "DISASTER"
                            print target_maf
                            print mafpair
                            print target_maf.seq
                            print mafpair.maf1.seq
                            sys.exit()
                elif target_maf.strand == "0":
                    start_coord = (target_maf.scaf_len - cds[1]) - (target_maf.scaf_len - target_maf.end)
                    end_coord = (target_maf.scaf_len - cds[0]) - (target_maf.scaf_len - target_maf.end)
                    if cds[1] > target_maf.start and cds[0] < target_maf.end:
                        start_coord = liftover(start_coord, target_maf.seq, start_coord)
                        end_coord = liftover(end_coord, target_maf.seq, start_coord)
                        if end_coord > seq_len:
                            end_coord = seq_len

                        insert_seq = target_maf.seq[start_coord:end_coord+1].upper()
                        insert_seq = insert_seq.replace("A", "a").replace("C", "c").replace("T", "t").replace("G","g")
                        target_maf.seq = target_maf.seq[:start_coord] + insert_seq + target_maf.seq[end_coord+1:]
                        if len(target_maf.seq) != seq_len:
                            print "DISASTER"
                            print target_maf
                            print mafpair
                            print target_maf.seq
                            print mafpair.maf1.seq
                            print start_coord
                            print end_coord
                            print insert_seq
                            sys.exit()

    return maf_list
            

def liftover(target_index, gapped_seq, i_start):
    i = i_start
    prev_i = i
    prev_gap = gapped_seq[:i].count("-")
    while i < len(gapped_seq):
        cur_gap_count = prev_gap + gapped_seq[prev_i:i].count("-")
        prev_gap = cur_gap_count
        prev_i = i
        if i - cur_gap_count == target_index:
            return i
        if i - cur_gap_count + 10000 < target_index:
            i += 10000
        elif i - cur_gap_count + 1000 < target_index:
            i += 1000
        elif i - cur_gap_count + 100 < target_index:
            i += 100
        elif i - cur_gap_count + 10 < target_index:
            i += 10
        else:
            i += 1
    return i

def write_mafs_to_file(maf_list, outdir, inspecies, outspecies):
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    if not os.path.exists("%s/maf_files" % outdir):
        os.mkdir("%s/maf_files" % outdir)
    for maf in maf_list:
        target_maf = maf.get_maf(inspecies)
        out_maf = maf.get_maf(outspecies)
        outfile = open("%s/maf_files/%s_%s_%s_%s_%s_%s.fa" % (outdir, target_maf.scaf, target_maf.start, target_maf.end, out_maf.scaf, out_maf.start, out_maf.end), 'w')
        outfile.write(">%s:%s:%s:%s:%s\n%s\n" % (inspecies, target_maf.scaf, target_maf.start, target_maf.end, target_maf.strand, target_maf.seq.replace("-", "-").swapcase()))
        outfile.write(">%s:%s:%s:%s:%s\n%s\n" % (outspecies, out_maf.scaf, out_maf.start, out_maf.end, out_maf.strand, out_maf.seq.replace("-", "").swapcase()))
        for other_maf in maf.mafothers:
            outfile.write(">%s:%s:%s:%s:%s\n%s\n" % (other_maf.species, other_maf.scaf, other_maf.start, other_maf.end, other_maf.strand, other_maf.seq.replace("N", "").replace("-", "").swapcase()))
        outfile.close()


def filter_small_mafs(maf_list, min_taxa):
    new_maf_list = []
    for maf in maf_list:
        species_list = []
        tuple_dic = {}
        if len(maf.mafothers) < min_taxa - 2:
            continue
        for other_maf in maf.mafothers:
            if other_maf.species not in species_list:
                tuple_dic[other_maf.species] = (other_maf.other_start, other_maf.other_end)
                species_list.append(other_maf.species)
            else:
                if other_maf.other_start < tuple_dic[other_maf.species][0]:
                    tuple_dic[other_maf.species] = (other_maf.other_start, tuple_dic[other_maf.species][1])
                else:
                    tuple_dic[other_maf.species] = (tuple_dic[other_maf.species][0], other_maf.other_end)
        merged_tuples, merged_count = merge_tuples(tuple_dic.values())

            
        if len(species_list) >= min_taxa-2 and merged_count >= min_taxa - 2:
            new_maf_list.append(maf)
    return new_maf_list

def merge_tuples(tuple_list):
    temp_list = copy.deepcopy(tuple_list)
    x = 0
    first_tuple = temp_list[0]
    first_index = 0
    merge_count = 0
    while first_index < len(temp_list):
        first_tuple = temp_list[first_index]
        while x < len(temp_list):
            if len(temp_list) == 1:
                break
            new_tuple = False
            if x == first_index:
                x += 1
                if x == len(temp_list):
                    break
                
            cur_tuple = temp_list[x]
            if first_tuple[0] >= cur_tuple[0] and first_tuple[1] <= cur_tuple[1]:
                new_tuple = cur_tuple
            elif first_tuple[0] <= cur_tuple[0] and first_tuple[1] >= cur_tuple[1]:
                new_tuple = first_tuple

            elif first_tuple[0] < cur_tuple[0] and first_tuple[1] >= cur_tuple[0]:
                new_tuple = (first_tuple[0], cur_tuple[1])
            elif first_tuple[0] <= cur_tuple[1] and first_tuple[1] > cur_tuple[1]:
                new_tuple = (cur_tuple[0], first_tuple[1])
            if new_tuple:

                temp_list.remove(cur_tuple)
                temp_list.remove(first_tuple)
                temp_list.append(new_tuple)
                x = 0
                first_index = 0
                first_tuple = temp_list[first_index]
                merge_count += 1
            else:
                x += 1
        if len(temp_list) == 1:
            break
        first_index += 1
        x = 0
    return temp_list, merge_count

                        
def realign_fsa(maf_scaf_dic, outdir, inspecies, outspecies, num_threads):
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    pool = multiprocessing.Pool(processes = num_threads)
    work_list = []
    for scaf, maf_list in maf_scaf_dic.items():
        for maf in maf_list:
            target_maf = maf.get_maf(inspecies)
            out_maf = maf.get_maf(outspecies)
            infile = "%s/maf_files/%s_%s_%s_%s_%s_%s.fa" % (outdir, target_maf.scaf, target_maf.start, target_maf.end, out_maf.scaf, out_maf.start, out_maf.end)
            outfile = "%s/maf_files/%s_%s_%s_%s_%s_%s.afa" % (outdir, target_maf.scaf, target_maf.start, target_maf.end, out_maf.scaf, out_maf.start, out_maf.end)
            if os.path.exists(outfile):
                num_needed = 0
                reader = SeqIO.parse(infile, format = 'fasta')
                for rec in reader:
                    num_needed += 1
                counter = 0
                reader = SeqIO.parse(outfile, format = 'fasta')
                for rec in reader:
                    counter += 1
                if counter == num_needed:
                    continue
            work_list.append([infile, outfile])
    pool.map_async(realign_fsa_worker, work_list).get(99999999)


def realign_fsa_worker(param_list):
    infile = param_list[0] 
    outfile = param_list[1]
    cmd = ["%s/fsa" % fsa_path, "--anchored", "--exonerate", "--softmasked", "--maxram", "18000", infile]
    FNULL = open(os.devnull, 'w')
    with open(outfile, 'w') as outwriter:
        subprocess.call(cmd, stdout = outwriter, stderr = FNULL)
    outwriter.close()


def make_taxa_dic(infile):
    #This reads in the minimum required taxonomy file. It returns two
    #dictionaries. manda_dic has tuples of taxon names as keys and 
    #the number of those taxa required as values. multi_dic has a tuple
    #of tuples as keys. The idea here is that the nested tuples represent
    #sets of species that must appear together. The values are the number
    #of nested tuples that must appear.
    reader = open(infile, 'rU')
    manda_dic = {}
    multi_dic = {}
    remove_list = []
    for line in reader:
        cur_line = line.split()
        if int(cur_line[0]) < 0:
            remove_list.append(cur_line[1])
        if "(" in cur_line[1]:
            taxa_list = cur_line[1].split("),(")
            first_pair = taxa_list[0].replace(")","").replace("(", "").split(",")
            first_tuple = (first_pair[0],)
            for taxon in first_pair[1:]:
                first_tuple = first_tuple + (taxon,)
            taxa_tuple = (first_tuple,)

            for taxon in taxa_list[1:]:
                next_pair = taxon.replace(")","").replace("(", "").split(",")
                next_tuple = (next_pair[0],)
                for next_taxon in next_pair[1:]:
                    next_tuple = next_tuple + (next_taxon,)
                taxa_tuple = taxa_tuple + (next_tuple,)
            multi_dic[taxa_tuple] = int(cur_line[0])
                
        else:
            taxa_list = cur_line[1].split(",")
            taxa_tuple = (taxa_list[0],)
            for taxon in taxa_list[1:]:
                taxa_tuple = taxa_tuple + (taxon,)
            manda_dic[taxa_tuple] = int(cur_line[0])
    return manda_dic, multi_dic, remove_list

def min_taxa_membership(manda_dic, multi_dic, remove_list, index_file, min_taxa):
    reader = open(index_file, 'rU')
    og_list = []
    for line in reader:
        if line.startswith("#"):
            continue
        include = True
        cur_line = line.split()
        taxa_to_remove = []
        taxa_list = [seq_name[0:4] for seq_name in cur_line[2].split(",")]
        for taxa in taxa_list:
            if taxa in remove_list:
                taxa_to_remove.append(taxa)
        if int(cur_line[1]) - len(remove_list) < min_taxa:
            continue
        for these_taxa, min_count in manda_dic.items():
            mycount = 0
            for this_taxon in these_taxa:
                if this_taxon in taxa_list:
                    mycount += 1
            if mycount < min_count:
                include = False
                break
        if not include:
            continue
        for these_sets, min_count in multi_dic.items():
            mycount = 0
            for this_set in these_sets:
                set_total = 0
                for this_taxon in this_set:
                    if this_taxon in taxa_list:
                        set_total += 1
                if set_total == len(this_set):
                    mycount += 1
            if mycount < min_count:
                include = False
        if include:
            og_list.append(cur_line[0])
    return og_list


def countseqs(fasta_file):
    counter = 0
    reader = SeqIO.parse(fasta_file, format = 'fasta')
    for rec in reader:
        counter += 1
    return counter

def genome_blast_dbs(genome_paths, outdir):
    if not os.path.exists("%s/blastdbs" % outdir):
        os.mkdir("%s/blastdbs" % outdir)
    blast_working = "%s/blastdbs" % outdir
    reader = open(genome_paths, 'rU')
    seq_dic = {}
    for line in reader:
        cur_line = line.split()
        cmd = ["makeblastdb", "-in", cur_line[1], "-out", "%s/%s_db" % (blast_working, cur_line[0]), "-dbtype", "nucl"]
        subprocess.call(cmd)

def genome_seqs(genome_paths):
    reader = open(genome_paths, 'rU')
    seq_dic = {}
    for line in reader:
        cur_line = line.split()
        seq_reader = SeqIO.parse(cur_line[1], format = 'fasta')
        seq_dic[cur_line[0]] = {}
        for rec in seq_reader:
            seq_dic[cur_line[0]][rec.id] = str(rec.seq)
    return seq_dic

def slide_baby_slide(maf_scaf_dic, outdir, inspecies, outspecies, window_size, window_step, min_taxa, constraint_tree, num_threads, genome_paths, gff_dir):
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
#    if not os.path.isdir(outdir + "/raxml_phylos"):
#        os.mkdir(outdir + "/raxml_phylos")
    if not os.path.isdir(outdir + "/trimal_windows"):
        os.mkdir(outdir + "/trimal_windows")
    if not os.path.isdir(outdir + "/filtered_loci"):
        os.mkdir(outdir + "/filtered_loci")
    if not os.path.isdir(outdir + "/blast_working"):
        os.mkdir(outdir + "/blast_working")
    if not os.path.isdir(outdir + "/ncar_annotations"):
        os.mkdir(outdir + "/ncar_annotations")

    pool = multiprocessing.Pool(processes = num_threads)
    work_list = []
    for scaf, maf_list in maf_scaf_dic.items():
#        if scaf != "NMEL_chr_6":
#            continue
        for maf in maf_list:
            target_maf = maf.get_maf(inspecies)
            out_maf = maf.get_maf(outspecies)
            infile = "%s/maf_files/%s_%s_%s_%s_%s_%s.afa" % (outdir, target_maf.scaf, target_maf.start, target_maf.end, out_maf.scaf, out_maf.start, out_maf.end)
            if os.path.exists(infile):
                if countseqs(infile) > 0:
                    work_list.append([maf, outdir, inspecies, outspecies, window_size, window_step, min_taxa, constraint_tree])
#                    slide_baby_slide_worker([maf, outdir, inspecies, outspecies, window_size, window_step, min_taxa, constraint_tree])
    print len(work_list)
    blens = pool.map_async(slide_baby_slide_worker, work_list).get(99999999)

    outfile = open("%s/filtered_loci.index" % outdir, 'w')
    gff_file_dic = {}
    ncar_file_dic = {}
    ncar_index = 10000
    genome_seq_dic = genome_seqs(genome_paths)
    genome_seq_dic = mask_unaligned_genomes(genome_seq_dic, gff_dir)
    for locus_dic in blens:
        for maf, blens_dic in locus_dic.items():
            for windex, trimal_dic in blens_dic.items():
                if inspecies not in [taxon[0:4] for taxon in trimal_dic.keys()]:
                    continue
                real_space_coord = windex - maf.realigned_seq[0:windex].count("-")
                windex_seq = maf.realigned_seq[windex:windex+500].replace("N","n")
                if windex_seq.startswith("n") or windex_seq.lstrip("-").startswith("n"):
                    n_count = len(windex_seq) - len(windex_seq.replace("-", "n").lstrip("n"))
                    real_space_coord = real_space_coord + n_count
                cur_locus = "%s_%s_%s_%s_%s_%s" % (maf.scaf, maf.start, maf.end, windex, real_space_coord, maf.strand)
                seqfile = open("%s/filtered_loci/NCAR_%s.afa" % (outdir, ncar_index), 'w')
                for rec_id, seq in trimal_dic.items():
                    seqfile.write(">%s\n%s\n" % (rec_id, str(seq)))
                    if rec_id[0:4] not in gff_file_dic.keys():
                        gff_file_dic[rec_id[0:4]] = open("%s/ncar_annotations/%s_ncars.gff" % (outdir, rec_id[0:4]), 'w')
                        ncar_file_dic[rec_id[0:4]] = open("%s/ncar_annotations/%s_ncars.fasta" % (outdir, rec_id[0:4]), 'w')
                    seq_features = rec_id.split(":")
                    cur_scaf = seq_features[1]
                    cur_start = int(seq_features[2])
                    cur_end = int(seq_features[3])
                    cur_strand = seq_features[4]
                    gff_file_dic[rec_id[0:4]].write("%s\tNCARpipe\tncar\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s;\n" % (cur_scaf, cur_start, cur_end, cur_strand, rec_id, rec_id))
                    cur_seq = genome_seq_dic[rec_id[0:4]][cur_scaf][cur_start:cur_end]
                    if cur_strand == "-":
                        cur_seq = str(Seq.Seq(cur_seq).reverse_complement())
                    ncar_file_dic[rec_id[0:4]].write(">%s\n%s\n" % (rec_id, cur_seq))
                seqfile.close()
                outfile.write("%s\t%s\t%s\t%s\n" % (cur_locus, ncar_index, len(trimal_dic.keys()), ",".join(trimal_dic.keys())))
                ncar_index += 1

    outfile.close()
    for gff_file in gff_file_dic.values():
        gff_file.close()
    for ncar_file in ncar_file_dic.values():
        ncar_file.close()


def slide_baby_slide_worker(param_list):
    mymaf = param_list[0]
    outdir = param_list[1]
    inspecies = param_list[2]
    outspecies = param_list[3]
    window_size = param_list[4]
    window_step = param_list[5]
    min_taxa = param_list[6]
    constraint_tree = param_list[7]
    blast_working = "%s/blast_working/" % outdir
    target_maf = mymaf.get_maf(inspecies)
    out_maf = mymaf.get_maf(outspecies)
    infile = "%s/maf_files/%s_%s_%s_%s_%s_%s.afa" % (outdir, target_maf.scaf, target_maf.start, target_maf.end, out_maf.scaf, out_maf.start, out_maf.end)

    reader = SeqIO.parse(infile, format = 'fasta')
    seq_dic ={}
    blens_dic = {}
    trimal_dic = {}
    seq_len = 0
    for rec in reader:
        seq_dic[rec.id] = str(rec.seq).swapcase()
        seq_len = len(rec.seq)
        if rec.id[0:4] == inspecies:
            realigned_ref = seq_dic[rec.id]
    x = 0
    while x < seq_len:
        cur_dic = {}
        for seq_id, seq in seq_dic.items():
            cur_seq = seq[x:x+window_size]
#            print seq_id
            cur_feature = seq_id.split(":")
            cur_species = cur_feature[0]
            cur_scaf = cur_feature[1]
#            cur_strand = cur_feature[4]
            cur_start, cur_end, cur_strand = get_ncar_coords(cur_seq, seq_id + ":" + str(x), outdir)
            cur_seq = cur_seq.replace("g", "n").replace("c", "n").replace("t", "n").replace("a", "n")
            new_seq_id = "%s:%s:%s:%s:%s" % (cur_species, cur_scaf, cur_start, cur_end, cur_strand)

            if cur_seq.count("G") + cur_seq.count("C") + cur_seq.count("A") + cur_seq.count("T") > window_size * 0.5:
                cur_dic[new_seq_id] = cur_seq
        if len(cur_dic) > min_taxa:
            trimal_dic[x] = trimal(cur_dic, outdir, infile.split(".afa")[0].split("/")[-1], x, window_size)
        x = x + window_step
    target_maf.realigned_seq = realigned_ref
    return {target_maf : trimal_dic}

def get_ncar_coords(ncar_seq, seq_id, blast_working):
    cur_seq = ncar_seq.replace("-","").upper()
    if len(cur_seq) < 20:
        return -1, -1, "+"
    seq_list = cur_seq.strip("N").split("N")
    blast_file = open("%s/blast_working/%s.query" % (blast_working, seq_id), 'w')
    seq_num = 0
    for seq in seq_list:
        if len(seq) == 0:
            continue
        blast_file.write(">%s_%s\n%s\n" % (seq_id, seq_num, seq))
        seq_num += 1
    blast_file.close()
    cur_species = seq_id.split(":")[0]
    cmd = ["blastn", "-query", "%s/blast_working/%s.query" % (blast_working, seq_id), "-db", "%s/blastdbs/%s_db" % (blast_working, cur_species), "-outfmt", "6", "-out", "%s/blast_working/%s.txt" % (blast_working, seq_id), "-max_target_seqs", "1", "-perc_identity", "95"]
    subprocess.call(cmd)
    reader = open("%s/blast_working/%s.txt" % (blast_working, seq_id), 'rU')
    smallest_start = -1
    biggest_end = -1
    for line in reader:
        cur_line = line.split()
        if float(cur_line[2]) < 90:
            continue
        cur_start = int(cur_line[8])
        cur_end = int(cur_line[9])
        if smallest_start < 0:
            smallest_start = cur_start
        elif cur_start < smallest_start:
            smallest_start = cur_start
        if cur_end > biggest_end:
            biggest_end = cur_end
        if abs(biggest_end - smallest_start) >= len(ncar_seq.strip("N"))-10:
            break
    if biggest_end - smallest_start > 10000:
        print "Huge deletion"
        print seq_id
        print ncar_seq
    cur_strand = "+"
    if biggest_end < smallest_start:
        temp_start = smallest_start
        smallest_start = biggest_end
        biggest_end = temp_start
        cur_strand = "-"
    os.remove("%s/blast_working/%s.query" % (blast_working, seq_id))
    os.remove("%s/blast_working/%s.txt" % (blast_working, seq_id))
#    print "%s\t%s\t%s" % (smallest_start, biggest_end, cur_strand)
    return smallest_start - 1, biggest_end, cur_strand

        
def blengths(seq_dic, outdir, rootname, windex, window_size, min_taxa, constraint_tree):
    seq_dic = trimal(seq_dic, outdir, rootname, windex, window_size)
    if len(seq_dic) < min_taxa:
        return
    seq_file = open("%s/raxml_phylos/%s_%s.afa" % (outdir, rootname, windex), 'w')
    for k, v in seq_dic.items():
        cur_length = len(v)
        break
    seq_file.write("%s %s\n" % (len(seq_dic), cur_length))
    species_list = []
    for k, v in seq_dic.items():
        seq_file.write("%s\n%s\n" % (k, v))
        species_list.append(k)
    seq_file.close()
    tree = PhyloTree(constraint_tree)
    
    tree.prune(species_list)
    tree_str = tree.write(format = 5)
    cons_file = open("%s/raxml_phylos/%s_%s.constraint" % (outdir, rootname, windex), 'w')
    cons_file.write(tree_str)
    cons_file.close()
    if os.path.exists("%s/raxml_phylos/RAxML_info.%s_%s.tree" % (outdir, rootname, windex)):
        os.remove("%s/raxml_phylos/RAxML_info.%s_%s.tree" % (outdir, rootname, windex))
    old_dir = os.getcwd()
    os.chdir("%s/raxml_phylos" % outdir)
    cmd = ["%s/raxmlHPC-SSE3" % raxml_path, '-f', 'a', '-x', '12345', '-p', '12345', '-m', 'GTRGAMMA', '-#', '10', '-s', "%s_%s.afa" % (rootname, windex), '-n', "%s_%s.tree" % (rootname, windex), "-g", "%s/raxml_phylos/%s_%s.constraint" % (outdir, rootname, windex)]
    subprocess.call(cmd)
    os.chdir(old_dir)
    os.remove("%s/raxml_phylos/%s_%s.afa" % (outdir, rootname, windex))
    if os.path.exists("%s/raxml_phylos/%s_%s.afa.reduced" % (outdir, rootname, windex)):
        os.remove("%s/raxml_phylos/%s_%s.afa.reduced" % (outdir, rootname, windex))
    os.remove("%s/raxml_phylos/%s_%s.constraint" % (outdir, rootname, windex))
    os.remove("%s/raxml_phylos/RAxML_bipartitions.%s_%s.tree" % (outdir, rootname, windex))
    os.remove("%s/raxml_phylos/RAxML_bipartitionsBranchLabels.%s_%s.tree" % (outdir, rootname, windex))
    os.remove("%s/raxml_phylos/RAxML_bootstrap.%s_%s.tree" % (outdir, rootname, windex))
    os.remove("%s/raxml_phylos/RAxML_info.%s_%s.tree" % (outdir, rootname, windex))
    tree = PhyloTree("%s/raxml_phylos/RAxML_bestTree.%s_%s.tree" % (outdir, rootname, windex))
    keep_tips = []
    for leaf in tree:
        if leaf.dist <= 0.4:
            keep_tips.append(leaf.name)
    if len(keep_tips) >= min_taxa: 
        tree.prune(keep_tips)
        return tree.write(format = 5)
    else:
        return None


def baseml_blengths(ncar_list, aligndir, outdir, treefile, num_threads, remove_list):
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    pool = multiprocessing.Pool(processes = num_threads)
    work_list = []
    for ncar in ncar_list:
        work_list.append([ncar, aligndir, outdir, treefile, remove_list])
    pool.map_async(baseml_worker, work_list).get(9999999)

def baseml_worker(param_list):
    ncar = param_list[0]
    aligndir = param_list[1]
    outdir = param_list[2]
    treefile = param_list[3]
    remove_list = param_list[4]
    reader = SeqIO.parse("%s/%s.afa" % (aligndir, ncar), format = 'fasta')
    seq_dic = {}
    for rec in reader:
        seq_dic[rec.id] = str(rec.seq)
    outfile = open("%s/%s.afa" % (outdir, ncar), 'w')
    taxa_list = []
    for seq_name, seq in seq_dic.items():
        cur_taxa = seq_name.split(":")[0]
        if cur_taxa not in remove_list:
            taxa_list.append(cur_taxa)
            outfile.write(">%s\n%s\n" % (seq_name, seq))
    outfile.close()
    tree = PhyloTree(treefile)
    tree.prune(taxa_list)
    tree.unroot()
    tree_str = tree.write(format = 5)
    cons_file = open("%s/%s.constraint" % (outdir, ncar), 'w')
    cons_file.write(tree_str)
    cons_file.close()
    cml = baseml.Baseml(alignment = "%s/%s.afa" % (outdir, ncar), tree = "%s/%s.constraint" % (outdir, ncar), out_file = "%s/%s.alt" % (outdir, ncar), working_dir = "%s/%s_working" % (outdir, ncar))
    cml.set_options(runmode=0,fix_blength=0,model=7, clock = 0, Mgene = 0, fix_kappa = 0, kappa = 2, getSE = 0, RateAncestor = 0, cleandata = 0, Small_Diff = .45e-6, verbose = True)
    cml.run(command = "%s/baseml" % paml_path, verbose = True)

def read_baseml_phylos(ncar_list, indir, outdir, outfile):
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    output = open("%s/%s" % (outdir, outfile), 'w')
    for ncar in ncar_list:
        cur_blens = read_baseml_blengths("%s/%s.alt" % (indir, ncar))
        output.write("%s\t%s\n" % (ncar, cur_blens))
    output.close()

def read_baseml_blengths(infile):
    og_file =open(infile, 'rU')
    aa_tree = False
    first_line = True
    line = og_file.readline()
    while True:
        if first_line:
            seq_len = int(line.strip().split()[1])
            first_line = False
        if aa_tree:
            aa = PhyloTree(line.strip().replace("#", ":"))
            return aa.write(format = 5)
        if line.strip().startswith("tree length = "):
            og_file.readline()
            og_file.readline()
            og_file.readline()
            line = og_file.readline()
            aa_tree = True
            continue
        line = og_file.readline()
        if not line:
            break

def aaml_blengths(seq_dic, outdir, rootname, windex, window_size, min_taxa, constraint_tree):
    seq_dic = trimal(seq_dic, outdir, rootname, windex, window_size)
    if len(seq_dic) < min_taxa:
        return
    
    seq_file = open("%s/raxml_phylos/%s_%s.afa" % (outdir, rootname, windex), 'w')
    for k, v in seq_dic.items():
        cur_length = len(v)
        break
    seq_file.write("%s %s\n" % (len(seq_dic), cur_length))
    species_list = []
    for k, v in seq_dic.items():
        seq_file.write("%s\n%s\n" % (k, v))
        species_list.append(k)
    seq_file.close()
    tree = PhyloTree(constraint_tree)    
    tree.prune(species_list)
    tree_str = tree.write(format = 5)
    cons_file = open("%s/raxml_phylos/%s_%s.constraint" % (outdir, rootname, windex), 'w')
    cons_file.write(tree_str)
    cons_file.close()

    cml = baseml.Baseml(alignment = "%s/raxml_phylos/%s_%s.afa" % (outdir, rootname, windex), tree = "%s/raxml_phylos/%s_%s.constraint" % (outdir, rootname, windex), out_file = "%s/raxml_phylos/%s_%s.alt" % (outdir, rootname, windex), working_dir = "%s/raxml_phylos/%s_%s_working" % (outdir, rootname, windex))
    cml.set_options(runmode=0,fix_blength=0,model=7, clock = 0, Mgene = 0, fix_kappa = 0, kappa = 2, getSE = 0, RateAncestor = 0, cleandata = 0, Small_Diff = .45e-6, verbose = True)
    cml.run(command = "%s/baseml" % paml_path, verbose = True)
    cur_blens = read_aaml_blengths("%s/raxml_phylos/%s_%s.alt" % (outdir, rootname, windex), min_taxa)
    return cur_blens


def read_aaml_blengths(aaml_file, min_taxa):
    og_file =open(aaml_file, 'rU')
    aa_tree = False
    first_line = True
    line = og_file.readline()
    while True:
        if first_line:
            seq_len = int(line.strip().split()[1])
            first_line = False
        if aa_tree:
            tree = PhyloTree(line.strip())

            keep_tips = []
            for leaf in tree:
                if leaf.dist <= 0.4:
                    keep_tips.append(leaf.name)
            if len(keep_tips) >= min_taxa: 
                tree.prune(keep_tips)
                return tree.write(format = 5)
            else:
                return None

            return aa.write(format = 5)
            aa_tree = False
            break
        if line.strip().startswith("tree length = "):
            og_file.readline()
            og_file.readline()
            og_file.readline()
            line = og_file.readline()
            aa_tree = True
            continue
        line = og_file.readline()
        if not line:
            break



def trimal(seq_dic, outdir, rootname, windex, window_size):
    seq_file = open("%s/trimal_windows/%s_%s.dirty" % (outdir, rootname, windex), 'w')
    for k, v in seq_dic.items():
#        seq_file.write(">%s\n%s\n" % (k.split(":")[0], v))
            
        seq_file.write(">%s\n%s\n" % (k.replace(":","."), v))
    seq_file.close()
    inseq = "%s/trimal_windows/%s_%s.dirty" % (outdir, rootname, windex)
    outseq = "%s/trimal_windows/%s_%s.trimal" % (outdir, rootname, windex)
    cmd = ["%s/trimal" % trimal_path, "-automated1", "-in", inseq, "-out", outseq]
    subprocess.call(cmd)
    reader = SeqIO.parse(outseq, format = 'fasta')
    trim_dic = {}
    for rec in reader:
        cur_seq = str(rec.seq)
        if cur_seq.count("G") + cur_seq.count("C") + cur_seq.count("A") + cur_seq.count("T") > 0.5* window_size:
            trim_dic[rec.id.replace(".",":")] = str(rec.seq)
    return trim_dic

def add_species_to_maf(inspecies, newspecies, old_maf_list, new_maf_list):
    species_merged_maf_list = []
    for new_maf in new_maf_list:
        
        added = False
        if not new_maf:
            sys.exit()
        for old_maf in old_maf_list:
            new_start = new_maf.maf1.start
            old_start = old_maf.maf1.start
            new_end = new_maf.maf1.end
            old_end = old_maf.maf1.end
            new_start_lift = 0
            new_end_lift = len(new_maf.maf1.seq)

            if new_start >= old_start and new_start < old_end:
                new_start_lift = liftover(new_start - old_start, new_maf.maf1.seq, 0)
                if new_end > old_end:
                    new_end_lift = liftover(old_end-new_start, new_maf.maf1.seq, 0)
                maf2 = new_maf.maf2
                new_start_coord = maf2.start + new_start_lift - maf2.seq[:new_start_lift].count("-")
                new_end_coord = maf2.start + new_end_lift - maf2.seq[:new_end_lift].count("-")
                if old_maf.maf1.strand == new_maf.maf1.strand:
                    merged_maf = MAFSeq(maf2.species, maf2.scaf, maf2.start, maf2.orig_len, FILL_CHAR*new_start_lift + maf2.seq[:new_end_lift], maf2.strand, maf2.scaf_len, maf2.end, new_maf.maf1.start, new_maf.maf1.end, new_maf.maf1.seq) 
                else:
                    print "STRAND PROBLEM"
                    sys.exit()
                old_maf.other_merge(merged_maf)
                species_merged_maf_list.append(old_maf)
                added = True

            elif new_start <= old_start and new_end > old_start:
                new_start_lift = liftover(old_start - new_start, new_maf.maf1.seq, 0)
                if new_end > old_end:
                    new_end_lift = liftover(old_end-new_start, new_maf.maf1.seq, 0)
                    
                maf2 = new_maf.maf2
                new_start_coord = maf2.start + new_start_lift - maf2.seq[:new_start_lift].count("-")
                new_end_coord = maf2.start + new_end_lift - maf2.seq[:new_end_lift].count("-")
                if old_maf.maf1.strand == new_maf.maf1.strand:
                    merged_maf = MAFSeq(maf2.species, maf2.scaf, maf2.start, maf2.orig_len, maf2.seq[new_start_lift:new_end_lift], maf2.strand, maf2.scaf_len, maf2.end, new_maf.maf1.start, new_maf.maf1.end, new_maf.maf1.seq)
                else:
                    print "STRAND PROBLEM"
                    sys.exit()

                    merged_maf = MAFSeq(maf2.species, maf2.scaf, new_start_coord, maf2.orig_len, maf2.seq[new_start_lift:new_end_lift], maf2.strand, maf2.scaf_len, new_end_coord, new_maf.maf1.start, new_maf.maf1.end)

                old_maf.other_merge(merged_maf)
                if not added:
                    species_merged_maf_list.append(old_maf)
                added = True
        if not added:
            species_merged_maf_list.append(new_maf)
    return species_merged_maf_list
                    
def mask_unaligned_genomes(genome_seqs, gff_dir):
    masked_genomes = copy.deepcopy(genome_seqs)
    for species, seq_dic in genome_seqs.items():
        reader = open("%s/%s.gff" % (gff_dir, species), 'rU')
        cds_tuples = {}
        for line in reader:
            cur_line = line.split()
            cur_type = cur_line[2]
            if cur_type != "CDS":
                continue
            cur_scaf = cur_line[0]
            cur_start = int(cur_line[3])
            cur_end = int(cur_line[4])
            if cds_tuples.get(cur_scaf, None) == None:
                cds_tuples[cur_scaf] = []
            cds_tuples[cur_scaf].append((cur_start, cur_end))
        for scaf, cds_list in cds_tuples.items():
            scaf_seq = masked_genomes[species][scaf]
            scaf_seq_list = list(scaf_seq)
            for cds_coords in cds_list:
                cur_index = cds_coords[0] - 1
                while cur_index < cds_coords[1]:
                    scaf_seq_list[cur_index] = "n"
                    cur_index += 1
            masked_scaf = "".join(scaf_seq_list)
            masked_genomes[species][scaf] = masked_scaf
    return masked_genomes
