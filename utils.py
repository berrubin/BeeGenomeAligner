import os
from MAFSeq import MAFSeq
from MAFPair import MAFPair
import numpy
import copy
from Bio import SeqIO

HIT_OVERLAP = 0
INTRON_SIZE = 10000

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

#        print cur_hsp
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
#        print i
#        print my_hsp
#        print cur_hsp
#        if my_hsp_target.equals(cur_hsp_target):
#            i += 1
#            continue
        if my_hsp_out.strand == "+":
            if my_hsp_target.start <= cur_hsp_target.start and my_hsp_target.end >= cur_hsp_target.start - border_size:
                if my_hsp_target.end >= cur_hsp_target.end:
#                if False:
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
                        target_seq = my_hsp_target.seq + "N"*blank_len + "-"*target_gaps + cur_hsp_target.seq
                        out_seq = my_hsp_out.seq + "N"*out_blank_len + "-"*out_gaps + cur_hsp_out.seq 
                        target_maf = MAFSeq(my_hsp_target.species, my_hsp_target.scaf, my_hsp_target.start, cur_hsp_target.end - my_hsp_target.start, target_seq, my_hsp_target.strand, my_hsp_target.scaf_len, cur_hsp_target.end)
#                        print target_maf
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
                    target_seq = cur_hsp_target.seq + "N"*blank_len + "-"*target_gaps + my_hsp_target.seq
                    out_seq = cur_hsp_out.seq + "N"*out_blank_len + "-"*out_gaps + my_hsp_out.seq 
                    target_maf = MAFSeq(my_hsp_target.species, my_hsp_target.scaf, cur_hsp_target.start, my_hsp_target.end - cur_hsp_target.start, target_seq, my_hsp_target.strand, my_hsp_target.scaf_len, my_hsp_target.end)
#                    print target_maf
                    out_maf = MAFSeq(my_hsp_out.species, my_hsp_out.scaf, cur_hsp_out.start, my_hsp_out.end - cur_hsp_out.start, out_seq, my_hsp_out.strand, my_hsp_out.scaf_len, my_hsp_out.end)
                    new_hsp = MAFPair(target_maf.species, target_maf, out_maf.species, out_maf)

                    my_hsp = new_hsp
                    add_hsp = True
                    insert_pos = i
                    del temp_hsp_list[i]
                    i = -1
#            else:
#                break
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
                        target_seq = my_hsp_target.seq + "N"*blank_len + "-"*target_gaps + cur_hsp_target.seq
                        out_seq = my_hsp_out.seq + "N"*out_blank_len + "-"*out_gaps + cur_hsp_out.seq 
                        target_maf = MAFSeq(my_hsp_target.species, my_hsp_target.scaf, my_hsp_target.start, cur_hsp_target.end - my_hsp_target.start, target_seq, my_hsp_target.strand, my_hsp_target.scaf_len, cur_hsp_target.end)
#                        print target_maf

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
                    target_seq = cur_hsp_target.seq + "N"*blank_len + "-"*target_gaps + my_hsp_target.seq
                    out_seq = cur_hsp_out.seq + "N"*out_blank_len + "-"*out_gaps + my_hsp_out.seq 
                    target_maf = MAFSeq(my_hsp_target.species, my_hsp_target.scaf, cur_hsp_target.start, my_hsp_target.end - cur_hsp_target.start, target_seq, my_hsp_target.strand, my_hsp_target.scaf_len, my_hsp_target.end)
#                    print target_maf
                    out_maf = MAFSeq(my_hsp_out.species, my_hsp_out.scaf, my_hsp_out.start, cur_hsp_out.end - my_hsp_out.start, out_seq, cur_hsp_out.strand, my_hsp_out.scaf_len, cur_hsp_out.end)
                    new_hsp = MAFPair(target_maf.species, target_maf, out_maf.species, out_maf)
                    my_hsp = new_hsp
                    add_hsp = True
                    insert_post = i
                    del temp_hsp_list[i]
                    i = -1                
#            else:
#                break
#        print i
#        if my_hsp_target.start > cur_hsp_target.start and my_hsp_target.end < cur_hsp_target.end:
#            included = True
#            break
        i += 1
#    if add_hsp:
#    print "merge"
#    print insert_pos
#    print my_hsp
#    if not included:
    temp_hsp_list.append(my_hsp)
#        add_hsp = False
    return temp_hsp_list


def overlap(mytuple, tuplelist):
    included = False
    i = 0
    merged = []
    while i < len(tuplelist):
        cur_tuple = tuplelist[i]
        if mytuple[0] <= cur_tuple[0] and mytuple[1] >= cur_tuple[0]-1000:
#            merged.append(mytuple)
#            merged.append(cur_tuple)
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
#        print x
#        print len(seq1)
#        print len(seq2)

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
#        target_maf = mafpair.get_maf(inspecies)
#        out_maf = mafpair.get_maf(outspecies)
#        if target_maf.scaf not in target_scaf_dists.keys():
#            target_scaf_dists[target_maf.scaf] = []
#        target_scaf_dists[target_maf.scaf] = target_scaf_dists[target_maf.scaf] + in_tuples
#        if out_maf.scaf not in out_scaf_dists.keys():
#            out_scaf_dists[out_maf.scaf] = []
#        out_scaf_dists[out_maf.scaf] = out_scaf_dists[out_maf.scaf] + out_tuples
#    return target_scaf_dists, out_scaf_dists

def sliding_window_transfer_frame(maf_dic, frame_maf_dic, window_size, window_step, inspecies, outspecies, framespecies, output_dir):
    target_scaf_dists = {}
    out_scaf_dists = {}
    outfile = open("%s/%s_vs_%s_%s_%s_distances_frame_%s.txt" % (output_dir, inspecies, outspecies, window_size, window_step, framespecies), 'w')    
    outfile.write("inspecies\tinscaf\tincoord\toutspecies\toutscaf\toutcoord\tdxy\tframescaf\tframecoord\n")
    for scaf, maf_list in maf_dic.items():
        for mafpair in maf_list:
            distance_tuples = mafpair.sliding_distance(window_size, window_step, inspecies, outspecies)
        for mafpair in maf_list:
#            inmaf = mafpair.get_maf(inspecies)
            mafpair.transfer_coords(frame_maf_dic[scaf], inspecies, framespecies)
        for mafpair in maf_list:
#            print mafpair
#            print mafpair.framed_distances
            inmaf = mafpair.get_maf(inspecies)
#            if inmaf.start == 3645435:
#                print inmaf
#                print inmaf.framed_distances
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
    reader = open(inmaf, 'rU')
    maf_blocks = []
    seq_count = 0
    for line in reader:
        if line.startswith("#"):
            continue
        if line.startswith("a"):
            if seq_count > 0:
#                print seq1.species
#                print seq2.species
#                if seq1.start > 1000000 or seq1.start < 950000:
#                    seq_count = 0
#                    continue
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
#            cur_species = cur_name[0]
#            cur_scaf = cur_name[1]
            if seq_count == 0:
                cur_species = inspecies
            else:
                cur_species = outspecies
            cur_scaf = cur_line[1]
            cur_start = int(cur_line[2])
            cur_len = int(cur_line[3])
            cur_strand = cur_line[4]
            cur_scaf_len = int(cur_line[5])
            cur_seq = cur_line[6]

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
    print len(maf_blocks)
    return maf_blocks

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
        print syntenics
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
#    print hsp_list
#    print hsp_count
    while True:
        if len(hsp_list) == 0:
            break
        smallest_hsp = hsp_list[0]
#        print smallest_hsp
        for hsp in hsp_list:
            if hsp.less_than(smallest_hsp):
                smallest_hsp = hsp
        new_hsp_list.append(smallest_hsp)
        hsp_list.remove(smallest_hsp)
#        print new_hsp_list
#        print len(new_hsp_list)
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
#                print "contained"
#                print hsp
#                print contained_hsp
                break
        if keep_hsp:
            new_hsp_list.append(hsp)
#    for hsp in new_hsp_list:
#        print hsp
    return new_hsp_list

def merge_overlord(maf_list, inspecies, outspecies):
    
    x = 0
    temp_maf_list = copy.deepcopy(maf_list)
    cur_maf_len = len(maf_list)
    while x < cur_maf_len:
        temp_maf_list = merge_HSPs(temp_maf_list[x], temp_maf_list, inspecies, outspecies)
        if len(temp_maf_list) < cur_maf_len:
            x = 0
#            temp_maf_list = sort_hsps(temp_maf_list)
            cur_maf_len = len(temp_maf_list)
#            print cur_maf_len
        else:
            x += 1

#    for hsp in temp_maf_list:
#        print hsp
    print len(maf_list)
    print len(temp_maf_list)
    temp_maf_list = remove_shorts(temp_maf_list, inspecies)
    temp_maf_list = sort_hsps(temp_maf_list)
    for block in temp_maf_list:
        print block
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
#                print cds
                if target_maf.strand == "+":
                    if cds[0] > target_maf.start and cds[0] < target_maf.end:
                        #                    print cds
                        start_coord = liftover(cds[0] - target_maf.start, target_maf.seq, cds[0] - target_maf.start)
                        end_coord = liftover(cds[1] - target_maf.start, target_maf.seq, start_coord)
#                        if target_maf.start == 4623570:
#                            print end_coord
                        if end_coord >= seq_len:
                            end_coord = seq_len -1
#                    print start_coord
                        insert_seq = target_maf.seq[start_coord:end_coord+1].upper()
#                        if target_maf.start == 4623570:
#                            print insert_seq

                        insert_seq = insert_seq.replace("A", "N").replace("C", "N").replace("T", "N").replace("G","N")
#                        if target_maf.start == 4623570:
#                            print insert_seq


#                        insert_seq = "N" * (end_coord - start_coord + 1)
                        target_maf.seq = target_maf.seq[:start_coord-1] + insert_seq + target_maf.seq[end_coord:]
                elif target_maf.strand == "-":
                    start_coord = (target_maf.scaf_len - cds[1]) - (target_maf.scaf_len - target_maf.end)
                    end_coord = (target_maf.scaf_len - cds[0]) - (target_maf.scaf_len - target_maf.end)
#                    print cds
#                    print start_coord
#                    print end_coord
#                    if start_coord < 0 or end_coord < 0:
#                        sys.exit()
                    if cds[0] > target_maf.start and cds[0] < target_maf.end:
                        #                    print cds
                        start_coord = liftover(start_coord, target_maf.seq, start_coord)
                        end_coord = liftover(end_coord, target_maf.seq, start_coord)
#                    print start_coord
                        if end_coord > seq_len:
                            end_coord = seq_len

                        insert_seq = target_maf.seq[start_coord:end_coord+1].upper()
                        insert_seq = insert_seq.replace("A", "N").replace("C", "N").replace("T", "N").replace("G","N")
#                        insert_seq = "N" * (end_coord - start_coord + 1)
#                        x = start_coord
#                        while x < end_coord:
#                            if target_maf.seq[x] == "-":
#                                insert_seq[x - start_coord] = "-"
                        target_maf.seq = target_maf.seq[:start_coord] + insert_seq + target_maf.seq[end_coord+1:]

#                else:
#                        start
#                print target_maf.seq
#        mafpair.replace_maf(species, target_maf)
    return maf_list
            

def liftover(target_index, gapped_seq, i_start):
    i = i_start
    prev_i = i
    prev_gap = gapped_seq[:i].count("-")
    while i < len(gapped_seq):
        cur_gap_count = prev_gap + gapped_seq[prev_i:i].count("-")
        prev_gap = cur_gap_count
        prev_i = i
#        cur_gap_count = gapped_seq[:i].count("-")
#        print i
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
    for maf in maf_list:
        target_maf = maf.get_maf(inspecies)
        out_maf = maf.get_maf(outspecies)
#        print target_maf
#        print out_maf
        outfile = open("%s/%s_%s_%s_%s_%s_%s.fa" % (outdir, target_maf.scaf, target_maf.start, target_maf.end, out_maf.scaf, out_maf.start, out_maf.end), 'w')
        outfile.write(">%s:%s:%s\n%s\n" % (target_maf.scaf, target_maf.start, target_maf.end, target_maf.seq))
        outfile.write(">%s:%s:%s\n%s\n" % (out_maf.scaf, out_maf.start, out_maf.end, out_maf.seq))
        outfile.close()

