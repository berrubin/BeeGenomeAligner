import numpy
import copy
from Bio import SeqIO
import utils
import os
import sys
from optparse import OptionParser

parser = OptionParser()

parser.add_option("-i", "--input_maf", dest = "input_maf", type = "str")
parser.add_option("-g", "--inspecies_genome", dest = "inspecies_genome", type = "str")
parser.add_option("-j", "--outspecies_genome", dest = "outspecies_genome", type = "str")
parser.add_option("-s", "--inspecies", dest = "inspecies", type = "str")
parser.add_option("-t", "--outspecies", dest = "outspecies", type = "str")
parser.add_option("-p", "--inspecies_gff", dest = "inspecies_gff", type = "str")
parser.add_option("-r", "--outspecies_gff", dest = "outspecies_gff", type = "str")
parser.add_option("-w", "--window_size", dest = "window_size", type = "int", default = 200)
parser.add_option("-z", "--step_size", dest = "step_size", type = "int", default = 50)
parser.add_option("-f", "--frame_species", dest = "frame_species", type = "str")
parser.add_option("-m", "--frame_maf", dest = "frame_maf", type = "str")
parser.add_option("-e", "--frame_gff", dest = "frame_gff", type = "str")
parser.add_option("-o", "--output_dir", dest = "output_dir", type = "str")

(options, args) = parser.parse_args()


def main():
    if not os.path.isdir(options.output_dir):
        os.mkdir(options.output_dir)
    maf_list = utils.read_maf(options.input_maf, options.inspecies, options.outspecies)
    print len(maf_list)
    maf_dic = utils.mafs_by_scaf(maf_list, options.inspecies)
    print len(maf_dic)
    maf_dic = utils.only_syntenic(maf_dic, options.outspecies)
    print len(maf_dic)
#    frame_maf_list = utils.read_maf(options.frame_maf, options.inspecies, options.frame_species)
#    frame_maf_dic = utils.mafs_by_scaf(frame_maf_list, options.inspecies)
#    frame_maf_dic = utils.only_syntenic(frame_maf_dic, options.frame_species)
    processed_maf_dic = {}

    for scaf, maf_list in maf_dic.items():
        if scaf != "3R":
            continue
        maf_list = utils.sort_hsps(maf_list)
        print scaf
        print "maf len: " + str(len(maf_list))
        merged_maf_list = utils.merge_overlord(maf_list, options.inspecies, options.outspecies)
        print "merged maf len: " + str(len(merged_maf_list))
        merged_maf_list = utils.remove_contained(merged_maf_list, options.inspecies)
        utils.write_mafs_to_file(merged_maf_list, options.output_dir, options.inspecies, options.outspecies)
        maf_list = utils.mask_proteins(merged_maf_list, options.inspecies, options.inspecies_gff)
        maf_list = utils.mask_proteins(merged_maf_list, options.outspecies, options.outspecies_gff)
        utils.write_mafs_to_file(maf_list, options.output_dir, options.inspecies, options.outspecies)
        processed_maf_dic[scaf] = maf_list
    
        utils.sliding_window(maf_list, options.window_size, options.step_size, options.inspecies, options.outspecies, options.output_dir, scaf)

    processed_frame_dic = {}
    for scaf, frame_list in frame_maf_dic.items():
        frame_list = utils.sort_hsps(frame_list)
        merged_frame_list = utils.merge_overlord(frame_list, options.inspecies, options.frame_species)
        merged_frame_list = utils.remove_contained(merged_frame_list, options.inspecies)
#        utils.write_mafs_to_file(merged_frame_list, options.output_dir, options.inspecies, options.frame_species)
        merged_frame_list = utils.mask_proteins(merged_frame_list, options.inspecies, options.inspecies_gff)
        merged_frame_list = utils.mask_proteins(merged_frame_list, options.frame_species, options.frame_gff)
        utils.write_mafs_to_file(merged_frame_list, options.output_dir, options.inspecies, options.frame_species)
        utils.sliding_window(merged_frame_list, options.window_size, options.step_size, options.inspecies, options.frame_species, options.output_dir, scaf)
        processed_frame_dic[scaf] = merged_frame_list

    print processed_frame_dic.keys()
    
    utils.sliding_window_transfer_frame(processed_maf_dic, processed_frame_dic, options.window_size, options.step_size, options.inspecies, options.outspecies, options.frame_species, options.output_dir)

if __name__ == '__main__':
    main()
