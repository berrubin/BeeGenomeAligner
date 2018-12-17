import numpy
import copy
from Bio import SeqIO
import utils
import os
import sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i", "--maf_dir", dest = "maf_dir", type = "str")
parser.add_option("-f", "--gff_dir", dest = "gff_dir", type = "str")
parser.add_option("-s", "--inspecies", dest = "inspecies", type = "str")
parser.add_option("-t", "--outspecies", dest = "outspecies", type = "str")
parser.add_option("-w", "--window_size", dest = "window_size", type = "int", default = 200)
parser.add_option("-z", "--step_size", dest = "step_size", type = "int", default = 50)
parser.add_option("-m", "--min_taxa", dest = "min_taxa", type = "int", default = 10)
parser.add_option("-c", "--constraint_tree", dest = "constraint_tree", type = "str")
parser.add_option("-o", "--output_dir", dest = "output_dir", type = "str")
parser.add_option("-p", "--num_threads", dest = "num_threads", type = "int", default = 1)
parser.add_option("-k", "--target_scafs", dest = "target_scafs", type = "str", default = "AMEL_all_scafs.txt")

(options, args) = parser.parse_args()
outspecies_list = options.outspecies.split(",")

def main():
    if not os.path.isdir(options.output_dir):
        os.mkdir(options.output_dir)
    pairs_dic = {}
    in_masked = False
    pairs_dic = utils.read_mafs_overlord(options.output_dir, options.num_threads, outspecies_list, options.maf_dir, options.inspecies, options.gff_dir, options.target_scafs)
    nearest_dic = pairs_dic[outspecies_list[0]]
    # for outspecies in outspecies_list:
    #     maf_file = "%s/%s.%s.sing.maf" % (options.maf_dir, options.inspecies, outspecies)
    #     maf_list = utils.read_maf(maf_file, options.inspecies, outspecies)
    #     print len(maf_list)
    #     maf_dic = utils.mafs_by_scaf(maf_list, options.inspecies)
    #     print len(maf_dic)
    #     maf_dic = utils.only_syntenic(maf_dic, outspecies)
    #     print len(maf_dic)
    #     processed_maf_dic = {}
    #     for scaf, maf_list in maf_dic.items():
    #         maf_list = utils.sort_hsps(maf_list)
    #         print scaf
    #         print "maf len: " + str(len(maf_list))
    #         merged_maf_list = utils.merge_overlord(maf_list, options.inspecies, outspecies)
    #         print "merged maf len: " + str(len(merged_maf_list))
    #         merged_maf_list = utils.remove_contained(merged_maf_list, options.inspecies)
    #         ingff = "%s/%s.gff" % (options.gff_dir, options.inspecies)
    #         if not in_masked:
    #             maf_list = utils.mask_proteins(merged_maf_list, options.inspecies, ingff)
    #             in_masked = True
    #         outgff = "%s/%s.gff" % (options.gff_dir, outspecies)
    #         print "masked"
    #         maf_list = utils.mask_proteins(maf_list, outspecies, outgff)
    #         maf_list = merged_maf_list
    #         print "masked"
    #         processed_maf_dic[scaf] = maf_list
    #     pairs_dic[outspecies] = processed_maf_dic
    nearest_dic = pairs_dic[outspecies_list[0]]

    extend_pairs_dic = {}
    new_maf_scaf_dic = {}
    for scaf, maf_list in nearest_dic.items():
        print scaf
        print len(maf_list)
        for outspecies in outspecies_list[1:]:
            if scaf in pairs_dic[outspecies].keys():
                utils.add_species_to_maf(options.inspecies, outspecies, maf_list, pairs_dic[outspecies][scaf])
        print "write to file"
        maf_list = utils.filter_small_mafs(maf_list, options.min_taxa)
        utils.write_mafs_to_file(maf_list, options.output_dir, options.inspecies, outspecies_list[0])
        new_maf_scaf_dic[scaf] = maf_list
    utils.realign_fsa(new_maf_scaf_dic, options.output_dir, options.inspecies, outspecies_list[0], options.num_threads)
    utils.slide_baby_slide(new_maf_scaf_dic, options.output_dir, options.inspecies, outspecies_list[0], options.window_size, options.step_size, options.min_taxa, options.constraint_tree, options.num_threads)
        
    
if __name__ == '__main__':
    main()
