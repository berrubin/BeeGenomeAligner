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
parser.add_option("-k", "--target_scafs", dest = "target_scafs", type = "str")

(options, args) = parser.parse_args()
outspecies_list = options.outspecies.split(",")

def main():
    if not os.path.isdir(options.output_dir):
        os.mkdir(options.output_dir)
    pairs_dic = {}
    in_masked = False
    pairs_dic = utils.read_mafs_overlord(options.output_dir, options.num_threads, outspecies_list, options.maf_dir, options.inspecies, options.gff_dir, options.target_scafs)
    nearest_dic = pairs_dic[outspecies_list[0]]
    new_maf_scaf_dic = {}
    for scaf, maf_list in nearest_dic.items():
        print scaf
        for outspecies in outspecies_list[1:]:
            if scaf in pairs_dic[outspecies].keys():
                utils.add_species_to_maf(options.inspecies, outspecies, maf_list, pairs_dic[outspecies][scaf])
        maf_list = utils.filter_small_mafs(maf_list, options.min_taxa)
        utils.write_mafs_to_file(maf_list, options.output_dir, options.inspecies, outspecies_list[0])
        new_maf_scaf_dic[scaf] = maf_list
    utils.realign_fsa(new_maf_scaf_dic, options.output_dir, options.inspecies, outspecies_list[0], options.num_threads)
    utils.slide_baby_slide(new_maf_scaf_dic, options.output_dir, options.inspecies, outspecies_list[0], options.window_size, options.step_size, options.min_taxa, options.constraint_tree, options.num_threads)
        
    
if __name__ == '__main__':
    main()
