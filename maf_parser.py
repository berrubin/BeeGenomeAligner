import numpy
import copy
from Bio import SeqIO
import utils
import os
import sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i", "--maf_dir", dest = "maf_dir", type = "str")
#parser.add_option("-g", "--genome_dir", dest = "genome_dir", type = "str")
parser.add_option("-f", "--gff_dir", dest = "gff_dir", type = "str")
parser.add_option("-s", "--inspecies", dest = "inspecies", type = "str")
parser.add_option("-t", "--outspecies", dest = "outspecies", type = "str")
#parser.add_option("-p", "--inspecies_gff", dest = "inspecies_gff", type = "str")
#parser.add_option("-r", "--outspecies_gff", dest = "outspecies_gff", type = "str")
parser.add_option("-w", "--window_size", dest = "window_size", type = "int", default = 200)
parser.add_option("-z", "--step_size", dest = "step_size", type = "int", default = 50)
#parser.add_option("-f", "--frame_species", dest = "frame_species", type = "str")
#parser.add_option("-m", "--frame_maf", dest = "frame_maf", type = "str")
#parser.add_option("-e", "--frame_gff", dest = "frame_gff", type = "str")
parser.add_option("-m", "--min_taxa", dest = "min_taxa", type = "int", default = 10)
parser.add_option("-c", "--constraint_tree", dest = "constraint_tree", type = "str")
parser.add_option("-o", "--output_dir", dest = "output_dir", type = "str")
parser.add_option("-p", "--num_threads", dest = "num_threads", type = "int", default = 1)
parser.add_option("-k", "--target_scafs", dest = "target_scafs", type = "str", default = "/Genomics/kocherlab/berubin/alignment/10bees/genomes/AMEL_all_scafs.txt")

(options, args) = parser.parse_args()
outspecies_list = options.outspecies.split(",")

def main():
    if not os.path.isdir(options.output_dir):
        os.mkdir(options.output_dir)
    pairs_dic = {}
    in_masked = False
    pairs_dic = utils.read_mafs_overlord(options.output_dir, options.num_threads, outspecies_list, options.maf_dir, options.inspecies, options.gff_dir, options.target_scafs)
    nearest_dic = pairs_dic[outspecies_list[0]]
#     for outspecies in outspecies_list:
#         maf_file = "%s/%s.%s.sing.maf" % (options.maf_dir, options.inspecies, outspecies)
#         maf_list = utils.read_maf(maf_file, options.inspecies, outspecies)
#         print len(maf_list)
# #        maf_list = utils.restrict_region(maf_list, options.inspecies)
#         maf_dic = utils.mafs_by_scaf(maf_list, options.inspecies)
#         print len(maf_dic)
#         maf_dic = utils.only_syntenic(maf_dic, outspecies)
#         print len(maf_dic)
# #    frame_maf_list = utils.read_maf(options.frame_maf, options.inspecies, options.frame_species)
# #    frame_maf_dic = utils.mafs_by_scaf(frame_maf_list, options.inspecies)
# #    frame_maf_dic = utils.only_syntenic(frame_maf_dic, options.frame_species)
#         processed_maf_dic = {}
#         for scaf, maf_list in maf_dic.items():
#             if scaf != "Group10.1":
#                 continue
#             maf_list = utils.sort_hsps(maf_list)
#             print scaf
#             print "maf len: " + str(len(maf_list))
#             merged_maf_list = utils.merge_overlord(maf_list, options.inspecies, outspecies)
#             print "merged maf len: " + str(len(merged_maf_list))
#             merged_maf_list = utils.remove_contained(merged_maf_list, options.inspecies)
# #            utils.write_mafs_to_file(merged_maf_list, options.output_dir, options.inspecies, outspecies)
#             ingff = "%s/%s.gff" % (options.gff_dir, options.inspecies)
#             if not in_masked:
#                 maf_list = utils.mask_proteins(merged_maf_list, options.inspecies, ingff)
#                 in_masked = True
#             outgff = "%s/%s.gff" % (options.gff_dir, outspecies)
#             print "masked"
#             maf_list = utils.mask_proteins(maf_list, outspecies, outgff)
#             maf_list = merged_maf_list
#             print "masked"
# #            utils.write_mafs_to_file(maf_list, options.output_dir, options.inspecies, outspecies)
#             processed_maf_dic[scaf] = maf_list
#         pairs_dic[outspecies] = processed_maf_dic
#     nearest_dic = pairs_dic[outspecies_list[0]]

    extend_pairs_dic = {}
    new_maf_scaf_dic = {}
    for scaf, maf_list in nearest_dic.items():
        print scaf
        print len(maf_list)
        for outspecies in outspecies_list[1:]:
            if scaf in pairs_dic[outspecies].keys():
                utils.add_species_to_maf(options.inspecies, outspecies, maf_list, pairs_dic[outspecies][scaf])
        # print len(pairs_dic["BIMP"][scaf])
        
        # new_maf_list = utils.add_species_to_maf(options.inspecies, "EMEX", maf_list, pairs_dic["EMEX"][scaf])
        # for new_maf in maf_list:
        #     print new_maf
        #     for other in new_maf.mafothers:
        #         print "\t" + str(other)

        # new_maf_list = utils.add_species_to_maf(options.inspecies, "HLAB", maf_list, pairs_dic["HLAB"][scaf])
#        for new_maf in maf_list:
#            print new_maf
#            for other in new_maf.mafothers:
#                print "\t" + str(other)

#        print "print to screen"
#        for maf in new_maf_list:
#            print maf
#        maf_list = utils.restrict_region(maf_list, options.inspecies)
        print "write to file"
        maf_list = utils.filter_small_mafs(maf_list, options.min_taxa)        

#uncomment to run from beginning
#        utils.write_mafs_to_file(maf_list, options.output_dir, options.inspecies, outspecies_list[0])
        new_maf_scaf_dic[scaf] = maf_list
###    utils.blast_filter(new_maf_scaf_dic, options.output_dir, options.inspecies, outspecies_list[0], options.num_threads)
#uncomment to run from beginning
#    utils.realign_fsa(new_maf_scaf_dic, options.output_dir, options.inspecies, outspecies_list[0], options.num_threads)
    utils.slide_baby_slide(new_maf_scaf_dic, options.output_dir, options.inspecies, outspecies_list[0], options.window_size, options.step_size, options.min_taxa, options.constraint_tree, options.num_threads)
        
    
#            utils.sliding_window(maf_list, options.window_size, options.step_size, options.inspecies, options.outspecies, options.output_dir, scaf)

#     processed_frame_dic = {}
#     for scaf, frame_list in frame_maf_dic.items():
#         frame_list = utils.sort_hsps(frame_list)
#         merged_frame_list = utils.merge_overlord(frame_list, options.inspecies, options.frame_species)
#         merged_frame_list = utils.remove_contained(merged_frame_list, options.inspecies)
# #        utils.write_mafs_to_file(merged_frame_list, options.output_dir, options.inspecies, options.frame_species)
#         merged_frame_list = utils.mask_proteins(merged_frame_list, options.inspecies, options.inspecies_gff)
#         merged_frame_list = utils.mask_proteins(merged_frame_list, options.frame_species, options.frame_gff)
#         utils.write_mafs_to_file(merged_frame_list, options.output_dir, options.inspecies, options.frame_species)
#         utils.sliding_window(merged_frame_list, options.window_size, options.step_size, options.inspecies, options.frame_species, options.output_dir, scaf)
#         processed_frame_dic[scaf] = merged_frame_list

#     print processed_frame_dic.keys()
    
#     utils.sliding_window_transfer_frame(processed_maf_dic, processed_frame_dic, options.window_size, options.step_size, options.inspecies, options.outspecies, options.frame_species, options.output_dir)

if __name__ == '__main__':
    main()
