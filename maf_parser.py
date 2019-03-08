import cPickle as pickle
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
parser.add_option("--taxa_inclusion", dest = "taxa_inclusion", type = str)
parser.add_option("--outputfile", dest = "outputfile", type = str)
parser.add_option("-a", "--action", dest = "action", type = str)
parser.add_option("-r", "--params", dest = "params", type = str)

(options, args) = parser.parse_args()
outspecies_list = options.outspecies.split(",")

def main():
    if options.action == "build":
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
        pickle_file = open("%s/maf_dic.pickle" % (options.output_dir), 'wb')
        pickle.dump(new_maf_scaf_dic, pickle_file)
        pickle_file.close()
        sys.exit()

    elif options.action == "realign":
        pickle_file = open("%s/maf_dic.pickle" % (options.output_dir), 'rb')
        new_maf_scaf_dic = pickle.load(pickle_file)
        pickle_file.close()
        utils.realign_fsa(new_maf_scaf_dic, options.output_dir, options.inspecies, outspecies_list[0], options.num_threads)
        sys.exit()

    elif options.action == "write_loci":
        pickle_file = open("%s/maf_dic.pickle" % (options.output_dir), 'rb')
        new_maf_scaf_dic = pickle.load(pickle_file)
        pickle_file.close()
        utils.genome_blast_dbs(options.params, options.output_dir)
        print "write_loci"
        utils.slide_baby_slide(new_maf_scaf_dic, options.output_dir, options.inspecies, outspecies_list[0], options.window_size, options.step_size, options.min_taxa, options.constraint_tree, options.num_threads, options.params, options.gff_dir)
        sys.exit()

    elif options.action == "rer_converge":
        manda_taxa, multi_taxa, remove_list = utils.make_taxa_dic(options.taxa_inclusion)
        ncar_list = utils.min_taxa_membership(manda_taxa, multi_taxa, remove_list, "%s/filtered_loci.index" % (options.output_dir), options.min_taxa)
        utils.baseml_blengths(ncar_list, "%s/filtered_loci" % (options.output_dir), "%s/baseml_%s" % (options.output_dir, options.outputfile.split(".")[0]), options.constraint_tree, options.num_threads, remove_list)
        utils.read_baseml_phylos(ncar_list, "%s/baseml_%s" % (options.output_dir, options.outputfile.split(".")[0]), "%s/baseml_compiled" % (options.output_dir), options.outputfile)
        sys.exit()
        
    
if __name__ == '__main__':
    main()
