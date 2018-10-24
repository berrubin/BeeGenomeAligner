import random
from glob import glob
from Bio import SeqIO
import copy
from ete3 import PhyloTree
from optparse import OptionParser

parser = OptionParser()

parser.add_option("-i", "--input_file", dest = "infile")
parser.add_option("-o", "--output_file", dest = "outfile")
parser.add_option("-r", "--ref_gff", dest = "ref_gff")
parser.add_option("-p", "--promoter_size", dest = "promoter_size", type = int)
(options, args) = parser.parse_args()

def main():
    assign_gene()    

def assign_gene():
    gff_dic = read_gff(options.ref_gff)
    outfile = open(options.outfile, 'w')
    outfile.write("chr\tcoord\tproxtype\tgene\n")
    reader = open(options.infile, 'rU')
    phylo_dic = {}
    for line in reader:
        cur_line = line.split()
        cur_scaf = cur_line[0]
        if cur_scaf not in coord_dic.keys():
            coord_dic[cur_scaf] = []
        cur_start = int(cur_line[1])
        cur_proxim = process_coord(cur_start, cur_scaf, gff_dic)
        outfile.write("%s\t%s\t%s\t%s\n" % (cur_scaf, cur_start, cur_proxim[0], cur_proxim[1]))
    outfile.close()
            
def process_coord(start_coord, scaf, gene_dic):
    mod_name_dic = get_old_og_dic()
    if scaf not in gene_dic.keys():
        return "intergenic", "NA", "NA", "NA", "NA", "NA"
    added = False
    coord_proxim_dic = {}
    for gene in gene_dic[scaf]:
        cur_proxim = gene.proximity(start_coord)
        if not cur_proxim:
            continue
        coord_proxim_dic[gene.name] = cur_proxim
        added = True
    if added:
        best_proxim = choose_best_proxim(coord_proxim_dic)
        for gene in gene_dic[scaf]:
            if gene.name == best_proxim[0]:
                best_gene = gene
                break
        express_bias = best_gene.expression_bias
        wq_express_bias = best_gene.wq_expression_bias
        he_wq_express_bias = best_gene.he_wq_expression_bias
        proxtype = best_proxim[1]
        added = True
        return proxtype, best_gene.name
    if not added:
        return "intergenic", "NA"


def read_gff(gff_file):
    reader = open(gff_file, 'rU')
    ingene = False
    first_gene = True
    scaf_dic = {}
    for line in reader:
        cur_line = line.split()
        cur_scaf = cur_line[0]
        if cur_scaf not in scaf_dic.keys():
            scaf_dic[cur_scaf] = []
        if cur_line[2] == "gene":
            if not first_gene:
                cur_gene.get_introns()
                scaf_dic[cur_gene.scaf].append(cur_gene)
            first_gene = False
            cur_id = cur_line[-1].split("ID=")[1]
            cur_gene = Gene(cur_id, int(cur_line[3]), int(cur_line[4]), cur_line[6], cur_line[0])
        elif cur_line[2] == "CDS":
            cur_gene.add_cds(int(cur_line[3]), int(cur_line[4]))
        elif cur_line[2] == "five_prime_UTR":
            cur_gene.add_five(int(cur_line[3]), int(cur_line[4]))
        elif cur_line[2] == "three_prime_UTR":
            cur_gene.add_three(int(cur_line[3]), int(cur_line[4]))

    scaf_dic[cur_scaf].append(cur_gene)
    return scaf_dic

def choose_best_proxim(coord_proxim_dic):
    
    proxtype_counts = {}
    proxtype_genes = {}
    for proxtype in ["intron", "three_utr", "five_utr", "promoter", "stream"]:
        
        proxtype_counts[proxtype] = []
    best_per_gene = {}
    for gene_name, proxim in coord_proxim_dic.items():
        if proxim[0] in ["downstream", "upstream"]:
            proxtype_counts["stream"].append((gene_name, proxim[0], proxim[1])) 
        else:
            proxtype_counts[proxim[0]].append((gene_name, proxim[0], proxim[1]))
    for proxtype in ["intron", "three_utr", "five_utr", "promoter", "stream"]:
        prox_list = proxtype_counts[proxtype]
        if len(prox_list) == 0:
            continue
        if proxtype == "intron":
            return random.choice(prox_list)
        elif proxtype == "three_utr":
            return random.choice(prox_list)
        elif proxtype == "five_utr":
            return random.choice(prox_list)
        elif proxtype == "promoter":
            return choose_closest(prox_list)
        else:
            return choose_closest(prox_list)

def choose_closest(prox_list):
    closest_coord = prox_list[0][2]
    closest_prox = prox_list[0]
    for prox in prox_list:
        if prox[2] < closest_coord:
            closest_coord = prox[2]
            closest_prox = prox
    return closest_prox


class Gene:
    def __init__(self, gene_name, gene_start, gene_end, gene_strand, gene_scaf):
        self.name = gene_name
        self.start = gene_start
        self.end = gene_end
        self.strand = gene_strand
        self.scaf = gene_scaf
        self.cds = []
        self.five_utrs = []
        self.three_utrs = []
        self.introns = []
        self.expression_bias = "NA"
        self.wq_expression_bias = "NA"
        self.he_wq_expression_bias = "he_NA"

    def add_expression_bias(self, bias_towards):
        self.expression_bias = bias_towards

    def add_wq_expression_bias(self, bias_towards):
        self.wq_expression_bias = bias_towards

    def add_he_wq_expression_bias(self, bias_towards):
        self.he_wq_expression_bias = bias_towards
        
    def add_cds(self, cur_start, cur_end):
        self.cds.append((cur_start+1, cur_end-1))

    def add_five(self, cur_start, cur_end):
        self.five_utrs.append((cur_start-1, cur_end+1))

    def add_three(self, cur_start, cur_end):
        self.three_utrs.append((cur_start-1, cur_end+1))

    def get_introns(self):
        intron_list = []
        for x in range(len(self.cds)-1):
            if self.strand == "+":
                cur_intron = (self.cds[x][1], self.cds[x+1][0])
            else:
                cur_intron = (self.cds[x+1][0], self.cds[x][1])
            intron_list.append(cur_intron)
        self.introns = intron_list
        

    def coord_in(self, target_coord, target_list):
        for coords in target_list:
            if target_coord >= coords[0] and target_coord <= coords[1]:
                return True
            else:
                if target_coord >= coords[1] and target_coord <= coords[0]:
                    return True

    def proximity(self, target_coord):
        if target_coord < self.end + 10000 and target_coord > self.start - 10000:

            if self.coord_in(target_coord, self.introns):
                return ("intron", -9)
            elif self.coord_in(target_coord, self.five_utrs):
                return ("five_utr", -9)
            elif self.coord_in(target_coord, self.three_utrs):
                return ("three_utr", -9)
            elif self.strand == "+":
                if target_coord <= self.start:
                    if target_coord >= self.start - 1500:
                        return ("promoter", self.start - target_coord)
                    else:
                        return ("upstream", self.start - target_coord)
                else:
                    return ("downstream", target_coord - self.end)
            elif self.strand == "-":
                if target_coord >= self.end:
                    if target_coord <= self.end + 1500:
                        return ("promoter", target_coord - self.end)
                    else:
                        return ("upstream", target_coord - self.end)
                else:
                    return ("downstream", self.start - target_coord)

            elif self.coord_in(target_coord, self.cds):
                print self.name
                print target_coord
                return ("cds", -9)
        else:
            return False

if __name__ == '__main__':
    main()



