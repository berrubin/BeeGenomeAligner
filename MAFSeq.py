class MAFSeq:
    def __init__(self, species, scaf, start, orig_len, seq, strand, scaf_len, end, other_start = -1, other_end = -1, other_seq = -1):

        self.species = species
        self.scaf = scaf
        self.strand = strand
        self.orig_len = orig_len
        self.seq = seq
        self.realigned_seq = ""
        self.scaf_len = scaf_len
        self.start = start
        self.end = end
        self.other_start = other_start
        self.other_end = other_end
        self.other_seq = other_seq
#        if strand == "+":
#            self.start = start
#            self.end = start + orig_len
##        elif strand == "-":
#            self.end = scaf_len - start
#            self.start = self.end - orig_len


    def __str__(self):
        return "%s\t%s\t%s\t%s\t%s\t%s" % (self.species, self.scaf, self.start, self.end, self.orig_len, self.strand)

    def equals(self, other_maf):
        if self.species == other_maf.species:
            if self.scaf == other_maf.scaf:
                if self.start == other_maf.start:
                    if self.end == other_maf.end:
                        return True
        return False
