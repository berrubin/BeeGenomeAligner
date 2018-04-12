import sys
import utils

FILL_CHAR ="N"

class MAFPair:
    def __init__(self, species1, maf1, species2, maf2):
        self.species1 = species1
        self.maf1 = maf1
        self.species2 = species2
        self.maf2 = maf2
        self.distances = []
        self.framed_distances = []
        self.mafothers = []

    def other_merge(self, new_other):
        merged = False
#        if len(self.mafothers) > 0:
#            print "OTHERLIST"
#            for other in self.mafothers:
#                print other
        for other in self.mafothers:
            if other.species == new_other.species and other.scaf == new_other.scaf and other.strand == new_other.strand:
#                print "OTHER"
#                print self.maf1
#                print self.maf2
#                print other
#                print new_other
                if new_other.start == other.start and new_other.end == other.end:
                    return
                #my instinct is that all of these gaps should be N's to help correctly justify the alignments. however, FSA doesn't align N's with other stuff so that doesn't work at all.
#                print new_other
#                print other
                if new_other.strand == "+":
                    if new_other.start >= other.end:
                        position = len(new_other.seq) - len(new_other.seq.lstrip("N"))
                        new_seq = other.seq + FILL_CHAR*(position - other.end) + new_other.seq.lstrip("N")
                        new_other_seq = other.other_seq + FILL_CHAR*(position - other.other_end) + new_other.other_seq.lstrip("N")
                        other.end = new_other.end
                        other.other_end = new_other.other_end
                    elif other.start >= new_other.end:
                        position = len(other.seq) - len(other.seq.lstrip("N"))
                        new_seq = new_other.seq + FILL_CHAR*(position - new_other.end) + other.seq.lstrip("N")
                        new_other_seq = new_other.other_seq + FILL_CHAR*(position - new_other.other_end) + other.other_seq.lstrip("N")
                        other.start = new_other.start
                        other.other_start = new_other.other_start
                    elif new_other.start > other.start and new_other.end < other.end:
                        new_seq = other.seq[:new_other.start] + new_other.seq + other.seq[new_other.end:]
                        new_other_seq = other.other_seq[:new_other.other_start] + new_other.other_seq + other.other_seq[new_other.other_end:]
                else:
                    if new_other.start >= other.end:
                        position = len(new_other.seq) - len(new_other.seq.lstrip("N"))
                        new_seq = other.seq + FILL_CHAR*(position - other.end) + new_other.seq.lstrip("N")
                        new_other_seq = other.other_seq + FILL_CHAR*(position - other.other_end) + new_other.other_seq.lstrip("N")
                        other.end = new_other.end
                        other.other_end = new_other.other_end
                    elif other.start >= new_other.end:
                        position = len(other.seq) - len(other.seq.lstrip("N"))
                        new_seq = new_other.seq + FILL_CHAR*(position - new_other.end) + other.seq.lstrip("N")
                        new_other_seq = new_other.other_seq + FILL_CHAR*(position - new_other.other_end) + other.other_seq.lstrip("N")
                        other.start = new_other.start
                        other.other_start = new_other.other_start
                    elif new_other.start > other.start and new_other.end < other.end:
                        new_seq = other.seq[:new_other.start] + new_other.seq + other.seq[new_other.end:]
                        new_other_seq = other.other_seq[:new_other.other_start] + new_other.other_seq + other.other_seq[new_other.other_end:]
                    
                other.seq = new_seq
                other.other_seq = new_other_seq
                merged = True
                break
        if not merged:
            self.mafothers.append(new_other)
                

    def get_maf(self, species):
        if species == self.species1:
            return self.maf1
        elif species == self.species2:
            return self.maf2
        else:
            return self.maf2
#            return False
    
    def replace_maf(self, species, new_maf):
        if self.species1 == species:
            self.maf1 = new_maf
        elif self.species2 == species:
            self.maf2 = new_maf

    def less_than(self, other_mafpair):
        if self.maf1.start < other_mafpair.get_maf(self.species1).start:
            return True
        return False
        
    def __str__(self):
        return str(self.maf1) + "\n" + str(self.maf2)

    def sliding_distance(self, window_size, window_step, inspecies, outspecies):
        x = 0
        target_seq = self.get_maf(inspecies).seq
        out_seq = self.get_maf(outspecies).seq
        if len(target_seq) != len(out_seq):
            print len(target_seq)
            print len(out_seq)
            print target_seq
            print out_seq
            print self.maf1
            print self.maf2
            return []
            sys.exit()
        target_dist_list = []
        while x < len(target_seq):
            target_window = target_seq[x:x+window_size]
            out_window = out_seq[x:x+window_size]
            cur_dist = utils.distance(target_window, out_window)
            if cur_dist[1] > 0.5*window_size:
                target_dist_list.append((x, cur_dist[0] / cur_dist[1]))
            x = x + window_step
        target_coords = []
        out_coords = []
        for coord_tuple in target_dist_list:
            target_coord, target_scaf = self.coordinate_in_frame(inspecies, coord_tuple[0])
            out_coord, out_scaf = self.coordinate_in_frame(outspecies, coord_tuple[0])
            target_coords.append((target_coord, out_coord, coord_tuple[1]))

#            out_coords.append(out_coord, coord_tuple[1])
        self.distances = target_coords
        return target_coords

    def coordinate_in_frame(self, frame_species, coordinate):
#        print frame_species
#        print str(self)
        frame_maf = self.get_maf(frame_species)
        if frame_maf.strand == "+":
            new_coord = frame_maf.start + coordinate - frame_maf.seq[:coordinate].count("-")
        elif frame_maf.strand == "-":
#            new_coord = frame_maf.scaf_len - frame_maf.end + coordinate - frame_maf.seq[:coordinate].count("-")
                
            new_coord = frame_maf.end - coordinate + frame_maf.seq[:coordinate].count("-")
#            new_coord = frame_maf.scaf_len - (frame_maf.start + coordinate - frame_maf.seq[:coordinate].count("-"))

        return new_coord, frame_maf.scaf

    def transfer_coords(self, frame_list, inspecies, framespecies):
        new_distances = []
        for dist_tuple in self.distances:
            cur_coord = dist_tuple[0]
#            print cur_coord
            for frame_maf in frame_list:
                target_maf = frame_maf.get_maf(inspecies)
#                print str(self)
#                print target_maf
                gap_count = target_maf.seq[:cur_coord].count("-")
                if cur_coord >= target_maf.start and cur_coord <= target_maf.end:
                    coordinate = utils.liftover(cur_coord - target_maf.start, target_maf.seq, cur_coord - target_maf.start)
                    new_coord, new_scaf = frame_maf.coordinate_in_frame(framespecies, coordinate)
#                    if target_maf.start == 8134507:
#                        print frame_maf
#                        print cur_coord
#                        print coordinate
#                        print new_coord
#                    print new_coord
#                    print new_scaf
                    
                    new_distances.append((cur_coord, dist_tuple[1], dist_tuple[2], new_scaf, new_coord))
                    break
#        print new_distances
        self.framed_distances = new_distances
        return new_distances

    def reconstitute_coord(inspecies, frame_species, coordinate):
        target_maf = self.get_maf(inspecies)
#        coordinate = utils.liftover(coordinate, target_maf)
        frame_maf = self.get_maf(frame_species)
        if frame_maf.strand == "+":
            new_coord = frame_maf.start + coordinate - frame_maf.seq[:coordinate].count("-")
        elif frame_maf.strand == "-":
#            new_coord = frame_maf.scaf_len - frame_maf.end + coordinate - frame_maf.seq[:coordinate].count("-")
                
            new_coord = frame_maf.end - coordinate + frame_maf.seq[:coordinate].count("-")
            if frame_maf.start == 3645435:
                print str(self)
                print coordinate
                print frame_maf.seq[:coordinate].count("-")
                print new_coord
#            new_coord = frame_maf.scaf_len - (frame_maf.start + coordinate - frame_maf.seq[:coordinate].count("-"))

        return new_coord, frame_maf.scaf
