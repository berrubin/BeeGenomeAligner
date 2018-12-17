# README

Scripts for generating genome alignments used in Rubin et al. (2019)

## Getting Started

### Repeat Masking

First, you'll need to have RepeatModeler (http://www.repeatmasker.org/RepeatModeler/) and RepeatMasker (http://www.repeatmasker.org/) installed. We used RepeatModeler version 1.0.4 and RepeatMasker version 4.0. We also used a local installation of NCBI BLAST version 2.5.0+. The processing of the repeats identified by RepeatModeler also depends on CD-HIT (http://weizhongli-lab.org/cd-hit/). We used version 4.6.

Our procedure for identifying repetitive elements in a genome and then masking those elements is given in repeats_run.sbatch. This was the SLURM script used for running repeat masking on *Apis florea* but the same process was used for each genome individually.

For comparing identified repetitive elements to known proteins, we downloaded all manually reviewed proteins in UniProtKB (http://www.uniprot.org/downloads) on Dec. 2, 2016. We used FlyBase release 5.50 for comparison to *Drosophila melanogaster* proteins. 

The newly identified repetitive elements were combined with all arthropod-derived repetitive sequences present in Repbase (https://www.girinst.org/) downloaded on March 8, 2017.