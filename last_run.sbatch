#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --time=23:59:59 --qos=1day
#SBATCH --mem=50000
#SBATCH --mail-user=berubin@princeton.edu
#SBATCH --job-name=MROT_MQUA
#SBATCH --output=/Genomics/kocherlab/berubin/alignment/10bees/outputs/last_MROT_MQUA-%j.out
#SBATCH --error=/Genomics/kocherlab/berubin/alignment/10bees/errors/last_MROT_MQUA-%j.err

module add boost

/Genomics/kocherlab/berubin/local/src/last-914/src/lastdb -P 4 -cR11 alignments/MROT_db repeats/MROT/MROT_rm_masked/MROT_genome.fasta.masked

/Genomics/kocherlab/berubin/local/src/last-914/scripts/last-train -P 4 --revsym --matsym --gapsym -C 2 -E 0.01 alignments/MROT_db repeats/MQUA/MQUA_rm_masked/MQUA_genome.fasta.masked > alignments/MQUA_MROT_E1.mat

/Genomics/kocherlab/berubin/local/src/last-914/src/lastal -P 4 -m 100 -E 0.01 -C 2 -p alignments/MQUA_MROT_E1.mat alignments/MROT_db repeats/MQUA/MQUA_rm_masked/MQUA_genome.fasta.masked | /Genomics/kocherlab/berubin/local/src/last-914/src/last-split -m 1 > alignments/MQUA_MROT_E1-1.maf

/Genomics/kocherlab/berubin/local/src/last-914/scripts/maf-swap alignments/MQUA_MROT_E1-1.maf | /Genomics/kocherlab/berubin/local/src/last-914/src/last-split -m1 > alignments/MQUA.MROT.sing.maf

sed -i 's/# LAST version 914/##maf version=1/g' alignments/MQUA.MROT.sing.maf