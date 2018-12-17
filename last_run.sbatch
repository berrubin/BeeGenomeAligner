#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --time=23:59:59 --qos=1day
#SBATCH --mem=50000
#SBATCH --mail-user=berubin@princeton.edu
#SBATCH --job-name=AFLO_AMEL
#SBATCH --output=last_AFLO_AMEL-%j.out
#SBATCH --error=last_AFLO_AMEL-%j.err

lastdb -P 4 -cR11 AFLO_db AFLO_genome.fasta.masked

last-train -P 4 --revsym --matsym --gapsym -C 2 -E 0.01 AMEL_genome.fasta.masked > AMEL_AFLO_E1.mat

lastal -P 4 -m 100 -E 0.01 -C 2 -p AMEL_AFLO_E1.mat AFLO_db AMEL_genome.fasta.masked | last-split -m 1 > AMEL_AFLO_E1-1.maf

maf-swap AMEL_AFLO_E1-1.maf | last-split -m1 > AMEL.AFLO.sing.maf

sed -i 's/# LAST version 914/##maf version=1/g' AMEL.AFLO.sing.maf