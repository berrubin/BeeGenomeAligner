#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=10
#SBATCH --time=6-23:00 --qos=1wk
#SBATCH --mem=20000
#SBATCH --mail-user=berubin@princeton.edu
#SBATCH --job-name=AFLOrep
#SBATCH --output=repeats_AFLO-%j.out
#SBATCH --error=repeats_AFLO-%j.err

module add RepeatModeler

export PATH=$PATH:/Genomics/kocherlab/berubin/local/src/ncbi-blast-2.5.0+/bin/

BuildDatabase -name AFLO_db -engine ncbi AFLO_genome.fasta

RepeatModeler -database AFLO_db -pa 3 -engine ncbi >& AFLO.out

cd-hit-est -i consensi.fa -o rm_noredun.fa -c 0.8 -M 20000 -n 5 -aS 0.8 -r 1 -T 10

blastx -query rm_noredun.fa -db uniprot_sprot_blastdb -outfmt '6 qseqid sseqid pident length evalue bitscore qlen slen' -out rm_noredun_against_uniprot.txt -num_threads 10 

blastx -query rm_noredun.fa -db dmel_blastdb -outfmt '6 qseqid sseqid pident length evalue bitscore qlen slen' -out rm_noredun_against_dmel.txt -num_threads 10

python filter_repeats.py rm_noredun.fa rm_noredun_noprot.fa rm_noredun_against_uniprot.txt rm_noredun_against_dmel.txt

 module add RepeatMasker

cat rm_noredun_noprot.fa arthropoda_repeatmasker_lib.fa > AFLO_arthropoda_consensi.fa

RepeatMasker -e ncbi -pa 10 -lib AFLO_arthropoda_consensi.fa -gff -xsmall -source -dir AFLO_rm_masked AFLO_genome.fasta

