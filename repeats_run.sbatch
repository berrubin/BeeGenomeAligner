#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=10
#SBATCH --time=6-23:00 --qos=1wk
#SBATCH --mem=20000
#SBATCH --mail-user=berubin@princeton.edu
#SBATCH --job-name=EMEXrep
#SBATCH --output=/Genomics/kocherlab/berubin/alignment/10bees/repeats/EMEX/repeats_EMEX-%j.out
#SBATCH --error=/Genomics/kocherlab/berubin/alignment/10bees/repeats/EMEX/repeats_EMEX-%j.err

module add RepeatModeler

export PATH=$PATH:/Genomics/kocherlab/berubin/local/src/ncbi-blast-2.5.0+/bin/

cd EMEX

BuildDatabase -name EMEX_db -engine ncbi /Genomics/kocherlab/berubin/alignment/10bees/genomes/EMEX_genome.fasta

RepeatModeler -database EMEX_db -pa 3 -engine ncbi >& EMEX.out

/Genomics/kocherlab/berubin/local/src/cdhit/cd-hit-est -i ./RM_*/consensi.fa -o rm_noredun.fa -c 0.8 -M 20000 -n 5 -aS 0.8 -r 1 -T 10

blastx -query rm_noredun.fa -db /Genomics/kocherlab/berubin/annotation/reference_genomes/uniprot/uniprot_sprot_blastdb -outfmt '6 qseqid sseqid pident length evalue bitscore qlen slen' -out rm_noredun_against_uniprot.txt -num_threads 10 

blastx -query rm_noredun.fa -db /Genomics/kocherlab/berubin/annotation/reference_genomes/drosophila/dmel_blastdb -outfmt '6 qseqid sseqid pident length evalue bitscore qlen slen' -out rm_noredun_against_dmel.txt -num_threads 10

python /Genomics/kocherlab/berubin/annotation/filter_repeats.py rm_noredun.fa rm_noredun_noprot.fa rm_noredun_against_uniprot.txt rm_noredun_against_dmel.txt

 module add RepeatMasker

cat rm_noredun_noprot.fa /Genomics/kocherlab/berubin/annotation/arthropoda_repeatmasker_lib.fa > EMEX_arthropoda_consensi.fa

mkdir EMEX_rm_masked

RepeatMasker -e ncbi -pa 10 -lib /Genomics/kocherlab/berubin/alignment/10bees/repeats/AFLO/AFLO_arthropoda_consensi.fa -gff -xsmall -source -dir EMEX_rm_masked /Genomics/kocherlab/berubin/alignment/10bees/genomes/EMEX_genome.fasta

