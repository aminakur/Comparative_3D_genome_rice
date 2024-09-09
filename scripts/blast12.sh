#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00
#SBATCH --mem=256GB
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=ak8725@nyu.edu

module purge
module load  blast+/2.13.0

blastn -query npb_tads.fna -subject ../genomes/IR64.fna -out npb_ir64_all_blast_results.txt -max_hsps 1 -max_target_seqs 1 \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"


