#!/bin/bash
#
##SBATCH --nodes=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=4:00:00
#SBATCH --mem=10GB
#SBATCH --mail-type=END
#SBATCH --mail-user=ak8725@nyu.edu

module load samtools/intel/1.14
samtools faidx azucena.fna
cut -f 1,2 azucena.fna.fai > azucena.genome

