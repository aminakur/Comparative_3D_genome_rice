#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=4:00:00
#SBATCH --mem=20GB
#SBATCH --mail-user=ak8725@nyu.edu
#SBATCH --mail-type=END,FAIL

module load bwa/intel/0.7.17
bwa index /scratch/ak8725/genomes/azucena.fna
