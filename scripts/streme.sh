#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=128GB
#SBATCH --mail-type=END
#SBATCH --mail-user=ak8725@nyu.edu

module load meme/openmpi/intel/5.3.0

# streme --oc hicexp_1kb --n streme_background_20k.fa --p /scratch/ak8725/az_mrg/hicexp_1kb_boundaries.fna
# streme --oc hicexp_2kb --n streme_background_20k.fa  --p /scratch/ak8725/az_mrg/hicexp_2kb_boundaries.fna
streme --oc hicexp_5kb --n streme_background_12k.fa --p /scratch/ak8725/az_mrg/hicexp_5kb_boundaries.fna
