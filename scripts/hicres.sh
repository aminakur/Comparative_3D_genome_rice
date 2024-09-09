#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --time=48:00:00
#SBATCH --mem=64GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ak8725@nyu.edu

singularity run --bind /scratch/ak8725/az_mrg:/home/input \
--bind /scratch/ak8725/az_mrg:/home/output ~/hicres/hicres.sif -m bam -t 40 \
-c azucena.chrom.sizes -b az_mrg_mapped.PT.bam 
