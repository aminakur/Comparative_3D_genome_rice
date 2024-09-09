#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=96:00:00
#SBATCH --mem=256GB
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=ak8725@nyu.edu

module purge

module load samtools/intel/1.12
bgzip az_new_mapped.pairs
module load pypairix/0.3.7
pairix az_new_mapped.pairs.gz

../coolerfile "cload pairix -p 16 azucena.genome:200 /scratch/ak8725/az_new/az_new_mapped.pairs.gz /scratch/ak8725/az_new/az_new.cool"

../coolerfile "zoomify -r 200,400,600,800,1000,2000,5000,10000,20000,40000,80000,160000,320000,640000,1280000 \
--balance -p 16 az_new.cool"

