#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=128GB
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=ak8725@nyu.edu

#parameters from Cantaloupe melon reveals 3D, 2022

module purge

/scratch/work/public/singularity/run-hicexplorer-3.7.2.bash hicFindTADs -m /scratch/ak8725/az_new/az_new.mcool::resolutions/5000 \
--outPrefix az_new_5kb_hicex -p 8 --correctForMultipleTesting fdr \
--minDepth 15000 --maxDepth 50000 --step 5000 --thresholdComparisons 0.05 --delta 0.15

/scratch/work/public/singularity/run-hicexplorer-3.7.2.bash hicFindTADs -m /scratch/ak8725/az_new/az_new.mcool::resolutions/2000 \
--outPrefix az_new_2kb_hicex -p 8 --correctForMultipleTesting fdr \
--minDepth 15000 --maxDepth 150000 --step 5000 --thresholdComparisons 0.05 --delta 0.01

/scratch/work/public/singularity/run-hicexplorer-3.7.2.bash hicFindTADs -m /scratch/ak8725/az_new/az_new.mcool::resolutions/1000 \
--outPrefix NPB_new_1kb_hicex -p 8 --correctForMultipleTesting fdr \
--minDepth 3000 --maxDepth 30000 --step 1000 --thresholdComparisons 0.05 --delta 0.01


