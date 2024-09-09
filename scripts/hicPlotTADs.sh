#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --time=4:00:00
#SBATCH --mem=64GB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ak8725@nyu.edu

module purge

singularity exec --overlay /home/ak8725/hicexplorer/hicexplorer.ext3:ro \
/scratch/work/public/singularity/cuda11.2.2-cudnn8-devel-ubuntu20.04.sif \
/bin/bash -c "source /ext3/env.sh; hicPlotTADs --tracks tracks_marks_bw.ini --region chr11:7000000-9500000 -o TAD_callers_bw_chr11.png"
