#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=168:00:00
#SBATCH --mem=64GB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ak8725@nyu.edu

module load jdk/11.0.9
java -Xmx48000m -Djava.awt.headless=true -jar /scratch/ak8725/az_mrg/arrowhead/juicer_tools_1.22.01.jar \
arrowhead -m 400 -r 2000 -k KR --threads 16 /scratch/ak8725/az_new/az_new.hic az_mrg_2kb_arrowhead






