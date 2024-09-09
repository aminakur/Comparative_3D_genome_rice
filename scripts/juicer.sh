#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=10:00:00
#SBATCH --mem=256GB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ak8725@nyu.edu

module load jdk/11.0.9
java -Xmx48000m -Djava.awt.headless=true -jar ../juicer_tools_1.22.01.jar \
pre -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,2000,1000 \
--threads 16 az-2_mapped.pairs az-2.hic azucena.genome






