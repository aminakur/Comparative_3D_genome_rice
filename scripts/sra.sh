#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=15:00:00
#SBATCH --mem=64GB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ak8725@nyu.edu

module load aspera/4.1.3.93

ascp -i /scratch/ak8725/aspera.openssh -QT -l600m -k1 -d /scratch/ak8725/sra_upload/ subasp@upload.ncbi.nlm.nih.gov:uploads/amina.kur_gmail.com_IxEhSWh1