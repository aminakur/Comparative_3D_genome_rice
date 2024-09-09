#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --time=12:00:00
#SBATCH --mem=256GB
#SBATCH --mail-user=ak8725@nyu.edu
#SBATCH --mail-type=END,FAIL

module load repeatmasker/4.1.2

#/scratch/work/public/singularity/run-RepeatMasker-4.1.2-p1.bash RepeatMasker -xsmall -pa 12 -lib TIGR_Oryza_Repeats_v_3_3.fna azucena.m.fna
/scratch/work/public/singularity/run-RepeatMasker-4.1.2-p1.bash RepeatMasker -xsmall -pa 12 -lib TIGR_Oryza_Repeats_v_3_3.fna IR64.m.fna
#/scratch/work/public/singularity/run-RepeatMasker-4.1.2-p1.bash RepeatMasker -xsmall -pa 12 -lib TIGR_Oryza_Repeats_v_3_3.fna NPB.m.fna
# /scratch/work/public/singularity/run-RepeatMasker-4.1.2-p1.bash RepeatMasker -xsmall -pa 12 -lib TIGR_Oryza_Repeats_v_3_3.fna omer.m.fna
# /scratch/work/public/singularity/run-RepeatMasker-4.1.2-p1.bash RepeatMasker -xsmall -pa 12 -lib TIGR_Oryza_Repeats_v_3_3.fna orufi.m.fna
