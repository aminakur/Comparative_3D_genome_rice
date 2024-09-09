#code failed, see hitad5kb.log for error message
#used az_mrg TADs instead of az_new

#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=128GB
#SBATCH --mail-user=ak8725@nyu.edu
#SBATCH --mail-type=BEGIN,END,FAIL,ARRAY_TASKS

bash /scratch/work/public/singularity/run-hitad-0.4.5-r1.bash \
hitad --output az_new_hitad5kb.txt -d meta_file5kb --logFile hitad5kb.log --maxsize 3000000 -p 8

bash /scratch/work/public/singularity/run-hitad-0.4.5-r1.bash \
hitad --output az_new_hitad2kb.txt -d meta_file2kb --logFile hitad2kb.log --maxsize 3000000 -p 8