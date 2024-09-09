#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=168:00:00
#SBATCH --mem=10GB
#SBATCH --mail-type=END
#SBATCH --mail-user=ak8725@nyu.edu

module load python/intel/3.8.6
python -m pip install tabulate
python3 get_qc.py -p az_new_stats.txt

