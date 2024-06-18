#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=64GB

module load bedtools/intel/2.29.2

# bedtools intersect -a NPB_500bp_windows.bed -b npb_del.bed -wao > NPB_del_in_windows.bed
# bedtools intersect -a NPB_500bp_windows.bed -b npb_dup.bed -wao > NPB_dup_in_windows.bed
# bedtools intersect -a NPB_500bp_windows.bed -b npb_ins.bed -wao > NPB_ins_in_windows.bed
# bedtools intersect -a NPB_500bp_windows.bed -b npb_inv.bed -wao > NPB_inv_in_windows.bed

bedtools intersect -a NPB_500bp_windows.bed -b NB_final_snp.bed -wao > NPB_snp_in_windows.bed