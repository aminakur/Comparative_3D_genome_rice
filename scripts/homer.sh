#!/bin/bash
#
##SBATCH --nodes=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=6:00:00
#SBATCH --mem=64GB
#SBATCH --mail-type=END
#SBATCH --mail-user=ak8725@nyu.edu

module load homer/4.11
#default parameters
findMotifsGenome.pl /scratch/ak8725/az_mrg/hicFindTADs/hicFindTADs2_out/az_1kb_boundaries.bed \
/scratch/ak8725/az_mrg/azucena.fna hicexp_1kb -size given -len 10 -p 8 \
-bg /scratch/ak8725/az_mrg/HOMER/homer_hicexp_background_seqs.bed

findMotifsGenome.pl /scratch/ak8725/az_mrg/hicFindTADs/hicFindTADs1_out/az_2kb_boundaries.bed \
/scratch/ak8725/az_mrg/azucena.fna hicexp_2kb -size given -len 10 -p 8 \
-bg /scratch/ak8725/az_mrg/HOMER/homer_hicexp_background_seqs.bed

findMotifsGenome.pl /scratch/ak8725/az_mrg/hicFindTADs/hicFindTADs13_out/az_5kb_boundaries.bed \
/scratch/ak8725/az_mrg/azucena.fna hicexp_5kb -size given -len 10 -p 8 \
-bg /scratch/ak8725/az_mrg/HOMER/homer_hicexp_background_seqs.bed

#parameters from cotton paper
# findMotifsGenome.pl /scratch/ak8725/az_mrg/hicFindTADs/hicFindTADs2_out/az_1kb_boundaries.bed \
# /scratch/ak8725/az_mrg/azucena.fna hicexp_1kb_cotton -size given -size 200 -p 8 \
# -bg /scratch/ak8725/az_mrg/HOMER/homer_hicexp_background_seqs.bed

# findMotifsGenome.pl /scratch/ak8725/az_mrg/hicFindTADs/hicFindTADs1_out/az_2kb_boundaries.bed \
# /scratch/ak8725/az_mrg/azucena.fna hicexp_2kb_cotton -size given -size 200 -p 8 \
# -bg /scratch/ak8725/az_mrg/HOMER/homer_hicexp_background_seqs.bed

# findMotifsGenome.pl /scratch/ak8725/az_mrg/hicFindTADs/hicFindTADs13_out/az_5kb_boundaries.bed \
# /scratch/ak8725/az_mrg/azucena.fna hicexp_5kb_cotton -size given -size 200 -p 8 \
# -bg /scratch/ak8725/az_mrg/HOMER/homer_hicexp_background_seqs.bed

