#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00
#SBATCH --mem=32GB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ak8725@nyu.edu

module purge
module load deeptools/3.5.0

computeMatrix scale-regions -S /scratch/ak8725/az_mrg/fitcons/PRO-seq_true.bw \
                -R /scratch/ak8725/az_mrg/arrowhead/az_mrg_2kb_arrowhead/2000_blocks.bedpe \
/scratch/ak8725/az_mrg/arrowhead/az_mrg_5kb_arrowhead/5000_blocks.bedpe \
/scratch/ak8725/az_mrg/hicFindTADs/hicFindTADs2_out/az_1kb_domains.bed \
/scratch/ak8725/az_mrg/hicFindTADs/hicFindTADs1_out/az_2kb_domains.bed \
/scratch/ak8725/az_mrg/hicFindTADs/hicFindTADs13_out/az_5kb_domains.bed \
/scratch/ak8725/az_mrg/hitad/az_hitad1kb_0lvl.txt \
/scratch/ak8725/az_mrg/hitad/az_hitad2kb_0lvl.txt \
/scratch/ak8725/az_mrg/hitad/az_hitad5kb_0lvl.txt \
                              --beforeRegionStartLength 20000 \
                              --regionBodyLength 20000 \
                              --afterRegionStartLength 20000 \
                              --binSize 200 \
                              --skipZeros -p 4 -o matrix5.mat.gz \


plotProfile -m matrix5.mat.gz \
            --colors palevioletred palevioletred palevioletred palevioletred palevioletred palevioletred palevioletred palevioletred \
            -out TADs_all_pro.pdf \
            --numPlotsPerRow 4 --startLabel 0 --endLabel 0 \
            --perGroup \
            --plotTitle 'PRO-seq' \
            --regionsLabel 'Arrowhead 2kb' 'Arrowhead 5kb' 'HiCExplorer 1kb' 'HiCExplorer 2kb' 'HiCExplorer 5kb' 'HiTAD 1kb' 'HiTAD 2kb' 'HiTAD 5kb' \
            --samplesLabel 'PRO-seq'
            # --yMin 0.6 0.5 0.4 1 7 0.6 \
            # --yMax 0.75 0.65 0.55 1.55 19 0.75 \
            #--plotTitle "80% overlapping TADs"