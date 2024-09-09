## Scripts for data analysis

Ipython/Jupyter notebooks

R, python, bash, batch scripts  

- TAD calling
    - arrowhead.sh - Arrowhead from juicertools
    - hicFindTADs.sh - HiCExplorer utility
    - hitad.sbatch - HiTAD
    - meta_file10kb - example meta file for HiTAD
      
- TAD plotting with HiCExplorer
    - hicPlotTADs.sh
    - tracks_marks_bw.ini

- Metagene profiles with deeptools
    - deeptools.sh

- Map resolution estimation with HiCRes
    - hicres.sh

- Motif analysis
    - HOMER.ipynb
    - homer.sh
    - STREME.ipynb
    - streme.sh
    
- Comparison of strong and weak TAD boundaries
    - Strong_vs_weak_boundaries_comparison.ipynb

- Analyzing TADs and boundaries
    - Overlapping_TADs_boundaries.ipynb
    - TADs vs boundaries vs non-TADs analysis.ipynb

- Micro-C reads processing, contact map generation and QC
    - bwa_index.sh
    - genome.sh
    - micro-C.bash
    - micro-C.sbatch - this script invokes micro-C.bash to generate .bam and .pairs files
    - qc.sh - invokes get_qc.py (Dovetail Genomics' script)
    - get_qc.py
    - cooler.sh
    - juicer.sh
    - coolerfile - to invoke cooler within cooler.sh

- Analyzing TAD conservation with BLAST
    - 1-2.Analyzing_TAD_conservation_with_BLAST.ipynb - Startegies 1 and 2
    - blast12.sh - Startegies 1 and 2
    - 3.Analyzing_TAD_conservation_with_BLAST.ipynb - Strategy 3 (used in the paper)
    - blast3.sh - Strategy 3 
    - 4.NPB-Az.ipynb, 4.NPB-IR64.ipynb, 4.NPB-omer.ipynb, 4.NPB-oruf.ipynb - Startegy 4
    - blast4.sh - Strategy 4
    - Conserved_TADs.ipynb - analysis of features in TADs conservation groups

- Repeatmasking genomes
    - repeatmasker.sh

- dNdS calculation with orthologR
    - dNdS.R
    - dNdS.sh
    - orthologr_dNdS_genomes.R

- Uploading fastq files to SRA
    - sra.sh

- Comparative analysis of TADs/boundaries identified with liftover tool
    - Comparative analysis boundaries 5kb.ipynb
    - Comparative_analysis_TADs_5kb.ipynb

- Assessing the reproducibility of Micro-C data with HiCRep
    - HiCRep.R

- Lifting over coordinates from one genome to another
    - liftOver.ipynb
    - liftover.sh

- Creating and visualizing dotplots of individual chromosome alignments
    - Dotplot_with_mummer.ipynb
    - dotplot.sh
    - dotplot-all.sh

- Creating dotplots for regions with conserved TADs
    - Dotplots_of_regions_with_conserved_tads.ipynb

- Visualizing TADs, genetic and epigenetic features in Coolbox
    - coolbox_conseved_tads.ipynb
    - TADs+bw+bed coolbox.ipynb

- Distance decay curves
    - Distance decay curves.ipynb

- Other scripts
    - Extracting features from GFF for plotting.ipynb
    - Gene orientation at boundaries.ipynb
    - Genes_per_TAD.ipynb
    - How many promoters or genes in boundaries.ipynb
    - Nextflow_example.ipynb

- Comparison of Micro-C maps with CHESS
    - Whole-genome-lifted-analysis-all.ipynb
    - NPB-Omer.ipynb
    - CHESS_example_analysis.ipynb
    - NPB-Omer_tails_analysis.ipynb
    - NPB-Oruf_tails_analysis.ipynb