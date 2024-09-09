#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=1:00:00
#SBATCH --mem=256GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ak8725@nyu.edu

module purge
module load perl/intel/5.32.0
module load mummer/intel/4.0.0rc1
module load gnuplot/gcc/5.4.1
module load imagemagick/intel/7.0.10

reference_dir="./oruf"
query_dir="./omer"

output_dir="./oruf-omer_alignments"
png_output="O.rufipogon-O.meridionalis.png"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Iterate through chromosomes and perform alignments
for chr_index in {01..12}; do
    chr_name="chr${chr_index}"
    
    nucmer --mum -l 100 -c 1000 -d 50 -p "${output_dir}/${chr_name}" \
        "${reference_dir}/${chr_name}.fna" "${query_dir}/${chr_name}.fna"
    
    mummerplot -postscript -p "${output_dir}/${chr_name}" \
        "${output_dir}/${chr_name}.delta"
    
    convert "${output_dir}/${chr_name}.ps" -rotate 90 \
        "${output_dir}/${chr_name}.png"
done

# Combine PNG images into a 4x3 grid
convert \( "${output_dir}/chr01.png" "${output_dir}/chr02.png" "${output_dir}/chr03.png" "${output_dir}/chr04.png" +append \) \
        \( "${output_dir}/chr05.png" "${output_dir}/chr06.png" "${output_dir}/chr07.png" "${output_dir}/chr08.png" +append \) \
        \( "${output_dir}/chr09.png" "${output_dir}/chr10.png" "${output_dir}/chr11.png" "${output_dir}/chr12.png" +append \) \
        -append -gravity Center "$png_output"