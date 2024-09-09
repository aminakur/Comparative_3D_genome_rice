#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH --mem=256GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ak8725@nyu.edu

module purge
module load perl/intel/5.32.0
module load mummer/intel/4.0.0rc1

# Directory list
directories=("NPB" "azucena" "IR64" "oruf" "omer")

# Iterate over the directories
for i in "${!directories[@]}"; do
    for ((j=i+1; j<${#directories[@]}; j++)); do
        reference_dir="./${directories[$i]}"
        query_dir="./${directories[$j]}"

        output_dir="./${directories[$i]}-${directories[$j]}_alignments"

        # Check if output directory exists
        if [ -d "$output_dir" ]; then
            echo "Output directory $output_dir already exists. Skipping to next pair."
            continue
        else
            # Create output directory if it doesn't exist
            mkdir -p "$output_dir"
        fi

        # Check if reference and query directories exist
        if [ ! -d "$reference_dir" ]; then
            echo "Reference directory $reference_dir does not exist."
            continue
        fi
        if [ ! -d "$query_dir" ]; then
            echo "Query directory $query_dir does not exist."
            continue
        fi

        # Iterate through chromosomes and perform alignments
        for chr_index in {01..12}; do
            chr_name="chr${chr_index}"

            ref_file="${reference_dir}/${chr_name}.fna"
            query_file="${query_dir}/${chr_name}.fna"

            # Check if reference and query files exist
            if [ ! -f "$ref_file" ]; then
                echo "Reference file $ref_file does not exist in $reference_dir"
                continue
            fi
            if [ ! -f "$query_file" ]; then
                echo "Query file $query_file does not exist in $query_dir"
                continue
            fi
            
            # Run nucmer alignment
            nucmer --mum -l 100 -c 1000 -d 50 -p "${output_dir}/${chr_name}" \
                "$ref_file" "$query_file"
            
            # Check if nucmer command was successful
            if [ $? -ne 0 ]; then
                echo "nucmer alignment failed for $ref_file and $query_file"
            else
                echo "nucmer alignment succeeded for $ref_file and $query_file"
            fi
        done
    done
done