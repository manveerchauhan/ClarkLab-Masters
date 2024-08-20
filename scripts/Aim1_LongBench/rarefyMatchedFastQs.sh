#!/bin/bash
# Script to create randomly subsampled fastq files from an original for a list of different sampling rates
# Usage- change sample_rates array to whatever you want to subsample to then run the following command in terminal:
# ./rarefyMatchedFastQs.sh <fast_q file> <output_prefix> <sampling_rate_1> <sampling_rate_2>...etc
# **Make sure bbmap is installed in the current conda environment**

# Define array of sampling rates
sample_rates=(0.75 0.50 0.25 0.10)
# Define the seed for reproducibility
seed=42

# Define global variables from user-entered parameters
fastq_file="$1"
output_prefix="$2"
readNumberTrackerFile="${output_prefix}_readRarePoints.txt"

# Shift the first two arguments so that "$@" contains only the sampling rates
#shift 2
# Define array of sampling rates from user input
#sample_rates=("$@")

# Function to calculate the total number of reads in any given FASTQ file
calculate_total_reads() {
    local fastq_file="$1"
    local total_lines
    local total_reads
    
    # Count the number of lines in the FASTQ file and divide by 4 to get the number of reads
    total_lines=$(wc -l < "$fastq_file")
    total_reads=$((total_lines / 4))
    
    # Return the total reads as the function's output
    echo "$total_reads"
}

total_reads=$(calculate_total_reads "$fastq_file")
echo "Original Total Reads: $total_reads reads" >> "$readNumberTrackerFile"

# Loop over the array of sample rates and perform subsampling, keeping track of total number of reads after
for rate in "${sample_rates[@]}"; do
    # Create an output file name based on the sampling rate
    output_file="${output_prefix}_${rate//./}percent.fastq"

    # Subsample the FASTQ file using reformat.sh with the current sample rate
    reformat.sh in="$fastq_file" out="$output_file" samplerate="$rate" sampleseed="$seed"

    # Calculate the number of reads in the subsampled FASTQ file
    reads_sample_rate=$(calculate_total_reads "$output_file")

    # Append the result to the read number tracker file
    echo "Sampling rate $rate: $reads_sample_rate reads" >> "$readNumberTrackerFile"
done