#!/bin/bash
#PBS -N meryl_error_filter
#PBS -l select=1:ncpus=8:mem=64gb:scratch_local=300gb
#PBS -l walltime=12:00:00
#PBS -j oe

# Clean scratch space on exit
trap 'clean_scratch' TERM EXIT

# Set paths and variables
Path=${PWD}  
Base="meryl_error_filter"  
TimeStamp=$(date +"%Y%m%d.%H%M")  
Outdir=${Path}/${Base}_${TimeStamp} 
SCRATCHDIR=${SCRATCHDIR:-"/scratch"}  

# Create output directory
mkdir -p "$Outdir"

# Input file 
KMER_DB="${1:-ERR11867203.fastq.gz.meryl}"
FILTERED_DB="without_errors.meryl"
CUTOFF=10  # Cut-off value from GenomeScope2.0

# Copy input data to scratch
cp -r "$Path/$KMER_DB" "$SCRATCHDIR/" || exit 1
cd "$SCRATCHDIR" || exit 2

# Load Meryl
export PATH=/storage/plzen1/home/jendrb00/meryl-1.4.1/bin:$PATH

# Run Meryl filtering
meryl greater-than $CUTOFF "$KMER_DB" output "$FILTERED_DB" || exit 3

# Move results back to output folder
cp -r "$FILTERED_DB" "$Outdir/" || exit 4

# Final message with result path
echo "Done! Filtered Meryl database saved in: $Outdir"
