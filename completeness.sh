#!/bin/bash
#PBS -N calculate_W_completeness
#PBS -l select=1:ncpus=2:mem=4gb:scratch_local=10gb
#PBS -l walltime=01:00:00
#PBS -j oe

# Clean scratch space on exit
trap 'clean_scratch' TERM EXIT

# Set paths and variables
Path=${PWD} 
Base="W_completeness_calc" 
TimeStamp=$(date +"%Y%m%d.%H%M")  
Outdir=${Path}/${Base}_${TimeStamp}  
SCRATCHDIR=${SCRATCHDIR:-"/scratch"}  

# Input files 
MISSING_FILE="${1:-W_missing_in_assembly_sum.txt}"
EXPECTED_FILE="${2:-expected_W_size_sum.txt}"

# Create output directory
mkdir -p "$Outdir"
cd "$SCRATCHDIR" || exit 1

# Copy input data to scratch
cp "$Path/$MISSING_FILE" "$SCRATCHDIR/" || exit 2
cp "$Path/$EXPECTED_FILE" "$SCRATCHDIR/" || exit 3

# Read values from files
expected=$(cat "$EXPECTED_FILE")
missing=$(cat "$MISSING_FILE")

# Calculate W chromosome completeness percentage
completeness=$(awk -v missing="$missing" -v expected="$expected" 'BEGIN { printf "%.5f\n", (1 - (missing / expected)) * 100 }') || exit 4
echo "$completeness" > W_FINAL_completeness.txt

# Print and store results
cat W_FINAL_completeness.txt || exit 5
cp W_FINAL_completeness.txt "$Outdir/" || { export CLEAN_SCRATCH=false; exit 6; }

# Final message 
echo "W chromosome completeness results saved in: $Outdir/W_FINAL_completeness.txt"
