#!/bin/bash
#PBS -N merqury_Errors_check
#PBS -l select=1:ncpus=8:mem=60gb:scratch_local=200gb
#PBS -l walltime=24:00:00
#PBS -j oe

# Clean scratch space on exit
trap 'clean_scratch' TERM EXIT

# Set paths and variables
Path=${PWD}  
Base="merqury_errors_check" 
TimeStamp=$(date +"%Y%m%d.%H%M") 
Outdir=${Path}/${Base}_${TimeStamp}  
SCRATCHDIR=${SCRATCHDIR:-"/scratch"}  

# Input file
MERYL_FILE="${1:-ERR11867203.fastq.gz.meryl}"

# Create output directory
mkdir -p "$Outdir"
cd "$SCRATCHDIR" || exit 1

# Copy input data to scratch
cp -r "$Path/$MERYL_FILE" "$SCRATCHDIR/" || exit 2

# Load required tools
module add mambaforge || exit 3
mamba activate /storage/plzen1/home/jendrb00/merqury-1.3/env_merqury || exit 4
module add r || exit 5
export R_LIBS_USER="/storage/brno2/home/jendrb00/Rpackages"
export PATH=/storage/plzen1/home/jendrb00/meryl-1.4.1/bin:$PATH

# Run Merqury error check
OUTPUT_NAME="merqury_ERRORS_PacBio_${TimeStamp}"
merqury.sh "$MERYL_FILE" "$MERYL_FILE" "$OUTPUT_NAME" || exit 6

# Move results back to output folder
cp -r "$SCRATCHDIR"/* "$Outdir" || { export CLEAN_SCRATCH=false; exit 7; }

# Final message with result path
echo "Analysis completed. Results are saved in: $Outdir"
