#!/bin/bash
#PBS -q interactive@pbs-m1.metacentrum.cz 
#PBS -l walltime=1:0:0
#PBS -l select=1:ncpus=16:mem=128gb:gpu_mem=128gb:scratch_local=200gb:vnode=kirke59
#PBS -N merqury_check_eroors_job

# Clean scratch space on exit
trap 'clean_scratch' TERM EXIT

# Set paths and variables
Path=${PWD}  
Base="merqury_errors_check" 
TimeStamp=$(date +"%Y%m%d.%H%M") 
Outdir=${Path}/${Base}_${TimeStamp}  
SCRATCHDIR=${SCRATCHDIR:-"/scratch"}  

# Input file
MERYL_FILE="${1:-fastq.gz.meryl}"

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

# Final message
echo "Analysis completed. Results are saved in: $Outdir"
