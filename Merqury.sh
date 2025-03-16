#!/bin/bash
#PBS -q interactive@pbs-m1.metacentrum.cz 
#PBS -l walltime=5:0:0
#PBS -l select=1:ncpus=16:mem=128gb:gpu_mem=128gb:scratch_local=200gb:vnode=kirke59
#PBS -N merqury_job

# Clean scratch space on exit
trap 'clean_scratch' TERM EXIT

# Set paths and variables
Path=${PWD}  
Base="merqury_analysis" 
TimeStamp=$(date +"%Y%m%d.%H%M") 
Outdir=${Path}/${Base}_${TimeStamp} 
SCRATCH_OUTDIR="$SCRATCHDIR/merqury_results_end"

# Input files 
MERYL_FILE="${1:-reads.meryl}"
ASSEMBLY_FILE="${2:-assembly.fasta}"

# Create output directory
mkdir -p "$Outdir"
cd "$SCRATCHDIR" || exit 1

# Copy input data to scratch
cp "$Path/$MERYL_FILE" "$SCRATCHDIR/" || exit 2
cp "$Path/$ASSEMBLY_FILE" "$SCRATCHDIR/" || exit 3

# Load necessary modules and activate Merqury environment
export PATH=/storage/plzen1/home/jendrb00/meryl-1.4.1/bin:$PATH
module add mambaforge || exit 4
mamba activate /storage/plzen1/home/jendrb00/merqury-1.3/env_merqury || exit 5
module add r || exit 6
export R_LIBS_USER="/storage/plzen1/home/jendrb00/Rpackages"

# Run Merqury
merqury.sh "$MERYL_FILE" "$ASSEMBLY_FILE" merqury_output || exit 7

# Move results back to output folder
mkdir -p "$Outdir/merqury_output"
cp -r merqury_output/* "$Outdir/merqury_output/" || exit 8

# Final message
echo "Analysis completed. Results are saved in: $Outdir/merqury_output"
