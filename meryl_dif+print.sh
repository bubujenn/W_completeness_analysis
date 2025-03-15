#!/bin/bash
#PBS -N meryl_difference_print
#PBS -l select=1:ncpus=8:mem=64gb:scratch_local=300gb
#PBS -l walltime=12:00:00
#PBS -j oe

# Clean scratch space on exit
trap 'clean_scratch' TERM EXIT

# Set paths and variables
Path=${PWD}  
Base="meryl_difference_analysis"  
TimeStamp=$(date +"%Y%m%d.%H%M")  
Outdir=${Path}/${Base}_${TimeStamp}  
SCRATCHDIR=${SCRATCHDIR:-"/scratch"}  

# Input files 
PACBIO_MERYL="${1:-without_errors.meryl}"
AA_Z_MERYL="${2:-AA_Z_chr.fasta.meryl}"
AA_Z_W_MERYL="${3:-AA_Z_Wasmbl.fasta.meryl}"

# Create output directory
mkdir -p "$Outdir"
cd "$SCRATCHDIR" || exit 1

# Copy input data to scratch
cp -r "$Path/$PACBIO_MERYL" "$SCRATCHDIR/" || exit 2
cp -r "$Path/$AA_Z_MERYL" "$SCRATCHDIR/" || exit 3
cp -r "$Path/$AA_Z_W_MERYL" "$SCRATCHDIR/" || exit 4

# Load Meryl
export PATH=/storage/plzen1/home/jendrb00/meryl-1.4.1/bin:$PATH

# Create directory for results
mkdir -p meryl_results

# Compare PacBio data with AA_Zchr
meryl difference "$PACBIO_MERYL" "$AA_Z_MERYL" output meryl_results/expected_W_size.meryl || exit 5

# Compare PacBio data with AA_Z_W
meryl difference "$PACBIO_MERYL" "$AA_Z_W_MERYL" output meryl_results/W_missing_in_assembly.meryl || exit 6

# Sum unique k-mers for expected W chromosome size
meryl print meryl_results/expected_W_size.meryl | awk '{sum += $2} END {print sum}' > meryl_results/expected_W_size_sum.txt || exit 7

# Sum missing unique k-mers of W in the assembled genome
meryl print meryl_results/W_missing_in_assembly.meryl | awk '{sum += $2} END {print sum}' > meryl_results/W_missing_in_assembly_sum.txt || exit 8

# Move results back to output folder
cp -r meryl_results/* "$Outdir" || { export CLEAN_SCRATCH=false; exit 9; }

# Final message 
echo "First part completed. Results saved in: $Outdir"
