#!/bin/bash
#PBS -N jellyfish_cutoff
#PBS -l select=1:ncpus=8:mem=64gb:scratch_local=300gb
#PBS -l walltime=12:00:00
#PBS -j oe

# Clean scratch space on exit
trap 'clean_scratch' TERM EXIT

# Set paths and variables
Path=${PWD}  
Base="jellyfish_cutoff"  
TimeStamp=$(date +"%Y%m%d.%H%M")  
Outdir="${Path}/${Base}_${TimeStamp}"  
SCRATCHDIR=${SCRATCHDIR:-"/scratch"}  

# Create output directory
mkdir -p "$Outdir"

# Input 
HISTO_FILE="${1:-reads.histo}"

# Copy input file to scratch
cp "$Path/$HISTO_FILE" "$SCRATCHDIR/" || exit 1
cd "$SCRATCHDIR" || exit 2

# Identify cutoff value 
CUTOFF=$(awk '
NR == 1 {prev_count = $2; prev_freq = $1; next}
{
    if ($2 > prev_count) {
        print prev_freq;
        exit;
    }
    prev_count = $2;
    prev_freq = $1;
}' "$HISTO_FILE")

# Output the determined cutoff value
echo "Determined cutoff value: $CUTOFF"

# Check if the cutoff value was found
if [[ -z "$CUTOFF" ]]; then
    echo "Error: No cutoff value found!"
    exit 3
fi

# Save the cutoff value to a text file
echo "$CUTOFF" > "cutoff_value.txt"

# Move results back to output directory
cp "cutoff_value.txt" "$Outdir/" || exit 4

# Final message
echo "Done! The determined cutoff value has been saved in: $Outdir/cutoff_value.txt"
