#!/bin/bash
#PBS -N meryl_creation
#PBS -l select=1:ncpus=8:mem=60gb:scratch_local=200gb
#PBS -l walltime=12:00:00
#PBS -j oe

# Clean scratch space on exit
trap 'clean_scratch' TERM EXIT

# Set paths and variables
Path=${PWD}  
Base="meryl_creation"  
TimeStamp=$(date +"%Y%m%d.%H%M")  
Outdir=${Path}/${Base}_${TimeStamp}  
SCRATCHDIR=${SCRATCHDIR:-"/scratch"} 

# Create output directory
mkdir -p "$Outdir"

# Input files
CHR_FASTA="${1:-AA_Z_chr.fasta}"
W_FASTA="${2:-AA_Z_W.fasta}"
FASTQ="${3:-*.fastq.gz}"  
KMER_SIZE="${4:-19}"  # (Default k-mer size is 19 unless specified)

# Copy input data to scratch
mkdir -p "$SCRATCHDIR/meryls"
cp "$Path/$CHR_FASTA" "$Path/$W_FASTA" "$Path/$FASTQ" "$SCRATCHDIR/meryls/" || exit 1
cd "$SCRATCHDIR/meryls" || exit 2

# Load Meryl
export PATH=/storage/plzen1/home/jendrb00/meryl-1.4.1/bin:$PATH

# Create Meryl databases with specified k-mer size
meryl k="$KMER_SIZE" memory=60G threads=8 count output AA_Z_chr.meryl "$CHR_FASTA"
meryl k="$KMER_SIZE" memory=60G threads=8 count output AA_Z_W.meryl "$W_FASTA"
meryl k="$KMER_SIZE" memory=60G threads=8 count output ERR_fastq.meryl "$FASTQ"

# Move results back to output folder
cp -r "$SCRATCHDIR/meryls"/* "$Outdir" || { export CLEAN_SCRATCH=false; exit 3; }

# Final message with result path
echo "Done! Meryl databases saved in: $Outdir (k-mer size: $KMER_SIZE)"





