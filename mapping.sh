#!/bin/bash
#PBS -l select=1:ncpus=8:mem=32gb:scratch_local=500gb
#PBS -l walltime=12:00:00
#PBS -N Minimap2_Samtools
#PBS -j oe

# Clean scratch space on exit
trap 'clean_scratch' TERM EXIT

# Set paths and variables
Path=${PWD} 
Base="Minimap2_Samtools"  
TimeStamp=$(date +"%Y%m%d.%H%M")  
Outdir=${Path}/${Base}_${TimeStamp}  
Outfile=${Outdir}/${Base}_job.pbs  
SCRATCHDIR=${SCRATCHDIR:-"/scratch"} 

# Create output directory
mkdir -p "$Outdir"

# Input files 
REFERENCE="${1:-reference.fasta}"
QUERY="${2:-query.fasta}"

# Copy input data to scratch
cp "$Path/$REFERENCE" "$SCRATCHDIR/" || exit 1
cp "$Path/$QUERY" "$SCRATCHDIR/" || exit 2
cd "$SCRATCHDIR" || exit 3

# Load required tools
module purge (just in case)
module load minimap2
module load samtools

# Run Minimap2 alignment
minimap2 -ax asm5 --secondary=no "$REFERENCE" "$QUERY" > alignment.sam || exit 4

# Process alignment file
samtools view -F 2048 -h alignment.sam > filtered.sam || exit 5
samtools sort -o filtered.bam filtered.sam || exit 6
samtools index filtered.bam || exit 7

# Extract scaffold-to-chromosome mapping
samtools view filtered.bam | awk '{print $3, $1}' | sort | uniq > scaffold_to_chr.tsv || exit 8

# Compute query coverage using awk
samtools view filtered.bam | awk '
    { 
        qname[$1] = $1;  # Query (scaffold) name
        qlen[$1] += $9;  # Aligned length
    }
    END {
        for (scaffold in qname)
            print scaffold, qlen[scaffold];  
    }' > query_coverage.tsv || exit 9

# Move results back to output folder
cp -r "$SCRATCHDIR"/* "$Outdir/" || { export CLEAN_SCRATCH=false; exit 10; }

# Final message 
echo "Done! Results saved in: $Outdir"
