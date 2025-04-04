#!/bin/bash
#PBS -q interactive@pbs-m1.metacentrum.cz (vnodes with samtools)
#PBS -l walltime=5:0:0
#PBS -l select=1:ncpus=16:mem=128gb:scratch_local=200gb
#PBS -N mapping_H2_H1

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

# Load required tools + module purge (just in case.. # Ensuring compatibility between minimap2 and samtools by loading compatible module versions and resolving Python conflicts.)
module purge
module load python/3.9.12-gcc-10.2.1
module load samtools/1.13-gcc-10.2.1
module load minimap2/2.22-gcc-10.2.1

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
