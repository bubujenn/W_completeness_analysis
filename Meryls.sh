#!/bin/bash
#PBS -N meryl_creation
#PBS -l select=1:ncpus=8:mem=60gb:scratch_local=200gb
#PBS -l walltime=12:00:00
#PBS -j oe

trap 'clean_scratch' TERM EXIT

DATADIR="/storage/brno2/home/jendrb00/my_analysis/Abraxas_sylvata/meryls"
OUTPUT_DIR="/storage/brno2/home/jendrb00/my_analysis/Abraxas_sylvata/meryls_results"
SCRATCHDIR=${SCRATCHDIR:-"/scratch"}

mkdir -p "$SCRATCHDIR/meryls"
cp "$DATADIR/AA_Z_chr.fasta" "$DATADIR/AA_Z_W.fasta" "$DATADIR/ERR12370385.fastq.gz" "$SCRATCHDIR/meryls/" || exit 1
cd "$SCRATCHDIR/meryls" || exit 2

export PATH=/storage/plzen1/home/jendrb00/meryl-1.4.1/bin:$PATH

meryl k=19 memory=60G threads=8 count output AA_Z_chr.meryl AA_Z_chr.fasta
meryl k=19 memory=60G threads=8 count output AA_Z_W.meryl AA_Z_W.fasta
meryl k=19 memory=60G threads=8 count output ERR12370385.meryl ERR12370385.fastq.gz

mkdir -p "$OUTPUT_DIR"
cp -r "$SCRATCHDIR/meryls"/* "$OUTPUT_DIR" || { export CLEAN_SCRATCH=false; exit 3; }

echo "Meryl databáze byly úspěšně vytvořeny a zkopírovány do $OUTPUT_DIR."




