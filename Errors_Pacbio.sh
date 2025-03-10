#!/bin/bash
#PBS -N meryl_error_filter
#PBS -l select=1:ncpus=8:mem=64gb:scratch_local=300gb
#PBS -l walltime=12:00:00
#PBS -j oe

trap 'clean_scratch' TERM EXIT

DATADIR="/storage/brno2/home/jendrb00/my_analysis/Ostrinia_nubilalis/creatig_my_datasets"
RESULTDIR="/storage/brno2/home/jendrb00/results/meryl_filtered"
SCRATCHDIR=${SCRATCHDIR:-"/scratch"}

cp -r $DATADIR/ERR11867203.fastq.gz.meryl $SCRATCHDIR/ || exit 1
cd $SCRATCHDIR || exit 2

export PATH=/storage/plzen1/home/jendrb00/meryl-1.4.1/bin:$PATH

KMER_DB="ERR11867203.fastq.gz.meryl"
FILTERED_DB="without_errors.meryl"
CUTOFF=10  #  (cut-off) z GenomeScope2.0

meryl greater-than $CUTOFF $KMER_DB output $FILTERED_DB || exit 3

OUTPUT_NAME="meryl_filtered_errors_$(date +"%Y%m%d")"
mkdir -p $RESULTDIR/$OUTPUT_NAME

cp -r $FILTERED_DB $RESULTDIR/$OUTPUT_NAME/ || exit 4

echo "Filtrace chyb dokončena. Výsledky jsou v: $RESULTDIR/$OUTPUT_NAME"
