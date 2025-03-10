#!/bin/bash
#PBS -N meryl_difference_print
#PBS -l select=1:ncpus=8:mem=64gb:scratch_local=300gb
#PBS -l walltime=12:00:00
#PBS -j oe

trap 'clean_scratch' TERM EXIT

DATADIR="/storage/brno2/home/jendrb00/my_analysis/Ostrinia_nubilalis/creatig_my_datasets"
OUTPUT_DIR="/storage/brno2/home/jendrb00/results/meryl_difference_print_results"
export PATH=/storage/plzen1/home/jendrb00/meryl-1.4.1/bin:$PATH

PACBIO_MERYL="without_errors.meryl"
AA_Z_MERYL="AA_Z_chr.fasta.meryl"
AA_Z_W_MERYL="AA_Z_Wasmbl.fasta.meryl"

cp -r $DATADIR/$PACBIO_MERYL $SCRATCHDIR/ || exit 1
cp -r $DATADIR/$AA_Z_MERYL $SCRATCHDIR/ || exit 2
cp -r $DATADIR/$AA_Z_W_MERYL $SCRATCHDIR/ || exit 3
cd $SCRATCHDIR || exit 4

export PATH=/storage/plzen1/home/jendrb00/meryl-1.4.1/bin:$PATH

mkdir -p meryl_results

echo "Porovnání PacBio dat s AA_Zchr..."
meryl difference $PACBIO_MERYL $AA_Z_MERYL output meryl_results/expected_W_size.meryl || exit 5

echo "Porovnání PacBio dat s AA_Z_W..."
meryl difference $PACBIO_MERYL $AA_Z_W_MERYL output meryl_results/W_missing_in_assembly.meryl || exit 6

echo "Počítám sumy unikátních k-merů pro očekávanou velikost W chromozomu..."
meryl print meryl_results/expected_W_size.meryl | awk '{sum += $2} END {print sum}' > meryl_results/expected_W_size_sum.txt || exit 7

echo "Počítám sumy chybějících unikátních k-merů W v sestaveném genomu..."
meryl print meryl_results/W_missing_in_assembly.meryl | awk '{sum += $2} END {print sum}' > meryl_results/W_missing_in_assembly_sum.txt || exit 8

mkdir -p $OUTPUT_DIR
cp -r meryl_results/* $OUTPUT_DIR || export CLEAN_SCRATCH=false || exit 9

echo "První část dokončena. Výsledky jsou v: $OUTPUT_DIR"
