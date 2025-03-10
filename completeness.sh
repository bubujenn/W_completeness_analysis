#!/bin/bash
#PBS -N calculate_W_completeness
#PBS -l select=1:ncpus=2:mem=4gb:scratch_local=10gb
#PBS -l walltime=01:00:00
#PBS -j oe

trap 'clean_scratch' TERM EXIT

DATADIR="/storage/brno2/home/jendrb00/my_analysis/Ostrinia_nubilalis/creatig_my_datasets"
RESULTDIR="/storage/brno2/home/jendrb00/my_analysis/Ostrinia_nubilalis/creatig_my_datasets"

cp $DATADIR/W_missing_in_assembly_sum.txt $SCRATCHDIR/ || exit 1
cp $DATADIR/expected_W_size_sum.txt $SCRATCHDIR/ || exit 2
cd $SCRATCHDIR || exit 3

expected=$(cat expected_W_size_sum.txt)
missing=$(cat W_missing_in_assembly_sum.txt)

# Výpočet procenta kompletnosti chromozomu W
completeness=$(awk -v missing="$missing" -v expected="$expected" 'BEGIN { printf "%.5f\n", (1 - (missing / expected)) * 100 }')
echo "$completeness" > W_FINAL_completeness.txt

echo "Procento kompletnosti chromozomu W:"
cat W_FINAL_completeness.txt

mkdir -p $RESULTDIR
cp W_FINAL_completeness.txt $RESULTDIR/ || export CLEAN_SCRATCH=false

echo "Výsledky kompletnosti chromozomu W jsou uloženy v: $RESULTDIR/W_FINAL_completeness.txt"
