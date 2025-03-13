#!/bin/bash
#PBS -N calculate_W_completeness
#PBS -l select=1:ncpus=2:mem=4gb:scratch_local=10gb
#PBS -l walltime=01:00:00
#PBS -j oe

trap 'clean_scratch' TERM EXIT

DATADIR="/storage/brno2/home/jendrb00/my_analysis/Ostrinia_nubilalis/creatig_my_datasets"
OUTPUT_DIR="/storage/brno2/home/jendrb00/my_analysis/Ostrinia_nubilalis/creatig_my_datasets"

cp $DATADIR/W_missing_in_assembly_sum.txt $SCRATCHDIR/ || exit 1
cp $DATADIR/expected_W_size_sum.txt $SCRATCHDIR/ || exit 2
cd $SCRATCHDIR || exit 3

expected=$(cat expected_W_size_sum.txt)
missing=$(cat W_missing_in_assembly_sum.txt)

completeness=$(awk -v missing="$missing" -v expected="$expected" 'BEGIN { printf "%.5f\n", (1 - (missing / expected)) * 100 }')  || exit 4
echo "$completeness" > W_FINAL_completeness.txt

cat W_FINAL_completeness.txt  || exit 5

mkdir -p $OUTPUT_DIR
cp W_FINAL_completeness.txt $OUTPUT_DIR/ || export CLEAN_SCRATCH=false  || exit 6

echo "Výsledky completeness chromozomu W jsou uloženy v: $OUTPUT_DIR/W_FINAL_completeness.txt" 
