#!/bin/bash
#PBS -N merqury_cytogenetic_analysis
#PBS -l select=1:ncpus=8:mem=64gb:scratch_local=200gb
#PBS -l walltime=24:00:00
#PBS -j oe

trap 'clean_scratch' TERM EXIT

MERYL_FILE="reads.meryl"
ASSEMBLY_FILE="assembly.fasta"
OUTPUT_DIR="/storage/brno2/home/jendrb00/merqury__analysis"
SCRATCH_OUTDIR="$SCRATCHDIR/merqury_results_end"

cd $SCRATCHDIR || exit 1

cp $OUTPUT_DIR/$MERYL_FILE $SCRATCHDIR/ || exit 2
cp $OUTPUT_DIR/$ASSEMBLY_FILE $SCRATCHDIR/ || exit 3

export PATH=/storage/plzen1/home/jendrb00/meryl-1.4.1/bin:$PATH
module add mambaforge || exit 4
mamba activate /storage/plzen1/home/jendrb00/merqury-1.3/env_merqury || exit 5
module add r || exit 5
export R_LIBS_USER="/storage/plzen1/home/jendrb00/Rpackages"

merqury.sh $MERYL_FILE $ASSEMBLY_FILE merqury_output || exit 6

mkdir -p $OUTPUT_DIR/merqury_output
cp -r merqury_output/* $OUTPUT_DIR/merqury_output/ || exit 7

clean_scratch

echo "Analysis completed. Results are saved in: $OUTPUT_DIR/merqury_output2"
