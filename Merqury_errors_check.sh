#!/bin/bash
#PBS -N merqury_Errors_check
#PBS -l select=1:ncpus=8:mem=60gb:scratch_local=200gb
#PBS -l walltime=24:00:00
#PBS -j oe

trap 'clean_scratch' TERM EXIT

DATADIR="/storage/brno2/home/jendrb00/my_analysis/Ostrinia_nubilalis/creatig_my_datasets/ERR11867203.fastq.gz.meryl"
OUTPUT_DIR="/storage/brno2/home/jendrb00/my_analysis/Ostrinia_nubilalis/creatig_my_datasets/merqury_errors_look"
SCRATCHDIR=${SCRATCHDIR:-"/scratch"}

cp -r "$DATADIR" "$SCRATCHDIR" || exit 1
cd "$SCRATCHDIR" || exit 2

module add mambaforge
mamba activate /storage/plzen1/home/jendrb00/merqury-1.3/env_merqury

module add r
export R_LIBS_USER="/storage/brno2/home/jendrb00/Rpackages"

export PATH=/storage/plzen1/home/jendrb00/meryl-1.4.1/bin:$PATH

OUTPUT_NAME="merqury_ERRORS_PacBio_$(date +"%Y%m%d")"

merqury.sh "$(basename "$DATADIR")" "$(basename "$DATADIR")" "$OUTPUT_NAME" || exit 3

cp -r "$SCRATCHDIR"/* "$OUTPUT_DIR" || { export CLEAN_SCRATCH=false; exit 4; }
