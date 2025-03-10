#!/bin/bash
#PBS -l select=1:ncpus=8:mem=32gb:scratch_local=500gb
#PBS -l walltime=12:00:00
#PBS -N Minimap2_Samtools
#PBS -j oe

trap 'clean_scratch' TERM EXIT

DATADIR="/storage/plzen1/home/jendrb00/data"
RESULTDIR="/storage/plzen1/home/jendrb00/results"
SCRATCHDIR=${SCRATCHDIR:-"/scratch"}

cp $DATADIR/reference.fasta $SCRATCHDIR || exit 1
cp $DATADIR/query.fasta $SCRATCHDIR || exit 2
cd $SCRATCHDIR || exit 3

module load minimap2
module load samtools

minimap2 -ax asm5 --secondary=no reference.fasta query.fasta > alignment.sam || exit 4

samtools view -F 2048 -h alignment.sam > filtered.sam || exit 5
samtools sort -o filtered.bam filtered.sam || exit 6
samtools index filtered.bam || exit 7
samtools view filtered.bam | awk '{print $3, $1}' | sort | uniq > scaffold_to_chr.tsv || exit 8

cp -r $SCRATCHDIR/* $RESULTDIR || export CLEAN_SCRATCH=false || exit 9
