#!/bin/bash
#PBS -N bcalm_missing_seq
#PBS -l select=1:ncpus=8:mem=60gb:scratch_local=200gb
#PBS -l walltime=12:00:00
#PBS -j oe

# Clean scratch on exit
trap 'clean_scratch' TERM EXIT

# PATHS
Path="/storage/brno2/home/jendrb00/DIPLOMOVA_PRACE_VYSLEDKY/my_analysis/Abraxas_sylvata_merfin/AA_Z_W/try2/bcalm"
Base="bcalm_missing_seq"
TimeStamp=$(date +"%Y%m%d.%H%M")
Outdir="${Path}/${Base}_${TimeStamp}"
SCRATCHDIR="${SCRATCHDIR:-"/scratch"}"

MERFIN_STATS="merfin_stats.tsv"
KMER_SIZE=19
BCALM_BIN="/storage/brno2/home/jendrb00/bcalm/build/bcalm"

#Prepare scratch space 
mkdir -p "$SCRATCHDIR/bcalm_run"
cp "$Path/$MERFIN_STATS" "$SCRATCHDIR/bcalm_run/" || exit 1
cd "$SCRATCHDIR/bcalm_run" || exit 2
mkdir -p "$Outdir"

#Filter problematic k-mers
echo ":small_blue_diamond: Filtering missing and collapsed k-mers..."
awk -F'\t' '$0 ~ /missing|collapsed/ {print $1}' "$MERFIN_STATS" > problem_kmers.txt
if [[ ! -s problem_kmers.txt ]]; then
  echo ":warning: No problematic k-mers found. Exiting."
  exit 0
fi

#Convert to FASTA 
echo ":small_blue_diamond: Converting k-mers to FASTA format..."
awk '{print ">kmer" NR "\n" $1}' problem_kmers.txt > problem_kmers.fasta

#Run BCALM2 
echo ":small_blue_diamond: Running BCALM2..."
"$BCALM_BIN" \
  -in problem_kmers.fasta \
  -kmer-size "$KMER_SIZE" \
  -abundance-min 1 \
  -out unitigs

UNITIG_FASTA="unitigs.unitigs.fa"
if [[ ! -f "$UNITIG_FASTA" ]]; then
  echo ":x: Error: BCALM2 failed to produce output!"
  exit 1
fi

#Summarize unitigs
echo ":small_blue_diamond: Summarizing unitig lengths..."
TOTAL_BP=$(awk '/^>/ {next} {sum += length($0)} END {print sum}' "$UNITIG_FASTA")
NUM_UNITIGS=$(grep -c "^>" "$UNITIG_FASTA")

#Copy results back
cp -r * "$Outdir" || { export CLEAN_SCRATCH=false; exit 3; }

#Final report 
echo ":white_tick: Done!"
echo "- Output directory: $Outdir"
echo "- Number of unitigs: $NUM_UNITIGS"
echo "- Total bp in unitigs: $TOTAL_BP"
