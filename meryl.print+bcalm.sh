#!/bin/bash
#PBS -N missing_kmers_extract
#PBS -l select=1:ncpus=8:mem=40gb:scratch_local=200gb
#PBS -l walltime=04:00:00
#PBS -j oe

trap 'clean_scratch' TERM EXIT

# PATHS
WORKDIR="/storage/brno2/home/jendrb00/DIPLOMOVA_PRACE_VYSLEDKY/my_analysis/Abraxas_sylvata_merfin/AA_Z_W/try2"
OUTDIR="${WORKDIR}/bcalm"
TIMESTAMP=$(date +"%Y%m%d_%H%M")
SCRATCH_SUBDIR="$SCRATCHDIR/missing_kmers_${TIMESTAMP}"


READMERS="without_errors.meryl"
SEQMERS="AA_Z_W.meryl"
KMER_SIZE=19

# Meryl
export PATH="/storage/plzen1/home/jendrb00/meryl-1.4.1/bin:$PATH"

#SCRATCHDIR
mkdir -p "$SCRATCH_SUBDIR"
cd "$SCRATCH_SUBDIR" || exit 1

#COPY
cp -r "$WORKDIR/$READMERS" "$SCRATCH_SUBDIR"
cp -r "$WORKDIR/$SEQMERS" "$SCRATCH_SUBDIR"

#Extract and compare k-mers
echo "🧪 Vytahuju a třídím k-mery z reads..."
meryl print "$READMERS" | cut -f1 | sort -T . > reads_kmers.txt

echo "🧪 Vytahuju a třídím k-mery z assembly..."
meryl print "$SEQMERS" | cut -f1 | sort -T . > asm_kmers.txt

echo "🔍 Hledám chybějící k-mery..."
comm -23 reads_kmers.txt asm_kmers.txt > missing_kmers.txt

# FASTA
echo "📄 Převádím na FASTA..."
awk '{print ">kmer" NR "\n" $1}' missing_kmers.txt > missing_kmers.fasta

#Copy back
echo "💾 Kopíruju výsledky zpět do: $OUTDIR"
cp missing_kmers.txt missing_kmers.fasta "$OUTDIR" || { export CLEAN_SCRATCH=false; exit 2; }

#Shrnutí
TOTAL_KMERS=$(wc -l < missing_kmers.txt)
echo "Done. Počet chybějících k-merů: $TOTAL_KMERS"
