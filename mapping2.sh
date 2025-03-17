#!/bin/bash
#PBS -l select=1:ncpus=8:mem=32gb:scratch_local=500gb
#PBS -l walltime=12:00:00
#PBS -N Minimap2_Samtools_Bedtools
#PBS -j oe

# Clean scratch space on exit
trap 'clean_scratch' TERM EXIT

# Nastavení cest a proměnných
Path=${PWD}
Base="Minimap2_Samtools"
TimeStamp=$(date +"%Y%m%d.%H%M")
Outdir=${Path}/${Base}_${TimeStamp}
SCRATCHDIR=${SCRATCHDIR:-"/scratch"}

# Vytvoření výstupního adresáře
mkdir -p "$Outdir"

# Vstupní soubory 
REFERENCE="${1:-reference.fasta}"
QUERY="${2:-query.fasta}"

# Kopírování vstupních souborů do scratch prostoru
cp "$Path/$REFERENCE" "$SCRATCHDIR/" || exit 1
cp "$Path/$QUERY" "$SCRATCHDIR/" || exit 2
cd "$SCRATCHDIR" || exit 3

# Načtení požadovaných modulů
module purge 
module load python/3.9.12-gcc-10.2.1
module load samtools/1.13-gcc-10.2.1
module load minimap2/2.22-gcc-10.2.1
module load bedtools

# Spuštění Minimap2 zarovnání
minimap2 -ax asm5 --secondary=no "$REFERENCE" "$QUERY" > alignment.sam || exit 4

# Zpracování zarovnání pomocí Samtools
samtools view -F 2048 -h alignment.sam > filtered.sam || exit 5
samtools sort -o filtered.bam filtered.sam || exit 6
samtools index filtered.bam || exit 7

# Vytvoření scaffold-to-chromosome mapování
samtools view filtered.bam | awk '{print $3, $1}' | sort | uniq > scaffold_to_chr.tsv || exit 8

# Výpočet coverage pomocí bedtools
bedtools genomecov -ibam filtered.bam -d > bedtools_coverage.tsv || exit 9

# Přesunutí výsledků do výstupního adresáře
cp -r "$SCRATCHDIR"/* "$Outdir/" || { export CLEAN_SCRATCH=false; exit 10; }

# Finální zpráva s umístěním výsledků
echo "Done! Results saved in: $Outdir"
