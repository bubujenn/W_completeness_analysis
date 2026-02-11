#!/bin/bash
#PBS -N RE_RUN_AAZ0_v2
#PBS -l select=1:ncpus=20:mem=1024gb:scratch_ssd=500gb:singularity=True
#PBS -l walltime=144:00:00
#PBS -q large_mem@pbs-m1.metacentrum.cz
#PBS -M blanka.jendriskova0000@seznam.cz
#PBS -m bea
#PBS -o /storage/brno2/home/jendrb00/20_10_consultation/RE_final_inputs/RE_RUN_AAZ0_v2_Output/stdout.txt
#PBS -e /storage/brno2/home/jendrb00/20_10_consultation/RE_final_inputs/RE_RUN_AAZ0_v2_Output/stderr.txt

# ============================================================
# RepeatExplorer - AAZ0 vs PBrw (složení repetic W)
# Abraxas sylvata - Blanka Jendrisková
# Opravená verze 2 (11.2.2026)
# Změna: scratch_ssd místo scratch_local (řeší disk I/O error)
# ============================================================

trap 'clean_scratch' TERM EXIT

INPUT_FILE="AAZ0_vs_PBrw.fa"
INPUT_DIR="/storage/brno2/home/jendrb00/20_10_consultation/RE_final_inputs"
OUTPUT_DIR="${INPUT_DIR}/RE_RUN_AAZ0_v2_Output"
REPO_IMAGE="library://repeatexplorer/default/repex_tarean:0.3.12-7a7dc9e"

# Vytvoř výstupní složku
mkdir -p "$OUTPUT_DIR"

# 1. Kopie vstupních dat na scratch SSD
echo "=== KROK 1: Kopíruji data na scratch ==="
echo "Scratch: $SCRATCHDIR"
df -h "$SCRATCHDIR"
cp "${INPUT_DIR}/${INPUT_FILE}" "$SCRATCHDIR" || exit 1
cd "$SCRATCHDIR"

# 2. Stažení Singularity image
echo "=== KROK 2: Stahuji Singularity image ==="
export SINGULARITY_CACHEDIR="$SCRATCHDIR"
export SINGULARITY_LOCALCACHEDIR="$SCRATCHDIR"
export SINGULARITY_TMPDIR="$SCRATCHDIR"
export TMPDIR="$SCRATCHDIR"
singularity pull --arch amd64 "$REPO_IMAGE"
echo "Stav disku po stažení image:"
df -h "$SCRATCHDIR"

# 3. Spuštění RepeatExplorer
echo "=== KROK 3: Spouštím RepeatExplorer ==="
echo "Vstup: ${INPUT_FILE}, prefix_length=5"
echo "Start: $(date)"

singularity exec --bind "$SCRATCHDIR:/data" "$SCRATCHDIR/repex_tarean_0.3.12-7a7dc9e.sif" \
    seqclust --paired \
    --automatic_filtering \
    --mincl 0.001 \
    --cpu 18 \
    --prefix_length 5 \
    --output_dir=/data/output/ \
    --max_memory 900000000 \
    --cleanup \
    --assembly_min 3 \
    --taxon METAZOA3.0 \
    "/data/${INPUT_FILE}"

SEQCLUST_EXIT=$?
echo "Seqclust skončil s kódem: $SEQCLUST_EXIT"
echo "Konec: $(date)"
df -h "$SCRATCHDIR"

# 4. Kopie výsledků zpět
echo "=== KROK 4: Kopíruji výsledky ==="
if [ -d "$SCRATCHDIR/output" ]; then
    cp -r "$SCRATCHDIR/output"/* "$OUTPUT_DIR" || export CLEAN_SCRATCH=false
    echo "Výsledky úspěšně zkopírovány do: $OUTPUT_DIR"
else
    echo "CHYBA: složka output nenalezena!"
    export CLEAN_SCRATCH=false
fi
echo "=== HOTOVO ==="
