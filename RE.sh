#!/bin/bash
# =============================================================================
# RepeatExplorer Pipeline for W Chromosome Completeness Analysis
# =============================================================================

# =============================================================================
# PART 1: SIMULATION (randomreads.sh)
# =============================================================================
# Parametry: len=100, reads=20M, seed=42, paired=t, adderrors=f

cd /storage/brno2/home/jendrb00/20_10_consultation

# --- AAZ (genom samce, BEZ W chromozomu) ---
randomreads.sh \
    ref=AA_Z_chr.fasta \
    out1=AAZ_R1.fastq.gz \
    out2=AAZ_R2.fastq.gz \
    len=100 \
    reads=20000000 \
    seed=42 \
    paired=t \
    adderrors=f

# --- AAZW (genom samice, S W chromozomem) ---
randomreads.sh \
    ref=AA_Z_Wasmbl.fasta \
    out1=AAZW_R1.fastq.gz \
    out2=AAZW_R2.fastq.gz \
    len=100 \
    reads=20000000 \
    seed=42 \
    paired=t \
    adderrors=f

# --- PBrw (PacBio data samice) ---
# Krok 1: 2× vzorek z PacBio (14M readů = ~1.4 Gb)
seqtk sample -s42 pb_100_se.fastq.gz 14000000 > PB_2xgenome_new.fastq
seqtk seq -a PB_2xgenome_new.fastq > PB_2xgenome_new.fasta

# Krok 2: Simulace
randomreads.sh \
    ref=PB_2xgenome_new.fasta \
    out1=PBrw_R1_new.fastq.gz \
    out2=PBrw_R2_new.fastq.gz \
    len=100 \
    reads=20000000 \
    seed=42 \
    paired=t \
    adderrors=f

# =============================================================================
# PART 2: SUBSAMPLE (normalizace na 0.4× coverage)
# =============================================================================
# Počty: AAZ 438646 párů, AAZW/PBrw 449857 párů

cd /storage/brno2/home/jendrb00/20_10_consultation/RE_final_inputs

# AAZ (438646 párů)
seqtk sample -s42 ../AAZ_R1.fastq.gz 438646 > AAZ_R1.sampled.fq
seqtk sample -s42 ../AAZ_R2.fastq.gz 438646 > AAZ_R2.sampled.fq

# AAZW (449857 párů)
seqtk sample -s42 ../AAZW_R1.fastq.gz 449857 > AAZW_R1.sampled.fq
seqtk sample -s42 ../AAZW_R2.fastq.gz 449857 > AAZW_R2.sampled.fq

# PBrw (449857 párů)
seqtk sample -s42 ../PBrw_R1_new.fastq.gz 449857 > PBrw_R1.sampled.fq
seqtk sample -s42 ../PBrw_R2_new.fastq.gz 449857 > PBrw_R2.sampled.fq

# =============================================================================
# PART 3: INTERLEAVE + FASTA
# =============================================================================

seqtk mergepe AAZ_R1.sampled.fq AAZ_R2.sampled.fq > AAZ_interleaved.fastq
seqtk seq -a AAZ_interleaved.fastq > AAZ_interleaved_new.fa

seqtk mergepe AAZW_R1.sampled.fq AAZW_R2.sampled.fq > AAZW_interleaved.fastq
seqtk seq -a AAZW_interleaved.fastq > AAZW_interleaved_new.fa

seqtk mergepe PBrw_R1.sampled.fq PBrw_R2.sampled.fq > PBrw_interleaved.fastq
seqtk seq -a PBrw_interleaved.fastq > PBrw_interleaved_new.fa

# =============================================================================
# PART 4: ZKRÁCENÍ NÁZVŮ + TAGOVÁNÍ (v jednom kroku)
# =============================================================================

awk '/^>/{split($0,a," "); print a[1]; next} {print}' AAZ_interleaved_new.fa | sed 's/>/>AAZ0_/' > AAZ0_tagged.fa
awk '/^>/{split($0,a," "); print a[1]; next} {print}' AAZW_interleaved_new.fa | sed 's/>/>AAZW_/' > AAZW_tagged.fa
awk '/^>/{split($0,a," "); print a[1]; next} {print}' PBrw_interleaved_new.fa | sed 's/>/>PBrw_/' > PBrw_tagged.fa

# =============================================================================
# PART 5: CONCATENACE
# =============================================================================

cat AAZ0_tagged.fa PBrw_tagged.fa > AAZ0_vs_PBrw.fa   # Složení repetic W
cat AAZW_tagged.fa PBrw_tagged.fa > AAZW_vs_PBrw.fa   # Kompletnost W

# Ověření
grep -c '^>' AAZ0_vs_PBrw.fa   # 1777006
grep -c '^>' AAZW_vs_PBrw.fa   # 1799428

# =============================================================================
# PART 6: SPUŠTĚNÍ RepeatExplorer
# =============================================================================

# Dataset 1: AAZ0 vs PBrw
bash /storage/plzen1/home/p817n421/Shared/blues/scripts/RE_clustering.sh AAZ0_vs_PBrw.fa 5

# Dataset 2: AAZW vs PBrw
bash /storage/plzen1/home/p817n421/Shared/blues/scripts/RE_clustering.sh AAZW_vs_PBrw.fa 5
