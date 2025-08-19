#!/bin/bash
#PBS -l select=1:ncpus=8:mem=64gb:scratch_local=200gb
#PBS -l walltime=12:00:00
#PBS -N repeats_kmer_lookup
#PBS -j oe

trap 'clean_scratch' TERM EXIT

#nastaveni vstupu a promennych
WORKDIR="${PWD}"
CLUSTERS="cluster_list.txt"
CONTIGS="contigs.fasta"
TADEAM_DIR="tandem_consensi"

MERYL_BIN="/storage/plzen1/home/jendrb00/meryl-1.4.1/bin"
export PATH="$MERYL_BIN:$PATH"

MERYL_READS="reads.meryl"
MERYL_AAZ="AAZ.meryl"
MERYL_AAZW="AAZW.meryl"

KMER=63
DO_DIMER=1
PB_PEAK=45
EPS=1e-9

TS=$(date +"%Y%m%d_%H%M")
OUTDIR="$WORKDIR/repeats_k${KMER}_$TS"
OUT_TABLE="$OUTDIR/repeats_k${KMER}_summary.tsv"
OUT_DIMER="$OUTDIR/repeats_dimer_profile.tsv"

#kopie na scratch
mkdir -p "$SCRATCHDIR/repeat_analysis"
cp "$WORKDIR/$CLUSTERS" "$SCRATCHDIR/repeat_analysis/"
cp "$WORKDIR/$CONTIGS" "$SCRATCHDIR/repeat_analysis/"
cp -r "$WORKDIR/$TADEAM_DIR" "$SCRATCHDIR/repeat_analysis/" || true
cp -r "$WORKDIR/$MERYL_READS" "$SCRATCHDIR/repeat_analysis/"
cp -r "$WORKDIR/$MERYL_AAZ" "$SCRATCHDIR/repeat_analysis/"
cp -r "$WORKDIR/$MERYL_AAZW" "$SCRATCHDIR/repeat_analysis/"
cd "$SCRATCHDIR/repeat_analysis"

TADEAM_DIR="$(basename "$TADEAM_DIR")"
MERYL_READS="$(basename "$MERYL_READS")"
MERYL_AAZ="$(basename "$MERYL_AAZ")"
MERYL_AAZW="$(basename "$MERYL_AAZW")"

#hlavicky vystupu
mkdir -p "$OUTDIR"
echo -e "Cluster\tNumKmers\tPB_sum\tAAZ_sum\tAAZW_sum\tPB_norm\tAAZW_minus_AAZ\tW_ratio" > "$OUT_TABLE"
if [[ "$DO_DIMER" -eq 1 ]]; then
  echo -e "Cluster\tA\tC\tG\tT\tAA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT" > "$OUT_DIMER"
fi

#pomocne funkce
sum_lookup_counts() {
  awk 'NF==2 {s+=$2; next} NF==1 {s+=$1; next} END{print (s+0)}'
}

pb_norm() {
  awk -v cov="$PB_PEAK" '{printf "%.6f", $1/cov}'
}

dimer_profile() {
  local fasta="$1"
  meryl count k=2 output tmp_dimer.meryl "$fasta"
  {
    for b in A C G T; do echo "$b"; done
    for x in A C G T; do for y in A C G T; do echo "${x}${y}"; done; done
  } | while read -r k; do
        echo -e ">$k\n$k" | meryl-lookup -dump -sequence - -mers tmp_dimer.meryl | sum_lookup_counts
      done | paste -sd$'\t' -
  rm -rf tmp_dimer.meryl
}

#for loop pres klastry
while read -r cid; do
  [[ -z "$cid" ]] && continue
  echo ">>> delam cluster${cid}"

  # 4.1 hledam TADEAM nebo contig
  if [[ -f "${TADEAM_DIR}/cluster${cid}_TADEAM.consensus.fasta" ]]; then
    seqfile="${TADEAM_DIR}/cluster${cid}_TADEAM.consensus.fasta"
    source_type="TADEAM"
  else
    seqfile="cluster${cid}.fa"
    grep -A1 -w ">cluster${cid}" "$CONTIGS" > "$seqfile" || true
    source_type="CONTIG"
  fi

  #kdyz neni sekvence tak preskocim
  if [[ ! -s "$seqfile" ]]; then
    echo "!! nic pro cluster${cid}"
    continue
  fi

  #spocitam kmery z celeho souboru
  meryl count k=$KMER memory=4 output "tmp_${cid}.meryl" "$seqfile"

  #vytvorim fasta se seznamem kmeru
  meryl print "tmp_${cid}.meryl" | awk '{print ">"NR"\n"$1}' > "kmers_${cid}.fa"
  nkmers=$(grep -c '^>' "kmers_${cid}.fa")

  #lookup v databazich
  pb_sum=$(meryl-lookup -dump -sequence "kmers_${cid}.fa" -mers "$MERYL_READS" | sum_lookup_counts)
  aaz_sum=$(meryl-lookup -dump -sequence "kmers_${cid}.fa" -mers "$MERYL_AAZ"   | sum_lookup_counts)
  aazw_sum=$(meryl-lookup -dump -sequence "kmers_${cid}.fa" -mers "$MERYL_AAZW" | sum_lookup_counts)

  #normalizace pacbio
  pb_norm_val=$(printf "%s\n" "$pb_sum" | pb_norm)

  #rozdil a pomer pro W
  delta_w=$(awk -v x="$aazw_sum" -v y="$aaz_sum" 'BEGIN{printf "%.6f", x - y}')
  w_ratio=$(awk -v x="$aazw_sum" -v y="$aaz_sum" -v e="$EPS" 'BEGIN{printf "%.6f", (x+e)/(y+e)}')

  #zapisu do tabulky
  echo -e "cluster${cid}\t${nkmers}\t${pb_sum}\t${aaz_sum}\t${aazw_sum}\t${pb_norm_val}\t${delta_w}\t${w_ratio}" >> "$OUT_TABLE"

  #kdyz mam TADEAM udelam dimery
  if [[ "$DO_DIMER" -eq 1 && "$source_type" == "TADEAM" ]]; then
    line=$(dimer_profile "$seqfile")
    echo -e "cluster${cid}\t${line}" >> "$OUT_DIMER"
  fi

  #uklid
  rm -rf "tmp_${cid}.meryl" "kmers_${cid}.fa"
  if [[ "$source_type" == "CONTIG" ]]; then rm -f "$seqfile"; fi

  echo "✓ cluster${cid} hotovo"
done < "$CLUSTERS"

#zkopiruju vysledky zpet
cp "$OUT_TABLE" "$OUTDIR/"
if [[ "$DO_DIMER" -eq 1 && -f "$OUT_DIMER" ]]; then
  cp "$OUT_DIMER" "$OUTDIR/"
fi

echo "HOTOVO vse ulozeno do $OUTDIR"
