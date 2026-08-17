#!/bin/bash
# Build a self-contained subset workspace from an existing full simulation.
#
# Nothing is re-simulated. Long reads and barcodes are already keyed by genome,
# so the subset is pure selection:
#
#   genome/      symlinks to the selected .fna
#   pacbio/      long reads filtered by header prefix
#   single_cell/ symlink to the parent (barcodes are filtered via truth/)
#   injected/    symlink to the parent
#   truth/       every table filtered to the selected genomes
#
# The result is a normal pipeline root, so every existing script runs against it
# unchanged with ROOT=<subset dir>.
#
# Usage: bash 11_build_subset.sh <subset_genomes.txt> <subset_dir>
set -euo pipefail

LIST=${1:?usage: 11_build_subset.sh <subset_genomes.txt> <subset_dir>}
DEST=${2:?usage: 11_build_subset.sh <subset_genomes.txt> <subset_dir>}

ROOT=${ROOT:-$(pwd)}
LIST=$(readlink -f "$LIST")
cd "$ROOT"
mkdir -p "$DEST"/{genome,pacbio,truth}
DEST=$(readlink -f "$DEST")

N=$(wc -l < "$LIST")
echo "[subset] $N genomes -> $DEST"

echo "[1/5] genome symlinks"
while read -r g; do
  [ -n "$g" ] || continue
  ln -sf "$ROOT/genome/${g}.fna" "$DEST/genome/${g}.fna"
done < "$LIST"
echo "  $(ls "$DEST/genome" | wc -l) linked"

echo "[2/5] long reads (filter by header prefix)"
zcat pacbio/pbsim_long_read.fastq.gz \
  | awk -v L="$LIST" '
      BEGIN { while ((getline g < L) > 0) if (g != "") keep_g[g]=1 }
      NR%4==1 { h=substr($1,2); sub(/_[0-9]+_[0-9]+$/,"",h); k=(h in keep_g) }
      k { print }' \
  | gzip > "$DEST/pacbio/pbsim_long_read.fastq.gz"
echo "  $(zcat "$DEST/pacbio/pbsim_long_read.fastq.gz" | awk 'END{print NR/4}') reads"

echo "[3/5] read symlinks"
ln -sfn "$ROOT/single_cell" "$DEST/single_cell"
[ -d "$ROOT/injected" ] && ln -sfn "$ROOT/injected" "$DEST/injected"

echo "[4/5] truth tables"
# barcode_truth: keep barcodes whose genome is in the subset
awk 'NR==FNR{g[$1];next} $3 in g' "$LIST" truth/barcode_truth.tsv \
  > "$DEST/truth/barcode_truth.tsv"
# contig_truth: MUST preserve per-genome order, so filter without sorting
awk -F'\t' 'NR==FNR{g[$1];next} $2 in g' "$LIST" truth/contig_truth.tsv \
  > "$DEST/truth/contig_truth.tsv"
# longread_truth: rebuild from the filtered reads
zcat "$DEST/pacbio/pbsim_long_read.fastq.gz" | awk 'NR%4==1' | sed 's/^@//' \
  | awk -v OFS='\t' '{acc=$1; sub(/_[0-9]+_[0-9]+$/,"",acc); print $1, acc}' \
  > "$DEST/truth/longread_truth.tsv"
# pair tables, strain groups, genome list
awk -F'\t' 'NR==FNR{g[$1];next} FNR==1 || ($1 in g && $2 in g)' \
  "$LIST" truth/genome_pairs.tsv > "$DEST/truth/genome_pairs.tsv"
awk -F'\t' 'NR==FNR{g[$1];next} FNR==1 || ($2 in g && $4 in g)' \
  "$LIST" truth/plasmid_pairs.tsv > "$DEST/truth/plasmid_pairs.tsv"
awk -F'\t' 'NR==FNR{g[$1];next} FNR==1 || $1 in g' \
  "$LIST" truth/strain_groups.tsv > "$DEST/truth/strain_groups.tsv"
awk -F'\t' 'NR==FNR{g[$1];next} FNR==1 || $1 in g' \
  "$LIST" truth/strain_nn.tsv > "$DEST/truth/strain_nn.tsv" 2>/dev/null || true
cp "$LIST" "$DEST/truth/genome_list.txt"
cp "$LIST" "$DEST/truth/pilot_genomes.txt"
# injected truth, if present, restricted to subset genomes
for f in truth/injected_*.tsv; do
  [ -e "$f" ] || continue
  awk -F'\t' 'NR==FNR{g[$1];next} FNR==1 || $5 in g' "$LIST" "$f" \
    > "$DEST/$f"
done

for f in "$DEST"/truth/*.tsv; do
  printf "  %-28s %s\n" "$(basename "$f")" "$(wc -l < "$f")"
done

echo "[5/5] contig truth round-trip"
python3 "$(dirname "$0")/07_longread_contigs.py" --truth "$DEST/truth"

cat <<EOF

Subset ready. Run the pipeline against it with:

  export ROOT=$DEST
  cd $DEST
  LR_SUBSAMPLE=0 BATCH=2000 bash $ROOT/scripts/02_link.sh 1x full
  python3 $ROOT/scripts/03_score.py --links links/1x --out links/1x/retained.tsv.gz --fdr 0.01
  python3 $ROOT/scripts/04_analyse.py --tiers 1x --linkdir links --truth truth --out results
  python3 $ROOT/scripts/08_plasmid_linkage.py --tier 1x --links links/1x/retained.tsv.gz --truth truth --out results
EOF
