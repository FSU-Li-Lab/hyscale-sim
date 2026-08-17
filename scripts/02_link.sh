#!/bin/bash
# Stage 3: link barcodes to long reads.
#
# Design notes:
#   - PAF not SAM. Only (query, target) is needed; alignment records are
#     discarded, so nothing but aggregated counts ever touches disk.
#   - minimap2 pipes straight into reduce_links.py.
#   - Barcodes are batched per index load: ~200 barcodes concatenated into one
#     minimap2 call, each read renamed <BARCODE>:<mate>:<n>. One index load
#     serves the whole batch instead of 200 loads.
#   - R1/R2 mapped as independent single-end reads. Pairing adds nothing to a
#     co-occurrence count and breaks the concatenation.
#   - Injected doublet/ambient barcodes (injected/<tier>/) are included
#     automatically. They are the negative class; without them FDR is
#     undefined.
#
# Usage: bash 02_link.sh <tier> [pilot|full]
set -euo pipefail

TIER=${1:?usage: 02_link.sh <tier> [pilot|full]}
SCOPE=${2:-full}

ROOT=${ROOT:-$(pwd)}
THREADS=${THREADS:-32}
BATCH=${BATCH:-200}
MIN_HITS=${MIN_HITS:-1}
LR_SUBSAMPLE=${LR_SUBSAMPLE:-0.30}   # 0 = index all long reads
MASK=${MASK:-}                       # optional BED of repeat intervals

cd "$ROOT"
SCRIPTS=$(cd "$(dirname "$0")" && pwd)
OUT="links/$TIER"
mkdir -p "$OUT" work

# ---------------------------------------------------------------- long reads
LR=work/longreads.fa
IDX=work/longreads.mmi

if [ ! -s "$IDX" ]; then
  echo "[index] preparing long-read reference"
  if [ "$LR_SUBSAMPLE" != "0" ]; then
    echo "  subsampling long reads to $LR_SUBSAMPLE (linkage does not need 50x;"
    echo "  the full set is still used for assembly in stage 7)"
    seqkit sample -p "$LR_SUBSAMPLE" -s 42 pacbio/pbsim_long_read.fastq.gz \
      | seqkit fq2fa > "$LR"
  else
    seqkit fq2fa pacbio/pbsim_long_read.fastq.gz > "$LR"
  fi
  grep -c '>' "$LR" | xargs echo "  long reads indexed:"
  minimap2 -x sr -t "$THREADS" -I 16G -d "$IDX" "$LR"
fi

# ---------------------------------------------------------------- barcodes
if [ "$SCOPE" = "pilot" ]; then
  [ -s truth/pilot_genomes.txt ] || { echo "missing truth/pilot_genomes.txt" >&2; exit 1; }
  awk -v t="$TIER" 'NR==FNR{g[$1];next} $1==t && $3 in g {print $2}' \
    truth/pilot_genomes.txt truth/barcode_truth.tsv > "$OUT/barcodes.txt"
else
  awk -v t="$TIER" '$1==t {print $2}' truth/barcode_truth.tsv > "$OUT/barcodes.txt"
fi
N_CLEAN=$(wc -l < "$OUT/barcodes.txt")

# injected negatives for this tier
N_INJ=0
INJ_TRUTH="truth/injected_${TIER}.tsv"
if [ -s "$INJ_TRUTH" ]; then
  if [ "$SCOPE" = "pilot" ]; then
    # keep only injections whose major genome is in the pilot set
    awk 'NR==FNR{g[$1];next} FNR>1 && $5 in g {print $1}' \
      truth/pilot_genomes.txt "$INJ_TRUTH" >> "$OUT/barcodes.txt"
  else
    awk 'FNR>1 {print $1}' "$INJ_TRUTH" >> "$OUT/barcodes.txt"
  fi
  N_INJ=$(( $(wc -l < "$OUT/barcodes.txt") - N_CLEAN ))
else
  echo "  WARNING: no $INJ_TRUTH -- running without a negative class."
  echo "           FDR will be undefined. Run 06_inject_negatives.py first."
fi

echo "[link] tier=$TIER scope=$SCOPE clean=$N_CLEAN injected=$N_INJ batch=$BATCH"

rm -f "$OUT"/batch_[0-9]*
split -l "$BATCH" -d -a 4 "$OUT/barcodes.txt" "$OUT/batch_"

# resolve a barcode+mate to its FASTQ, injected first
resolve() {
  local bc=$1 mate=$2
  if [ -s "injected/$TIER/${bc}${mate}.fq.gz" ]; then
    echo "injected/$TIER/${bc}${mate}.fq.gz"
  elif [ -s "single_cell/${bc}${mate}.fq.gz" ]; then
    echo "single_cell/${bc}${mate}.fq.gz"
  fi
}
export -f resolve

REDUCE_ARGS="--min-hits $MIN_HITS"
[ -n "$MASK" ] && REDUCE_ARGS="$REDUCE_ARGS --mask $MASK"

# ---------------------------------------------------------------- run
for B in "$OUT"/batch_[0-9]*; do
  TAG=$(basename "$B")
  DEST="$OUT/${TAG}.links.tsv.gz"
  [ -s "$DEST" ] && { echo "  $TAG done, skipping"; continue; }
  echo "  $TAG"

  {
    while read -r bc; do
      for mate in 1 2; do
        f=$(TIER="$TIER" resolve "$bc" "$mate")
        [ -n "$f" ] || continue
        zcat "$f" | awk -v b="$bc" -v m="$mate" '
          NR%4==1 { n++; print "@" b ":" m ":" n; next }
          { print }'
      done
    done < "$B"
  } | minimap2 -x sr -t "$THREADS" --secondary=yes -N 300 -p 0.6 "$IDX" - 2>/dev/null \
    | python3 "$SCRIPTS/reduce_links.py" $REDUCE_ARGS \
    | gzip > "$DEST.tmp" && mv "$DEST.tmp" "$DEST"
done

echo "[link] done. per-batch tables in $OUT/"
zcat "$OUT"/batch_*.links.tsv.gz | wc -l | xargs echo "  total (bc,longread) pairs:"
echo
echo "Next: python3 scripts/03_score.py --links $OUT --out $OUT/retained.tsv.gz"
