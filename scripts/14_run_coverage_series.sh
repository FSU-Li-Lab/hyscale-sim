#!/bin/bash
# Complete the coverage series: inject negatives, link, score and analyse the
# remaining tiers with parameters identical to the tuned 0.01x/0.1x/1x runs.
#
# Parameters are fixed here deliberately. The coverage series is only
# interpretable if every tier is scored the same way -- a per-tier threshold
# would confound coverage with calibration.
#
#   --fdr 1e-8 --min-reads 5     tuned threshold (see threshold sweep)
#   LR_SUBSAMPLE=0               index the full HiFi set
#   BATCH=2000                   amortise index loading
#
# Resumable: any tier whose retained.tsv.gz already exists is skipped.
#
# Usage:
#   bash 14_run_coverage_series.sh                      # 0.2x 0.4x 0.8x
#   bash 14_run_coverage_series.sh 0.2x 0.4x            # named tiers
set -euo pipefail

TIERS=${@:-"0.2x 0.4x 0.8x"}

ROOT=${ROOT:-$(pwd)}
SCRIPTS=${SCRIPTS:-$(cd "$(dirname "$0")" && pwd)}
FDR=${FDR:-1e-8}
MINREADS=${MINREADS:-5}
export LR_SUBSAMPLE=${LR_SUBSAMPLE:-0}
export BATCH=${BATCH:-2000}
export THREADS=${THREADS:-32}

cd "$ROOT"
echo "[series] root=$ROOT"
echo "[series] tiers: $TIERS"
echo "[series] threshold: --fdr $FDR --min-reads $MINREADS"

for T in $TIERS; do
  echo
  echo "================== tier $T =================="

  # barcodes for this tier must exist in the subset truth table
  NBC=$(awk -v t="$T" '$1==t' truth/barcode_truth.tsv | wc -l)
  if [ "$NBC" -eq 0 ]; then
    echo "  no barcodes for tier $T in truth/barcode_truth.tsv -- skipping" >&2
    continue
  fi
  echo "  $NBC clean barcodes"

  # 1. negatives
  if [ ! -s "truth/injected_${T}.tsv" ]; then
    python3 "$SCRIPTS/06_inject_negatives.py" --root . --tier "$T" \
      --doublet-rate 0.05 --ambient-rate 0.05 --ambient-fracs 0.01 0.05 0.10
  else
    echo "  [inject] truth/injected_${T}.tsv exists, skipping"
  fi

  # 2. link
  if ! ls "links/$T"/*.links.tsv.gz >/dev/null 2>&1; then
    bash "$SCRIPTS/02_link.sh" "$T" full
  else
    echo "  [link] batch tables exist, skipping"
  fi

  # 3. score
  if [ ! -s "links/$T/retained.tsv.gz" ]; then
    python3 "$SCRIPTS/03_score.py" --links "links/$T" \
      --out "links/$T/retained.tsv.gz" --fdr "$FDR" --min-reads "$MINREADS"
  else
    echo "  [score] retained.tsv.gz exists, skipping"
  fi

  # 4. M5.1 for this tier
  python3 "$SCRIPTS/08_plasmid_linkage.py" --tier "$T" \
    --links "links/$T/retained.tsv.gz" --truth truth --out results_series
done

# ---------------------------------------------------------------- combined
ALL=""
for T in 0.01x 0.1x 0.2x 0.4x 0.8x 1x; do
  [ -s "links/$T/retained.tsv.gz" ] && ALL="$ALL $T"
done
echo
echo "================== combined analysis =================="
echo "  tiers with results:$ALL"
python3 "$SCRIPTS/04_analyse.py" --tiers $ALL \
  --linkdir links --truth truth --out results_series

# ---------------------------------------------------------------- summary
echo
echo "=== coverage series: plasmid-to-host assignment ==="
printf "%-8s %10s %10s %10s %12s\n" tier assigned pct acc_adj chrom_acc
for T in $ALL; do
  f="results_series/M5.1_plasmid_linkage_${T}.tsv"
  [ -s "$f" ] || continue
  awk -F'\t' -v t="$T" '
    $2=="all_reads" && $3=="plasmid"    { pa=$4; acc=$9 }
    $2=="all_reads" && $3=="chromosome" { ca=$9 }
    END { printf "%-8s %10d %10s %10.4f %12.4f\n", t, pa, "-", acc, ca }' "$f"
done
echo
echo "Full tables in results_series/"
