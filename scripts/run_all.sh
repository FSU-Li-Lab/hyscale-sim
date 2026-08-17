#!/bin/bash
# Driver. Run the pilot first; only scale once E1/E2 look sane on 10 genomes.
#
#   bash run_all.sh pilot     # ~1 hour, 10 genomes, 1x tier
#   bash run_all.sh full      # the coverage ladder on all 473
#
set -euo pipefail

MODE=${1:-pilot}
ROOT=${ROOT:-$(pwd)}
SCRIPTS=$(cd "$(dirname "$0")" && pwd); [ -d "$SCRIPTS/scripts" ] && SCRIPTS="$SCRIPTS/scripts"
export ROOT

cd "$ROOT"

# ---------------------------------------------------------------- truth
if [ ! -s truth/strain_groups.tsv ]; then
  bash "$SCRIPTS/01_build_truth.sh"
else
  echo "[truth] already built, skipping (rm truth/strain_groups.tsv to redo)"
fi

# ---------------------------------------------------------------- tiers
if [ "$MODE" = "pilot" ]; then
  TIERS="1x"
  SCOPE=pilot
else
  # 10x is the saturated control and is ~400 GB of the 500 GB; run it on a
  # subset separately rather than as part of the ladder
  TIERS="0.01x 0.1x 0.2x 0.4x 0.8x 1x"
  SCOPE=full
fi

for T in $TIERS; do
  echo
  echo "================ tier $T ================"

  if [ ! -s "truth/injected_${T}.tsv" ]; then
    python3 "$SCRIPTS/06_inject_negatives.py" --root . --tier "$T" \
      --doublet-rate 0.05 --ambient-rate 0.05 --ambient-fracs 0.01 0.05 0.10
  else
    echo "[inject] truth/injected_${T}.tsv exists, skipping"
  fi

  bash "$SCRIPTS/02_link.sh" "$T" "$SCOPE"

  python3 "$SCRIPTS/03_score.py" \
    --links "links/$T" --out "links/$T/retained.tsv.gz" --fdr 0.01
done

# ---------------------------------------------------------------- analyse
python3 "$SCRIPTS/04_analyse.py" --tiers $TIERS \
  --linkdir links --truth truth --out results

echo
echo "results/ written. Key files:"
ls -la results/
