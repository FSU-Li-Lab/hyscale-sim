#!/bin/bash
# Characterise the validation community: ANI, distance tree, plasmid clusters,
# and genome quality. Everything the community figure needs.
#
# Outputs (all under <outdir>/):
#   ani_matrix.tsv        36x36 skani ANI, full matrix
#   mash_dist.tsv         pairwise mash distances (used for the tree)
#   plasmid_pairs.tsv     skani all-vs-all over individual plasmid contigs
#   plasmid_groups.tsv    genome, contig, plasmid_group, length
#   checkm2/quality_report.tsv
#   organisms.tsv         genome, organism name, size, n_plasmids
#
# Usage: bash 12_characterise_subset.sh <genome_list> <genome_dir> <outdir>
set -euo pipefail

LIST=${1:?usage: 12_characterise_subset.sh <genome_list> <genome_dir> <outdir>}
GDIR=${2:?}
OUT=${3:?}

THREADS=${THREADS:-16}
PLASMID_ANI=${PLASMID_ANI:-99}      # ANI to call two plasmids the same element
PLASMID_AF=${PLASMID_AF:-80}        # aligned fraction required

mkdir -p "$OUT"
LIST=$(readlink -f "$LIST"); GDIR=$(readlink -f "$GDIR"); OUT=$(readlink -f "$OUT")

# a directory holding ONLY the subset genomes, for tools that take a folder
WORKG="$OUT/genomes"
mkdir -p "$WORKG"
while read -r g; do
  [ -n "$g" ] || continue
  ln -sf "$GDIR/${g}.fna" "$WORKG/${g}.fna"
done < "$LIST"
N=$(ls "$WORKG" | wc -l)
echo "[0] $N genomes linked into $WORKG"

# ---------------------------------------------------------------- 1. ANI
echo "[1] skani ANI matrix"
if [ ! -s "$OUT/ani_matrix.tsv" ]; then
  ( cd "$OUT" && skani triangle "$WORKG"/*.fna -t "$THREADS" --full-matrix \
      -o ani_matrix.tsv > skani.log 2>&1 )
fi
head -1 "$OUT/ani_matrix.tsv"

# ------------------------------------------------------- 2. mash distances
# skani reports nothing below ~80% ANI, so distant genera would have no
# distance at all. mash is defined across the whole range and is what the
# tree is built from; skani supplies the heatmap values where they exist.
echo "[2] mash distances"
if [ ! -s "$OUT/mash_dist.tsv" ]; then
  mash sketch -p "$THREADS" -s 20000 -k 21 -o "$OUT/subset" "$WORKG"/*.fna \
    > "$OUT/mash.log" 2>&1
  mash dist -p "$THREADS" "$OUT/subset.msh" "$OUT/subset.msh" \
    > "$OUT/mash_dist.tsv"
fi
wc -l < "$OUT/mash_dist.tsv" | xargs echo "  pairs:"

# ------------------------------------------------------- 3. plasmid clusters
echo "[3] plasmid extraction and clustering"
: > "$OUT/plasmids.fa"
while read -r g; do
  [ -n "$g" ] || continue
  awk -v g="$g" '
    /^>/ { p = (tolower($0) ~ /plasmid|extrachromosomal|megaplasmid/)
           if (p) { sub(/^>/, ">" g "|") } }
    p { print }' "$GDIR/${g}.fna"
done < "$LIST" >> "$OUT/plasmids.fa"
NP=$(grep -c '>' "$OUT/plasmids.fa" || true)
echo "  $NP plasmid contigs"

if [ "$NP" -gt 1 ] && [ ! -s "$OUT/plasmid_pairs.tsv" ]; then
  # -E / --sparse gives an edge list; without it skani writes a lower-triangle
  # matrix, which the parser below also handles for pre-existing files
  skani triangle -i "$OUT/plasmids.fa" -E -t "$THREADS" \
    -o "$OUT/plasmid_pairs.tsv" > "$OUT/skani_plasmid.log" 2>&1 || \
  skani triangle -i "$OUT/plasmids.fa" -t "$THREADS" \
    -o "$OUT/plasmid_pairs.tsv" >> "$OUT/skani_plasmid.log" 2>&1
fi

# single-linkage cluster plasmids at the chosen ANI/AF, then name groups
python3 - "$OUT" "$PLASMID_ANI" "$PLASMID_AF" <<'PY'
import sys, os
out, min_ani, min_af = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])

lens, order = {}, []
name = None
with open(os.path.join(out, "plasmids.fa")) as fh:
    for line in fh:
        if line.startswith(">"):
            name = line[1:].split()[0]
            order.append(name); lens[name] = 0
        elif name:
            lens[name] += len(line.strip())

parent = {n: n for n in order}
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]; x = parent[x]
    return x
def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb: parent[rb] = ra

def read_triangle(path):
    """skani lower-triangle matrix -> {(name_i, name_j): value}."""
    vals, names = {}, []
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]
    if not lines:
        return vals, names
    body = lines[1:] if lines[0].strip().isdigit() else lines
    for row in body:
        f = row.split("\t")
        names.append(f[0].split()[0])
        i = len(names) - 1
        for j, v in enumerate(f[1:]):
            if j >= i:
                break
            try:
                vals[(names[i], names[j])] = float(v)
            except ValueError:
                pass
    return vals, names


pp = os.path.join(out, "plasmid_pairs.tsv")
npair = 0
if os.path.exists(pp) and os.path.getsize(pp) > 0:
    with open(pp) as fh:
        first = fh.readline().rstrip("\n")
    is_matrix = first.strip().isdigit() or ("\t" in first and
                                            not first.lower().startswith(("ref", "query")))
    if is_matrix:
        # lower-triangle matrix, with aligned fractions in the sidecar .af file
        ani_v, _ = read_triangle(pp)
        af_v = {}
        if os.path.exists(pp + ".af"):
            af_v, _ = read_triangle(pp + ".af")
        for (a, b), ani in ani_v.items():
            af = max(af_v.get((a, b), 0.0), af_v.get((b, a), 0.0))
            if ani >= min_ani and (af >= min_af or not af_v):
                if a in parent and b in parent:
                    union(a, b); npair += 1
        print(f"  parsed lower-triangle matrix ({len(ani_v)} pairs scored)")
    else:
        with open(pp) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            def col(*names):
                for nm in names:
                    for i, h in enumerate(header):
                        if h.strip().lower().startswith(nm):
                            return i
                return None
            ci  = col("ani")
            cr  = col("align_fraction_ref", "align_fraction")
            cq  = col("align_fraction_query")
            cn1 = col("ref_name")
            cn2 = col("query_name")
            if ci is None or cn1 is None or cn2 is None:
                ci, cr, cq, cn1, cn2 = 2, 3, 4, 5, 6
                fh.seek(0)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) <= max(ci, cn1, cn2):
                    continue
                try:
                    ani = float(f[ci])
                    af = max(float(f[cr]) if cr is not None and cr < len(f) else 0.0,
                             float(f[cq]) if cq is not None and cq < len(f) else 0.0)
                except ValueError:
                    continue
                if ani >= min_ani and af >= min_af:
                    a = f[cn1].split()[0]; b = f[cn2].split()[0]
                    if a in parent and b in parent:
                        union(a, b); npair += 1
    print(f"  {npair} plasmid pairs at >={min_ani}% ANI / >={min_af}% AF")
else:
    print("  no plasmid comparison available -- every plasmid is its own group",
          file=sys.stderr)

groups = {}
for n in order: groups.setdefault(find(n), []).append(n)
# name each group after its longest member, ordered by size then length
ordered = sorted(groups.values(),
                 key=lambda m: (-len(m), -max(lens[x] for x in m)))
with open(os.path.join(out, "plasmid_groups.tsv"), "w") as fh:
    fh.write("genome\tcontig\tplasmid_group\tgroup_size\tlength\n")
    for i, members in enumerate(ordered, 1):
        gid = f"PG{i:03d}"
        hosts = len({m.split("|")[0] for m in members})
        for m in members:
            g, c = m.split("|", 1)
            fh.write(f"{g}\t{c}\t{gid}\t{hosts}\t{lens[m]}\n")
print(f"  {len(ordered)} plasmid groups; "
      f"{sum(1 for m in ordered if len({x.split('|')[0] for x in m})>1)} "
      f"shared across >1 host")
PY

# ---------------------------------------------------------------- 4. CheckM2
echo "[4] CheckM2"
if [ ! -s "$OUT/checkm2/quality_report.tsv" ]; then
  if command -v checkm2 >/dev/null; then
    checkm2 predict --input "$WORKG" -x fna \
      --output-directory "$OUT/checkm2" --threads "$THREADS" --force \
      > "$OUT/checkm2.log" 2>&1 || {
        echo "  checkm2 failed -- see $OUT/checkm2.log" >&2; }
  else
    echo "  checkm2 not on PATH; activate its env and rerun this step" >&2
  fi
fi
[ -s "$OUT/checkm2/quality_report.tsv" ] && \
  echo "  $(( $(wc -l < "$OUT/checkm2/quality_report.tsv") - 1 )) genomes scored"

# ---------------------------------------------------------------- 5. metadata
echo "[5] organism names and sizes"
{
  printf "genome\torganism\tsize_mb\tn_contigs\tn_plasmids\n"
  while read -r g; do
    [ -n "$g" ] || continue
    org=$(head -1 "$GDIR/${g}.fna" | sed 's/^>[^ ]* //; s/ chromosome.*//; s/,.*//' \
          | cut -c1-60)
    awk -v g="$g" -v org="$org" '
      /^>/ { n++; pl += (tolower($0) ~ /plasmid|extrachromosomal|megaplasmid/); next }
      { L += length($0) }
      END { printf "%s\t%s\t%.2f\t%d\t%d\n", g, org, L/1e6, n, pl }' "$GDIR/${g}.fna"
  done < "$LIST"
} > "$OUT/organisms.tsv"
echo "  $(( $(wc -l < "$OUT/organisms.tsv") - 1 )) rows"

echo
echo "Done. Next:"
echo "  python3 $(dirname "$0")/13_plot_community.py --data $OUT --out community_figure"
