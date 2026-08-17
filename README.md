# hyscale-sim

Simulation and benchmarking framework for barcode-guided linkage of
single-cell short reads to bulk PacBio HiFi long reads.

This repository contains the in silico development pipeline for Task 5 of the
Hy-SCALE project: assigning plasmids to their host cells by propagating
cell-specific barcodes from a short-read stream onto raw long reads, and
measuring how accurately that can be done when ground truth is known exactly.

\---

## What problem this addresses

Metagenomic binning assigns contigs to genomes using sequence composition and
coverage depth. Both properties differ systematically between a plasmid and its
host chromosome, so plasmid–host linkage is exactly the case where binning is
least reliable. Droplet or hydrogel single-cell methods supply a physical
constraint instead — two sequences in the same compartment came from the same
cell — but the barcode is present only on the short reads, while the long reads
that span mobile elements carry none.

This framework simulates both streams from known reference genomes, propagates
the barcodes, and scores the result against ground truth. Because every read's
origin is known, false discovery is measured rather than estimated.

\---

## Findings

These emerged from running the framework and are documented here because
several are non-obvious and cost significant compute to establish.

**Read accuracy is the binding constraint, not read length.** An initial run
used PacBio CLR reads (88.5% identity, Q10). Short reads could not find exact
anchors, so the mean number of supporting short reads per long read stayed near
1.0 regardless of coverage — flat across a hundred-fold change. Regenerating as
HiFi (99.5% identity via multi-pass CCS) restored the expected scaling. Any
implementation of this approach requires HiFi-grade long reads; CLR does not
work at all.

**Reference redundancy defeats per-read mapping.** Mapping short reads to raw
long reads at 13× depth means every locus appears in \~13 targets, so MAPQ is
uninformative — it measures reference redundancy, not strain ambiguity. Only
1,577 of 3.6M alignments had MAPQ ≥ 1. Aggregate assignment across a barcode's
whole read set works; per-read assignment does not.

**Community composition determines feasibility.** With 473 congeneric genomes
at 96–99% ANI, a core-genome short read has \~6,000 valid targets and the
`-N` cap admits a near-random sample of them; own-genome hits are crowded out
and linkage fails entirely. Restricting to a realistically structured community
— 8 co-resident *E. coli* at 98–99.98% ANI within 28 genomes across 16 genera —
reduces this to \~104 targets and linkage succeeds. The failure is a property of
the benchmark, not of the method, but it defines the density of near-identical
strains the approach tolerates.

**Coverage has an optimum.** Assignment completeness and accuracy respond to
per-cell coverage in opposite directions. Below λ = c·L/l ≈ 10 too few plasmid
fragments recruit supporting reads; above λ ≈ 40 completeness saturates while
accuracy declines, because additional reads recruit long reads from other
strains sharing the same sequence. Overall association precision falls
monotonically from 99.6% at 0.01× to 86.0% at 1×.

**The residual error is structural.** Varying the significance threshold over
two orders of magnitude and doubling the minimum supporting-read requirement
changed accuracy by under 0.4 percentage points. The misassignments are
well-supported associations to genuinely shared sequence, not marginal calls.

**Error scales with genome similarity.** Per-genome false discovery correlates
with nearest-neighbour ANI at every coverage level (Spearman ρ = +0.46 to
+0.50, all p < 5 × 10⁻³, n = 36).

**Artifact detection needs normalization against strain diffusion.** A raw
minor-strain fraction does not separate multi-cell compartments from
single-cell ones in a community containing near-identical strains, because
profiles are diffused both by a genuine second cell and by mapping spread. The
raw statistic ranked near-identical doublets *above* distant ones — backwards
for a mixing statistic, and a sign it was partly classifying taxonomy. Scoring
each barcode against the per-strain median profile cancels the diffusion term
and restores the expected monotone decline with ANI.

\---

## Pipeline

```
09  simulate HiFi reads        pbsim3 multi-pass -> CCS -> traceable headers
01  build truth tables         barcode/long-read/contig maps, strain groups
05  strain groups              collapse indistinguishable assemblies
10  select subset              choose a realistic community by ANI structure
11  build subset               filter reads and truth to that community
06  inject negatives           doublets and ambient DNA at controlled rates
02  link                       minimap2 -> reduce\_links.py -> incidence table
03  score                      hypergeometric test + Benjamini-Hochberg
07  long-read contigs          resolve reads to plasmid vs chromosome
08  plasmid linkage            host assignment accuracy (the M5.1 metric)
04  analyse                    precision/recall, error vs similarity
16  detect artifacts           profile-normalized doublet and ambient detection
12  characterise               ANI, mash, plasmid clusters, CheckM2
13  plot community             composite community figure
14  coverage series            driver for the full coverage ladder
```

Scripts are numbered by the order they were developed, not strictly by
execution order; see `docs/PIPELINE\_DESIGN.md` for the intended sequence.

`15\_detect\_artifacts.py` is superseded by `16\_detect\_artifacts\_normalized.py`
and is retained only because the comparison between them documents the
diffusion confound described above.

\---

## Requirements

```bash
mamba env create -f environment.yml
conda activate hyscale-sim
```

External tools: `pbsim3`, `pbccs`, `samtools`, `minimap2` (≥2.17), `seqkit`,
`skani`, `mash`, `checkm2`, `parallel`. Python: `numpy`, `scipy`,
`matplotlib`.

Hardware used for the reference run: 80 cores, 754 GB RAM. The HiFi simulation
across 473 genomes takes roughly 9 h at `JOBS=8`; the linking stage is the
other significant cost and scales with per-cell coverage.

\---

Short-read simulation. 

Barcode-linked short reads are generated by

scripts/00\_simulate\_short\_reads.sh using ART. Each simulated cell is an

independent ART run against its source genome at the target per-cell

coverage, so a barcode's reads sample that genome independently of every

other cell — the property the linkage analysis depends on.



```

art\_illumina -ss HS25 -p -na -l 150 -f <per-cell coverage> -m 350 -s 50

```

HiSeq 2500 error profile, 150 bp paired-end, 350 bp mean fragment (SD 50),

30 cells per genome, per-cell coverages 0.01, 0.1, 0.2, 0.4, 0.8, 1 and 10×.

Barcode identifiers are assigned sequentially across the whole run

(BC000001 onward), so a barcode ID is globally unique but does not encode

its coverage tier or source genome; that mapping is recovered from the

directory structure by 01\_build\_truth.sh.





## Reproducing

Simulation inputs (reference genomes) and outputs (reads) are far too large for
version control — the HiFi set alone is 28 GB and the short-read set 500 GB.
The genomes are public RefSeq/GenBank assemblies; accessions are listed in
`docs/genome\_accessions.txt`. Simulated reads can be regenerated from those
accessions with `scripts/09\_simulate\_hifi\_all.sh` and the short-read simulator
parameters recorded in `docs/PIPELINE\_DESIGN.md`.

```bash
# 1. simulate HiFi long reads from reference genomes
TMPDIR=/path/to/scratch JOBS=8 bash scripts/09\_simulate\_hifi\_all.sh

# 2. build truth tables and strain groups
bash scripts/01\_build\_truth.sh

# 3. choose and build a benchmark community
python3 scripts/10\_select\_subset.py --truth truth \\
  --out truth/subset\_genomes.txt --n-close 5 --n-diverse 25
bash scripts/11\_build\_subset.sh truth/subset\_genomes.txt subset\_v1

# 4. run the coverage series
cd subset\_v1 \&\& export ROOT=$PWD
bash ../scripts/14\_run\_coverage\_series.sh 0.01x 0.1x 0.2x 0.4x 0.8x 1x

# 5. artifact detection
python3 ../scripts/16\_detect\_artifacts\_normalized.py \\
  --tiers 0.1x 0.2x 0.4x 0.8x --linkdir links --truth truth \\
  --out results\_detect\_norm
```

\---

## Status and limitations

This is research code from an active project, not a packaged tool. Specifically:

* **Only the association-filtering layer is implemented.** Barcode Jaccard
clustering into strain bins and allelic phasing concordance are designed but
not built. Reported accuracy figures are for association filtering alone.
* **Paths and environment assumptions are partly hard-coded** to the
development system. Most are overridable by environment variable; some are
not.
* **`03\_score.py` holds all tested pairs in memory** and needs roughly 10 GB at
30M pairs. It would need chunking for substantially larger runs.
* **The short-read simulation is not included here.** Barcode-linked short
reads were generated with ART in a separate step; only the parameters are
documented.
* **No test suite.** Individual components were validated against synthetic
fixtures during development, but those checks are not packaged as tests.

Issues and questions are welcome, but responses may be slow while the project
is in progress.

\---

## Citation

If this framework is useful in your work, please cite the associated
publication (in preparation) and link to this repository.

## License

MIT — see `LICENSE`.

