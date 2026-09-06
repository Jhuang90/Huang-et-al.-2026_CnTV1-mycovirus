# Changelog

All notable changes to RIBBON RUNNER. Versions are git tags; there is only
ever one `synteny_pipeline.py` and one `synteny_plot.py` on `main`.

## v1.3.0

### Fixed

- **`synteny_pipeline.py` on `main` was still the v1.0.0 code.** The README
  told users to run it, but it ignored `--min_len` / `--min_idy` at the
  filtering step, hardcoded `conda activate syri2`, and created the SLURM
  log directory too late for `#SBATCH --output` to resolve. `main` now
  carries the current code, and the older `*_v1.x.y.py` copies are gone —
  version history lives in git tags.
- **All rows shared one width instead of one scale.** `compute_layout`
  recomputed `scale` per strain from that strain's own total length, so
  every row was stretched to the full drawing width and genome-size
  differences became invisible. The scale bar was drawn from the first
  strain's scale only, so it was wrong for every row below it. There is now
  a single scale derived from the longest genome; shorter assemblies end
  short of the right margin, which is the point.
- **"All chromosomes" only ever meant Chr1–Chr14.** `DEFAULT_CHROMS` was
  hardcoded, so `synteny_all_chr.svg` silently dropped MT and anything past
  Chr14 for other species. The chromosome list is now read from the `.fai`
  files (union across strains, natural-sorted, organelles last), with
  `--min_chrom_len` for draft assemblies full of small contigs.
- **Malformed `show-coords` output was silently skipped**, producing an
  empty figure with no error. Unparseable lines are now counted and
  reported, zero parsed alignments is a hard error, and zero drawn ribbons
  prints a warning that names the likely cause.
- **Nothing was shell-quoted.** Every path and strain name is now passed
  through `shlex.quote()`, so paths containing spaces work and a strain
  name can no longer inject commands into the generated script.
- Merged blocks reported an unweighted running mean identity; it is now
  length-weighted.
- `numpy` was imported but never used, and is no longer a dependency.

### Added

- `examples/make_example_data.py` — a seeded 3-strain synthetic dataset
  (~375 kb per genome, deliberate inversion, translocation and deletion,
  three different header styles). The whole pipeline runs on it in about
  two seconds, which makes it usable as a smoke test.
- `environment.yml`, `.gitignore`, this changelog.
- `--version` on both scripts.
- `--threads`, passed to nucmer only when non-zero (MUMmer 4 only) and
  mirrored into `#SBATCH --cpus-per-task`. Previously 4 CPUs were requested
  and 3 sat idle.
- Validation before any work starts: duplicate strain names, unreadable
  FASTA paths, fewer than two strains, and two sequences that normalise to
  the same chromosome name.

### Changed

- **samtools and seqkit are no longer required.** Renaming, rewrapping and
  `.fai` indexing happen in Python, in one place, instead of in a bash
  heredoc that duplicated the same logic. Dependencies are now Python,
  matplotlib and MUMmer.
- Inversions are drawn as twisted ribbons rather than parallel bands that
  only differ in colour.
- `--conda_init` is auto-detected from the environment instead of
  defaulting to one specific Duke HPC path.
- Default `--time` is `24:00:00`; `240:00:00` is rejected by most cluster
  QOS policies.
- Generated scripts use `set -euo pipefail`.
- The strain colour palette holds 12 entries instead of 8.

## v1.2.1

- Fixed cosmetic formatting in the generated SLURM script where the conda
  activation line and `set -e` were not cleanly separated.

## v1.2.0

- Removed dead functions (`detect_chr_format`, `rename_fasta`,
  `reformat_fasta`, `index_cmd`) that duplicated the inline heredoc logic.
- Fixed SLURM log directory race so `#SBATCH --output` resolves on a fresh
  output directory.
- Switched title quoting from double to single quotes.
- `parse_strains` requires tab-separated input, with the line number in the
  error message.
- Added `--conda_init` and `--conda_env`.
- Print the alignment order at runtime.

## v1.1.0

- `filter_cmd` passes `--min_idy` and `--min_len` through to `delta-filter`;
  previously filtering was hardcoded at `-i 90 -l 10000`.
- Added contact email to the header banner.

## v1.0.0

- Initial pipeline: FASTA input → chromosome renaming → nucmer alignment →
  show-coords → SVG ribbon figure, with SLURM job generation.
