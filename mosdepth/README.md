# Genome-wide coverage plots with mosdepth

[Open the notebook](mosdepth_coverage_plot.ipynb).

This notebook maps paired-end whole-genome DNA reads with minimap2, filters alignments at MAPQ >= 30, computes mean depth in 500 bp windows with mosdepth, and plots coverage normalized to each sample's median window depth.

The code cells are unchanged from Jun Huang's original notebook. Citation information and a brief clarification of coverage interpretation are included in the notebook text.

## Requirements

- Command line: `minimap2`, `samtools`, `mosdepth`, and Bash.
- R: `dplyr`, `readr`, and `ggplot2` >= 3.4.0 (with its `tibble` dependency).
- To run the notebook: Jupyter with a Python 3 kernel, `rpy2`, and R.
- Optional rasterized points: `ggrastr` and `ragg`.

Alternatively, run the Bash cells in a terminal without the `%%bash` line and the R plotting cell in RStudio without the `%%R` line.

## Usage

1. Place the reference FASTA and paired reads (for example, `Sample1_S1_R1_001.fastq.gz` / `Sample1_S1_R2_001.fastq.gz`) in your working directory.
2. Set `REF` in Step 1 and run the mapping, filtering, sorting, and indexing commands.
3. Run Step 2 to generate mosdepth window coverage files.
4. Edit `fai_path`, `samples`, `colors`, `y_max`, and `output_file` in Step 3, then run the R cell to generate the coverage summary and PDF.

If the reference `.fai` and `<sample>.regions.bed.gz` files already exist, you can run just the plotting step. Sample names in the notebook are placeholders.

Coverage values are relative to each sample's median window depth. They are not absolute ploidy measurements, and zero coverage alone does not establish a deletion.

## Citation

If you use this notebook or its analyses in your research, please cite [Huang et al. (2026), bioRxiv](https://www.biorxiv.org/content/10.64898/2026.09.15.751489v1), DOI: [10.64898/2026.09.15.751489](https://doi.org/10.64898/2026.09.15.751489). The full reference is in the [repository README](../README.md#preprint-and-citation).

Tool documentation and citation guidance: [minimap2](https://github.com/lh3/minimap2), [samtools](https://www.htslib.org/), and [mosdepth](https://github.com/brentp/mosdepth).
