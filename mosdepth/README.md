# Genome-wide coverage plots with mosdepth

[Open the notebook](mosdepth_coverage_plot.ipynb).

This notebook maps paired-end whole-genome DNA reads with minimap2, filters alignments at MAPQ >= 30, computes mean depth in 500 bp windows with mosdepth, and plots each sample relative to its median window depth across the retained chromosomes.

## Requirements

- Command line: `minimap2`, `samtools`, `mosdepth`, and Bash.
- R: `dplyr`, `readr`, `tibble`, and `ggplot2` >= 3.4.0.
- Notebook execution: Jupyter with a Python 3 kernel and `rpy2`, with R installed and available to `rpy2`.
- Optional rasterized points: `ggrastr` and `ragg`.

The R dependencies can be installed with:

```r
install.packages(c("dplyr", "readr", "tibble", "ggplot2"))
```

For notebook execution, install Jupyter and `rpy2` in the Python environment used by the notebook kernel:

```bash
python -m pip install jupyterlab rpy2
```

Alternatively, run the Bash cells in a terminal without the `%%bash` line, and run the R plotting cell in RStudio without its `%%R` line. This route does not require Jupyter or `rpy2`.

## Usage

1. Work in a directory containing a reference FASTA and paired reads named `Sample1_S1_R1_001.fastq.gz` and `Sample1_S1_R2_001.fastq.gz` (and equivalent pairs for other samples).
2. Set `REF` in Step 1 to your reference filename. Run Step 1 to produce sorted, indexed BAM files and a reference `.fai` index.
3. Run Step 2 to produce `<sample>.regions.bed.gz` files. Unused per-base coverage output is disabled; standard mosdepth overlapping-mate handling is retained.
4. Edit `fai_path`, `samples`, `colors`, `y_max`, and `output_file` in Step 3. Run the R cell to print a coverage summary and save `Mosdepth_normalized_coverage.pdf`.

All samples must use the same reference and window scheme. If `.fai` and mosdepth window files already exist, you can run only the plotting step. Rerunning a step overwrites outputs with the same names.

The example names are placeholders, not bundled study data. No sequencing data or executed notebook outputs are included.

## Interpretation

Normalized coverage is relative to each sample's median window depth, not an absolute copy-number measurement. A value of 2 indicates twice the sample baseline; interpreting it as disomy assumes a predominantly haploid baseline. Zero coverage can reflect a deletion, mapping filters, or poor mappability. Widespread copy-number changes can shift the median. Excluding a contig from the chromosome layout also excludes it from normalization. Each window contributes equally to the median, including shorter terminal windows.

## Validation

The R plotting code was exercised with synthetic coverage files using R 4.6.0, dplyr 1.2.1, readr 2.2.0, ggplot2 4.0.3, and tibble 3.3.1. Checks covered relative coverage of 0/1/2, chromosome ordering, values above the plot ceiling, PDF generation, and rejection of zero-median, negative, non-finite, out-of-bounds, or missing-chromosome inputs. The generated plot was visually inspected.

Both Bash cells passed syntax checks; stub-command tests verified missing-input errors and upstream failure propagation. The full minimap2/samtools/mosdepth workflow and the Jupyter/rpy2 integration were not executed in this validation environment.

## Citation and tool documentation

If you use this notebook or its analyses in your research, please cite [Huang et al. (2026), bioRxiv](https://www.biorxiv.org/content/10.64898/2026.09.15.751489v1), DOI: [10.64898/2026.09.15.751489](https://doi.org/10.64898/2026.09.15.751489). The full reference is in the [repository README](../README.md#preprint-and-citation).

Also follow the citation guidance for the tools used: [minimap2](https://github.com/lh3/minimap2), [samtools](https://www.htslib.org/), and [mosdepth](https://github.com/brentp/mosdepth).
