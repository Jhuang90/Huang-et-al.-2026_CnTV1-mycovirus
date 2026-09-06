# sRNA_Viewer

Visualization of small RNA-seq alignments, emphasizing RNA size classes and strandedness.

This is a fork of [MikeAxtell/sRNA_Viewer](https://github.com/MikeAxtell/sRNA_Viewer) (original author: Michael J. Axtell, Penn State University). Version 0.4.0 fixes a critical strand-separation bug present in the original, extends annotation support, and adds several new options. See [Changes from the original](#changes-from-the-original) for a full breakdown.

---

## Installation

### Dependencies

- Linux or macOS (Windows not supported)
- [R](https://www.r-project.org) with the following packages:
  - [tidyverse](https://www.tidyverse.org)
  - [cowplot](https://cran.r-project.org/web/packages/cowplot/vignettes/introduction.html)
  - [IRanges](https://bioconductor.org/packages/release/bioc/html/IRanges.html)
  - [docopt](https://github.com/docopt/docopt.R)
- [samtools](http://www.htslib.org) >= 1.10 in your PATH
- [tabix](http://www.htslib.org/doc/tabix.html) in your PATH (required even without `-g`, checked at startup)

Install R packages:

```r
install.packages(c("tidyverse", "cowplot", "docopt"))
BiocManager::install("IRanges")
```

### Script setup

```bash
chmod +x sRNA_Viewer
# optionally move to somewhere on your PATH
mv sRNA_Viewer ~/bin/
```

---

## Usage

```
sRNA_Viewer [-g TABIXGFF -l VLINE -r RGLIST -y YLIMIT --transcript TRANSCRIPT] -c COORDINATES -b BAMLIST -p OUTPUTFILE
sRNA_Viewer -h | --help
sRNA_Viewer -v | --version
```

### Required arguments

| Flag | Description |
|------|-------------|
| `-c COORDINATES` | Genomic region in `Chr:Start-Stop` format (1-based, inclusive). Commas in numbers are accepted, e.g. `6:1,371,693-1,372,731`. Maximum interval size is 100,000 nt. |
| `-b BAMLIST` | CSV file: BAM paths in column 1, display names in column 2. Each BAM must be sorted and indexed. |
| `-p OUTPUTFILE` | Output PDF path. |

### Optional arguments

| Flag | Description |
|------|-------------|
| `-g TABIXGFF` | bgzip-compressed, tabix-indexed GFF3 file. Enables the annotation track. See [Preparing a GFF3](#preparing-a-tabix-indexed-gff3) below. |
| `-l VLINE` | Genomic coordinate at which to draw a vertical dashed line on all panels. |
| `-r RGLIST` | CSV file: read group IDs in column 1, display names in column 2. Plots each read group as a separate track from the first BAM in `-b`. |
| `-y YLIMIT` | RPM threshold for Y-axis splitting. Samples whose maximum absolute RPM exceeds this value get a free Y scale; all others are clamped to ±YLIMIT RPM. Useful when one sample has much higher coverage than others. |
| `--transcript ID` | Restrict the annotation track to a single transcript ID (e.g. `--transcript CNAG_06757.t01`). |

---

## Input file formats

### BAMLIST (`-b`)

A headerless CSV with two columns:

```
/path/to/sample1.bam,Sample 1 label
/path/to/sample2.bam,Sample 2 label
```

Every BAM must be positionally sorted and indexed (`samtools index`).

### RGLIST (`-r`)

A headerless CSV with two columns. Applies to the first BAM in BAMLIST only:

```
RG001,Condition A
RG002,Condition B
```

Read group IDs must match `@RG ID` tags in the BAM header.

---

## Preparing a tabix-indexed GFF3

Only the feature types `mRNA`, `transcript`, `ncRNA`, `lncRNA`, `pseudogenic_transcript`, `snoRNA`, `snRNA`, `tRNA`, `rRNA`, `miRNA`, `pre_miRNA`, `exon`, and `CDS` are used. Everything else is ignored.

```bash
grep -e '\smRNA\s' -e '\stranscript\s' -e '\sexon\s' -e '\sCDS\s' \
     -e '\sncRNA\s' -e '\slncRNA\s' -e '\stRNA\s' -e '\srRNA\s' \
     mygff.gff \
  | sort -k1,1 -k4,4n \
  | bgzip > mygff_sorted.gff.gz

tabix -p gff mygff_sorted.gff.gz
```

Pass the `.gz` file to `-g`. The `.tbi` index must be in the same directory.

---

## Output

A PDF where each BAM file (or read group) has its own panel. Plus-strand coverage plots above zero; minus-strand coverage plots below zero. Bars are coloured by RNA size class:

| Colour | Size class |
|--------|------------|
| Light gray | < 21 nt |
| Blue | 21 nt |
| Medium sea green | 22 nt |
| Orange | 23 nt |
| Tomato | 24 nt |
| Dark gray | > 24 nt |

All values are in RPM (reads per million), normalised against the total read count of each BAM or read group.

When `-g` is provided, an annotation track is drawn above the coverage panels showing transcript models (arrow indicating strand direction), exons (salmon), and CDS (light blue).

---

## Changes from the original

This fork diverges from [MikeAxtell/sRNA_Viewer](https://github.com/MikeAxtell/sRNA_Viewer) in the following ways.

### Bug fixes

**Strand separation (critical correctness fix)**

The original script uses `samtools depth -g 16` and `-G 16` to separate plus and minus strand reads. These flags do not reliably filter by strand in all samtools versions and contexts. This version replaces that approach by pre-filtering with `samtools view -F 16` (exclude flag 16 → plus-strand reads) and `samtools view -f 16` (require flag 16 → minus-strand reads) into separate temporary BAMs, then running `samtools depth` on each. This is the correct and reliable method.

**Annotation panel X-axis misalignment**

The original annotation panel uses `theme_void()`, which can cause the plot area to be sized differently from the coverage panels when assembled with `cowplot::plot_grid`, resulting in the annotation track not aligning with the coverage tracks below it. This version uses `theme_classic()` with a `geom_blank()` dummy layer to force the annotation panel's X axis to lock precisely to the requested coordinate range.

**Arrow overshoot on annotation track**

The original draws transcript arrows as a single `geom_segment(..., arrow = arrow())` spanning the full transcript length. When a transcript extends beyond the requested coordinate window, the arrowhead can render outside the plot area or appear distorted. This version clips arrow start/end to `userStart`/`userEnd` using `pmax`/`pmin`, and separates the backbone (a plain segment) from a short arrowhead-only segment at the 3′ end (5% of the visible span).

**CLI argument parsing for absent numeric options**

The original calls `as.numeric(opts$l)` directly, which returns `NA` (not `NULL`) when `-l` is not supplied. Downstream `!is.null(vl)` checks then incorrectly pass, causing errors. This version wraps all optional numeric arguments in a `safe_numeric()` helper that returns `NULL` for absent options.

**Misleading error messages**

The original's PATH-check error messages reference "Function sRNA_depth_by_size", which does not exist in the script. Error messages now correctly reference "sRNA_Viewer".

### New features

**`-y YLIMIT` — per-sample Y-axis split**

When one sample has much higher coverage than others, a shared Y axis either compresses low-coverage samples into near-flat lines or clips high-coverage peaks. `-y YLIMIT` splits samples into two panels: those whose maximum absolute RPM exceeds the threshold get a free Y scale; the rest are clamped to ±YLIMIT. Both panels share the same X axis and are vertically aligned.

**`--transcript ID` — single-transcript annotation filter**

Restricts the annotation track to one transcript ID. Useful for loci with many overlapping isoforms where the annotation track becomes cluttered.

**Extended GFF3 feature type support**

The original recognises only `mRNA`, `transcript`, `exon`, and `CDS`. This version additionally displays `ncRNA`, `lncRNA`, `pseudogenic_transcript`, `snoRNA`, `snRNA`, `tRNA`, `rRNA`, `miRNA`, and `pre_miRNA` as transcript-level features. This matters for any organism or GFF where non-coding genes are annotated under these types.

**Orphan transcript synthesis**

If the parent `mRNA`/`transcript` feature lies outside the query window but its child `exon` features fall within it (common when zooming into the middle of a long gene), the original plots nothing. This version detects this situation and synthesises a minimal parent row from the exon coordinate extents so the annotation still renders correctly.

**Comma-tolerant coordinate input**

The original fails to parse coordinates written with thousand-separator commas (e.g. `6:1,371,693-1,372,731`). This version strips commas before parsing, so both formats are accepted.

**Shared legend and improved multi-panel layout**

The legend is extracted from one sRNA panel, removed from all individual panels, and added back as a single shared legend to the right of the assembled figure. This prevents the legend from interfering with `plot_grid` vertical alignment across panels.

**`build_raw_cov` refactoring**

The original duplicates the entire samtools depth loop for plus and minus strands (~30 lines each). This version extracts the shared logic into an inner helper function `build_raw_cov(strand_bam)`, reducing repetition and making the strand-separation logic easier to follow.

---

## Performance note

Using multiple BAM files via `-b` is faster than splitting a single merged BAM by read groups via `-r`, because read-group filtering requires scanning the full BAM for each group.

---

## License

MIT — same as the original. See [LICENSE.txt](LICENSE.txt).

## Fork author

Jun Huang, Duke University (Jun.Huang@duke.edu)

## Original author

Michael J. Axtell, The Pennsylvania State University ([MikeAxtell/sRNA_Viewer](https://github.com/MikeAxtell/sRNA_Viewer))
