# Huang-etal-2026-CnTV1-mycovirus
Data, analyses, and resources for Huang et al. (2026) on the CnTV1 mycovirus.

## Preprint and citation

This repository accompanies our [bioRxiv preprint](https://www.biorxiv.org/content/10.64898/2026.09.15.751489v1).

**If you use the data, code, analyses, or other resources in this repository in your research, please cite our preprint:**

> Huang, J., Larmore, C. J., Davenport, T. C., Averette, A. F., Choi, Y., Debat, H., Gupta, P., Babaian, A., Meneghini, M. D., Sun, S., & Heitman, J. (2026). **An ancestral mitochondrial DNA insertion disrupts RNAi and enables persistence of a novel mycovirus in *Cryptococcus neoformans*.** bioRxiv [Preprint]. https://doi.org/10.64898/2026.09.15.751489

## Analysis tools

- [sRNA_Viewer](sRNA_Viewer/README.md): small RNA-seq coverage visualization by RNA length and strand, with optional GFF3 annotation tracks. Includes the script, documentation, MIT license, and example outputs.
- [smallRNA SAM read-length counter](smallRNA/README.md): Perl utility for counting primary mapped SAM records by read length and biological 5′ nucleotide, with forward/reverse summaries, alignment filtering, and optional contig exclusion.
- [RIBBON RUNNER](ribbonrunner/README.md): automated whole-genome synteny pipeline that runs adjacent-pair MUMmer alignments and produces publication-quality SVG ribbon figures.
- [mosdepth coverage notebook](mosdepth/README.md): paired-end DNA read mapping, 500 bp window depth, and genome-wide coverage plots normalized to each sample's median window depth.
- [Sources and versions](SOURCE.md): origins, pinned commits, authorship, and licenses for all included analysis tools.
