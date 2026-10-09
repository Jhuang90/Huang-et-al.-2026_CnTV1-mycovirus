# NRHc5010 mitochondrial annotation and Circos NUMT plot

This directory contains the data and configuration needed to reproduce the NRHc5010 mitochondrial-to-nuclear Circos plot.

## Files

- `NRHc5010_mt_vs_core.blast`: original mitochondrial-versus-nuclear BLAST results (23 records, standard 12-column tabular format).
- `NRHc5010_MT.annotated`: MFannot mitochondrial annotation, including the 27,977 bp mitochondrial sequence.
- `NRHc5010_mt_vs_core.filtered.blast`: 15 alignments with identity ≥90% and alignment length ≥100 bp.
- `karyotype.txt`: chromosome lengths, order, and colors.
- `NRHc5010_links2.txt`: 15 colored mitochondrial-to-nuclear links.
- `NRHc5010_highlight_v2.txt`: 39 mitochondrial feature highlights.
- `NRHc5010_mt_annot_gene_only_v2.txt`: 18 mitochondrial feature labels derived from the MFannot annotation.
- `NRHc5010_circos_v2.conf`: Circos V2 configuration.
- `5010_mtcircos_V2.png` and `5010_mtcircos_V2.svg`: regenerated plot output.

## Run

The configuration was tested with Circos 0.69-8. Use the Perl environment associated with your Circos installation, including its required modules, configuration files, and fonts.

In `NRHc5010_circos_v2.conf`, replace `/path/to/circos/etc/` with your Circos configuration directory and `/path/to/data/` with the absolute path to this directory. Then run:

```bash
mkdir -p output
/path/to/circos -conf NRHc5010_circos_v2.conf
```

The command creates `output/5010_mtcircos_V2.png` and `output/5010_mtcircos_V2.svg`.

The included PNG and SVG were generated with the supplied configuration and plotting inputs.

These files reproduce the Circos plot from the existing BLAST results. Rerunning BLAST additionally requires the mitochondrial and nuclear genome FASTA files and BLAST+.

## Citation

**If you use the code, data, or plotting files from this repository in your research, please cite the following preprint:**

> Huang, J., Larmore, C. J., Davenport, T. C., Averette, A. F., Choi, Y., Debat, H., Gupta, P., Babaian, A., Meneghini, M. D., Sun, S., & Heitman, J. (2026). **An ancestral mitochondrial DNA insertion disrupts RNAi and enables persistence of a novel mycovirus in *Cryptococcus neoformans*.** bioRxiv [Preprint]. https://doi.org/10.64898/2026.09.15.751489

[View the bioRxiv preprint](https://www.biorxiv.org/content/10.64898/2026.09.15.751489v1).
