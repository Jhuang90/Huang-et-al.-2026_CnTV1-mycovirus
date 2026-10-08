# Sources and versions

This repository contains self-contained copies of the analysis tools listed below. Imported repository snapshots are pinned to the commits listed here and are not automatically synchronized with their standalone repositories. Locally contributed tools are identified separately below.

## sRNA_Viewer

- Included directory: [`sRNA_Viewer/`](sRNA_Viewer/)
- Standalone repository: [Jhuang90/sRNA_Viewer](https://github.com/Jhuang90/sRNA_Viewer), a fork of [MikeAxtell/sRNA_Viewer](https://github.com/MikeAxtell/sRNA_Viewer)
- Version: 0.4.0
- Source commit: [`7953c9e1e96e8d42b740350039634e2c71246b2b`](https://github.com/Jhuang90/sRNA_Viewer/commit/7953c9e1e96e8d42b740350039634e2c71246b2b)
- Original author: Michael J. Axtell
- Fork maintainer: Jun Huang
- License: [MIT](sRNA_Viewer/LICENSE.txt)

The six source files were copied unchanged from the pinned commit. `Example.pdf`, `Example.png`, and `README.html` are inherited source-repository assets, not newly generated results for the CnTV1 study.

## smallRNA SAM read-length counter

- Included script: [`smallRNA/count-read-length-from-sam.pl`](smallRNA/count-read-length-from-sam.pl)
- Standalone repository: [Jhuang90/smallRNA](https://github.com/Jhuang90/smallRNA), a fork of [timdahlmann/smallRNA](https://github.com/timdahlmann/smallRNA)
- Source commit: [`ddbadd0e0cce17c3c898a65f165d16790419e1d0`](https://github.com/Jhuang90/smallRNA/commit/ddbadd0e0cce17c3c898a65f165d16790419e1d0)
- Original author: Tim A. Dahlmann
- Updated by: Jun Huang
- License: [GNU General Public License v3.0](https://github.com/Jhuang90/smallRNA/blob/main/LICENSE)

Only the updated script and its project-specific [`smallRNA/README.md`](smallRNA/README.md) are included here. The filename does not carry an internal version suffix.

## RIBBON RUNNER

- Included directory: [`ribbonrunner/`](ribbonrunner/)
- Standalone repository: [Jhuang90/ribbonrunner](https://github.com/Jhuang90/ribbonrunner)
- Version: 1.3.0
- Source commit: [`3e74f37b3a665293fabc74c2a0c9896243fecd76`](https://github.com/Jhuang90/ribbonrunner/commit/3e74f37b3a665293fabc74c2a0c9896243fecd76)
- Author and maintainer: Jun Huang
- License: [MIT](ribbonrunner/LICENSE)

The included snapshot contains the two pipeline scripts, Conda environment, documentation, changelog, license, and bundled synthetic example. Generated alignment intermediates are not included.

## mosdepth coverage notebook

- Included directory: [`mosdepth/`](mosdepth/)
- Source: `mosdepth_coverage_plot.ipynb`, supplied by Jun Huang on 2026-10-08.
- Original file SHA-256: `2b61d54e348dd3d6e151b5e3ed61d9c3b845b895b614c45a16266ef3017c9972`
- Repository copy: all code cells match the supplied original notebook exactly. Only citation text and relative-coverage interpretation wording were added or clarified.
- Dependencies, usage, and validation scope: [`mosdepth/README.md`](mosdepth/README.md).

The original notebook was preserved outside this repository. This contribution contains placeholder sample names and no sequencing data or executed outputs.
