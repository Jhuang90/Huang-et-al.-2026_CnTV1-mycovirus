# SAM read-length and 5′-nucleotide counter

Perl utility used for small-RNA sequencing analysis in the CnTV1 project.

The original `count-read-length-from-sam.pl` was written by Dr. Tim A. Dahlmann and distributed in [timdahlmann/smallRNA](https://github.com/timdahlmann/smallRNA). This updated version retains the original 14–73 nt reporting range while adding strand-aware SAM parsing, filtering, richer output, and safer operation.

## Updated SAM read-length counter

`count-read-length-from-sam.pl` counts primary mapped SAM records by read length and 5′ nucleotide. Results are reported for forward-oriented records, reverse-oriented records, and both strands combined.

### Requirements

- Perl 5 (no non-core modules required)
- A SAM-format alignment file

### Usage

Exclude one contig:

```bash
perl smallRNA/count-read-length-from-sam.pl input.sam output.txt Mt
```

Exclude multiple contigs with a comma-separated list:

```bash
perl smallRNA/count-read-length-from-sam.pl input.sam output.txt Mt,CP045660.1
```

Run without the optional third argument to choose exclusions interactively. Press Enter at the prompt to include all contigs:

```bash
perl smallRNA/count-read-length-from-sam.pl input.sam output.txt
```

### What the updated script does

- Parses standard SAM fields and bitwise flags.
- Counts primary mapped records only.
- Skips unmapped (`0x4`), secondary (`0x100`), and supplementary (`0x800`) records.
- Allows selected contigs to be excluded and reports how many records were removed from each.
- Separates forward- and reverse-oriented records using flag `0x10`.
- Reverse-complements reverse-oriented sequences, including IUPAC ambiguity codes, before determining their biological 5′ nucleotide.
- Reports counts for lengths 14–73 nt and records how many reads fall outside that range.
- Produces one tab-delimited table containing total, 5′ U, 5′ A, 5′ G, and 5′ C counts for all, forward, and reverse records.
- Prefixes metadata and filtering summaries with `#`, making the output easy to import into R.
- Refuses to use the same path for input and output, protecting the SAM file from accidental overwrite.

The output columns are:

```text
length  all_total  all_U  all_A  all_G  all_C  fwd_total  fwd_U  fwd_A  fwd_G  fwd_C  rev_total  rev_U  rev_A  rev_G  rev_C
```

Import the table into R while ignoring the metadata lines:

```r
df <- read.delim("output.txt", comment.char = "#")
```

> The script counts primary mapped SAM records, not deduplicated read names. If a SAM file contains multiple primary alignments for the same read name, each record is counted.

## Improvements over the original version

| Area | Original version | Updated version |
| --- | --- | --- |
| SAM parsing | Regular-expression matching designed around a limited record pattern | Field-based parsing of standard SAM records and bitwise flags |
| Strand support | Counted only the original matched orientation | Counts forward and reverse orientations separately and together |
| Reverse-strand 5′ base | Not handled | Reverse-complements the sequence before extracting the 5′ base |
| Alignment filtering | No explicit primary/unmapped filtering | Removes unmapped, secondary, and supplementary records |
| Contig filtering | Not available | Optional command-line or interactive exclusion list |
| Output structure | Separate unheaded blocks for total and each 5′ nucleotide | One labeled, R-readable wide table |
| Nucleotide summary | Separate U/A/G/C sections | U/A/G/C counts for all, forward, and reverse records |
| Range reporting | Values outside 14–73 nt were not clearly summarized | Explicit out-of-range counter |
| Safety | Output could overwrite the input path | Input/output same-path check |
| Audit information | Minimal run summary | Filtering and record totals saved as `#` metadata lines |

## License

This derivative remains available under the [GNU General Public License v3.0](https://github.com/timdahlmann/smallRNA/blob/master/LICENSE).
