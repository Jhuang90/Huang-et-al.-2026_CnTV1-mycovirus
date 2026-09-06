#!/usr/bin/perl
use strict;
use warnings;

=head1 AUTHOR
Tim A. Dahlmann (original)
Updated by JH in 2026

=head2 CHANGE HISTORY
Updated in 2026:
  1. Counts both forward- and reverse-oriented primary mapped records.
  2. Reports length and 5-prime U/A/G/C composition for each strand and
     for both strands combined.
  3. Supports optional exclusion of one or more contigs.
  4. Skips unmapped, secondary, and supplementary alignments.
  5. Prefixes output metadata with "#" for clean R import using
     comment.char = "#".
  6. Ignores empty contig names in exclusion lists and reports records
     outside the 14-73 nt table range.
  7. Prevents the input SAM file from being used as the output path.

=head3 DESCRIPTION
Counts primary mapped SAM records after removing unmapped, secondary,
supplementary, and user-excluded contig alignments. For reverse-oriented
records (bit 0x10 set), SEQ is reverse-complemented before extracting
the 5-prime nucleotide. Output summarizes length distribution and
5-prime nucleotide composition by strand and combined across strands.

NOTE: This script counts primary mapped SAM records, not necessarily
unique reads. If the SAM file contains multiple primary records for
the same read name, that read will be counted multiple times. For
bowtie1 with -M 1 this is unlikely but should be kept in mind.

For formal library-level sRNA analysis, also generate length
distribution and 5-prime bias from trimmed FASTQ before mapping.

=head4 USAGE
# Exclude contigs on command line (comma-separated, no spaces):
perl count-read-length-from-sam.pl input.sam output.txt Mt
perl count-read-length-from-sam.pl input.sam output.txt Mt,CP045660.1

# Interactive mode (prompted at runtime):
perl count-read-length-from-sam.pl input.sam output.txt

# Import in R:
# df <- read.delim("output.txt", comment.char = "#")

=cut

# ---------------------------------------------------------------
# Subroutine: reverse complement — full IUPAC ambiguity codes
# ---------------------------------------------------------------
sub revcomp {
    my ($seq) = @_;
    $seq = reverse $seq;
    $seq =~ tr/ACGTRYKMSWBDHVNacgtrykmswbdhvn/TGCAYRMKSWVHDBNtgcayrmkswvhdbn/;
    return $seq;
}

# ---------------------------------------------------------------
# Argument check
# ---------------------------------------------------------------
my $num_args = $#ARGV + 1;
if ($num_args < 2) {
    print "\nUsage: perl count-read-length-from-sam.pl inputfile.sam outputfile.txt [contigs_to_exclude]\n";
    print "Example: perl count-read-length-from-sam.pl input.sam output.txt Mt,Pt\n\n";
    exit;
}

my $inputfile   = $ARGV[0];
my $outputfile  = $ARGV[1];
my $exclude_arg = $ARGV[2] // '';

# Safety check: prevent overwriting input file
if ($inputfile eq $outputfile) {
    die "Error: input file and output file cannot be the same path.\n";
}

open INPUT,  "< $inputfile"  or die "Can't open inputfile: $!\n";
open OUTPUT, "> $outputfile" or die "Can't open outputfile: $!\n";

print "\n===========================================\n";
print " SAM Read Length Counter - Both Strands\n";
print "===========================================\n\n";

# ---------------------------------------------------------------
# First pass: collect contig names and counts
# Only primary mapped records
# ---------------------------------------------------------------
my %contig_counts;
while (my $line = <INPUT>) {
    chomp $line;
    next if $line =~ /^@/;
    my @fields = split(/\t/, $line);
    next unless defined $fields[1] && defined $fields[2];
    my $flag = $fields[1];
    next if $flag & 4;
    next if $flag & 256;
    next if $flag & 2048;
    $contig_counts{$fields[2]}++;
}
close INPUT;

# Print contig summary to screen only
print "Contigs found in SAM file (primary mapped records only):\n";
print "------------------------------------------\n";
printf "%-30s %s\n", "Contig", "Record count";
print "------------------------------------------\n";
foreach my $contig (sort { $contig_counts{$b} <=> $contig_counts{$a} } keys %contig_counts) {
    printf "%-30s %d\n", $contig, $contig_counts{$contig};
}
print "------------------------------------------\n\n";

# ---------------------------------------------------------------
# Get exclusion list: command line OR interactive
# ---------------------------------------------------------------
my %exclude;

if ($exclude_arg =~ /\S/) {
    my @to_exclude = split(/\s*,\s*/, $exclude_arg);
    foreach my $c (@to_exclude) {
        $c =~ s/^\s+|\s+$//g;
        next if $c eq '';          # skip empty strings
        $exclude{$c} = 1;
    }
    print "Contigs to exclude (from command line): " . join(", ", sort keys %exclude) . "\n";
} else {
    print "Enter contig names to EXCLUDE, separated by commas.\n";
    print "Example: Mt,Pt,scaffold_1\n";
    print "Press ENTER with no input to include all contigs.\n\n";
    print "Contigs to exclude: ";
    my $user_input = <STDIN>;
    chomp $user_input;
    if ($user_input =~ /\S/) {
        my @to_exclude = split(/\s*,\s*/, $user_input);
        foreach my $c (@to_exclude) {
            $c =~ s/^\s+|\s+$//g;
            next if $c eq '';      # skip empty strings
            $exclude{$c} = 1;
        }
        print "\nExcluding contigs: " . join(", ", sort keys %exclude) . "\n";
    } else {
        print "\nNo contigs excluded. Analysing all records.\n";
    }
}

# Write exclusion info to output as comments
if (%exclude) {
    print OUTPUT "# Excluded contigs: " . join(", ", sort keys %exclude) . "\n";
} else {
    print OUTPUT "# No contigs excluded. Analysing all records.\n";
}

# ---------------------------------------------------------------
# Second pass: load and filter primary mapped records only
# ---------------------------------------------------------------
open INPUT, "< $inputfile" or die "Can't open inputfile: $!\n";

my @fwd_reads;
my @rev_reads;
my %excluded_counts;
my $skipped_unmapped      = 0;
my $skipped_secondary     = 0;
my $skipped_supplementary = 0;

while (my $line = <INPUT>) {
    chomp $line;
    next if $line =~ /^@/;
    my @fields = split(/\t/, $line);
    next unless defined $fields[1] && defined $fields[9];

    my $flag  = $fields[1];
    my $seq   = $fields[9];
    my $rname = $fields[2] // '*';

    if ($flag & 4)    { $skipped_unmapped++;      next; }
    if ($flag & 256)  { $skipped_secondary++;     next; }
    if ($flag & 2048) { $skipped_supplementary++; next; }

    if (exists $exclude{$rname}) {
        $excluded_counts{$rname}++;
        next;
    }

    next if $seq eq '*';

    if ($flag & 16) {
        push @rev_reads, revcomp($seq);
    } else {
        push @fwd_reads, $seq;
    }
}
close INPUT;

my $n_fwd   = scalar(@fwd_reads);
my $n_rev   = scalar(@rev_reads);
my $n_total = $n_fwd + $n_rev;

# Print filtering summary to screen
print "\n--- Record filtering summary ---\n";
printf "  Unmapped records skipped:          %d\n", $skipped_unmapped;
printf "  Secondary alignments skipped:      %d\n", $skipped_secondary;
printf "  Supplementary alignments skipped:  %d\n", $skipped_supplementary;
if (%excluded_counts) {
    foreach my $contig (sort keys %excluded_counts) {
        printf "  %-30s %d records excluded\n", $contig, $excluded_counts{$contig};
    }
}
print "\nForward-oriented records (0x10 not set): $n_fwd\n";
print "Reverse-oriented records (0x10 set):     $n_rev\n";
print "Total primary mapped records:            $n_total\n";

# Write filtering summary to output as comments
print OUTPUT "# --- Record filtering summary ---\n";
printf OUTPUT "# Unmapped records skipped:\t%d\n",          $skipped_unmapped;
printf OUTPUT "# Secondary alignments skipped:\t%d\n",      $skipped_secondary;
printf OUTPUT "# Supplementary alignments skipped:\t%d\n",  $skipped_supplementary;
if (%excluded_counts) {
    foreach my $contig (sort keys %excluded_counts) {
        printf OUTPUT "# Excluded %s:\t%d records\n", $contig, $excluded_counts{$contig};
    }
}
printf OUTPUT "# Forward-oriented records (0x10 not set):\t%d\n", $n_fwd;
printf OUTPUT "# Reverse-oriented records (0x10 set):\t%d\n",     $n_rev;
printf OUTPUT "# Total primary mapped records:\t%d\n",             $n_total;

# ---------------------------------------------------------------
# Count lengths and 5' nt
# ---------------------------------------------------------------
my %counts;
for my $len (14..73) {
    for my $strand ('all', 'fwd', 'rev') {
        for my $nt ('total', 'U', 'A', 'G', 'C') {
            $counts{$len}{$strand}{$nt} = 0;
        }
    }
}

my $out_of_range = 0;

sub tally {
    my ($seq, $strand) = @_;
    my $len = length($seq);
    if ($len < 14 || $len > 73) {
        $out_of_range++;
        return;
    }
    $counts{$len}{$strand}{'total'}++;
    $counts{$len}{'all'}{'total'}++;
    my $nt = uc(substr($seq, 0, 1));
    $nt = 'U' if $nt eq 'T';
    if ($nt eq 'U' || $nt eq 'A' || $nt eq 'G' || $nt eq 'C') {
        $counts{$len}{$strand}{$nt}++;
        $counts{$len}{'all'}{$nt}++;
    }
}

foreach my $seq (@fwd_reads) { tally($seq, 'fwd'); }
foreach my $seq (@rev_reads) { tally($seq, 'rev'); }

print        "\nRecords outside 14-73 nt ignored in table: $out_of_range\n\n";
printf OUTPUT "# Records outside 14-73 nt ignored in table:\t%d\n", $out_of_range;

# ---------------------------------------------------------------
# Output: wide table — R-readable, comments prefixed with #
# ---------------------------------------------------------------
print OUTPUT "length\tall_total\tall_U\tall_A\tall_G\tall_C\tfwd_total\tfwd_U\tfwd_A\tfwd_G\tfwd_C\trev_total\trev_U\trev_A\trev_G\trev_C\n";

for my $len (14..73) {
    printf OUTPUT "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
        $len,
        $counts{$len}{'all'}{'total'},
        $counts{$len}{'all'}{'U'},
        $counts{$len}{'all'}{'A'},
        $counts{$len}{'all'}{'G'},
        $counts{$len}{'all'}{'C'},
        $counts{$len}{'fwd'}{'total'},
        $counts{$len}{'fwd'}{'U'},
        $counts{$len}{'fwd'}{'A'},
        $counts{$len}{'fwd'}{'G'},
        $counts{$len}{'fwd'}{'C'},
        $counts{$len}{'rev'}{'total'},
        $counts{$len}{'rev'}{'U'},
        $counts{$len}{'rev'}{'A'},
        $counts{$len}{'rev'}{'G'},
        $counts{$len}{'rev'}{'C'};
}

# Screen summary
print "length\tall_total\tfwd_total\trev_total\n";
for my $len (14..73) {
    next if $counts{$len}{'all'}{'total'} == 0;
    printf "%d\t%d\t%d\t%d\n",
        $len,
        $counts{$len}{'all'}{'total'},
        $counts{$len}{'fwd'}{'total'},
        $counts{$len}{'rev'}{'total'};
}

print "\nDone! Results written to $outputfile\n";
close OUTPUT;
