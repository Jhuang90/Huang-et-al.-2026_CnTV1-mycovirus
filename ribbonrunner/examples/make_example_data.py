#!/usr/bin/env python3
"""
make_example_data.py — generate a tiny synthetic 3-strain dataset for
testing and demonstrating RIBBON RUNNER.

Produces three small "genomes" (~375 kb each) with known, deliberate
rearrangements, so the resulting figure must show all three ribbon types:

    strainA   reference. Headers use strain-prefixed style:  A_Chr1 ...
    strainB   inversion on Chr2; 30 kb of Chr3 translocated onto Chr1.
              Headers use lowercase style:                   chr1 ...
    strainC   further inversion on Chr1; 20 kb deletion on Chr3.
              Headers use numeric (H99) style:               1, 2, 3, Mt

The three header styles are intentional: they exercise the chromosome
name normalisation code. The Chr3 deletion in strainC is also
intentional: on a correctly scaled figure strainC's Chr3 must be visibly
shorter than the others.

Everything is seeded, so the data is byte-identical on every machine.

    python3 make_example_data.py --outdir .
"""

import argparse
import os
import random

LINE_WIDTH = 60

# name -> length (bp)
CHROM_SIZES = {
    "Chr1": 150_000,
    "Chr2": 120_000,
    "Chr3": 80_000,
    "MT":    25_000,
}

# How each strain spells its sequence names in the FASTA header.
HEADER_STYLE = {
    "strainA": lambda c: f"A_{c}",                              # strain-prefixed
    "strainB": lambda c: c.lower() if c != "MT" else "Mt",      # lowercase
    "strainC": lambda c: c.replace("Chr", "") if c != "MT" else "Mt",  # numeric
}


def random_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def revcomp(s):
    return s.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def mutate(rng, seq, rate):
    """Introduce point substitutions at the given per-base rate."""
    out = list(seq)
    for i in range(len(out)):
        if rng.random() < rate:
            out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
    return "".join(out)


def write_fasta(path, records, header_fn):
    with open(path, "w") as f:
        for name, seq in records:
            f.write(f">{header_fn(name)}\n")
            for i in range(0, len(seq), LINE_WIDTH):
                f.write(seq[i:i + LINE_WIDTH] + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default=".", help="Where to write the FASTA files")
    p.add_argument("--seed", type=int, default=20260830)
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    rng = random.Random(args.seed)

    # ---- strainA: the reference ----
    A = {c: random_seq(rng, n) for c, n in CHROM_SIZES.items()}

    # ---- strainB: inversion on Chr2 + translocation Chr3 -> Chr1 ----
    B = {c: mutate(rng, s, 0.004) for c, s in A.items()}

    # Chr2: invert 40 kb (40,000 - 80,000)
    c2 = B["Chr2"]
    B["Chr2"] = c2[:40_000] + revcomp(c2[40_000:80_000]) + c2[80_000:]

    # Chr3: move its last 30 kb onto the end of Chr1 (a translocation)
    c3 = B["Chr3"]
    moved, B["Chr3"] = c3[50_000:], c3[:50_000]
    B["Chr1"] = B["Chr1"] + moved

    # ---- strainC: derived from B, extra inversion + a real deletion ----
    C = {c: mutate(rng, s, 0.004) for c, s in B.items()}

    # Chr1: invert 50 kb (20,000 - 70,000)
    c1 = C["Chr1"]
    C["Chr1"] = c1[:20_000] + revcomp(c1[20_000:70_000]) + c1[70_000:]

    # Chr3: delete 20 kb — this is what makes the scale bar bug visible
    c3 = C["Chr3"]
    C["Chr3"] = c3[:15_000] + c3[35_000:]

    order = ["Chr1", "Chr2", "Chr3", "MT"]
    for strain, genome in (("strainA", A), ("strainB", B), ("strainC", C)):
        recs = [(c, genome[c]) for c in order]
        path = os.path.join(args.outdir, f"{strain}.fasta")
        write_fasta(path, recs, HEADER_STYLE[strain])
        total = sum(len(s) for _, s in recs)
        sizes = ", ".join(f"{c}={len(genome[c]):,}" for c in order)
        print(f"{path}  total={total:,} bp  ({sizes})")

    strains_txt = os.path.join(args.outdir, "strains.txt")
    with open(strains_txt, "w") as f:
        f.write("# strain_name<TAB>/path/to/assembly.fasta\n")
        for strain in ("strainA", "strainB", "strainC"):
            f.write(f"{strain}\t{os.path.abspath(os.path.join(args.outdir, strain + '.fasta'))}\n")
    print(f"{strains_txt}")


if __name__ == "__main__":
    main()
