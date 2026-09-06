#!/usr/bin/env python3
"""
synteny_plot.py — publication-quality synteny ribbon figure for multiple
genome assemblies, drawn from nucmer show-coords output.

Part of RIBBON RUNNER. Version 1.3.0.
Jun Huang, Heitman Lab, Duke University — jun.huang@duke.edu

Usage:
    python3 synteny_plot.py --config plot_config.txt \\
        --output synteny.svg \\
        --title "My title" \\
        [--chroms Chr1,Chr2,Chr3,Chr4]

Config file (tab-separated, one line per strain, top row first):
    strain_name <TAB> /path/to/fai <TAB> /path/to/coords_vs_next_strain
    ...the last strain has "-" for coords.

Coords must be `show-coords -THrd` output (11 columns).

If --chroms is not given, the chromosome list is taken from the .fai
files themselves (union across strains, natural-sorted, organelles last).
"""

import argparse
import os
import re
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
from matplotlib.path import Path

__version__ = "1.3.0"

# -------------------------------------------------------
# Defaults
# -------------------------------------------------------
STRAIN_COLORS = [
    "#4fa3e0", "#e07b4f", "#4caf50", "#ab47bc",
    "#e91e8c", "#ff9800", "#00bcd4", "#8bc34a",
    "#3f51b5", "#795548", "#607d8b", "#c2185b",
]
SYNTENY_COLOR   = "#b8cfe8"
INVERSION_COLOR = "#f5a623"
TRANSLOC_COLOR  = "#e05c4a"
CHR_HEIGHT      = 0.016
MARGIN          = 0.10
MIN_LEN         = 10000
MIN_IDY         = 90.0
MERGE_GAP_FRAC  = 0.015   # merge blocks within this fraction of chr length

ORGANELLE_NAMES = {"MT", "MITO", "MITOCHONDRIA", "CHLOROPLAST", "PT"}


def die(msg):
    sys.exit(f"ERROR: {msg}")


# -------------------------------------------------------
# I/O
# -------------------------------------------------------

def read_fai(path):
    """Read a .fai index; returns {sequence_name: length}."""
    if not os.path.exists(path):
        die(f"fai index not found: {path}\n"
            f"       Run the pipeline's preparation step first, or fix the path "
            f"in the plot config.")
    sizes = {}
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                try:
                    sizes[p[0]] = int(p[1])
                except ValueError:
                    continue
    if not sizes:
        die(f"no sequences found in {path} — is it a valid .fai index?")
    return sizes


def strip_chr_prefix(name):
    """NRHc5010_Chr1 -> Chr1, A2-102-5_Chr3 -> Chr3, MT -> MT"""
    for tag in ("_Chr", "_MT"):
        if tag in name:
            return name[name.index(tag) + 1:]
    return name


def read_coords(path):
    """
    Parse `show-coords -THrd` output:
        S1 E1 S2 E2 LEN1 LEN2 %IDY FRM1 FRM2 TAG_REF TAG_QRY

    Unlike earlier versions this counts and reports malformed lines instead
    of silently skipping them, so a format mismatch cannot quietly produce
    an empty figure.
    """
    if not os.path.exists(path):
        die(f"coords file not found: {path}")

    alns, skipped, total = [], 0, 0
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            total += 1
            p = line.split("\t")
            if len(p) < 11:
                skipped += 1
                continue
            try:
                s1, e1 = int(p[0]), int(p[1])
                s2, e2 = int(p[2]), int(p[3])
                len1   = int(p[4])
                idy    = float(p[6])
                strand_q = int(p[8])
                ref = strip_chr_prefix(p[9])
                qry = strip_chr_prefix(p[10])
            except (ValueError, IndexError):
                skipped += 1
                continue

            inverted = (strand_q == -1)
            if s2 > e2:
                s2, e2 = e2, s2
            if s1 > e1:
                s1, e1 = e1, s1

            alns.append({
                "s1": s1, "e1": e1, "s2": s2, "e2": e2,
                "len1": len1, "idy": idy,
                "ref": ref, "qry": qry,
                "inverted": inverted,
                "is_transloc": (ref != qry),
            })

    if total and not alns:
        die(f"{path}: {total} lines read but none could be parsed.\n"
            f"       Expected 11 tab-separated columns from `show-coords -THrd`.\n"
            f"       Check that -T -H -r -d were all passed to show-coords.")
    if skipped:
        print(f"  WARNING: {path}: skipped {skipped}/{total} unparseable lines")
    return alns


def parse_config(path):
    if not os.path.exists(path):
        die(f"plot config not found: {path}")
    strains = []
    with open(path) as f:
        for lineno, line in enumerate(f, 1):
            line = line.rstrip("\n").strip()
            if not line or line.startswith("#"):
                continue
            p = line.split("\t")
            if len(p) < 2:
                die(f"{path} line {lineno}: expected at least "
                    f"'strain<TAB>fai', got {line!r}")
            coords = p[2] if len(p) > 2 and p[2].strip() != "-" else None
            strains.append({"strain": p[0], "fai": p[1], "coords": coords})
    if not strains:
        die(f"{path}: no strains found")
    names = [s["strain"] for s in strains]
    dup = {n for n in names if names.count(n) > 1}
    if dup:
        die(f"{path}: duplicate strain name(s): {', '.join(sorted(dup))}")
    return strains


# -------------------------------------------------------
# Chromosome selection
# -------------------------------------------------------

def chrom_sort_key(name):
    """Natural sort: Chr1 < Chr2 < Chr10; organelles last; others in between."""
    if name.upper() in ORGANELLE_NAMES:
        return (2, 0, name)
    m = re.match(r"^Chr(\d+)$", name)
    if m:
        return (0, int(m.group(1)), name)
    m = re.search(r"(\d+)$", name)
    if m:
        return (1, int(m.group(1)), name)
    return (1, 10 ** 9, name)


def auto_chroms(strain_data, min_chrom_len=0):
    """
    Chromosome list derived from the .fai files instead of a hardcoded
    Chr1..Chr14. A sequence is kept if its longest occurrence across the
    strains is at least min_chrom_len.
    """
    longest = defaultdict(int)
    for sd in strain_data:
        for name, length in sd["sizes"].items():
            longest[name] = max(longest[name], length)
    names = [n for n, ln in longest.items() if ln >= min_chrom_len]
    if not names:
        die("no sequences passed --min_chrom_len; lower it or set --chroms")
    return sorted(names, key=chrom_sort_key)


# -------------------------------------------------------
# Block merging — collapse adjacent alignments
# -------------------------------------------------------

def merge_blocks(alns, ref_sizes, gap_frac=MERGE_GAP_FRAC):
    """
    Merge adjacent alignments sharing (ref_chr, qry_chr, strand) into single
    larger blocks, for cleaner ribbons. Identity of a merged block is the
    length-weighted mean of its parts.
    """
    groups = defaultdict(list)
    for a in alns:
        groups[(a["ref"], a["qry"], a["inverted"])].append(a)

    merged = []
    for (rc, _qc, _inv), group in groups.items():
        group = sorted(group, key=lambda x: x["s1"])
        max_gap = ref_sizes.get(rc, 1) * gap_frac

        cur = dict(group[0])
        cur_idy_w = cur["idy"] * max(cur["len1"], 1)
        cur_w     = max(cur["len1"], 1)

        for nxt in group[1:]:
            gap = nxt["s1"] - cur["e1"]
            if 0 <= gap <= max_gap:
                cur["e1"] = max(cur["e1"], nxt["e1"])
                cur["s2"] = min(cur["s2"], nxt["s2"])
                cur["e2"] = max(cur["e2"], nxt["e2"])
                cur["len1"] = cur["e1"] - cur["s1"]
                w = max(nxt["len1"], 1)
                cur_idy_w += nxt["idy"] * w
                cur_w     += w
                cur["idy"] = cur_idy_w / cur_w
            else:
                if cur["len1"] > 0:
                    merged.append(cur)
                cur = dict(nxt)
                cur_w = max(cur["len1"], 1)
                cur_idy_w = cur["idy"] * cur_w
        if cur["len1"] > 0:
            merged.append(cur)

    return merged


# -------------------------------------------------------
# Layout
# -------------------------------------------------------

def compute_layout(strain_data, chroms):
    """
    Lay chromosomes out left to right, one row per strain.

    All strains share ONE scale (figure units per bp), derived from the
    longest genome. A strain whose sequences total less than the longest
    genome ends short of the right margin — that white space is a real
    length difference and is meant to be visible.
    """
    spans = [sum(sd["sizes"].get(c, 0) for c in chroms) for sd in strain_data]
    max_span = max(spans) if spans else 0
    if max_span == 0:
        die("no chromosome lengths found — check the fai paths and that the "
            "chromosome names in --chroms match the assemblies")

    n_gaps  = max(len(chroms) - 1, 0)
    gap_bp  = max_span * 0.006
    # One scale for every row: longest genome + its gaps fills the drawing area.
    scale   = (1.0 - 2 * MARGIN) / (max_span + n_gaps * gap_bp)

    layout = {}
    n = len(strain_data)
    for i, sd in enumerate(strain_data):
        y = 1.0 - (i + 0.5) / n
        layout[sd["strain"]] = {}
        x = MARGIN
        for c in chroms:
            clen = sd["sizes"].get(c, 0)
            layout[sd["strain"]][c] = (x, x + clen * scale, y)
            x += clen * scale + gap_bp * scale
    return layout, scale


# -------------------------------------------------------
# Drawing
# -------------------------------------------------------

def draw_chrom(ax, x0, x1, y, height, color):
    w = x1 - x0
    if w <= 0:
        return
    ax.add_patch(FancyBboxPatch(
        (x0, y - height / 2), w, height,
        boxstyle="round,pad=0.002",
        linewidth=0.5, edgecolor="white", facecolor=color,
        zorder=5, clip_on=False))


def draw_ribbon(ax, x1s, x1e, y1, x2s, x2e, y2, color, alpha, height,
                twist=False):
    """
    Smooth cubic-bezier ribbon between two chromosome rows.

    twist=True crosses the two edges, so an inverted block is drawn as a
    visible twist rather than a parallel band that merely differs in colour.
    """
    if x1s > x1e:
        x1s, x1e = x1e, x1s
    if x2s > x2e:
        x2s, x2e = x2e, x2s
    if twist:
        x2s, x2e = x2e, x2s

    y1b = y1 - height / 2      # bottom of upper chromosome
    y2t = y2 + height / 2      # top of lower chromosome
    y_cp1 = y1b + (y2t - y1b) * 0.33
    y_cp2 = y1b + (y2t - y1b) * 0.67

    verts = [
        (x1s, y1b),
        (x1s, y_cp1), (x2s, y_cp2), (x2s, y2t),
        (x2e, y2t),
        (x2e, y_cp2), (x1e, y_cp1), (x1e, y1b),
        (x1s, y1b),
    ]
    codes = [
        Path.MOVETO,
        Path.CURVE4, Path.CURVE4, Path.CURVE4,
        Path.LINETO,
        Path.CURVE4, Path.CURVE4, Path.CURVE4,
        Path.CLOSEPOLY,
    ]
    ax.add_patch(mpatches.PathPatch(
        Path(verts, codes), facecolor=color, edgecolor="none",
        alpha=alpha, zorder=1, clip_on=False))


# -------------------------------------------------------
# Main plot
# -------------------------------------------------------

def make_plot(strain_data, chroms, output, title=None,
              min_len=MIN_LEN, min_idy=MIN_IDY, figsize=(16, None),
              legacy_scale=False):

    n_strains = len(strain_data)
    fig_h = figsize[1] if figsize[1] else max(4, 3.0 * n_strains)
    fig, ax = plt.subplots(figsize=(figsize[0], fig_h))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    if legacy_scale:
        layout, bp_per_unit = _legacy_layout(strain_data, chroms)
    else:
        layout, bp_per_unit = compute_layout(strain_data, chroms)

    # ---- Chromosomes ----
    for i, sd in enumerate(strain_data):
        strain, color = sd["strain"], sd["color"]
        ly = layout[strain]
        for c in chroms:
            if c not in ly:
                continue
            x0, x1, y = ly[c]
            draw_chrom(ax, x0, x1, y, CHR_HEIGHT, color)
            if i == 0 and x1 > x0:
                ax.text((x0 + x1) / 2, y + CHR_HEIGHT / 2 + 0.010,
                        c.replace("Chr", ""), ha="center", va="bottom",
                        fontsize=7, color="#444444", fontfamily="DejaVu Sans")

        first = next((c for c in chroms if c in ly), None)
        if first:
            x0, _, y = ly[first]
            ax.text(x0 - 0.012, y, strain, ha="right", va="center",
                    fontsize=9, fontweight="bold", color=color, clip_on=False)

    # ---- Ribbons ----
    n_drawn = 0
    for i, sd in enumerate(strain_data[:-1]):
        ref_strain = sd["strain"]
        qry_strain = strain_data[i + 1]["strain"]
        ref_sizes  = sd["sizes"]
        qry_sizes  = strain_data[i + 1]["sizes"]

        alns = [a for a in sd.get("alns", [])
                if a["len1"] >= min_len and a["idy"] >= min_idy]
        alns = merge_blocks(alns, ref_sizes)
        alns = sorted(alns, key=lambda a: (a["is_transloc"], a["inverted"]))

        ly_ref, ly_qry = layout[ref_strain], layout[qry_strain]
        for aln in alns:
            rc, qc = aln["ref"], aln["qry"]
            if rc not in ly_ref or qc not in ly_qry:
                continue
            rx0, rx1, ry = ly_ref[rc]
            qx0, qx1, qy = ly_qry[qc]
            r_len = ref_sizes.get(rc, 0)
            q_len = qry_sizes.get(qc, 0)
            if r_len <= 0 or q_len <= 0:
                continue
            r_sc = (rx1 - rx0) / r_len
            q_sc = (qx1 - qx0) / q_len

            if aln["is_transloc"]:
                color, alpha = TRANSLOC_COLOR, 0.70
            elif aln["inverted"]:
                color, alpha = INVERSION_COLOR, 0.60
            else:
                color, alpha = SYNTENY_COLOR, 0.38

            draw_ribbon(ax,
                        rx0 + aln["s1"] * r_sc, rx0 + aln["e1"] * r_sc, ry,
                        qx0 + aln["s2"] * q_sc, qx0 + aln["e2"] * q_sc, qy,
                        color, alpha, CHR_HEIGHT,
                        twist=aln["inverted"])
            n_drawn += 1

    if n_drawn == 0:
        print("  WARNING: no ribbons drawn. Either the assemblies share no "
              "alignments passing the filters, or the chromosome names in the "
              "coords files do not match those in the fai files.")

    # ---- Scale bar (valid for every row: all rows share one scale) ----
    first_strain = strain_data[0]["strain"]
    first_chr = next((c for c in chroms if c in layout[first_strain]), None)
    if first_chr and bp_per_unit:
        _, _, top_y = layout[first_strain][first_chr]
        genome_bp = max(sum(sd["sizes"].get(c, 0) for c in chroms)
                        for sd in strain_data)
        bar_bp, label = _pick_scale_bar(genome_bp)
        bar_len = bar_bp * bp_per_unit
        bar_x = MARGIN
        bar_y = top_y - CHR_HEIGHT / 2 - 0.04

        ax.plot([bar_x, bar_x + bar_len], [bar_y, bar_y],
                color="#444444", linewidth=1.2, solid_capstyle="butt", zorder=6)
        tick_h = 0.008
        for xt in (bar_x, bar_x + bar_len):
            ax.plot([xt, xt], [bar_y - tick_h, bar_y + tick_h],
                    color="#444444", linewidth=1.0, zorder=6)
        ax.text(bar_x + bar_len / 2, bar_y - 0.018, label,
                ha="center", va="top", fontsize=8, color="#444444",
                fontfamily="DejaVu Sans")

    # ---- Legend ----
    handles = [mpatches.Patch(facecolor=sd["color"], label=sd["strain"])
               for sd in strain_data]
    handles += [
        mpatches.Patch(facecolor=SYNTENY_COLOR,   alpha=0.8, label="Syntenic"),
        mpatches.Patch(facecolor=INVERSION_COLOR, alpha=0.8, label="Inversion"),
        mpatches.Patch(facecolor=TRANSLOC_COLOR,  alpha=0.8, label="Translocation"),
    ]
    ax.legend(handles=handles, loc="lower right", frameon=True, framealpha=0.9,
              fontsize=8, ncol=2, borderpad=0.8)

    if title:
        ax.set_title(title, fontsize=11, fontweight="bold", pad=12,
                     color="#1a1a1a")

    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close()
    print(f"Saved: {output}")


def _pick_scale_bar(genome_bp):
    """Round scale-bar length appropriate to the genome size."""
    for bp, label in ((10_000_000, "10 Mb"), (1_000_000, "1 Mb"),
                      (100_000, "100 kb"), (10_000, "10 kb"),
                      (1_000, "1 kb")):
        if genome_bp >= bp * 5:
            return bp, label
    return 1_000, "1 kb"


def _legacy_layout(strain_data, chroms):
    """
    The pre-1.3.0 layout: every row stretched to the same width, scale bar
    taken from the first strain only. Kept ONLY so the documentation can
    show what the bug looked like (--legacy_scale). Do not use for figures.
    """
    spans = [sum(sd["sizes"].get(c, 0) for c in chroms) for sd in strain_data]
    max_span = max(spans) if spans else 0
    if max_span == 0:
        die("no chromosome lengths found")
    layout, bp_per_unit = {}, None
    n, n_gaps = len(strain_data), len(chroms) - 1
    for i, sd in enumerate(strain_data):
        y = 1.0 - (i + 0.5) / n
        total_len = sum(sd["sizes"].get(c, 0) for c in chroms)
        gap_bp = max_span * 0.006
        total_span = total_len + n_gaps * gap_bp
        scale = (1.0 - 2 * MARGIN) / total_span if total_span > 0 else 1.0
        if bp_per_unit is None:
            bp_per_unit = scale
        layout[sd["strain"]] = {}
        x = MARGIN
        for c in chroms:
            clen = sd["sizes"].get(c, 0)
            layout[sd["strain"]][c] = (x, x + clen * scale, y)
            x += clen * scale + gap_bp * scale
    return layout, bp_per_unit


# -------------------------------------------------------
# CLI
# -------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description="Synteny ribbon figure from nucmer show-coords output")
    p.add_argument("--config", required=True)
    p.add_argument("--chroms", default=None,
                   help="Comma-separated subset, e.g. Chr1,Chr2. "
                        "Default: all sequences found in the fai files.")
    p.add_argument("--min_chrom_len", type=int, default=0,
                   help="When --chroms is not given, ignore sequences shorter "
                        "than this (bp). Useful for draft assemblies with many "
                        "small contigs.")
    p.add_argument("--output",  default="synteny.svg")
    p.add_argument("--title",   default=None)
    p.add_argument("--min_len", type=int,   default=MIN_LEN)
    p.add_argument("--min_idy", type=float, default=MIN_IDY)
    p.add_argument("--width",   type=float, default=16)
    p.add_argument("--height",  type=float, default=0)
    p.add_argument("--legacy_scale", action="store_true",
                   help=argparse.SUPPRESS)  # reproduces the pre-1.3.0 bug
    p.add_argument("--version", action="version",
                   version=f"synteny_plot.py {__version__}")
    args = p.parse_args()

    config = parse_config(args.config)
    strain_data = []
    for i, entry in enumerate(config):
        strain_data.append({
            "strain": entry["strain"],
            "sizes":  read_fai(entry["fai"]),
            "alns":   read_coords(entry["coords"]) if entry["coords"] else [],
            "color":  STRAIN_COLORS[i % len(STRAIN_COLORS)],
        })

    if args.chroms:
        chroms = [c.strip() for c in args.chroms.split(",") if c.strip()]
        known = set().union(*[set(sd["sizes"]) for sd in strain_data])
        missing = [c for c in chroms if c not in known]
        if missing:
            print(f"  WARNING: requested chromosome(s) not present in any "
                  f"assembly: {', '.join(missing)}")
            print(f"           available: {', '.join(sorted(known, key=chrom_sort_key))}")
    else:
        chroms = auto_chroms(strain_data, args.min_chrom_len)
        print(f"Chromosomes (auto-detected from fai): {', '.join(chroms)}")
        if len(chroms) > 30:
            print(f"  WARNING: {len(chroms)} sequences will be drawn. Consider "
                  f"--min_chrom_len or an explicit --chroms.")

    make_plot(strain_data, chroms,
              output=args.output, title=args.title,
              min_len=args.min_len, min_idy=args.min_idy,
              figsize=(args.width, args.height if args.height > 0 else None),
              legacy_scale=args.legacy_scale)


if __name__ == "__main__":
    main()
