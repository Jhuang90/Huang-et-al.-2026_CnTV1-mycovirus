#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════╗
║    RIBBON RUNNER                                                         ║
║    Automated whole-genome synteny ribbon figure pipeline                 ║
║    Version 1.3.0  |  Jun Huang, Heitman Lab, Duke University             ║
║    Contact: jun.huang@duke.edu                                           ║
╚══════════════════════════════════════════════════════════════════════════╝

RIBBON RUNNER takes a list of genome assemblies, normalises their
chromosome names, runs adjacent-pair whole-genome alignments with nucmer,
and produces publication-quality SVG synteny ribbon figures in one command.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
REQUIRED SOFTWARE
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

    conda env create -f environment.yml
    conda activate ribbonrunner

  Tool          Version tested    Purpose
  ──────────    ──────────────    ──────────────────────────────────────
  Python        3.9+              Pipeline orchestration and plotting
  matplotlib    3.7+              SVG/PDF/PNG figure rendering
  nucmer        3.23 or 4.x       Whole-genome pairwise alignment
  delta-filter  (MUMmer suite)    Alignment quality filtering
  show-coords   (MUMmer suite)    Alignment coordinate extraction

  Since 1.3.0 samtools and seqkit are NO LONGER required: FASTA rewrapping
  and .fai indexing are done in Python by this script.

  synteny_plot.py must sit in the same directory as synteny_pipeline.py.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
CHROMOSOME NAME AUTO-DETECTION
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

    Format               Example          Normalised to
    ────────────────     ───────────      ─────────────
    Already standard     Chr1, MT         Chr1, MT
    Strain-prefixed      NRHc5010_Chr1    Chr1
    Numeric (H99-style)  1, 2, Mt         Chr1, Chr2, MT
    Lowercase            chr1             Chr1
    Long form            chromosome_1     Chr1

  Anything unrecognised is kept unchanged and reported. Two sequences that
  would normalise to the same name is a hard error, not a silent overwrite.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
USAGE
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  1. Create a TAB-separated strains file (order matters — strains are
     aligned as adjacent pairs, strain1 → strain2 → strain3 → ...):

       # strain_name    /path/to/assembly.fasta
       NRHc5010         /path/to/NRHc5010.fasta
       NRHc5028         /path/to/NRHc5028.fasta
       H99              /path/to/H99.fasta

  2. Run locally:

       python3 synteny_pipeline.py --input strains.txt --outdir out

     or submit to SLURM:

       python3 synteny_pipeline.py --input strains.txt --outdir out \\
           --title "My synteny figure" --chroms Chr1,Chr2,Chr3,Chr4 \\
           --email your@email.edu --slurm

  3. Output:
       synteny_all_chr.svg              — every chromosome found
       synteny_Chr1_Chr2_Chr3_Chr4.svg  — subset, if --chroms was given

  Try it on the bundled example first (runs in under a minute):

       cd examples && python3 make_example_data.py --outdir .
       python3 ../synteny_pipeline.py --input strains.txt --outdir demo
"""

import argparse
import os
import re
import shlex
import subprocess
import sys
import textwrap

__version__ = "1.3.0"

LINE_WIDTH = 60   # FASTA output line width; also written into the .fai


def die(msg):
    sys.exit(f"ERROR: {msg}")


# -------------------------------------------------------
# FASTA preparation: normalise names, rewrap, write .fai
#
# Done in Python rather than by shelling out to seqkit + samtools. That
# removes two dependencies, and — more importantly — keeps this logic in
# one testable place instead of duplicating it inside a bash heredoc.
# -------------------------------------------------------

def normalize_chrom_name(h):
    """Map an assembly's sequence name onto the ChrN / MT convention."""
    hl = h.lower()
    if re.match(r"^Chr\d+$", h) or h == "MT":
        return h, True
    m = re.match(r"^.+_Chr(\d+)$", h)
    if m:
        return f"Chr{m.group(1)}", True
    if re.match(r"^.+_MT$", h, re.IGNORECASE):
        return "MT", True
    if re.match(r"^\d+$", h):
        return f"Chr{h}", True
    if hl in ("mt", "mito", "mitochondria"):
        return "MT", True
    m = re.match(r"^chromosome[_\-]?(\d+)$", hl)
    if m:
        return f"Chr{m.group(1)}", True
    m = re.match(r"^chr(\d+)$", hl)
    if m:
        return f"Chr{m.group(1)}", True
    return h, False           # unrecognised — kept as is


def read_fasta(path):
    """Yield (name, sequence). Sequence is uppercased and unwrapped."""
    name, chunks = None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                header = line[1:].strip()
                if not header:
                    die(f"{path}: empty FASTA header")
                name, chunks = header.split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def prepare_fasta(in_path, out_path, line_width=LINE_WIDTH):
    """
    Write a renamed, uniformly wrapped copy of in_path, plus a valid
    5-column .fai index (name, length, offset, linebases, linewidth).
    """
    records, unrecognised = [], []
    seen = {}
    for orig, seq in read_fasta(in_path):
        new, ok = normalize_chrom_name(orig)
        if not ok:
            unrecognised.append(orig)
        if new in seen:
            die(f"{in_path}: '{orig}' and '{seen[new]}' both normalise to "
                f"'{new}'. Rename one of them in the assembly first — "
                f"continuing would silently merge two sequences.")
        seen[new] = orig
        records.append((new, seq))

    if not records:
        die(f"{in_path}: no sequences found")

    pos = 0
    with open(out_path, "w") as fo, open(out_path + ".fai", "w") as fi:
        for name, seq in records:
            header = f">{name}\n"
            fo.write(header)
            pos += len(header)
            offset = pos
            for i in range(0, len(seq), line_width):
                chunk = seq[i:i + line_width]
                fo.write(chunk + "\n")
                pos += len(chunk) + 1
            fi.write(f"{name}\t{len(seq)}\t{offset}\t{line_width}\t{line_width + 1}\n")

    renamed = {o: n for n, o in ((n, seen[n]) for n, _ in records) if o != n}
    print(f"  {os.path.basename(in_path)}: {len(records)} sequences"
          + (f", renamed {len(renamed)}" if renamed else ""))
    for o, n in sorted(renamed.items()):
        print(f"    {o} -> {n}")
    for u in unrecognised:
        print(f"    WARNING: unrecognised sequence name '{u}' — kept as is")


# -------------------------------------------------------
# Shell command builders. Everything user-supplied is quoted with
# shlex.quote(): paths with spaces work, and a strain name can no longer
# inject shell commands into the generated script.
# -------------------------------------------------------

def q(s):
    return shlex.quote(str(s))


def nucmer_cmd(ref, qry, prefix, threads=0):
    thread_arg = f"-t {int(threads)} " if threads else ""
    return f"nucmer --maxmatch {thread_arg}-p {q(prefix)} {q(ref)} {q(qry)}"


def filter_cmd(prefix, min_idy, min_len):
    return (f"delta-filter -m -i {min_idy} -l {min_len} "
            f"{q(prefix + '.delta')} > {q(prefix + '.filter.delta')}")


def coords_cmd(prefix):
    return (f"show-coords -THrd {q(prefix + '.filter.delta')} "
            f"> {q(prefix + '.coords')}")


def plot_cmd(python, script, config, output, title, chroms,
             min_len, min_idy, height, width, min_chrom_len):
    parts = [f"{q(python)} {q(script)}",
             f"--config {q(config)}",
             f"--output {q(output)}",
             f"--min_len {min_len}",
             f"--min_idy {min_idy}",
             f"--width {width}"]
    if title:
        parts.append(f"--title {q(title)}")
    if chroms:
        parts.append(f"--chroms {q(chroms)}")
    elif min_chrom_len:
        parts.append(f"--min_chrom_len {min_chrom_len}")
    if height:
        parts.append(f"--height {height}")
    return " \\\n    ".join(parts)


# -------------------------------------------------------
# Input parsing and validation
# -------------------------------------------------------

def parse_strains(path):
    if not os.path.exists(path):
        die(f"strains file not found: {path}")
    strains, names = [], set()
    with open(path) as f:
        for lineno, line in enumerate(f, 1):
            line = line.rstrip("\n").strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                die(f"line {lineno} of {path} is not tab-separated: {line!r}\n"
                    f"       Expected: strain_name<TAB>/path/to/assembly.fasta")
            strain, fasta = parts[0].strip(), parts[1].strip()
            if not strain or not fasta:
                die(f"line {lineno} of {path}: empty strain name or path")
            # Strain names become file name prefixes that are handed to nucmer.
            # nucmer builds its own shell commands internally and cannot cope
            # with spaces or metacharacters in a prefix, so reject them here
            # with a clear message rather than letting nucmer fail obscurely.
            if not re.match(r"^[A-Za-z0-9._+-]+$", strain):
                die(f"line {lineno} of {path}: strain name {strain!r} contains "
                    f"characters that cannot be used in a file name prefix.\n"
                    f"       Use only letters, digits, dot, underscore, plus "
                    f"and hyphen (e.g. A2-102-5, NRHc5010, H99).")
            if strain in names:
                die(f"line {lineno} of {path}: duplicate strain name "
                    f"{strain!r}. Strain names must be unique — they are used "
                    f"as output file names.")
            if not os.path.exists(fasta):
                die(f"line {lineno} of {path}: fasta not found: {fasta}")
            names.add(strain)
            strains.append({"strain": strain, "fasta": os.path.abspath(fasta)})
    if len(strains) < 2:
        die(f"{path}: need at least 2 strains, found {len(strains)}")
    return strains


# -------------------------------------------------------
# Build the pipeline shell script
# -------------------------------------------------------

def build_pipeline(strains, outdir, this_script, plot_script, args):
    py = "python3"
    lines = ["set -euo pipefail",
             f"mkdir -p {q(outdir)}",
             f"cd {q(outdir)}",
             ""]

    lines.append("echo '=== Step 1: Normalising chromosome names and indexing ==='")
    # Work with names relative to outdir. The script has already cd'd there,
    # and nucmer cannot handle spaces anywhere in the paths it is given —
    # keeping them bare file names means an output directory whose path
    # contains a space still works.
    renamed_fastas, fai_files = {}, {}
    for sd in strains:
        strain, fasta = sd["strain"], sd["fasta"]
        renamed = f"{strain}_renamed.fasta"
        renamed_fastas[strain] = renamed
        fai_files[strain] = os.path.join(outdir, renamed + ".fai")
        lines.append(f"echo 'Processing {strain}...'")
        lines.append(f"{q(py)} {q(this_script)} --prepare-fasta "
                     f"{q(fasta)} {q(renamed)}")
    lines.append("")

    lines.append("echo '=== Step 2: Adjacent-pair nucmer alignments ==='")
    coords_files = {}
    for i in range(len(strains) - 1):
        ref, qry = strains[i]["strain"], strains[i + 1]["strain"]
        prefix = f"{ref}_{qry}"                     # relative to outdir
        coords_files[ref] = os.path.join(outdir, prefix + ".coords")
        lines.append(f"echo 'Aligning {ref} vs {qry}...'")
        lines.append(nucmer_cmd(renamed_fastas[ref], renamed_fastas[qry],
                                prefix, args.threads))
        lines.append(filter_cmd(prefix, args.min_idy, args.min_len))
        lines.append(coords_cmd(prefix))
        lines.append(f"echo '  Alignment blocks: '$(wc -l < {q(prefix + '.coords')})")
        lines.append(
            f"echo '  Translocations:'; "
            f"awk 'NF>=11 && $10!=$11 {{print \"   \"$10\"->\"$11\"  \"$5\"bp\"}}' "
            f"{q(prefix + '.coords')} | head -20 || true")
        lines.append("")
    coords_files[strains[-1]["strain"]] = "-"

    lines.append("echo '=== Step 3: Writing plot config ==='")
    config_path = os.path.join(outdir, "plot_config.txt")
    lines.append(f"cat > {q(config_path)} << 'RIBBONRUNNER_EOF'")
    for sd in strains:
        s = sd["strain"]
        lines.append(f"{s}\t{fai_files[s]}\t{coords_files[s]}")
    lines.append("RIBBONRUNNER_EOF")
    lines.append("")

    lines.append("echo '=== Step 4: Generating figures ==='")
    full_output = os.path.join(outdir, "synteny_all_chr.svg")
    lines.append("echo 'Full genome plot...'")
    lines.append(plot_cmd(
        py, plot_script, config_path, full_output,
        args.title or "Whole genome synteny", None,
        args.min_len, args.min_idy,
        args.height or len(strains) * 2.5, args.width, args.min_chrom_len))
    lines.append("")

    if args.chroms:
        subset = os.path.join(
            outdir, f"synteny_{args.chroms.replace(',', '_')}.svg")
        lines.append(f"echo 'Subset plot ({args.chroms})...'")
        lines.append(plot_cmd(
            py, plot_script, config_path, subset,
            args.title or f"Synteny — {args.chroms}", args.chroms,
            args.min_len, args.min_idy,
            args.height or max(6, len(strains) * 2.0), args.width,
            args.min_chrom_len))
        lines.append("")

    lines.append("echo '=== Done ==='")
    lines.append(f"ls -lh {q(outdir)}/*.svg")
    return "\n".join(lines)


# -------------------------------------------------------
# SLURM wrapper
# -------------------------------------------------------

def default_conda_init():
    """Locate conda.sh instead of hardcoding one cluster's path."""
    base = os.environ.get("CONDA_PREFIX_1") or os.environ.get("CONDA_EXE", "")
    if base.endswith("/bin/conda"):
        base = base[:-len("/bin/conda")]
    for cand in ([os.path.join(base, "etc/profile.d/conda.sh")] if base else []) + [
            os.path.expanduser("~/miniforge3/etc/profile.d/conda.sh"),
            os.path.expanduser("~/miniconda3/etc/profile.d/conda.sh"),
            os.path.expanduser("~/anaconda3/etc/profile.d/conda.sh")]:
        if os.path.exists(cand):
            return cand
    return ""


def wrap_slurm(pipeline_body, args, outdir):
    log_dir = os.path.join(outdir, "logs")
    os.makedirs(log_dir, exist_ok=True)   # must exist before sbatch reads #SBATCH --output

    email_lines = ""
    if args.email:
        email_lines = (f"#SBATCH --mail-user={args.email}\n"
                       f"#SBATCH --mail-type=BEGIN,END,FAIL\n")

    conda_lines = ""
    if args.conda_env:
        conda_init = args.conda_init or default_conda_init()
        if not conda_init:
            die("could not locate conda.sh — pass --conda_init "
                "/path/to/etc/profile.d/conda.sh, or --conda_env '' to skip "
                "conda activation entirely.")
        conda_lines = (f"source {q(conda_init)}\n"
                       f"conda activate {q(args.conda_env)}\n")

    header = textwrap.dedent(f"""\
        #!/bin/bash
        #SBATCH --ntasks=1
        #SBATCH --cpus-per-task={max(1, args.threads) if args.threads else 1}
        #SBATCH --job-name=ribbonrunner
        #SBATCH --time={args.time}
        #SBATCH --mem={args.mem}
        {email_lines}#SBATCH --output={log_dir}/synteny_%j.out
        #SBATCH --error={log_dir}/synteny_%j.err
        #SBATCH --partition={args.partition}

        """) + conda_lines + "\n"
    return header + pipeline_body + "\n"


# -------------------------------------------------------
# Main
# -------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description="RIBBON RUNNER: assemblies → synteny ribbon figures",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    p.add_argument("--input",   help="Strains file (strain<TAB>fasta per line)")
    p.add_argument("--outdir",  help="Output directory")
    p.add_argument("--title",   default=None, help="Figure title")
    p.add_argument("--chroms",  default=None,
                   help="Chromosome subset for the second figure, "
                        "e.g. Chr1,Chr2,Chr3,Chr4")
    p.add_argument("--min_chrom_len", type=int, default=0,
                   help="Ignore sequences shorter than this (bp) in the "
                        "full-genome figure. Useful for draft assemblies.")
    p.add_argument("--min_len", type=int,   default=10000,
                   help="Minimum alignment block length, bp (default: 10000)")
    p.add_argument("--min_idy", type=float, default=90.0,
                   help="Minimum alignment identity, %% (default: 90)")
    p.add_argument("--threads", type=int, default=0,
                   help="Threads for nucmer (MUMmer 4 only; 0 = single-threaded, "
                        "required for nucmer 3.23)")
    p.add_argument("--height",  type=float, default=0, help="Figure height (0=auto)")
    p.add_argument("--width",   type=float, default=16, help="Figure width")
    p.add_argument("--slurm",      action="store_true", help="Submit as a SLURM job")
    p.add_argument("--slurm_only", action="store_true",
                   help="Write the SLURM script but do not submit it")
    p.add_argument("--mem",       default="32G")
    p.add_argument("--time",      default="24:00:00",
                   help="SLURM time limit (default: 24:00:00)")
    p.add_argument("--partition", default="common",
                   help="SLURM partition (default: common — change for your cluster)")
    p.add_argument("--conda_init", default=None,
                   help="Path to conda.sh (default: auto-detected)")
    p.add_argument("--conda_env",  default="ribbonrunner",
                   help="Conda env to activate in the SLURM job "
                        "(empty string = do not activate anything)")
    p.add_argument("--email", default=None, help="Email for SLURM notifications")
    p.add_argument("--prepare-fasta", nargs=2, metavar=("IN", "OUT"),
                   help=argparse.SUPPRESS)   # internal: used by the generated script
    p.add_argument("--version", action="version",
                   version=f"RIBBON RUNNER {__version__}")
    args = p.parse_args()

    # Internal mode: normalise one FASTA and write its index, then exit.
    if args.prepare_fasta:
        prepare_fasta(args.prepare_fasta[0], args.prepare_fasta[1])
        return

    if not args.input or not args.outdir:
        die("--input and --outdir are required")

    this_script = os.path.abspath(__file__)
    plot_script = os.path.join(os.path.dirname(this_script), "synteny_plot.py")
    if not os.path.exists(plot_script):
        die(f"synteny_plot.py not found next to this script "
            f"(expected {plot_script}). Both files are part of RIBBON RUNNER.")

    outdir = os.path.abspath(args.outdir)
    # nucmer interpolates the working directory into its own internal shell
    # commands, so it fails with an unhelpful "prenuc returned non-zero" if
    # that path contains a space or a shell metacharacter. Nothing this
    # script does can quote its way around that — catch it before a long job
    # dies for a reason nobody can diagnose from the log.
    if re.search(r"""[\s;|&$<>()'"`\\]""", outdir):
        die(f"--outdir contains a space or a shell metacharacter:\n"
            f"         {outdir}\n"
            f"       nucmer cannot run in such a directory. Choose an output "
            f"path made of letters, digits, '/', '.', '_', '+' and '-'.")

    strains = parse_strains(args.input)

    print(f"RIBBON RUNNER {__version__}")
    print(f"Strains ({len(strains)}):")
    for sd in strains:
        print(f"  {sd['strain']}: {sd['fasta']}")
    print("\nAlignment order (adjacent pairs only):")
    for i in range(len(strains) - 1):
        print(f"  {strains[i]['strain']} → {strains[i + 1]['strain']}")
    print()

    pipeline_body = build_pipeline(strains, outdir, this_script, plot_script, args)
    os.makedirs(outdir, exist_ok=True)

    if args.slurm or args.slurm_only:
        script_path = os.path.join(outdir, "ribbonrunner.sb")
        with open(script_path, "w") as f:
            f.write(wrap_slurm(pipeline_body, args, outdir))
        print(f"SLURM script written: {script_path}")
        if args.slurm and not args.slurm_only:
            r = subprocess.run(["sbatch", script_path],
                               capture_output=True, text=True)
            print(r.stdout.strip())
            if r.returncode != 0:
                sys.exit(r.stderr.strip() or "sbatch failed")
        else:
            print(f"Submit with: sbatch {script_path}")
    else:
        script_path = os.path.join(outdir, "run_pipeline.sh")
        with open(script_path, "w") as f:
            f.write("#!/bin/bash\n" + pipeline_body + "\n")
        os.chmod(script_path, 0o755)
        r = subprocess.run(["bash", script_path])
        if r.returncode != 0:
            sys.exit(f"pipeline failed (exit {r.returncode}); "
                     f"see the output above. Script: {script_path}")


if __name__ == "__main__":
    main()
