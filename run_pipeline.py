#!/usr/bin/env python3
"""HBB sequence pipeline: QC -> stats -> pairwise identity -> tree -> report.

Run the whole thing with one command:

    python run_pipeline.py

Every step reads the previous step's file from disk and writes its own
output file, so each intermediate result is kept in the results folder:

    results/01_qc/clean.fasta            sequences that passed QC
    results/01_qc/qc_log.tsv             pass/fail and reason for every record
    results/02_stats/seq_stats.tsv       length, GC %, base counts per sequence
    results/03_identity/identity_matrix.tsv   pairwise % identity (global alignment)
    results/04_tree/tree.nwk             UPGMA tree in Newick format
    results/05_report/report.md          human-readable summary

Only the Python standard library is used.  Same input + same parameters
gives byte-identical output files.
"""
import argparse
import hashlib
import os
import sys
from collections import Counter

# ---------------------------------------------------------------------------
# PARAMETERS - all settings live here.  Change them here, nowhere else.
# ---------------------------------------------------------------------------
INPUT_FASTA = "data/hbb_primates.fasta"   # input sequences (nucleotide FASTA)
OUTPUT_DIR = "results"                    # where all step folders are written
MIN_LENGTH = 200                          # QC: drop sequences shorter than this (nt)
MAX_N_FRACTION = 0.05                     # QC: drop sequences with more than 5% N
REFERENCE_ID = "Homo_sapiens"             # report: identities are ranked against this
MATCH_SCORE = 1                           # alignment: score for a matching base
MISMATCH_SCORE = -1                       # alignment: score for a mismatch
GAP_SCORE = -2                            # alignment: score for each gap position
BAND_WIDTH = 40                           # alignment: band half-width (speed-up)
# ---------------------------------------------------------------------------


# ----------------------------- helpers -------------------------------------
def read_fasta(path):
    records, name, chunks = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name, chunks = line[1:].strip(), []
            else:
                if name is None:
                    sys.exit(f"{path}: sequence data before first '>' header")
                chunks.append(line.upper())
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def write_fasta(path, records):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")


def short_id(name):
    """Label used in tables/tree: the part of the header before the first '|' or space."""
    return name.split("|")[0].split()[0]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        h.update(fh.read())
    return h.hexdigest()[:12]


# ----------------------------- step 1: QC ----------------------------------
def step_qc(fasta_in, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    records = read_fasta(fasta_in)
    if not records:
        sys.exit("No sequences found in the input file.")
    kept, log = [], ["id\tlength\tn_fraction\tstatus\treason"]
    for name, seq in records:
        n_frac = seq.count("N") / len(seq) if seq else 1.0
        reason = ""
        if len(seq) < MIN_LENGTH:
            reason = f"shorter than {MIN_LENGTH} nt"
        elif n_frac > MAX_N_FRACTION:
            reason = f"more than {MAX_N_FRACTION:.0%} N"
        elif set(seq) - set("ACGTNRYKMSWBDHVU"):
            reason = "invalid characters"
        status = "fail" if reason else "pass"
        if not reason:
            kept.append((name, seq))
        log.append(f"{short_id(name)}\t{len(seq)}\t{n_frac:.4f}\t{status}\t{reason}")
    clean_path = os.path.join(out_dir, "clean.fasta")
    write_fasta(clean_path, kept)
    with open(os.path.join(out_dir, "qc_log.tsv"), "w") as fh:
        fh.write("\n".join(log) + "\n")
    return clean_path, len(records), len(kept)


# ----------------------------- step 2: stats -------------------------------
def step_stats(clean_fasta, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    rows = ["id\tlength\tgc_percent\tA\tC\tG\tT\tN\tother"]
    for name, seq in read_fasta(clean_fasta):
        c = Counter(seq)
        acgt = c["A"] + c["C"] + c["G"] + c["T"]
        gc = 100.0 * (c["G"] + c["C"]) / acgt if acgt else float("nan")
        other = len(seq) - acgt - c["N"]
        rows.append(f"{short_id(name)}\t{len(seq)}\t{gc:.2f}\t{c['A']}\t{c['C']}\t{c['G']}\t{c['T']}\t{c['N']}\t{other}")
    path = os.path.join(out_dir, "seq_stats.tsv")
    with open(path, "w") as fh:
        fh.write("\n".join(rows) + "\n")
    return path


# ----------------------------- step 3: identity ----------------------------
def global_identity(a, b):
    """Banded Needleman-Wunsch. Returns percent identity = matches / alignment columns."""
    n, m = len(a), len(b)
    d = m - n
    NEG = -10**9
    H = [[NEG] * (m + 1) for _ in range(n + 1)]
    H[0][0] = 0
    for j in range(1, min(m, BAND_WIDTH + max(d, 0)) + 1):
        H[0][j] = j * GAP_SCORE
    for i in range(1, n + 1):
        lo = max(0, i + min(d, 0) - BAND_WIDTH)
        hi = min(m, i + max(d, 0) + BAND_WIDTH)
        if lo == 0:
            H[i][0] = i * GAP_SCORE
        row, prev, ai = H[i], H[i - 1], a[i - 1]
        for j in range(max(lo, 1), hi + 1):
            s = MATCH_SCORE if ai == b[j - 1] else MISMATCH_SCORE
            best = prev[j - 1] + s
            up = prev[j] + GAP_SCORE
            left = row[j - 1] + GAP_SCORE
            if up > best:
                best = up
            if left > best:
                best = left
            row[j] = best
    if H[n][m] <= NEG // 2:
        sys.exit("Sequence lengths differ by more than BAND_WIDTH; raise BAND_WIDTH.")
    i, j, matches, cols = n, m, 0, 0
    while i > 0 or j > 0:
        cur = H[i][j]
        if i > 0 and j > 0:
            s = MATCH_SCORE if a[i - 1] == b[j - 1] else MISMATCH_SCORE
            if cur == H[i - 1][j - 1] + s:
                matches += a[i - 1] == b[j - 1]
                cols += 1
                i, j = i - 1, j - 1
                continue
        if i > 0 and cur == H[i - 1][j] + GAP_SCORE:
            cols += 1
            i -= 1
        else:
            cols += 1
            j -= 1
    return 100.0 * matches / cols


def step_identity(clean_fasta, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    recs = read_fasta(clean_fasta)
    ids = [short_id(n) for n, _ in recs]
    k = len(recs)
    mat = [[100.0] * k for _ in range(k)]
    for x in range(k):
        for y in range(x + 1, k):
            mat[x][y] = mat[y][x] = global_identity(recs[x][1], recs[y][1])
    path = os.path.join(out_dir, "identity_matrix.tsv")
    with open(path, "w") as fh:
        fh.write("id\t" + "\t".join(ids) + "\n")
        for i, row in enumerate(mat):
            fh.write(ids[i] + "\t" + "\t".join(f"{v:.2f}" for v in row) + "\n")
    return path


# ----------------------------- step 4: tree --------------------------------
def read_matrix(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")[1:]
        mat = [[float(v) for v in line.rstrip("\n").split("\t")[1:]] for line in fh]
    return header, mat


def upgma(ids, dist):
    """UPGMA. dist is a full symmetric matrix of distances. Returns Newick string."""
    clusters = {i: (ids[i], 1, 0.0) for i in range(len(ids))}  # key -> (newick, size, height)
    d = {(i, j): dist[i][j] for i in range(len(ids)) for j in range(i + 1, len(ids))}
    nxt = len(ids)
    while len(clusters) > 1:
        (a, b), dab = min(d.items(), key=lambda kv: (kv[1], kv[0]))  # ties broken by index: deterministic
        na, sa, ha = clusters[a]
        nb, sb, hb = clusters[b]
        h = dab / 2.0
        newick = f"({na}:{h - ha:.5f},{nb}:{h - hb:.5f})"
        for c in list(clusters):
            if c in (a, b):
                continue
            dac = d[(min(a, c), max(a, c))]
            dbc = d[(min(b, c), max(b, c))]
            d[(c, nxt)] = (dac * sa + dbc * sb) / (sa + sb)
        d = {k: v for k, v in d.items() if a not in k and b not in k}
        del clusters[a], clusters[b]
        clusters[nxt] = (newick, sa + sb, h)
        nxt += 1
    return next(iter(clusters.values()))[0] + ";"


def step_tree(identity_tsv, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    ids, mat = read_matrix(identity_tsv)
    dist = [[(100.0 - v) / 100.0 for v in row] for row in mat]  # p-distance
    path = os.path.join(out_dir, "tree.nwk")
    with open(path, "w") as fh:
        fh.write(upgma(ids, dist) + "\n")
    return path


# ----------------------------- step 5: report ------------------------------
def step_report(paths, counts, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    n_in, n_kept = counts
    ids, mat = read_matrix(paths["identity"])
    stats = [line.rstrip("\n").split("\t") for line in open(paths["stats"])][1:]
    gcs = [float(r[2]) for r in stats]
    lens = [int(r[1]) for r in stats]
    pairs = [(mat[i][j], ids[i], ids[j]) for i in range(len(ids)) for j in range(i + 1, len(ids))]
    closest, farthest = max(pairs), min(pairs)
    ref = ids.index(REFERENCE_ID) if REFERENCE_ID in ids else 0
    ranked = sorted(((mat[ref][j], ids[j]) for j in range(len(ids)) if j != ref), reverse=True)
    lines = [
        "# Pipeline report",
        "",
        "## Parameters used",
        f"- Input: `{INPUT_FASTA}`",
        f"- QC: minimum length {MIN_LENGTH} nt, maximum N fraction {MAX_N_FRACTION}",
        f"- Alignment: match {MATCH_SCORE}, mismatch {MISMATCH_SCORE}, gap {GAP_SCORE}, band width {BAND_WIDTH}",
        f"- Reference for ranking: `{ids[ref]}`",
        "",
        "## Results",
        f"- Sequences read: {n_in}; passed QC: {n_kept}; removed: {n_in - n_kept}",
        f"- Length range: {min(lens)} to {max(lens)} nt",
        f"- GC content range: {min(gcs):.2f}% to {max(gcs):.2f}%",
        f"- Most similar pair: {closest[1]} and {closest[2]} ({closest[0]:.2f}% identity)",
        f"- Most different pair: {farthest[1]} and {farthest[2]} ({farthest[0]:.2f}% identity)",
        "",
        f"## Identity to {ids[ref]}",
        "",
        "| Sequence | % identity |",
        "|---|---|",
    ]
    lines += [f"| {name} | {val:.2f} |" for val, name in ranked]
    lines += [
        "",
        "## Tree (UPGMA, Newick)",
        "",
        "```",
        open(paths["tree"]).read().strip(),
        "```",
        "",
        "## Output files (with checksums for reproducibility)",
        "",
    ]
    for key in ("clean", "qc_log", "stats", "identity", "tree"):
        lines.append(f"- `{os.path.relpath(paths[key], os.path.dirname(out_dir))}`  (sha256 {sha256(paths[key])})")
    path = os.path.join(out_dir, "report.md")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


# ----------------------------- main ----------------------------------------
def main():
    global INPUT_FASTA
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--input", default=INPUT_FASTA, help=f"input FASTA (default: {INPUT_FASTA})")
    ap.add_argument("--outdir", default=OUTPUT_DIR, help=f"output folder (default: {OUTPUT_DIR})")
    args = ap.parse_args()
    INPUT_FASTA = args.input
    out = args.outdir

    paths = {}
    print(f"[1/5] QC            {INPUT_FASTA}")
    paths["clean"], n_in, n_kept = step_qc(INPUT_FASTA, os.path.join(out, "01_qc"))
    paths["qc_log"] = os.path.join(out, "01_qc", "qc_log.tsv")
    print(f"      {n_kept}/{n_in} sequences passed")
    print("[2/5] Sequence stats")
    paths["stats"] = step_stats(paths["clean"], os.path.join(out, "02_stats"))
    print("[3/5] Pairwise identity (global alignment)")
    paths["identity"] = step_identity(paths["clean"], os.path.join(out, "03_identity"))
    print("[4/5] UPGMA tree")
    paths["tree"] = step_tree(paths["identity"], os.path.join(out, "04_tree"))
    print("[5/5] Report")
    paths["report"] = step_report(paths, (n_in, n_kept), os.path.join(out, "05_report"))
    print(f"Done. See {paths['report']}")


if __name__ == "__main__":
    main()
