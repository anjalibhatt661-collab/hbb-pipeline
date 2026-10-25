#!/usr/bin/env python3
"""Provenance helper: build data/hbb_primates.fasta from a BLAST text report.

The pipeline does NOT need this script to run; the FASTA it produces is
already committed in data/. It is included so the origin of the input data
is transparent and repeatable.

The BLAST report is the Task 02 search (RID DVTA26TX014) of human HBB
(NM_000518.5) against core_nt.  For each chosen accession the aligned
subject region is taken from the first alignment block (gap characters
removed), which gives the HBB mRNA region of that species.

Usage:
    python tools/extract_from_blast.py DVTA26TX014-Alignment.txt data/hbb_primates.fasta
"""
import re
import sys

# accession -> label used in the FASTA header
SELECTED = {
    "NM_000518.5": "Homo_sapiens",
    "XM_508242.5": "Pan_troglodytes",
    "XM_003819029.5": "Pan_paniscus",
    "XM_019036164.3": "Gorilla_gorilla",
    "XM_054440830.2": "Pongo_pygmaeus",
    "XM_002822127.6": "Pongo_abelii",
    "XM_055282707.1": "Symphalangus_syndactylus",
    "XM_004090649.3": "Nomascus_leucogenys",
    "XM_032166808.2": "Hylobates_moloch",
    "XM_010363344.2": "Rhinopithecus_roxellana",
    "XM_008052701.1": "Carlito_syrichta",
    "XM_053561585.1": "Nycticebus_coucang",
}


def main(report, out):
    seqs, current, in_first_hsp = {}, None, False
    for line in open(report):
        m = re.match(r"Sequence ID:\s*(\S+)", line)
        if m:
            current = m.group(1) if m.group(1) in SELECTED and m.group(1) not in seqs else None
            in_first_hsp = current is not None
            if current:
                seqs[current] = []
            continue
        if current and line.startswith(("Score", "Score:")) and seqs[current]:
            in_first_hsp = False  # second HSP of the same subject starts
        if current and in_first_hsp and line.startswith("Sbjct"):
            seqs[current].append(line.split()[2].upper().replace("-", ""))
    missing = [a for a in SELECTED if a not in seqs or not seqs[a]]
    if missing:
        sys.exit(f"Not found in report: {missing}")
    with open(out, "w") as f:
        for acc, label in SELECTED.items():
            seq = "".join(seqs[acc])
            f.write(f">{label}|{acc}\n")
            for i in range(0, len(seq), 70):
                f.write(seq[i:i + 70] + "\n")
    print(f"Wrote {len(SELECTED)} sequences to {out}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
