# HBB Scripted Analysis Pipeline

A small pipeline that takes a set of nucleotide sequences through five chained steps and finishes with a report. It runs from **one command** and keeps the output of every step.

```bash
python run_pipeline.py
```

Python 3.6+ and nothing else: the standard library only, nothing to install.

## What the pipeline does

| Step | What happens | Output (kept in `results/`) |
|---|---|---|
| 1. QC | Drops sequences that are too short, have too many `N`, or contain invalid characters | `01_qc/clean.fasta`, `01_qc/qc_log.tsv` |
| 2. Stats | Length, GC % and base counts for each sequence | `02_stats/seq_stats.tsv` |
| 3. Identity | Pairwise percent identity from a global (Needleman-Wunsch) alignment | `03_identity/identity_matrix.tsv` |
| 4. Tree | UPGMA tree built from the identity matrix | `04_tree/tree.nwk` |
| 5. Report | Summary of the parameters, results, the tree and file checksums | `05_report/report.md` |

Each step reads the file written by the previous step, so any intermediate result can be opened and checked on its own.

## Input data

`data/hbb_primates.fasta` holds 12 primate beta-globin (HBB) mRNA sequences: human, chimpanzee, bonobo, gorilla, two orangutans, three gibbons, a snub-nosed monkey, a tarsier and a slow loris. They are the aligned HBB regions from my BLAST search of human HBB (`NM_000518.5`) against NCBI `core_nt` (job RID DVTA26TX014), where the identity percentages can be compared with this pipeline's results. `tools/extract_from_blast.py` shows how the file was made. It is not needed to run the pipeline.

To use your own data, point the pipeline at any nucleotide FASTA:

```bash
python run_pipeline.py --input my_sequences.fasta --outdir my_results
```

## Parameters

All settings are in one block at the top of `run_pipeline.py`:

| Parameter | Default | Meaning |
|---|---|---|
| `INPUT_FASTA` | `data/hbb_primates.fasta` | input file |
| `OUTPUT_DIR` | `results` | output folder |
| `MIN_LENGTH` | 200 | QC: minimum sequence length (nt) |
| `MAX_N_FRACTION` | 0.05 | QC: maximum fraction of `N` |
| `REFERENCE_ID` | `Homo_sapiens` | sequence the report ranks the others against |
| `MATCH_SCORE` / `MISMATCH_SCORE` / `GAP_SCORE` | 1 / -1 / -2 | alignment scoring |
| `BAND_WIDTH` | 40 | alignment band half-width (speed-up) |

## Reproducibility

- Run it from a fresh checkout with the single command above. There are no external downloads and no random numbers.
- I ran the pipeline twice from scratch and the output files were byte-identical. The report lists a checksum for each output file so you can check this yourself.
- Running it takes about one second.

## Results on the included data

- All 12 sequences passed QC.
- Identity to human: chimpanzee 99.84%, bonobo 99.68%, gorilla 99.20%, orangutans 98.41%, gibbons 97.6 to 98.3%, snub-nosed monkey 96.35%, tarsier 88.38%, slow loris 86.78%.
- The human, chimpanzee and bonobo values match the identities BLAST reported, which is a useful check that the alignment step is correct.
- The tree groups human with chimpanzee and bonobo, then gorilla, then the orangutans, then the gibbons, with the monkey, tarsier and slow loris further out. This matches the accepted primate relationships.
- The two orangutan sequences are 100% identical over this region, which is why their branch length is zero.

Full numbers are in `results/05_report/report.md`.

## Decisions and limitations

- **Banded alignment.** The sequences are close relatives of similar length, so a band of 40 positions around the diagonal is enough and keeps the run fast. If sequences differ in length by more than the band, the script stops with a message asking you to raise `BAND_WIDTH`.
- **Identity** is matches divided by alignment columns, counting gap columns, as BLAST does.
- **UPGMA** assumes a constant rate of change across lineages, so it is a quick illustration of relationships, not a rigorous phylogeny. For a publication-quality tree use a maximum-likelihood method.
- **Input is assumed to be nucleotide FASTA.** There is no protein mode.
- The input contains one gene per species. The pipeline does not check that it is the same gene.
