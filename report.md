# Pipeline report

## Parameters used
- Input: `data/hbb_primates.fasta`
- QC: minimum length 200 nt, maximum N fraction 0.05
- Alignment: match 1, mismatch -1, gap -2, band width 40
- Reference for ranking: `Homo_sapiens`

## Results
- Sequences read: 12; passed QC: 12; removed: 0
- Length range: 622 to 629 nt
- GC content range: 51.27% to 53.70%
- Most similar pair: Pongo_pygmaeus and Pongo_abelii (100.00% identity)
- Most different pair: Gorilla_gorilla and Nycticebus_coucang (86.31% identity)

## Identity to Homo_sapiens

| Sequence | % identity |
|---|---|
| Pan_troglodytes | 99.84 |
| Pan_paniscus | 99.68 |
| Gorilla_gorilla | 99.20 |
| Pongo_pygmaeus | 98.41 |
| Pongo_abelii | 98.41 |
| Symphalangus_syndactylus | 98.25 |
| Nomascus_leucogenys | 97.93 |
| Hylobates_moloch | 97.61 |
| Rhinopithecus_roxellana | 96.35 |
| Carlito_syrichta | 88.38 |
| Nycticebus_coucang | 86.78 |

## Tree (UPGMA, Newick)

```
(Nycticebus_coucang:0.06555,(Carlito_syrichta:0.05866,(Rhinopithecus_roxellana:0.01960,((Hylobates_moloch:0.00598,(Symphalangus_syndactylus:0.00400,Nomascus_leucogenys:0.00400):0.00198):0.00402,((Pongo_pygmaeus:0.00000,Pongo_abelii:0.00000):0.00795,(Gorilla_gorilla:0.00373,(Pan_paniscus:0.00120,(Homo_sapiens:0.00080,Pan_troglodytes:0.00080):0.00040):0.00253):0.00422):0.00204):0.00961):0.03906):0.00689);
```

## Output files (with checksums for reproducibility)

- `01_qc/clean.fasta`  (sha256 e46a8754203d)
- `01_qc/qc_log.tsv`  (sha256 2a26794cba7d)
- `02_stats/seq_stats.tsv`  (sha256 fe3d20a03bae)
- `03_identity/identity_matrix.tsv`  (sha256 b50748a6dcd7)
- `04_tree/tree.nwk`  (sha256 70c6ad81e12c)
