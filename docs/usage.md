# Installation and usage

## Install

Default: the container image built from the `Dockerfile`
(`gaiusaugustus/scribonia`), run with Singularity/Apptainer:

```bash
singularity pull scribonia.sif docker://gaiusaugustus/scribonia:0.1.1
singularity run scribonia.sif classify genome.fa -o result.tsv
```

Docker works the same way (see the README). Without a container:

```bash
pip install .            # numpy only: rule-based classifier
pip install ".[train]"   # + scikit-learn, joblib: shipped model and training
```

The shipped model was trained with scikit-learn 1.9.1 (the version pinned in
the `Dockerfile`); other versions load it with a version warning.

## Classify

```bash
scribonia classify genome.fa -o result.tsv
scribonia classify bin07.fa --sample lake_2020 --bag Tetrahymena -o result.tsv
scribonia classify *.fa -o all.tsv --per-contig votes.tsv
scribonia classify genome.fa --rules -v  # rule-based fallback, no model
```

`--sample` and `--bag` are optional labels for the first two output columns
(see below); they do not affect the classification.

The trained model is loaded from `scribonia/models/` (file name
`DEFAULT_MODEL` in `scribonia/model.py`). Without it, or with `--rules`, the
rule-based classifier (`scribonia/rules.py`) is used. `scribonia --help`
lists the other subcommands (`features`, `simulate`, `recode`, `train`, `evaluate`).

## Output

One row per input, tab separated. Column names are stable, since downstream
pipelines read them by name (`OUTPUT_COLUMNS` in `scribonia/cli.py`).

| column | meaning |
| --- | --- |
| `sample` | label given with `--sample`; empty if not given |
| `bag` | label given with `--bag`; the input file name if not given |
| `size_bp`, `n_contigs`, `n50` | size of the input in bp, number of contigs, contig N50 |
| `gc`, `gc3` | GC content of the input and at third codon positions of the long open runs |
| `n_orf_eval` | number of long open runs the call rests on |
| `p_taa_stop`, `p_tag_stop`, `p_tga_stop` | probability that the codon terminates translation |
| `d_taa`, `d_tag`, `d_tga` | depletion D of the codon inside long open runs (about 0 for a stop, about 1 for a sense codon; `nan` without evidence) |
| `stop_set` | most probable stop set among those with an NCBI table, e.g. `TGA` for most ciliates |
| `table` | canonical NCBI table of the stop set: 1, 6, 10, 15 or 33 |
| `confidence` | `high`, `medium` or `low` (small bins, little evidence) |
| `ambiguous_codons` | codons with intermediate probability despite much evidence; `none` otherwise |
| `contamination_flag` | `yes` when contigs disagree about the stop set; not reliable, see [discarded ideas](discarded_ideas.md#contamination-flag-from-contig-votes) |
| `minority_vote_frac` | largest length-weighted share of contigs voting against the call of a codon |
| `classifier_version` | model version (`scribonia_v3`) or `rules-0.4` for the rule-based fallback |

The target is the stop set, not the NCBI table. Tables that share a stop set
give the same gene structures: 6, 27, 29 and 30 all stop only at TGA and
differ only in the amino acid that TAA/TAG encode. `table` is the canonical
table of the stop set; the amino acid of a reassigned codon has to come from
taxonomy.

## Tests

```bash
python3 -m pytest tests
```

The tests use synthetic genomes only (`scribonia/synthetic.py`). The model
tests, including those of the shipped model (`tests/test_shipped_model.py`),
need scikit-learn and are skipped without it. They run in the release image
(which has no pytest, hence the install):

```bash
docker run --rm -v "$PWD":/repo -w /repo --entrypoint sh \
    gaiusaugustus/scribonia:0.1.1 \
    -c "pip install -q pytest && python -m pytest tests -rs -p no:cacheprovider"
```

The GitHub workflow `.github/workflows/tests.yml` does the same on every
push, and fails if a test is skipped there. What Scribonia gets wrong is
listed in the [README](../README.md#limitations).
