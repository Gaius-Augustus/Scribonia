# Tested and discarded ideas

Each entry gives what was tried, what happened and what the code does
instead. The protocol and terms (held-out species, leave-genus-out) are
explained in [experiments.md](experiments.md). Most tests date from
2026-09-26/27.

## Model

### Isotonic calibration of the probabilities

* **Tried:** each head was calibrated on leave-genus-out predictions (v1).
* **Result:** held-out exact match fell from 0.987 to 0.910. The calibrators
  were fitted on too few genera per non-standard class and flattened good
  probabilities.
* **Now:** no calibration by default. `scribonia train --calibrate` keeps it
  as an option.

### One model per codon

* **Tried:** a single gradient-boosting model per codon.
* **Result:** calls swung with the simulation seed. Stentor, for example,
  scored anywhere from 0.44 to 0.96 across seeds.
* **Now:** 10 models per codon, each fitted on a bootstrap sample of genera,
  and their outputs averaged.

### Independent thresholds at 0.5

* **Tried:** each codon called a stop when p ≥ 0.5.
* **Result:** this can produce stop sets that have no NCBI table, and it was
  slightly worse in a threshold sweep (0.9961 against 0.9972 exact on
  real-length held-out bins).
* **Now:** joint decoding over the stop sets that have a table
  (`decode_joint` in `scribonia/model.py`).

### Genome-composition features (GC, GC3, codon frequencies)

* **Tried:** keeping GC, GC3 and the three stop-codon frequencies as
  features.
* **Result:** once the {TGA} recodes were removed (next section), the model
  used composition as a shortcut. It called the 16 %-GC Ichthyophthirius
  {TAA,TGA}.
* **Now:** these five features are excluded (`EXCLUDED_FEATURES`).

### Contamination flag from contig votes

* **Tried:** flagging a bin as contaminated when contigs of at least 5 kb vote
  for different stop sets.
* **Result:** on bins of at least 0.3 Mb the flag fired on 91 % of clean bins
  (96 % and 99 % of the contaminated ones). On 10–300 kb bins it fired on
  33 % of clean bins. It does not separate contaminated from clean bins.
* **Now:** the `contamination_flag` column is still written, because
  downstream pipelines read the columns by name, but it should not be
  trusted. The per-contig vote
  features also added nothing measurable to accuracy (0.993 with and 0.994
  without them). They are kept, since they cost nothing.

## Training data

### Recoding standard genomes into {TGA}

* **Tried:** recoded CAA/CAG in annotated CDS into TAA/TAG, to add {TGA}
  training genomes.
* **Result:** even at a recoding fraction of 0.8, depletion of TAA/TAG over
  coding-like runs stayed at 0.1–0.4. Real {TGA} ciliates reach 0.45–0.8.
  These genomes taught the model to call AT-rich standard genomes
  (Dictyostelium, C. elegans) {TGA}.
* **Now:** only {TAA,TAG} and {TAA,TGA} are recoded. {TGA} comes from real
  genomes: ciliates, diplomonads and aphelids.

### Recoding fraction 0.05

* **Tried:** recoding only 5 % of the donor codons.
* **Result:** the recoded genomes are indistinguishable from standard ones and
  only add label noise.
* **Now:** fractions 0.25, 0.5 and 0.8 for {TAA,TAG}; 0.5, 0.8 and 1.0 for
  {TAA,TGA}.

### Intron-rich recoding sources

* **Tried:** C. elegans, Drosophila and Arabidopsis as sources for recoded
  genomes.
* **Result:** a 2–20 kb contig of these genomes holds little coding sequence,
  so the recoded copies still looked standard. The model then called Danio
  {TGA} or {TAA,TGA}.
* **Now:** recoding sources are compact genomes only (`RECODE_SOURCES`).

### Training on real contig lengths

* **Tried:** training bins with the real contig lengths of metagenomic bins
  instead of 2–20 kb.
* **Result:** leave-genus-out exact on protist bins of at least 1 Mb fell from
  0.973 to 0.956. The largest loss was Euplotes focardii, from 0.98 to 0.36.
  A model trained on 2–20 kb contigs does better on real-length bins than one
  trained on them.
* **Now:** training uses 2–20 kb; real lengths are for testing only.

### Training on small bins

* **Tried:** training bins of 10 kb–30 Mb instead of 0.3–30 Mb.
* **Result:** on held-out species, 10–300 kb bins improved from 0.87 to 0.93,
  and false non-standard calls on small standard bins fell from 15 % to 6 %.
  But leave-genus-out exact collapsed to 0.867: Euplotes focardii 0.06,
  Blepharisma 0.16, Tetrahymena 0.76.
* **Now:** bins of 0.3–30 Mb. Bins below 1 Mb are outside the intended use;
  a pipeline should take their code from taxonomy instead.

### Votes only from contigs of at least 20 kb

* **Tried:** per-contig votes only from contigs of at least 20 kb.
* **Result:** only 28 % of the 150 real metagenomic bins have such a contig,
  while all have one of at least 5 kb.
* **Now:** contigs of at least 5 kb vote (`VOTE_MIN_CONTIG`). This needed an
  overflow-safe sigmoid, because D of a short contig can be far above 1.

### Acetabularia as a {TGA} training genome

* **Tried:** Acetabularia acetabulum (GCA_963931725.1) as a training genome,
  labelled {TGA} after Cocquyt et al. 2010.
* **Result:** every feature says standard, and v2 called it standard in
  150 of 150 bins. Training on it would repeat the recoding error above.
* **Now:** held out; see [experiments.md](experiments.md#acetabularia).

### Other excluded genomes

These and three further candidates are listed with their reasons in
`data/excluded_genomes.tsv`.

* **Streblomastix strix (oxymonad, {TGA}):** a single-cell metagenome from a
  termite hindgut with contig N50 5 kb. Its bacterial contamination has not
  been assessed.
* **Blastophysa rhizopus (ulvophyte, {TGA}):** contig N50 277 bp, shorter than
  the 2 kb contig minimum.
* **Amoebophrya A25/A120:** labelled {TGA} earlier, but no paper states their
  code. Only the strain ex Karlodinium veneficum is shown to be non-standard.
  Details in [genetic_codes.md](genetic_codes.md).

## Compute

### GPU nodes

* **Considered:** running training on GPU nodes.
* **Why not:** simulation and training are CPU-only (numpy, scikit-learn).
  The available GPU nodes have 16 cores against 48 on the CPU nodes, and the
  job would block GPUs.
* **Now:** the speed-up comes from splitting the simulation over an array job
  (`scripts/submit_train.sh`): about 30 min instead of about 3.5 h.
