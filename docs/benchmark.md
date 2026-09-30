# Comparison with Codetta and PhyloFisher

**Status: final (2026-09-30).** All three tools have finished all samples of
both sets: 102 samples of the held-out set and 108 of the stopCR ciliates.

| set | | Scribonia | Codetta | PhyloFisher |
| --- | --- | --- | --- | --- |
| held-out set, 11 species | genomes exact | 10/11 | 8/11 | 10/11 |
| | bins exact | 50/55 | 33/55 | 39/55 |
| stopCR ciliates, 11 species | genomes exact | 8/11 | 10/11 | 9/11 |
| | bins exact | 36/55 | 41/55 | 26/55 |
| both sets | CPU time, all 210 samples | 0.6 h | 3,378 h | 3.4 h |

The accuracy rows count only species with a clear stop set. The other 13
species are reported separately below.

For how Scribonia itself is evaluated, see [experiments.md](experiments.md).

## Tools

Only tools that work on eukaryotic nuclear genomes are compared.
Prokaryote and phage tools (gTranslate, KACI, MgCod, prodigal-gv) are out of
scope. FACIL (Dutilh et al. 2011) is not run: Codetta works on the same
principle and was shown to outperform it (Shulgina & Eddy 2023), and FACIL's
web server no longer responds.

| tool | principle | training | output | stop set used here |
| --- | --- | --- | --- | --- |
| Scribonia (model scribonia_v3) | depletion of stop codons inside long open reading frames, per stop-set hypothesis | boosted trees on simulated bins | stop set | as output |
| Codetta 2.0 | Pfam profile HMMs aligned to the six-frame translation; the amino acid of each codon is decoded from the aligned columns | none (fixed profiles) | amino acid or `?` for each of the 64 codons | a stop codon counts as reassigned when Codetta assigns it an amino acid |
| PhyloFisher 1.2.14, `genetic_code_examiner` | tblastn of 240 curated orthologs; codons are counted at alignment columns conserved in > 70 % of the database | none (fixed database, v1.0) | counts of conserved amino acids per codon (bar plots) | a stop codon counts as reassigned at >= 5 conserved positions and >= 0.1 % of all conserved positions |

Codetta and PhyloFisher do not name a stop set. The rules in the last column
are ours. For Codetta the rule follows its own output: it does not infer that
a codon is a stop, so a stop codon without an inferred amino acid stays a
stop. For PhyloFisher the thresholds are a choice; a sweep over them is given
below. Both tools ran with their default parameters. PhyloFisher's queries
were Homo sapiens, Arabidopsis, Dictyostelium, Naegleria, Plasmodium and
Tetrahymena (for each ortholog, the first of these in the database).

## Samples

* **Held-out set:** the held-out species of Scribonia (never used in
  training), except Danio (a 1.4 Gb animal genome, outside the protist scope).
  17 species: 11 with a clear stop set, 6 with an ambiguous code (stop codons
  that also encode amino acids), scored separately.
* **stopCR ciliates:** 18 ciliate assemblies analysed by Chen et al. 2023,
  none used in training. 16 are spirotrich assemblies of gene-sized
  chromosomes (contig N50 0.9–3.3 kb); Stentor roeselii (141 kb) and
  Mesodinium rubrum (27 kb) have longer contigs.
  * 7 species of genera new to the model, with a clear stop set.
  * 4 Euplotes. Their genus is in the training set (E. focardii), so they are
    listed apart. The same holds for Euplotes crassus in the held-out set.
  * 7 species in which TGA reads Arg but still ends 49–98 % of the genes
    (Chen et al. 2023). They are scored separately against {TGA}, the codon
    that terminates.
* **Per species:** the whole assembly and 5 simulated bins, 102 + 108
  samples. Bins follow the held-out evaluation of
  [experiments.md](experiments.md#evaluation-protocol): real contig lengths of
  metagenomic bins, 1–30 Mb, and in half of the bins up to 30 % of sequence
  from a genome with another stop set.
* **Hardware:** Codetta and PhyloFisher ran on 16 cores per sample,
  Scribonia on one (it is single-threaded). Costs are compared in CPU time.

## Results on the held-out set

### All three tools

Species with a clear stop set, all samples.

| tool | genomes exact | bins exact | standard genomes called non-standard | standard bins called non-standard |
| --- | --- | --- | --- | --- |
| Scribonia | 10/11 | 50/55 | 0/6 | 0/30 |
| Codetta | 8/11 | 33/55 | 2/6 | 13/30 |
| PhyloFisher | 10/11 | 39/55 | 0/6 | 1/30 |

All three tools miss Acetabularia, the genome and every bin (see
[Acetabularia](#acetabularia)). Without Acetabularia, the bins are Scribonia
50/50, Codetta 33/50 and PhyloFisher 39/50.

### Per species

Genome calls ("ok" = correct), and correct bins out of 5.

| species | stop set | genome: Scribonia | Codetta | PhyloFisher | bins: Scribonia | Codetta | PhyloFisher |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Euplotes crassus | TAA,TAG | ok | ok | ok | 5/5 | 3/5 | 3/5 |
| Hexamita inflata | TGA | ok | ok | ok | 5/5 | 5/5 | 3/5 |
| Ichthyophthirius multifiliis | TGA | ok | ok | ok | 5/5 | 4/5 | 3/5 |
| Pseudocohnilembus persalinus | TGA | ok | ok | ok | 5/5 | 4/5 | 1/5 |
| Acetabularia acetabulum | TGA | standard | standard | standard | 0/5 | 0/5 | 0/5 |
| Stentor coeruleus | standard | ok | ok | ok | 5/5 | 3/5 | 5/5 |
| Entodinium bursa | standard | ok | TAA,TAG | ok | 5/5 | 1/5 | 5/5 |
| Entodinium caudatum | standard | ok | ok | ok | 5/5 | 3/5 | 5/5 |
| Isotricha intestinalis | standard | ok | ok | ok | 5/5 | 5/5 | 5/5 |
| Plasmodium falciparum | standard | ok | ok | ok | 5/5 | 2/5 | 4/5 |
| Blastocystis hominis | standard | ok | TAA,TGA | ok | 5/5 | 3/5 | 5/5 |

### Where each tool goes wrong

**Scribonia.** Only Acetabularia (1 genome, 5 bins). It never called a
standard genome or bin non-standard.

**Codetta: contamination.** Codetta pools the aligned codons of all sequences
in a sample. A contaminant with a reassigned code then decides the call for
the whole bin. 15 of its 17 bin errors outside Acetabularia are of this kind:

| bin of | contaminant (share of bp) | contaminant's code | Codetta's call |
| --- | --- | --- | --- |
| Stentor coeruleus | Paramecium tetraurelia (0.6 %) | TAA/TAG = Gln | TAA read as Gln |
| Stentor coeruleus | Stylonychia lemnae (1.7 %) | TAA/TAG = Gln | TAA and TAG read as Gln |
| Entodinium caudatum | Mycoplasma genitalium (0.3 %) | TGA = Trp | TGA read as Trp |
| Entodinium caudatum | Loxodes magnus (18.9 %) | TAA/TAG = Gln, TGA context-dependent | TAA read as Gln, TGA as Trp |
| Entodinium bursa | Stylonychia lemnae (9.6 %) | TAA/TAG = Gln | all three read as amino acids |
| Entodinium bursa | Pseudocohnilembus persalinus (14.4 %) | TAA/TAG = Gln | TAA read as Gln |
| Entodinium bursa | Euplotes crassus (5.8 %) | TGA = Cys | TGA read as Cys |
| Plasmodium falciparum | Mycoplasma genitalium (12.4 %) | TGA = Trp | TGA read as Trp |
| Plasmodium falciparum | Blastocrithidia nonstop (3.2 %) | all three sense, TAG = Glu, TGA = Trp | TAG read as Glu, TGA as Trp |
| Plasmodium falciparum | Paramecium tetraurelia (4.1 %) | TAA/TAG = Gln | TAA and TAG read as Gln |
| Blastocystis hominis | Spironucleus salmonicida (4.8 %) | TAA/TAG = Gln | TAA and TAG read as Gln |
| Euplotes crassus | Loxodes magnus (6.7 %) | TAA/TAG = Gln, TGA context-dependent | TAA read as Gln |
| Euplotes crassus | Oxytricha trifallax (6.8 %) | TAA/TAG = Gln | TAA and TAG read as Gln, no stop left |
| Ichthyophthirius multifiliis | Blepharisma stoltei (8.6 %) | TGA = Trp | TGA read as Trp, no stop left |
| Pseudocohnilembus persalinus | Blastocrithidia triatomae (6.8 %) | all three sense, TGA = Trp | TGA read as Trp, no stop left |

Scribonia called all 15 bins correctly. It was trained on bins with the
same kind of contamination, and its depletion statistics are dominated by the
bin's main genome.

**Codetta: Blastocystis and Entodinium bursa.** The other two bin errors and
two genome errors have no contaminant behind them. Codetta reads TAG as Gln
in the Blastocystis hominis genome (66 aligned positions) and in one bin
(36), and TGA as Gly in the Entodinium bursa genome (256) and in one bin
(104). Both genomes are labelled standard, and Scribonia and PhyloFisher
agree. We have not checked this further.

**PhyloFisher: few conserved positions.** A bin of 1–30 Mb yields only a few
hundred to a few thousand conserved positions. A reassigned codon that is
rare in the genome then gets 0 or 1 hits and stays a stop: 8 bins of {TGA}
species were called {TAG,TGA}, 2 Euplotes bins standard. One Plasmodium bin
was called {TAA,TAG}. Lowering the threshold to >= 1 hit (and >= 0.1 % of
conserved positions) is the best setting on this set: 45/55 bins, still below
Scribonia, and chosen on the test data itself.

### Genomes without a dedicated stop codon

In these species, stop codons also encode amino acids and termination
depends on context. Codetta calls all three codons sense in every sample of
Blastocrithidia (12), Condylostoma (6), Loxodes (6) and Plagiopyla (6), and
in 4 of 6 of Amoebophrya. PhyloFisher does so for Blastocrithidia (11 of 12
samples)
but calls Condylostoma and Plagiopyla mostly {TGA}. Scribonia calls {TGA}:
its output always keeps at least one stop codon (see
[genetic_codes.md](genetic_codes.md)). Neither answer is a stop set that a
gene finder can use directly, so these species are not scored.

### Cost

CPU time per sample, median over all 17 species (range in brackets).

| tool | genome | bin | peak memory | total CPU time |
| --- | --- | --- | --- | --- |
| Scribonia | 39 s | 4.7 s | 9.1 GB (genomes), 0.7 GB (bins) | 0.4 h |
| PhyloFisher | 111 s | 31 s | 0.8 GB | 1.7 h |
| Codetta | 22 h (5.7–164 h) | 78 min (11 min–20 h) | 2.0 GB | 1,223 h |

Per sample, Codetta took about 2,800 times the CPU time of Scribonia on the
genomes and about 1,100 times on the bins (median ratio).

## Results on the stopCR ciliates

### By group

Genomes and bins called exactly.

| group | species | Scribonia genomes | bins | Codetta genomes | bins | PhyloFisher genomes | bins |
| --- | --- | --- | --- | --- | --- | --- | --- |
| genus new to Scribonia | 7 | 5/7 | 20/35 | 6/7 | 23/35 | 5/7 | 16/35 |
| Euplotes (genus in training) | 4 | 3/4 | 16/20 | 4/4 | 18/20 | 4/4 | 10/20 |
| clear stop set, total | 11 | 8/11 | 36/55 | 10/11 | 41/55 | 9/11 | 26/55 |
| TGA = Arg, scored against {TGA} | 7 | 5/7 | 27/35 | 0/7 | 22/35 | 5/7 | 15/35 |

On this set Codetta is the most accurate tool on the species with a clear
stop set. Scribonia is the most accurate on the TGA = Arg species, where
Codetta leaves no stop codon at all.

### Per species, stopCR ciliates

Genome calls ("ok" = correct, "none" = no stop codon left), and correct bins
out of 5.

| species | stop set | genome: Scribonia | Codetta | PhyloFisher | bins: Scribonia | Codetta | PhyloFisher |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Laurentiella sp. | TGA | ok | none | ok | 5/5 | 2/5 | 2/5 |
| Halteria grandinella | TGA | ok | ok | TAA,TGA | 3/5 | 4/5 | 0/5 |
| Pseudokeronopsis carnea | TGA | standard | ok | TAA,TGA | 0/5 | 2/5 | 0/5 |
| Diophrys appendiculata | TGA | ok | ok | ok | 4/5 | 3/5 | 2/5 |
| Diophrys sp. 3 | TGA | ok | ok | ok | 3/5 | 5/5 | 4/5 |
| Mesodinium rubrum | TGA | TAA,TGA | ok | ok | 0/5 | 5/5 | 4/5 |
| Stentor roeselii | standard | ok | ok | ok | 5/5 | 2/5 | 4/5 |
| Euplotes woodruffi | TAA,TAG | ok | ok | ok | 5/5 | 4/5 | 2/5 |
| Euplotes parawoodruffi | TAA,TAG | ok | ok | ok | 5/5 | 5/5 | 3/5 |
| Euplotes weissei | TAA,TAG | standard | ok | ok | 1/5 | 4/5 | 0/5 |
| Euplotes cf. woodruffi 2 | TAA,TAG | ok | ok | ok | 5/5 | 5/5 | 5/5 |
| Sterkiella histriomuscorum | TGA (= Arg) | ok | none | ok | 5/5 | 5/5 | 2/5 |
| Tetmemena sp. | TGA (= Arg) | ok | none | ok | 5/5 | 4/5 | 1/5 |
| Paraurostyla sp. | TGA (= Arg) | ok | none | ok | 5/5 | 2/5 | 5/5 |
| Urostyla sp. | TGA (= Arg) | ok | none | ok | 5/5 | 3/5 | 4/5 |
| Uronychia setigera | TGA (= Arg) | ok | none | ok | 4/5 | 3/5 | 2/5 |
| Uronychia binucleata | TGA (= Arg) | TAA,TGA | none | TAA,TGA | 3/5 | 1/5 | 1/5 |
| Pseudokeronopsis flava | TGA (= Arg) | standard | none | TAA,TGA | 0/5 | 4/5 | 0/5 |

### Where each tool goes wrong, stopCR ciliates

**Scribonia** misses 3 of the 11 genomes with a clear stop set, and 2 of the
7 TGA = Arg genomes.

* **Pseudokeronopsis carnea and P. flava** are called standard in all 12
  samples. The assemblies are contaminated (see
  [Pseudokeronopsis](#pseudokeronopsis-contaminated-assemblies)).
* **Mesodinium rubrum** (6 samples) and **Uronychia binucleata** (genome and
  2 bins) are called {TAA,TGA}. In both, TAA encodes an amino acid and also
  ends many genes: 42 % and 49 % of terminations (Chen et al. 2023).
  Scribonia reports the codons that terminate. The label {TGA} counts this
  as an error.
* **Euplotes weissei** is called standard in the genome and 4 bins. TGA is
  rarer inside long open reading frames than in the other three Euplotes
  (D 0.86 against 1.03–1.08), the call is marked low confidence, and TGA is
  flagged as ambiguous (p_stop 0.71). The assembly has the shortest contigs
  of the set (N50 0.9 kb). We have not checked the cause.
* **Bins:** the remaining 6 errors are bins called {TAA,TGA} instead of
  {TGA} (Halteria 2, Diophrys 3, Uronychia setigera 1). All 6 contain
  6–18 % of sequence from a genome in which TAA is a stop.

Scribonia called no standard sample non-standard (Stentor roeselii 6/6).

**Codetta: TGA = Arg.** Codetta reads TGA as Arg in all 7 genomes in which
stopCR does (58–1,142 aligned positions). TAA and TAG are Gln there, so no
stop codon is left. The codon meaning is right, but the output names no
codon that ends a gene, although TGA ends 49–98 % of them. In 22 of the 35
bins of these species Codetta finds too few positions to decode TGA (0–32).
TGA then stays a stop, and the bin counts as correct.

**Codetta: Laurentiella.** Codetta also reads TGA as Arg in the Laurentiella
genome (116 positions) and in 3 of its bins (16–32), and in one Halteria bin
(9). stopCR calls TGA a stop in both. Laurentiella is a stichotrich like
Sterkiella and Tetmemena, so a rare Arg reading may be real. We have not
checked this. On all other codons of the 18 genomes, Codetta and stopCR
agree.

**Codetta: bins.** Of the other 10 bin errors in species with a clear stop
set, 6 follow a contaminant's code (Stentor roeselii 3, Diophrys
appendiculata 2, Euplotes weissei 1), as on the held-out set. In 3 small bins
(1.3–2.6 Mb) Codetta left a reassigned codon undecoded (Pseudokeronopsis
carnea 2, Euplotes woodruffi 1). One Pseudokeronopsis carnea bin has TGA read
as Trp from 7 positions.

**PhyloFisher: rare codons and few positions.** The 4 genome errors are all
{TAA,TGA}: TAA is read as Gln at 5–21 conserved positions, which is below
0.1 % of all conserved positions. In the bins the counts are smaller still:
half of the Euplotes bins are called standard (TGA read as an amino acid at
0–5 conserved positions), and Halteria, Pseudokeronopsis and Euplotes weissei have no
correct bin. In the 7 TGA = Arg genomes PhyloFisher finds TGA at one
conserved position in total.

**PhyloFisher: thresholds.** Requiring >= 1 hit and >= 0.03 % of conserved
positions gives 11/11 genomes and 39/55 bins here. On the held-out set the
same setting gives 8/11 genomes and 42/55 bins, and calls 7 of 30 standard
bins non-standard. The best setting of the held-out set (>= 1 hit, >= 0.1 %)
gives 9/11 genomes and 34/55 bins here. No single threshold works on both
sets.

### Cost, stopCR ciliates

CPU time per sample, median over all 18 species (range in brackets).

| tool | genome | bin | peak memory | total CPU time |
| --- | --- | --- | --- | --- |
| Scribonia | 18 s | 4.1 s | 0.5 GB | 0.2 h |
| PhyloFisher | 135 s | 32 s | 0.3 GB | 1.6 h |
| Codetta | 98 h (8.5–192 h) | 1.9 h (10 min–20 h) | 1.3 GB | 2,156 h |

Per sample, Codetta took about 18,000 times the CPU time of Scribonia on the
genomes and about 2,000 times on the bins (median ratio). On 16 cores, a
Codetta genome run took 20 h of wall time in the median.

## Acetabularia

All three tools call Acetabularia standard, although the literature labels it
{TGA} (TAA/TAG = Gln; Schneider et al. 1989; Cocquyt et al. 2010). Two checks
explain this.

* **The label is right, but the reassigned codons are rare.** PhyloFisher on
  the published transcriptome (ENA TSA HCDB01, project PRJEB40460, 246,083
  transcripts; Andresen et al. 2021b) finds
  TAA at 32 and TAG at 31 conserved positions, always at Gln. At the same
  positions, CAG occurs 463 times and CAA 206 times. TAA/TAG thus carry about
  9 % of Gln. In the {TGA} ciliates they carry most of it: 66 % in
  Pseudocohnilembus, 75 % in Ichthyophthirius.
* **The genome assembly holds few Acetabularia genes.** GCA_963931725.1 is a
  meta-assembly of DNA amplified from five single cells (Andresen et al.
  2021a). It covers an estimated 7–11 % of the genome, about 1.8 % of it is
  transcribed, and 11 % of BUSCO genes are found. PhyloFisher finds only 53 of
  its 240 genes there, and only 6 of the 53 hit contigs match an Acetabularia
  transcript (>= 97 % identity over >= 200 bp). At its 131 conserved Gln
  positions the assembly shows CAA or CAG only.

Scribonia also calls the transcriptome standard (depletion D of TAA 0.13, of
TAG 0.23). A reassigned codon used for a tenth of one amino acid depletes
too few stop codons from open reading frames to be seen. This is a real limit
of the method, not only of the assembly.

## Comparison with stopCR

stopCR (Chen et al. 2023) aligns six-frame translations to Pfam and counts,
for each stop codon, how often it occurs inside conserved domains (codon
meaning) and how often it follows a conserved C-terminus (termination
share). A codon can have both: stichotrich TGA reads Arg, yet about 95 % of
genes end with it. Scribonia predicts the codons that terminate, so it
should call such codons stop.

`benchmark/compare_stopcr.py` puts Scribonia's calls next to the per-species
stopCR numbers (`benchmark/chen2023_stopcr.tsv`) for the 20 held-out
genomes that stopCR also analysed. It uses the whole genomes and 5 bins per
genome. Codons are classed as stop (stopCR: '*'), bifunctional (amino-acid
meaning, termination share >= 5 %) or sense.

| stopCR class | codons | called stop, genomes | called stop, bins | median p_stop, genomes |
| --- | --- | --- | --- | --- |
| stop | 21 | 100 % | 100 % | 0.999 |
| bifunctional | 13 | 62 % | 63 % | 0.980 |
| sense | 26 | 15 % | 17 % | 0.059 |

Spearman's rho between p_stop and the termination share is 0.68 on genomes
and 0.66 on bins.

* **Stichotrich TGA** (Arg per Chen et al., 93–98 % of terminations):
  called stop in all samples. Codetta infers TGA = Arg in all seven of these
  genomes and in Laurentiella (see
  [Results on the stopCR ciliates](#results-on-the-stopcr-ciliates)).
* **Bifunctional TAA** is called stop only where it ends many genes:
  Uronychia binucleata (49 % of terminations) is called {TAA,TGA}, while
  Uronychia setigera (30 %) and Diophrys (7–17 %) are called {TGA}.
* **Euplotes TGA = Cys** (up to 4 % of terminations): never called stop.
* **Pseudokeronopsis carnea and P. flava** are called standard in all 12
  samples, although stopCR reads TAA/TAG as Gln (at most 3.6 % of
  terminations) and Codetta infers Gln in both genomes. These account for 4
  of the 26 sense codons called stop. Halteria (TAA in 2 of 5 bins) accounts
  for the rest.

### Pseudokeronopsis: contaminated assemblies

Both assemblies contain long contigs that do not fit gene-sized ciliate
macronuclear chromosomes. Scribonia was run separately on contigs below and
above 5 kb:

| genome | contigs | Mb | N50 | call |
| --- | --- | --- | --- | --- |
| P. carnea | < 5 kb | 49.1 | 1.7 kb | {TAA,TGA} |
| P. carnea | >= 5 kb | 27.7 | 48.7 kb | standard |
| P. flava | < 5 kb | 41.8 | 1.4 kb | {TAA,TGA} |
| P. flava | >= 5 kb | 10.7 | 48.2 kb | standard |
| Laurentiella (control) | < 5 kb | 37.1 | 2.4 kb | {TGA} |
| Laurentiella (control) | >= 5 kb | 12.0 | 6.6 kb | {TGA} |

The longest contigs are 0.59 Mb (P. carnea) and 1.18 Mb (P. flava). In
P. carnea, 177 of the 200 longest contigs vote standard. Only contigs of at
least 5 kb vote, so these contigs set the whole-genome call although they
make up only 36 % and 20 % of the bp. Scribonia raises its contamination
flag, but still reports the contaminant's stop set. The flag does not single
these two out: it is raised for 17 of the 18 stopCR assemblies. stopCR avoids
the problem by keeping only contigs with ciliate protein hits.

On the short contigs TAA still looks like a stop (depletion D 0.22–0.27).
Either short contaminant contigs remain, or TAA is too rare in these genes
to be seen, as in Acetabularia. Telling the two apart needs protein
homology.

## Reproducing

The scripts are in `benchmark/` of the repository. `setup.sh` fetches the
tools and databases, `setup.sbatch` builds them and runs a test (Codetta must
infer TGA = Trp on its own example genome), and `submit_bench.sh` runs the
inputs, the three tools and the scoring as batch jobs. Paths are set through
environment variables listed in each script's header. `benchmark/score.py`
writes `summary.tsv`, `per_sample.tsv` and the PhyloFisher threshold sweep
`pf_sweep.tsv`. `benchmark/compare_stopcr.py` reads the same results and puts
Scribonia's and Codetta's per-codon calls next to the stopCR numbers of Chen
et al. 2023 (`benchmark/chen2023_stopcr.tsv`, doi:10.1093/molbev/msad064),
with each codon classed as stop, bifunctional (it has an amino-acid meaning
but ends at least 5 % of genes) or sense.

### Versions behind the numbers

Both runs were recorded with a modified working tree. The modifications
were checked file by file against the repository afterwards:

| | held-out set | stopCR ciliates |
| --- | --- | --- |
| Scribonia package and model | 0.1.0 source (contigs scanned one by one), model scribonia_v3 | 0.1.1 source (contigs scanned in batches), model scribonia_v3 |
| modified at run time | five scripts in `benchmark/` | the stopCR rows of `data/manifest.tsv`, `benchmark/chen2023_stopcr.tsv`, the species list in `scripts/train.sbatch` |
| state of those files | the code in this repository (0.1.1 changed only comments and the version record of the score stage) | as in this repository in every field the benchmark reads; the notes of three rows and the `table` field of the Mesodinium row were edited later |
| Codetta | commit 863359e | commit 863359e |
| PhyloFisher | bioconda package 1.2.14 (reports itself as 1.2.13) | same |

The Scribonia package differs between the two releases only in how contigs
are scanned (0.1.1 scans them in batches). The features are the same up to
float rounding (`test_batched_runs_match_contig_by_contig_scan`). The 102
held-out samples were classified again with 0.1.1: every column of the
output, including the probabilities and D values as printed, is the same as
in the recorded run. `benchmark/bench.sbatch` now saves the difference to
the commit whenever the working tree is not clean.

## References

* Andresen IJ, Orr RJS, Krabberød AK, Shalchian-Tabrizi K, Bråte J (2021a).
  Genome sequencing and de novo assembly of the giant unicellular alga
  Acetabularia acetabulum using droplet MDA. Scientific Reports 11:12820.
  doi:10.1038/s41598-021-92092-4
* Andresen IJ, Orr RJS, Shalchian-Tabrizi K, Bråte J (2021b).
  Compartmentalization of mRNAs in the giant, unicellular green alga
  Acetabularia acetabulum. Algal Research 59:102440.
  doi:10.1016/j.algal.2021.102440
* Chen et al. (2023). Stop or Not: Genome-Wide Profiling of Reassigned Stop
  Codons in Ciliates. Molecular Biology and Evolution 40:msad064.
  doi:10.1093/molbev/msad064
* Cocquyt E, Gile GH, Leliaert F, Verbruggen H, Keeling PJ, De Clerck O
  (2010). Complex phylogenetic distribution of a non-canonical genetic code in
  green algae. BMC Evolutionary Biology 10:327. doi:10.1186/1471-2148-10-327
* Dutilh BE, et al. (2011). FACIL: Fast and Accurate genetic Code Inference
  and Logo. Bioinformatics 27:1929–1933. doi:10.1093/bioinformatics/btr316
* Schneider SU, Leible MB, Yang XP (1989). Strong homology between the small
  subunit of ribulose-1,5-bisphosphate carboxylase/oxygenase of two species
  of Acetabularia and the occurrence of unusual codon usage. Molecular and
  General Genetics 218:445–452. doi:10.1007/BF00332408
* Shulgina Y, Eddy SR (2023). Codetta: predicting the genetic code from
  nucleotide sequence. Bioinformatics 39:btac802.
  doi:10.1093/bioinformatics/btac802
* Tice AK, et al. (2021). PhyloFisher: a phylogenomic package for resolving
  eukaryotic relationships. PLoS Biology 19:e3001365.
  doi:10.1371/journal.pbio.3001365
