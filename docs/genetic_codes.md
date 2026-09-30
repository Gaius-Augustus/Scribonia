# Genetic codes: literature basis and training genomes

This file covers which eukaryotic lineages change the meaning of a nuclear stop
codon, which genomes Scribonia is trained and tested on, and which cases are
open. Every assignment below comes from a publication whose DOI and claim were
checked against the paper. `data/manifest.tsv` holds the per-genome labels
and `data/stop_codon_clades.tsv` the clade list. Candidate genomes that were
looked at and not used are in `data/excluded_genomes.tsv`, with the reason.

## Clades with nuclear stop-codon reassignments

Only bins from these clades, and bins without taxonomy, need to be
classified. Outside them no reassignment has been reported. That is absence of
reports rather than proof, but a screen of 265 non-ciliate transcriptomes from
the Marine Microbial Eukaryote Transcriptome Sequencing Project (MMETSP) found
no new code (Swart et al. 2016, doi:10.1016/j.cell.2016.06.020).

| clade (NCBI taxid) | stop sets found | standard relatives | key references |
| --- | --- | --- | --- |
| Ciliophora (5878) | {TGA}, {TAA,TAG}, {TAA,TGA}, ambiguous | Stentor, rumen Litostomatea, Fabrea, Chilodonella uncinata | 10.1016/j.cell.2016.06.020, 10.1093/molbev/msw166, 10.1371/journal.pgen.1010913, 10.1371/journal.pgen.1011512, 10.24072/pcjournal.141 |
| Fornicata (207245) | {TGA} in Hexamitinae; {TAA,TGA} in Iotanema | Giardia | 10.1093/oxfordjournals.molbev.a025832, 10.1371/journal.pgen.1004053, 10.1186/s12915-017-0353-y |
| Oxymonadida (66288) | {TGA} (Streblomastix) | Saccinobaculus | 10.1016/s0022-2836(03)00057-3 |
| Heterolobosea (5752) | {TAA,TGA} (Dactylomonas) | most Heterolobosea | 10.1016/j.ympev.2025.108289 |
| Trypanosomatidae (5654) | ambiguous (Blastocrithidia) | Obscuromonas, Trypanosoma, Leishmania | 10.1128/mbio.00885-25 |
| Ulvophyceae (33103) | {TGA} (Dasycladales, Cladophorales, Trentepohliales, Blastophysa); {TAA,TGA} (Scotinosphaerales) | Bryopsidales, Ulvales | 10.1186/1471-2148-10-327, 10.1093/gbe/evab101 |
| Cercozoa (136419) | {TAA,TGA} (one uncultured Sainouroidea lineage) | other Sainouroidea | 10.1186/s12915-017-0353-y |
| Aphelidiomycota (2316435) | {TGA} (Amoeboaphelidium protococcarum) | A. occidentale | 10.1016/j.cub.2022.08.071 |
| Syndiniales (88547) | ambiguous (Amoebophrya ex Karlodinium) | other Amoebophrya strains | 10.1371/journal.pone.0212912 |

Within ciliates, class membership does not predict the code. There are at
least five independent origins of TAA/TAG → Gln and at least three of
TGA → Trp (McGowan et al. 2023). Colpodea, Heterotrichea, Phyllopharyngea and
Litostomatea each contain several codes.

{TAA,TGA} (TAG read as sense) has arisen independently at least five times:
Phyllopharyngea, Scotinosphaerales, Iotanema, Dactylomonas and the
Sainouroidea lineage. NCBI table 15 covers this stop set, and Scribonia can
call it. No real genome of this class has a public assembly, so training uses
only recoded genomes.

## Training and test genomes

`data/manifest.tsv` has 81 rows, all with a checked citation: 60 primary
genomes for training and the held-out test (table below), three prokaryotic
contaminants, and the 18 stopCR ciliates used only in the
[benchmark](benchmark.md#samples). Recoded genomes are added at training
time.

| class | training genomes | held-out genomes |
| --- | --- | --- |
| {TGA} | Tetrahymena, Paramecium, Oxytricha, Stylonychia, Spironucleus, Amoeboaphelidium protococcarum (2 strains) | Ichthyophthirius, Pseudocohnilembus, Hexamita, Acetabularia |
| {TAA,TAG} | Euplotes focardii, Blepharisma stoltei; recoded genomes | Euplotes crassus |
| {TAA,TGA} | recoded genomes only | – |
| standard, relatives of non-standard lineages | 7 rumen ciliates, Fabrea, Giardia, Obscuromonas (2), A. occidentale | Stentor, Entodinium (2), Isotricha |
| standard, other | 21 protists, fungi, animals and plants | Plasmodium, Danio, Blastocystis |
| ambiguous (evaluation only) | – | Condylostoma, Blastocrithidia (2), Amoebophrya ex Karlodinium, Loxodes, Plagiopyla |

Standard relatives inside the non-standard clades matter most. They are the
negatives a bin from these clades is most likely to resemble.

## Corrections to NCBI and earlier assumptions

* **Blepharisma:** TGA encodes Trp and TAA/TAG are stops, giving stop set
  {TAA,TAG}. This is not NCBI table 15 (Singh et al. 2023,
  doi:10.1073/pnas.2213887120).
* **Fabrea salina:** standard code. NCBI had TGA = Cys (Heaphy et al. 2016,
  Table 1).
* **Amoebophrya:** only the strain ex Karlodinium veneficum is shown to read
  all three stops as sense (Bachvaroff 2019). The ex-Akashiwo control strain
  was canonical. No paper states the code of strains A25 and A120, so they are
  not in the manifest (`data/excluded_genomes.tsv`).
* **Platyophrya macrostoma and Colpoda aspera:** standard code (Heaphy et al.
  2016). NCBI lists Platyophrya as table 6.

## Open cases

* **Stylonychia and other stichotrichs.** Chen et al. 2023
  (doi:10.1093/molbev/msad064) report TGA = Arg, which would leave no
  dedicated stop codon. This contradicts the established {TGA} label, and
  Stylonychia lemnae is a training genome with that label. It has not been
  independently confirmed, so the label is unchanged. The same paper's own
  quantification (stopCR, Fig. 4c) puts about 98 % of Stylonychia's
  terminations at TGA, so TGA would be bifunctional and {TGA} remains the
  stop set a gene finder needs. The paper's per-species numbers are in
  `benchmark/chen2023_stopcr.tsv`. Codetta, a separate tool that also works
  from protein alignments, infers TGA = Arg in Sterkiella, Tetmemena and
  Uronychia binucleata, and Scribonia calls TGA a stop in all these genomes
  ([benchmark.md](benchmark.md#comparison-with-stopcr)).
* **Acetabularia.** Labelled {TGA}, but the assembly looks standard; see
  [experiments.md](experiments.md#acetabularia).
* **Palmarella salina (Armophorea).** TAA/TAG are reassigned, but the amino
  acid is not confirmed (the full text was inaccessible).
* **Chilodonella.** Reported as standard by McGowan et al. 2024. Another
  report of a "Chilodonella code" could not be traced to a verified source.
* **Genomes with no dedicated stop codon** (Condylostoma, karyorelicts,
  Blastocrithidia, Plagiopyla, Amoebophrya ex Karlodinium) are called {TGA}
  in all but 2 of 168 held-out bins. Detecting such codes would need a
  separate output.
