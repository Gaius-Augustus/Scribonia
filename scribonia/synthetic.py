"""Synthetic genomes with a chosen stop set, for tests and augmentation.

Genes are random protein-like codon strings under a codon-usage model derived
from the requested genetic code: amino acids are drawn from a typical
frequency table, and each amino acid is spelled with one of its synonymous
codons (including the reassigned stop codons where the code says so).
Intergenic sequence is random with a chosen GC content, and both strands are
used.  Optional short introns mimic ciliate macronuclear genomes.

This is a toy: there is no codon-usage bias beyond the genetic code itself,
so a rare reassigned codon has to be requested explicitly with `reassigned_
usage`.
"""

import numpy as np

from .tables import CODONS

# Amino-acid frequencies (rough eukaryotic proteome averages).
AA_FREQ = {
    "A": 7.4, "R": 5.2, "N": 4.5, "D": 5.3, "C": 1.8, "Q": 4.1, "E": 6.7,
    "G": 6.6, "H": 2.3, "I": 5.3, "L": 9.3, "K": 6.0, "M": 2.3, "F": 3.9,
    "P": 5.1, "S": 8.1, "T": 5.6, "W": 1.2, "Y": 3.1, "V": 6.4,
}

STANDARD = {
    "A": ["GCT", "GCC", "GCA", "GCG"],
    "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "N": ["AAT", "AAC"], "D": ["GAT", "GAC"], "C": ["TGT", "TGC"],
    "Q": ["CAA", "CAG"], "E": ["GAA", "GAG"],
    "G": ["GGT", "GGC", "GGA", "GGG"], "H": ["CAT", "CAC"],
    "I": ["ATT", "ATC", "ATA"],
    "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "K": ["AAA", "AAG"], "M": ["ATG"], "F": ["TTT", "TTC"],
    "P": ["CCT", "CCC", "CCA", "CCG"],
    "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "T": ["ACT", "ACC", "ACA", "ACG"], "W": ["TGG"], "Y": ["TAT", "TAC"],
    "V": ["GTT", "GTC", "GTA", "GTG"],
}

# Which amino acid a reassigned codon encodes, per stop set (canonical
# nuclear choices: table 6 for {TGA}, 10 for {TAA,TAG}, 15 for {TAA,TGA}).
REASSIGNMENT = {
    frozenset({"TAA", "TAG", "TGA"}): {},
    frozenset({"TGA"}): {"TAA": "Q", "TAG": "Q"},
    frozenset({"TAA", "TAG"}): {"TGA": "C"},
    frozenset({"TAA", "TGA"}): {"TAG": "Q"},
    frozenset({"TAG"}): {"TAA": "Y", "TGA": "W"},
}

_COMP = str.maketrans("ACGT", "TGCA")


def _random_dna(rng, n, gc):
    p = np.array([(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2])
    return "".join(rng.choice(list("ACGT"), size=n, p=p))


def make_genome(stop_set, size_bp=300_000, gc=0.45, coding_fraction=0.6,
                gene_len_codons=(150, 800), intron_len=0, introns_per_gene=0,
                reassigned_usage=0.5, n_contigs=1, seed=0):
    """Return list of (name, sequence str) contigs.

    reassigned_usage: share of the synonymous codons that the reassigned
    codons take among their amino acid (0.5 = as frequent as the canonical
    codons combined; 0.05 = rare).
    """
    rng = np.random.default_rng(seed)
    stop_set = frozenset(stop_set)
    stops = [c for c in CODONS if c in stop_set]
    if not stops:
        raise ValueError("synthetic genomes need at least one stop codon")
    table = {aa: list(cods) for aa, cods in STANDARD.items()}
    weights = {}
    for aa, cods in table.items():
        weights[aa] = np.ones(len(cods)) / len(cods)
    for codon, aa in REASSIGNMENT[stop_set].items():
        cods = table[aa]
        w = weights[aa] * (1.0 - reassigned_usage)
        cods.append(codon)
        weights[aa] = np.concatenate((w, [reassigned_usage]))
    # Normalise in case several codons were added to one amino acid.
    for aa in weights:
        weights[aa] = weights[aa] / weights[aa].sum()
    aas = list(AA_FREQ)
    aa_p = np.array([AA_FREQ[a] for a in aas])
    aa_p = aa_p / aa_p.sum()

    def gene():
        n = int(rng.integers(gene_len_codons[0], gene_len_codons[1] + 1))
        picks = rng.choice(aas, size=n, p=aa_p)
        codons = ["ATG"]
        for aa in picks:
            codons.append(table[aa][rng.choice(len(table[aa]), p=weights[aa])])
        codons.append(stops[rng.integers(len(stops))])
        cds = "".join(codons)
        if intron_len and introns_per_gene:
            pieces, last = [], 0
            cuts = sorted(rng.integers(3, len(cds) - 3, size=introns_per_gene))
            for cut in cuts:
                pieces.append(cds[last:cut])
                pieces.append("GT" + _random_dna(rng, intron_len - 4, gc) + "AG")
                last = cut
            pieces.append(cds[last:])
            cds = "".join(pieces)
        return cds

    contig_size = size_bp // n_contigs
    contigs = []
    for k in range(n_contigs):
        parts, total = [], 0
        while total < contig_size:
            g = gene()
            if rng.random() < 0.5:
                g = g.translate(_COMP)[::-1]
            inter_len = int(len(g) * (1 - coding_fraction) / coding_fraction)
            inter = _random_dna(rng, max(20, int(rng.exponential(inter_len))), gc)
            parts.append(inter)
            parts.append(g)
            total += len(inter) + len(g)
        contigs.append((f"synthetic_{k}", "".join(parts)[:contig_size]))
    return contigs


def fragment(contigs, rng, median_len=25_000, sigma=1.0, min_len=5_000,
             max_len=2_000_000):
    """Cut sequences into log-normally distributed pieces."""
    out = []
    for name, seq in contigs:
        pos, k = 0, 0
        while pos < len(seq):
            L = int(np.clip(rng.lognormal(np.log(median_len), sigma),
                            min_len, max_len))
            piece = seq[pos:pos + L]
            if len(piece) >= min_len or pos == 0:
                out.append((f"{name}_{k}", piece))
            pos += L
            k += 1
    return out
