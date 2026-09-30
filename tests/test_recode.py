import csv

import numpy as np
import pytest

from scribonia.cli import main
from scribonia.recode import (_cds_positions, read_cds, read_fasta_bytes,
                              recode_cds, write_fasta_bytes)
from scribonia.synthetic import AA_FREQ, STANDARD

_COMP = str.maketrans("ACGT", "TGCA")
STOPS = ("TAA", "TAG", "TGA")
CODE = {c: aa for aa, cods in STANDARD.items() for c in cods}


def _standard_genome_with_gff(tmp_path, size_bp, seed=0, gc=0.4):
    """Standard-code genome with genes on both strands, every third gene with
    an intron, and a matching GFF3.  Returns (fasta, gff)."""
    rng = np.random.default_rng(seed)
    aas = list(AA_FREQ)
    p = np.array([AA_FREQ[a] for a in aas]); p = p / p.sum()
    nt = np.array([(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2])

    def dna(n):
        return "".join(rng.choice(list("ACGT"), size=n, p=nt))

    parts, gff, pos, g = [], [], 0, 0
    while pos < size_bp:
        inter = dna(int(rng.integers(200, 800)))
        parts.append(inter); pos += len(inter)
        n = int(rng.integers(150, 500))
        cds = "ATG" + "".join(
            STANDARD[a][rng.integers(len(STANDARD[a]))] for a in rng.choice(aas, size=n, p=p)
        ) + STOPS[rng.integers(3)]
        segs = [(0, len(cds))]          # sense coordinates of the CDS pieces
        gene = cds
        if g % 3 == 0:
            cut = int(rng.integers(10, len(cds) - 10))
            intron = "GT" + dna(80) + "AG"
            gene = cds[:cut] + intron + cds[cut:]
            segs = [(0, cut), (cut + len(intron), len(gene))]
        strand = "+" if g % 2 == 0 else "-"
        L = len(gene)
        done = 0
        for a, b in segs:
            phase = (3 - done % 3) % 3
            done += b - a
            s0, e0 = (pos + a, pos + b) if strand == "+" else (pos + L - b, pos + L - a)
            gff.append(f"chr1\ttest\tCDS\t{s0 + 1}\t{e0}\t.\t{strand}\t{phase}\tID=cds-{g};Parent=rna-{g}\n")
        parts.append(gene if strand == "+" else gene.translate(_COMP)[::-1])
        pos += L; g += 1
    fa, gf = tmp_path / "std.fna", tmp_path / "std.gff"
    write_fasta_bytes(fa, {"chr1": bytearray("".join(parts), "ascii")})
    gf.write_text("##gff-version 3\n" + "".join(gff))
    return fa, gf


def _cds_codons(seqs, cds):
    out = {}
    for key, (seqid, strand, segs) in cds.items():
        pos = _cds_positions(strand, segs)
        s = "".join(chr(seqs[seqid][i]) for i in pos)
        if strand == "-":
            s = s.translate(_COMP)
        out[key] = [s[i:i + 3] for i in range(0, len(s), 3)]
    return out


def _translate(codons, extra):
    return "".join(extra.get(c) or CODE.get(c, "*") for c in codons)


@pytest.mark.parametrize("stop_set,extra,donors", [
    ({"TGA"}, {"TAA": "Q", "TAG": "Q"}, {"CAA", "CAG"}),
    ({"TAA", "TAG"}, {"TGA": "C"}, {"TGT", "TGC"}),
])
def test_recode_full(tmp_path, stop_set, extra, donors):
    fa, gff = _standard_genome_with_gff(tmp_path, 60_000, seed=1)
    before = read_fasta_bytes(fa)
    after = read_fasta_bytes(fa)
    cds = read_cds(gff)
    stats = recode_cds(after, cds, stop_set, fraction=1.0, seed=0)
    assert stats["cds_used"] == len(cds) and stats["cds_skipped"] == 0
    old, new = _cds_codons(before, cds), _cds_codons(after, cds)
    for k in cds:
        # same protein under the target code, stop from the target set
        assert _translate(new[k][:-1], extra) == _translate(old[k][:-1], {})
        assert new[k][-1] in stop_set
        # every donor codon was changed (fraction 1)
        assert not donors.intersection(new[k][:-1])
    # only CDS positions changed
    in_cds = np.zeros(len(before["chr1"]), dtype=bool)
    for seqid, strand, segs in cds.values():
        in_cds[_cds_positions(strand, segs)] = True
    diff = np.frombuffer(bytes(before["chr1"]), np.uint8) != np.frombuffer(bytes(after["chr1"]), np.uint8)
    assert diff.any() and not (diff & ~in_cds).any()


def test_recode_skips_non_standard_cds(tmp_path):
    fa, gff = _standard_genome_with_gff(tmp_path, 20_000, seed=2)
    seqs = read_fasta_bytes(fa)
    cds = read_cds(gff)
    # plant an internal stop in the first CDS
    key = sorted(cds)[0]
    seqid, strand, segs = cds[key]
    pos = _cds_positions(strand, segs)
    for x, ch in zip(pos[3:6], b"TAA" if strand == "+" else b"TTA"[::-1]):
        seqs[seqid][int(x)] = ch
    stats = recode_cds(seqs, cds, {"TGA"}, fraction=0.5)
    assert stats["cds_skipped"] == 1


def test_recoded_genome_is_called_tga(tmp_path):
    fa, gff = _standard_genome_with_gff(tmp_path, 400_000, seed=3, gc=0.35)
    out = tmp_path / "recoded.fna"
    main(["recode", str(fa), str(gff), "--stop-set", "TGA", "--fraction", "0.5",
          "-o", str(out)])
    tsv = tmp_path / "call.tsv"
    main(["classify", str(out), "--rules", "-o", str(tsv)])
    with open(tsv) as fh:
        row = next(csv.DictReader(fh, delimiter="\t"))
    assert row["stop_set"] == "TGA"
