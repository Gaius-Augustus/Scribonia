"""Synthetic recoding of standard-code genomes into another stop set.

Training augmentation for thin classes: the {TAA,TAG} class has two real
genomes, {TAA,TGA} none.  A standard-code genome with a CDS annotation (NCBI
GFF3) is rewritten inside its annotated coding sequences only:

  * reassigned codons appear in the CDS interior: for each codon c that is
    a stop in the standard code but not in the target set, a fraction
    `fraction` of the synonymous codons one substitution away from c
    (amino acid from synthetic.REASSIGNMENT, e.g. CAA->TAA, CAG->TAG for
    Gln in table 6, TGT/TGC->TGA for Cys in table 10) are changed to c;
  * a terminal stop that is no longer a stop becomes the nearest codon of
    the target set.

Intergenic sequence, introns and unannotated genes are left alone.  CDS
that do not look like standard-code genes (length not a multiple of 3,
internal stop, no terminal stop) are skipped, as are codons overlapping a
position another CDS already changed.  `fraction` is an augmentation
parameter, not a measured usage: the evaluation runs on real held-out
genomes only.
"""

import gzip

import numpy as np

from .synthetic import REASSIGNMENT, STANDARD
from .tables import CODONS

_COMP = bytes.maketrans(b"ACGTacgt", b"TGCAtgca")


def _open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_fasta_bytes(path):
    """{name: bytearray} (upper case), names are the first header word."""
    seqs, name, chunks = {}, None, []
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = bytearray("".join(chunks).upper(), "ascii")
                name = line[1:].split()[0] if len(line) > 1 else ""
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = bytearray("".join(chunks).upper(), "ascii")
    return seqs


def read_cds(gff):
    """{cds id: (seqid, strand, [(start0, end0, phase), ...])} from GFF3 CDS
    lines.  Segments are grouped by ID (NCBI gives all segments of one
    protein the same ID), else by Parent.  Pseudogenes are skipped."""
    cds = {}
    with _open_text(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS" or f[6] not in "+-":
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            if attrs.get("pseudo") == "true":
                continue
            key = attrs.get("ID") or attrs.get("Parent")
            if key is None:
                continue
            phase = int(f[7]) if f[7] in "012" else 0
            seg = (int(f[3]) - 1, int(f[4]), phase)
            entry = cds.setdefault(key, (f[0], f[6], []))
            if entry[0] != f[0] or entry[1] != f[6]:
                continue   # trans-spliced across contigs/strands: ignore
            entry[2].append(seg)
    return cds


def _cds_positions(strand, segments):
    """Genome positions of the CDS in reading order, phase trimmed."""
    segs = sorted(segments, reverse=(strand == "-"))
    parts = [np.arange(s, e) if strand == "+" else np.arange(e - 1, s - 1, -1)
             for s, e, _ in segs]
    pos = np.concatenate(parts)
    return pos[segs[0][2]:]


def _donors(codon, aa):
    """Synonymous codons of aa one substitution away from codon."""
    near = [c for c in STANDARD[aa]
            if sum(a != b for a, b in zip(c, codon)) == 1]
    return near or list(STANDARD[aa])


def _nearest(codon, targets, rng):
    d = [sum(a != b for a, b in zip(codon, t)) for t in targets]
    best = [t for t, x in zip(targets, d) if x == min(d)]
    return best[int(rng.integers(len(best)))]


def recode_cds(seqs, cds, stop_set, fraction, seed=0):
    """Recode `seqs` ({name: bytearray}, changed in place) to `stop_set`.

    Returns counts: cds_used, cds_skipped, stops_changed and one entry per
    reassigned codon with the number of interior codons changed to it."""
    rng = np.random.default_rng(seed)
    stop_set = frozenset(stop_set)
    reassigned = REASSIGNMENT[stop_set]
    targets = sorted(c for c in CODONS if c in stop_set)
    # donor codon -> reassigned codon it becomes
    change = {}
    for c, aa in reassigned.items():
        for d in _donors(c, aa):
            change.setdefault(d, c)
    std_stops = set(CODONS)
    stats = {"cds_used": 0, "cds_skipped": 0, "stops_changed": 0}
    stats.update({c: 0 for c in reassigned})
    touched_by_seq = {}

    for key in sorted(cds):
        seqid, strand, segments = cds[key]
        seq = seqs.get(seqid)
        if seq is None:
            stats["cds_skipped"] += 1
            continue
        pos = _cds_positions(strand, segments)
        if len(pos) < 6 or len(pos) % 3 or pos.max() >= len(seq):
            stats["cds_skipped"] += 1
            continue
        raw = bytes(seq[i] for i in pos)
        if strand == "-":
            raw = raw.translate(_COMP)
        codons = [raw[i:i + 3].decode() for i in range(0, len(raw), 3)]
        if codons[-1] not in std_stops or std_stops.intersection(codons[:-1]):
            stats["cds_skipped"] += 1
            continue
        stats["cds_used"] += 1
        touched = touched_by_seq.setdefault(seqid, set())
        new = list(codons)
        for i, c in enumerate(codons[:-1]):
            if c in change and rng.random() < fraction:
                new[i] = change[c]
        if codons[-1] not in stop_set:
            new[-1] = _nearest(codons[-1], targets, rng)
        for i, (old, cod) in enumerate(zip(codons, new)):
            if old == cod:
                continue
            p = pos[3 * i:3 * i + 3]
            if any(int(x) in touched for x in p):
                continue
            b = cod.encode()
            if strand == "-":
                b = b.translate(_COMP)
            for x, ch in zip(p, b):
                seq[int(x)] = ch
                touched.add(int(x))
            if i == len(codons) - 1:
                stats["stops_changed"] += 1
            else:
                stats[cod] += 1
    return stats


def write_fasta_bytes(path, seqs, width=80):
    with open(path, "w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            s = seq.decode("ascii")
            for i in range(0, len(s), width):
                fh.write(s[i:i + width] + "\n")
