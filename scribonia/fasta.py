"""Minimal streaming FASTA reader producing integer-encoded sequences.

A, C, G, T map to 0..3 (case-insensitive); everything else (N, IUPAC codes,
gaps) maps to 4.  Sequences are returned as numpy uint8 arrays, one contig at
a time, so memory stays proportional to the largest contig.
"""

import gzip
import io

import numpy as np

_LUT = np.full(256, 4, dtype=np.uint8)
for _ch, _code in (("A", 0), ("C", 1), ("G", 2), ("T", 3), ("U", 3)):
    _LUT[ord(_ch)] = _code
    _LUT[ord(_ch.lower())] = _code

COMPLEMENT = np.array([3, 2, 1, 0, 4], dtype=np.uint8)


def encode(seq_bytes):
    """bytes/str -> uint8 array with A,C,G,T = 0..3 and N = 4."""
    if isinstance(seq_bytes, str):
        seq_bytes = seq_bytes.encode()
    return _LUT[np.frombuffer(seq_bytes, dtype=np.uint8)]


def reverse_complement(arr):
    return COMPLEMENT[arr[::-1]]


def _open(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def read_fasta(path):
    """Yield (name, uint8 array) for every record in a (gzipped) FASTA file."""
    name = None
    chunks = []
    with _open(path) as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None:
                    yield name, encode(b"".join(chunks))
                name = line[1:].split()[0].decode() if len(line) > 1 else ""
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            yield name, encode(b"".join(chunks))


def read_fasta_str(text):
    """Same as read_fasta but on an in-memory string (tests)."""
    name = None
    chunks = []
    for line in io.StringIO(text):
        if line.startswith(">"):
            if name is not None:
                yield name, encode("".join(chunks))
            name = line[1:].split()[0] if len(line) > 1 else ""
            chunks = []
        else:
            chunks.append(line.strip())
    if name is not None:
        yield name, encode("".join(chunks))


def write_fasta(path, records, width=80):
    """Write (name, str) records to a FASTA file."""
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")
