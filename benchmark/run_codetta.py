#!/usr/bin/env python3
"""Run Codetta 2.0 on one FASTA with its hmmscan chunks in parallel.

Codetta runs its hmmscan chunk scripts one after another, or submits them
as its own SLURM array (--parallelize_hmmscan s), which does not fit a
benchmark job that must account for its own CPU time.  This shim calls the
same three steps as codetta.py (align, summary, infer; default parameters:
Pfam-A_enone.hmm, e-value 1e-10, probability 0.9999, max fraction 0.01,
all problematic Pfam groups excluded) but lets processing_genome() only
write the chunk scripts, then runs them --threads at a time with one
hmmscan thread each.  The inference is unchanged.

Output: <prefix>.genetic_code.out (Codetta's inference report); the
alignment temp files are removed unless --keep.

Usage: python3 run_codetta.py --codetta-dir DIR --threads N FASTA PREFIX
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("fasta")
    p.add_argument("prefix", help="output prefix (path allowed)")
    p.add_argument("--codetta-dir", required=True)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--keep", action="store_true", help="keep alignment temp files")
    a = p.parse_args(argv)

    sys.path.insert(0, a.codetta_dir)
    import codetta   # noqa: E402  (needs helper_functions from codetta_dir)

    resources = os.path.join(a.codetta_dir, "resources")
    ns = argparse.Namespace(
        sequence_file=a.fasta, align_output=a.prefix,
        inference_output=a.prefix + ".genetic_code.out",
        profiles="Pfam-A_enone.hmm", resource_directory=resources, bad_profiles=None,
        results_summary=None, identifier=None, download_type=None,
        # Any value other than None or 's': write the chunk scripts, run nothing.
        parallelize_hmmscan="local",
        evalue=1e-10, probability_threshold=0.9999, max_fraction=0.01,
        mito_pfams=False, transposon_pfams=False, viral_pfams=False,
        selenocysteine_pfams=False, pyrrolysine_pfams=False)
    codetta.initialize_globals()
    codetta.initialize_emissions_dict(resources, ns.profiles, None)
    gc = codetta.GeneticCode(ns)
    gc.processing_genome()

    scripts = sorted(glob.glob(os.path.join(gc.scratch_dir, "hmmscan_*.sh")))
    env = dict(os.environ, HMMER_NCPU="1")

    def run(script):
        return subprocess.call(["bash", script], env=env)

    print(f"[codetta] {len(scripts)} hmmscan chunks on {a.threads} threads", flush=True)
    with ThreadPoolExecutor(max_workers=max(a.threads, 1)) as ex:
        codes = list(ex.map(run, scripts))
    if any(codes):
        sys.exit(f"[codetta] {sum(1 for c in codes if c)} hmmscan chunks failed")

    gc.process_hmmscan_results()
    gc.compute_decoding_probabilities()
    if not a.keep:
        shutil.rmtree(gc.scratch_dir, ignore_errors=True)
        for ext in (".preliminary_translation.faa", ".preliminary_translation.faa.ssi",
                    ".sequence_pieces.fna"):
            if os.path.exists(a.prefix + ext):
                os.remove(a.prefix + ext)


if __name__ == "__main__":
    main()
