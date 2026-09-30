"""Rule-based stop-set classifier: thresholds on the statistics in features.py.

The fallback when the trained model cannot be loaded (`--rules`, missing
model file, no scikit-learn) and a cross-check of the trained model.  The
thresholds were set on synthetic genomes (see tests/) and are deliberately
conservative: when the evidence is thin the probabilities stay near 0.5 and
the confidence is low, so that a downstream pipeline falls back to its
taxonomy prior.

Evidence for a codon c, in order of preference:

  1. the two-codon hypothesis without c, when it has an excess of long runs
     (its runs have the shortest non-coding flanks, so D_c is cleanest);
  2. otherwise the better of the two single-codon hypotheses without c;
  3. otherwise (every hypothesis without c chops the genes, which happens
     when both other codons are sense codons) the coverage ratio of c alone:
     a stop codon leaves long c-free runs over every gene, a sense codon
     does not.

Scoring lives in features.codon_scores, shared with the per-contig votes.
"""

from .features import CODONS, codon_scores, evidence_for as _evidence_for
from .tables import canonical_table, format_stop_set

RULES_VERSION = "rules-0.4"

# When contigs disagree and the voting contigs cover at least this share of
# the bin, the length-weighted majority of per-contig calls replaces the
# pooled call: pooled statistics follow whichever component supplies the
# long runs, the majority follows the dominant organism.
MAJORITY_MIN_VOTING_FRAC = 0.5


def evidence_for(summary, c):
    return _evidence_for(summary["hyp"], c)


def codon_probabilities(summary):
    """Return {codon: (p_stop, evidence_runs, note)} from a feature summary."""
    return codon_scores(summary["hyp"], summary.get("covratio", {}))


def confidence_of(probs, n_orf_eval):
    margins = [max(p, 1 - p) for p, _, _ in probs.values()]
    if min(margins) >= 0.95 and n_orf_eval >= 200:
        return "high"
    if min(margins) >= 0.80 and n_orf_eval >= 50:
        return "medium"
    return "low"


def classify_summary(summary, contamination_thresholds=(0.15, 0.4)):
    """Rule-based call from a feature summary.  Returns a result dict."""
    probs = codon_probabilities(summary)
    # Evidence actually used per codon (the feature n_orf_eval only counts
    # hypotheses without the codon, which is 0 when a codon is judged by its
    # own runs).
    n_orf_eval = min(n for _, n, _ in probs.values())
    stop_set = frozenset(c for c, (p, _, _) in probs.items() if p >= 0.5)
    ambiguous = [c for c, (p, n, _) in probs.items()
                 if 0.25 <= p <= 0.75 and n >= 200]
    conf = confidence_of(probs, n_orf_eval)
    if summary.get("size_bp", 0) < 500_000 or n_orf_eval < 50:
        conf = "low"

    minority, iqr_max = 0.0, 0.0
    for c in CODONS:
        v = summary["votes"][c]
        if v["n"] >= 2:
            minority = max(minority, min(v["stop_frac"], v["sense_frac"]))
            iqr_max = max(iqr_max, v["iqr"])
    contamination = (minority > contamination_thresholds[0]
                     or iqr_max > contamination_thresholds[1])

    majority = summary.get("majority")
    decision = "pooled"
    if (contamination and majority
            and majority["voting_bp_frac"] >= MAJORITY_MIN_VOTING_FRAC
            and majority["stop_set"] != "NA"):
        maj_set = frozenset(majority["stop_set"].split(","))
        if maj_set != stop_set:
            stop_set = maj_set
            decision = f"majority({majority['weight_frac']:.2f})"
            # Probabilities follow the majority call so that the table
            # stays consistent with stop_set.
            probs = {c: (0.75 if c in maj_set else 0.25, probs[c][1],
                         probs[c][2] + ";majority") for c in CODONS}
            conf = "medium" if conf == "high" else conf

    return {
        "decision": decision,
        "p_stop": {c: probs[c][0] for c in CODONS},
        "evidence": {c: probs[c][1] for c in CODONS},
        "notes": {c: probs[c][2] for c in CODONS},
        "stop_set": stop_set,
        "stop_set_str": format_stop_set(stop_set),
        "table": canonical_table(stop_set),
        "confidence": conf,
        "ambiguous": ambiguous,
        "contamination_flag": contamination,
        "minority_vote_frac": minority,
        "classifier_version": RULES_VERSION,
    }
