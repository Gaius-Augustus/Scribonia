import csv
import itertools

from benchmark.make_inputs import make_inputs
from benchmark.score import codetta_call, pf_call, read_codetta, score
from scribonia.fasta import read_fasta, write_fasta
from scribonia.synthetic import make_genome

HEADER = ("#accession\tspecies\ttaxid\tstop_set\ttable\tambiguous\tlabel_source\t"
          "citation\tverified\tnote\trole\n")


def _manifest(tmp_path):
    rows = [("ACC_T", "Tetra testis", "1", "TGA", "6", "no", "lit", "", "yes", "", "primary"),
            ("ACC_S", "Stand ardis", "2", "TAA,TAG,TGA", "1", "no", "lit", "", "yes", "", "primary")]
    for i, (acc, *_r) in enumerate(rows):
        d = tmp_path / "genomes" / acc
        d.mkdir(parents=True)
        stops = {"TGA"} if acc == "ACC_T" else {"TAA", "TAG", "TGA"}
        write_fasta(d / "g.fna", make_genome(stops, size_bp=200_000, n_contigs=4, seed=i))
    m = tmp_path / "manifest.tsv"
    m.write_text(HEADER + "".join("\t".join(r) + "\n" for r in rows))
    return m


def test_make_inputs_writes_genomes_and_bins(tmp_path):
    m = _manifest(tmp_path)
    rows = make_inputs(m, tmp_path / "genomes", tmp_path / "bench", species="Tetra",
                       n_bins=3, min_mb=0.05, max_mb=0.1, contig_range=(2_000, 20_000),
                       contamination_prob=1.0, seed=1)
    assert [r["kind"] for r in rows] == ["genome", "bin", "bin", "bin"]
    assert all(r["stop_set"] == "TGA" for r in rows)
    assert all(r["contaminant"] == "ACC_S" for r in rows[1:])
    for r in rows:
        lens = [len(a) for _n, a in read_fasta(r["fasta"])]
        assert sum(lens) == r["size_bp"]
    assert rows[0]["size_bp"] == sum(len(a) for _n, a in
                                     read_fasta(tmp_path / "genomes" / "ACC_T" / "g.fna"))
    with open(tmp_path / "bench" / "samples.tsv") as fh:
        assert len(list(csv.DictReader(fh, delimiter="\t"))) == 4
    # bench.sbatch reads the fasta column with cut: no CRLF line ends.
    assert b"\r" not in (tmp_path / "bench" / "samples.tsv").read_bytes()


def _codetta_out(path, inferred):
    """Minimal Codetta report: codon rows plus the final code line."""
    lines = ["# codon   inference  std code  diff?    N aligned  N used"]
    code = ""
    for c in ("".join(t) for t in itertools.product("TCAG", repeat=3)):
        inf = inferred.get(c, "?" if c in ("TAA", "TAG", "TGA") else "A")
        code += inf
        lines.append(f"{c:<10}{inf:<11}{'A':<10}{'N':<9}{100:<11}{50:<10}")
    lines += ["# codon      logP(A)", "TTT    " + " ".join(["-1.0"] * 21),
              "#", "# Final genetic code inference", code]
    path.write_text("\n".join(lines))


def test_codetta_stop_set(tmp_path):
    f = tmp_path / "x.genetic_code.out"
    _codetta_out(f, {"TAA": "Q", "TAG": "Q"})
    stop, evidence, detail = codetta_call(read_codetta(f))
    assert stop == frozenset({"TGA"})
    assert evidence == 64 * 50
    assert detail.startswith("TAA:Q/50")


def test_phylofisher_rule():
    counts = {"CAA": {"Q": 9000}, "TAA": {"Q": 40, "K": 2}, "TGA": {"W": 1}}
    stop, total, _ = pf_call(counts)
    assert total == 9043
    assert stop == frozenset({"TAG", "TGA"})     # one TGA hit is noise
    assert pf_call(counts, min_n=100)[0] == frozenset({"TAA", "TAG", "TGA"})


def test_score_end_to_end(tmp_path):
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "sample\tkind\taccession\tspecies\tstop_set\tambiguous\tsize_bp\t"
        "contamination\tcontaminant\tfasta\n"
        "G_T\tgenome\tACC_T\tTetra testis\tTGA\tno\t100\t0\t\tx\n"
        "G_S\tgenome\tACC_S\tStand ardis\tTAA,TAG,TGA\tno\t100\t0\t\tx\n"
        "G_U\tgenome\tACC_U\tUnfin ished\tTAA,TAG,TGA\tno\t100\t0\t\tx\n")
    res = tmp_path / "results"
    for sample, call in (("G_T", "TGA"), ("G_S", "TAA,TAG,TGA"), ("G_U", "TAA,TAG,TGA")):
        d = res / "scribonia" / sample
        d.mkdir(parents=True)
        (d / "out.tsv").write_text(f"sample\tstop_set\tn_orf_eval\n\t{call}\t10\n")
        (d / "time.tsv").write_text("2.0 1.5 0.5 1000\n")
        d = res / "codetta" / sample
        d.mkdir(parents=True)
        if sample != "G_U":
            _codetta_out(d / "out.genetic_code.out", {"TAA": "Q", "TAG": "Q"})   # always {TGA}
        else:                    # still running: GNU time made an empty time.tsv
            (d / "time.tsv").write_text("")
        d = res / "phylofisher" / sample
        d.mkdir(parents=True)
        if sample == "G_T":
            (d / "counts.tsv").write_text("codon\taa\tn\nCAA\tQ\t900\nTAA\tQ\t30\nTAG\tQ\t20\n")
        elif sample == "G_S":    # ended without counts: failed
            (d / "time.tsv").write_text("Command exited with non-zero status 1\n9.0 1.0 0.1 50\n")
        else:
            (d / "counts.tsv").write_text("codon\taa\tn\nCAA\tQ\t900\n")
    per, summary, sweep = score(samples, res, tmp_path / "report")
    s = {(r["subset"], r["tool"]): r for r in summary}
    assert s["all", "scribonia"]["exact"] == 1.0 and s["all", "scribonia"]["cpu_s_median"] == 2.0
    assert s["all", "codetta"]["exact"] == 0.5 and s["all", "codetta"]["pending"] == 1
    assert s["all", "codetta"]["std_called_nonstd"] == 1.0
    assert s["all", "phylofisher"]["failed"] == 1 and s["all", "phylofisher"]["exact"] == 1.0
    # Only G_T is finished by all three tools.
    assert all(s["common", t]["n"] == 1 for t in ("scribonia", "codetta", "phylofisher"))
    assert s["common", "codetta"]["exact"] == 1.0
    assert b"\r" not in (tmp_path / "report" / "per_sample.tsv").read_bytes()
    assert sweep and all(r["kind"] == "genome" for r in sweep)


def test_compare_stopcr_joins_by_accession_and_classes(tmp_path):
    from benchmark.compare_stopcr import acc_key, compare
    assert acc_key("GCF_000165425.1") == acc_key("GCA_000165425.2") == "000165425"
    assert acc_key("SRR1296700") is None
    ref = tmp_path / "stopcr.tsv"
    ref.write_text(
        "# comment\n"
        "species\taccession\tfig4c_TAA\tfig4c_TAG\tfig4c_TGA\tterm_pct_TAA\tterm_pct_TAG\t"
        "term_pct_TGA\tcds_n_TAA\tcds_n_TAG\tcds_n_TGA\tcds_n_source\tnote\n"
        "Stich\tGCA_000000001.2\tQ\tQ\tR\t0.0\t1.8\t98.2\t12512\t6848\t116\tFig. S5\t\n"
        "Trans\tSRR1\t*\t*\t*\t50\t40\t10\t\t\t\t\t\n")
    samples = tmp_path / "samples.tsv"
    samples.write_text("sample\tkind\taccession\tspecies\tstop_set\n"
                       "G\tgenome\tGCA_000000001.1\tStich\tTGA\n"
                       "B1\tbin\tGCA_000000001.1\tStich\tTGA\n"
                       "B2\tbin\tGCA_000000001.1\tStich\tTGA\n"
                       "X\tgenome\tGCA_999.1\tOther\tTGA\n")
    res = tmp_path / "res"
    head = "stop_set\tp_taa_stop\tp_tag_stop\tp_tga_stop\td_taa\td_tag\td_tga\n"
    for sample, line in (("G", "TGA\t0.01\t0.02\t0.99\t1.0\t1.1\tnan"),
                         ("B1", "TGA\t0.1\t0.1\t0.9\t1.0\t1.0\tnan"),
                         ("B2", "TAG,TGA\t0.1\t0.6\t0.9\t1.0\t0.3\tnan")):
        d = res / "scribonia" / sample
        d.mkdir(parents=True)
        (d / "out.tsv").write_text(head + line + "\n")
    per, summary = compare(samples, res, tmp_path / "out", ref)
    assert len(per) == 6                                      # 1 accession x 2 kinds x 3 codons
    r = {(p["kind"], p["codon"]): p for p in per}
    assert r["genome", "TGA"]["class"] == "bifunctional" and r["genome", "TAA"]["class"] == "sense"
    assert r["genome", "TGA"]["scribonia_stop_frac"] == 1.0 and r["genome", "TGA"]["cds_n"] == "116"
    assert r["bin", "TAG"]["scribonia_stop_frac"] == 0.5 and r["bin", "TAG"]["n_samples"] == 2
    assert r["bin", "TGA"]["d_median"] == "NA" and r["genome", "TGA"]["codetta"] == "NA"
    s = {(x["kind"], x["class"]): x for x in summary}
    assert s["genome", "stop"]["n_codons"] == 0
    assert s["genome", "all"]["spearman_p_stop_vs_term"] == 1.0
