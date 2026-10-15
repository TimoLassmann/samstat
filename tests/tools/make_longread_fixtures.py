#!/usr/bin/env python3
"""Generate the SAMStat long-read test fixtures and their ground truth.

Standard library only. Fully deterministic: a private SplitMix64 generator is
used instead of `random`, so the same --seed gives the same bytes on every
Python version. Run it once with `uv run python`; the outputs are committed, so
CI and ctest never run Python.

Two modes:

  make_longread_fixtures.py [--outdir DIR] [--seed N] [--test-view PATH]
      Write the committed fixtures and README.md (the expected values) into
      DIR (default: tests/data/longread next to this script's repository).
      BAM files are made with htslib's test/test_view, as tests/data/README.md
      describes. test_view is looked up in --test-view, $TEST_VIEW and PATH; if
      none is found, `samtools view -b --no-PG` is used instead (the BAM bytes
      may then differ from the committed ones).

  make_longread_fixtures.py --reads N [--profile ont|short] [--seed N] [--out F] [--force]
      Write N reads as SAM (stdout, or F) for tests/perf/bench.sh. A pool of
      simulated templates is cycled with fresh names and positions, so a 1 Gb
      file takes seconds. Not byte-stable across --seed values. N is capped at
      200000 (about 75 MB short, 3.5 GB ONT-like) unless --force is given, so
      that a typo cannot fill the disk.

Ground truth: every expected value in README.md is computed from the
simulation (the generator knows the true alignment of each read) and from the
definitions in docs/plan-v2.4/standard-plots.md section 3. Nothing is ever
taken from a samstat run.
"""

import argparse
import math
import os
import shutil
import statistics
import subprocess
import sys
import tempfile
from bisect import bisect_right

MASK = (1 << 64) - 1
ACGT = "ACGT"

# Storage grid of standard-plots.md section 3.1: exact lengths below 1000,
# then g(1000) (about +1% per step).
_G = [1000]
while _G[-1] < 10_000_000:
    _G.append(max(_G[-1] + 1, (_G[-1] * 101 + 99) // 100))

BASE_TABLE = bytes(ord(ACGT[b & 3]) for b in range(256))


class Rng:
    """SplitMix64. Stable by construction (integer arithmetic only)."""

    def __init__(self, seed):
        self.s = seed & MASK

    def u64(self):
        self.s = (self.s + 0x9E3779B97F4A7C15) & MASK
        z = self.s
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK
        return z ^ (z >> 31)

    def random(self):
        return (self.u64() >> 11) * (1.0 / 9007199254740992.0)

    def below(self, n):
        return self.u64() % n

    def between(self, a, b):
        return a + self.below(b - a + 1)

    def choice(self, seq):
        return seq[self.below(len(seq))]

    def shuffle(self, lst):
        for i in range(len(lst) - 1, 0, -1):
            j = self.below(i + 1)
            lst[i], lst[j] = lst[j], lst[i]

    def expo(self, rate):
        return -math.log(1.0 - self.random()) / rate

    def raw(self, n):
        out = bytearray()
        while len(out) < n:
            out += self.u64().to_bytes(8, "little")
        return bytes(out[:n])

    def bases(self, n):
        if n <= 0:
            return ""
        return self.raw(n).translate(BASE_TABLE).decode()

    def other_base(self, b):
        return ACGT[(ACGT.index(b) + 1 + self.below(3)) & 3]


_NORMAL = statistics.NormalDist()
_QUANTILES = [_NORMAL.inv_cdf((i + 0.5) / 256.0) for i in range(256)]


def qual_string(rng, n, mu, sd, lo, hi):
    """n Phred+33 characters, Gaussian(mu, sd) clamped to [lo, hi]."""
    if n <= 0:
        return ""
    table = [chr(min(hi, max(lo, int(round(mu + sd * z)))) + 33) for z in _QUANTILES]
    return "".join(table[b] for b in rng.raw(n))


# ---------------------------------------------------------------------------
# Alignment simulation
# ---------------------------------------------------------------------------

class Aln:
    """An alignment of a query against a random reference, with exact truth."""

    def __init__(self):
        self.segs = []      # ('=', n) ('X', 1) ('I', n) ('D', n)
        self.qseq = ""
        self.md = ""
        self.M = 0          # bases in = and X (and M) operations
        self.X = 0
        self.I_ev = self.I_b = self.D_ev = self.D_b = 0

    @property
    def nm(self):
        return self.X + self.I_b + self.D_b

    @property
    def ref_span(self):
        return self.M + self.D_b

    def cigar_m(self):
        out = []
        for op, n in self.segs:
            op = "M" if op in "=X" else op
            if out and out[-1][0] == op:
                out[-1][1] += n
            else:
                out.append([op, n])
        return "".join("%d%s" % (n, op) for op, n in out)

    def cigar_eqx(self):
        out = []
        for op, n in self.segs:
            if out and out[-1][0] == op:
                out[-1][1] += n
            else:
                out.append([op, n])
        return "".join("%d%s" % (n, op) for op, n in out)


def _event_len(rng):
    k = 1
    while k < 8 and rng.random() < 0.35:
        k += 1
    return k


def build_alignment(rng, n_query, p_sub, p_ins, p_del):
    """Simulate an alignment consuming exactly n_query query bases.

    p_* are event probabilities per aligned base. Every event is preceded by at
    least one matching base and the alignment ends in a match, so MD is always
    well formed (a number between any two letters) and no two indels touch.
    """
    a = Aln()
    rate = p_sub + p_ins + p_del
    read = []
    md = []
    run = 0
    q = 0
    while True:
        remain = n_query - q
        if rate <= 0:
            gap = remain
        else:
            gap = 1 + int(rng.expo(rate))
        pick = rng.random() * rate if rate > 0 else 0.0
        if pick < p_sub:
            kind, k = "X", 1
        elif pick < p_sub + p_ins:
            kind, k = "I", _event_len(rng)
        else:
            kind, k = "D", _event_len(rng)
        consumed = k if kind in "XI" else 0
        if rate <= 0 or q + gap + consumed + 1 > n_query:
            gap = remain
            if gap > 0:
                read.append(rng.bases(gap))
                a.segs.append(("=", gap))
                a.M += gap
                run += gap
            break
        read.append(rng.bases(gap))
        a.segs.append(("=", gap))
        a.M += gap
        run += gap
        q += gap
        if kind == "X":
            ref = rng.choice(ACGT)
            read.append(rng.other_base(ref))
            a.segs.append(("X", 1))
            a.M += 1
            a.X += 1
            md.append("%d%s" % (run, ref))
            run = 0
            q += 1
        elif kind == "I":
            read.append(rng.bases(k))
            a.segs.append(("I", k))
            a.I_ev += 1
            a.I_b += k
            q += k
        else:
            md.append("%d^%s" % (run, rng.bases(k)))
            run = 0
            a.segs.append(("D", k))
            a.D_ev += 1
            a.D_b += k
    md.append(str(run))
    a.qseq = "".join(read)
    a.md = "".join(md)
    assert len(a.qseq) == n_query, (len(a.qseq), n_query)
    return a


# ---------------------------------------------------------------------------
# Records
# ---------------------------------------------------------------------------

class Rec:
    def __init__(self, qname, flag, rname="*", pos=0, mapq=0, cigar="*",
                 rnext="*", pnext=0, tlen=0, seq="*", qual="*", tags=(), tr=None):
        self.qname, self.flag, self.rname, self.pos = qname, flag, rname, pos
        self.mapq, self.cigar, self.rnext, self.pnext = mapq, cigar, rnext, pnext
        self.tlen, self.seq, self.qual, self.tags = tlen, seq, qual, list(tags)
        self.tr = tr or {"mapped": False}

    def line(self):
        return "\t".join([self.qname, str(self.flag), self.rname, str(self.pos),
                          str(self.mapq), self.cigar, self.rnext, str(self.pnext),
                          str(self.tlen), self.seq, self.qual] + self.tags)


def unmapped_rec(qname, seq, qual, flag=4, tags=()):
    return Rec(qname, flag, seq=seq, qual=qual, tags=tags)


def mapped_tr(aln, soft, has_nm, has_md, form):
    return {"mapped": True, "M": aln.M, "X": aln.X, "I_ev": aln.I_ev,
            "I_b": aln.I_b, "D_ev": aln.D_ev, "D_b": aln.D_b, "S": soft,
            "has_nm": has_nm, "has_md": has_md, "eqx": form == "eqx"}


def aligned_rec(rng, qname, length, clip5, clip3, rates, flag, rname, pos, mapq,
                qual, form="M", nm=True, md=True, extra_tags=(), mate=None):
    """A mapped record built from a simulated alignment; returns (Rec, Aln)."""
    aln = build_alignment(rng, length - clip5 - clip3, *rates)
    if form == "eqx":
        assert aln.X > 0, "an =/X read must contain at least one X"
    cig = aln.cigar_eqx() if form == "eqx" else aln.cigar_m()
    cigar = ("%dS" % clip5 if clip5 else "") + cig + ("%dS" % clip3 if clip3 else "")
    seq = rng.bases(clip5) + aln.qseq + rng.bases(clip3)
    assert len(seq) == length
    tags = []
    if nm:
        tags.append("NM:i:%d" % aln.nm)
    if md:
        tags.append("MD:Z:%s" % aln.md)
    tags.extend(extra_tags)
    rec = Rec(qname, flag, rname, pos, mapq, cigar, seq=seq, qual=qual, tags=tags,
              tr=mapped_tr(aln, clip5 + clip3, nm, md, form))
    if mate:
        rec.rnext, rec.pnext, rec.tlen = mate
    return rec, aln


def perfect_rec(qname, flag, rname, pos, mapq, seq, cigar, qual, *, mismatch_md=None,
                nm=0, soft=0, extra_tags=(), has_md=True, has_nm=True, M=None, X=0,
                I_ev=0, I_b=0, D_ev=0, D_b=0):
    """Hand-built record; the caller states the truth explicitly."""
    tags = []
    if has_nm:
        tags.append("NM:i:%d" % nm)
    if has_md:
        tags.append("MD:Z:%s" % mismatch_md)
    tags.extend(extra_tags)
    tr = {"mapped": True, "M": M, "X": X, "I_ev": I_ev, "I_b": I_b, "D_ev": D_ev,
          "D_b": D_b, "S": soft, "has_nm": has_nm, "has_md": has_md, "eqx": False}
    return Rec(qname, flag, rname, pos, mapq, cigar, seq=seq, qual=qual, tags=tags, tr=tr)


def header(sq=(), hd="@HD\tVN:1.6", extra=()):
    lines = [hd] if hd else []
    lines += ["@SQ\tSN:%s\tLN:%d" % (n, ln) for n, ln in sq]
    lines += list(extra)
    return lines


# ---------------------------------------------------------------------------
# Truth (definitions: standard-plots.md section 3)
# ---------------------------------------------------------------------------

def edge_of(length):
    if length < 1000:
        return length
    return _G[bisect_right(_G, length) - 1]


def mapq_class(q):
    if q >= 30:
        return "mapq_ge30"
    if q >= 10:
        return "mapq_10_29"
    return "mapq_lt10"


def truth(recs, has_sq):
    prim = [r for r in recs if not r.flag & 0x900]
    t = {}
    t["total_reads"] = len(prim)
    t["read2_reads"] = sum(1 for r in prim if r.flag & 1 and r.flag & 0x80)
    t["read1_reads"] = t["total_reads"] - t["read2_reads"]
    t["paired"] = sum(1 for r in prim if r.flag & 1)
    t["properly_paired"] = sum(1 for r in prim if r.flag & 2)
    t["secondary"] = sum(1 for r in recs if r.flag & 0x100)
    t["supplementary"] = sum(1 for r in recs if r.flag & 0x800 and not r.flag & 0x100)
    lens = [0 if r.seq == "*" else len(r.seq) for r in prim]
    total = sum(lens)
    t["total_bases"] = total
    t["min_length"] = min(lens) if lens else 0
    t["max_length"] = max(lens) if lens else 0
    cum = 0
    n50 = 0
    for n in sorted(lens, reverse=True):
        cum += n
        if cum * 2 >= total:
            n50 = n
            break
    t["n50_exact"] = n50
    t["n50"] = edge_of(n50)
    asc = sorted(lens)
    med = asc[(len(asc) + 1) // 2 - 1] if asc else 0   # lower median
    t["median_length_exact"] = med
    t["median_length"] = edge_of(med)
    t["layout"] = "long" if t["median_length"] >= 1000 else "short"
    qsum = 0.0
    qn = 0
    phred_sum = 0
    maxq = 0
    for r in prim:
        if r.qual != "*":
            for ch in r.qual:
                q = ord(ch) - 33
                qsum += 10.0 ** (-q / 10.0)
                phred_sum += q
                maxq = max(maxq, q)
            qn += len(r.qual)
    t["qual_bases"] = qn
    t["mean_q"] = -10.0 * math.log10(qsum / qn) if qn else None
    t["mean_phred_arithmetic"] = phred_sum / qn if qn else None
    t["max_q"] = maxq if qn else None
    mapped = [r for r in prim if r.tr["mapped"]]
    cnt = {"mapq_ge30": 0, "mapq_10_29": 0, "mapq_lt10": 0, "unmapped": 0}
    for r in prim:
        cnt["unmapped" if not r.tr["mapped"] else mapq_class(r.mapq)] += 1
    t.update(cnt)
    t["is_aligned"] = 1 if (has_sq or mapped) else 0
    n = len(prim)
    t["pct_mapped"] = 100.0 * len(mapped) / n if n else None
    t["pct_mapq_ge30"] = 100.0 * cnt["mapq_ge30"] / n if n else None
    mb = sum(len(r.seq) for r in mapped)
    t["pct_clipped"] = 100.0 * sum(r.tr["S"] for r in mapped) / mb if mb else None
    bins = []
    exact = []
    na = 0
    wo_md = 0
    for r in mapped:
        x = r.tr
        if not x["has_md"] and not (x["eqx"] and x["X"] > 0):
            wo_md += 1
        if x["has_nm"] or x["eqx"]:
            den = x["M"] + x["I_ev"] + x["D_ev"]
            num = x["M"] - x["X"]
            bins.append((1000 * num) // den)
            exact.append(100.0 * num / den)
        else:
            na += 1
    t["identity_reads"] = len(bins)
    t["identity_na"] = na
    t["reads_without_md"] = wo_md
    sb = sorted(bins)
    t["median_identity"] = sb[(len(sb) + 1) // 2 - 1] / 10.0 if sb else None
    t["median_identity_exact"] = (statistics.median_low(exact) if exact else None)
    md_reads = [r.tr for r in mapped if r.tr["has_md"] or (r.tr["eqx"] and r.tr["X"] > 0)]
    t["mismatch_bases"] = sum(x["X"] for x in md_reads)
    t["insertion_events"] = sum(r.tr["I_ev"] for r in mapped)
    t["deletion_events"] = sum(r.tr["D_ev"] for r in mapped)
    t["aligned_bases"] = sum(r.tr["M"] + r.tr["I_b"] for r in mapped)
    t["aligned_md_bases"] = sum(x["M"] + x["I_b"] for x in md_reads)
    t["clipped_bases"] = sum(r.tr["S"] for r in mapped)
    return t


# ---------------------------------------------------------------------------
# Fixture builders. Each returns (lines, records, has_sq).
# ---------------------------------------------------------------------------

SQ3 = [("chr1", 250000), ("chr2", 180000), ("chr3", 60000)]
ONT_RATES = (0.035, 0.020, 0.025)      # sub, ins, del events per aligned base: ~8% error
HIFI_RATES = (0.0012, 0.0006, 0.0006)
SHORT_RATES = (0.006, 0.0005, 0.0005)
CLIP_CHOICES = [0, 0, 0, 0, 7, 23, 61, 140]


def place(rng, sq, span):
    name, ln = rng.choice(sq)
    return name, rng.between(1, ln - span - 1)


def ont_reads(rng, lengths, mapqs, unmapped_idx, sq, prefix="ont"):
    """ONT-like records for the given read lengths (list). Returns records."""
    recs = []
    mi = 0
    for i, length in enumerate(lengths):
        name = "%s_%03d" % (prefix, i + 1)
        mu = 8 + rng.random() * 10
        qual = qual_string(rng, length, mu, 4.0, 2, 40)
        if i in unmapped_idx:
            recs.append(unmapped_rec(name, rng.bases(length), qual))
            continue
        c5 = rng.choice(CLIP_CHOICES) if length >= 500 else 0
        c3 = rng.choice(CLIP_CHOICES) if length >= 500 else 0
        flag = 16 if rng.below(2) else 0
        # the reference span is known only after simulation; place afterwards
        rec, aln = aligned_rec(rng, name, length, c5, c3, ONT_RATES, flag, "*", 0,
                               mapqs[mi], qual)
        mi += 1
        rname, pos = place(rng, sq, aln.ref_span)
        rec.rname, rec.pos = rname, pos
        recs.append(rec)
    return recs


def sort_recs(recs, sq):
    order = {n: i for i, (n, _) in enumerate(sq)}
    return sorted(recs, key=lambda r: (1, 0) if not r.tr["mapped"] else (0, order[r.rname] * 10**9 + r.pos))


def log_lengths(rng, n, lo, hi):
    out = [int(round(lo * (hi / lo) ** ((i + rng.random()) / n))) for i in range(n)]
    rng.shuffle(out)
    return out


ONT_MAPQ = [60] * 9 + [50, 44, 30, 29, 20, 15, 10, 9, 7, 5, 1, 0, 0]


def build_ont_like(seed):
    rng = Rng(seed)
    lengths = log_lengths(rng, 24, 500, 3500)
    unmapped_idx = {5, 17}
    mapqs = list(ONT_MAPQ)
    rng.shuffle(mapqs)
    recs = ont_reads(rng, lengths, mapqs, unmapped_idx, SQ3)
    recs = sort_recs(recs, SQ3)
    hdr = header(SQ3, "@HD\tVN:1.6\tSO:coordinate")
    # the extra read of the BAM: 70000 CIGAR operations (> 65535), so that BAM
    # stores the CIGAR in a CG:B:I tag.
    unit = 17500
    pat = rng.bases(1000)
    seq = (pat * 53)[:3 * unit]
    cigar = "1M1I1M1D" * unit
    md = "2^A" * unit + "0"
    aln = Aln()
    aln.M, aln.X = 2 * unit, 0
    aln.I_ev = aln.I_b = unit
    aln.D_ev = aln.D_b = unit
    long_rec = Rec("ont_longcigar", 0, "chr1", 1, 60, cigar, seq=seq, qual="5" * len(seq),
                   tags=["NM:i:%d" % aln.nm, "MD:Z:%s" % md], tr=mapped_tr(aln, 0, True, True, "M"))
    bam_recs = sort_recs(recs + [long_rec], SQ3)
    return hdr, recs, bam_recs


def build_supp_chain(seed):
    rng = Rng(seed)
    read = rng.bases(2000)
    q = qual_string(rng, 2000, 15, 3, 5, 30)
    recs = []
    sa = "SA:Z:chr2,5001,+,1000H600M400H,40,0;chr3,7001,+,1600H400M,35,0;"
    recs.append(perfect_rec("chain", 0, "chr1", 10001, 60, read, "1000M1000S", q,
                            mismatch_md="1000", M=1000, soft=1000, extra_tags=[sa]))
    recs.append(perfect_rec("chain", 2048, "chr2", 5001, 40, read[1000:1600], "1000H600M400H",
                            q[1000:1600], mismatch_md="600", M=600,
                            extra_tags=["SA:Z:chr1,10001,+,1000M1000S,60,0;"]))
    recs.append(perfect_rec("chain", 2048, "chr3", 7001, 35, read[1600:], "1600H400M",
                            q[1600:], mismatch_md="400", M=400,
                            extra_tags=["SA:Z:chr1,10001,+,1000M1000S,60,0;"]))
    recs.append(perfect_rec("chain", 256, "chr1", 40001, 0, read, "2000M", q,
                            mismatch_md="2000", M=2000))
    solo = rng.bases(1000)
    solo_md = "500%s499" % ("A" if solo[500] != "A" else "C")
    recs.append(perfect_rec("solo", 16, "chr1", 20001, 60, solo, "1000M",
                            qual_string(rng, 1000, 20, 3, 5, 30), mismatch_md=solo_md,
                            nm=1, M=1000, X=1))
    hdr = header([("chr1", 250000), ("chr2", 180000), ("chr3", 60000)])
    return hdr, recs


def make_mm_ml(rng, seq):
    skips = []
    count = 0
    for ch in seq:
        if ch == "C":
            if rng.random() < 0.08:
                skips.append(count)
                count = 0
            else:
                count += 1
    if not skips:
        skips = [0]
    mm = "MM:Z:C+m," + ",".join(str(s) for s in skips) + ";"
    ml = "ML:B:C," + ",".join(str(rng.below(256)) for _ in skips)
    return [mm, ml]


def build_hifi_like(seed):
    rng = Rng(seed)
    lengths = [3000, 4200, 4800, 6000]
    mapqs = [60, 60, 40, 20]
    recs = []
    for i, length in enumerate(lengths):
        mu = 46 + 4 * rng.random()
        qual = list(qual_string(rng, length, mu, 10.0, 12, 60))
        qual[100:140] = ["%c" % (60 + 33)] * 40      # a run of Q60 in every read
        qual = "".join(qual)
        # the sequence is needed before MM/ML, so build the alignment first
        aln = build_alignment(rng, length, *HIFI_RATES)
        seq = aln.qseq
        mm = make_mm_ml(rng, seq)
        cig = aln.cigar_m()
        rname, pos = place(rng, SQ3, aln.ref_span)
        recs.append(Rec("hifi_%03d" % (i + 1), 0, rname, pos, mapqs[i], cig, seq=seq,
                        qual=qual, tags=["NM:i:%d" % aln.nm, "MD:Z:%s" % aln.md] + mm,
                        tr=mapped_tr(aln, 0, True, True, "M")))
    recs = sort_recs(recs, SQ3)
    return header(SQ3, "@HD\tVN:1.6\tSO:coordinate"), recs


def build_long_fq(seed):
    """FASTQ records as (name, seq, qual). Last record is not newline-terminated."""
    rng = Rng(seed)
    lengths = [1, 40, 250, 800, 1000, 1200, 1500, 2000, 2000, 2600, 3000, 5000]
    recs = []
    for i, n in enumerate(lengths):
        name = "lr_%02d" % (i + 1)
        seq = list(rng.bases(n))
        qual = qual_string(rng, n, 10 + rng.random() * 8, 4.0, 2, 40)
        recs.append([name, "".join(seq), qual])
    # reads 8 and 9 (2000 bp each): identical Phred mean (20), different
    # error-probability mean.
    recs[7][0] = "q20_flat"
    recs[7][2] = chr(20 + 33) * 2000
    recs[8][0] = "q20_split"
    recs[8][2] = (chr(10 + 33) * 100 + chr(30 + 33) * 100) * 10
    # a few N bases in read 7
    s = recs[6][1]
    recs[6][1] = s[:300] + "NNN" + s[303:]
    return recs


def fq_truth(recs):
    fake = []
    for name, seq, qual in recs:
        fake.append(Rec(name, 4, seq=seq, qual=qual))
    return truth(fake, False)


def build_eqx_nomd(seed):
    rng = Rng(seed)
    lengths = [300, 420, 510, 600, 700]
    mapqs = [60, 60, 45, 12, 60]
    recs = []
    for i, n in enumerate(lengths):
        has_nm = i < 3
        c5 = rng.choice([0, 9, 21])
        c3 = rng.choice([0, 0, 15])
        rec, aln = aligned_rec(rng, "eqx_%d" % (i + 1), n, c5, c3, (0.03, 0.01, 0.01), 0, "*", 0,
                               mapqs[i], qual_string(rng, n, 15, 3, 5, 30), form="eqx",
                               nm=has_nm, md=False)
        rec.rname, rec.pos = place(rng, SQ3, aln.ref_span)
        recs.append(rec)
    return header(SQ3), recs


def build_nonm(seed):
    rng = Rng(seed)
    lengths = [300, 380, 450, 520, 650]
    mapqs = [60, 33, 8, 60, 60]
    recs = []
    for i, n in enumerate(lengths):
        rec, aln = aligned_rec(rng, "nonm_%d" % (i + 1), n, 0, 0, (0.02, 0.01, 0.01), 0, "*", 0,
                               mapqs[i], qual_string(rng, n, 16, 3, 5, 30), nm=False, md=True)
        rec.rname, rec.pos = place(rng, SQ3, aln.ref_span)
        recs.append(rec)
    return header(SQ3), recs


def build_ubam_like(seed):
    rng = Rng(seed)
    lengths = [300, 450, 600, 750, 900, 1100]
    recs = []
    for i, n in enumerate(lengths):
        qual = qual_string(rng, n, 14, 4.0, 2, 40)
        recs.append(unmapped_rec("ubam_%d" % (i + 1), rng.bases(n), qual,
                                 tags=["RG:Z:run1", "qs:f:%.1f" % (14 + i * 0.5), "ns:i:%d" % (n * 9)]))
    hdr = ["@HD\tVN:1.6\tSO:unsorted", "@RG\tID:run1\tPL:ONT\tSM:sample",
           "@PG\tID:basecaller\tPN:synthetic"]
    return hdr, recs


def build_sq_nomap(seed):
    rng = Rng(seed)
    lengths = [300, 420, 550, 700, 850, 1000]
    recs = []
    for i, n in enumerate(lengths):
        recs.append(unmapped_rec("nomap_%d" % (i + 1), rng.bases(n),
                                 qual_string(rng, n, 14, 4.0, 2, 40)))
    return header(SQ3, "@HD\tVN:1.6\tSO:unsorted"), recs


def build_mapq255(seed):
    rng = Rng(seed)
    lengths = [320, 400, 470, 560, 640, 700, 760]
    mapqs = [255, 255, 255, 255, 255, 60, 5]
    recs = []
    for i, n in enumerate(lengths):
        rec, aln = aligned_rec(rng, "star_%d" % (i + 1), n, 0, 0, (0.02, 0.005, 0.005), 0, "*", 0,
                               mapqs[i], qual_string(rng, n, 30, 4.0, 10, 40))
        rec.rname, rec.pos = place(rng, SQ3, aln.ref_span)
        recs.append(rec)
    return header(SQ3), recs


def build_single_len(seed):
    rng = Rng(seed)
    recs = []
    mapqs = [60, 60, 60, 60, 60, 60, 25, 25, 3, 60]
    for i in range(10):
        rname, p1 = "chr1", 1000 + 400 * i
        p2 = p1 + 150
        name = "pair_%02d" % (i + 1)
        r1, a1 = aligned_rec(rng, name, 150, 0, 0, SHORT_RATES, 99, rname, p1, mapqs[i],
                             qual_string(rng, 150, 34, 3, 10, 40), mate=("=", p2, 300))
        r2, a2 = aligned_rec(rng, name, 150, 0, 0, SHORT_RATES, 147, rname, p2, mapqs[i],
                             qual_string(rng, 150, 30, 4, 10, 40), mate=("=", p1, -300))
        recs += [r1, r2]
    return header([("chr1", 20000)], "@HD\tVN:1.6\tSO:coordinate"), recs


def build_n50_five(seed):
    rng = Rng(seed)
    recs = []
    for i, n in enumerate([100, 200, 300, 400, 600]):
        seq = rng.bases(n)
        recs.append(perfect_rec("n50_%d" % n, 0, "chr1", 1 + 1000 * i, 60, seq, "%dM" % n,
                                "?" * n, mismatch_md=str(n), M=n))
    return header([("chr1", 20000)]), recs


def build_identity_indel(seed):
    rng = Rng(seed)
    recs = []
    # A: 50M1I49M, no mismatch.   NM = 1
    seq = rng.bases(100)
    recs.append(perfect_rec("id_ins1", 0, "chr1", 1001, 60, seq, "50M1I49M", "?" * 100,
                            mismatch_md="99", nm=1, M=99, I_ev=1, I_b=1))
    # B: 50M2D50M with one mismatch (ref A, read not A) at offset 30. NM = 1 + 2 = 3
    seq = rng.bases(100)
    seq = seq[:30] + ("C" if seq[30] == "A" else seq[30]) + seq[31:]
    recs.append(perfect_rec("id_del2_mm1", 0, "chr1", 2001, 60, seq, "50M2D50M", "?" * 100,
                            mismatch_md="30A19^GT50", nm=3, M=100, X=1, D_ev=1, D_b=2))
    # C: perfect 100M. NM = 0
    seq = rng.bases(100)
    recs.append(perfect_rec("id_perfect", 0, "chr1", 3001, 60, seq, "100M", "?" * 100,
                            mismatch_md="100", nm=0, M=100))
    # D: 40M3I57M: one insertion EVENT of 3 bases. NM = 3
    seq = rng.bases(100)
    recs.append(perfect_rec("id_ins3", 0, "chr1", 4001, 60, seq, "40M3I57M", "?" * 100,
                            mismatch_md="97", nm=3, M=97, I_ev=1, I_b=3))
    return header([("chr1", 20000)]), recs


def build_softclip_pct(seed):
    rng = Rng(seed)
    recs = []
    spec = [("clip_both", "10S80M10S", 10, 10, 80), ("clip_5p", "20S80M", 20, 0, 80),
            ("clip_3p", "80M20S", 0, 20, 80), ("clip_none", "100M", 0, 0, 100)]
    for i, (name, cig, c5, c3, m) in enumerate(spec):
        recs.append(perfect_rec(name, 0, "chr1", 1 + 500 * i, 60, rng.bases(100), cig, "?" * 100,
                                mismatch_md=str(m), M=m, soft=c5 + c3))
    recs.append(unmapped_rec("clip_unmapped", rng.bases(100), "?" * 100))
    return header([("chr1", 20000)]), recs


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------

def sam_text(hdr, recs):
    return "\n".join(hdr + [r.line() for r in recs]) + "\n"


def fq_text(recs):
    # last record deliberately has no trailing newline
    return "\n".join("@%s\n%s\n+\n%s" % (n, s, q) for n, s, q in recs)


def find_test_view(explicit):
    for cand in (explicit, os.environ.get("TEST_VIEW")):
        if cand:
            if not os.access(cand, os.X_OK):
                sys.exit("test_view not executable: %s" % cand)
            return [cand, "-b"]
    found = shutil.which("test_view")
    if found:
        return [found, "-b"]
    st = shutil.which("samtools")
    if st:
        sys.stderr.write("warning: test_view not found; using samtools view -b --no-PG\n")
        return [st, "view", "-b", "--no-PG"]
    sys.exit("neither test_view nor samtools found; pass --test-view PATH")


def make_bam(cmd, sam_path, bam_path):
    with open(bam_path, "wb") as out:
        subprocess.run(cmd + [sam_path], stdout=out, check=True)


# ---------------------------------------------------------------------------
# Bulk mode for tests/perf/bench.sh
# ---------------------------------------------------------------------------

def bulk(n_reads, profile, seed, out):
    rng = Rng(seed)
    bsq = [("chr1", 100_000_000), ("chr2", 80_000_000), ("chr3", 60_000_000)]
    if profile == "ont":
        pool_n = 256
        lengths = log_lengths(rng, pool_n, 500, 30000)
        mapqs = [rng.choice(ONT_MAPQ) for _ in range(pool_n)]
        pool = ont_reads(rng, lengths, mapqs, {i for i in range(pool_n) if i % 40 == 7}, bsq,
                         prefix="bench")
    else:
        pool_n = 2048
        pool = []
        for i in range(pool_n):
            mq = rng.choice([60, 60, 60, 60, 40, 22, 5, 0])
            flag = 16 if rng.below(2) else 0
            rec, aln = aligned_rec(rng, "bench", 150, 0, 0, SHORT_RATES, flag, "*", 0, mq,
                                   qual_string(rng, 150, 33, 4, 5, 40))
            rec.rname, rec.pos = place(rng, bsq, aln.ref_span)
            pool.append(rec)
    own = out is not sys.stdout.buffer
    fh = out
    fh.write(("\n".join(header(bsq, "@HD\tVN:1.6\tSO:unsorted")) + "\n").encode())
    lens = {n: ln for n, ln in bsq}
    batch = []
    for i in range(n_reads):
        t = pool[i % pool_n]
        f = t.line().split("\t")
        f[0] = "bench_%d" % (i + 1)
        if t.tr["mapped"]:
            span = t.tr["M"] + t.tr["D_b"]
            f[3] = str(1 + (i * 7919) % (lens[t.rname] - span - 1))
        batch.append("\t".join(f))
        if len(batch) >= 20000:
            fh.write(("\n".join(batch) + "\n").encode())
            batch = []
    if batch:
        fh.write(("\n".join(batch) + "\n").encode())
    if own:
        fh.close()


# ---------------------------------------------------------------------------
# README
# ---------------------------------------------------------------------------

KEYS_SUMMARY = ["total_reads", "read1_reads", "read2_reads", "paired", "properly_paired",
                "secondary", "supplementary", "is_aligned", "total_bases", "min_length",
                "max_length", "n50_exact", "n50", "median_length_exact", "median_length",
                "layout", "qual_bases", "mean_q", "mean_phred_arithmetic", "max_q",
                "mapq_ge30", "mapq_10_29", "mapq_lt10", "unmapped", "pct_mapped",
                "pct_mapq_ge30", "pct_clipped", "median_identity", "median_identity_exact",
                "identity_reads", "identity_na", "reads_without_md", "mismatch_bases",
                "insertion_events", "deletion_events", "aligned_bases", "aligned_md_bases",
                "clipped_bases"]


def fmt(v):
    if v is None:
        return "null (not applicable)"
    if isinstance(v, float):
        return "%.4f" % v
    return str(v)


def table(t, keys=None, skip_zero_paired=True):
    rows = ["| key | expected |", "|---|---|"]
    for k in (keys or KEYS_SUMMARY):
        if k not in t:
            continue
        rows.append("| `%s` | %s |" % (k, fmt(t[k])))
    return "\n".join(rows)


def build_readme(info):
    L = []
    a = L.append
    a("# SAMStat long-read fixtures")
    a("")
    a("Generated by `tests/tools/make_longread_fixtures.py` (seed %d); do not edit by hand. "
      "Regenerate with `uv run python tests/tools/make_longread_fixtures.py` from the repository root "
      "(BAM files need htslib `test/test_view`, see the script header)." % info["seed"])
    a("")
    a("Every expected value below comes from the generator's own ground truth (the true alignment "
      "of each simulated read) and from the definitions in `docs/plan-v2.4/standard-plots.md` "
      "section 3. None was read from a samstat run. Where a number can be checked by hand, the "
      "derivation is shown.")
    a("")
    a("## Conventions")
    a("")
    a("- Keys are the `summary` keys of `standard-plots.md` section 7. Under the transitional "
      "`v2.` prefix (waves 3-4) they live in `v2.summary.*`.")
    a("- Statistics are primary-only: secondary (0x100) and supplementary (0x800) records are "
      "excluded; a record with both flags counts as secondary.")
    a("- Read length is the full `l_qseq`, soft-clipped bases included. A record with `SEQ *` has length 0.")
    a("- `n50_exact` and `median_length_exact` are the exact values (median: lower median). `n50` and "
      "`median_length` are the lower edge of the storage length bin containing them "
      "(exact below 1000 bp, about +1% steps above; section 3.1). Compare samstat output "
      "to the binned value exactly, or to the exact value within 1%.")
    a("- `mean_q` is `-10 log10(mean of 10^(-q/10))` over all bases with a quality, clipped bases "
      "included. `mean_phred_arithmetic` is the plain mean of the Phred values, for contrast.")
    a("- `median_identity` is the lower edge (percent) of the identity bin (0.1% wide) holding the "
      "lower median read; `median_identity_exact` is the exact lower-median identity. Identity is "
      "`(M - X) / (M + I_events + D_events)` (section 3.4); samstat should agree within 0.5 "
      "percentage points.")
    a("- `mismatch_bases`, `insertion_events`, `deletion_events`, `aligned_bases`, `aligned_md_bases` "
      "are the sums of the `relative` position arrays `mismatches`, `insertions`, `deletions`, "
      "`aligned`, `aligned_md` over all MAPQ classes (the relative axis loses nothing; the absolute "
      "axis stops at position 1000). Position axes cover the aligned part of the read only.")
    a("- `clipped_bases` is the number of soft-clipped bases of mapped primary reads.")
    a("- SEQ and QUAL are written as stored in the file (reference orientation); flag 16 marks "
      "reverse-strand reads but SEQ is not reverse-complemented by the generator, so ground truth never "
      "depends on strand.")
    a("")
    a("## Size and scope notes")
    a("")
    a("The package asked for about 40 reads of 500 bp to 30 kb in `ont_like`, a 15 kb `hifi_like` "
      "and a 1-50 kb `long.fq`. SAM and FASTQ are uncompressed text (about 2.4 bytes per base for "
      "ONT-like reads) and each file is capped at 100 KB (300 KB in total), so the committed files "
      "are scaled down: `ont_like` has 24 reads of 500 bp to 3.5 kb, `hifi_like` 4 reads of 3.0 to "
      "6.0 kb, `long.fq` 12 reads of 1 bp to 5 kb. `ont_like.bam` additionally holds one 52.5 kb read. "
      "The 1 Mb case is built by CMake at test time (`tests/cmake/run_ultralong.cmake`), and "
      "`make_longread_fixtures.py --reads N` makes large files for `tests/perf/bench.sh`.")
    a("")
    for sec in info["sections"]:
        a("## %s" % sec["title"])
        a("")
        for p in sec["text"]:
            a(p)
            a("")
        if sec.get("truth"):
            a(table(sec["truth"], sec.get("keys")))
            a("")
        for p in sec.get("after", []):
            a(p)
            a("")
    a("## Checking the files")
    a("")
    a("`samtools quickcheck` and `samtools view` parse every SAM and BAM here. "
      "`samtools stats` and `samtools view -q` agree with the read, base, mapped, unmapped and MAPQ "
      "counts above for `ont_like.sam`, `supp_chain.sam` and `mapq255.sam` (a manual cross-check, not a "
      "test). Every MD tag, CIGAR, NM and SEQ of each SAM was also checked for mutual consistency "
      "(CIGAR query length equals SEQ length, MD matches and deletions equal the CIGAR, NM equals "
      "mismatches plus indel bases, each MD mismatch differs from the read base).")
    a("")
    return "\n".join(L)


def sec(title, text, t=None, keys=None, after=()):
    return {"title": title, "text": text, "truth": t, "keys": keys, "after": list(after)}


# ---------------------------------------------------------------------------

def default_outdir():
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.normpath(os.path.join(here, "..", "data", "longread"))


def generate(outdir, seed, test_view):
    os.makedirs(outdir, exist_ok=True)
    sections = []
    written = []

    def put(name, text):
        path = os.path.join(outdir, name)
        with open(path, "w", newline="\n") as fh:
            fh.write(text)
        written.append(path)

    cmd = find_test_view(test_view)

    # ---- ont_like -------------------------------------------------------
    hdr, recs, bam_recs = build_ont_like(seed + 1)
    put("ont_like.sam", sam_text(hdr, recs))
    with tempfile.TemporaryDirectory() as tmp:
        tmp_sam = os.path.join(tmp, "ont_like_bam.sam")
        with open(tmp_sam, "w", newline="\n") as fh:
            fh.write(sam_text(hdr, bam_recs))
        bam_path = os.path.join(outdir, "ont_like.bam")
        make_bam(cmd, tmp_sam, bam_path)
        written.append(bam_path)
    t_sam = truth(recs, True)
    t_bam = truth(bam_recs, True)
    mapq_list = sorted(r.mapq for r in recs if r.tr["mapped"])
    sections.append(sec(
        "ont_like.sam",
        ["24 reads simulated ONT-style: lengths log-uniform from 500 bp to 3.5 kb, about 8% error "
         "(3.5% substitutions, 2.0% insertion and 2.5% deletion events per aligned base, indel "
         "length 1-8), Q2-Q40 qualities with a per-read mean of Q8-Q18, 0-140 soft-clipped bases at "
         "either end, NM and MD on every mapped read, coordinate-sorted over three contigs "
         "(`@SQ` chr1 250000, chr2 180000, chr3 60000). Two reads are unmapped (flag 4). 22 reads are "
         "mapped with MAPQ values " + ", ".join(str(m) for m in mapq_list) + ".",
         "Median length is %d, so the layout is long (relative position axis)." % t_sam["median_length"]],
        t_sam))
    sections.append(sec(
        "ont_like.bam",
        ["The 24 reads of `ont_like.sam` plus one read, `ont_longcigar`: chr1:1, MAPQ 60, flag 0, "
         "CIGAR `1M1I1M1D` repeated 17500 times (70000 operations, more than the 65535 a BAM record "
         "can hold, so htslib writes the CIGAR to a `CG:B:I` tag), SEQ 52500 bp (a repeated random "
         "1 kb block), all qualities `5` (Q20), `NM:i:35000`, `MD:Z:2^A` repeated 17500 times then `0`. "
         "Its identity is (35000 - 0) / (35000 + 17500 + 17500) = 50.0%; its X count is "
         "NM - I_bases - D_bases = 35000 - 17500 - 17500 = 0.",
         "Hand check of the extra read: M = 2 x 17500 = 35000, I events = 17500, D events = 17500, "
         "aligned bases (M + I) = 52500."],
        t_bam))
    # ---- supp_chain -----------------------------------------------------
    hdr, recs = build_supp_chain(seed + 2)
    put("supp_chain.sam", sam_text(hdr, recs))
    t = truth(recs, True)
    sections.append(sec(
        "supp_chain.sam",
        ["Hand-computable. Read `chain` (2000 bp) has a primary alignment (flag 0, chr1:10001, "
         "`1000M1000S`, MAPQ 60), two hard-clipped supplementary alignments (flag 2048: chr2:5001 "
         "`1000H600M400H` MAPQ 40, chr3:7001 `1600H400M` MAPQ 35, SEQ only the aligned part) and one "
         "secondary alignment (flag 256, chr1:40001, `2000M`, MAPQ 0, full SEQ). Read `solo` (1000 bp, "
         "flag 16, `1000M`, `NM:i:1`, `MD:Z:500x499`) is the second primary.",
         "Only the two primaries count: 2 reads, 2000 + 1000 = 3000 bases, N50 = 2000 (2000 >= half of "
         "3000; the binned `n50` is 1993). Secondary = 1, supplementary = 2 (not double-counted in reads, bases or lengths). "
         "Clipped bases: 1000 of 3000 mapped bases = 33.3333%. "
         "Identity: `chain` 1000/1000 = 100.0%, `solo` (1000 - 1)/1000 = 99.9%; the lower median is 99.9."],
        t))
    # ---- hifi_like ------------------------------------------------------
    hdr, recs = build_hifi_like(seed + 3)
    put("hifi_like.sam", sam_text(hdr, recs))
    t = truth(recs, True)
    sections.append(sec(
        "hifi_like.sam",
        ["4 HiFi-style reads of 3000, 4200, 4800 and 6000 bp (the package asked for 15 kb; see the "
         "size note), about 0.2% error, per-base qualities Gaussian around Q46-Q50 (sd 10) clamped to "
         "Q12-Q60, plus a 40-base run of exactly Q60 at positions 100-139 of every read. MAPQ 60, 60, 40, "
         "20. Every read carries `MM:Z:C+m,...;` and `ML:B:C,...` tags (about 8% of C bases "
         "modified). Quality must not saturate at Q41: `max_q` is 60."],
        t))
    # ---- long.fq --------------------------------------------------------
    fq = build_long_fq(seed + 4)
    put("long.fq", fq_text(fq))
    t = fq_truth(fq)
    flat = [r for r in fq if r[0] == "q20_flat"][0]
    split = [r for r in fq if r[0] == "q20_split"][0]
    eflat = sum(10 ** (-(ord(c) - 33) / 10) for c in flat[2]) / len(flat[2])
    esplit = sum(10 ** (-(ord(c) - 33) / 10) for c in split[2]) / len(split[2])
    sections.append(sec(
        "long.fq",
        ["12 reads of 1, 40, 250, 800, 1000, 1200, 1500, 2000, 2000, 2600, 3000 and 5000 bp "
         "(the package asked for 1-50 kb; see the size note). The last record has no trailing newline. "
         "`lr_07` (1500 bp) contains 3 N bases. FASTQ has no alignment: unaligned layout, every read "
         "unmapped.",
         "Reads `q20_flat` and `q20_split` (2000 bp each) have the same arithmetic Phred mean, Q20 "
         "(flat: all Q20; split: alternating blocks of 100 bases at Q10 and Q30), but different "
         "mean error probability: flat %.4f, split %.4f (hand check: (0.1 + 0.001)/2 = 0.0505). "
         "A mean taken over Phred values would give the same value for both; a mean over error "
         "probabilities does not. Over the whole file the arithmetic Phred mean is %.3f and `mean_q` is %.3f."
         % (eflat, esplit, t["mean_phred_arithmetic"], t["mean_q"])],
        t, keys=["total_reads", "total_bases", "min_length", "max_length", "n50_exact", "n50",
                 "median_length_exact", "median_length", "layout", "qual_bases", "mean_q",
                 "mean_phred_arithmetic", "max_q", "unmapped", "is_aligned"]))
    # ---- eqx_nomd -------------------------------------------------------
    hdr, recs = build_eqx_nomd(seed + 5)
    put("eqx_nomd.sam", sam_text(hdr, recs))
    t = truth(recs, True)
    sections.append(sec(
        "eqx_nomd.sam",
        ["5 mapped reads whose CIGARs use `=` and `X` operations (minimap2 `--eqx` style), with no MD "
         "tag. Reads 1-3 carry NM; reads 4-5 have no NM, so their identity comes from the X "
         "operations (section 3.4 rule 2). Every read has at least one X operation, so "
         "`reads_without_md` is 0 and mismatches are counted from the X operations."],
        t))
    # ---- nonm -----------------------------------------------------------
    hdr, recs = build_nonm(seed + 6)
    put("nonm.sam", sam_text(hdr, recs))
    t = truth(recs, True)
    sections.append(sec(
        "nonm.sam",
        ["5 mapped reads with MD but no NM tag and plain `M` CIGARs: identity is not available for "
         "any of them (`identity_na` = 5, `identity_reads` = 0, `median_identity` null)."],
        t))
    # ---- ubam_like ------------------------------------------------------
    hdr, recs = build_ubam_like(seed + 7)
    put("ubam_like.sam", sam_text(hdr, recs))
    t = truth(recs, False)
    sections.append(sec(
        "ubam_like.sam",
        ["6 unaligned reads as Dorado or PacBio write them: no `@SQ` line (`@RG` and `@PG` only), every "
         "record flag 4, RG/qs/ns tags. `is_aligned` is 0: unaligned layout, no MAPQ table."],
        t, keys=["total_reads", "total_bases", "min_length", "max_length", "n50_exact", "n50",
                 "median_length_exact", "median_length", "layout", "qual_bases", "mean_q", "unmapped",
                 "is_aligned"]))
    # ---- sq_nomap -------------------------------------------------------
    hdr, recs = build_sq_nomap(seed + 8)
    put("sq_nomap.sam", sam_text(hdr, recs))
    t = truth(recs, True)
    sections.append(sec(
        "sq_nomap.sam",
        ["6 reads, every record flag 4, but the header has `@SQ` lines (chr1-chr3). `is_aligned` is 1 "
         "(header rule): aligned layout, the MAPQ table shows 100% unmapped, plots 5 and 6 show "
         "`No mapped reads.`"],
        t))
    # ---- mapq255 --------------------------------------------------------
    hdr, recs = build_mapq255(seed + 9)
    put("mapq255.sam", sam_text(hdr, recs))
    t = truth(recs, True)
    sections.append(sec(
        "mapq255.sam",
        ["7 mapped reads: five with MAPQ 255 (STAR-style unique hits), one with MAPQ 60, one with MAPQ 5. "
         "MAPQ 255 counts under `mapq_ge30`, so `mapq_ge30` = 6, `mapq_lt10` = 1."],
        t))
    # ---- single_len -----------------------------------------------------
    hdr, recs = build_single_len(seed + 10)
    put("single_len.sam", sam_text(hdr, recs))
    t = truth(recs, True)
    sections.append(sec(
        "single_len.sam",
        ["10 proper pairs (flags 99/147), every read 150 bp (`median_length` = 150, `min_length` = "
         "`max_length` = 150, short layout), about 0.7% substitutions, MAPQ 60 except pairs 7-8 (25) "
         "and 9 (3). R1 and R2 each have 10 reads."],
        t))
    # ---- hand cases -----------------------------------------------------
    hdr, recs = build_n50_five(seed + 11)
    put("n50_five.sam", sam_text(hdr, recs))
    t = truth(recs, True)
    assert t["n50"] == 400 and t["median_length"] == 300 and t["total_bases"] == 1600
    sections.append(sec(
        "n50_five.sam",
        ["Hand-computable N50. Five perfectly matching reads of 100, 200, 300, 400 and 600 bp. Total "
         "1600 bases; half is 800. From the longest: 600, then 600 + 400 = 1000 >= 800, so N50 = 400. "
         "Median length (3rd of 5) = 300."],
        t, keys=["total_reads", "total_bases", "n50", "median_length", "layout", "mapq_ge30",
                 "median_identity", "mismatch_bases"]))
    hdr, recs = build_identity_indel(seed + 12)
    put("identity_indel.sam", sam_text(hdr, recs))
    t = truth(recs, True)
    assert t["median_identity"] == 98.9 and t["identity_reads"] == 4
    sections.append(sec(
        "identity_indel.sam",
        ["Hand-computable gap-compressed identity, `(M - X) / (M + I_events + D_events)`, four 100 bp reads:",
         "- `id_ins1` `50M1I49M`, NM 1: M = 99, X = NM - 1 = 0, one insertion event: 99/100 = 99.0% (bin 990).",
         "- `id_del2_mm1` `50M2D50M`, NM 3, `MD:Z:30A19^GT50`: M = 100, X = 3 - 0 - 2 = 1, one deletion event "
         "(the 2 deleted bases count once): 99/101 = 98.0198% (bin 980).",
         "- `id_perfect` `100M`, NM 0: 100.0% (bin 1000).",
         "- `id_ins3` `40M3I57M`, NM 3: M = 97, X = 3 - 3 = 0, one insertion event of 3 bases: 97/98 = "
         "98.9796% (bin 989).",
         "Sorted bins 980, 989, 990, 1000: the lower median is bin 989, i.e. 98.9%."],
        t, keys=["total_reads", "total_bases", "identity_reads", "identity_na", "median_identity",
                 "median_identity_exact", "mismatch_bases", "insertion_events", "deletion_events",
                 "aligned_bases"]))
    hdr, recs = build_softclip_pct(seed + 13)
    put("softclip_pct.sam", sam_text(hdr, recs))
    t = truth(recs, True)
    assert abs(t["pct_clipped"] - 15.0) < 1e-9 and abs(t["pct_mapped"] - 80.0) < 1e-9
    sections.append(sec(
        "softclip_pct.sam",
        ["Hand-computable soft-clip percentage. Four mapped 100 bp reads (`10S80M10S`, `20S80M`, "
         "`80M20S`, `100M`) and one unmapped 100 bp read. Soft-clipped bases: 20 + 20 + 20 + 0 = 60 of "
         "the 400 bases of mapped primary reads: `pct_clipped` = 15.0%. The unmapped read is excluded "
         "from that denominator but not from `total_bases` (500); `pct_mapped` = 4/5 = 80.0%."],
        t, keys=["total_reads", "total_bases", "pct_mapped", "pct_mapq_ge30", "pct_clipped",
                 "clipped_bases", "unmapped", "mapq_ge30"]))

    readme = build_readme({"seed": seed, "sections": sections})
    put("README.md", readme)
    return written


def check_sizes(paths):
    total = 0
    bad = False
    for p in sorted(paths):
        n = os.path.getsize(p)
        total += n
        sys.stderr.write("%8d  %s\n" % (n, os.path.basename(p)))
        if n > 100 * 1024:
            sys.stderr.write("error: %s exceeds 100 KB\n" % p)
            bad = True
    sys.stderr.write("%8d  total\n" % total)
    if total > 300 * 1024:
        sys.stderr.write("error: total exceeds 300 KB\n")
        bad = True
    return not bad


# bulk mode: more reads than this need --force (200000 ONT-like reads are about 3.5 GB)
MAX_BULK_READS = 200_000


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--outdir", default=None, help="fixture directory (default: tests/data/longread)")
    ap.add_argument("--seed", type=int, default=240, help="random seed (default 240)")
    ap.add_argument("--test-view", default=None, help="path to htslib test/test_view")
    ap.add_argument("--reads", type=int, default=None,
                    help="bulk mode: write N reads as SAM for bench.sh instead of the fixtures")
    ap.add_argument("--profile", choices=["ont", "short"], default="ont",
                    help="bulk mode read profile (default ont: 500 bp-30 kb; short: 150 bp)")
    ap.add_argument("--out", default=None, help="bulk mode output file (default stdout)")
    ap.add_argument("--force", action="store_true",
                    help="bulk mode: allow more than %d reads" % MAX_BULK_READS)
    args = ap.parse_args()
    if args.reads is not None:
        if args.reads < 1:
            ap.error("--reads must be >= 1")
        if args.reads > MAX_BULK_READS and not args.force:
            ap.error("--reads %d is above the cap of %d; add --force to write it anyway"
                     % (args.reads, MAX_BULK_READS))
        out = open(args.out, "wb") if args.out else sys.stdout.buffer
        bulk(args.reads, args.profile, args.seed, out)
        return 0
    outdir = args.outdir or default_outdir()
    written = generate(outdir, args.seed, args.test_view)
    return 0 if check_sizes(written) else 1


if __name__ == "__main__":
    sys.exit(main())
