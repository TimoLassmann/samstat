#!/bin/sh
# SAMStat benchmark: not a CTest, no gates. Timing is noisy; read the ratios.
#
# usage: bench.sh SAMSTAT [OLD_SAMSTAT]
#   SAMSTAT      samstat binary to measure (build dir, outside the repo)
#   OLD_SAMSTAT  optional 2.3.0 binary, timed on the same files
#
# environment:
#   BENCH_DIR          work directory, REQUIRED; must not be inside a git work tree
#                      (the inputs are hundreds of MB; never let them land in a repo
#                      or a backed-up source checkout)
#   BENCH_SHORT_READS  150 bp reads                  (default 200000,  at most 200000)
#   BENCH_ONT_READS    ONT-like reads, 500 bp-30 kb  (default 10000,   at most 20000)
#   BENCH_HUGE_READS   reads of 1 Mb                 (default 8,       at most 8)
#   BENCH_THREADS      thread counts to run          (default "1 4 8", each at most 8)
#   BENCH_FORCE        1 lifts the caps above (for the large published runs, e.g.
#                      BENCH_SHORT_READS=2000000 BENCH_ONT_READS=150000, about 3 GB)
#   BENCH_PYTHON       python launcher               (default "uv run python")
#
# The expected disk use is printed before anything is generated. Generation time
# is reported separately from the timed runs. Inputs are SAM plus a FASTQ.gz made
# from the short-read SAM; they are reused when present (delete BENCH_DIR to
# regenerate).
set -eu

if [ $# -lt 1 ] || [ ! -x "$1" ]; then
    echo "usage: $0 SAMSTAT [OLD_SAMSTAT]" >&2
    exit 2
fi
NEW=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
OLD=
if [ $# -ge 2 ]; then
    [ -x "$2" ] || { echo "not executable: $2" >&2; exit 2; }
    OLD=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
fi

HERE=$(cd "$(dirname "$0")" && pwd)
GEN="$HERE/../tools/make_longread_fixtures.py"
[ -f "$GEN" ] || { echo "generator not found: $GEN" >&2; exit 2; }

# No default work directory: a default inside the current directory once put
# gigabytes of inputs into a source checkout.
if [ -z "${BENCH_DIR:-}" ]; then
    echo "BENCH_DIR is not set: set it to a scratch directory outside any repository," >&2
    echo "e.g. BENCH_DIR=/tmp/samstat-bench sh $0 $*" >&2
    exit 2
fi
: "${BENCH_SHORT_READS:=200000}"
: "${BENCH_ONT_READS:=10000}"
: "${BENCH_HUGE_READS:=8}"
: "${BENCH_THREADS:=1 4 8}"
: "${BENCH_FORCE:=0}"
: "${BENCH_PYTHON:=uv run python}"

# is_count VALUE: a positive decimal integer
is_count() {
    case $1 in
        ''|*[!0-9]*|0*) return 1 ;;
    esac
    return 0
}
# check_cap NAME VALUE MAX: refuse VALUE > MAX unless BENCH_FORCE=1
check_cap() {
    is_count "$2" || { echo "$1 must be a positive integer, got '$2'" >&2; exit 2; }
    if [ "$BENCH_FORCE" != 1 ] && [ "$2" -gt "$3" ]; then
        echo "$1=$2 is above the cap of $3; set BENCH_FORCE=1 to run it anyway" >&2
        exit 2
    fi
}
check_cap BENCH_SHORT_READS "$BENCH_SHORT_READS" 200000
check_cap BENCH_ONT_READS "$BENCH_ONT_READS" 20000
check_cap BENCH_HUGE_READS "$BENCH_HUGE_READS" 8
for t in $BENCH_THREADS; do
    check_cap BENCH_THREADS "$t" 8
done

# Refuse a work directory inside a git work tree: the nearest existing ancestor
# decides, since BENCH_DIR may not exist yet.
case $BENCH_DIR in
    /*) bench_abs=$BENCH_DIR ;;
    *)  bench_abs=$PWD/$BENCH_DIR ;;
esac
probe=$bench_abs
while [ ! -d "$probe" ]; do
    probe=$(dirname "$probe")
done
if command -v git >/dev/null 2>&1 &&
   [ "$(git -C "$probe" rev-parse --is-inside-work-tree 2>/dev/null || true)" = true ]; then
    echo "BENCH_DIR=$BENCH_DIR is inside the git work tree $(git -C "$probe" rev-parse --show-toplevel);" >&2
    echo "choose a directory outside any repository" >&2
    exit 2
fi

# Expected disk use (rounded up): about 380 bytes per 150 bp SAM record plus
# 120 for its FASTQ.gz copy, 17.5 kB per ONT-like record (mean length 7.2 kb,
# about 2.4 bytes per base) and 2 MB per 1 Mb record.
est_mb=$(( (BENCH_SHORT_READS * 500 + BENCH_ONT_READS * 17500 + BENCH_HUGE_READS * 2000100) / 1000000 + 1 ))
echo "bench.sh: inputs in $bench_abs need about $est_mb MB of disk (less when they are reused)"

case "$(uname -s)" in
    Darwin) TIME_KIND=bsd ;;
    *)      TIME_KIND=gnu ;;
esac
[ -x /usr/bin/time ] || { echo "/usr/bin/time not found" >&2; exit 2; }

mkdir -p "$BENCH_DIR"
BENCH_DIR=$(cd "$BENCH_DIR" && pwd)
OUT="$BENCH_DIR/out"
rm -rf "$OUT"
mkdir -p "$OUT"

now() { date +%s; }

# --- input generation ------------------------------------------------------
gen_t0=$(now)

SHORT="$BENCH_DIR/short_${BENCH_SHORT_READS}.sam"
ONT="$BENCH_DIR/ont_${BENCH_ONT_READS}.sam"
HUGE="$BENCH_DIR/huge_${BENCH_HUGE_READS}x1Mb.sam"
FQGZ="$BENCH_DIR/short_${BENCH_SHORT_READS}.fastq.gz"

# the generator has its own cap on --reads; BENCH_FORCE=1 lifts both
GEN_FORCE=
if [ "$BENCH_FORCE" = 1 ]; then GEN_FORCE=--force; fi
if [ ! -s "$SHORT" ]; then
    # shellcheck disable=SC2086
    $BENCH_PYTHON "$GEN" --reads "$BENCH_SHORT_READS" --profile short --out "$SHORT" $GEN_FORCE
fi
if [ ! -s "$ONT" ]; then
    # shellcheck disable=SC2086
    $BENCH_PYTHON "$GEN" --reads "$BENCH_ONT_READS" --profile ont --out "$ONT" $GEN_FORCE
fi
if [ ! -s "$HUGE" ]; then
    awk -v n="$BENCH_HUGE_READS" 'BEGIN {
        srand(240); L = 1000000; split("A C G T", b, " ");
        print "@HD\tVN:1.6\tSO:unsorted";
        print "@SQ\tSN:chr1\tLN:" (L * 3);
        for (r = 1; r <= n; r++) {
            s = ""; q = "";
            for (i = 0; i < L / 1000; i++) {
                chunk = ""; qc = "";
                for (j = 0; j < 1000; j++) { chunk = chunk b[int(rand() * 4) + 1]; qc = qc sprintf("%c", 41 + int(rand() * 25)); }
                s = s chunk; q = q qc;
            }
            printf "huge_%d\t0\tchr1\t%d\t60\t%dM\t*\t0\t0\t%s\t%s\n", r, 1 + (r - 1) * 1000, L, s, q;
        }
    }' > "$HUGE"
fi
if [ ! -s "$FQGZ" ]; then
    awk -F '\t' '!/^@/ { print "@" $1; print $10; print "+"; print $11 }' "$SHORT" | gzip -1 > "$FQGZ"
fi

gen_t1=$(now)

# --- helpers ---------------------------------------------------------------
# count_sam FILE / count_fq FILE -> "reads bases"
count_sam() {
    awk -F '\t' '!/^@/ { n++; if ($10 != "*") b += length($10) } END { printf "%d %d\n", n, b }' "$1"
}
count_fq() {
    gzip -dc "$1" | awk 'NR % 4 == 2 { n++; b += length($0) } END { printf "%d %d\n", n, b }'
}

# run_one BIN THREADS FILE -> "wall_s rss_mb"; exits the script if samstat fails
run_one() {
    bin=$1; thr=$2; file=$3
    tlog="$OUT/time.log"
    rc=0
    rm -rf "$OUT/rep"
    if [ "$TIME_KIND" = bsd ]; then
        /usr/bin/time -l "$bin" -t "$thr" -d "$OUT/rep" "$file" >"$OUT/stdout.log" 2>"$tlog" || rc=$?
    else
        /usr/bin/time -v "$bin" -t "$thr" -d "$OUT/rep" "$file" >"$OUT/stdout.log" 2>"$tlog" || rc=$?
    fi
    if [ "$rc" -ne 0 ]; then
        echo "samstat failed (exit $rc): $bin -t $thr $file" >&2
        cat "$tlog" >&2
        exit 1
    fi
    if [ "$TIME_KIND" = bsd ]; then
        awk '/ real / { w = $1 } /maximum resident set size/ { r = $1 / 1048576 } END { printf "%.2f %.0f\n", w, r }' "$tlog"
    else
        awk -F': ' '/Elapsed \(wall clock\)/ { n = split($2, p, ":"); w = 0; for (i = 1; i <= n; i++) w = w * 60 + p[i] }
                    /Maximum resident set size/ { r = $2 / 1024 }
                    END { printf "%.2f %.0f\n", w, r }' "$tlog"
    fi
}

# bench LABEL BIN FILE READS BASES
bench() {
    label=$1; bin=$2; file=$3; reads=$4; bases=$5
    for t in $BENCH_THREADS; do
        res=$(run_one "$bin" "$t" "$file")
        awk -v l="$label" -v t="$t" -v f="$(basename "$file")" -v res="$res" -v n="$reads" -v b="$bases" 'BEGIN {
            split(res, a, " "); w = a[1]; r = a[2];
            if (w <= 0) w = 0.01;
            printf "%-6s %-30s -t %-2d %8.2f s %6d MB %9.3f Mreads/s %9.2f Mbases/s\n", l, f, t, w, r, n / w / 1e6, b / w / 1e6 }'
    done
}

# --- report ----------------------------------------------------------------
echo "SAMStat benchmark"
echo "date:       $(date -u +%Y-%m-%dT%H:%MZ)"
echo "machine:    $(uname -srm)"
if [ "$(uname -s)" = Darwin ]; then
    echo "cpu:        $(sysctl -n machdep.cpu.brand_string 2>/dev/null || echo unknown), $(sysctl -n hw.ncpu) logical cores"
else
    echo "cpu:        $(grep -m1 'model name' /proc/cpuinfo 2>/dev/null | cut -d: -f2- | sed 's/^ //'), $(nproc 2>/dev/null || echo ?) logical cores"
fi
echo "samstat:    $("$NEW" --version 2>&1 | head -1)"
if [ -n "$OLD" ]; then echo "reference:  $("$OLD" --version 2>&1 | head -1)"; fi

# htslib version and libdeflate: from the build tree above the binary
# (bundled htslib), else from pkg-config (system htslib).
BUILD_ROOT=$(dirname "$(dirname "$NEW")")
HTS_VER=unknown; DEFL=unknown
for h in "$BUILD_ROOT/htslib/src/htslib/version.h" "$BUILD_ROOT/htslib/include/htslib/hts.h"; do
    if [ -f "$h" ]; then
        HTS_VER=$(sed -n 's/^#define HTS_VERSION_TEXT "\(.*\)"/\1/p' "$h" | head -1)
        if [ -z "$HTS_VER" ]; then HTS_VER=unknown; fi
        break
    fi
done
if [ -f "$BUILD_ROOT/htslib/src/htslib/config.h" ]; then
    if grep -q '^#define HAVE_LIBDEFLATE 1' "$BUILD_ROOT/htslib/src/htslib/config.h"; then DEFL=yes; else DEFL=no; fi
fi
if [ "$HTS_VER" = unknown ] && command -v pkg-config >/dev/null 2>&1; then
    HTS_VER=$(pkg-config --modversion htslib 2>/dev/null || echo unknown)
    if pkg-config --libs --static htslib 2>/dev/null | grep -q deflate; then DEFL=yes; fi
fi
echo "htslib:     $HTS_VER (libdeflate: $DEFL)"
echo "generation: $((gen_t1 - gen_t0)) s (not part of the timed runs; 0 when inputs were reused)"
echo "threads:    $BENCH_THREADS"
echo

for spec in "short:$SHORT:sam" "ont:$ONT:sam" "huge:$HUGE:sam" "fastqgz:$FQGZ:fq"; do
    name=${spec%%:*}; rest=${spec#*:}; file=${rest%:*}; kind=${rest##*:}
    if [ "$kind" = fq ]; then counts=$(count_fq "$file"); else counts=$(count_sam "$file"); fi
    reads=${counts% *}; bases=${counts#* }
    echo "# $name: $reads reads, $bases bases"
    bench new "$NEW" "$file" "$reads" "$bases"
    if [ -n "$OLD" ]; then bench old "$OLD" "$file" "$reads" "$bases"; fi
done
rm -rf "$OUT"
