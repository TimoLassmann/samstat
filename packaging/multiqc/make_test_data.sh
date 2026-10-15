#!/usr/bin/env bash
# Write the example SAMStat JSON reports for MultiQC's test-data repository
# (data/modules/samstat/) from samstat's own test fixtures.
# Usage: packaging/multiqc/make_test_data.sh <samstat binary> <output dir>
set -euo pipefail
if [ $# -ne 2 ]; then
        echo "usage: $0 <samstat binary> <output dir>" >&2
        exit 2
fi
SAMSTAT="$1"
OUT="$2"
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
F="$ROOT/tests/data"
mkdir -p "$OUT"
# aligned short reads, paired reads, ONT-like and HiFi-like long reads, FASTQ, FASTA
"$SAMSTAT" -t 2 --json -d "$OUT" \
        "$F/small.sam" "$F/paired_small.sam" "$F/longread/ont_like.bam" \
        "$F/longread/hifi_like.sam" "$F/longread/long.fq" "$F/reads.fa.gz"
rm -f "$OUT"/*.samstat.html
ls -l "$OUT"
