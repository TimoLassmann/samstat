# CRAM fixtures (WP5.2)

The reference is the 5000 bp sequence of `chr1` (a 50 bp unit repeated 100 times), stored the
way htslib looks references up with `REF_PATH=<dir>/%s`: one file `ref/<md5>` holding the bare
sequence, no header, no newline. The CRAM headers carry `M5` but no `UR`, so the reference can
only be found through `REF_PATH` / `REF_CACHE`.

- `small.cram`  tests/data/small.sam (20 reads, MD tags) with `M5` added to the @SQ line.
- `lzma.cram`   300 random 100 bp reads from the reference (3% substitutions) plus a 3000 bp
                `XX:Z` tag; written with `-O cram,version=3.0,use_lzma,level=9`, so its blocks
                are LZMA coded (`samtools cram-size -v` shows `lzma`). 109 reads have
                MAPQ >= 30.

Regenerate (samtools 1.24; fixture generation only, never at ctest time), in a scratch
directory with an empty `ref/`:

    uv run python mkref.py            # writes ref.fa and ref/<md5>
    uv run python mklzma.py           # writes lz.sam and lz2.sam
    samtools view -C -T ref.fa --no-PG -o small.cram small_with_m5.sam
    samtools view -T ref.fa --no-PG -O cram,version=3.0,use_lzma,level=9 -o lzma.cram lz2.sam

then `samtools reheader` each CRAM with a header file holding only the @HD and @SQ (with M5)
lines, which drops the `UR:` field samtools adds.
