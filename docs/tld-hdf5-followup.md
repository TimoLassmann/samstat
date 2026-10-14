# tld: making HDF5 optional (follow-up)

## Problem

tld (pinned at f745dd8) defines `BOOT` publicly when built with
`TLD_USE_HDF5=OFF`. `BOOT` guards both the hdf5 module and the seq module
(`seq/tld-seq.h`), so samstat loses `FILE_TYPE_FASTA`, the sequence reader
API and (transitively) `PRId64`. samstat therefore has to build with HDF5
ON even though it never uses it.

## The patch

Two commits on top of f745dd8, branch `fix/hdf5-optional-keeps-seq`:

1. `Keep seq module when TLD_USE_HDF5=OFF; add TLD_NO_HDF5` - a new PUBLIC
   define `TLD_NO_HDF5` guards only the hdf5 headers/sources (`tld.h`,
   `tld.c`). `BOOT` is kept for the gentable bootstrap build only. The
   `unit-hdf5` test is registered only when `TLD_USE_HDF5` is ON.
2. `Include inttypes.h from tld.h` - `PRId64` and friends were only
   available through the hdf5 headers.

Files (git format-patch output):

- `docs/tld-patches/0001-Keep-seq-module-when-TLD_USE_HDF5-OFF-add-TLD_NO_HDF.patch`
- `docs/tld-patches/0002-Include-inttypes.h-from-tld.h.patch`

Apply them in a tld checkout at f745dd8 with `git am docs/tld-patches/*.patch`.

Verified: tld ctest passes with HDF5 ON (36 tests) and OFF (35 tests, no
unit-hdf5). samstat built against the patched tld with `-DTLD_USE_HDF5=OFF`
passes ctest (29 smoke tests) and `otool -L` lists only libSystem and libz.

## Steps once the patch is upstream

1. Apply to tld (`git am` the two patches) and push to
   github.com/timolassmann/tld main.
2. In samstat: bump the `src/tld` submodule to the new commit.
3. In samstat's CMake, before `add_subdirectory(tld)`, add
   `set(TLD_USE_HDF5 OFF CACHE BOOL "" FORCE)`.
4. Remove HDF5 from CI (install steps for libhdf5-dev / `brew install hdf5`)
   and from the build documentation.
5. Do a clean configure and build, run ctest, and check that `otool -L`
   (macOS) or `ldd` (Linux) shows no hdf5.
