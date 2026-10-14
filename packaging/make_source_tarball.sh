#!/usr/bin/env bash
# Build a source tarball that includes the tld submodule (git archive does not).
# Usage: packaging/make_source_tarball.sh <version> [output-dir] [ref]
#   version     e.g. 2.3.0 (tarball: samstat-<version>.tar.gz)
#   output-dir  default: dist
#   ref         git ref to archive, default: HEAD
# Writes <tarball> and <tarball>.sha256 into output-dir.
# Requires the submodules to be checked out (git submodule update --init --recursive).
set -euo pipefail

if [ $# -lt 1 ]; then
        echo "usage: $0 <version> [output-dir] [ref]" >&2
        exit 2
fi
VERSION="$1"
OUT="${2:-dist}"
REF="${3:-HEAD}"

ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"
mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"
NAME="samstat-${VERSION}"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

# Main repository content (no submodule content).
git archive --format=tar --prefix="${NAME}/" "$REF" | tar -x -C "$TMP"

# Add each submodule at the commit recorded in REF.
git ls-tree -r "$REF" | awk '$2 == "commit" {print $3 " " $4}' | while read -r sha path; do
        if [ ! -e "$ROOT/$path/.git" ]; then
                echo "submodule $path not initialised; run: git submodule update --init --recursive" >&2
                exit 1
        fi
        mkdir -p "$TMP/${NAME}/$path"
        git -C "$ROOT/$path" archive --format=tar "$sha" | tar -x -C "$TMP/${NAME}/$path"
done

tar -C "$TMP" -czf "$OUT/${NAME}.tar.gz" "$NAME"

cd "$OUT"
if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "${NAME}.tar.gz" > "${NAME}.tar.gz.sha256"
else
        shasum -a 256 "${NAME}.tar.gz" > "${NAME}.tar.gz.sha256"
fi
echo "wrote $OUT/${NAME}.tar.gz"
cat "${NAME}.tar.gz.sha256"
