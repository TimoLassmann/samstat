#!/usr/bin/env bash
set -euo pipefail

cmake -S . -B build \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="${PREFIX}"
cmake --build build -j"${CPU_COUNT:-2}"
cmake --install build
