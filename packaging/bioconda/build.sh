#!/usr/bin/env bash
set -euo pipefail

cmake -S . -B build \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
        -DGIT_SUBMODULE=OFF
cmake --build build -j"${CPU_COUNT:-2}"
cmake --install build
