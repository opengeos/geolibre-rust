#!/usr/bin/env bash
# Build the WASI tool runner and stage it into npm/ for publishing.
#
#   rustup target add wasm32-wasip1
#   cargo install wasm-opt   # or use the binaryen system package
#   ./build.sh
set -euo pipefail

cd "$(dirname "$0")"

TARGET=wasm32-wasip1
ARTIFACT="target/${TARGET}/release/geolibre.wasm"

echo "==> cargo build (${TARGET}, release)"
cargo build -p geolibre-cli --release --target "${TARGET}"

echo "==> wasm-opt -Oz"
# Match whitebox-wasm's feature overrides so older wasm-opt accepts the output.
wasm-opt -Oz \
  --enable-bulk-memory \
  --enable-nontrapping-float-to-int \
  --enable-sign-ext \
  --enable-mutable-globals \
  "${ARTIFACT}" -o npm/geolibre-cli.wasm

echo "==> staged npm/geolibre-cli.wasm"
ls -lh npm/geolibre-cli.wasm
