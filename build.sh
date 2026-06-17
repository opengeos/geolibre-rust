#!/usr/bin/env bash
# Build both WASM artifacts and stage them into npm/ for publishing:
#   1. the wasm-bindgen browser library (geolibre_wasm*.{js,wasm,d.ts})
#   2. the WASI tool runner (geolibre-cli.wasm)
#
#   rustup target add wasm32-unknown-unknown wasm32-wasip1
#   cargo install wasm-pack
#   # optional, to shrink the WASI runner: apt-get install -y binaryen
#   ./build.sh
set -euo pipefail

cd "$(dirname "$0")"

# ── 1. Browser library (wasm-bindgen) ───────────────────────────────────────
echo "==> wasm-pack build geolibre-wasm (browser library)"
wasm-pack build crates/geolibre-wasm --release --target web \
  --out-dir "$PWD/target/gl-wasm-pkg"
cp target/gl-wasm-pkg/geolibre_wasm.js          npm/
cp target/gl-wasm-pkg/geolibre_wasm_bg.wasm     npm/
cp target/gl-wasm-pkg/geolibre_wasm.d.ts        npm/
cp target/gl-wasm-pkg/geolibre_wasm_bg.wasm.d.ts npm/

# ── 2. WASI tool runner ─────────────────────────────────────────────────────
echo "==> cargo build geolibre-cli (wasm32-wasip1, release)"
cargo build -p geolibre-cli --release --target wasm32-wasip1
CLI_WASM="target/wasm32-wasip1/release/geolibre.wasm"

if command -v wasm-opt >/dev/null 2>&1; then
  echo "==> wasm-opt -Oz the WASI runner"
  wasm-opt -Oz \
    --enable-bulk-memory \
    --enable-nontrapping-float-to-int \
    --enable-sign-ext \
    --enable-mutable-globals \
    "$CLI_WASM" -o npm/geolibre-cli.wasm
else
  echo "==> wasm-opt not found; shipping the WASI runner unoptimized"
  cp "$CLI_WASM" npm/geolibre-cli.wasm
fi

echo "==> staged npm/ artifacts:"
ls -lh npm/geolibre_wasm_bg.wasm npm/geolibre-cli.wasm | awk '{print $5"\t"$9}'
