#!/usr/bin/env bash
# Serve the geolibre-rust demo locally.
#
# The page (index.html) needs sibling runtime files that live elsewhere in the
# repo: the JS runner, the WASI binary, the wasm-bindgen browser library from
# npm/, and a sample DEM from examples/. This script stages them into a temp dir
# alongside a copy of index.html (with the __BUILD__ cache-busting placeholder
# filled in), serves that dir over HTTP, and cleans up on exit. The repo's demo/
# stays clean.
#
#   ./demo/serve.sh            # serve on http://localhost:8000
#   ./demo/serve.sh 8731       # serve on a different port
set -euo pipefail

PORT="${1:-8000}"
DEMO_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$DEMO_DIR/.." && pwd)"

CLI_WASM="$ROOT/npm/geolibre-cli.wasm"
TOOLS_MJS="$ROOT/npm/tools.mjs"
LIB_JS="$ROOT/npm/geolibre_wasm.js"
LIB_WASM="$ROOT/npm/geolibre_wasm_bg.wasm"
SAMPLE="$ROOT/examples/sample.tif"

for f in "$CLI_WASM" "$TOOLS_MJS" "$LIB_JS" "$LIB_WASM" "$SAMPLE"; do
  if [ ! -f "$f" ]; then
    echo "Missing $f" >&2
    echo "Run ./build.sh first to produce the npm/ WASM artifacts." >&2
    exit 1
  fi
done

STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT

sed 's/__BUILD__/dev/g' "$DEMO_DIR/index.html" > "$STAGE/index.html"
cp "$TOOLS_MJS" "$CLI_WASM" "$LIB_JS" "$LIB_WASM" "$SAMPLE" "$STAGE/"

echo "Serving demo at http://localhost:$PORT/  (Ctrl-C to stop)"
cd "$STAGE"
exec python3 -m http.server "$PORT"
