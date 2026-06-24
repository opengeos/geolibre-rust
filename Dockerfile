# syntax=docker/dockerfile:1
#
# Self-hostable image for the geolibre-rust browser demo.
#
# Stage 1 compiles the two WASM artifacts (wasm-bindgen library + WASI tool
# runner) with the same toolchain CI uses and assembles the static site exactly
# like .github/workflows/pages.yml. Stage 2 serves that site with nginx — no
# Rust, no GDAL, no server-side compute; every tool still runs in the visitor's
# browser.
#
#   docker build -t geolibre-rust .
#   docker run --rm -p 8080:80 geolibre-rust   # http://localhost:8080

# ---- Stage 1: build the WASM artifacts and assemble the site ----------------
FROM rust:1-bookworm AS builder

# wasm targets + tooling. binaryen provides wasm-opt (shrinks the WASI runner);
# clang/llvm/lld let cc-rs cross-compile C deps to wasm32-unknown-unknown (the
# browser-library target — GitHub's CI runners ship these preinstalled); wasi-sdk
# is the separate C toolchain for the wasm32-wasip1 WASI runner (zstd-sys).
RUN rustup target add wasm32-unknown-unknown wasm32-wasip1 \
 && apt-get update \
 && apt-get install -y --no-install-recommends binaryen curl ca-certificates clang llvm lld \
 && rm -rf /var/lib/apt/lists/* \
 && curl https://rustwasm.github.io/wasm-pack/installer/init.sh -sSf | sh

ENV WASI_SDK_VERSION=33.0
RUN curl -sL "https://github.com/WebAssembly/wasi-sdk/releases/download/wasi-sdk-33/wasi-sdk-${WASI_SDK_VERSION}-x86_64-linux.tar.gz" \
      | tar xz -C /opt
ENV WASI_SDK=/opt/wasi-sdk-${WASI_SDK_VERSION}-x86_64-linux \
    CC_wasm32_wasip1=/opt/wasi-sdk-${WASI_SDK_VERSION}-x86_64-linux/bin/clang \
    AR_wasm32_wasip1=/opt/wasi-sdk-${WASI_SDK_VERSION}-x86_64-linux/bin/llvm-ar \
    CFLAGS_wasm32_wasip1=--sysroot=/opt/wasi-sdk-${WASI_SDK_VERSION}-x86_64-linux/share/wasi-sysroot

WORKDIR /src
COPY . .

# Cache-bust token baked into index.html (mirrors pages.yml's __BUILD__ swap so
# returning visitors never mix an old cached loader with a fresh module). The
# publish workflow passes the commit SHA; local builds default to "docker".
ARG BUILD_REF=docker
RUN ./build.sh \
 && mkdir -p site \
 && cp demo/index.html           site/index.html \
 && cp npm/tools.mjs             site/tools.mjs \
 && cp npm/geolibre-cli.wasm     site/geolibre-cli.wasm \
 && cp npm/geolibre_wasm.js      site/geolibre_wasm.js \
 && cp npm/geolibre_wasm_bg.wasm site/geolibre_wasm_bg.wasm \
 && cp examples/sample.tif       site/sample.tif \
 && sed -i "s/__BUILD__/${BUILD_REF}/g" site/index.html

# ---- Stage 2: serve the static site -----------------------------------------
FROM nginx:alpine AS runtime

# Ensure .wasm/.mjs are served with the MIME types ES-module + WebAssembly
# streaming compilation require.
COPY docker/nginx.conf /etc/nginx/conf.d/default.conf
COPY --from=builder /src/site /usr/share/nginx/html

EXPOSE 80
