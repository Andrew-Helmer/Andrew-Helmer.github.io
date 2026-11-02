#!/usr/bin/env bash
set -euo pipefail
export PATH="/Library/TeX/texbin:/usr/local/bin:/opt/homebrew/bin:$PATH"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
POST_ASSETS_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BUILD_DIR="$(mktemp -d)"
trap 'rm -rf "$BUILD_DIR"' EXIT

echo "1. Compiling TikZ diagram and building sz-matrix-header.svg..."
cp "$SCRIPT_DIR/sz-matrix-header.tex" "$BUILD_DIR/"
(
  cd "$BUILD_DIR"
  latex -interaction=nonstopmode sz-matrix-header.tex > latex.log 2>&1 || (cat latex.log && exit 1)
  dvisvgm --no-fonts -o raw_header.svg sz-matrix-header.dvi >/dev/null 2>&1
  python3 "$SCRIPT_DIR/build_header_svg.py" \
    raw_header.svg \
    "$POST_ASSETS_DIR/sz-matrix-header.svg"
)

echo "2. Rendering sz-matrix-header.png and sz-matrix-header-meta.jpg..."
(
  cd "$BUILD_DIR"
  pdflatex -interaction=nonstopmode sz-matrix-header.tex >/dev/null 2>&1
  gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r300 \
    -dTextAlphaBits=4 -dGraphicsAlphaBits=4 \
    -sOutputFile="$POST_ASSETS_DIR/sz-matrix-header.png" \
    sz-matrix-header.pdf >/dev/null 2>&1
  magick "$POST_ASSETS_DIR/sz-matrix-header.png" \
    -resize 1200x -quality 92 "$POST_ASSETS_DIR/sz-matrix-header-meta.jpg"
)

echo "Done! Updated:"
echo "  - $POST_ASSETS_DIR/sz-matrix-header.svg"
echo "  - $POST_ASSETS_DIR/sz-matrix-header.png"
echo "  - $POST_ASSETS_DIR/sz-matrix-header-meta.jpg"

