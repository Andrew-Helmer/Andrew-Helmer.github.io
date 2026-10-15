#!/usr/bin/env bash
set -euo pipefail
export PATH="/Library/TeX/texbin:/usr/local/bin:/opt/homebrew/bin:$PATH"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
POST_ASSETS_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BUILD_DIR="$(mktemp -d)"
trap 'rm -rf "$BUILD_DIR"' EXIT

echo "1. Compiling TikZ diagram and building permute-header.svg..."
cp "$SCRIPT_DIR/permute-header.tex" "$BUILD_DIR/"
(
  cd "$BUILD_DIR"
  latex -interaction=nonstopmode permute-header.tex >/dev/null
  dvisvgm --no-fonts -o raw_header.svg permute-header.dvi >/dev/null 2>&1
  python3 "$SCRIPT_DIR/build_header_svg.py" \
    raw_header.svg \
    "$POST_ASSETS_DIR/permute-header.svg"
)

echo "2. Rendering permute-header.png and permute-header-meta.jpg..."
(
  cd "$BUILD_DIR"
  pdflatex -interaction=nonstopmode permute-header.tex >/dev/null
  gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r300 \
    -dTextAlphaBits=4 -dGraphicsAlphaBits=4 \
    -sOutputFile="$POST_ASSETS_DIR/permute-header.png" \
    permute-header.pdf >/dev/null
  magick "$POST_ASSETS_DIR/permute-header.png" \
    -resize 1200x -quality 92 "$POST_ASSETS_DIR/permute-header-meta.jpg"
)

echo "Done! Updated:"
echo "  - $POST_ASSETS_DIR/permute-header.svg"
echo "  - $POST_ASSETS_DIR/permute-header.png"
echo "  - $POST_ASSETS_DIR/permute-header-meta.jpg"

