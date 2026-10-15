#!/usr/bin/env bash
set -euo pipefail
export PATH="/Library/TeX/texbin:/usr/local/bin:/opt/homebrew/bin:$PATH"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
IMAGES_DIR="$REPO_ROOT/images"
BUILD_DIR="$(mktemp -d)"
trap 'rm -rf "$BUILD_DIR"' EXIT

echo "1. Decoding assets/pantheon-hdr.webp..."
dwebp "$REPO_ROOT/assets/pantheon-hdr.webp" -ppm -o "$BUILD_DIR/pantheon.ppm" >/dev/null 2>&1

echo "2. Compiling and running 3D sphere raytracer..."
clang++ -O3 "$SCRIPT_DIR/render_spheres.cpp" -o "$BUILD_DIR/render_spheres"
"$BUILD_DIR/render_spheres" "$BUILD_DIR/pantheon.ppm" "$BUILD_DIR/sphere"

for m in 0 1 2; do
  magick "$BUILD_DIR/sphere_mode${m}.ppm" "$BUILD_DIR/sphere_mode${m}_alpha.pgm" \
    -alpha off -compose CopyOpacity -composite "PNG32:$BUILD_DIR/sphere_mode${m}.png"
done

echo "3. Compiling LaTeX labels and building spherical-harmonics-header.svg..."
cp "$IMAGES_DIR/spherical-harmonics-header.tex" "$BUILD_DIR/"
(
  cd "$BUILD_DIR"
  latex -interaction=nonstopmode spherical-harmonics-header.tex >/dev/null
  dvisvgm --no-fonts -o raw_header.svg spherical-harmonics-header.dvi >/dev/null 2>&1
  python3 "$SCRIPT_DIR/build_header_svg.py" \
    raw_header.svg \
    "$IMAGES_DIR/spherical-harmonics-header.svg" \
    sphere_mode0.png sphere_mode1.png sphere_mode2.png
)

echo "4. Rendering spherical-harmonics-header.png and spherical-harmonics-header-meta.jpg..."
(
  cd "$BUILD_DIR"
  pdflatex -interaction=nonstopmode spherical-harmonics-header.tex >/dev/null
  gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r300 \
    -dTextAlphaBits=4 -dGraphicsAlphaBits=4 \
    -sOutputFile="$IMAGES_DIR/spherical-harmonics-header.png" \
    spherical-harmonics-header.pdf >/dev/null
  magick "$IMAGES_DIR/spherical-harmonics-header.png" \
    -resize 1200x -quality 92 "$IMAGES_DIR/spherical-harmonics-header-meta.jpg"
)

echo "Done! Updated:"
echo "  - $IMAGES_DIR/spherical-harmonics-header.svg"
echo "  - $IMAGES_DIR/spherical-harmonics-header.png"
echo "  - $IMAGES_DIR/spherical-harmonics-header-meta.jpg"
