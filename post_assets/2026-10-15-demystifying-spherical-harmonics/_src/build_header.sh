#!/usr/bin/env bash
set -euo pipefail
export PATH="/Library/TeX/texbin:/usr/local/bin:/opt/homebrew/bin:$PATH"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
POST_ASSETS_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BUILD_DIR="$(mktemp -d)"
trap 'rm -rf "$BUILD_DIR"' EXIT

echo "1. Decoding and un-pinching pantheon-hdr.webp..."
if [ -f "$SCRIPT_DIR/pantheon-hdr-raw.webp" ] && [ -f "$SCRIPT_DIR/cube_top_gen.webp" ] && [ -f "$SCRIPT_DIR/cube_bot_gen.webp" ]; then
  dwebp "$SCRIPT_DIR/pantheon-hdr-raw.webp" -ppm -o "$BUILD_DIR/pantheon_raw.ppm" >/dev/null 2>&1
  dwebp "$SCRIPT_DIR/cube_top_gen.webp" -ppm -o "$BUILD_DIR/cube_top_gen.ppm" >/dev/null 2>&1
  dwebp "$SCRIPT_DIR/cube_bot_gen.webp" -ppm -o "$BUILD_DIR/cube_bot_gen.ppm" >/dev/null 2>&1
  clang++ -O3 -std=c++17 "$SCRIPT_DIR/fix_poles.cpp" -o "$BUILD_DIR/fix_poles"
  "$BUILD_DIR/fix_poles" "$BUILD_DIR/pantheon_raw.ppm" "$BUILD_DIR/cube_top_gen.ppm" "$BUILD_DIR/cube_bot_gen.ppm" "$BUILD_DIR/pantheon.ppm" >/dev/null
  cwebp -lossless -m 6 "$BUILD_DIR/pantheon.ppm" -o "$POST_ASSETS_DIR/pantheon-hdr.webp" >/dev/null 2>&1
else
  dwebp "$POST_ASSETS_DIR/pantheon-hdr.webp" -ppm -o "$BUILD_DIR/pantheon.ppm" >/dev/null 2>&1
fi

echo "2. Compiling and running 3D sphere raytracer..."
clang++ -O3 "$SCRIPT_DIR/render_spheres.cpp" -o "$BUILD_DIR/render_spheres"
"$BUILD_DIR/render_spheres" "$BUILD_DIR/pantheon.ppm" "$BUILD_DIR/sphere"

for m in 0 1 2; do
  magick "$BUILD_DIR/sphere_mode${m}.ppm" "$BUILD_DIR/sphere_mode${m}_alpha.pgm" \
    -alpha off -compose CopyOpacity -composite "PNG32:$BUILD_DIR/sphere_mode${m}.png"
  cwebp -q 94 -alpha_q 100 -m 6 "$BUILD_DIR/sphere_mode${m}.png" -o "$BUILD_DIR/sphere_mode${m}.webp" >/dev/null 2>&1
done

echo "3. Compiling LaTeX labels and building spherical-harmonics-header.svg..."
cp "$SCRIPT_DIR/spherical-harmonics-header.tex" "$BUILD_DIR/"
(
  cd "$BUILD_DIR"
  latex -interaction=nonstopmode spherical-harmonics-header.tex >/dev/null
  dvisvgm --no-fonts -o raw_header.svg spherical-harmonics-header.dvi >/dev/null 2>&1
  python3 "$SCRIPT_DIR/build_header_svg.py" \
    raw_header.svg \
    "$POST_ASSETS_DIR/spherical-harmonics-header.svg" \
    sphere_mode0.webp sphere_mode1.webp sphere_mode2.webp
)

echo "4. Rendering spherical-harmonics-header.png and spherical-harmonics-header-meta.jpg..."
(
  cd "$BUILD_DIR"
  pdflatex -interaction=nonstopmode spherical-harmonics-header.tex >/dev/null
  gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r300 \
    -dTextAlphaBits=4 -dGraphicsAlphaBits=4 \
    -sOutputFile="$POST_ASSETS_DIR/spherical-harmonics-header.png" \
    spherical-harmonics-header.pdf >/dev/null
  magick "$POST_ASSETS_DIR/spherical-harmonics-header.png" \
    -resize 1200x -quality 92 "$POST_ASSETS_DIR/spherical-harmonics-header-meta.jpg"
)

echo "Done! Updated:"
echo "  - $POST_ASSETS_DIR/spherical-harmonics-header.svg"
echo "  - $POST_ASSETS_DIR/spherical-harmonics-header.png"
echo "  - $POST_ASSETS_DIR/spherical-harmonics-header-meta.jpg"
