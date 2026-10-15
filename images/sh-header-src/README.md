# Spherical Harmonics Header Image & Diagrams Build Guide

This directory contains the source code and build pipeline for the header image of `_posts/2026-10-19-demystifying-spherical-harmonics.md`, as well as documentation for a future AI or developer on how to regenerate the header image or restyle the post's visuals if the site's color scheme or typography changes.

---

## 1. Quick Regeneration (One Command)

From anywhere in the repository, run:

```bash
bash images/sh-header-src/build_header.sh
```

**Prerequisites (macOS):**
- TeX Live (`latex`, `pdflatex`, `dvisvgm`, `gs` — typically in `/Library/TeX/texbin` and `/usr/local/bin`)
- `clang++` (Xcode Command Line Tools)
- `dwebp` (`libwebp`) and `magick` (`ImageMagick`)
- `python3`

Running `build_header.sh` automatically updates all three output files in `images/`:
1. `images/spherical-harmonics-header.svg` — Used as the main post `image:` in the blog post front-matter. Contains crisp vector LaTeX math paths (`dvisvgm --no-fonts`) + embedded base64 32-bit RGBA PNGs of the three 3D spheres.
2. `images/spherical-harmonics-header.png` — High-resolution 300 DPI PNG rendered via `pdflatex` + Ghostscript (`gs`).
3. `images/spherical-harmonics-header-meta.jpg` — Social/OpenGraph preview image (`meta_image:` in the post front-matter), resized to 1200px wide.

---

## 2. How the Header Image Pipeline Works

1. **HDR Probe Decoding (`assets/pantheon-hdr.webp`):**
   - `assets/pantheon-hdr.webp` is a 1024×512 lossless 8-bit WebP storing log-encoded luminance values (`LOG_EPS = 0.002`, `LOG_MAX_LUM = 12.15625`).
   - `build_header.sh` decodes `assets/pantheon-hdr.webp` into a temporary PPM via `dwebp`.

2. **3D Sphere Raytracing (`images/sh-header-src/render_spheres.cpp`):**
   - Matches the exact WebGL2 fragment shader in `initPantheonFitDiagram` (`assets/sh-diagrams.js`), including:
     - Camera orientation (`uYaw = 0.38`, `uPitch = 0.32`, `uDist = 2.85`, `tanFov = 0.382`)
     - Exact L1 SH coefficients of the Pantheon probe: `[0.317151, 0.064818, 0.355999, -0.029964]`
     - Frostbite photographic display mapper (`rangeCompress(val, 0.25)`) and `1 / 2.2` gamma correction
     - Analytical silhouette anti-aliasing (`smoothstep`) + 4×4 supersampling per pixel
   - Outputs three 560×560 unpremultiplied RGB PPMs (`sphere_mode{0,1,2}.ppm`) and exact coverage alpha masks (`sphere_mode{0,1,2}_alpha.pgm`), which `magick` combines into 32-bit transparent RGBA PNGs (`sphere_mode{0,1,2}.png`):
     - `mode 0`: Original spherical function $f(\omega)$
     - `mode 1`: Clamped L1 SH approximation $\hat{f}(\omega)$
     - `mode 2`: Absolute error $\lvert\hat{f}(\omega) - f(\omega)\rvert$, tinted red (`vec3(ldr, ldr * 0.16, ldr * 0.12)`)

3. **LaTeX Layout & Vector Math (`images/spherical-harmonics-header.tex`):**
   - Uses `standalone` + `tikz` + `amsmath` + `mathtools`.
   - When compiled with `pdflatex` (`\pdfoutput > 0`), it embeds `sphere_mode{0,1,2}.png` directly via `\includegraphics`.
   - When compiled with `latex` (`\pdfoutput = 0`), it uses `\def\pgfsysdriver{pgfsys-dvisvgm.def}` and draws three placeholder marker squares (`#010101`, `#020202`, `#030303`).
   - `dvisvgm --no-fonts` converts the DVI into an SVG with pure `<path>` vector outlines for all math glyphs.
   - `images/sh-header-src/build_header_svg.py` replaces the three marker paths with `<image xlink:href="data:image/png;base64,...">` elements (self-contained so `<img src="/images/spherical-harmonics-header.svg">` works in all browsers) and ensures the background `<rect>` is a standard SVG `<rect>`.

---

## 3. What to Change If the Site's Style Changes

### A. Changing Colors (e.g., Light Mode or New Dark Palette)
1. **Header Image Background & Label Colors (`images/spherical-harmonics-header.tex`):**
   - Edit the `\definecolor` declarations at the top of `images/spherical-harmonics-header.tex`:
     ```latex
     \definecolor{bg}{HTML}{1E1E1E}        % Match $bg in _sass/_variables.scss
     \definecolor{textlight}{HTML}{ECECEC} % Match $text / $heading in _sass/_variables.scss
     \definecolor{errred}{HTML}{E84435}    % Red accent matching the error sphere
     ```
   - Because the three rendered spheres (`sphere_mode{0,1,2}.png`) have **transparent RGBA backgrounds**, you do *not* need to modify `render_spheres.cpp` if only the page background color changes (unless you want to change the red tint of the error sphere on line ~182 of `render_spheres.cpp`).
   - Re-run `bash images/sh-header-src/build_header.sh`.

2. **Interactive WebGL Diagrams (`assets/sh-diagrams.js` and `style.scss`):**
   - All 5 WebGL2 canvases in `assets/sh-diagrams.js` use transparent backgrounds (`gl.clearColor(0, 0, 0, 0)` with `alpha: true, premultipliedAlpha: true`), so they automatically adapt to any CSS page background color.
   - If switching to a light theme, check:
     - `.callout-box` and `.sh-diagram` button/slider styles in `style.scss` (lines 308–569).
     - The SVG $x, y, z$ axes gizmo colors in `createAxesSVG()` at the top of `assets/sh-diagrams.js`.
     - The incoming light arrow colors in Diagrams 2, 4, and 5 in `assets/sh-diagrams.js` (currently light gray/white to contrast against the dark `#1e1e1e` background).

### B. Changing Math or Display Fonts
- The math in the post body is rendered by **KaTeX** (configured in `_layouts/default.html`), which uses Computer Modern (`KaTeX_Main` / `KaTeX_Math`).
- The math in `images/spherical-harmonics-header.tex` uses LaTeX's default **Computer Modern** font to match KaTeX.
- If you change the math font in the future, add the corresponding font package (e.g., `\usepackage{newpxtext,newpxmath}` or `\usepackage{libertinust1math}`) to the preamble of `images/spherical-harmonics-header.tex` and run `bash images/sh-header-src/build_header.sh`.

