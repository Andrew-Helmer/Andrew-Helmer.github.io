# Permute Header Image Build Guide

This directory (`post_assets/2021-03-22-permute/_src/`) contains the TikZ source code and build script for the header image of `_posts/2021-03-22-permute.md`.

> **Note:** Because this folder is named `_src` (with a leading underscore), Jekyll automatically excludes it from the generated `_site/` build output while keeping the source files co-located right next to the post's published assets.

---

## 1. Quick Regeneration (One Command)

From the repository root, run:

```bash
bash post_assets/2021-03-22-permute/_src/build_header.sh
```

**Prerequisites (macOS):**
- TeX Live (`latex`, `pdflatex`, `dvisvgm`, `gs` — typically in `/Library/TeX/texbin` and `/usr/local/bin`)
- `magick` (`ImageMagick`)
- `python3`

Running `build_header.sh` automatically updates all three output files in `post_assets/2021-03-22-permute/`:
1. `permute-header.svg` — Used as the main post `image:` in the blog post front-matter. Contains pure vector paths (`dvisvgm --no-fonts`) scaled to 600pt × 234pt (800px × 312px).
2. `permute-header.png` — High-resolution 300 DPI PNG rendered via `pdflatex` + Ghostscript (`gs`).
3. `permute-header-meta.jpg` — Social/OpenGraph preview image (`meta_image:` in the post front-matter), resized to 1200px wide.

---

## 2. Changing Colors or Typography

Edit the `\definecolor` declarations at the top of `permute-header.tex`:

```latex
\definecolor{bg}{HTML}{1E1E1E}         % Match $bg in _sass/_variables.scss
\definecolor{treeline}{HTML}{BEBEBE}   % Highlighted boxes and target array borders
\definecolor{mutedline}{HTML}{555555}  % Background bijection arrows and lane dividers
\definecolor{mutedbox}{HTML}{6E6E6E}   % Unselected input and function box borders
\definecolor{digitcolor}{HTML}{D9D9D9} % Highlighted index numbers
\definecolor{muteddigit}{HTML}{808080} % Unselected index numbers
\definecolor{swapcyan}{HTML}{19D4FA}   % Highlighted stateless query paths
\definecolor{boxbg}{HTML}{252526}      % Function box fill color ($codeBg)
```

Then re-run `bash post_assets/2021-03-22-permute/_src/build_header.sh`.

