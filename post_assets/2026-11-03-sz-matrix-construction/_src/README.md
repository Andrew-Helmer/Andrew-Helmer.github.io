# SZ Matrix Construction Header Assets

This directory contains the reproducible LaTeX/TikZ source and build script for the header diagram of `_posts/2026-11-03-sz-matrix-construction.md`.

## Files
- `sz-matrix-header.tex`: TikZ source visualizing the top product $M = \mathcal{P}(E_0)\mathcal{P}(E_1)\mathcal{P}(E_2)\mathcal{P}(E_3)\cdots$ and the lower grid of basis and companion matrices ($E_0, E_1, E_2, E_3$, diagonal repeats $E'_i$, and companion matrices $\alpha_1, \alpha_2, \alpha_3$).
- `build_header_svg.py`: Post-processes the `dvisvgm --no-fonts` output to set display dimensions and ensure crisp dark-theme rendering.
- `build_header.sh`: Builds `../sz-matrix-header.svg`, `../sz-matrix-header.png`, and `../sz-matrix-header-meta.jpg`.

