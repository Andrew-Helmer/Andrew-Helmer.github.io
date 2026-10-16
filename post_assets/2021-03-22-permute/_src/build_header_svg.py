import re
import sys

raw_svg_path = sys.argv[1]
out_svg_path = sys.argv[2]

with open(raw_svg_path, "r", encoding="utf-8") as f:
    svg = f.read()

# Set display height to 234pt (312px in CSS pixels, matching owen-shuffling.svg)
# and proportional width 405.60pt (541px) for the 13.38 x 7.72 bounding box.
svg = re.sub(
    r"width='([0-9.]+)pt' height='([0-9.]+)pt'",
    "width='405.60pt' height='234.00pt'",
    svg,
    count=1,
)

# Extract the page1 transform so we can map the exact TikZ bounding-box rectangle
# to root SVG coordinates and set a tight, perfectly centered viewBox (ignoring
# dvisvgm's internal origin excursion for scaled text nodes).
tx_match = re.search(
    r"transform='translate\(([-+0-9.]+),([-+0-9.]+)\)scale\(([-+0-9.]+),([-+0-9.]+)\)'",
    svg,
)
assert tx_match is not None, "Failed to find page1 transform"
tx, ty, sx, sy = (float(tx_match.group(i)) for i in range(1, 5))

bg_pattern = r"<g fill='(#[0-9a-fA-F]{6})' stroke='\1'>\s*<path d='([^']+)'[^>]*/>\s*</g>"
bbox_local = {}


def replace_bg(m):
    color = m.group(1)
    nums = [float(v) for v in re.findall(r"[-+]?\d*\.?\d+", m.group(2))]
    xs = nums[0::2]
    ys = nums[1::2]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    bbox_local["min_x"] = min_x
    bbox_local["max_x"] = max_x
    bbox_local["min_y"] = min_y
    bbox_local["max_y"] = max_y
    rx = min_x - 1.0
    ry = min_y - 1.0
    rw = (max_x - min_x) + 2.0
    rh = (max_y - min_y) + 2.0
    return (
        f"<g fill='{color}' stroke='none'>\n"
        f"<rect x='{rx:.2f}' y='{ry:.2f}' width='{rw:.2f}' height='{rh:.2f}' fill='{color}'/>\n"
        f"</g>"
    )


svg, n_bg = re.subn(bg_pattern, replace_bg, svg, count=1)
assert n_bg == 1, "Failed to replace background rect"

# Compute exact viewBox in root SVG coordinates from the TikZ bounding box
vb_x1 = tx + sx * bbox_local["min_x"]
vb_x2 = tx + sx * bbox_local["max_x"]
vb_y1 = ty + sy * bbox_local["max_y"]
vb_y2 = ty + sy * bbox_local["min_y"]
vb_x = min(vb_x1, vb_x2)
vb_y = min(vb_y1, vb_y2)
vb_w = abs(vb_x2 - vb_x1)
vb_h = abs(vb_y2 - vb_y1)

svg = re.sub(
    r"viewBox='[^']+'",
    f"viewBox='{vb_x:.5f} {vb_y:.5f} {vb_w:.5f} {vb_h:.5f}'",
    svg,
    count=1,
)

# Ensure all font glyph <path> and <use> elements have stroke='none' so scaled TikZ nodes
# never inherit a stroke outline from enclosing <g> groups.
svg = re.sub(r"<path id='(g[^']+)' d='([^']+)'/>", r"<path id='\1' stroke='none' d='\2'/>", svg)
svg = re.sub(r"<use x='([^']+)' y='([^']+)' xlink:href='#(g[^']+)'/>", r"<use x='\1' y='\2' stroke='none' xlink:href='#\3'/>", svg)

with open(out_svg_path, "w", encoding="utf-8") as f:
    f.write(svg)

print("Successfully wrote", out_svg_path)

