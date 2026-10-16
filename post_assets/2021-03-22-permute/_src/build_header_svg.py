import re
import sys

raw_svg_path = sys.argv[1]
out_svg_path = sys.argv[2]

with open(raw_svg_path, "r", encoding="utf-8") as f:
    svg = f.read()

# Set display width/height to 600pt x 234pt (800px x 312px in CSS pixels),
# matching the tree-shuffling header image's intrinsic layout size.
svg = re.sub(
    r"width='([0-9.]+)pt' height='([0-9.]+)pt'",
    "width='600.00pt' height='234.00pt'",
    svg,
    count=1,
)

# Replace the dvisvgm background path with a standard SVG <rect> using the exact path bounds
# so all renderers reliably fill the #1e1e1e dark background.
bg_pattern = r"<g fill='(#[0-9a-fA-F]{6})' stroke='\1'>\s*<path d='([^']+)'[^>]*/>\s*</g>"


def replace_bg(m):
    color = m.group(1)
    nums = [float(v) for v in re.findall(r"[-+]?\d*\.?\d+", m.group(2))]
    xs = nums[0::2]
    ys = nums[1::2]
    rx = min(xs) - 1.0
    ry = min(ys) - 1.0
    rw = max(xs) - min(xs) + 2.0
    rh = max(ys) - min(ys) + 2.0
    return (
        f"<g fill='{color}' stroke='none'>\n"
        f"<rect x='{rx:.2f}' y='{ry:.2f}' width='{rw:.2f}' height='{rh:.2f}' fill='{color}'/>\n"
        f"</g>"
    )


svg, n_bg = re.subn(bg_pattern, replace_bg, svg, count=1)
assert n_bg == 1, "Failed to replace background rect"

# Ensure all font glyph <path> and <use> elements have stroke='none' so scaled TikZ nodes
# never inherit a stroke outline from enclosing <g> groups.
svg = re.sub(r"<path id='(g[^']+)' d='([^']+)'/>", r"<path id='\1' stroke='none' d='\2'/>", svg)
svg = re.sub(r"<use x='([^']+)' y='([^']+)' xlink:href='#(g[^']+)'/>", r"<use x='\1' y='\2' stroke='none' xlink:href='#\3'/>", svg)

with open(out_svg_path, "w", encoding="utf-8") as f:
    f.write(svg)

print("Successfully wrote", out_svg_path)

