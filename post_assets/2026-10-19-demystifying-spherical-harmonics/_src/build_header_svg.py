import base64
import re
import sys

raw_svg_path = sys.argv[1]
out_svg_path = sys.argv[2]
sphere_paths = [sys.argv[3], sys.argv[4], sys.argv[5]]

with open(raw_svg_path, "r", encoding="utf-8") as f:
    svg = f.read()

# Scale up the display width/height (preserving exact viewBox and aspect ratio)
scale_factor = 1.52
svg = re.sub(
    r"width='([0-9.]+)pt' height='([0-9.]+)pt'",
    lambda m: f"width='{float(m.group(1)) * scale_factor:.2f}pt' height='{float(m.group(2)) * scale_factor:.2f}pt'",
    svg,
    count=1,
)

# Replace the dvisvgm background path with a standard SVG <rect> so all renderers fill #1e1e1e
bg_pattern = r"<g fill='(#[0-9a-fA-F]{6})' stroke='\1'>\s*<path[^>]+/>\s*</g>"
svg, n_bg = re.subn(
    bg_pattern,
    lambda m: (
        f"<g fill='{m.group(1)}' stroke='none'>\n"
        f"<rect x='-158' y='-74' width='316' height='143' fill='{m.group(1)}'/>\n"
        f"</g>"
    ),
    svg,
    count=1,
)
assert n_bg == 1, "Failed to replace background rect"

markers = [
    ("#010101", -145.53552, -43.19101, 86.38276, 86.38202, sphere_paths[0]),
    ("#020202", -43.19101, -43.19101, 86.38202, 86.38202, sphere_paths[1]),
    ("#030303", 59.15276, -43.19101, 86.38276, 86.38202, sphere_paths[2]),
]

for color, x, y, w, h, img_path in markers:
    with open(img_path, "rb") as imgf:
        b64 = base64.b64encode(imgf.read()).decode("ascii")
    img_tag = (
        f"<g transform='scale(1,-1)'>\n"
        f"<image x='{x}' y='{y}' width='{w}' height='{h}' "
        f"preserveAspectRatio='none' xlink:href='data:image/png;base64,{b64}'/>\n"
        f"</g>"
    )
    pattern = rf"<g fill='{color}' stroke='{color}'>\s*<path[^>]+/>\s*</g>"
    svg, n = re.subn(pattern, img_tag, svg, count=1)
    assert n == 1, f"Failed to replace marker {color}"

with open(out_svg_path, "w", encoding="utf-8") as f:
    f.write(svg)

print("Successfully wrote", out_svg_path)

