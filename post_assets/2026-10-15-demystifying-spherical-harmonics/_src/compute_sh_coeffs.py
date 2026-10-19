#!/usr/bin/env python3
"""Compute exact L1 and L2 Spherical Harmonics projection coefficients
for the Pantheon HDR equirectangular light probe (pantheon-hdr.webp).

Zero-dependency Python 3 script (uses `dwebp` to decode the WebP image).
Matches the coordinate system and inverse log-luminance decoding in
post_assets/2026-10-15-demystifying-spherical-harmonics/sh-diagrams.js.
"""

from pathlib import Path
import math
import subprocess
import tempfile


LOG_EPS = 0.002
LOG_MAX_LUM = 12.15625


def load_pantheon_lum(webp_path: Path):
    """Decode pantheon-hdr.webp and convert 8-bit log-encoded values to linear luminance."""
    with tempfile.NamedTemporaryFile(suffix=".ppm") as tmp:
        subprocess.run(
            ["dwebp", str(webp_path), "-ppm", "-o", tmp.name],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        data = Path(tmp.name).read_bytes()

    parts = data.split(maxsplit=4)
    width = int(parts[1])
    height = int(parts[2])
    raw = parts[4]

    lo = math.log2(LOG_EPS)
    hi = math.log2(LOG_MAX_LUM + LOG_EPS)
    lut = [0.0] * 256
    for i in range(1, 256):
        t = i / 255.0
        lut[i] = max(0.0, (2.0 ** (lo + t * (hi - lo))) - LOG_EPS)

    # Extract red channel and map through LUT
    lum = [lut[raw[i * 3]] for i in range(width * height)]
    return width, height, lum


def compute_sh_coefficients(width: int, height: int, lum: list) -> dict:
    """Integrate over the equirectangular grid using exact solid-angle weights per pixel."""
    # Precompute longitude (phi) trig tables
    cos_phi = [0.0] * width
    sin_phi = [0.0] * width
    for x in range(width):
        u = (x + 0.5) / width
        phi = (u - 0.5) * (2.0 * math.pi)
        cos_phi[x] = math.cos(phi)
        sin_phi[x] = math.sin(phi)

    m_1 = 0.0
    m_x = 0.0
    m_y = 0.0
    m_z = 0.0
    m_xy = 0.0
    m_xz = 0.0
    m_yz = 0.0
    m_x2 = 0.0
    m_y2 = 0.0
    m_z2 = 0.0

    for y in range(height):
        v = (y + 0.5) / height
        theta = (0.5 - v) * math.pi  # elevation angle in [-pi/2, pi/2]
        cos_t = math.cos(theta)
        sin_t = math.sin(theta)

        # Exact solid-angle weight for pixel (d(sin(theta)) * dphi / (4*pi))
        theta_top = (0.5 - y / height) * math.pi
        theta_bot = (0.5 - (y + 1.0) / height) * math.pi
        w_row = (math.sin(theta_top) - math.sin(theta_bot)) / (2.0 * width)

        row_offset = y * width
        sum_f = 0.0
        sum_fx = 0.0
        sum_fz = 0.0
        sum_fxz = 0.0
        sum_fx2 = 0.0
        sum_fz2 = 0.0

        for x in range(width):
            f = lum[row_offset + x]
            cp = cos_phi[x]
            sp = sin_phi[x]
            sum_f += f
            sum_fx += f * cp
            sum_fz += f * sp
            sum_fxz += f * cp * sp
            sum_fx2 += f * cp * cp
            sum_fz2 += f * sp * sp

        m_1 += w_row * sum_f
        m_x += w_row * cos_t * sum_fx
        m_y += w_row * sin_t * sum_f
        m_z += w_row * cos_t * sum_fz
        m_xy += w_row * cos_t * sin_t * sum_fx
        m_xz += w_row * cos_t * cos_t * sum_fxz
        m_yz += w_row * cos_t * sin_t * sum_fz
        m_x2 += w_row * cos_t * cos_t * sum_fx2
        m_y2 += w_row * sin_t * sin_t * sum_f
        m_z2 += w_row * cos_t * cos_t * sum_fz2

    # Orthogonal L2 form:
    # f_hat(w) = c_{0,0}
    #          + c_{1,-1} w_y + c_{1,0} w_z + c_{1,1} w_x
    #          + c_{2,-2} w_x w_y + c_{2,-1} w_y w_z + c_{2,0} (3 w_z^2 - 1)
    #          + c_{2,1} w_x w_z + c_{2,2} (w_x^2 - w_y^2)
    ortho = {
        "c00": m_1,
        "c1_1": 3.0 * m_y,
        "c10": 3.0 * m_z,
        "c11": 3.0 * m_x,
        "c2_2": 15.0 * m_xy,
        "c2_1": 15.0 * m_yz,
        "c20": 1.25 * (3.0 * m_z2 - m_1),
        "c21": 15.0 * m_xz,
        "c22": 3.75 * (m_x2 - m_y2),
    }

    # Fast Reconstruction Form (FRF):
    # f_hat(w) = c_0 + c_x w_x + c_y w_y + c_z w_z
    #          + c_{xy} w_x w_y + c_{xz} w_x w_z + c_{yz} w_y w_z
    #          + c_{x^2} w_x^2 + c_{y^2} w_y^2
    frf = {
        "c0": 6.0 * m_1 - 7.5 * m_x2 - 7.5 * m_y2,
        "cx": 3.0 * m_x,
        "cy": 3.0 * m_y,
        "cz": 3.0 * m_z,
        "cxy": 15.0 * m_xy,
        "cxz": 15.0 * m_xz,
        "cyz": 15.0 * m_yz,
        "cx2": -7.5 * m_1 + 15.0 * m_x2 + 7.5 * m_y2,
        "cy2": -7.5 * m_1 + 7.5 * m_x2 + 15.0 * m_y2,
    }

    return {"ortho": ortho, "frf": frf}


def main():
    asset_dir = Path(__file__).resolve().parent.parent
    webp_path = asset_dir / "pantheon-hdr.webp"
    width, height, lum = load_pantheon_lum(webp_path)
    coeffs = compute_sh_coefficients(width, height, lum)
    ortho = coeffs["ortho"]
    frf = coeffs["frf"]

    print("=== L1 Coefficients [c0, cx, cy, cz] ===")
    print(
        f"[{ortho['c00']:.6f}, {ortho['c11']:.6f}, {ortho['c1_1']:.6f}, {ortho['c10']:.6f}]"
    )

    print("\n=== L2 Orthogonal Coefficients ===")
    for k, v in ortho.items():
        print(f"  {k:5s} = {v:+.6f}")

    print("\n=== L2 Fast Reconstruction Form (FRF) Coefficients ===")
    for k, v in frf.items():
        print(f"  {k:5s} = {v:+.6f}")


if __name__ == "__main__":
    main()
