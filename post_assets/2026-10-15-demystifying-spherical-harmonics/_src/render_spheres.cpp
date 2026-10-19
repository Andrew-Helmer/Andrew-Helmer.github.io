#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>
#include <string>

struct Vec3 {
  double x, y, z;
  Vec3(double x_ = 0, double y_ = 0, double z_ = 0) : x(x_), y(y_), z(z_) {}
  Vec3 operator+(const Vec3& o) const { return Vec3(x + o.x, y + o.y, z + o.z); }
  Vec3 operator-(const Vec3& o) const { return Vec3(x - o.x, y - o.y, z - o.z); }
  Vec3 operator*(double s) const { return Vec3(x * s, y * s, z * s); }
};

double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
Vec3 cross(const Vec3& a, const Vec3& b) {
  return Vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}
Vec3 normalize(const Vec3& v) {
  double len = std::sqrt(dot(v, v));
  return len > 0 ? v * (1.0 / len) : Vec3(0, 0, 0);
}
double clampd(double x, double lo, double hi) {
  return std::max(lo, std::min(hi, x));
}
double smoothstep(double edge0, double edge1, double x) {
  double t = clampd((x - edge0) / (edge1 - edge0), 0.0, 1.0);
  return t * t * (3.0 - 2.0 * t);
}

// Frostbite exponential ("Photographic") luma range compression
double rangeCompress(double x) {
  return 1.0 - std::exp(-x);
}
double rangeCompress(double val, double threshold) {
  double v1 = val;
  double v2 = threshold + (1.0 - threshold) * rangeCompress((val - threshold) / (1.0 - threshold));
  return val < threshold ? v1 : v2;
}

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr, "Usage: %s <pantheon.ppm> <out_prefix>\n", argv[0]);
    return 1;
  }

  FILE* f = std::fopen(argv[1], "rb");
  if (!f) {
    std::perror("fopen");
    return 1;
  }
  char magic[16];
  int texW = 0, texH = 0, maxVal = 0;
  if (std::fscanf(f, "%15s %d %d %d", magic, &texW, &texH, &maxVal) != 4) {
    std::fprintf(stderr, "Invalid PPM header\n");
    return 1;
  }
  std::fgetc(f); // single whitespace

  std::vector<unsigned char> rgb(texW * texH * 3);
  std::fread(rgb.data(), 1, rgb.size(), f);
  std::fclose(f);

  // Inverse log-luminance LUT matching loadPantheonHDR() in assets/sh-diagrams.js
  const double LOG_EPS = 0.002;
  const double LOG_MAX_LUM = 12.15625;
  const double lo = std::log2(LOG_EPS);
  const double hi = std::log2(LOG_MAX_LUM + LOG_EPS);
  double lut[256];
  lut[0] = 0.0;
  for (int i = 1; i < 256; ++i) {
    double t = i / 255.0;
    lut[i] = std::max(0.0, std::pow(2.0, lo + t * (hi - lo)) - LOG_EPS);
  }

  std::vector<double> lumData(texW * texH);
  for (int i = 0; i < texW * texH; ++i) {
    lumData[i] = lut[rgb[i * 3]];
  }

  auto sampleEquirectLum = [&](const Vec3& dir) -> double {
    Vec3 d = normalize(dir);
    const double PI = 3.14159265358979323846;
    double phi = std::atan2(d.z, d.x);
    double theta = std::asin(clampd(d.y, -1.0, 1.0));
    double u = phi / (2.0 * PI) + 0.5;
    double v = 0.5 - theta / PI;

    // Bilinear filtering: REPEAT in U, CLAMP_TO_EDGE in V
    double px = u * texW - 0.5;
    double py = v * texH - 0.5;
    int x0 = (int)std::floor(px);
    int y0 = (int)std::floor(py);
    double fx = px - x0;
    double fy = py - y0;

    auto getTex = [&](int x, int y) -> double {
      int wx = ((x % texW) + texW) % texW;
      int cy = std::max(0, std::min(texH - 1, y));
      return lumData[cy * texW + wx];
    };

    double v00 = getTex(x0, y0);
    double v10 = getTex(x0 + 1, y0);
    double v01 = getTex(x0, y0 + 1);
    double v11 = getTex(x0 + 1, y0 + 1);
    return (v00 * (1.0 - fx) + v10 * fx) * (1.0 - fy) +
           (v01 * (1.0 - fx) + v11 * fx) * fy;
  };

  // Precomputed L1 SH coefficients of the Pantheon HDR probe (from assets/sh-diagrams.js)
  const double c0 = 0.095252;
  const Vec3 cVec(0.087969, 0.184707, 0.049513);

  // Camera parameters matching initPantheonFitDiagram in assets/sh-diagrams.js
  const double uYaw = 0.38;
  const double uPitch = 0.32;
  const double uDist = 2.85;
  const double tanFov = 0.382;

  double cp = std::cos(uPitch), sp = std::sin(uPitch);
  double cy = std::cos(uYaw),   sy = std::sin(uYaw);
  Vec3 camForward = normalize(Vec3(-sy * cp, -sp, -cy * cp));
  Vec3 worldUp(0.0, 1.0, 0.0);
  Vec3 camRight = normalize(cross(camForward, worldUp));
  Vec3 camUp = cross(camRight, camForward);
  Vec3 ro = camForward * (-uDist);

  const int imgSize = 560; // High-res per sphere
  const int ss = 4;        // 4x4 supersampling on top of analytical AA

  for (int mode = 0; mode < 3; ++mode) {
    std::vector<unsigned char> outRgb(imgSize * imgSize * 3);
    std::vector<unsigned char> outAlpha(imgSize * imgSize);

    for (int py = 0; py < imgSize; ++py) {
      for (int px = 0; px < imgSize; ++px) {
        Vec3 accumCol(0, 0, 0);
        double accumAlpha = 0.0;
        for (int sy_i = 0; sy_i < ss; ++sy_i) {
          for (int sx_i = 0; sx_i < ss; ++sx_i) {
            double u = (px + (sx_i + 0.5) / ss) / imgSize;
            double v = 1.0 - (py + (sy_i + 0.5) / ss) / imgSize;
            double ndcX = u * 2.0 - 1.0;
            double ndcY = v * 2.0 - 1.0;

            Vec3 rd = normalize(camForward + camRight * (ndcX * tanFov) + camUp * (ndcY * tanFov));
            double b = dot(ro, rd);
            double c = dot(ro, ro) - 1.0;
            double h = b * b - c;

            double distFromCenter = std::sqrt(std::max(0.0, 1.0 - h));
            double pixelRadius = (2.0 * tanFov * uDist) / imgSize;
            double coverage = 1.0 - smoothstep(1.0 - pixelRadius * 0.75, 1.0 + pixelRadius * 0.75, distFromCenter);

            if (coverage <= 0.0) {
              continue;
            }

            double clampedH = std::max(0.0, h);
            double t = -b - std::sqrt(clampedH);
            Vec3 pos = ro + rd * t;
            Vec3 n = normalize(pos);

            double rawLum = sampleEquirectLum(n);
            double shRaw = c0 + dot(cVec, n);

            Vec3 col;
            if (mode == 0) {
              // Mode 0: Original spherical function f(w)
              double mapped = rangeCompress(rawLum, 0.25);
              double ldr = std::pow(clampd(mapped, 0.0, 1.0), 1.0 / 2.2);
              col = Vec3(ldr, ldr, ldr);
            } else if (mode == 1) {
              // Mode 1: L1 SH Approximation f_hat(w)
              double mapped = rangeCompress(std::max(0.0, shRaw), 0.25);
              double ldr = std::pow(clampd(mapped, 0.0, 1.0), 1.0 / 2.2);
              col = Vec3(ldr, ldr, ldr);
            } else {
              // Mode 2: Absolute Error |f_hat(w) - f(w)|
              double absErr = std::abs(rawLum - shRaw);
              double mapped = rangeCompress(absErr, 0.25);
              double ldr = std::pow(clampd(mapped, 0.0, 1.0), 1.0 / 2.2);
              col = Vec3(ldr, ldr * 0.16, ldr * 0.12);
            }

            accumCol = accumCol + col * coverage;
            accumAlpha += coverage;
          }
        }
        double alpha = accumAlpha / (ss * ss);
        Vec3 unpremult = (accumAlpha > 1e-6) ? (accumCol * (1.0 / accumAlpha)) : Vec3(0, 0, 0);
        int pIdx = py * imgSize + px;
        outRgb[pIdx * 3 + 0] = (unsigned char)clampd(std::round(unpremult.x * 255.0), 0.0, 255.0);
        outRgb[pIdx * 3 + 1] = (unsigned char)clampd(std::round(unpremult.y * 255.0), 0.0, 255.0);
        outRgb[pIdx * 3 + 2] = (unsigned char)clampd(std::round(unpremult.z * 255.0), 0.0, 255.0);
        outAlpha[pIdx] = (unsigned char)clampd(std::round(alpha * 255.0), 0.0, 255.0);
      }
    }

    char fnameRgb[512], fnameAlpha[512];
    std::snprintf(fnameRgb, sizeof(fnameRgb), "%s_mode%d.ppm", argv[2], mode);
    std::snprintf(fnameAlpha, sizeof(fnameAlpha), "%s_mode%d_alpha.pgm", argv[2], mode);

    FILE* out = std::fopen(fnameRgb, "wb");
    std::fprintf(out, "P6\n%d %d\n255\n", imgSize, imgSize);
    std::fwrite(outRgb.data(), 1, outRgb.size(), out);
    std::fclose(out);

    FILE* outA = std::fopen(fnameAlpha, "wb");
    std::fprintf(outA, "P5\n%d %d\n255\n", imgSize, imgSize);
    std::fwrite(outAlpha.data(), 1, outAlpha.size(), outA);
    std::fclose(outA);
  }
  return 0;
}

