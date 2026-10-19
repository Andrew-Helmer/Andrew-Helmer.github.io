#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <string>

static inline double clampd(double x, double lo, double hi) {
  return std::max(lo, std::min(hi, x));
}

static inline double smoothstep(double a, double b, double x) {
  double t = clampd((x - a) / (b - a), 0.0, 1.0);
  return t * t * (3.0 - 2.0 * t);
}

const double PI = 3.14159265358979323846;
const double LOG_EPS = 0.002;
const double LOG_MAX_LUM = 12.15625;
const double LO_LOG2 = std::log2(LOG_EPS);
const double HI_LOG2 = std::log2(LOG_MAX_LUM + LOG_EPS);

double byteNormToLum(double t) {
  if (t <= 0.0) return 0.0;
  return std::max(0.0, std::pow(2.0, LO_LOG2 + t * (HI_LOG2 - LO_LOG2)) - LOG_EPS);
}

double lumToByteNorm(double L) {
  if (L <= 0.0) return 0.0;
  double t = (std::log2(L + LOG_EPS) - LO_LOG2) / (HI_LOG2 - LO_LOG2);
  return clampd(t, 0.0, 1.0);
}

double rangeCompress(double x) {
  return 1.0 - std::exp(-x);
}
double rangeCompress(double val, double threshold) {
  double v1 = val;
  double v2 = threshold + (1.0 - threshold) * rangeCompress((val - threshold) / (1.0 - threshold));
  return val < threshold ? v1 : v2;
}

double inverseTonemap(double Y) {
  Y = clampd(Y, 0.0, 1.0);
  double M = std::pow(Y, 2.2);
  const double thresh = 0.25;
  if (M <= thresh) {
    return M;
  }
  const double maxCompressed = rangeCompress(LOG_MAX_LUM, thresh);
  double M_scaled = thresh + (M - thresh) * ((maxCompressed - thresh) / (1.0 - thresh));
  double norm = clampd((M_scaled - thresh) / (1.0 - thresh), 0.0, 0.999999);
  double L = thresh - (1.0 - thresh) * std::log(1.0 - norm);
  return clampd(L, 0.0, LOG_MAX_LUM);
}

struct Image {
  int W, H;
  std::vector<double> luma;
};

Image loadLumaPPM(const char* path) {
  FILE* f = std::fopen(path, "rb");
  if (!f) { std::perror(path); std::exit(1); }
  char magic[16]; int maxVal = 0;
  Image img;
  if (std::fscanf(f, "%15s %d %d %d", magic, &img.W, &img.H, &maxVal) != 4) {
    std::fprintf(stderr, "Invalid PPM: %s\n", path);
    std::exit(1);
  }
  std::fgetc(f);
  std::vector<unsigned char> rgb(img.W * img.H * 3);
  std::fread(rgb.data(), 1, rgb.size(), f);
  std::fclose(f);
  img.luma.resize(img.W * img.H);
  for (int i = 0; i < img.W * img.H; ++i) {
    img.luma[i] = (0.2126 * rgb[i * 3 + 0] + 0.7152 * rgb[i * 3 + 1] + 0.0722 * rgb[i * 3 + 2]) / 255.0;
  }
  return img;
}

void saveGrayPPM(const char* path, int W, int H, const std::vector<unsigned char>& gray) {
  FILE* f = std::fopen(path, "wb");
  std::fprintf(f, "P6\n%d %d\n255\n", W, H);
  std::vector<unsigned char> rgb(W * H * 3);
  for (int i = 0; i < W * H; ++i) {
    rgb[i * 3 + 0] = gray[i];
    rgb[i * 3 + 1] = gray[i];
    rgb[i * 3 + 2] = gray[i];
  }
  std::fwrite(rgb.data(), 1, rgb.size(), f);
  std::fclose(f);
}

double sampleBilinearClamp(const Image& img, double u, double v) {
  double px = u * img.W - 0.5;
  double py = v * img.H - 0.5;
  int x0 = (int)std::floor(px);
  int y0 = (int)std::floor(py);
  double fx = px - x0;
  double fy = py - y0;
  auto getP = [&](int x, int y) -> double {
    int cx = std::max(0, std::min(img.W - 1, x));
    int cy = std::max(0, std::min(img.H - 1, y));
    return img.luma[cy * img.W + cx];
  };
  return (getP(x0, y0) * (1 - fx) + getP(x0 + 1, y0) * fx) * (1 - fy) +
         (getP(x0, y0 + 1) * (1 - fx) + getP(x0 + 1, y0 + 1) * fx) * fy;
}

double sampleBilinearRepeatU(const Image& img, double u, double v) {
  double px = u * img.W - 0.5;
  double py = v * img.H - 0.5;
  int x0 = (int)std::floor(px);
  int y0 = (int)std::floor(py);
  double fx = px - x0;
  double fy = py - y0;
  auto getP = [&](int x, int y) -> double {
    int wx = ((x % img.W) + img.W) % img.W;
    int cy = std::max(0, std::min(img.H - 1, y));
    return img.luma[cy * img.W + wx];
  };
  return (getP(x0, y0) * (1 - fx) + getP(x0 + 1, y0) * fx) * (1 - fy) +
         (getP(x0, y0 + 1) * (1 - fx) + getP(x0 + 1, y0 + 1) * fx) * fy;
}

int main(int argc, char** argv) {
  if (argc < 5) {
    std::fprintf(stderr, "Usage: %s <dome_top.ppm> <equirect_raw.ppm> <floor_bot.ppm> <out_pantheon.ppm> [preview.ppm]\n", argv[0]);
    return 1;
  }

  Image domeTop  = loadLumaPPM(argv[1]);
  Image wallRaw  = loadLumaPPM(argv[2]);
  Image floorBot = loadLumaPPM(argv[3]);
  const char* outHdrPath = argv[4];
  const char* outPreviewPath = (argc >= 6) ? argv[5] : nullptr;

  // Normalize out the built-in radial vignette of floorBot so that mirror-tiling across tile centers has zero brightness seams
  {
    const int S = floorBot.W;
    const int G = 8;
    std::vector<double> gridMean(G * G, 0.40);
    for (int gy = 0; gy < G; ++gy) {
      for (int gx = 0; gx < G; ++gx) {
        double sum = 0.0, cnt = 0.0;
        for (int y = (gy * S) / G; y < ((gy + 1) * S) / G; ++y) {
          for (int x = (gx * S) / G; x < ((gx + 1) * S) / G; ++x) {
            sum += floorBot.luma[y * S + x];
            cnt += 1.0;
          }
        }
        gridMean[gy * G + gx] = sum / cnt;
      }
    }
    for (int pass = 0; pass < 12; ++pass) {
      std::vector<double> next = gridMean;
      for (int gy = 0; gy < G; ++gy) {
        for (int gx = 0; gx < G; ++gx) {
          double s = gridMean[gy * G + gx] * 4.0;
          double w = 4.0;
          if (gy > 0)     { s += gridMean[(gy - 1) * G + gx]; w += 1.0; }
          if (gy < G - 1) { s += gridMean[(gy + 1) * G + gx]; w += 1.0; }
          if (gx > 0)     { s += gridMean[gy * G + (gx - 1)]; w += 1.0; }
          if (gx < G - 1) { s += gridMean[gy * G + (gx + 1)]; w += 1.0; }
          next[gy * G + gx] = s / w;
        }
      }
      gridMean = next;
    }
    for (int y = 0; y < S; ++y) {
      double fy = ((y + 0.5) / S) * G - 0.5;
      int gy0 = clampd(std::floor(fy), 0.0, G - 2.0);
      double ty = clampd(fy - gy0, 0.0, 1.0);
      for (int x = 0; x < S; ++x) {
        double fx = ((x + 0.5) / S) * G - 0.5;
        int gx0 = clampd(std::floor(fx), 0.0, G - 2.0);
        double tx = clampd(fx - gx0, 0.0, 1.0);
        double env = (gridMean[gy0 * G + gx0] * (1 - tx) + gridMean[gy0 * G + gx0 + 1] * tx) * (1 - ty) +
                     (gridMean[(gy0 + 1) * G + gx0] * (1 - tx) + gridMean[(gy0 + 1) * G + gx0 + 1] * tx) * ty;
        floorBot.luma[y * S + x] = clampd(floorBot.luma[y * S + x] * (0.42 / std::max(0.15, env)), 0.0, 1.0);
      }
    }
  }

  // 1. Pre-build the leveled 2048x1024 cylindrical wall from wallRaw, with a LOWER EFFECTIVE CAMERA HEIGHT:
  //    - Upper cornice (base of 5-ring dome): y_cornice_tgt = 384 (theta = +22.50 deg)
  //    - Middle entablature (above columns):  y_ent_tgt     = 496 (theta = +2.81 deg)
  //    - Column bases / marble floor edge:    y_base_tgt    = 624 (theta = -19.69 deg)
  const int W_WALL = wallRaw.W;
  const int H_WALL = wallRaw.H;
  Image wallLeveled{W_WALL, H_WALL, std::vector<double>(W_WALL * H_WALL, 0.0)};

  const double y_cornice_tgt = 384.0; // theta = +22.50 deg (base of 5-ring dome)
  const double y_ent_tgt     = 496.0; // theta = +2.81 deg  (entablature above columns)
  const double y_base_tgt    = 624.0; // theta = -19.69 deg (base of columns / floor edge)

  for (int y = 0; y < H_WALL; ++y) {
    for (int x = 0; x < W_WALL; ++x) {
      double u = (x + 0.5) / W_WALL;
      double d2 = (u - 0.5) * (u - 0.5);

      double y_cornice_raw = 526.0 - 488.0 * d2;
      double y_ent_raw     = 638.0 + 115.0 * d2;
      double y_base_raw    = 772.0 + 532.0 * d2;

      double y_raw;
      if (y <= y_ent_tgt) {
        double t = (y - y_cornice_tgt) / (y_ent_tgt - y_cornice_tgt);
        y_raw = y_cornice_raw + t * (y_ent_raw - y_cornice_raw);
      } else if (y <= y_base_tgt) {
        double t = (y - y_ent_tgt) / (y_base_tgt - y_ent_tgt);
        y_raw = y_ent_raw + t * (y_base_raw - y_ent_raw);
      } else {
        double t = (y - y_base_tgt) / (860.0 - y_base_tgt);
        y_raw = y_base_raw + t * (1020.0 - y_base_raw);
      }

      wallLeveled.luma[y * W_WALL + x] = sampleBilinearRepeatU(wallRaw, u, y_raw / H_WALL);
    }
  }

  // Smooth the wrap-around seam at x = 0 / W_WALL - 1 over +-12 pixels
  Image wallSmoothed = wallLeveled;
  const int W_SEAM = 12;
  for (int y = 0; y < H_WALL; ++y) {
    for (int dx = -W_SEAM; dx < W_SEAM; ++dx) {
      int wx = (dx + W_WALL) % W_WALL;
      double t = (dx + W_SEAM + 0.5) / (2.0 * W_SEAM);
      int xL = (W_WALL - 1 - std::max(0, dx));
      int xR = std::max(0, -dx);
      double vOrig = wallLeveled.luma[y * W_WALL + wx];
      double vMirror = (1.0 - t) * wallLeveled.luma[y * W_WALL + xL] + t * wallLeveled.luma[y * W_WALL + xR];
      double wBlend = (1.0 - std::abs(dx + 0.5) / W_SEAM) * 0.75;
      wallSmoothed.luma[y * W_WALL + wx] = vOrig * (1.0 - wBlend) + vMirror * wBlend;
    }
  }

  const double theta_cornice = (0.5 - y_cornice_tgt / H_WALL) * PI; // +0.3927 rad (+22.50 deg)
  const double theta_floor   = (0.5 - y_base_tgt / H_WALL) * PI;    // -0.3436 rad (-19.69 deg)
  const double alpha_cornice = 0.5 * PI - theta_cornice;            // 1.1781 rad (67.50 deg)
  const double u_wall_shift = 0.19;
  const double dome_rot_phi = 0.25;

  // Map the 5-ring hemispherical dome (alpha in [0, alpha_cornice]) to dome_5rings.ppm around its exact center (515.27, 515.52):
  auto evalDomeLuma = [&](double dx, double dy, double dz) -> double {
    double alpha = std::acos(clampd(dy, -1.0, 1.0));
    double phi = std::atan2(dz, dx) + dome_rot_phi;
    double r_tgt = clampd(alpha / alpha_cornice, 0.0, 1.0) * 0.96;

    const double cx = 515.27 / 1024.0;
    const double cy = 515.52 / 1024.0;
    const double c = std::cos(phi);
    const double s = std::sin(phi);

    auto sampleAtR = [&](double r_src, double scale) -> double {
      double r_corr = r_src + 0.008 * clampd((r_src - 0.30) / 0.50, 0.0, 1.0) * (-c);
      double u_dome = cx + 0.5 * r_corr * c;
      double v_dome = cy + 0.5 * r_corr * s;
      return sampleBilinearClamp(domeTop, u_dome, v_dome) * scale;
    };

    if (r_tgt <= 0.74) {
      double r_src = r_tgt * (0.865 / 0.74);
      return sampleAtR(r_src, 1.0);
    } else if (r_tgt <= 0.78) {
      double t = smoothstep(0.74, 0.78, r_tgt);
      double r_a = r_tgt * (0.865 / 0.74);
      double r_b = 0.740 + (r_tgt - 0.74) * ((0.928 - 0.740) / (0.96 - 0.74));
      return sampleAtR(r_a, 1.0) * (1.0 - t) + sampleAtR(r_b, 0.92) * t;
    } else {
      double r_src = 0.740 + (r_tgt - 0.74) * ((0.928 - 0.740) / (0.96 - 0.74));
      return sampleAtR(r_src, 0.92);
    }
  };

  auto evalWallLuma = [&](double dx, double dy, double dz) -> double {
    double phi = std::atan2(dz, dx);
    double theta = std::asin(clampd(dy, -1.0, 1.0));
    double u = phi / (2.0 * PI) + 0.5 - u_wall_shift;
    u = u - std::floor(u);
    double v = 0.5 - theta / PI;
    double rawVal = sampleBilinearRepeatU(wallSmoothed, u, v);
    return 0.54 * std::pow(clampd(rawVal, 0.0, 1.0), 1.18);
  };

  // Rectilinear opus sectile marble floor across the circular floor plane (theta <= theta_floor = -19.69 deg):
  auto evalFloorLuma = [&](double dx, double dy, double dz) -> double {
    double ay = std::max(1e-4, -dy);
    double cotMax = std::cos(theta_floor) / (-std::sin(theta_floor)); // 2.794 at theta_floor = -19.69 deg
    double nx = (dx / ay) / cotMax;
    double ny = (dz / ay) / cotMax;
    double r_floor = std::hypot(nx, ny);

    const double cos45 = 0.70710678;
    double rx = (nx * cos45 - ny * cos45) * 1.45;
    double ry = (nx * cos45 + ny * cos45) * 1.45;

    auto mirrorBetween = [](double t_val, double c0, double c1) -> double {
      double u = (t_val * 0.5 + 0.5 - c0) / (c1 - c0);
      u = std::abs(u);
      int ip = (int)std::floor(u);
      double fp = u - ip;
      double m = (ip & 1) ? (1.0 - fp) : fp;
      return c0 + m * (c1 - c0);
    };

    double u = mirrorBetween(rx, 0.1060, 0.8830);
    double v = mirrorBetween(ry, 0.0980, 0.8790);
    double tileVal = sampleBilinearClamp(floorBot, u, v);

    double vignette = 0.36 - 0.13 * clampd(r_floor * r_floor, 0.0, 1.2);
    return clampd(tileVal * vignette, 0.0, 1.0);
  };

  auto evalPantheonFilmic = [&](double dx, double dy, double dz) -> double {
    double len = std::sqrt(dx * dx + dy * dy + dz * dz);
    dx /= len; dy /= len; dz /= len;
    double theta = std::asin(clampd(dy, -1.0, 1.0));

    // 1. Upper cornice fascia transition: theta in [theta_cornice, theta_cornice + 0.035]
    if (theta >= theta_cornice + 0.035) {
      return evalDomeLuma(dx, dy, dz);
    } else if (theta > theta_cornice) {
      double wDome = smoothstep(theta_cornice, theta_cornice + 0.035, theta);
      return evalWallLuma(dx, dy, dz) * (1.0 - wDome) + evalDomeLuma(dx, dy, dz) * wDome;
    }

    // 2. Column base / floor step transition: theta in [theta_floor - 0.040, theta_floor + 0.015]
    if (theta >= theta_floor + 0.015) {
      return evalWallLuma(dx, dy, dz);
    } else if (theta > theta_floor - 0.040) {
      double wWall = smoothstep(theta_floor - 0.040, theta_floor + 0.015, theta);
      return evalFloorLuma(dx, dy, dz) * (1.0 - wWall) + evalWallLuma(dx, dy, dz) * wWall;
    } else {
      return evalFloorLuma(dx, dy, dz);
    }
  };

  // 3. Render the 1024x512 equirectangular HDR probe (with 4x4 supersampling)
  const int OUT_W = 1024;
  const int OUT_H = 512;
  const int SS = 4;
  std::vector<unsigned char> hdrBytes(OUT_W * OUT_H);
  std::vector<unsigned char> previewBytes(OUT_W * OUT_H);

  for (int py = 0; py < OUT_H; ++py) {
    for (int px = 0; px < OUT_W; ++px) {
      double accumHdr = 0.0;
      double accumDisp = 0.0;
      for (int sy = 0; sy < SS; ++sy) {
        for (int sx = 0; sx < SS; ++sx) {
          double u = (px + (sx + 0.5) / SS) / OUT_W;
          double v = (py + (sy + 0.5) / SS) / OUT_H;
          double phi = (u - 0.5) * 2.0 * PI;
          double theta = (0.5 - v) * PI;
          double cosT = std::cos(theta);
          double dx = cosT * std::cos(phi);
          double dy = std::sin(theta);
          double dz = cosT * std::sin(phi);

          double Y = evalPantheonFilmic(dx, dy, dz);
          double L = clampd(inverseTonemap(Y) * 2.0, 0.0, LOG_MAX_LUM);
          double disp = std::pow(clampd(rangeCompress(L, 0.25), 0.0, 1.0), 1.0 / 2.2);
          accumDisp += disp;
          accumHdr += L;
        }
      }
      double avgDisp = accumDisp / (SS * SS);
      double avgHdr  = accumHdr / (SS * SS);
      hdrBytes[py * OUT_W + px] = (unsigned char)clampd(std::round(lumToByteNorm(avgHdr) * 255.0), 0.0, 255.0);
      previewBytes[py * OUT_W + px] = (unsigned char)clampd(std::round(avgDisp * 255.0), 0.0, 255.0);
    }
  }

  saveGrayPPM(outHdrPath, OUT_W, OUT_H, hdrBytes);
  if (outPreviewPath) {
    saveGrayPPM(outPreviewPath, OUT_W, OUT_H, previewBytes);
  }
  return 0;
}
