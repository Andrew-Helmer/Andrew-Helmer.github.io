(function () {
  "use strict";

  const HDR_WEBP_URL =
    "/post_assets/2026-10-15-demystifying-spherical-harmonics/pantheon-hdr.webp";
  const LOG_EPS = 0.002;
  const LOG_MAX_LUM = 12.15625;
  let hdrPromise = null;

  function loadPantheonHDR() {
    if (!hdrPromise) {
      hdrPromise = new Promise((resolve, reject) => {
        const img = new Image();
        img.decoding = "async";
        img.onload = () => {
          const width = img.naturalWidth;
          const height = img.naturalHeight;
          const offscreen = document.createElement("canvas");
          offscreen.width = width;
          offscreen.height = height;
          const ctx = offscreen.getContext("2d", { willReadFrequently: true });
          ctx.drawImage(img, 0, 0);
          const rgba = ctx.getImageData(0, 0, width, height).data;

          // Build 256-entry inverse log-luminance lookup table
          const lo = Math.log2(LOG_EPS);
          const hi = Math.log2(LOG_MAX_LUM + LOG_EPS);
          const lut = new Float32Array(256);
          for (let i = 1; i < 256; i++) {
            const t = i / 255.0;
            lut[i] = Math.max(0.0, Math.pow(2.0, lo + t * (hi - lo)) - LOG_EPS);
          }

          const numPixels = width * height;
          const lumData = new Float32Array(numPixels);
          for (let i = 0; i < numPixels; i++) {
            lumData[i] = lut[rgba[i * 4]];
          }

          resolve({ width, height, data: lumData });
        };
        img.onerror = () => reject(new Error("Failed to load " + HDR_WEBP_URL));
        img.src = HDR_WEBP_URL;
      });
    }
    return hdrPromise;
  }

  function createAxesSVG() {
    const svgNS = "http://www.w3.org/2000/svg";
    const svg = document.createElementNS(svgNS, "svg");
    svg.setAttribute("class", "axes-gizmo");
    svg.setAttribute("viewBox", "-36 -36 72 72");
    svg.setAttribute("aria-hidden", "true");

    const axes = ["x", "y", "z"].map((label) => {
      const g = document.createElementNS(svgNS, "g");
      const line = document.createElementNS(svgNS, "line");
      line.setAttribute("x1", "0");
      line.setAttribute("y1", "0");
      line.setAttribute("stroke-width", "1.15");
      line.setAttribute("stroke-linecap", "round");

      const arrow = document.createElementNS(svgNS, "polygon");

      const text = document.createElementNS(svgNS, "text");
      text.textContent = label;
      text.setAttribute("text-anchor", "middle");
      text.setAttribute("dominant-baseline", "central");
      text.setAttribute("class", "axis-label");

      g.appendChild(line);
      g.appendChild(arrow);
      g.appendChild(text);
      svg.appendChild(g);
      return { label, g, line, arrow, text };
    });

    const originDot = document.createElementNS(svgNS, "circle");
    originDot.setAttribute("cx", "0");
    originDot.setAttribute("cy", "0");
    originDot.setAttribute("r", "1.6");
    originDot.setAttribute("fill", "#8a8a8a");
    svg.appendChild(originDot);

    function update(yaw, pitch) {
      const cp = Math.cos(pitch),
        sp = Math.sin(pitch);
      const cy = Math.cos(yaw),
        sy = Math.sin(yaw);

      // Camera basis matching shader
      const fx = -sy * cp,
        fy = -sp,
        fz = -cy * cp;
      // camRight = normalize(cross(camForward, worldUp(0,1,0))) = (-fz, 0, fx) / cp
      const rx = cy,
        ry = 0,
        rz = -sy;
      // camUp = cross(camRight, camForward)
      const ux = -sy * sp,
        uy = cp,
        uz = -cy * sp;

      const worldAxes = [
        { item: axes[0], wx: 1, wy: 0, wz: 0 },
        { item: axes[1], wx: 0, wy: 1, wz: 0 },
        { item: axes[2], wx: 0, wy: 0, wz: 1 }
      ];

      const axisLen = 21;
      const labelDist = 29;

      const projected = worldAxes.map(({ item, wx, wy, wz }) => {
        const vx = wx * rx + wy * ry + wz * rz;
        const vy = wx * ux + wy * uy + wz * uz;
        const vz = -(wx * fx + wy * fy + wz * fz); // +vz is toward viewer
        return { item, vx, vy, vz };
      });

      // Depth-sort back-to-front
      projected.sort((a, b) => a.vz - b.vz);

      projected.forEach(({ item, vx, vy, vz }) => {
        svg.appendChild(item.g);
        const sx = vx * axisLen;
        const syScreen = -vy * axisLen; // SVG +Y is down

        const depthFactor = 0.5 * (vz + 1.0); // 0 (away) to 1 (toward)
        const grayVal = Math.round(95 + 75 * depthFactor);
        const labelGray = Math.round(135 + 80 * depthFactor);
        const strokeColor = `rgb(${grayVal}, ${grayVal}, ${grayVal})`;
        const textColor = `rgb(${labelGray}, ${labelGray}, ${labelGray})`;

        item.line.setAttribute("x2", sx.toFixed(2));
        item.line.setAttribute("y2", syScreen.toFixed(2));
        item.line.setAttribute("stroke", strokeColor);

        const len2D = Math.hypot(sx, syScreen);
        if (len2D > 2.5) {
          const ux2 = sx / len2D;
          const uy2 = syScreen / len2D;
          const px2 = -uy2;
          const py2 = ux2;
          const headLen = 4.2;
          const headWidth = 1.8;
          const bx = sx - ux2 * headLen;
          const by = syScreen - uy2 * headLen;
          const p1 = `${sx.toFixed(2)},${syScreen.toFixed(2)}`;
          const p2 = `${(bx + px2 * headWidth).toFixed(2)},${(by + py2 * headWidth).toFixed(2)}`;
          const p3 = `${(bx - px2 * headWidth).toFixed(2)},${(by - py2 * headWidth).toFixed(2)}`;
          item.arrow.setAttribute("points", `${p1} ${p2} ${p3}`);
          item.arrow.setAttribute("fill", strokeColor);
          item.arrow.style.display = "";
        } else {
          item.arrow.style.display = "none";
        }

        const lx = vx * labelDist;
        const ly = -vy * labelDist;
        item.text.setAttribute("x", lx.toFixed(2));
        item.text.setAttribute("y", ly.toFixed(2));
        item.text.setAttribute("fill", textColor);
      });

      svg.appendChild(originDot);
    }

    return { svg, update };
  }

  function initReflectiveSphereDiagram(container) {
    const stage = document.createElement("div");
    stage.className = "diagram-stage";

    const canvas = document.createElement("canvas");
    canvas.className = "diagram-canvas";
    stage.appendChild(canvas);

    const axesGizmo = createAxesSVG();
    stage.appendChild(axesGizmo.svg);

    const controls = document.createElement("div");
    controls.className = "diagram-controls";

    const toggleGroup = document.createElement("div");
    toggleGroup.className = "segmented-toggle";
    toggleGroup.setAttribute("role", "group");
    toggleGroup.setAttribute("aria-label", "Sphere visualization mode");

    const btnProj = document.createElement("button");
    btnProj.type = "button";
    btnProj.className = "toggle-btn active";
    btnProj.textContent = "Spherical Projection";

    const btnChrome = document.createElement("button");
    btnChrome.type = "button";
    btnChrome.className = "toggle-btn";
    btnChrome.textContent = "Chrome Ball";

    toggleGroup.appendChild(btnProj);
    toggleGroup.appendChild(btnChrome);
    controls.appendChild(toggleGroup);

    container.innerHTML = "";
    container.appendChild(stage);
    container.appendChild(controls);

    const gl = canvas.getContext("webgl2", {
      alpha: true,
      premultipliedAlpha: true,
      antialias: true,
      powerPreference: "low-power"
    });

    if (!gl) {
      stage.innerHTML =
        '<p class="webgl-fallback">WebGL2 is required to view this interactive diagram.</p>';
      return;
    }

    gl.getExtension("OES_texture_float_linear");

    const vsSource = `#version 300 es
    out vec2 vUV;
    void main() {
      vec2 pos = vec2((gl_VertexID << 1) & 2, gl_VertexID & 2);
      vUV = pos;
      gl_Position = vec4(pos * 2.0 - 1.0, 0.0, 1.0);
    }`;

    const fsSource = `#version 300 es
    precision highp float;
    in vec2 vUV;
    out vec4 fragColor;

    uniform sampler2D uEnvMap;
    uniform vec2 uResolution;
    uniform float uYaw;
    uniform float uPitch;
    uniform float uDist;
    uniform int uProjMode;

    const float PI = 3.14159265358979323846;
    const float F0 = 0.55;

    vec2 dirToEquirectUV(vec3 dir) {
      vec3 d = normalize(dir);
      float phi = atan(d.z, d.x);
      float theta = asin(clamp(d.y, -1.0, 1.0));
      return vec2(phi / (2.0 * PI) + 0.5, 0.5 - theta / PI);
    }

    float sampleEquirectLum(vec3 dir) {
      return textureLod(uEnvMap, dirToEquirectUV(dir), 0.0).r;
    }

    // Frostbite exponential ("Photographic") luma range compression
    float rangeCompress(float x) {
      return 1.0 - exp(-x);
    }

    float rangeCompress(float val, float threshold) {
      float v1 = val;
      float v2 = threshold + (1.0 - threshold) * rangeCompress((val - threshold) / (1.0 - threshold));
      return val < threshold ? v1 : v2;
    }

    void main() {
      vec2 ndc = vUV * 2.0 - 1.0;
      float aspect = uResolution.x / uResolution.y;
      ndc.x *= aspect;

      const float tanFov = 0.40;

      float cp = cos(uPitch), sp = sin(uPitch);
      float cy = cos(uYaw),   sy = sin(uYaw);
      vec3 camForward = normalize(vec3(-sy * cp, -sp, -cy * cp));
      vec3 worldUp = vec3(0.0, 1.0, 0.0);
      vec3 camRight = normalize(cross(camForward, worldUp));
      vec3 camUp = cross(camRight, camForward);

      vec3 ro = -camForward * uDist;
      vec3 rd = normalize(camForward + ndc.x * tanFov * camRight + ndc.y * tanFov * camUp);

      // Ray-sphere intersection with analytical silhouette anti-aliasing
      float b = dot(ro, rd);
      float c = dot(ro, ro) - 1.0;
      float h = b * b - c;

      // Closest approach distance from sphere center to ray: d^2 = dot(ro,ro) - b^2 = 1.0 - h
      float distFromCenter = sqrt(max(0.0, 1.0 - h));
      float pixelRadius = (2.0 * tanFov * uDist) / uResolution.y;
      float coverage = 1.0 - smoothstep(1.0 - pixelRadius * 0.75, 1.0 + pixelRadius * 0.75, distFromCenter);

      if (coverage <= 0.0) {
        fragColor = vec4(0.0);
        return;
      }

      float clampedH = max(0.0, h);
      float t = -b - sqrt(clampedH);
      vec3 pos = ro + t * rd;
      vec3 n = normalize(pos);

      float hdrLum;
      if (uProjMode == 1) {
        // Spherical projection: sample HDR probe via surface normal
        float rawLum = sampleEquirectLum(n);

        // Subtle view-locked Lambertian shading from upper-right:
        vec3 viewLightDir = normalize(0.55 * camRight + 0.55 * camUp - 0.2 * camForward);
        float ndotl = max(dot(n, viewLightDir), 0.0);
        float shading = 0.3 + 0.7 * ndotl;

        hdrLum = rawLum * shading;
      } else {
        vec3 v = -rd;
        vec3 refl = reflect(rd, n);
        float envLum = sampleEquirectLum(refl);
        float ndotv = clamp(dot(n, v), 0.0, 1.0);
        float fresnel = F0 + (1.0 - F0) * pow(1.0 - ndotv, 5.0);
        hdrLum = envLum * fresnel;
      }

      const float linearSegmentEnd = 0.25;
      float linearVal = rangeCompress(hdrLum, linearSegmentEnd);
      float ldr = pow(clamp(linearVal, 0.0, 1.0), 1.0 / 2.2);
      fragColor = vec4(vec3(ldr) * coverage, coverage);
    }`;

    function compileShader(type, src) {
      const s = gl.createShader(type);
      gl.shaderSource(s, src);
      gl.compileShader(s);
      if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) {
        console.error(gl.getShaderInfoLog(s));
      }
      return s;
    }

    const prog = gl.createProgram();
    gl.attachShader(prog, compileShader(gl.VERTEX_SHADER, vsSource));
    gl.attachShader(prog, compileShader(gl.FRAGMENT_SHADER, fsSource));
    gl.linkProgram(prog);

    const vao = gl.createVertexArray();
    const uLocs = {
      envMap: gl.getUniformLocation(prog, "uEnvMap"),
      resolution: gl.getUniformLocation(prog, "uResolution"),
      yaw: gl.getUniformLocation(prog, "uYaw"),
      pitch: gl.getUniformLocation(prog, "uPitch"),
      dist: gl.getUniformLocation(prog, "uDist"),
      projMode: gl.getUniformLocation(prog, "uProjMode")
    };

    const DEFAULT_YAW = 0.38;
    const DEFAULT_PITCH = 0.32;
    let envTexture = null;
    let yaw = DEFAULT_YAW;
    let pitch = DEFAULT_PITCH;
    const dist = 2.85;
    let projMode = 1; // 0 = Chrome Ball, 1 = Spherical Projection
    let isVisible = false;
    let rafId = null;

    axesGizmo.update(yaw, pitch);

    function renderFrame() {
      rafId = null;
      if (!isVisible || !envTexture) return;

      const dpr = Math.min(window.devicePixelRatio || 1, 2);
      const rect = stage.getBoundingClientRect();
      const w = Math.max(1, Math.round(rect.width * dpr));
      const h = Math.max(1, Math.round(rect.height * dpr));

      if (canvas.width !== w || canvas.height !== h) {
        canvas.width = w;
        canvas.height = h;
      }

      gl.viewport(0, 0, w, h);
      gl.clearColor(0.0, 0.0, 0.0, 0.0);
      gl.clear(gl.COLOR_BUFFER_BIT);

      gl.useProgram(prog);
      gl.bindVertexArray(vao);
      gl.activeTexture(gl.TEXTURE0);
      gl.bindTexture(gl.TEXTURE_2D, envTexture);

      gl.uniform1i(uLocs.envMap, 0);
      gl.uniform2f(uLocs.resolution, w, h);
      gl.uniform1f(uLocs.yaw, yaw);
      gl.uniform1f(uLocs.pitch, pitch);
      gl.uniform1f(uLocs.dist, dist);
      gl.uniform1i(uLocs.projMode, projMode);

      gl.drawArrays(gl.TRIANGLES, 0, 3);
    }

    function requestRender() {
      if (rafId === null && isVisible && envTexture) {
        rafId = requestAnimationFrame(renderFrame);
      }
    }

    // Load shared HDR luminance probe
    loadPantheonHDR()
      .then((hdr) => {
        const t = gl.createTexture();
        gl.bindTexture(gl.TEXTURE_2D, t);
        gl.texImage2D(
          gl.TEXTURE_2D,
          0,
          gl.R16F,
          hdr.width,
          hdr.height,
          0,
          gl.RED,
          gl.FLOAT,
          hdr.data
        );
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.REPEAT);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
        envTexture = t;
        requestRender();
      })
      .catch((err) => {
        console.error("Error loading HDR light probe:", err);
      });

    // Drag to orbit
    let dragging = false;
    let lastX = 0;
    let lastY = 0;

    canvas.addEventListener("pointerdown", (e) => {
      if (e.button !== 0 && e.pointerType === "mouse") return;
      dragging = true;
      lastX = e.clientX;
      lastY = e.clientY;
      canvas.setPointerCapture(e.pointerId);
    });

    canvas.addEventListener("pointermove", (e) => {
      if (!dragging) return;
      const dx = e.clientX - lastX;
      const dy = e.clientY - lastY;
      lastX = e.clientX;
      lastY = e.clientY;

      yaw -= dx * 0.007;
      pitch = Math.max(-1.45, Math.min(1.45, pitch + dy * 0.007));

      axesGizmo.update(yaw, pitch);
      requestRender();
    });

    function endDrag(e) {
      if (!dragging) return;
      dragging = false;
      if (e && e.pointerId !== undefined && canvas.hasPointerCapture(e.pointerId)) {
        canvas.releasePointerCapture(e.pointerId);
      }
    }

    canvas.addEventListener("pointerup", endDrag);
    canvas.addEventListener("pointercancel", endDrag);

    canvas.addEventListener("dblclick", () => {
      yaw = DEFAULT_YAW;
      pitch = DEFAULT_PITCH;
      axesGizmo.update(yaw, pitch);
      requestRender();
    });

    function setMode(mode) {
      if (projMode === mode) return;
      projMode = mode;
      btnChrome.classList.toggle("active", mode === 0);
      btnProj.classList.toggle("active", mode === 1);
      requestRender();
    }

    btnChrome.addEventListener("click", () => setMode(0));
    btnProj.addEventListener("click", () => setMode(1));

    // Only render when on screen or resized
    if ("IntersectionObserver" in window) {
      const io = new IntersectionObserver(
        (entries) => {
          for (const entry of entries) {
            isVisible = entry.isIntersecting;
            if (isVisible) requestRender();
          }
        },
        { threshold: 0.01 }
      );
      io.observe(container);
    } else {
      isVisible = true;
      requestRender();
    }

    if ("ResizeObserver" in window) {
      const ro = new ResizeObserver(() => {
        requestRender();
      });
      ro.observe(stage);
    } else {
      window.addEventListener("resize", requestRender, { passive: true });
    }
  }

  function initLinearApproxDiagram(container) {
    const stage = document.createElement("div");
    stage.className = "diagram-stage";

    const canvas = document.createElement("canvas");
    canvas.className = "diagram-canvas";
    stage.appendChild(canvas);

    const axesGizmo = createAxesSVG();
    stage.appendChild(axesGizmo.svg);

    const slidersWrap = document.createElement("div");
    slidersWrap.className = "diagram-sliders";

    const DEFAULT_C0 = 0.5;
    const dirNorm = Math.hypot(1.0, 1.75, -1.75);
    const DEFAULT_CX = (1.0 / dirNorm) * DEFAULT_C0;
    const DEFAULT_CY = (1.75 / dirNorm) * DEFAULT_C0;
    const DEFAULT_CZ = (-1.75 / dirNorm) * DEFAULT_C0;

    const sliderDefs = [
      { key: "c0", sub: "0", min: 0.0, max: 1.0, def: DEFAULT_C0, centered: false },
      { key: "cx", sub: "x", min: -1.0, max: 1.0, def: DEFAULT_CX, centered: true },
      { key: "cy", sub: "y", min: -1.0, max: 1.0, def: DEFAULT_CY, centered: true },
      { key: "cz", sub: "z", min: -1.0, max: 1.0, def: DEFAULT_CZ, centered: true }
    ];

    const coeffs = {
      c0: DEFAULT_C0,
      cx: DEFAULT_CX,
      cy: DEFAULT_CY,
      cz: DEFAULT_CZ
    };

    let requestRender = () => {};

    sliderDefs.forEach((def) => {
      const row = document.createElement("div");
      row.className = "slider-row";

      const label = document.createElement("span");
      label.className = "slider-label";
      label.innerHTML = `<i>c</i><sub>${def.sub}</sub>`;

      const trackWrap = document.createElement("div");
      trackWrap.className = "slider-track-wrap";

      const track = document.createElement("div");
      track.className = "slider-track";

      const fill = document.createElement("div");
      fill.className = "slider-fill";
      track.appendChild(fill);

      const input = document.createElement("input");
      input.type = "range";
      input.className = "slider-input";
      input.min = String(def.min);
      input.max = String(def.max);
      input.step = "any";
      input.value = String(def.def);
      input.setAttribute("aria-label", `c_${def.sub}`);

      function updateFill(val) {
        if (def.centered) {
          const pct = ((val - def.min) / (def.max - def.min)) * 100;
          const left = Math.min(50, pct);
          const width = Math.abs(pct - 50);
          fill.style.left = `${left}%`;
          fill.style.width = `${width}%`;
        } else {
          const pct = ((val - def.min) / (def.max - def.min)) * 100;
          fill.style.left = "0%";
          fill.style.width = `${pct}%`;
        }
      }

      updateFill(def.def);

      input.addEventListener("input", () => {
        const val = parseFloat(input.value);
        coeffs[def.key] = val;
        updateFill(val);
        requestRender();
      });

      input.addEventListener("dblclick", () => {
        input.value = String(def.def);
        coeffs[def.key] = def.def;
        updateFill(def.def);
        requestRender();
      });

      trackWrap.appendChild(track);
      trackWrap.appendChild(input);
      row.appendChild(label);
      row.appendChild(trackWrap);
      slidersWrap.appendChild(row);
    });

    container.innerHTML = "";
    container.appendChild(stage);
    container.appendChild(slidersWrap);

    const gl = canvas.getContext("webgl2", {
      alpha: true,
      premultipliedAlpha: true,
      antialias: true,
      powerPreference: "low-power"
    });

    if (!gl) {
      stage.innerHTML =
        '<p class="webgl-fallback">WebGL2 is required to view this interactive diagram.</p>';
      return;
    }

    const vsSource = `#version 300 es
    out vec2 vUV;
    void main() {
      vec2 pos = vec2((gl_VertexID << 1) & 2, gl_VertexID & 2);
      vUV = pos;
      gl_Position = vec4(pos * 2.0 - 1.0, 0.0, 1.0);
    }`;

    const fsSource = `#version 300 es
    precision highp float;
    in vec2 vUV;
    out vec4 fragColor;

    uniform vec2 uResolution;
    uniform float uYaw;
    uniform float uPitch;
    uniform float uDist;
    uniform vec4 uCoeffs; // (c0, cx, cy, cz)

    // Exact 2D signed distance to an arrow along +x from x=0 to x=L, centered on y=0.
    float sdArrow2D(vec2 p, float shaftLen, float headLen, float shaftHalfW, float headHalfW) {
      float v = abs(p.y);
      float u = p.x;

      // Rectangular shaft: [0, shaftLen] x [0, shaftHalfW]
      vec2 dShaft = abs(vec2(v, u - 0.5 * shaftLen)) - vec2(shaftHalfW, 0.5 * shaftLen);
      float sdShaft = length(max(dShaft, 0.0)) + min(max(dShaft.x, dShaft.y), 0.0);

      // Triangular head: base at u = shaftLen with half-width headHalfW, tip at u = shaftLen + headLen
      float uh = u - shaftLen;
      vec2 e = vec2(-headHalfW, headLen);
      vec2 w = vec2(v - headHalfW, uh);
      vec2 dSlant = w - e * clamp(dot(w, e) / max(dot(e, e), 1e-6), 0.0, 1.0);
      vec2 dBase = vec2(v - clamp(v, 0.0, headHalfW), uh);
      float s = (uh >= 0.0 && (v * headLen + uh * headHalfW <= headHalfW * headLen)) ? -1.0 : 1.0;
      float sdHead = sqrt(min(dot(dSlant, dSlant), dot(dBase, dBase))) * s;

      return min(sdShaft, sdHead);
    }

    void main() {
      vec2 ndc = vUV * 2.0 - 1.0;
      float aspect = uResolution.x / uResolution.y;
      ndc.x *= aspect;

      const float tanFov = 0.40;

      float cp = cos(uPitch), sp = sin(uPitch);
      float cy = cos(uYaw),   sy = sin(uYaw);
      vec3 camForward = normalize(vec3(-sy * cp, -sp, -cy * cp));
      vec3 worldUp = vec3(0.0, 1.0, 0.0);
      vec3 camRight = normalize(cross(camForward, worldUp));
      vec3 camUp = cross(camRight, camForward);

      vec3 ro = -camForward * uDist;
      vec3 rd = normalize(camForward + ndc.x * tanFov * camRight + ndc.y * tanFov * camUp);

      float b = dot(ro, rd);
      float c = dot(ro, ro) - 1.0;
      float h = b * b - c;

      float distFromCenter = sqrt(max(0.0, 1.0 - h));
      float pixelRadius = (2.0 * tanFov * uDist) / uResolution.y;
      float sphereCoverage = 1.0 - smoothstep(1.0 - pixelRadius * 0.75, 1.0 + pixelRadius * 0.75, distFromCenter);

      vec3 sphereCol = vec3(0.0);
      if (sphereCoverage > 0.0) {
        float clampedH = max(0.0, h);
        float tSphere = -b - sqrt(clampedH);
        vec3 pos = ro + tSphere * rd;
        vec3 n = normalize(pos);

        // Evaluate L1 linear model: f_hat(w) = c0 + cx*w_x + cy*w_y + cz*w_z
        float val = uCoeffs.x + dot(uCoeffs.yzw, n);
        float ldr = clamp(abs(val), 0.0, 1.0);

        sphereCol = (val >= 0.0)
          ? vec3(ldr)
          : vec3(ldr, ldr * 0.22, ldr * 0.18);
      }

      // Render 3D vector arrow for (cx, cy, cz) with a clean, flat 2D diagram style.
      // Since the unit sphere represents radius 0.5, scale vector by 2.0 so length 0.5 -> 1.0.
      vec3 cVec = uCoeffs.yzw * 2.0;
      float vecLen = length(cVec);
      vec3 arrowCol = vec3(0.0);
      float arrowAlpha = 0.0;

      if (vecLen > 0.005) {
        vec3 uDir = cVec / vecLen;

        // Project 3D arrow endpoint cVec into 2D ndc space (origin projects to (0, 0))
        vec3 vTip = cVec - ro;
        float z0 = uDist;
        float z1 = dot(vTip, camForward);
        vec2 s1 = vec2(dot(vTip, camRight), dot(vTip, camUp)) / (z1 * tanFov);
        float L2D = length(s1);

        float scale3D = clamp(vecLen / 0.10, 0.0, 1.0);
        float shaftHalfW = 0.0075 * scale3D;
        float headHalfW = 0.023 * scale3D;

        if (L2D > 1e-4) {
          vec2 dir2D = s1 / L2D;
          vec2 pLocal = vec2(
            dot(ndc, dir2D),
            ndc.x * dir2D.y - ndc.y * dir2D.x
          );

          float foreshorten = clamp(L2D / max(vecLen * 0.75, 1e-4), 0.25, 1.0);
          float headLen = min(0.068 * scale3D * foreshorten, 0.65 * L2D);
          float shaftLen = max(0.0, L2D - headLen);

          float sd = sdArrow2D(pLocal, shaftLen, headLen, shaftHalfW, headHalfW);
          float aaNdc = (2.0 / uResolution.y) * 0.75;
          float arrowCoverage = 1.0 - smoothstep(-aaNdc, aaNdc, sd);

          if (arrowCoverage > 0.0) {
            // Perspective-correct fractional position alpha in [0, 1] along the 3D vector
            float tau = clamp(pLocal.x / L2D, 0.0, 1.0);
            float alpha = (tau * z1) / max(z0 + tau * (z1 - z0), 1e-4);
            float rArrow = alpha * vecLen;

            // A point on the arrow is in front of the sphere only if it extends past
            // the unit sphere radius (rArrow > 1.0) on the camera-facing hemisphere.
            float aa3D = pixelRadius * 1.25;
            float outsideRadius = smoothstep(1.0, 1.0 + aa3D, rArrow);
            float inFrontHemisphere = smoothstep(-0.02, 0.02, -dot(uDir, camForward));
            float inFrontOfSphere = outsideRadius * inFrontHemisphere;

            // Any part behind the sphere's front surface within the sphere's silhouette is occluded
            float occluded = sphereCoverage * (1.0 - inFrontOfSphere);

            vec3 occludedCol = vec3(0.48);
            float occludedAlpha = 0.20;

            vec3 unoccludedCol = vec3(0.62);
            float unoccludedAlpha = 1.0;

            arrowCol = mix(unoccludedCol, occludedCol, occluded);
            arrowAlpha = arrowCoverage * mix(unoccludedAlpha, occludedAlpha, occluded);
          }
        }
      }

      // Composite arrow over sphere (premultiplied alpha)
      vec4 dst = vec4(sphereCol * sphereCoverage, sphereCoverage);
      vec4 src = vec4(arrowCol * arrowAlpha, arrowAlpha);
      fragColor = src + dst * (1.0 - src.a);
    }`;

    function compileShader(type, src) {
      const s = gl.createShader(type);
      gl.shaderSource(s, src);
      gl.compileShader(s);
      if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) {
        console.error(gl.getShaderInfoLog(s));
      }
      return s;
    }

    const prog = gl.createProgram();
    gl.attachShader(prog, compileShader(gl.VERTEX_SHADER, vsSource));
    gl.attachShader(prog, compileShader(gl.FRAGMENT_SHADER, fsSource));
    gl.linkProgram(prog);

    const vao = gl.createVertexArray();
    const uLocs = {
      resolution: gl.getUniformLocation(prog, "uResolution"),
      yaw: gl.getUniformLocation(prog, "uYaw"),
      pitch: gl.getUniformLocation(prog, "uPitch"),
      dist: gl.getUniformLocation(prog, "uDist"),
      coeffs: gl.getUniformLocation(prog, "uCoeffs")
    };

    const DEFAULT_YAW = 0.38;
    const DEFAULT_PITCH = 0.32;
    let yaw = DEFAULT_YAW;
    let pitch = DEFAULT_PITCH;
    const dist = 2.85;
    let isVisible = false;
    let rafId = null;

    axesGizmo.update(yaw, pitch);

    function renderFrame() {
      rafId = null;
      if (!isVisible) return;

      const dpr = Math.min(window.devicePixelRatio || 1, 2);
      const rect = stage.getBoundingClientRect();
      const w = Math.max(1, Math.round(rect.width * dpr));
      const h = Math.max(1, Math.round(rect.height * dpr));

      if (canvas.width !== w || canvas.height !== h) {
        canvas.width = w;
        canvas.height = h;
      }

      gl.viewport(0, 0, w, h);
      gl.clearColor(0.0, 0.0, 0.0, 0.0);
      gl.clear(gl.COLOR_BUFFER_BIT);

      gl.useProgram(prog);
      gl.bindVertexArray(vao);

      gl.uniform2f(uLocs.resolution, w, h);
      gl.uniform1f(uLocs.yaw, yaw);
      gl.uniform1f(uLocs.pitch, pitch);
      gl.uniform1f(uLocs.dist, dist);
      gl.uniform4f(uLocs.coeffs, coeffs.c0, coeffs.cx, coeffs.cy, coeffs.cz);

      gl.drawArrays(gl.TRIANGLES, 0, 3);
    }

    requestRender = function () {
      if (rafId === null && isVisible) {
        rafId = requestAnimationFrame(renderFrame);
      }
    };

    let dragging = false;
    let lastX = 0;
    let lastY = 0;

    canvas.addEventListener("pointerdown", (e) => {
      if (e.button !== 0 && e.pointerType === "mouse") return;
      dragging = true;
      lastX = e.clientX;
      lastY = e.clientY;
      canvas.setPointerCapture(e.pointerId);
    });

    canvas.addEventListener("pointermove", (e) => {
      if (!dragging) return;
      const dx = e.clientX - lastX;
      const dy = e.clientY - lastY;
      lastX = e.clientX;
      lastY = e.clientY;

      yaw -= dx * 0.007;
      pitch = Math.max(-1.45, Math.min(1.45, pitch + dy * 0.007));

      axesGizmo.update(yaw, pitch);
      requestRender();
    });

    function endDrag(e) {
      if (!dragging) return;
      dragging = false;
      if (e && e.pointerId !== undefined && canvas.hasPointerCapture(e.pointerId)) {
        canvas.releasePointerCapture(e.pointerId);
      }
    }

    canvas.addEventListener("pointerup", endDrag);
    canvas.addEventListener("pointercancel", endDrag);

    canvas.addEventListener("dblclick", () => {
      yaw = DEFAULT_YAW;
      pitch = DEFAULT_PITCH;
      axesGizmo.update(yaw, pitch);
      requestRender();
    });

    if ("IntersectionObserver" in window) {
      const io = new IntersectionObserver(
        (entries) => {
          for (const entry of entries) {
            isVisible = entry.isIntersecting;
            if (isVisible) requestRender();
          }
        },
        { threshold: 0.01 }
      );
      io.observe(container);
    } else {
      isVisible = true;
      requestRender();
    }

    if ("ResizeObserver" in window) {
      const ro = new ResizeObserver(() => {
        requestRender();
      });
      ro.observe(stage);
    } else {
      window.addEventListener("resize", requestRender, { passive: true });
    }
  }

  function initPantheonFitDiagram(container) {
    const stage = document.createElement("div");
    stage.className = "diagram-stage";

    const canvas = document.createElement("canvas");
    canvas.className = "diagram-canvas";
    stage.appendChild(canvas);

    const axesGizmo = createAxesSVG();
    stage.appendChild(axesGizmo.svg);

    const controls = document.createElement("div");
    controls.className = "diagram-controls";

    const toggleGroup = document.createElement("div");
    toggleGroup.className = "segmented-toggle";
    toggleGroup.setAttribute("role", "group");
    toggleGroup.setAttribute("aria-label", "Spherical harmonics comparison view");

    const btnProj = document.createElement("button");
    btnProj.type = "button";
    btnProj.className = "toggle-btn active";
    btnProj.textContent = "Original";

    const btnSH = document.createElement("button");
    btnSH.type = "button";
    btnSH.className = "toggle-btn";
    btnSH.textContent = "Approximation";

    const btnErr = document.createElement("button");
    btnErr.type = "button";
    btnErr.className = "toggle-btn";
    btnErr.textContent = "Error";

    toggleGroup.appendChild(btnProj);
    toggleGroup.appendChild(btnSH);
    toggleGroup.appendChild(btnErr);
    controls.appendChild(toggleGroup);

    container.innerHTML = "";
    container.appendChild(stage);
    container.appendChild(controls);

    const gl = canvas.getContext("webgl2", {
      alpha: true,
      premultipliedAlpha: true,
      antialias: true,
      powerPreference: "low-power"
    });

    if (!gl) {
      stage.innerHTML =
        '<p class="webgl-fallback">WebGL2 is required to view this interactive diagram.</p>';
      return;
    }

    gl.getExtension("OES_texture_float_linear");

    const vsSource = `#version 300 es
    out vec2 vUV;
    void main() {
      vec2 pos = vec2((gl_VertexID << 1) & 2, gl_VertexID & 2);
      vUV = pos;
      gl_Position = vec4(pos * 2.0 - 1.0, 0.0, 1.0);
    }`;

    const fsSource = `#version 300 es
    precision highp float;
    in vec2 vUV;
    out vec4 fragColor;

    uniform sampler2D uEnvMap;
    uniform vec2 uResolution;
    uniform float uYaw;
    uniform float uPitch;
    uniform float uDist;
    uniform int uViewMode; // 0 = Spherical Projection, 1 = SH Approximation, 2 = Absolute Error
    uniform vec4 uCoeffs;  // (c0, cx, cy, cz)

    const float PI = 3.14159265358979323846;

    vec2 dirToEquirectUV(vec3 dir) {
      vec3 d = normalize(dir);
      float phi = atan(d.z, d.x);
      float theta = asin(clamp(d.y, -1.0, 1.0));
      return vec2(phi / (2.0 * PI) + 0.5, 0.5 - theta / PI);
    }

    float sampleEquirectLum(vec3 dir) {
      return textureLod(uEnvMap, dirToEquirectUV(dir), 0.0).r;
    }

    // Frostbite exponential ("Photographic") luma range compression
    float rangeCompress(float x) {
      return 1.0 - exp(-x);
    }

    float rangeCompress(float val, float threshold) {
      float v1 = val;
      float v2 = threshold + (1.0 - threshold) * rangeCompress((val - threshold) / (1.0 - threshold));
      return val < threshold ? v1 : v2;
    }

    void main() {
      vec2 ndc = vUV * 2.0 - 1.0;
      float aspect = uResolution.x / uResolution.y;
      ndc.x *= aspect;

      const float tanFov = 0.40;

      float cp = cos(uPitch), sp = sin(uPitch);
      float cy = cos(uYaw),   sy = sin(uYaw);
      vec3 camForward = normalize(vec3(-sy * cp, -sp, -cy * cp));
      vec3 worldUp = vec3(0.0, 1.0, 0.0);
      vec3 camRight = normalize(cross(camForward, worldUp));
      vec3 camUp = cross(camRight, camForward);

      vec3 ro = -camForward * uDist;
      vec3 rd = normalize(camForward + ndc.x * tanFov * camRight + ndc.y * tanFov * camUp);

      float b = dot(ro, rd);
      float c = dot(ro, ro) - 1.0;
      float h = b * b - c;

      float distFromCenter = sqrt(max(0.0, 1.0 - h));
      float pixelRadius = (2.0 * tanFov * uDist) / uResolution.y;
      float coverage = 1.0 - smoothstep(1.0 - pixelRadius * 0.75, 1.0 + pixelRadius * 0.75, distFromCenter);

      if (coverage <= 0.0) {
        fragColor = vec4(0.0);
        return;
      }

      float clampedH = max(0.0, h);
      float t = -b - sqrt(clampedH);
      vec3 pos = ro + t * rd;
      vec3 n = normalize(pos);

      float rawLum = sampleEquirectLum(n);
      float shRaw = uCoeffs.x + dot(uCoeffs.yzw, n);

      vec3 col;
      if (uViewMode == 0) {
        // Pure spherical projection (no Lambertian shading)
        float mapped = rangeCompress(rawLum, 0.25);
        float ldr = pow(clamp(mapped, 0.0, 1.0), 1.0 / 2.2);
        col = vec3(ldr);
      } else if (uViewMode == 1) {
        // L1 SH approximation with negative values clamped to black
        float mapped = rangeCompress(max(0.0, shRaw), 0.25);
        float ldr = pow(clamp(mapped, 0.0, 1.0), 1.0 / 2.2);
        col = vec3(ldr);
      } else {
        // Absolute error |f(w) - f_hat(w)| compressed with Frostbite rangeCompress -> red
        float absErr = abs(rawLum - shRaw);
        float mapped = rangeCompress(absErr, 0.25);
        float ldr = pow(clamp(mapped, 0.0, 1.0), 1.0 / 2.2);
        col = vec3(ldr, ldr * 0.16, ldr * 0.12);
      }

      fragColor = vec4(col * coverage, coverage);
    }`;

    function compileShader(type, src) {
      const s = gl.createShader(type);
      gl.shaderSource(s, src);
      gl.compileShader(s);
      if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) {
        console.error(gl.getShaderInfoLog(s));
      }
      return s;
    }

    const prog = gl.createProgram();
    gl.attachShader(prog, compileShader(gl.VERTEX_SHADER, vsSource));
    gl.attachShader(prog, compileShader(gl.FRAGMENT_SHADER, fsSource));
    gl.linkProgram(prog);

    const vao = gl.createVertexArray();
    const uLocs = {
      envMap: gl.getUniformLocation(prog, "uEnvMap"),
      resolution: gl.getUniformLocation(prog, "uResolution"),
      yaw: gl.getUniformLocation(prog, "uYaw"),
      pitch: gl.getUniformLocation(prog, "uPitch"),
      dist: gl.getUniformLocation(prog, "uDist"),
      viewMode: gl.getUniformLocation(prog, "uViewMode"),
      coeffs: gl.getUniformLocation(prog, "uCoeffs")
    };

    // Precomputed exact L1 SH projection of the Pantheon HDR probe
    const PANTHEON_L1_COEFFS = [0.317151, 0.064818, 0.355999, -0.029964];

    const DEFAULT_YAW = 0.38;
    const DEFAULT_PITCH = 0.32;
    let envTexture = null;
    let yaw = DEFAULT_YAW;
    let pitch = DEFAULT_PITCH;
    const dist = 2.85;
    let viewMode = 0; // 0 = Spherical Projection, 1 = SH Approximation, 2 = Absolute Error
    let isVisible = false;
    let rafId = null;

    axesGizmo.update(yaw, pitch);

    function renderFrame() {
      rafId = null;
      if (!isVisible || !envTexture) return;

      const dpr = Math.min(window.devicePixelRatio || 1, 2);
      const rect = stage.getBoundingClientRect();
      const w = Math.max(1, Math.round(rect.width * dpr));
      const h = Math.max(1, Math.round(rect.height * dpr));

      if (canvas.width !== w || canvas.height !== h) {
        canvas.width = w;
        canvas.height = h;
      }

      gl.viewport(0, 0, w, h);
      gl.clearColor(0.0, 0.0, 0.0, 0.0);
      gl.clear(gl.COLOR_BUFFER_BIT);

      gl.useProgram(prog);
      gl.bindVertexArray(vao);
      gl.activeTexture(gl.TEXTURE0);
      gl.bindTexture(gl.TEXTURE_2D, envTexture);

      gl.uniform1i(uLocs.envMap, 0);
      gl.uniform2f(uLocs.resolution, w, h);
      gl.uniform1f(uLocs.yaw, yaw);
      gl.uniform1f(uLocs.pitch, pitch);
      gl.uniform1f(uLocs.dist, dist);
      gl.uniform1i(uLocs.viewMode, viewMode);
      gl.uniform4f(
        uLocs.coeffs,
        PANTHEON_L1_COEFFS[0],
        PANTHEON_L1_COEFFS[1],
        PANTHEON_L1_COEFFS[2],
        PANTHEON_L1_COEFFS[3]
      );

      gl.drawArrays(gl.TRIANGLES, 0, 3);
    }

    function requestRender() {
      if (rafId === null && isVisible && envTexture) {
        rafId = requestAnimationFrame(renderFrame);
      }
    }

    loadPantheonHDR()
      .then((hdr) => {
        const t = gl.createTexture();
        gl.bindTexture(gl.TEXTURE_2D, t);
        gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
        gl.texImage2D(
          gl.TEXTURE_2D,
          0,
          gl.R16F,
          hdr.width,
          hdr.height,
          0,
          gl.RED,
          gl.FLOAT,
          hdr.data
        );
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.REPEAT);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
        envTexture = t;
        requestRender();
      })
      .catch((err) => {
        console.error("Error loading HDR light probe:", err);
      });

    let dragging = false;
    let lastX = 0;
    let lastY = 0;

    canvas.addEventListener("pointerdown", (e) => {
      if (e.button !== 0 && e.pointerType === "mouse") return;
      dragging = true;
      lastX = e.clientX;
      lastY = e.clientY;
      canvas.setPointerCapture(e.pointerId);
    });

    canvas.addEventListener("pointermove", (e) => {
      if (!dragging) return;
      const dx = e.clientX - lastX;
      const dy = e.clientY - lastY;
      lastX = e.clientX;
      lastY = e.clientY;

      yaw -= dx * 0.007;
      pitch = Math.max(-1.45, Math.min(1.45, pitch + dy * 0.007));

      axesGizmo.update(yaw, pitch);
      requestRender();
    });

    function endDrag(e) {
      if (!dragging) return;
      dragging = false;
      if (e && e.pointerId !== undefined && canvas.hasPointerCapture(e.pointerId)) {
        canvas.releasePointerCapture(e.pointerId);
      }
    }

    canvas.addEventListener("pointerup", endDrag);
    canvas.addEventListener("pointercancel", endDrag);

    canvas.addEventListener("dblclick", () => {
      yaw = DEFAULT_YAW;
      pitch = DEFAULT_PITCH;
      axesGizmo.update(yaw, pitch);
      requestRender();
    });

    function setViewMode(mode) {
      if (viewMode === mode) return;
      viewMode = mode;
      btnProj.classList.toggle("active", mode === 0);
      btnSH.classList.toggle("active", mode === 1);
      btnErr.classList.toggle("active", mode === 2);
      requestRender();
    }

    btnProj.addEventListener("click", () => setViewMode(0));
    btnSH.addEventListener("click", () => setViewMode(1));
    btnErr.addEventListener("click", () => setViewMode(2));

    if ("IntersectionObserver" in window) {
      const io = new IntersectionObserver(
        (entries) => {
          for (const entry of entries) {
            isVisible = entry.isIntersecting;
            if (isVisible) requestRender();
          }
        },
        { threshold: 0.01 }
      );
      io.observe(container);
    } else {
      isVisible = true;
      requestRender();
    }

    if ("ResizeObserver" in window) {
      const ro = new ResizeObserver(() => {
        requestRender();
      });
      ro.observe(stage);
    } else {
      window.addEventListener("resize", requestRender, { passive: true });
    }
  }

  function initClampedDirectionalDiagram(container) {
    const stage = document.createElement("div");
    stage.className = "diagram-stage";

    const canvas = document.createElement("canvas");
    canvas.className = "diagram-canvas";
    stage.appendChild(canvas);

    const axesGizmo = createAxesSVG();
    stage.appendChild(axesGizmo.svg);

    const controls = document.createElement("div");
    controls.className = "diagram-controls";

    const toggleGroup = document.createElement("div");
    toggleGroup.className = "segmented-toggle";
    toggleGroup.setAttribute("role", "group");
    toggleGroup.setAttribute("aria-label", "Directional light SH projection comparison");

    const btnNaive = document.createElement("button");
    btnNaive.type = "button";
    btnNaive.className = "toggle-btn";
    btnNaive.textContent = "Naive";

    const btnSplit = document.createElement("button");
    btnSplit.type = "button";
    btnSplit.className = "toggle-btn active";
    btnSplit.textContent = "Split";

    const btnCorrect = document.createElement("button");
    btnCorrect.type = "button";
    btnCorrect.className = "toggle-btn";
    btnCorrect.textContent = "Correct";

    toggleGroup.appendChild(btnNaive);
    toggleGroup.appendChild(btnSplit);
    toggleGroup.appendChild(btnCorrect);
    controls.appendChild(toggleGroup);

    container.innerHTML = "";
    container.appendChild(stage);
    container.appendChild(controls);

    const gl = canvas.getContext("webgl2", {
      alpha: true,
      premultipliedAlpha: true,
      antialias: true,
      powerPreference: "low-power"
    });

    if (!gl) {
      stage.innerHTML =
        '<p class="webgl-fallback">WebGL2 is required to view this interactive diagram.</p>';
      return;
    }

    const vsSource = `#version 300 es
    out vec2 vUV;
    void main() {
      vec2 pos = vec2((gl_VertexID << 1) & 2, gl_VertexID & 2);
      vUV = pos;
      gl_Position = vec4(pos * 2.0 - 1.0, 0.0, 1.0);
    }`;

    const fsSource = `#version 300 es
    precision highp float;
    in vec2 vUV;
    out vec4 fragColor;

    uniform vec2 uResolution;
    uniform float uYaw;
    uniform float uPitch;
    uniform float uDist;
    uniform int uViewMode; // 0 = Naive, 1 = Split, 2 = Correct
    uniform vec3 uLightDir;

    // Exact 2D signed distance to an arrow along +x from x=0 to x=shaftLen+headLen, centered on y=0.
    float sdArrow2D(vec2 p, float shaftLen, float headLen, float shaftHalfW, float headHalfW) {
      float v = abs(p.y);
      float u = p.x;

      vec2 dShaft = abs(vec2(v, u - 0.5 * shaftLen)) - vec2(shaftHalfW, 0.5 * shaftLen);
      float sdShaft = length(max(dShaft, 0.0)) + min(max(dShaft.x, dShaft.y), 0.0);

      float uh = u - shaftLen;
      vec2 e = vec2(-headHalfW, headLen);
      vec2 w = vec2(v - headHalfW, uh);
      vec2 dSlant = w - e * clamp(dot(w, e) / max(dot(e, e), 1e-6), 0.0, 1.0);
      vec2 dBase = vec2(v - clamp(v, 0.0, headHalfW), uh);
      float s = (uh >= 0.0 && (v * headLen + uh * headHalfW <= headHalfW * headLen)) ? -1.0 : 1.0;
      float sdHead = sqrt(min(dot(dSlant, dSlant), dot(dBase, dBase))) * s;

      return min(sdShaft, sdHead);
    }

    void main() {
      vec2 ndc = vUV * 2.0 - 1.0;
      float aspect = uResolution.x / uResolution.y;
      ndc.x *= aspect;

      const float tanFov = 0.40;

      float cp = cos(uPitch), sp = sin(uPitch);
      float cy = cos(uYaw),   sy = sin(uYaw);
      vec3 camForward = normalize(vec3(-sy * cp, -sp, -cy * cp));
      vec3 worldUp = vec3(0.0, 1.0, 0.0);
      vec3 camRight = normalize(cross(camForward, worldUp));
      vec3 camUp = cross(camRight, camForward);

      vec3 ro = -camForward * uDist;
      vec3 rd = normalize(camForward + ndc.x * tanFov * camRight + ndc.y * tanFov * camUp);

      float b = dot(ro, rd);
      float c = dot(ro, ro) - 1.0;
      float h = b * b - c;

      float distFromCenter = sqrt(max(0.0, 1.0 - h));
      float pixelRadius = (2.0 * tanFov * uDist) / uResolution.y;
      float sphereCoverage = 1.0 - smoothstep(1.0 - pixelRadius * 0.75, 1.0 + pixelRadius * 0.75, distFromCenter);

      vec3 sphereCol = vec3(0.0);
      if (sphereCoverage > 0.0) {
        float clampedH = max(0.0, h);
        float t = -b - sqrt(clampedH);
        vec3 pos = ro + t * rd;
        vec3 n = normalize(pos);

        float ndotl = dot(n, uLightDir);

        // With I = 0.25 so that the peak across both views is 1.0:
        // - Naive L1 projection + Lambertian convolution: c0 = 0.25, c_xyz = 0.50 * L -> max(0.25 + 0.50 * (n . L), 0)
        // - Exact Clamped L1 representation:              c0 = 0.00, c_xyz = 1.00 * L -> max(n . L, 0)
        float naiveVal = max(0.25 + 0.50 * ndotl, 0.0);
        float clampedVal = max(ndotl, 0.0);

        float linearVal;
        if (uViewMode == 0) {
          linearVal = naiveVal;
        } else if (uViewMode == 2) {
          linearVal = clampedVal;
        } else {
          linearVal = (ndc.x < 0.0) ? naiveVal : clampedVal;
        }

        float ldr = pow(clamp(linearVal, 0.0, 1.0), 1.0 / 2.2);
        sphereCol = vec3(ldr);
      }

      // Incoming directional light arrow outside the sphere, pointing towards the sphere center
      // from pStart (r = 1.42) to pEnd (r = 1.05) along uLightDir.
      vec3 pStart = uLightDir * 1.42;
      vec3 pEnd   = uLightDir * 1.05;

      vec3 v0 = pStart - ro;
      vec3 v1 = pEnd - ro;
      float z0 = dot(v0, camForward);
      float z1 = dot(v1, camForward);

      vec2 s0 = vec2(dot(v0, camRight), dot(v0, camUp)) / (z0 * tanFov);
      vec2 s1 = vec2(dot(v1, camRight), dot(v1, camUp)) / (z1 * tanFov);

      vec2 seg = s1 - s0;
      float L2D = length(seg);
      vec3 arrowCol = vec3(0.0);
      float arrowAlpha = 0.0;

      if (L2D > 1e-4) {
        vec2 dir2D = seg / L2D;
        vec2 rel = ndc - s0;
        vec2 pLocal = vec2(
          dot(rel, dir2D),
          rel.x * dir2D.y - rel.y * dir2D.x
        );

        float foreshorten = clamp(L2D / (0.37 * 0.75), 0.25, 1.0);
        float shaftHalfW = 0.0075;
        float headHalfW = 0.023;
        float headLen = min(0.068 * foreshorten, 0.65 * L2D);
        float shaftLen = max(0.0, L2D - headLen);

        float sd = sdArrow2D(pLocal, shaftLen, headLen, shaftHalfW, headHalfW);
        float aaNdc = (2.0 / uResolution.y) * 0.75;
        float arrowCoverage = 1.0 - smoothstep(-aaNdc, aaNdc, sd);

        if (arrowCoverage > 0.0) {
          // Since the entire arrow is outside the sphere (r >= 1.05), it is only occluded
          // when it lies on the back hemisphere (dot(uLightDir, camForward) > 0) behind the sphere disc.
          float inFrontHemisphere = smoothstep(-0.02, 0.02, -dot(uLightDir, camForward));
          float occluded = sphereCoverage * (1.0 - inFrontHemisphere);

          arrowCol = vec3(0.92);
          arrowAlpha = arrowCoverage * mix(1.0, 0.20, occluded);
        }
      }

      vec4 dst = vec4(sphereCol * sphereCoverage, sphereCoverage);
      vec4 src = vec4(arrowCol * arrowAlpha, arrowAlpha);
      fragColor = src + dst * (1.0 - src.a);
    }`;

    function compileShader(type, src) {
      const s = gl.createShader(type);
      gl.shaderSource(s, src);
      gl.compileShader(s);
      if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) {
        console.error(gl.getShaderInfoLog(s));
      }
      return s;
    }

    const prog = gl.createProgram();
    gl.attachShader(prog, compileShader(gl.VERTEX_SHADER, vsSource));
    gl.attachShader(prog, compileShader(gl.FRAGMENT_SHADER, fsSource));
    gl.linkProgram(prog);

    const vao = gl.createVertexArray();
    const uLocs = {
      resolution: gl.getUniformLocation(prog, "uResolution"),
      yaw: gl.getUniformLocation(prog, "uYaw"),
      pitch: gl.getUniformLocation(prog, "uPitch"),
      dist: gl.getUniformLocation(prog, "uDist"),
      viewMode: gl.getUniformLocation(prog, "uViewMode"),
      lightDir: gl.getUniformLocation(prog, "uLightDir")
    };

    const dirNorm = Math.hypot(0.705, 0.708, -0.033);
    const LIGHT_DIR = [0.705 / dirNorm, 0.708 / dirNorm, -0.033 / dirNorm];

    const DEFAULT_YAW = 0.38;
    const DEFAULT_PITCH = 0.32;
    let yaw = DEFAULT_YAW;
    let pitch = DEFAULT_PITCH;
    const dist = 2.85;
    let viewMode = 1; // 0 = Naive, 1 = Split, 2 = Correct
    let isVisible = false;
    let rafId = null;

    axesGizmo.update(yaw, pitch);

    function renderFrame() {
      rafId = null;
      if (!isVisible) return;

      const dpr = Math.min(window.devicePixelRatio || 1, 2);
      const rect = stage.getBoundingClientRect();
      const w = Math.max(1, Math.round(rect.width * dpr));
      const h = Math.max(1, Math.round(rect.height * dpr));

      if (canvas.width !== w || canvas.height !== h) {
        canvas.width = w;
        canvas.height = h;
      }

      gl.viewport(0, 0, w, h);
      gl.clearColor(0.0, 0.0, 0.0, 0.0);
      gl.clear(gl.COLOR_BUFFER_BIT);

      gl.useProgram(prog);
      gl.bindVertexArray(vao);

      gl.uniform2f(uLocs.resolution, w, h);
      gl.uniform1f(uLocs.yaw, yaw);
      gl.uniform1f(uLocs.pitch, pitch);
      gl.uniform1f(uLocs.dist, dist);
      gl.uniform1i(uLocs.viewMode, viewMode);
      gl.uniform3f(uLocs.lightDir, LIGHT_DIR[0], LIGHT_DIR[1], LIGHT_DIR[2]);

      gl.drawArrays(gl.TRIANGLES, 0, 3);
    }

    function requestRender() {
      if (rafId === null && isVisible) {
        rafId = requestAnimationFrame(renderFrame);
      }
    }

    let dragging = false;
    let lastX = 0;
    let lastY = 0;

    canvas.addEventListener("pointerdown", (e) => {
      if (e.button !== 0 && e.pointerType === "mouse") return;
      dragging = true;
      lastX = e.clientX;
      lastY = e.clientY;
      canvas.setPointerCapture(e.pointerId);
    });

    canvas.addEventListener("pointermove", (e) => {
      if (!dragging) return;
      const dx = e.clientX - lastX;
      const dy = e.clientY - lastY;
      lastX = e.clientX;
      lastY = e.clientY;

      yaw -= dx * 0.007;
      pitch = Math.max(-1.45, Math.min(1.45, pitch + dy * 0.007));

      axesGizmo.update(yaw, pitch);
      requestRender();
    });

    function endDrag(e) {
      if (!dragging) return;
      dragging = false;
      if (e && e.pointerId !== undefined && canvas.hasPointerCapture(e.pointerId)) {
        canvas.releasePointerCapture(e.pointerId);
      }
    }

    canvas.addEventListener("pointerup", endDrag);
    canvas.addEventListener("pointercancel", endDrag);

    canvas.addEventListener("dblclick", () => {
      yaw = DEFAULT_YAW;
      pitch = DEFAULT_PITCH;
      axesGizmo.update(yaw, pitch);
      requestRender();
    });

    function setViewMode(mode) {
      if (viewMode === mode) return;
      viewMode = mode;
      btnNaive.classList.toggle("active", mode === 0);
      btnSplit.classList.toggle("active", mode === 1);
      btnCorrect.classList.toggle("active", mode === 2);
      requestRender();
    }

    btnNaive.addEventListener("click", () => setViewMode(0));
    btnSplit.addEventListener("click", () => setViewMode(1));
    btnCorrect.addEventListener("click", () => setViewMode(2));

    if ("IntersectionObserver" in window) {
      const io = new IntersectionObserver(
        (entries) => {
          for (const entry of entries) {
            isVisible = entry.isIntersecting;
            if (isVisible) requestRender();
          }
        },
        { threshold: 0.01 }
      );
      io.observe(container);
    } else {
      isVisible = true;
      requestRender();
    }

    if ("ResizeObserver" in window) {
      const ro = new ResizeObserver(() => {
        requestRender();
      });
      ro.observe(stage);
    } else {
      window.addEventListener("resize", requestRender, { passive: true });
    }
  }

  function initHemisphereBiasDiagram(container) {
    const stage = document.createElement("div");
    stage.className = "diagram-stage";

    const canvas = document.createElement("canvas");
    canvas.className = "diagram-canvas";
    stage.appendChild(canvas);

    const axesGizmo = createAxesSVG();
    stage.appendChild(axesGizmo.svg);

    const controls = document.createElement("div");
    controls.className = "diagram-controls";

    const toggleGroup = document.createElement("div");
    toggleGroup.className = "segmented-toggle";
    toggleGroup.setAttribute("role", "group");
    toggleGroup.setAttribute("aria-label", "Hemispherical SH projection comparison");

    const btnNaive = document.createElement("button");
    btnNaive.type = "button";
    btnNaive.className = "toggle-btn active";
    btnNaive.textContent = "Naive";

    const btnCorrect = document.createElement("button");
    btnCorrect.type = "button";
    btnCorrect.className = "toggle-btn";
    btnCorrect.textContent = "Correct";

    toggleGroup.appendChild(btnNaive);
    toggleGroup.appendChild(btnCorrect);
    controls.appendChild(toggleGroup);

    container.innerHTML = "";
    container.appendChild(stage);
    container.appendChild(controls);

    const gl = canvas.getContext("webgl2", {
      alpha: true,
      premultipliedAlpha: true,
      antialias: true,
      powerPreference: "low-power"
    });

    if (!gl) {
      stage.innerHTML =
        '<p class="webgl-fallback">WebGL2 is required to view this interactive diagram.</p>';
      return;
    }

    const vsSource = `#version 300 es
    out vec2 vUV;
    void main() {
      vec2 pos = vec2((gl_VertexID << 1) & 2, gl_VertexID & 2);
      vUV = pos;
      gl_Position = vec4(pos * 2.0 - 1.0, 0.0, 1.0);
    }`;

    const fsSource = `#version 300 es
    precision highp float;
    in vec2 vUV;
    out vec4 fragColor;

    uniform vec2 uResolution;
    uniform float uYaw;
    uniform float uPitch;
    uniform float uDist;
    uniform int uViewMode; // 0 = Naive, 1 = Correct

    // Exact 2D signed distance to an arrow along +x from x=0 to x=shaftLen+headLen, centered on y=0.
    float sdArrow2D(vec2 p, float shaftLen, float headLen, float shaftHalfW, float headHalfW) {
      float v = abs(p.y);
      float u = p.x;

      vec2 dShaft = abs(vec2(v, u - 0.5 * shaftLen)) - vec2(shaftHalfW, 0.5 * shaftLen);
      float sdShaft = length(max(dShaft, 0.0)) + min(max(dShaft.x, dShaft.y), 0.0);

      float uh = u - shaftLen;
      vec2 e = vec2(-headHalfW, headLen);
      vec2 w = vec2(v - headHalfW, uh);
      vec2 dSlant = w - e * clamp(dot(w, e) / max(dot(e, e), 1e-6), 0.0, 1.0);
      vec2 dBase = vec2(v - clamp(v, 0.0, headHalfW), uh);
      float s = (uh >= 0.0 && (v * headLen + uh * headHalfW <= headHalfW * headLen)) ? -1.0 : 1.0;
      float sdHead = sqrt(min(dot(dSlant, dSlant), dot(dBase, dBase))) * s;

      return min(sdShaft, sdHead);
    }

    // Evaluate L1 SH reconstruction for outgoing Lambertian radiance under straight-down light L = (0, 1, 0).
    // For true outgoing radiance f(n) = E * n_y on the upper hemisphere (n_y >= 0):
    // - Correct least-squares over upper hemisphere: c0 = 0,     cy = E   -> f_hat(n) = E * n_y
    // - Naive projection averaged over hemisphere:   c0 = 0.5*E, cy = E   -> f_hat(n) = E * (0.5 + n_y)
    // With E = 2/3, the maximum value across both views (Naive at n_y = 1) is 1.5 * E = 1.0.
    float evalSH(vec3 n) {
      const float E = 2.0 / 3.0;
      float c0 = (uViewMode == 0) ? (0.5 * E) : 0.0;
      float cy = E;
      return max(0.0, c0 + cy * n.y);
    }

    void main() {
      vec2 ndc = vUV * 2.0 - 1.0;
      float aspect = uResolution.x / uResolution.y;
      ndc.x *= aspect;

      const float tanFov = 0.40;

      float cp = cos(uPitch), sp = sin(uPitch);
      float cy = cos(uYaw),   sy = sin(uYaw);
      vec3 camForward = normalize(vec3(-sy * cp, -sp, -cy * cp));
      vec3 worldUp = vec3(0.0, 1.0, 0.0);
      vec3 camRight = normalize(cross(camForward, worldUp));
      vec3 camUp = cross(camRight, camForward);

      vec3 target = vec3(0.0, 0.18, 0.0);
      vec3 ro = target - camForward * uDist;
      vec3 rd = normalize(camForward + ndc.x * tanFov * camRight + ndc.y * tanFov * camUp);

      float pixelRadius = (2.0 * tanFov * uDist) / uResolution.y;

      // Square side length = 2 * S = 1.68; hemisphere radius R = (1/3) * (2 * S) = 0.56
      const float S = 0.84;
      const float R = 0.56;

      // 1. Flat square in the y = 0 plane: [-S, S] x {0} x [-S, S], normal = (0, 1, 0)
      float squareCoverage = 0.0;
      vec3 squareCol = vec3(0.0);
      float safeRdY = (abs(rd.y) < 1e-5) ? -1e-5 : rd.y;
      float tPlane = -ro.y / safeRdY;
      vec3 posP = ro + clamp(tPlane, 0.0, 50.0) * rd;
      vec2 dBox = abs(posP.xz) - vec2(S);
      float sdSquare = length(max(dBox, 0.0)) + min(max(dBox.x, dBox.y), 0.0);
      float fwSquare = fwidth(sdSquare);
      if (rd.y < -1e-4 && tPlane > 0.0) {
        float aaSquare = max(fwSquare * 0.75, pixelRadius * 0.75);
        squareCoverage = 1.0 - smoothstep(-aaSquare, aaSquare, sdSquare);

        float valP = evalSH(vec3(0.0, 1.0, 0.0));
        float ldrP = pow(clamp(valP, 0.0, 1.0), 1.0 / 2.2);
        squareCol = vec3(ldrP);
      }

      // 2. Upper hemisphere of radius R sitting at the center of the square (y >= 0)
      float bH = dot(ro, rd);
      float cH = dot(ro, ro) - R * R;
      float hH = bH * bH - cH;

      float distFromCenter = sqrt(max(0.0, R * R - hH));
      float domeCoverage = 1.0 - smoothstep(R - pixelRadius * 0.75, R + pixelRadius * 0.75, distFromCenter);

      float clampedH = max(0.0, hH);
      float tH = -bH - sqrt(clampedH);
      vec3 posH = ro + tH * rd;

      float hemiCut = smoothstep(-pixelRadius * 0.85, pixelRadius * 0.85, posH.y);
      float hemiCoverage = domeCoverage * hemiCut;

      vec3 hemiCol = vec3(0.0);
      if (hemiCoverage > 0.0) {
        vec3 nH = normalize(posH);
        float valH = evalSH(nH);
        float ldrH = pow(clamp(valH, 0.0, 1.0), 1.0 / 2.2);
        hemiCol = vec3(ldrH);
      }

      // Composite hemisphere over the flat square
      vec4 sqRgba = vec4(squareCol * squareCoverage, squareCoverage);
      vec4 hemiRgba = vec4(hemiCol * hemiCoverage, hemiCoverage);
      vec4 sceneRgba = hemiRgba + sqRgba * (1.0 - hemiRgba.a);

      // 3. Downward directional light arrow floating above the hemisphere
      vec3 pStart = vec3(0.0, 1.04, 0.0);
      vec3 pEnd   = vec3(0.0, 0.66, 0.0);

      vec3 v0 = pStart - ro;
      vec3 v1 = pEnd - ro;
      float z0 = dot(v0, camForward);
      float z1 = dot(v1, camForward);

      vec2 s0 = vec2(dot(v0, camRight), dot(v0, camUp)) / (z0 * tanFov);
      vec2 s1 = vec2(dot(v1, camRight), dot(v1, camUp)) / (z1 * tanFov);

      vec2 seg = s1 - s0;
      float L2D = length(seg);
      vec3 arrowCol = vec3(0.0);
      float arrowAlpha = 0.0;

      if (L2D > 1e-4) {
        vec2 dir2D = seg / L2D;
        vec2 rel = ndc - s0;
        vec2 pLocal = vec2(
          dot(rel, dir2D),
          rel.x * dir2D.y - rel.y * dir2D.x
        );

        float foreshorten = clamp(L2D / (0.38 * 0.75), 0.25, 1.0);
        float shaftHalfW = 0.0075;
        float headHalfW = 0.023;
        float headLen = min(0.068 * foreshorten, 0.65 * L2D);
        float shaftLen = max(0.0, L2D - headLen);

        float sd = sdArrow2D(pLocal, shaftLen, headLen, shaftHalfW, headHalfW);
        float aaNdc = (2.0 / uResolution.y) * 0.75;
        float arrowCoverage = 1.0 - smoothstep(-aaNdc, aaNdc, sd);

        if (arrowCoverage > 0.0) {
          arrowCol = vec3(0.92);
          arrowAlpha = arrowCoverage;
        }
      }

      vec4 arrowRgba = vec4(arrowCol * arrowAlpha, arrowAlpha);
      fragColor = arrowRgba + sceneRgba * (1.0 - arrowRgba.a);
    }`;

    function compileShader(type, src) {
      const s = gl.createShader(type);
      gl.shaderSource(s, src);
      gl.compileShader(s);
      if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) {
        console.error(gl.getShaderInfoLog(s));
      }
      return s;
    }

    const prog = gl.createProgram();
    gl.attachShader(prog, compileShader(gl.VERTEX_SHADER, vsSource));
    gl.attachShader(prog, compileShader(gl.FRAGMENT_SHADER, fsSource));
    gl.linkProgram(prog);

    const vao = gl.createVertexArray();
    const uLocs = {
      resolution: gl.getUniformLocation(prog, "uResolution"),
      yaw: gl.getUniformLocation(prog, "uYaw"),
      pitch: gl.getUniformLocation(prog, "uPitch"),
      dist: gl.getUniformLocation(prog, "uDist"),
      viewMode: gl.getUniformLocation(prog, "uViewMode")
    };

    // Look down on the hemisphere and square by default
    const DEFAULT_YAW = 0.38;
    const DEFAULT_PITCH = 0.48;
    let yaw = DEFAULT_YAW;
    let pitch = DEFAULT_PITCH;
    const dist = 2.85;
    let viewMode = 0; // 0 = Naive, 1 = Correct
    let isVisible = false;
    let rafId = null;

    axesGizmo.update(yaw, pitch);

    function renderFrame() {
      rafId = null;
      if (!isVisible) return;

      const dpr = Math.min(window.devicePixelRatio || 1, 2);
      const rect = stage.getBoundingClientRect();
      const w = Math.max(1, Math.round(rect.width * dpr));
      const h = Math.max(1, Math.round(rect.height * dpr));

      if (canvas.width !== w || canvas.height !== h) {
        canvas.width = w;
        canvas.height = h;
      }

      gl.viewport(0, 0, w, h);
      gl.clearColor(0.0, 0.0, 0.0, 0.0);
      gl.clear(gl.COLOR_BUFFER_BIT);

      gl.useProgram(prog);
      gl.bindVertexArray(vao);

      gl.uniform2f(uLocs.resolution, w, h);
      gl.uniform1f(uLocs.yaw, yaw);
      gl.uniform1f(uLocs.pitch, pitch);
      gl.uniform1f(uLocs.dist, dist);
      gl.uniform1i(uLocs.viewMode, viewMode);

      gl.drawArrays(gl.TRIANGLES, 0, 3);
    }

    function requestRender() {
      if (rafId === null && isVisible) {
        rafId = requestAnimationFrame(renderFrame);
      }
    }

    let dragging = false;
    let lastX = 0;
    let lastY = 0;

    canvas.addEventListener("pointerdown", (e) => {
      if (e.button !== 0 && e.pointerType === "mouse") return;
      dragging = true;
      lastX = e.clientX;
      lastY = e.clientY;
      canvas.setPointerCapture(e.pointerId);
    });

    canvas.addEventListener("pointermove", (e) => {
      if (!dragging) return;
      const dx = e.clientX - lastX;
      const dy = e.clientY - lastY;
      lastX = e.clientX;
      lastY = e.clientY;

      yaw -= dx * 0.007;
      // Keep camera above the square plane so we always view the upper surface
      pitch = Math.max(0.08, Math.min(1.45, pitch + dy * 0.007));

      axesGizmo.update(yaw, pitch);
      requestRender();
    });

    function endDrag(e) {
      if (!dragging) return;
      dragging = false;
      if (e && e.pointerId !== undefined && canvas.hasPointerCapture(e.pointerId)) {
        canvas.releasePointerCapture(e.pointerId);
      }
    }

    canvas.addEventListener("pointerup", endDrag);
    canvas.addEventListener("pointercancel", endDrag);

    canvas.addEventListener("dblclick", () => {
      yaw = DEFAULT_YAW;
      pitch = DEFAULT_PITCH;
      axesGizmo.update(yaw, pitch);
      requestRender();
    });

    function setViewMode(mode) {
      if (viewMode === mode) return;
      viewMode = mode;
      btnNaive.classList.toggle("active", mode === 0);
      btnCorrect.classList.toggle("active", mode === 1);
      requestRender();
    }

    btnNaive.addEventListener("click", () => setViewMode(0));
    btnCorrect.addEventListener("click", () => setViewMode(1));

    if ("IntersectionObserver" in window) {
      const io = new IntersectionObserver(
        (entries) => {
          for (const entry of entries) {
            isVisible = entry.isIntersecting;
            if (isVisible) requestRender();
          }
        },
        { threshold: 0.01 }
      );
      io.observe(container);
    } else {
      isVisible = true;
      requestRender();
    }

    if ("ResizeObserver" in window) {
      const ro = new ResizeObserver(() => {
        requestRender();
      });
      ro.observe(stage);
    } else {
      window.addEventListener("resize", requestRender, { passive: true });
    }
  }

  function initClamped2DDiagram(container) {
    const svgNS = "http://www.w3.org/2000/svg";
    const svg = document.createElementNS(svgNS, "svg");
    svg.setAttribute("class", "diagram-svg");
    svg.setAttribute("viewBox", "0 -4 420 266");
    svg.setAttribute("role", "img");
    svg.setAttribute(
      "aria-label",
      "2D plot comparing f(x) = x^4 with a Naive Fit, Naive Clamped, and Optimal Clamped"
    );

    // Coordinate mapping:
    // x in [0, 1] -> [x0, x0 + plotW] = [26, 370] (width = 344)
    // y in [0, 1] -> [y0, y0 - plotH] = [206, 26] (height = 180)
    const x0 = 26;
    const y0 = 206;
    const plotW = 344;
    const plotH = 180;

    function sx(x) {
      return (x0 + x * plotW).toFixed(2);
    }
    function sy(y) {
      return (y0 - y * plotH).toFixed(2);
    }

    // Axes lines (purely positive domain with arrows for x and y)
    const axisColor = "#68686e";
    const labelColor = "#b8b8be";

    // Subtle negative-y dashed extension to give visual context when Naive Fit dips below 0
    const negYAxis = document.createElementNS(svgNS, "line");
    negYAxis.setAttribute("x1", sx(0));
    negYAxis.setAttribute("y1", sy(0));
    negYAxis.setAttribute("x2", sx(0));
    negYAxis.setAttribute("y2", sy(-0.24));
    negYAxis.setAttribute("stroke", "#424246");
    negYAxis.setAttribute("stroke-width", "1.2");
    negYAxis.setAttribute("stroke-dasharray", "3 3");
    svg.appendChild(negYAxis);

    // X axis: (0,0) to xEnd, with arrow pointing directly at the "x" label
    const xEnd = 386;
    const xAxis = document.createElementNS(svgNS, "line");
    xAxis.setAttribute("x1", sx(0));
    xAxis.setAttribute("y1", sy(0));
    xAxis.setAttribute("x2", String(xEnd));
    xAxis.setAttribute("y2", sy(0));
    xAxis.setAttribute("stroke", axisColor);
    xAxis.setAttribute("stroke-width", "1.4");
    xAxis.setAttribute("stroke-linecap", "round");
    svg.appendChild(xAxis);

    const xArrow = document.createElementNS(svgNS, "polygon");
    xArrow.setAttribute(
      "points",
      `${xEnd + 5},${y0} ${xEnd - 2},${y0 - 3.2} ${xEnd - 2},${y0 + 3.2}`
    );
    xArrow.setAttribute("fill", axisColor);
    svg.appendChild(xArrow);

    const xLabel = document.createElementNS(svgNS, "text");
    xLabel.textContent = "x";
    xLabel.setAttribute("x", String(xEnd + 14));
    xLabel.setAttribute("y", String(y0 - 1.5));
    xLabel.setAttribute("fill", labelColor);
    xLabel.setAttribute("font-family", "Georgia, serif");
    xLabel.setAttribute("font-style", "italic");
    xLabel.setAttribute("font-size", "14");
    xLabel.setAttribute("text-anchor", "middle");
    xLabel.setAttribute("dominant-baseline", "central");
    svg.appendChild(xLabel);

    // Y axis: (0,0) to yEnd, with arrow pointing directly at the "y" label above it
    const yEnd = 22;
    const yAxis = document.createElementNS(svgNS, "line");
    yAxis.setAttribute("x1", sx(0));
    yAxis.setAttribute("y1", sy(0));
    yAxis.setAttribute("x2", sx(0));
    yAxis.setAttribute("y2", String(yEnd));
    yAxis.setAttribute("stroke", axisColor);
    yAxis.setAttribute("stroke-width", "1.4");
    yAxis.setAttribute("stroke-linecap", "round");
    svg.appendChild(yAxis);

    const yArrow = document.createElementNS(svgNS, "polygon");
    yArrow.setAttribute(
      "points",
      `${x0},${yEnd - 5} ${x0 - 3.2},${yEnd + 2} ${x0 + 3.2},${yEnd + 2}`
    );
    yArrow.setAttribute("fill", axisColor);
    svg.appendChild(yArrow);

    const yLabel = document.createElementNS(svgNS, "text");
    yLabel.textContent = "y";
    yLabel.setAttribute("x", String(x0));
    yLabel.setAttribute("y", String(yEnd - 17.5));
    yLabel.setAttribute("fill", labelColor);
    yLabel.setAttribute("font-family", "Georgia, serif");
    yLabel.setAttribute("font-style", "italic");
    yLabel.setAttribute("font-size", "14");
    yLabel.setAttribute("text-anchor", "middle");
    yLabel.setAttribute("dominant-baseline", "central");
    svg.appendChild(yLabel);

    // Original function: f(x) = x^4 on [0, 1]
    // Analytical least-squares fits on [0, 1]:
    // 1. Naive Fit: y = a*x + b
    //    a = 12 * (1/6 - 1/10) = 4/5 = 0.8
    //    b = 1/5 - 0.8/2 = -1/5 = -0.2
    //    Zero-crossing at x0 = 0.25, y(0) = -0.2, y(1) = 0.6
    //
    // 2. Naive Clamped: y = max(0.8*x - 0.2, 0)
    //    Flat 0 on [0, 0.25], then rises linearly from (0.25, 0) to (1, 0.6)
    //
    // 3. Optimal Clamped: min_{a,b} \int_0^1 (max(a*x + b, 0) - x^4)^2 dx
    //    Exact algebraic root 9*x_k^4 + 12*x_k^3 + 6*x_k^2 + 4*x_k - 4 = 0:
    //    x_k = 0.4905977, a = 1.6730801, b = -0.8208093 (y(1) = 0.8522708)

    const aLin = 0.8;
    const bLin = -0.2;
    const xKinkLin = 0.25;

    const xKinkOpt = 0.4905977;
    const aOpt = 1.6730801;
    const bOpt = -0.8208093;

    // Shaded negative lobe between Naive Fit and y=0 on [0, 0.25] (shown when Naive Fit is active)
    const negFill = document.createElementNS(svgNS, "polygon");
    negFill.setAttribute(
      "points",
      `${sx(0)},${sy(0)} ${sx(0)},${sy(bLin)} ${sx(xKinkLin)},${sy(0)}`
    );
    negFill.setAttribute("fill", "rgba(232, 68, 53, 0.14)");
    svg.appendChild(negFill);

    // Build path for f(x) = x^4
    const N = 100;
    let dOrig = "";
    for (let i = 0; i <= N; i++) {
      const x = i / N;
      const y = x * x * x * x;
      dOrig += (i === 0 ? "M" : "L") + `${sx(x)},${sy(y)} `;
    }
    const pathOrig = document.createElementNS(svgNS, "path");
    pathOrig.setAttribute("d", dOrig.trim());
    pathOrig.setAttribute("fill", "none");
    pathOrig.setAttribute("stroke", "#ececec");
    pathOrig.setAttribute("stroke-width", "2.4");
    pathOrig.setAttribute("stroke-linecap", "round");
    pathOrig.setAttribute("stroke-linejoin", "round");
    svg.appendChild(pathOrig);

    // Label "f(x)" right next to the endpoint of the plotted curve
    const endX = Number(sx(1));
    const endY = Number(sy(1));
    const fLabel = document.createElementNS(svgNS, "text");
    fLabel.textContent = "f(x)";
    fLabel.setAttribute("x", String(endX - 8));
    fLabel.setAttribute("y", String(endY + 3));
    fLabel.setAttribute("fill", "#ececec");
    fLabel.setAttribute("font-family", "Georgia, serif");
    fLabel.setAttribute("font-style", "italic");
    fLabel.setAttribute("font-size", "13.5");
    fLabel.setAttribute("text-anchor", "end");
    svg.appendChild(fLabel);

    // 1. Naive Fit path: y = 0.8*x - 0.2
    const pathLin = document.createElementNS(svgNS, "line");
    pathLin.setAttribute("x1", sx(0));
    pathLin.setAttribute("y1", sy(bLin));
    pathLin.setAttribute("x2", sx(1));
    pathLin.setAttribute("y2", sy(aLin + bLin));
    pathLin.setAttribute("stroke", "#569cd6");
    pathLin.setAttribute("stroke-opacity", "0.65");
    pathLin.setAttribute("stroke-width", "2.2");
    pathLin.setAttribute("stroke-dasharray", "5 4");
    pathLin.setAttribute("stroke-linecap", "round");
    svg.appendChild(pathLin);

    // 2. Naive Clamped path: y = max(0.8*x - 0.2, 0)
    const pathClampedNaive = document.createElementNS(svgNS, "path");
    pathClampedNaive.setAttribute(
      "d",
      `M${sx(0)},${sy(0)} L${sx(xKinkLin)},${sy(0)} L${sx(1)},${sy(aLin + bLin)}`
    );
    pathClampedNaive.setAttribute("fill", "none");
    pathClampedNaive.setAttribute("stroke", "#e5a646");
    pathClampedNaive.setAttribute("stroke-opacity", "0.65");
    pathClampedNaive.setAttribute("stroke-width", "2.3");
    pathClampedNaive.setAttribute("stroke-linecap", "round");
    pathClampedNaive.setAttribute("stroke-linejoin", "round");
    svg.appendChild(pathClampedNaive);

    // 3. Optimal Clamped path: y = max(1.6731*x - 0.8208, 0)
    const pathClampedOpt = document.createElementNS(svgNS, "path");
    pathClampedOpt.setAttribute(
      "d",
      `M${sx(0)},${sy(0)} L${sx(xKinkOpt)},${sy(0)} L${sx(1)},${sy(aOpt + bOpt)}`
    );
    pathClampedOpt.setAttribute("fill", "none");
    pathClampedOpt.setAttribute("stroke", "#4ec9b0");
    pathClampedOpt.setAttribute("stroke-opacity", "0.65");
    pathClampedOpt.setAttribute("stroke-width", "2.4");
    pathClampedOpt.setAttribute("stroke-linecap", "round");
    pathClampedOpt.setAttribute("stroke-linejoin", "round");
    svg.appendChild(pathClampedOpt);

    // Controls: 3 independent toggle buttons
    const controls = document.createElement("div");
    controls.className = "diagram-controls";

    const toggleGroup = document.createElement("div");
    toggleGroup.className = "segmented-toggle";
    toggleGroup.setAttribute("role", "group");
    toggleGroup.setAttribute("aria-label", "Toggle approximation curves");

    const series = [
      {
        label: "Naive Fit",
        color: "#569cd6",
        active: true,
        els: [pathLin, negFill]
      },
      {
        label: "Naive Clamped",
        color: "#e5a646",
        active: true,
        els: [pathClampedNaive]
      },
      {
        label: "Optimal Clamped",
        color: "#4ec9b0",
        active: true,
        els: [pathClampedOpt]
      }
    ];

    function syncVisibility() {
      series.forEach((s) => {
        s.els.forEach((el) => {
          el.style.display = s.active ? "" : "none";
        });
        s.btn.classList.toggle("active", s.active);
        s.btn.setAttribute("aria-pressed", String(s.active));
      });
    }

    series.forEach((s) => {
      const btn = document.createElement("button");
      btn.type = "button";
      btn.className = "toggle-btn" + (s.active ? " active" : "");
      btn.setAttribute("aria-pressed", String(s.active));

      const dot = document.createElement("span");
      dot.className = "toggle-dot";
      dot.style.backgroundColor = s.color;

      const txt = document.createTextNode(s.label);
      btn.appendChild(dot);
      btn.appendChild(txt);

      btn.addEventListener("click", () => {
        s.active = !s.active;
        syncVisibility();
      });

      s.btn = btn;
      toggleGroup.appendChild(btn);
    });

    controls.appendChild(toggleGroup);

    container.innerHTML = "";
    container.appendChild(svg);
    container.appendChild(controls);

    syncVisibility();
  }

  function initAllDiagrams() {
    const el1 = document.getElementById("sh-diagram-reflective-sphere");
    if (el1) initReflectiveSphereDiagram(el1);

    const el2 = document.getElementById("sh-diagram-linear-approx");
    if (el2) initLinearApproxDiagram(el2);

    const el3 = document.getElementById("sh-diagram-pantheon-fit");
    if (el3) initPantheonFitDiagram(el3);

    const el4 = document.getElementById("sh-diagram-clamped-2d");
    if (el4) initClamped2DDiagram(el4);

    const el5 = document.getElementById("sh-diagram-clamped-directional");
    if (el5) initClampedDirectionalDiagram(el5);

    const el6 = document.getElementById("sh-diagram-hemisphere-bias");
    if (el6) initHemisphereBiasDiagram(el6);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initAllDiagrams);
  } else {
    initAllDiagrams();
  }
})();

