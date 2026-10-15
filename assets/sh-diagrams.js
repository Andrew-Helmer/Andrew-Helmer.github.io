(function () {
  "use strict";

  const HDR_WEBP_URL = "/assets/pantheon-hdr.webp";
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
    svg.setAttribute("class", "sh-axes-gizmo");
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
      text.setAttribute("class", "sh-axis-label");

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
    stage.className = "sh-diagram-stage";

    const canvas = document.createElement("canvas");
    canvas.className = "sh-diagram-canvas";
    stage.appendChild(canvas);

    const axesGizmo = createAxesSVG();
    stage.appendChild(axesGizmo.svg);

    const controls = document.createElement("div");
    controls.className = "sh-diagram-controls";

    const toggleGroup = document.createElement("div");
    toggleGroup.className = "sh-segmented-toggle";
    toggleGroup.setAttribute("role", "group");
    toggleGroup.setAttribute("aria-label", "Sphere visualization mode");

    const btnChrome = document.createElement("button");
    btnChrome.type = "button";
    btnChrome.className = "sh-toggle-btn active";
    btnChrome.textContent = "Chrome Ball";

    const btnProj = document.createElement("button");
    btnProj.type = "button";
    btnProj.className = "sh-toggle-btn";
    btnProj.textContent = "Spherical Projection";

    toggleGroup.appendChild(btnChrome);
    toggleGroup.appendChild(btnProj);
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
        '<p class="sh-webgl-fallback">WebGL2 is required to view this interactive diagram.</p>';
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
    let projMode = 0; // 0 = Chrome Ball, 1 = Spherical Projection
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
    stage.className = "sh-diagram-stage";

    const canvas = document.createElement("canvas");
    canvas.className = "sh-diagram-canvas";
    stage.appendChild(canvas);

    const axesGizmo = createAxesSVG();
    stage.appendChild(axesGizmo.svg);

    const slidersWrap = document.createElement("div");
    slidersWrap.className = "sh-sliders";

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
      row.className = "sh-slider-row";

      const label = document.createElement("span");
      label.className = "sh-slider-label";
      label.innerHTML = `<i>c</i><sub>${def.sub}</sub>`;

      const trackWrap = document.createElement("div");
      trackWrap.className = "sh-slider-track-wrap";

      const track = document.createElement("div");
      track.className = "sh-slider-track";

      const fill = document.createElement("div");
      fill.className = "sh-slider-fill";
      track.appendChild(fill);

      const input = document.createElement("input");
      input.type = "range";
      input.className = "sh-slider-input";
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
        '<p class="sh-webgl-fallback">WebGL2 is required to view this interactive diagram.</p>';
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
    stage.className = "sh-diagram-stage";

    const canvas = document.createElement("canvas");
    canvas.className = "sh-diagram-canvas";
    stage.appendChild(canvas);

    const axesGizmo = createAxesSVG();
    stage.appendChild(axesGizmo.svg);

    const controls = document.createElement("div");
    controls.className = "sh-diagram-controls";

    const toggleGroup = document.createElement("div");
    toggleGroup.className = "sh-segmented-toggle";
    toggleGroup.setAttribute("role", "group");
    toggleGroup.setAttribute("aria-label", "Spherical harmonics comparison view");

    const btnProj = document.createElement("button");
    btnProj.type = "button";
    btnProj.className = "sh-toggle-btn active";
    btnProj.textContent = "Original";

    const btnSH = document.createElement("button");
    btnSH.type = "button";
    btnSH.className = "sh-toggle-btn";
    btnSH.textContent = "Approximation";

    const btnErr = document.createElement("button");
    btnErr.type = "button";
    btnErr.className = "sh-toggle-btn";
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
        '<p class="sh-webgl-fallback">WebGL2 is required to view this interactive diagram.</p>';
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

  function initAllDiagrams() {
    const el1 = document.getElementById("sh-diagram-reflective-sphere");
    if (el1) initReflectiveSphereDiagram(el1);

    const el2 = document.getElementById("sh-diagram-linear-approx");
    if (el2) initLinearApproxDiagram(el2);

    const el3 = document.getElementById("sh-diagram-pantheon-fit");
    if (el3) initPantheonFitDiagram(el3);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initAllDiagrams);
  } else {
    initAllDiagrams();
  }
})();

