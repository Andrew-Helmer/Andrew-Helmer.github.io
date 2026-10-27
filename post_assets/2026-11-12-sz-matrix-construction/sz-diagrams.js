(function () {
  "use strict";

  // First 4 dimensions of `hierarchical_b2_xors` (for up to 256 points, m = 0..7)
  // from "Stochastic Generation of (t,s) Sample Sequences" (Helmer, Christensen, Kensler 2021).
  // Dim 0 (x) and Dim 1 (y) form a base-2 (0,2)-sequence;
  // Dim 2 (z) and Dim 3 (w) form a base-2 (0,2)-sequence;
  // Together (x, y, z, w) form the canonical 4D SZ sequence.
  const HIERARCHICAL_B2_XORS_4D = [
    [0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
    [0x0, 0x1, 0x1, 0x7, 0x1, 0x13, 0x15, 0x7f],
    [0x0, 0x0, 0x2, 0x3, 0x3, 0x1, 0x2d, 0x36],
    [0x0, 0x1, 0x3, 0x6, 0x2, 0x11, 0x39, 0x67]
  ];

  const SOBOL_02_Y_XORS = HIERARCHICAL_B2_XORS_4D[1];

  function randomUint32() {
    if (typeof crypto !== "undefined" && crypto.getRandomValues) {
      const buf = new Uint32Array(1);
      crypto.getRandomValues(buf);
      return buf[0] >>> 0;
    }
    return (Math.random() * 0x100000000) >>> 0;
  }

  // Generate `n` points (up to 256) of a stochastic base-2 Sobol (0,2)-sequence.
  function generateStochasticSobol02(n) {
    const xInts = new Uint32Array(n);
    const yInts = new Uint32Array(n);

    xInts[0] = randomUint32();
    yInts[0] = randomUint32();

    let m = 0;
    for (let prevLen = 1; prevLen < n; prevLen <<= 1, m++) {
      const bitToFlip = (1 << (31 - m)) >>> 0;
      const randMask = (bitToFlip - 1) >>> 0;
      const yXor = SOBOL_02_Y_XORS[m];

      for (let k = 0; k < prevLen && prevLen + k < n; k++) {
        const rx = randomUint32() & randMask;
        const ry = randomUint32() & randMask;

        xInts[prevLen + k] = (xInts[k] ^ bitToFlip ^ rx) >>> 0;
        yInts[prevLen + k] = (yInts[k ^ yXor] ^ bitToFlip ^ ry) >>> 0;
      }
    }

    const inv2p32 = 1.0 / 4294967296.0;
    const pts = new Array(n);
    for (let i = 0; i < n; i++) {
      pts[i] = {
        x: xInts[i] * inv2p32,
        y: yInts[i] * inv2p32
      };
    }
    return pts;
  }

  // Generate `n` 4D points (up to 256) using the first 4 dimensions of `hierarchical_b2_xors`.
  function generateStochasticSZ4D(n) {
    const ints = [
      new Uint32Array(n),
      new Uint32Array(n),
      new Uint32Array(n),
      new Uint32Array(n)
    ];

    for (let d = 0; d < 4; d++) {
      ints[d][0] = randomUint32();
    }

    let m = 0;
    for (let prevLen = 1; prevLen < n; prevLen <<= 1, m++) {
      const bitToFlip = (1 << (31 - m)) >>> 0;
      const randMask = (bitToFlip - 1) >>> 0;

      for (let d = 0; d < 4; d++) {
        const xorVal = HIERARCHICAL_B2_XORS_4D[d][m];
        const arr = ints[d];
        for (let k = 0; k < prevLen && prevLen + k < n; k++) {
          const r = randomUint32() & randMask;
          arr[prevLen + k] = (arr[k ^ xorVal] ^ bitToFlip ^ r) >>> 0;
        }
      }
    }

    const inv2p32 = 1.0 / 4294967296.0;
    const pts = new Array(n);
    for (let i = 0; i < n; i++) {
      pts[i] = [
        ints[0][i] * inv2p32,
        ints[1][i] * inv2p32,
        ints[2][i] * inv2p32,
        ints[3][i] * inv2p32
      ];
    }
    return pts;
  }

  // Shared helper: calculate effective log_b(x_strata) from total exponent `k` and `offset`
  // without mutating the stored `offset`.
  function computeLogXStrata(k, offset) {
    if (k <= 0) return 0;
    let effectiveOffset = offset;
    const isOddPower = (k & 1) === 1;
    const isWholeOffset = Math.abs(offset - Math.round(offset)) < 1e-6;

    if (isOddPower && isWholeOffset) {
      // Whole-number offset with an odd power rounds down to the lower 0.5 effective value
      // (e.g. offset 0 -> -0.5, giving 4x2 rather than 2x4 for 8 points).
      effectiveOffset = offset - 0.5;
    } else if (!isOddPower && !isWholeOffset) {
      // x.5 offset with an even power rounds towards 0
      // (e.g. -0.5 or +0.5 for 8 points both round to 0 -> 4x4 when moving to 16 points).
      effectiveOffset = Math.trunc(offset);
    }

    const rawLogX = k / 2.0 - effectiveOffset;
    return Math.max(0, Math.min(k, Math.round(rawLogX)));
  }

  function createDiscreteSliderRow(labelText, minVal, maxVal, initialVal, ariaLabel) {
    const row = document.createElement("div");
    row.className = "slider-row";

    const label = document.createElement("span");
    label.className = "slider-label slider-label-wide";
    label.textContent = labelText;

    const trackWrap = document.createElement("div");
    trackWrap.className = "slider-track-wrap";

    const track = document.createElement("div");
    track.className = "slider-track";

    const fill = document.createElement("div");
    fill.className = "slider-fill";
    track.appendChild(fill);

    const ticksContainer = document.createElement("div");
    ticksContainer.className = "slider-ticks";

    const input = document.createElement("input");
    input.type = "range";
    input.className = "slider-input";
    input.min = String(minVal);
    input.max = String(maxVal);
    input.step = "1";
    input.value = String(initialVal);
    input.setAttribute("aria-label", ariaLabel);

    const valLabel = document.createElement("span");
    valLabel.className = "slider-value";

    trackWrap.appendChild(track);
    trackWrap.appendChild(ticksContainer);
    trackWrap.appendChild(input);

    row.appendChild(label);
    row.appendChild(trackWrap);
    row.appendChild(valLabel);

    function renderTicks(count) {
      ticksContainer.innerHTML = "";
      if (count <= 1) return;
      for (let i = 0; i < count; i++) {
        const tick = document.createElement("span");
        tick.className = "slider-tick";
        const pct = (i / (count - 1)) * 100;
        tick.style.left = `${pct}%`;
        ticksContainer.appendChild(tick);
      }
    }

    return { row, input, fill, ticksContainer, valLabel, renderTicks };
  }

  // ============================================================================
  // DIAGRAM 1: 2D Base-2 (0,2)-Sequence Stratification
  // ============================================================================
  function initStratification02Diagram(container) {
    const MAX_POINTS = 256;
    let points = generateStochasticSobol02(MAX_POINTS);

    const svgNS = "http://www.w3.org/2000/svg";
    const svg = document.createElementNS(svgNS, "svg");
    svg.setAttribute("class", "diagram-svg");
    svg.setAttribute("viewBox", "0 0 340 318");
    svg.setAttribute("role", "img");
    svg.setAttribute(
      "aria-label",
      "Interactive 2D stratification diagram of a randomized base-2 (0,2)-sequence"
    );

    const boxX = 20;
    const boxY = 12;
    const boxSize = 300;

    const bgRect = document.createElementNS(svgNS, "rect");
    bgRect.setAttribute("x", String(boxX));
    bgRect.setAttribute("y", String(boxY));
    bgRect.setAttribute("width", String(boxSize));
    bgRect.setAttribute("height", String(boxSize));
    bgRect.setAttribute("fill", "#1e1e1e");
    svg.appendChild(bgRect);

    const gridGroup = document.createElementNS(svgNS, "g");
    svg.appendChild(gridGroup);

    const pointsGroup = document.createElementNS(svgNS, "g");
    svg.appendChild(pointsGroup);

    const borderRect = document.createElementNS(svgNS, "rect");
    borderRect.setAttribute("x", String(boxX));
    borderRect.setAttribute("y", String(boxY));
    borderRect.setAttribute("width", String(boxSize));
    borderRect.setAttribute("height", String(boxSize));
    borderRect.setAttribute("fill", "none");
    borderRect.setAttribute("stroke", "#bebebe");
    borderRect.setAttribute("stroke-width", "1.35");
    borderRect.setAttribute("stroke-linejoin", "miter");
    svg.appendChild(borderRect);

    const pointEls = new Array(MAX_POINTS);
    for (let i = 0; i < MAX_POINTS; i++) {
      const c = document.createElementNS(svgNS, "circle");
      c.setAttribute("fill", "#19d4fa");
      c.setAttribute("stroke", "#1e1e1e");
      c.setAttribute("stroke-width", "0.75");
      pointsGroup.appendChild(c);
      pointEls[i] = c;
    }

    function syncPointPositions() {
      for (let i = 0; i < MAX_POINTS; i++) {
        const cx = boxX + points[i].x * boxSize;
        const cy = boxY + (1.0 - points[i].y) * boxSize;
        pointEls[i].setAttribute("cx", cx.toFixed(2));
        pointEls[i].setAttribute("cy", cy.toFixed(2));
      }
    }

    syncPointPositions();

    let currentK = 4;
    let internalOffset = 0.0;

    const slidersWrap = document.createElement("div");
    slidersWrap.className = "diagram-sliders";
    slidersWrap.style.maxWidth = "340px";

    const pointsSlider = createDiscreteSliderRow(
      "Samples",
      1,
      8,
      currentK,
      "Power-of-two sample count"
    );
    pointsSlider.renderTicks(8);

    const strataSlider = createDiscreteSliderRow(
      "Strata",
      0,
      currentK,
      2,
      "Stratification grid aspect"
    );

    const controls = document.createElement("div");
    controls.className = "diagram-controls";
    controls.style.marginTop = "2px";
    controls.style.marginBottom = "2px";

    const regenBtn = document.createElement("button");
    regenBtn.type = "button";
    regenBtn.className = "action-btn";
    regenBtn.textContent = "Regenerate";
    regenBtn.addEventListener("click", () => {
      points = generateStochasticSobol02(MAX_POINTS);
      syncPointPositions();
    });
    controls.appendChild(regenBtn);

    slidersWrap.appendChild(pointsSlider.row);
    slidersWrap.appendChild(strataSlider.row);

    container.innerHTML = "";
    container.appendChild(svg);
    container.appendChild(controls);
    container.appendChild(slidersWrap);

    function updateDiagram() {
      const n = 1 << currentK;
      const log2X = computeLogXStrata(currentK, internalOffset);
      const log2Y = currentK - log2X;
      const xStrata = 1 << log2X;
      const yStrata = 1 << log2Y;

      const pointsPct = ((currentK - 1) / (8 - 1)) * 100;
      pointsSlider.fill.style.left = "0%";
      pointsSlider.fill.style.width = `${pointsPct}%`;
      pointsSlider.valLabel.innerHTML = `${n}`;

      const sVal = currentK - log2X;
      strataSlider.input.max = String(currentK);
      strataSlider.input.value = String(sVal);
      strataSlider.renderTicks(currentK + 1);

      const strataPct = currentK > 0 ? (sVal / currentK) * 100 : 0;
      strataSlider.fill.style.left = "0%";
      strataSlider.fill.style.width = `${strataPct}%`;
      strataSlider.valLabel.innerHTML = `${xStrata}&thinsp;&times;&thinsp;${yStrata}`;

      gridGroup.innerHTML = "";
      const strokeColor = currentK >= 7 ? "#4a4a50" : "#5a5a60";
      const strokeWidth = currentK >= 7 ? "0.75" : "1.0";

      for (let ix = 1; ix < xStrata; ix++) {
        const xPos = boxX + (ix / xStrata) * boxSize;
        const line = document.createElementNS(svgNS, "line");
        line.setAttribute("x1", xPos.toFixed(2));
        line.setAttribute("y1", String(boxY));
        line.setAttribute("x2", xPos.toFixed(2));
        line.setAttribute("y2", String(boxY + boxSize));
        line.setAttribute("stroke", strokeColor);
        line.setAttribute("stroke-width", strokeWidth);
        gridGroup.appendChild(line);
      }

      for (let iy = 1; iy < yStrata; iy++) {
        const yPos = boxY + (iy / yStrata) * boxSize;
        const line = document.createElementNS(svgNS, "line");
        line.setAttribute("x1", String(boxX));
        line.setAttribute("y1", yPos.toFixed(2));
        line.setAttribute("x2", String(boxX + boxSize));
        line.setAttribute("y2", yPos.toFixed(2));
        line.setAttribute("stroke", strokeColor);
        line.setAttribute("stroke-width", strokeWidth);
        gridGroup.appendChild(line);
      }

      const radius =
        currentK <= 4 ? "3.4" : currentK <= 6 ? "2.8" : currentK === 7 ? "2.3" : "1.9";

      for (let i = 0; i < MAX_POINTS; i++) {
        const el = pointEls[i];
        if (i < n) {
          el.style.display = "";
          el.setAttribute("r", radius);
        } else {
          el.style.display = "none";
        }
      }
    }

    pointsSlider.input.addEventListener("input", () => {
      currentK = parseInt(pointsSlider.input.value, 10);
      updateDiagram();
    });

    pointsSlider.input.addEventListener("dblclick", () => {
      currentK = 4;
      pointsSlider.input.value = "4";
      updateDiagram();
    });

    strataSlider.input.addEventListener("input", () => {
      const sVal = parseInt(strataSlider.input.value, 10);
      const log2X = currentK - sVal;
      internalOffset = currentK / 2.0 - log2X;
      updateDiagram();
    });

    strataSlider.input.addEventListener("dblclick", () => {
      internalOffset = 0.0;
      updateDiagram();
    });

    updateDiagram();
  }

  // ============================================================================
  // DIAGRAM 2: Multi-Dimensional (2D / 3D / 4D) SZ Sequence Visualization
  // ============================================================================
  function createCustomAxesGizmo(axisNames) {
    const svgNS = "http://www.w3.org/2000/svg";
    const svg = document.createElementNS(svgNS, "svg");
    svg.setAttribute("class", "cube-axes-gizmo");
    svg.setAttribute("viewBox", "-36 -36 72 72");
    svg.setAttribute("aria-hidden", "true");

    const axes = [0, 1, 2].map((idx) => {
      const g = document.createElementNS(svgNS, "g");
      const line = document.createElementNS(svgNS, "line");
      line.setAttribute("x1", "0");
      line.setAttribute("y1", "0");
      line.setAttribute("stroke-width", "1.15");
      line.setAttribute("stroke-linecap", "round");

      const arrow = document.createElementNS(svgNS, "polygon");

      const text = document.createElementNS(svgNS, "text");
      text.textContent = axisNames[idx];
      text.setAttribute("text-anchor", "middle");
      text.setAttribute("dominant-baseline", "central");
      text.setAttribute("class", "axis-label");

      g.appendChild(line);
      g.appendChild(arrow);
      g.appendChild(text);
      svg.appendChild(g);
      return { g, line, arrow, text };
    });

    const originDot = document.createElementNS(svgNS, "circle");
    originDot.setAttribute("cx", "0");
    originDot.setAttribute("cy", "0");
    originDot.setAttribute("r", "1.5");
    originDot.setAttribute("fill", "#8a8a8a");
    svg.appendChild(originDot);

    function setAxisNames(names) {
      for (let i = 0; i < 3; i++) {
        axes[i].text.textContent = names[i];
      }
    }

    function update(yaw, pitch) {
      const cp = Math.cos(pitch),
        sp = Math.sin(pitch);
      const cy = Math.cos(yaw),
        sy = Math.sin(yaw);

      const fx = -sy * cp,
        fy = -sp,
        fz = -cy * cp;
      const rx = cy,
        ry = 0,
        rz = -sy;
      const ux = -sy * sp,
        uy = cp,
        uz = -cy * sp;

      const worldAxes = [
        { item: axes[0], wx: 1, wy: 0, wz: 0 },
        { item: axes[1], wx: 0, wy: 1, wz: 0 },
        { item: axes[2], wx: 0, wy: 0, wz: 1 }
      ];

      const axisLen = 19;
      const labelDist = 27;

      const projected = worldAxes.map(({ item, wx, wy, wz }) => {
        const vx = wx * rx + wy * ry + wz * rz;
        const vy = wx * ux + wy * uy + wz * uz;
        const vz = -(wx * fx + wy * fy + wz * fz);
        return { item, vx, vy, vz };
      });

      projected.sort((a, b) => a.vz - b.vz);

      projected.forEach(({ item, vx, vy, vz }) => {
        svg.appendChild(item.g);
        const sx = vx * axisLen;
        const syScreen = -vy * axisLen;

        const depthFactor = 0.5 * (vz + 1.0);
        const grayVal = Math.round(95 + 75 * depthFactor);
        const labelGray = Math.round(140 + 80 * depthFactor);
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
          const headLen = 3.8;
          const headWidth = 1.65;
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

    return { svg, update, setAxisNames };
  }

  function initSZ4DProjectionsDiagram(container) {
    const MAX_POINTS = 256;
    let points4D = generateStochasticSZ4D(MAX_POINTS);

    let viewMode = "2D"; // "2D" | "3D" | "4D"
    let currentK = 4; // log2(N) in 1..8, defaults to 4 (16 points)

    // 2D stratification internal offset (same behavior as Diagram 1)
    let internalOffset2D = 0.0;

    // 3D stratification internal state:
    // For N = 2^k, let m4 = ceil(k / 2) be the base-4 budget so 4^m4 >= N (at most 1 point per cell).
    // The slider has T = m4 + 2 ticks (indices 0 .. m4 + 1) and two non-overlapping knobs
    // at integer tick indices 0 <= t1 < t2 <= m4 + 1.
    // Exponents are p1 = t1, p2 = t2 - t1 - 1, p3 = (m4 + 1) - t2 (with p1 + p2 + p3 = m4).
    // Internally we store two offsets (offset3D_1, offset3D_2) from the 1/3 and 2/3 cumulative splits:
    //   p1_ideal = m4 / 3 - offset3D_1
    //   (p1 + p2)_ideal = 2 * m4 / 3 - offset3D_2
    // Starting at (0, 0) gives:
    //   m4 = 1 (N=2,4):    O-X-X     -> 4x1x1
    //   m4 = 2 (N=8,16):   O-X-O-X   -> 4x4x1
    //   m4 = 3 (N=32,64):  O-X-O-X-O -> 4x4x4
    //   m4 = 4 (N=128,256): O-X-O-X-O-O -> 4x4x16
    let offset3D_1 = 0.0;
    let offset3D_2 = 0.0;

    function compute3DKnobs(m4) {
      if (m4 <= 0) {
        return { t1: 0, t2: 1, p1: 0, p2: 0, p3: 0, numTicks: 2 };
      }
      const numTicks = m4 + 2;
      const c1Raw = m4 / 3.0 - offset3D_1;
      const c2Raw = (2.0 * m4) / 3.0 - offset3D_2;

      let p1 = Math.max(0, Math.min(m4, Math.round(c1Raw)));
      let c2 = Math.max(p1, Math.min(m4, Math.round(c2Raw)));
      let p2 = c2 - p1;
      let p3 = m4 - c2;

      return {
        t1: p1,
        t2: p1 + p2 + 1,
        p1,
        p2,
        p3,
        numTicks
      };
    }

    // Top mode selector bar: 2D | 3D | 4D
    const topControls = document.createElement("div");
    topControls.className = "diagram-controls";
    topControls.style.marginBottom = "10px";

    const modeGroup = document.createElement("div");
    modeGroup.className = "segmented-toggle";
    modeGroup.setAttribute("role", "group");
    modeGroup.setAttribute("aria-label", "Visualization dimension mode");

    const modeButtons = {};
    ["2D", "3D", "4D"].forEach((m) => {
      const btn = document.createElement("button");
      btn.type = "button";
      btn.className = "toggle-btn" + (m === viewMode ? " active" : "");
      btn.textContent = m;
      btn.addEventListener("click", () => {
        if (viewMode === m) return;
        viewMode = m;
        Object.keys(modeButtons).forEach((k) => {
          modeButtons[k].classList.toggle("active", k === viewMode);
        });
        updateAll();
      });
      modeButtons[m] = btn;
      modeGroup.appendChild(btn);
    });
    topControls.appendChild(modeGroup);

    // ==========================================================================
    // 2D Stage: 2 rows of 3 2D projections (xy, xz, xw, yz, yw, zw)
    // ==========================================================================
    const svgNS = "http://www.w3.org/2000/svg";
    const svg2D = document.createElementNS(svgNS, "svg");
    svg2D.setAttribute("class", "diagram-svg");
    svg2D.setAttribute("viewBox", "0 0 560 390");
    svg2D.setAttribute("role", "img");
    svg2D.setAttribute(
      "aria-label",
      "6 pairwise 2D projections (xy, xz, xw, yz, yw, zw) of the 4D SZ sequence"
    );

    const PROJECTIONS_2D = [
      { d0: 0, d1: 1, name: "xy", isBase2: true, col: 0, row: 0 },
      { d0: 0, d1: 2, name: "xz", isBase2: false, col: 1, row: 0 },
      { d0: 0, d1: 3, name: "xw", isBase2: false, col: 2, row: 0 },
      { d0: 1, d1: 2, name: "yz", isBase2: false, col: 0, row: 1 },
      { d0: 1, d1: 3, name: "yw", isBase2: false, col: 1, row: 1 },
      { d0: 2, d1: 3, name: "zw", isBase2: true, col: 2, row: 1 }
    ];

    const subSize = 158;
    const colGap = 28;
    const rowGap = 36;
    const padX = (560 - (3 * subSize + 2 * colGap)) / 2; // 15
    const padY = 24;

    const subplots2D = PROJECTIONS_2D.map((proj) => {
      const bx = padX + proj.col * (subSize + colGap);
      const by = padY + proj.row * (subSize + rowGap);

      const g = document.createElementNS(svgNS, "g");

      // Header above subplot: e.g. "xy   4 × 4"
      const titleText = document.createElementNS(svgNS, "text");
      titleText.setAttribute("x", (bx + subSize / 2).toFixed(1));
      titleText.setAttribute("y", (by - 7).toFixed(1));
      titleText.setAttribute("text-anchor", "middle");
      titleText.setAttribute("fill", "#d9d9d9");
      titleText.setAttribute("font-family", "Georgia, serif");
      titleText.setAttribute("font-size", "13.5");

      const nameSpan = document.createElementNS(svgNS, "tspan");
      nameSpan.setAttribute("font-style", "italic");
      nameSpan.textContent = proj.name;

      const gridSpan = document.createElementNS(svgNS, "tspan");
      gridSpan.setAttribute("fill", "#8e8e96");
      gridSpan.setAttribute("font-family", '"Helvetica Neue", Helvetica, Arial, sans-serif');
      gridSpan.setAttribute("font-size", "11.5");
      gridSpan.setAttribute("dx", "6");

      titleText.appendChild(nameSpan);
      titleText.appendChild(gridSpan);
      g.appendChild(titleText);

      const bg = document.createElementNS(svgNS, "rect");
      bg.setAttribute("x", String(bx));
      bg.setAttribute("y", String(by));
      bg.setAttribute("width", String(subSize));
      bg.setAttribute("height", String(subSize));
      bg.setAttribute("fill", "#1e1e1e");
      g.appendChild(bg);

      const gridG = document.createElementNS(svgNS, "g");
      g.appendChild(gridG);

      const ptsG = document.createElementNS(svgNS, "g");
      g.appendChild(ptsG);

      const border = document.createElementNS(svgNS, "rect");
      border.setAttribute("x", String(bx));
      border.setAttribute("y", String(by));
      border.setAttribute("width", String(subSize));
      border.setAttribute("height", String(subSize));
      border.setAttribute("fill", "none");
      border.setAttribute("stroke", "#bebebe");
      border.setAttribute("stroke-width", "1.25");
      border.setAttribute("stroke-linejoin", "miter");
      g.appendChild(border);

      const circles = new Array(MAX_POINTS);
      for (let i = 0; i < MAX_POINTS; i++) {
        const c = document.createElementNS(svgNS, "circle");
        c.setAttribute("fill", "#19d4fa");
        c.setAttribute("stroke", "#1e1e1e");
        c.setAttribute("stroke-width", "0.6");
        ptsG.appendChild(c);
        circles[i] = c;
      }

      svg2D.appendChild(g);
      return { proj, bx, by, gridG, circles, gridSpan };
    });

    function sync2DPointPositions() {
      subplots2D.forEach(({ proj, bx, by, circles }) => {
        const d0 = proj.d0;
        const d1 = proj.d1;
        for (let i = 0; i < MAX_POINTS; i++) {
          const cx = bx + points4D[i][d0] * subSize;
          const cy = by + (1.0 - points4D[i][d1]) * subSize;
          circles[i].setAttribute("cx", cx.toFixed(2));
          circles[i].setAttribute("cy", cy.toFixed(2));
        }
      });
    }

    sync2DPointPositions();

    // ==========================================================================
    // 3D & 4D Stage: Orthographic WebGL2 Ray-Cast Cubes + HTML/SVG Overlay
    // ==========================================================================
    const cubeStageWrap = document.createElement("div");
    cubeStageWrap.className = "cube-stage-wrap";
    cubeStageWrap.style.display = "none";

    const glCanvas = document.createElement("canvas");
    glCanvas.className = "diagram-canvas";
    cubeStageWrap.appendChild(glCanvas);

    const cubeOverlay = document.createElement("div");
    cubeOverlay.className = "cube-grid-overlay";
    cubeStageWrap.appendChild(cubeOverlay);

    const cubeCells = [0, 1, 2, 3].map(() => {
      const cell = document.createElement("div");
      cell.className = "cube-cell-overlay";

      const header = document.createElement("div");
      header.className = "cube-cell-header";
      cell.appendChild(header);

      const gizmo = createCustomAxesGizmo(["x", "y", "z"]);
      cell.appendChild(gizmo.svg);

      cubeOverlay.appendChild(cell);
      return { cell, header, gizmo };
    });

    // Mid controls: Regenerate button
    const midControls = document.createElement("div");
    midControls.className = "diagram-controls";
    midControls.style.marginTop = "4px";
    midControls.style.marginBottom = "2px";

    const regenBtn = document.createElement("button");
    regenBtn.type = "button";
    regenBtn.className = "action-btn";
    regenBtn.textContent = "Regenerate";
    midControls.appendChild(regenBtn);

    // Bottom Sliders:
    // 1. Samples slider (shared across 2D, 3D, 4D)
    // 2. 2D Strata slider (shown in 2D mode)
    // 3. 3D Dual-Knob Strata slider (shown in 3D mode)
    const slidersWrap = document.createElement("div");
    slidersWrap.className = "diagram-sliders";
    slidersWrap.style.maxWidth = "370px";

    const samplesSlider = createDiscreteSliderRow(
      "Samples",
      1,
      8,
      currentK,
      "Power-of-two sample count"
    );
    samplesSlider.renderTicks(8);

    const strata2DSlider = createDiscreteSliderRow(
      "Strata",
      0,
      currentK,
      2,
      "2D stratification granularity"
    );

    // Custom Dual-Knob 3D Base-4 Strata Slider
    function createDualKnob3DSliderRow() {
      const row = document.createElement("div");
      row.className = "slider-row";
      row.style.display = "none";

      const label = document.createElement("span");
      label.className = "slider-label slider-label-wide";
      label.textContent = "Strata";

      const trackWrap = document.createElement("div");
      trackWrap.className = "slider-track-wrap dual-slider-interactive";
      trackWrap.setAttribute("role", "group");
      trackWrap.setAttribute("aria-label", "3D base-4 stratification dual slider");

      const track = document.createElement("div");
      track.className = "slider-track";

      const fillMid = document.createElement("div");
      fillMid.className = "slider-fill";
      track.appendChild(fillMid);

      const ticksContainer = document.createElement("div");
      ticksContainer.className = "slider-ticks";

      const thumbsContainer = document.createElement("div");
      thumbsContainer.className = "dual-slider-thumbs";

      const thumb1 = document.createElement("div");
      thumb1.className = "slider-thumb";
      const thumb2 = document.createElement("div");
      thumb2.className = "slider-thumb";
      thumbsContainer.appendChild(thumb1);
      thumbsContainer.appendChild(thumb2);

      const valLabel = document.createElement("span");
      valLabel.className = "slider-value";

      trackWrap.appendChild(track);
      trackWrap.appendChild(ticksContainer);
      trackWrap.appendChild(thumbsContainer);

      row.appendChild(label);
      row.appendChild(trackWrap);
      row.appendChild(valLabel);

      let activeKnob = null; // 1 or 2

      function getFractionalTick(clientX, numTicks) {
        const rect = trackWrap.getBoundingClientRect();
        const usableW = Math.max(1, rect.width - 14);
        const relX = clientX - (rect.left + 7);
        const t = (relX / usableW) * (numTicks - 1);
        return Math.max(0, Math.min(numTicks - 1, t));
      }

      function applyKnobTick(knobIdx, targetTick) {
        const m4 = Math.ceil(currentK / 2);
        if (m4 <= 0) return;
        const state = compute3DKnobs(m4);
        const numTicks = state.numTicks;
        let t1 = state.t1;
        let t2 = state.t2;

        const clamped = Math.max(0, Math.min(numTicks - 1, Math.round(targetTick)));
        if (knobIdx === 1) {
          if (clamped >= t2) {
            if (clamped < numTicks - 1) {
              t1 = clamped;
              t2 = clamped + 1;
            } else {
              t1 = numTicks - 2;
              t2 = numTicks - 1;
            }
          } else {
            t1 = clamped;
          }
        } else {
          if (clamped <= t1) {
            if (clamped > 0) {
              t2 = clamped;
              t1 = clamped - 1;
            } else {
              t1 = 0;
              t2 = 1;
            }
          } else {
            t2 = clamped;
          }
        }

        const p1 = t1;
        const p2 = t2 - t1 - 1;
        offset3D_1 = m4 / 3.0 - p1;
        offset3D_2 = (2.0 * m4) / 3.0 - (p1 + p2);
        updateAll();
      }

      trackWrap.addEventListener("pointerdown", (e) => {
        if (e.button !== 0 && e.pointerType === "mouse") return;
        const m4 = Math.ceil(currentK / 2);
        if (m4 <= 0) return;
        const state = compute3DKnobs(m4);
        const u = getFractionalTick(e.clientX, state.numTicks);

        const d1 = Math.abs(u - state.t1);
        const d2 = Math.abs(u - state.t2);
        if (d1 < d2) {
          activeKnob = 1;
        } else if (d2 < d1) {
          activeKnob = 2;
        } else {
          activeKnob = u < state.t1 ? 1 : 2;
        }

        thumb1.classList.toggle("active", activeKnob === 1);
        thumb2.classList.toggle("active", activeKnob === 2);
        trackWrap.setPointerCapture(e.pointerId);
        applyKnobTick(activeKnob, u);
      });

      trackWrap.addEventListener("pointermove", (e) => {
        if (activeKnob === null) return;
        const m4 = Math.ceil(currentK / 2);
        if (m4 <= 0) return;
        const state = compute3DKnobs(m4);
        const u = getFractionalTick(e.clientX, state.numTicks);
        applyKnobTick(activeKnob, u);
      });

      function endKnobDrag(e) {
        if (activeKnob === null) return;
        activeKnob = null;
        thumb1.classList.remove("active");
        thumb2.classList.remove("active");
        if (
          e &&
          e.pointerId !== undefined &&
          trackWrap.hasPointerCapture(e.pointerId)
        ) {
          trackWrap.releasePointerCapture(e.pointerId);
        }
      }

      trackWrap.addEventListener("pointerup", endKnobDrag);
      trackWrap.addEventListener("pointercancel", endKnobDrag);

      trackWrap.addEventListener("dblclick", () => {
        offset3D_1 = 0.0;
        offset3D_2 = 0.0;
        updateAll();
      });

      function renderState(m4) {
        const { t1, t2, p1, p2, p3, numTicks } = compute3DKnobs(m4);

        ticksContainer.innerHTML = "";
        for (let i = 0; i < numTicks; i++) {
          const tick = document.createElement("span");
          tick.className = "slider-tick";
          const pct = (i / (numTicks - 1)) * 100;
          tick.style.left = `${pct}%`;
          ticksContainer.appendChild(tick);
        }

        const pct1 = (t1 / (numTicks - 1)) * 100;
        const pct2 = (t2 / (numTicks - 1)) * 100;

        thumb1.style.left = `${pct1}%`;
        thumb2.style.left = `${pct2}%`;

        fillMid.style.left = `${pct1}%`;
        fillMid.style.width = `${Math.max(0, pct2 - pct1)}%`;

        const s1 = 1 << (2 * p1);
        const s2 = 1 << (2 * p2);
        const s3 = 1 << (2 * p3);
        valLabel.innerHTML = `${s1}&thinsp;&times;&thinsp;${s2}&thinsp;&times;&thinsp;${s3}`;
      }

      return { row, renderState };
    }

    const strata3DSlider = createDualKnob3DSliderRow();

    slidersWrap.appendChild(samplesSlider.row);
    slidersWrap.appendChild(strata2DSlider.row);
    slidersWrap.appendChild(strata3DSlider.row);

    const multiStageWrap = document.createElement("div");
    multiStageWrap.className = "multi-stage-wrap";
    multiStageWrap.appendChild(svg2D);
    multiStageWrap.appendChild(cubeStageWrap);

    container.innerHTML = "";
    container.appendChild(topControls);
    container.appendChild(multiStageWrap);
    container.appendChild(midControls);
    container.appendChild(slidersWrap);

    // ==========================================================================
    // WebGL2 Orthographic Cube Ray-Caster Setup
    // ==========================================================================
    const gl = glCanvas.getContext("webgl2", {
      alpha: true,
      premultipliedAlpha: true,
      antialias: true,
      powerPreference: "low-power"
    });

    let glProg = null;
    let glVao = null;
    let uLocs = null;

    if (gl) {
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
      uniform vec3 uCamRight;
      uniform vec3 uCamUp;
      uniform vec3 uCamForward;
      uniform vec3 uStrata; // (Sx, Sy, Sz)
      uniform float uYCenterOffset;
      uniform int uNumPoints;
      uniform vec4 uPoints[256]; // (u_cam, v_cam, depth_t, 0.0), sorted front-to-back

      // Compute anti-aliased grid line coverage on a plane with normal along 'faceAxis' (0, 1, or 2),
      // for subdivisions along 'subAxis' (with S strata), where 'lineAxis' is the 3rd axis.
      void evalFaceGrid(
        float coordIn01,
        float S,
        float absRdFace,
        float rdLineAxis,
        float pxWorld,
        out float outerEdgeCov,
        out float innerGridCov
      ) {
        outerEdgeCov = 0.0;
        innerGridCov = 0.0;
        if (absRdFace < 1e-4) return;

        float scalePerp = absRdFace / max(1e-4, sqrt(max(0.0, 1.0 - rdLineAxis * rdLineAxis)));

        // Distance to outer border (0 or 1) in [0, 0.5]
        float distBorder3D = 0.5 - abs(coordIn01 - 0.5);
        float distBorderPx = (distBorder3D * scalePerp) / pxWorld;
        outerEdgeCov = 1.0 - smoothstep(0.0, 1.25, distBorderPx);

        if (S > 1.5) {
          float scaled = coordIn01 * S;
          float nearestIdx = floor(scaled + 0.5);
          // Interior grid lines only (1 <= nearestIdx <= S - 1)
          if (nearestIdx >= 0.5 && nearestIdx <= S - 0.5) {
            float distGrid3D = abs(scaled - nearestIdx) / S;
            float distGridPx = (distGrid3D * scalePerp) / pxWorld;
            float faceFade = smoothstep(0.02, 0.09, absRdFace);
            innerGridCov = (1.0 - smoothstep(0.0, 1.12, distGridPx)) * faceFade;
          }
        }
      }

      // Trace interior stratification planes along 'axis' (0, 1, or 2) to render solid lines
      // at every interior plane-plane intersection inside the cube volume.
      void accumInteriorPlaneLines(
        int axis,
        float S_axis,
        float S_b,
        float S_c,
        vec3 ro,
        vec3 rd,
        vec3 absRd,
        float invD_axis,
        float cEn,
        float cEx,
        float tEn,
        float tEx,
        float tSphereHit,
        float pxWorld,
        inout float frontLineCov,
        inout float backLineCov
      ) {
        float absRdAx = (axis == 0) ? absRd.x : ((axis == 1) ? absRd.y : absRd.z);
        if (S_axis <= 1.5 || (S_b <= 1.5 && S_c <= 1.5) || absRdAx < 0.025) return;

        float lo = min(cEn, cEx) * S_axis;
        float hi = max(cEn, cEx) * S_axis;
        int kMin = max(1, int(floor(lo)));
        int kMax = min(int(S_axis) - 1, int(ceil(hi)));

        float roAxis = (axis == 0) ? ro.x : ((axis == 1) ? ro.y : ro.z);
        float rdLineB = (axis == 0) ? rd.z : ((axis == 1) ? rd.z : rd.y);
        float rdLineC = (axis == 0) ? rd.y : ((axis == 1) ? rd.x : rd.x);

        for (int k = 1; k < 64; k++) {
          if (k < kMin) continue;
          if (k > kMax) break;

          float planeCoord = -0.5 + float(k) / S_axis;
          float tP = (planeCoord - roAxis) * invD_axis;
          if (tP < tEn - pxWorld || tP > tEx + pxWorld) continue;

          vec3 pos = ro + tP * rd + 0.5;
          float uB = (axis == 0) ? pos.y : pos.x;
          float uC = (axis == 2) ? pos.y : pos.z;

          float outB = 0.0, inB = 0.0;
          float outC = 0.0, inC = 0.0;
          if (S_b > 1.5 && uC >= -0.002 && uC <= 1.002) {
            evalFaceGrid(clamp(uB, 0.0, 1.0), S_b, absRdAx, rdLineB, pxWorld, outB, inB);
          }
          if (S_c > 1.5 && uB >= -0.002 && uB <= 1.002) {
            evalFaceGrid(clamp(uC, 0.0, 1.0), S_c, absRdAx, rdLineC, pxWorld, outC, inC);
          }

          float cov = max(inB, inC);
          if (cov > 0.0) {
            if (tP <= tSphereHit) {
              frontLineCov = max(frontLineCov, cov);
            } else {
              backLineCov = max(backLineCov, cov);
            }
          }
        }
      }

      // Compute exact number of interior planes (1 .. S-1) crossed along axis a between c0 and c1 in [0, 1]
      float countPlanesCrossed(float c0, float c1, float S) {
        if (S <= 1.5) return 0.0;
        float lo = clamp(min(c0, c1), 0.0, 1.0) * S;
        float hi = clamp(max(c0, c1), 0.0, 1.0) * S;
        float kMin = max(1.0, floor(lo) + 1.0);
        float kMax = min(S - 1.0, floor(hi));
        return max(0.0, kMax - kMin + 1.0);
      }

      void main() {
        float aspect = uResolution.x / uResolution.y;
        const float orthoScale = 0.88;

        vec2 ndc = vUV * 2.0 - 1.0;
        vec2 pCam = vec2(ndc.x * aspect, ndc.y - uYCenterOffset) * orthoScale;
        float pxWorld = (2.0 * orthoScale) / uResolution.y;

        vec3 ro = pCam.x * uCamRight + pCam.y * uCamUp;
        vec3 rd = uCamForward;

        // Safe reciprocal ray direction for slab intersection with [-0.5, 0.5]^3
        vec3 safeRd = sign(rd) * max(abs(rd), vec3(1e-6));
        vec3 invD = 1.0 / safeRd;

        vec3 t0s = (-0.5 - ro) * invD;
        vec3 t1s = ( 0.5 - ro) * invD;
        vec3 tsmaller = min(t0s, t1s);
        vec3 tbigger  = max(t0s, t1s);

        float tEnter = max(max(tsmaller.x, tsmaller.y), tsmaller.z);
        float tExit  = min(min(tbigger.x, tbigger.y), tbigger.z);

        // Smooth silhouette anti-aliasing on the hull of the unit cube
        float hullMargin = pxWorld * 1.35;
        vec3 t0Pad = (-0.5 - hullMargin - ro) * invD;
        vec3 t1Pad = ( 0.5 + hullMargin - ro) * invD;
        vec3 tsPad = min(t0Pad, t1Pad);
        vec3 tbPad = max(t0Pad, t1Pad);
        float tEnterPad = max(max(tsPad.x, tsPad.y), tsPad.z);
        float tExitPad  = min(min(tbPad.x, tbPad.y), tbPad.z);

        if (tEnterPad > tExitPad) {
          fragColor = vec4(0.0);
          return;
        }

        float cubeCoverage = clamp((tExitPad - tEnterPad) / (2.0 * hullMargin), 0.0, 1.0);
        float tEn = min(tEnter, tExit);
        float tEx = max(tEnter, tExit);

        vec3 posEnter = clamp(ro + tEn * rd + 0.5, 0.0, 1.0);
        vec3 posExit  = clamp(ro + tEx * rd + 0.5, 0.0, 1.0);

        // Identify which face was entered and which face was exited
        int enterAxis = (tsmaller.x >= tsmaller.y && tsmaller.x >= tsmaller.z) ? 0 : ((tsmaller.y >= tsmaller.z) ? 1 : 2);
        int exitAxis  = (tbigger.x <= tbigger.y && tbigger.x <= tbigger.z) ? 0 : ((tbigger.y <= tbigger.z) ? 1 : 2);

        vec3 absRd = abs(rd);

        // Evaluate back-face grid lines (plane-wall intersections) and outer edges (at tEx)
        float backOuter1 = 0.0, backInner1 = 0.0;
        float backOuter2 = 0.0, backInner2 = 0.0;
        if (exitAxis == 0) {
          evalFaceGrid(posExit.y, uStrata.y, absRd.x, rd.z, pxWorld, backOuter1, backInner1);
          evalFaceGrid(posExit.z, uStrata.z, absRd.x, rd.y, pxWorld, backOuter2, backInner2);
        } else if (exitAxis == 1) {
          evalFaceGrid(posExit.x, uStrata.x, absRd.y, rd.z, pxWorld, backOuter1, backInner1);
          evalFaceGrid(posExit.z, uStrata.z, absRd.y, rd.x, pxWorld, backOuter2, backInner2);
        } else {
          evalFaceGrid(posExit.x, uStrata.x, absRd.z, rd.y, pxWorld, backOuter1, backInner1);
          evalFaceGrid(posExit.y, uStrata.y, absRd.z, rd.x, pxWorld, backOuter2, backInner2);
        }
        float backOuterCov = max(backOuter1, backOuter2);
        float backGridCov  = max(backInner1, backInner2);

        // Evaluate front-face grid lines (plane-wall intersections) and outer edges (at tEn)
        float frontOuter1 = 0.0, frontInner1 = 0.0;
        float frontOuter2 = 0.0, frontInner2 = 0.0;
        if (enterAxis == 0) {
          evalFaceGrid(posEnter.y, uStrata.y, absRd.x, rd.z, pxWorld, frontOuter1, frontInner1);
          evalFaceGrid(posEnter.z, uStrata.z, absRd.x, rd.y, pxWorld, frontOuter2, frontInner2);
        } else if (enterAxis == 1) {
          evalFaceGrid(posEnter.x, uStrata.x, absRd.y, rd.z, pxWorld, frontOuter1, frontInner1);
          evalFaceGrid(posEnter.z, uStrata.z, absRd.y, rd.x, pxWorld, frontOuter2, frontInner2);
        } else {
          evalFaceGrid(posEnter.x, uStrata.x, absRd.z, rd.y, pxWorld, frontOuter1, frontInner1);
          evalFaceGrid(posEnter.y, uStrata.y, absRd.z, rd.x, pxWorld, frontOuter2, frontInner2);
        }
        float frontOuterCov = max(frontOuter1, frontOuter2);
        float frontGridCov  = max(frontInner1, frontInner2);

        // Base back-face surface color and solid intersection lines on the back walls
        float backShade = (exitAxis == 1) ? 0.148 : ((exitAxis == 0) ? 0.135 : 0.122);
        vec3 backCol = vec3(backShade);
        backCol = mix(backCol, vec3(0.35, 0.36, 0.39), backGridCov * 0.85);
        backCol = mix(backCol, vec3(0.45, 0.46, 0.49), backOuterCov * 0.88);

        // Per-plane opacity scaled gently with total plane count
        float totalPlanes = max(1.0, (uStrata.x - 1.0) + (uStrata.y - 1.0) + (uStrata.z - 1.0));
        float basePlaneAlpha = clamp(0.18 / sqrt(totalPlanes), 0.015, 0.095);

        // Dark neutral smoked-glass plane colors (distinct from the bright cyan points)
        vec3 planeColX = vec3(0.08, 0.09, 0.10);
        vec3 planeColY = vec3(0.11, 0.12, 0.13);
        vec3 planeColZ = vec3(0.09, 0.10, 0.11);

        // Find closest 3D point sphere hit along the ray (uPoints are sorted front-to-back)
        float ptRadius = (uNumPoints <= 16) ? 0.030 : ((uNumPoints <= 64) ? 0.024 : 0.0185);
        float aaRadius = pxWorld * 0.85;
        float outerR = ptRadius + aaRadius;
        float outerR2 = outerR * outerR;

        float sphereCov = 0.0;
        float tSphereHit = tEx;
        vec3 sphereColor = vec3(0.0);

        for (int i = 0; i < 256; i++) {
          if (i >= uNumPoints) break;
          vec3 pt = uPoints[i].xyz;
          vec2 dp = pCam - pt.xy;
          float d2 = dot(dp, dp);
          if (d2 < outerR2) {
            float d = sqrt(d2);
            float cov = 1.0 - smoothstep(ptRadius - aaRadius, ptRadius + aaRadius, d);
            if (cov > 0.001) {
              float zLocal = sqrt(max(0.0, ptRadius * ptRadius - d2));
              float tHit = pt.z - zLocal;

              // Crisp dark rim + vivid cyan (#19d4fa) core with 3D highlight
              float rim = smoothstep(ptRadius - pxWorld * 1.15, ptRadius, d);
              vec2 n2 = dp / max(ptRadius, 1e-4);
              float nz = zLocal / max(ptRadius, 1e-4);
              float highlight = max(0.0, -0.35 * n2.x + 0.45 * n2.y + 0.82 * nz);
              vec3 coreCol = mix(vec3(0.10, 0.78, 0.95), vec3(0.32, 0.94, 1.0), highlight);
              vec3 cCol = mix(coreCol, vec3(0.07, 0.07, 0.09), rim * 0.88);

              // Composite front-to-back if first sphere only partially covers pixel
              if (sphereCov == 0.0) {
                tSphereHit = clamp(tHit, tEn, tEx);
              }
              float rem = 1.0 - sphereCov;
              sphereColor += cCol * (cov * rem);
              sphereCov += cov * rem;
              if (sphereCov > 0.99) break;
            }
          }
        }
        if (sphereCov > 0.0) {
          sphereColor /= sphereCov;
        }

        // Evaluate all interior plane-plane intersection lines along X, Y, Z planes
        float frontInteriorLineCov = 0.0;
        float backInteriorLineCov  = 0.0;
        accumInteriorPlaneLines(
          0, uStrata.x, uStrata.y, uStrata.z,
          ro, rd, absRd, invD.x, posEnter.x, posExit.x,
          tEn, tEx, tSphereHit, pxWorld,
          frontInteriorLineCov, backInteriorLineCov
        );
        accumInteriorPlaneLines(
          1, uStrata.y, uStrata.x, uStrata.z,
          ro, rd, absRd, invD.y, posEnter.y, posExit.y,
          tEn, tEx, tSphereHit, pxWorld,
          frontInteriorLineCov, backInteriorLineCov
        );
        accumInteriorPlaneLines(
          2, uStrata.z, uStrata.x, uStrata.y,
          ro, rd, absRd, invD.z, posEnter.z, posExit.z,
          tEn, tEx, tSphereHit, pxWorld,
          frontInteriorLineCov, backInteriorLineCov
        );

        // Compute interior planes crossed along the full ray [tEn, tEx]
        float nxFull = countPlanesCrossed(posEnter.x, posExit.x, uStrata.x);
        float nyFull = countPlanesCrossed(posEnter.y, posExit.y, uStrata.y);
        float nzFull = countPlanesCrossed(posEnter.z, posExit.z, uStrata.z);

        float transFull = pow(1.0 - basePlaneAlpha, nxFull + nyFull + nzFull);
        float wSumFull = max(1e-4, nxFull + nyFull + nzFull);
        vec3 avgPlaneColFull = (nxFull * planeColX + nyFull * planeColY + nzFull * planeColZ) / wSumFull;
        vec3 bgWithPlanes = mix(avgPlaneColFull, backCol, transFull);

        // Solid interior plane-plane intersection lines behind any hit sphere
        float interiorLineAlpha = (totalPlanes > 12.0) ? 0.52 : 0.68;
        bgWithPlanes = mix(bgWithPlanes, vec3(0.40, 0.41, 0.44), backInteriorLineCov * interiorLineAlpha);

        // If a sphere was hit, attenuate it clearly for each dark plane strictly in front of it [tEn, tSphereHit]
        vec3 finalInterior = bgWithPlanes;
        if (sphereCov > 0.0) {
          vec3 posSphere = clamp(ro + tSphereHit * rd + 0.5, 0.0, 1.0);
          float nxFront = countPlanesCrossed(posEnter.x, posSphere.x, uStrata.x);
          float nyFront = countPlanesCrossed(posEnter.y, posSphere.y, uStrata.y);
          float nzFront = countPlanesCrossed(posEnter.z, posSphere.z, uStrata.z);

          float spherePlaneAtten = clamp(0.30 / sqrt(max(1.0, totalPlanes / 3.0)), 0.11, 0.26);
          float transFront = pow(1.0 - spherePlaneAtten, nxFront + nyFront + nzFront);
          float wSumFront = max(1e-4, nxFront + nyFront + nzFront);
          vec3 avgPlaneColFront = (nxFront * planeColX + nyFront * planeColY + nzFront * planeColZ) / wSumFront;
          vec3 sphereWithFrontPlanes = mix(avgPlaneColFront, sphereColor, transFront);

          finalInterior = mix(bgWithPlanes, sphereWithFrontPlanes, sphereCov);
        }

        // Composite solid interior plane-plane intersection lines in front of the sphere
        vec3 col = mix(finalInterior, vec3(0.42, 0.43, 0.46), frontInteriorLineCov * interiorLineAlpha);

        // Composite front-face stratification grid lines (plane-wall intersections) & crisp front outer edges
        float frontGridAlpha = (totalPlanes > 12.0) ? 0.62 : 0.78;
        col = mix(col, vec3(0.50, 0.51, 0.54), frontGridCov * frontGridAlpha);
        col = mix(col, vec3(0.76, 0.76, 0.78), frontOuterCov * 0.95);

        fragColor = vec4(col * cubeCoverage, cubeCoverage);
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

      glProg = gl.createProgram();
      gl.attachShader(glProg, compileShader(gl.VERTEX_SHADER, vsSource));
      gl.attachShader(glProg, compileShader(gl.FRAGMENT_SHADER, fsSource));
      gl.linkProgram(glProg);

      glVao = gl.createVertexArray();
      uLocs = {
        resolution: gl.getUniformLocation(glProg, "uResolution"),
        camRight: gl.getUniformLocation(glProg, "uCamRight"),
        camUp: gl.getUniformLocation(glProg, "uCamUp"),
        camForward: gl.getUniformLocation(glProg, "uCamForward"),
        strata: gl.getUniformLocation(glProg, "uStrata"),
        yCenterOffset: gl.getUniformLocation(glProg, "uYCenterOffset"),
        numPoints: gl.getUniformLocation(glProg, "uNumPoints"),
        points: gl.getUniformLocation(glProg, "uPoints")
      };
    }

    const DEFAULT_YAW = 0.54;
    const DEFAULT_PITCH = 0.42;
    let yaw = DEFAULT_YAW;
    let pitch = DEFAULT_PITCH;
    let isVisible = true;
    let rafId = null;

    const pointUniformBuf = new Float32Array(256 * 4);

    function render3DAnd4DCubes() {
      rafId = null;
      if (!gl || !glProg || viewMode === "2D" || !isVisible) return;

      const dpr = Math.min(window.devicePixelRatio || 1, 2);
      const rect = cubeStageWrap.getBoundingClientRect();
      const w = Math.max(1, Math.round(rect.width * dpr));
      const h = Math.max(1, Math.round(rect.height * dpr));

      if (glCanvas.width !== w || glCanvas.height !== h) {
        glCanvas.width = w;
        glCanvas.height = h;
      }

      gl.viewport(0, 0, w, h);
      gl.clearColor(0.0, 0.0, 0.0, 0.0);
      gl.clear(gl.COLOR_BUFFER_BIT);

      gl.useProgram(glProg);
      gl.bindVertexArray(glVao);

      const cp = Math.cos(pitch),
        sp = Math.sin(pitch);
      const cy = Math.cos(yaw),
        sy = Math.sin(yaw);

      const fx = -sy * cp,
        fy = -sp,
        fz = -cy * cp;
      const rx = cy,
        ry = 0,
        rz = -sy;
      const ux = -sy * sp,
        uy = cp,
        uz = -cy * sp;

      gl.uniform3f(uLocs.camRight, rx, ry, rz);
      gl.uniform3f(uLocs.camUp, ux, uy, uz);
      gl.uniform3f(uLocs.camForward, fx, fy, fz);
      gl.uniform1f(uLocs.yCenterOffset, viewMode === "3D" ? 0.0 : -0.05);

      const n = 1 << currentK;
      const cellW = Math.floor(w / 2);
      const cellH = Math.floor(h / 2);

      // Build viewport descriptors for active mode
      let cellsToDraw = [];
      if (viewMode === "3D") {
        const m4 = Math.ceil(currentK / 2);
        const { p1, p2, p3 } = compute3DKnobs(m4);
        const s1 = 1 << (2 * p1);
        const s2 = 1 << (2 * p2);
        const s3 = 1 << (2 * p3);

        const triplets = [
          { dims: [0, 1, 2], col: 0, row: 0 }, // top-left: xyz
          { dims: [0, 1, 3], col: 1, row: 0 }, // top-right: xyw
          { dims: [0, 2, 3], col: 0, row: 1 }, // bottom-left: xzw
          { dims: [1, 2, 3], col: 1, row: 1 }  // bottom-right: yzw
        ];

        cellsToDraw = triplets.map((t) => ({
          vx: t.col === 0 ? 0 : cellW,
          vy: t.row === 0 ? h - cellH : 0,
          vw: t.col === 0 ? cellW : w - cellW,
          vh: t.row === 0 ? cellH : h - cellH,
          dims: t.dims,
          strata: [s1, s2, s3],
          wFilter: null
        }));
      } else {
        // 4D Mode
        if (n <= 16) {
          // 2 rows of 1 centered cube: w in [0, 0.5) on top, w in [0.5, 1.0) on bottom, each 2x2x2 in xyz
          const centerX = Math.floor((w - cellW) / 2);
          cellsToDraw = [
            {
              vx: centerX,
              vy: h - cellH,
              vw: cellW,
              vh: cellH,
              dims: [0, 1, 2],
              strata: [2, 2, 2],
              wFilter: [0.0, 0.5]
            },
            {
              vx: centerX,
              vy: 0,
              vw: cellW,
              vh: h - cellH,
              dims: [0, 1, 2],
              strata: [2, 2, 2],
              wFilter: [0.5, 1.01]
            }
          ];
        } else {
          // 2 rows of 2 cubes: 4 quarter-intervals of w, each 4x4x4 in xyz
          const quarters = [
            { col: 0, row: 0, wMin: 0.0, wMax: 0.25 },
            { col: 1, row: 0, wMin: 0.25, wMax: 0.5 },
            { col: 0, row: 1, wMin: 0.5, wMax: 0.75 },
            { col: 1, row: 1, wMin: 0.75, wMax: 1.01 }
          ];
          cellsToDraw = quarters.map((q) => ({
            vx: q.col === 0 ? 0 : cellW,
            vy: q.row === 0 ? h - cellH : 0,
            vw: q.col === 0 ? cellW : w - cellW,
            vh: q.row === 0 ? cellH : h - cellH,
            dims: [0, 1, 2],
            strata: [4, 4, 4],
            wFilter: [q.wMin, q.wMax]
          }));
        }
      }

      cellsToDraw.forEach((cell) => {
        gl.viewport(cell.vx, cell.vy, cell.vw, cell.vh);
        gl.uniform2f(uLocs.resolution, cell.vw, cell.vh);
        gl.uniform3f(
          uLocs.strata,
          cell.strata[0],
          cell.strata[1],
          cell.strata[2]
        );

        const [d0, d1, d2] = cell.dims;
        const projectedPts = [];

        for (let i = 0; i < n; i++) {
          const p = points4D[i];
          if (cell.wFilter !== null) {
            const wVal = p[3];
            if (wVal < cell.wFilter[0] || wVal >= cell.wFilter[1]) continue;
          }
          const px = p[d0] - 0.5;
          const py = p[d1] - 0.5;
          const pz = p[d2] - 0.5;

          const uCam = px * rx + py * ry + pz * rz;
          const vCam = px * ux + py * uy + pz * uz;
          const tDepth = px * fx + py * fy + pz * fz;
          projectedPts.push({ uCam, vCam, tDepth });
        }

        // Sort front-to-back along the camera ray (smallest tDepth first)
        projectedPts.sort((a, b) => a.tDepth - b.tDepth);

        const count = projectedPts.length;
        for (let i = 0; i < count; i++) {
          const base = i * 4;
          pointUniformBuf[base] = projectedPts[i].uCam;
          pointUniformBuf[base + 1] = projectedPts[i].vCam;
          pointUniformBuf[base + 2] = projectedPts[i].tDepth;
          pointUniformBuf[base + 3] = 0.0;
        }

        gl.uniform1i(uLocs.numPoints, count);
        gl.uniform4fv(uLocs.points, pointUniformBuf);
        gl.drawArrays(gl.TRIANGLES, 0, 3);
      });
    }

    function request3DRender() {
      if (rafId === null && isVisible && viewMode !== "2D") {
        rafId = requestAnimationFrame(render3DAnd4DCubes);
      }
    }

    // Pointer rotation on the 3D/4D canvas
    let draggingCanvas = false;
    let lastX = 0;
    let lastY = 0;

    glCanvas.addEventListener("pointerdown", (e) => {
      if (e.button !== 0 && e.pointerType === "mouse") return;
      draggingCanvas = true;
      lastX = e.clientX;
      lastY = e.clientY;
      glCanvas.setPointerCapture(e.pointerId);
    });

    glCanvas.addEventListener("pointermove", (e) => {
      if (!draggingCanvas) return;
      const dx = e.clientX - lastX;
      const dy = e.clientY - lastY;
      lastX = e.clientX;
      lastY = e.clientY;

      yaw -= dx * 0.007;
      pitch = Math.max(-1.45, Math.min(1.45, pitch + dy * 0.007));

      cubeCells.forEach((c) => c.gizmo.update(yaw, pitch));
      request3DRender();
    });

    function endCanvasDrag(e) {
      if (!draggingCanvas) return;
      draggingCanvas = false;
      if (
        e &&
        e.pointerId !== undefined &&
        glCanvas.hasPointerCapture(e.pointerId)
      ) {
        glCanvas.releasePointerCapture(e.pointerId);
      }
    }

    glCanvas.addEventListener("pointerup", endCanvasDrag);
    glCanvas.addEventListener("pointercancel", endCanvasDrag);

    glCanvas.addEventListener("dblclick", () => {
      yaw = DEFAULT_YAW;
      pitch = DEFAULT_PITCH;
      cubeCells.forEach((c) => c.gizmo.update(yaw, pitch));
      request3DRender();
    });

    function update2DSubplots() {
      const n = 1 << currentK;
      const log2X = computeLogXStrata(currentK, internalOffset2D);
      const log2Y = currentK - log2X;
      const b2X = 1 << log2X;
      const b2Y = 1 << log2Y;

      // Round up base-4 strata so two points never share the same stratum:
      // - When currentK is odd (N = 2, 8, 32, 128), one of (log2X, log2Y) is odd and the other is even;
      //   rounding the odd exponent up by +1 to the next even integer maps e.g. 8x1 -> 16x1, 4x2 -> 4x4, 8x4 -> 16x4.
      // - When currentK is even (N = 4, 16, 64, 256), m4 = currentK / 2 has 4^m4 = N strata (1 point per cell).
      let log4X, log4Y;
      if ((currentK & 1) === 1) {
        log4X = (log2X + (log2X & 1)) >> 1;
        log4Y = (log2Y + (log2Y & 1)) >> 1;
      } else {
        const m4 = currentK >> 1;
        log4X = computeLogXStrata(m4, internalOffset2D / 2.0);
        log4Y = m4 - log4X;
      }
      const b4X = 1 << (2 * log4X);
      const b4Y = 1 << (2 * log4Y);

      const sVal = currentK - log2X;
      strata2DSlider.input.max = String(currentK);
      strata2DSlider.input.value = String(sVal);
      strata2DSlider.renderTicks(currentK + 1);

      const strataPct = currentK > 0 ? (sVal / currentK) * 100 : 0;
      strata2DSlider.fill.style.left = "0%";
      strata2DSlider.fill.style.width = `${strataPct}%`;
      strata2DSlider.valLabel.innerHTML = `${b2X}&thinsp;&times;&thinsp;${b2Y}`;

      const radius =
        currentK <= 4 ? "2.7" : currentK <= 6 ? "2.2" : currentK === 7 ? "1.8" : "1.5";

      subplots2D.forEach(({ proj, bx, by, gridG, circles, gridSpan }) => {
        const sx = proj.isBase2 ? b2X : b4X;
        const sy = proj.isBase2 ? b2Y : b4Y;

        gridSpan.textContent = `(${sx}\u2009\u00d7\u2009${sy})`;

        gridG.innerHTML = "";
        const maxS = Math.max(sx, sy);
        const strokeColor = maxS >= 64 ? "#44444a" : "#56565c";
        const strokeWidth = maxS >= 64 ? "0.65" : "0.9";

        for (let ix = 1; ix < sx; ix++) {
          const xPos = bx + (ix / sx) * subSize;
          const line = document.createElementNS(svgNS, "line");
          line.setAttribute("x1", xPos.toFixed(2));
          line.setAttribute("y1", String(by));
          line.setAttribute("x2", xPos.toFixed(2));
          line.setAttribute("y2", String(by + subSize));
          line.setAttribute("stroke", strokeColor);
          line.setAttribute("stroke-width", strokeWidth);
          gridG.appendChild(line);
        }

        for (let iy = 1; iy < sy; iy++) {
          const yPos = by + (iy / sy) * subSize;
          const line = document.createElementNS(svgNS, "line");
          line.setAttribute("x1", String(bx));
          line.setAttribute("y1", yPos.toFixed(2));
          line.setAttribute("x2", String(bx + subSize));
          line.setAttribute("y2", yPos.toFixed(2));
          line.setAttribute("stroke", strokeColor);
          line.setAttribute("stroke-width", strokeWidth);
          gridG.appendChild(line);
        }

        for (let i = 0; i < MAX_POINTS; i++) {
          const c = circles[i];
          if (i < n) {
            c.style.display = "";
            c.setAttribute("r", radius);
          } else {
            c.style.display = "none";
          }
        }
      });
    }

    function updateCubeOverlayAndControls() {
      const n = 1 << currentK;
      const m4 = Math.ceil(currentK / 2);

      if (viewMode === "3D") {
        cubeOverlay.classList.remove("single-col");
        strata3DSlider.renderState(m4);

        const triplets = [
          { axes: ["x", "y", "z"] },
          { axes: ["x", "y", "w"] },
          { axes: ["x", "z", "w"] },
          { axes: ["y", "z", "w"] }
        ];

        cubeCells.forEach((c, idx) => {
          c.cell.style.display = "";
          c.header.innerHTML = "";
          c.gizmo.setAxisNames(triplets[idx].axes);
          c.gizmo.update(yaw, pitch);
        });
      } else if (viewMode === "4D") {
        if (n <= 16) {
          cubeOverlay.classList.add("single-col");

          const labels = [
            `<i>w</i> &lt; 0.5 <span class="cube-cell-sub">(2&thinsp;&times;&thinsp;2&thinsp;&times;&thinsp;2)</span>`,
            `<i>w</i> &ge; 0.5 <span class="cube-cell-sub">(2&thinsp;&times;&thinsp;2&thinsp;&times;&thinsp;2)</span>`
          ];

          cubeCells.forEach((c, idx) => {
            if (idx < 2) {
              c.cell.style.display = "";
              c.header.innerHTML = labels[idx];
              c.gizmo.setAxisNames(["x", "y", "z"]);
              c.gizmo.update(yaw, pitch);
            } else {
              c.cell.style.display = "none";
            }
          });
        } else {
          cubeOverlay.classList.remove("single-col");

          const labels = [
            `0 &le; <i>w</i> &lt; 0.25 <span class="cube-cell-sub">(4&thinsp;&times;&thinsp;4&thinsp;&times;&thinsp;4)</span>`,
            `0.25 &le; <i>w</i> &lt; 0.5 <span class="cube-cell-sub">(4&thinsp;&times;&thinsp;4&thinsp;&times;&thinsp;4)</span>`,
            `0.5 &le; <i>w</i> &lt; 0.75 <span class="cube-cell-sub">(4&thinsp;&times;&thinsp;4&thinsp;&times;&thinsp;4)</span>`,
            `0.75 &le; <i>w</i> &lt; 1 <span class="cube-cell-sub">(4&thinsp;&times;&thinsp;4&thinsp;&times;&thinsp;4)</span>`
          ];

          cubeCells.forEach((c, idx) => {
            c.cell.style.display = "";
            c.header.innerHTML = labels[idx];
            c.gizmo.setAxisNames(["x", "y", "z"]);
            c.gizmo.update(yaw, pitch);
          });
        }
      }
    }

    function updateAll() {
      const n = 1 << currentK;
      const samplesPct = ((currentK - 1) / (8 - 1)) * 100;
      samplesSlider.fill.style.left = "0%";
      samplesSlider.fill.style.width = `${samplesPct}%`;
      samplesSlider.valLabel.innerHTML = `${n}`;

      svg2D.style.display = viewMode === "2D" ? "" : "none";
      cubeStageWrap.style.display = viewMode === "2D" ? "none" : "";

      strata2DSlider.row.style.display = viewMode === "2D" ? "" : "none";
      strata3DSlider.row.style.display = viewMode === "3D" ? "" : "none";

      if (viewMode === "2D") {
        update2DSubplots();
      } else {
        updateCubeOverlayAndControls();
        request3DRender();
      }
    }

    regenBtn.addEventListener("click", () => {
      points4D = generateStochasticSZ4D(MAX_POINTS);
      sync2DPointPositions();
      updateAll();
    });

    samplesSlider.input.addEventListener("input", () => {
      currentK = parseInt(samplesSlider.input.value, 10);
      updateAll();
    });

    samplesSlider.input.addEventListener("dblclick", () => {
      currentK = 4;
      samplesSlider.input.value = "4";
      updateAll();
    });

    strata2DSlider.input.addEventListener("input", () => {
      const sVal = parseInt(strata2DSlider.input.value, 10);
      const log2X = currentK - sVal;
      internalOffset2D = currentK / 2.0 - log2X;
      updateAll();
    });

    strata2DSlider.input.addEventListener("dblclick", () => {
      internalOffset2D = 0.0;
      updateAll();
    });

    if ("IntersectionObserver" in window) {
      const io = new IntersectionObserver(
        (entries) => {
          for (const entry of entries) {
            isVisible = entry.isIntersecting;
            if (isVisible) request3DRender();
          }
        },
        { threshold: 0.01 }
      );
      io.observe(container);
    }

    if ("ResizeObserver" in window) {
      const ro = new ResizeObserver(() => {
        request3DRender();
      });
      ro.observe(cubeStageWrap);
    } else {
      window.addEventListener("resize", request3DRender, { passive: true });
    }

    updateAll();
  }

  function initAllDiagrams() {
    const elStrat = document.getElementById("sz-diagram-02-stratification");
    if (elStrat) initStratification02Diagram(elStrat);

    const el4D = document.getElementById("sz-diagram-4d-projections");
    if (el4D) initSZ4DProjectionsDiagram(el4D);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initAllDiagrams);
  } else {
    initAllDiagrams();
  }
})();
