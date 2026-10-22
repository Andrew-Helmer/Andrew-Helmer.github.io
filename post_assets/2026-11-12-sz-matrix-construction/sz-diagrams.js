(function () {
  "use strict";

  // Embedded y_xor values for the simple stochastic Sobol (0,2) algorithm
  // from "Stochastic Generation of (t,s) Sample Sequences" (Helmer, Christensen, Kensler 2021).
  // Entry m corresponds to column m of P^{-1} (strictly above the main diagonal) packed into bits 0..m-1.
  const SOBOL_02_Y_XORS = [0x0, 0x1, 0x1, 0x7, 0x1, 0x13, 0x15, 0x7f];

  function randomUint32() {
    if (typeof crypto !== "undefined" && crypto.getRandomValues) {
      const buf = new Uint32Array(1);
      crypto.getRandomValues(buf);
      return buf[0] >>> 0;
    }
    return (Math.random() * 0x100000000) >>> 0;
  }

  // Generate `n` points (up to 256) of a stochastic base-2 Sobol (0,2)-sequence.
  // Because the x generator matrix is I (all x_xor = 0), x is a 1D van der Corput
  // sequence and y pairs each new point with previous point `k ^ y_xor[m]`.
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

    // Background square (#1e1e1e)
    const bgRect = document.createElementNS(svgNS, "rect");
    bgRect.setAttribute("x", String(boxX));
    bgRect.setAttribute("y", String(boxY));
    bgRect.setAttribute("width", String(boxSize));
    bgRect.setAttribute("height", String(boxSize));
    bgRect.setAttribute("fill", "#1e1e1e");
    svg.appendChild(bgRect);

    // Group for interior stratification grid lines
    const gridGroup = document.createElementNS(svgNS, "g");
    svg.appendChild(gridGroup);

    // Group for sample points
    const pointsGroup = document.createElementNS(svgNS, "g");
    svg.appendChild(pointsGroup);

    // Sharp mitered outer border (#bebebe, 1.35px) drawn on top of interior grid lines
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

    // Pre-create 256 SVG circle elements so slider updates only toggle visibility and radius
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

    // State variables matching the specification:
    // Slider 1: power-of-two sample count from 2^1 = 2 to 2^8 = 256 (8 positions, default 4th position -> k = 4, N = 16).
    // Slider 2: stratification grid lines, internally storing `offset` from the middle stratification, starting at 0.
    let currentK = 4; // log2(n) in [1..8], position 4 is 2^4 = 16 points
    let internalOffset = 0.0; // offset = log2(n)/2.0 - log2(x_strata)

    // Calculate effective log2(x_strata) from currentK and internalOffset without mutating internalOffset
    function computeLog2XStrata(k, offset) {
      let effectiveOffset = offset;
      const isOddPower = (k & 1) === 1;
      const isWholeOffset = Math.abs(offset - Math.round(offset)) < 1e-6;

      if (isOddPower && isWholeOffset) {
        // Whole-number offset with an odd power of two rounds down to the lower 0.5 effective value
        // (e.g. offset 0 -> -0.5, giving 4x2 rather than 2x4 for 8 points).
        effectiveOffset = offset - 0.5;
      } else if (!isOddPower && !isWholeOffset) {
        // x.5 offset with an even power of two rounds towards 0
        // (e.g. -0.5 or +0.5 for 8 points both round to 0 -> 4x4 when moving to 16 points).
        effectiveOffset = Math.trunc(offset);
      }

      const rawLog2X = k / 2.0 - effectiveOffset;
      return Math.max(0, Math.min(k, Math.round(rawLog2X)));
    }

    // Build sliders UI
    const slidersWrap = document.createElement("div");
    slidersWrap.className = "diagram-sliders";
    slidersWrap.style.maxWidth = "340px";

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

    // Slider 1: Sample count (k = 1..8 -> N = 2..256, 8 discrete steps)
    const pointsSlider = createDiscreteSliderRow(
      "Samples",
      1,
      8,
      currentK,
      "Power-of-two sample count"
    );
    pointsSlider.renderTicks(8);

    // Slider 2: Stratification (index = 0..k, where index = k - log2(x_strata) so moving left->right goes from Nx1 to 1xN, or log2(x_strata) = k - s)
    // Wait: let's check slider direction!
    const strataSlider = createDiscreteSliderRow(
      "Strata",
      0,
      currentK,
      2,
      "Stratification grid aspect"
    );

    // Regenerate button placed compactly between the plot box and the sliders
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
      const log2X = computeLog2XStrata(currentK, internalOffset);
      const log2Y = currentK - log2X;
      const xStrata = 1 << log2X;
      const yStrata = 1 << log2Y;

      // Update points slider UI
      const pointsPct = ((currentK - 1) / (8 - 1)) * 100;
      pointsSlider.fill.style.left = "0%";
      pointsSlider.fill.style.width = `${pointsPct}%`;
      pointsSlider.valLabel.innerHTML = `${n}`;

      // Update stratification slider range, ticks, and thumb position
      // We map slider value s in [0, currentK] to log2(x_strata) = currentK - s
      // so s = 0 is (N x 1), middle s = currentK/2 is square (sqrt(N) x sqrt(N)), and s = currentK is (1 x N).
      // Notice that offset = currentK / 2.0 - log2(x_strata) = s - currentK / 2.0!
      const sVal = currentK - log2X;
      strataSlider.input.max = String(currentK);
      strataSlider.input.value = String(sVal);
      strataSlider.renderTicks(currentK + 1);

      const strataPct = currentK > 0 ? (sVal / currentK) * 100 : 0;
      strataSlider.fill.style.left = "0%";
      strataSlider.fill.style.width = `${strataPct}%`;
      strataSlider.valLabel.innerHTML = `${xStrata}&thinsp;&times;&thinsp;${yStrata}`;

      // Redraw interior stratification grid lines
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

      // Scale point radius slightly with sample count so 256 points remain crisp and uncrowded
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

    // When moving the power-of-two slider, the internal offset value is NOT recalculated.
    pointsSlider.input.addEventListener("input", () => {
      currentK = parseInt(pointsSlider.input.value, 10);
      updateDiagram();
    });

    pointsSlider.input.addEventListener("dblclick", () => {
      currentK = 4;
      pointsSlider.input.value = "4";
      updateDiagram();
    });

    // When moving the stratification slider, the internal offset IS recalculated as:
    // offset = log2(n) / 2.0 - log2(x_strata)
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

  function initAllDiagrams() {
    const elStrat = document.getElementById("sz-diagram-02-stratification");
    if (elStrat) initStratification02Diagram(elStrat);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initAllDiagrams);
  } else {
    initAllDiagrams();
  }
})();

