let etaSSlider, etaSValDisplay, lcToggle, canvasOriginal, canvasBoosted, canvasProjective, canvasTime, mathMatrix;
let etaRSlider, etaRValDisplay, ampSlider, ampValDisplay, lemniToggle;
let mmSlider, mmValDisplay;
let ctxOrig, ctxBoost, ctxProj, ctxTime;
let etaS = 0;
let time = 0;
let spacelike = false;
let lemniscate = false;

// --- Motion of the point ------------------------------------------------
let etaR = 0.2;   // centre of the oscillation, eta_r
let etaAmp = 0.2; // amplitude, i.e. Delta_max
const WINDOW_T = 10; // seconds shown in the time chart

// The radial mode obeys
//     (1/2) etaDot^2 = -M cosh 4 Delta + alpha_1D^2 / 4,   Delta = eta - eta_r,
// whose turning points satisfy cosh 4 Delta_max = alpha_1D^2 / (4 M). Fixing M
// and reading the amplitude off the slider therefore fixes alpha_1D^2, and
// differentiating gives the second-order form integrated below.
let MM = 0.1; // the invariant M; sets the overall time scale

function alphaSq() { return 4 * MM * Math.cosh(4 * etaAmp); }

// Complete elliptic integral of the first kind by the arithmetic-geometric mean
function ellipticK(k) {
    let a = 1, b = Math.sqrt(Math.max(1 - k * k, 0));
    for (let i = 0; i < 40 && Math.abs(a - b) > 1e-15; i++) {
        const an = 0.5 * (a + b);
        b = Math.sqrt(a * b);
        a = an;
    }
    return Math.PI / (2 * a);
}

// With alpha^2 = 4 M cosh 4 eta_0 the period collapses to a function of the
// amplitude alone, up to the overall scale sqrt(M):
//     T = K(tanh 2 eta_0) / (sqrt(M) cosh 2 eta_0),  eta_0 = Delta_max.
function period() {
    return ellipticK(Math.tanh(2 * etaAmp)) / (Math.sqrt(MM) * Math.cosh(2 * etaAmp));
}
function accel(eta) { return -4 * MM * Math.sinh(4 * (eta - etaR)); }

let etaState = etaR + etaAmp;
let etaDot = 0;
let history = [];

function resetMotion() {
    etaState = etaR + etaAmp; // start at a turning point
    etaDot = 0;
    history = [{ t: time, eta: etaState }];
}

function stepMotion(dt) {
    const sub = 16, h = dt / sub;
    for (let i = 0; i < sub; i++) {
        const a1 = accel(etaState);
        etaState += etaDot * h + 0.5 * a1 * h * h;
        const a2 = accel(etaState);
        etaDot += 0.5 * (a1 + a2) * h;
    }
}

// --- Sectors ------------------------------------------------------------
// A point is given by (s, u): s = spatial component, u = temporal component.
// Timelike sector : drawn as is, the particle runs on u^2 - s^2 = 1.
// Spacelike sector: everything is reflected in the null line u = s, i.e.
//                   (s, u) -> (u, s), so the particle runs on s^2 - u^2 = 1.
// The reflection fixes the light cone pointwise on u = s and as a set on
// u = -s, and it commutes with the boost, which is why the same Lambda
// generates the motion in either sector.

function project(s, u) {
    return spacelike ? [u, s] : [s, u];
}

function viewOrigin(width, height) {
    return [width / 2, height / 2];
}

// Tied to the canvas size, so raising the drawing resolution enlarges the
// picture instead of padding it with empty margin.
function viewScale(width, height) {
    return Math.min(width, height) / 8;
}

// Update UI and Math
function updateUI() {
    etaSValDisplay.textContent = etaS.toFixed(2);
    if (etaRValDisplay) etaRValDisplay.textContent = etaR.toFixed(2);
    if (ampValDisplay) ampValDisplay.textContent = etaAmp.toFixed(2);
    if (mmValDisplay) mmValDisplay.textContent = MM.toFixed(3);
    const invBox = document.getElementById('invariants');
    if (invBox && typeof katex !== 'undefined') {
        try {
            katex.render(
                `\\begin{aligned} \\alpha_{\\mathrm{1D}}^{2} &= 4\\mathcal{M}\\cosh 4\\eta_0 = ${alphaSq().toFixed(3)} \\\\[0.5ex] \\mathcal{T} &= \\frac{K(\\tanh 2\\eta_0)}{\\sqrt{\\mathcal{M}}\\,\\cosh 2\\eta_0} = ${period().toFixed(3)} \\end{aligned}`,
                invBox, { throwOnError: false, displayMode: true });
        } catch (e) { /* ignore */ }
    }

    const ch = Math.cosh(etaS).toFixed(2);
    const sh = Math.sinh(etaS).toFixed(2);

    const matrixTex = `\\begin{aligned} x' &= \\Lambda(\\eta_s)\\,x \\\\[0.3ex] e'_a &= \\Lambda(-\\eta_s)\\,e_a \\\\[0.8ex] \\Lambda(\\eta_s) &= \\begin{pmatrix} \\cosh \\eta_s & \\sinh \\eta_s \\\\ \\sinh \\eta_s & \\cosh \\eta_s \\end{pmatrix} \\\\[0.5ex] &\\approx \\begin{pmatrix} ${ch} & ${sh} \\\\ ${sh} & ${ch} \\end{pmatrix} \\end{aligned}`;
    if (typeof katex !== 'undefined' && mathMatrix) {
        try {
            katex.render(matrixTex, mathMatrix, {
                throwOnError: false,
                displayMode: true
            });
        } catch (e) {
            console.error("KaTeX rendering error:", e);
        }
    }
}

// --- Drawing utilities --------------------------------------------------

function makePlotter(ctx, width, height) {
    const [cx, cy] = viewOrigin(width, height);
    const scale = viewScale(width, height);
    return {
        cx, cy, scale,
        // screen position of the point with components (s, u)
        pt(s, u) {
            const [h, v] = project(s, u);
            return [cx + h * scale, cy - v * scale];
        }
    };
}

function strokePath(ctx, points, color, lineWidth, dash) {
    if (!points.length) return;
    ctx.strokeStyle = color;
    ctx.lineWidth = lineWidth;
    if (dash) ctx.setLineDash(dash);
    ctx.beginPath();
    points.forEach(([x, y], i) => (i === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y)));
    ctx.stroke();
    if (dash) ctx.setLineDash([]);
}

function drawGrid(ctx, P, reach) {
    // The axes of the rest frame: the time axis (s = 0) and the space axis (u = 0).
    const faint = "rgba(255, 255, 255, 0.2)";
    strokePath(ctx, [P.pt(-reach, 0), P.pt(reach, 0)], faint, 1);
    strokePath(ctx, [P.pt(0, -reach), P.pt(0, reach)], faint, 1);
}

function drawNullLines(ctx, P, reach) {
    const yellow = "#eab308";
    strokePath(ctx, [P.pt(-reach, -reach), P.pt(reach, reach)], yellow, 1.5, [5, 5]);
    strokePath(ctx, [P.pt(reach, -reach), P.pt(-reach, reach)], yellow, 1.5, [5, 5]);
}

function drawHyperbola(ctx, P) {
    // u^2 - s^2 = 1, upper branch, parametrized by rapidity
    const pts = [];
    for (let e = -3.6; e <= 3.6; e += 0.05) {
        pts.push(P.pt(Math.sinh(e), Math.cosh(e)));
    }
    strokePath(ctx, pts, "rgba(255, 255, 255, 0.5)", 2);
}

function drawVector(ctx, P, s, u, color, label) {
    const [x0, y0] = P.pt(0, 0);
    const [x1, y1] = P.pt(s, u);
    ctx.strokeStyle = color;
    ctx.fillStyle = color;
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(x0, y0);
    ctx.lineTo(x1, y1);
    ctx.stroke();

    const angle = Math.atan2(y1 - y0, x1 - x0);
    ctx.beginPath();
    ctx.moveTo(x1, y1);
    ctx.lineTo(x1 - 10 * Math.cos(angle - Math.PI / 6), y1 - 10 * Math.sin(angle - Math.PI / 6));
    ctx.lineTo(x1 - 10 * Math.cos(angle + Math.PI / 6), y1 - 10 * Math.sin(angle + Math.PI / 6));
    ctx.fill();

    if (label) {
        ctx.font = "14px Inter";
        ctx.fillText(label, x1 + 6, y1 - 6);
    }
}

function drawDot(ctx, P, s, u, color) {
    const [x, y] = P.pt(s, u);
    ctx.fillStyle = color;
    ctx.beginPath();
    ctx.arc(x, y, 6, 0, Math.PI * 2);
    ctx.fill();
}

function drawArc(ctx, P, etaFrom, etaTo, color, lineWidth) {
    const pts = [];
    const steps = 80;
    for (let i = 0; i <= steps; i++) {
        const e = etaFrom + (etaTo - etaFrom) * (i / steps);
        pts.push(P.pt(Math.sinh(e), Math.cosh(e)));
    }
    strokePath(ctx, pts, color, lineWidth);
}

function drawTick(ctx, P, eta, color, half) {
    // short mark across the hyperbola, along its normal-ish direction
    const [x, y] = P.pt(Math.sinh(eta), Math.cosh(eta));
    const [xa, ya] = P.pt(Math.sinh(eta - 0.02), Math.cosh(eta - 0.02));
    const [xb, yb] = P.pt(Math.sinh(eta + 0.02), Math.cosh(eta + 0.02));
    const dx = xb - xa, dy = yb - ya;
    const len = Math.hypot(dx, dy) || 1;
    const nx = -dy / len, ny = dx / len;
    strokePath(ctx, [[x - nx * half, y - ny * half], [x + nx * half, y + ny * half]], color, 2);
}

function drawCentre(ctx, P, eta, color) {
    const [x, y] = P.pt(Math.sinh(eta), Math.cosh(eta));
    ctx.strokeStyle = color;
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.arc(x, y, 5.5, 0, Math.PI * 2);
    ctx.stroke();
}

// --- Lemniscate overlay -------------------------------------------------
// The projective coordinate of a state is w = tanh(eta) = X/Y, so the point of
// the lemniscate lies on the same ray from the origin as the point of the
// hyperbola. With rho^2 = X^2 + Y^2 = w/(1+w^2) this gives the rational
// parametrization X = w^{3/2}/(1+w^2), Y = w^{1/2}/(1+w^2).

function lemniPoint(eta) {
    // w = tanh(eta) must be non-negative: w = x^2 on the lemniscate.
    // For eta < 0 there is no point of the curve on that ray, so we clamp to
    // the node and report the state as unmappable.
    const w = Math.tanh(eta);
    if (w < 0) return [0, 0, false];
    const d = 1 + w * w;
    return [Math.pow(w, 1.5) / d, Math.sqrt(w) / d, true];
}

function drawLemniscate(ctx, P) {
    const pts = [];
    for (let x = -6; x <= 6; x += 0.02) {
        const d = 1 + x * x * x * x;
        pts.push(P.pt(x * x * x / d, x / d));
    }
    strokePath(ctx, pts, "rgba(148, 163, 184, 0.55)", 1.5);
}

function drawLemniArc(ctx, P, etaFrom, etaTo, color, lineWidth) {
    const pts = [];
    const steps = 90;
    for (let i = 0; i <= steps; i++) {
        const e = etaFrom + (etaTo - etaFrom) * (i / steps);
        const [X, Y] = lemniPoint(e);
        pts.push(P.pt(X, Y));
    }
    strokePath(ctx, pts, color, lineWidth);
}

function drawLemniMark(ctx, P, eta, color, filled) {
    const [X, Y, ok] = lemniPoint(eta);
    const [x, y] = P.pt(X, Y);
    const c = ok ? color : "#ef4444"; // red when the ray misses the curve
    ctx.beginPath();
    ctx.arc(x, y, 4.5, 0, Math.PI * 2);
    if (filled) { ctx.fillStyle = c; ctx.fill(); }
    else { ctx.strokeStyle = c; ctx.lineWidth = 2; ctx.stroke(); }
    if (!ok) {
        // small cross at the node to say the projection is undefined
        ctx.strokeStyle = c;
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.moveTo(x - 7, y - 7); ctx.lineTo(x + 7, y + 7);
        ctx.moveTo(x + 7, y - 7); ctx.lineTo(x - 7, y + 7);
        ctx.stroke();
    }
}

function drawProjectionRay(ctx, P, eta) {
    strokePath(ctx,
        [P.pt(0, 0), P.pt(Math.sinh(eta), Math.cosh(eta))],
        "rgba(16, 185, 129, 0.35)", 1, [3, 4]);
}

// --- The two frame canvases --------------------------------------------

function drawFrame(ctx, canvas, boosted) {
    if (!canvas) return;
    const width = canvas.width;
    const height = canvas.height;
    const P = makePlotter(ctx, width, height);
    const reach = 12;

    ctx.clearRect(0, 0, width, height);
    drawGrid(ctx, P, reach);
    drawNullLines(ctx, P, reach);
    drawHyperbola(ctx, P);

    // Basis vectors. The state transforms as eta -> eta + eta_s, so coordinates
    // go with Lambda(eta_s) and the basis with its inverse Lambda(-eta_s):
    // the primed legs drawn in the rest-frame canvas are (cosh, -sinh) and
    // (-sinh, cosh). In the boosted canvas the roles are swapped.
    const th = boosted ? etaS : -etaS;
    const grey = "#94a3b8", blue = "#3b82f6", red = "#ef4444";

    if (boosted) {
        drawVector(ctx, P, 1, 0, blue, "W'₁");
        drawVector(ctx, P, 0, 1, red, "W'₂");
        drawVector(ctx, P, Math.cosh(th), Math.sinh(th), grey, "W₁");
        drawVector(ctx, P, Math.sinh(th), Math.cosh(th), grey, "W₂");
    } else {
        drawVector(ctx, P, 1, 0, grey, "W₁");
        drawVector(ctx, P, 0, 1, grey, "W₂");
        drawVector(ctx, P, Math.cosh(th), Math.sinh(th), blue, "W'₁");
        drawVector(ctx, P, Math.sinh(th), Math.cosh(th), red, "W'₂");
    }

    // The band swept by the motion, its two turning points, and the fixed
    // centre of oscillation. In the boosted canvas everything is shifted by
    // +eta_s, so the band keeps its width: the boost moves the centre only.
    const shift = boosted ? -etaS : 0;
    const cEta = etaR - shift;
    const amber = "#f59e0b";

    drawArc(ctx, P, cEta - etaAmp, cEta + etaAmp, amber, 5);
    drawTick(ctx, P, cEta - etaAmp, amber, 9);
    drawTick(ctx, P, cEta + etaAmp, amber, 9);
    drawCentre(ctx, P, cEta, amber);

    const eta = etaState - shift;

    // Optional projection onto the lemniscate along the same rays
    if (lemniscate) {
        drawLemniscate(ctx, P);
        if (cEta + etaAmp > 0) {
            const clipped = cEta - etaAmp < 0;
            drawLemniArc(ctx, P, Math.max(cEta - etaAmp, 0), cEta + etaAmp,
                         clipped ? "#ef4444" : amber, 4);
            drawLemniMark(ctx, P, cEta, amber, false);
        }
        drawProjectionRay(ctx, P, eta);
        drawLemniMark(ctx, P, eta, "#10b981", true);
    }

    // The particle
    drawDot(ctx, P, Math.sinh(eta), Math.cosh(eta), "#10b981");
}

// --- Projective coordinate chart ---------------------------------------

function drawProjectiveFrame(ctx) {
    if (!canvasProjective) return;
    const width = canvasProjective.width;
    const height = canvasProjective.height;
    const cx = width / 2;
    const cy = height / 2;

    ctx.clearRect(0, 0, width, height);

    // The projective coordinate is the ratio of the components as drawn, and
    // project() swaps them in the spacelike sector. So w = tanh(eta) there
    // becomes w = coth(eta): the sector is exactly the outside of the unit
    // interval, r^2 > 1, which is what the lemniscate overlay shows as the
    // outer branch. The boost acts the same way on both, since
    // coth(eta + etaS) = (w + V)/(1 + wV) with w = coth(eta) as well.
    const wOf = (eta) => (spacelike ? 1 / Math.tanh(eta) : Math.tanh(eta));

    // Both axes carry w itself, at one and the same scale, so the Mobius map
    // is drawn undistorted. The spacelike sector then runs off the frame
    // whenever coth(eta) is large, which is the honest thing to show: a
    // compressed axis would fit it in at the cost of bending the curve.
    const wMax = 5 / 3;                       // half-range of both axes
    const scale = Math.min(cx, cy) / wMax;
    const V = Math.tanh(etaS);
    const toX = (w) => cx + w * scale;
    const toY = (w) => cy - w * scale;
    const clampX = (px) => Math.max(cx, Math.min(width, px));
    const clampY = (py) => Math.max(0, Math.min(cy, py));

    // The oscillation band. In the spacelike sector an interval of eta
    // straddling zero maps to two pieces, one running out to each end of the
    // line; they are clipped by the frame like everything else.
    const band = (shift) => {
        const lo = etaR - etaAmp + shift, hi = etaR + etaAmp + shift;
        if (!spacelike) return [[Math.tanh(lo), Math.tanh(hi)]];
        if (lo < 0 && hi > 0) return [[-Infinity, 1 / Math.tanh(lo)],
                                      [1 / Math.tanh(hi), Infinity]];
        const a = 1 / Math.tanh(lo), b = 1 / Math.tanh(hi);
        return [[Math.min(a, b), Math.max(a, b)]];
    };

    // The bands are drawn only in the first quadrant, where both coordinates
    // are non-negative and the state has an image on the lemniscate.
    ctx.fillStyle = "rgba(245, 158, 11, 0.14)";
    for (const [a, b] of band(0)) {
        const x0 = clampX(toX(a)), x1 = clampX(toX(b));
        if (x1 > x0) ctx.fillRect(x0, 0, x1 - x0, cy);
    }
    for (const [a, b] of band(etaS)) {
        const y0 = clampY(toY(b)), y1 = clampY(toY(a));
        if (y1 > y0) ctx.fillRect(cx, y0, width - cx, y1 - y0);
    }

    const xC = toX(wOf(etaR));
    const yC = toY(wOf(etaR + etaS));
    ctx.strokeStyle = "rgba(245, 158, 11, 0.55)";
    ctx.lineWidth = 1;
    ctx.setLineDash([4, 4]);
    ctx.beginPath();
    if (xC >= cx) { ctx.moveTo(xC, 0); ctx.lineTo(xC, cy); }
    if (yC <= cy) { ctx.moveTo(cx, yC); ctx.lineTo(width, yC); }
    ctx.stroke();
    ctx.setLineDash([]);

    ctx.strokeStyle = "rgba(255, 255, 255, 0.2)";
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(0, cy); ctx.lineTo(width, cy);
    ctx.moveTo(cx, 0); ctx.lineTo(cx, height);
    ctx.stroke();

    // the unit square: the light cone, and the wall between the two sectors
    ctx.strokeStyle = "rgba(148, 163, 184, 0.30)";
    ctx.setLineDash([2, 4]);
    ctx.strokeRect(toX(-1), toY(1), 2 * scale, 2 * scale);
    ctx.setLineDash([]);

    // the diagonal: the identity at zero boost
    ctx.strokeStyle = "rgba(148, 163, 184, 0.45)";
    ctx.setLineDash([5, 5]);
    ctx.beginPath();
    ctx.moveTo(toX(-wMax), toY(-wMax)); ctx.lineTo(toX(wMax), toY(wMax));
    ctx.stroke();
    ctx.setLineDash([]);

    ctx.fillStyle = "#94a3b8";
    ctx.font = "17px Inter";
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    ctx.fillText("1", toX(1), cy + 6);
    ctx.fillText("\u22121", toX(-1), cy + 6);
    ctx.textAlign = "right";
    ctx.textBaseline = "middle";
    ctx.fillText("1", cx - 6, toY(1));
    ctx.fillText("\u22121", cx - 6, toY(-1));
    ctx.textAlign = "left";
    ctx.textBaseline = "alphabetic";
    ctx.font = "18px Inter";
    const nm = spacelike ? "coth" : "tanh";
    ctx.fillText(nm + "(\u03b7)", width - 96, cy - 10);
    ctx.fillText(nm + "(\u03b7')", cx + 10, 24);

    ctx.strokeStyle = "rgba(239, 68, 68, 0.4)";
    ctx.setLineDash([4, 4]);
    ctx.beginPath();
    if (Math.abs(V) > 0.01) {
        const asymX = toX(-1 / V);
        ctx.moveTo(asymX, 0); ctx.lineTo(asymX, height);
        const asymY = toY(1 / V);
        ctx.moveTo(0, asymY); ctx.lineTo(width, asymY);
    }
    ctx.stroke();
    ctx.setLineDash([]);

    ctx.strokeStyle = "#a855f7";
    ctx.lineWidth = 2;
    ctx.beginPath();
    let first = true, side = 0;
    for (let i = 0; i <= width; i++) {
        const w = (i - cx) / scale;
        const den = 1 + w * V;
        if (Math.sign(den) !== side) { side = Math.sign(den); first = true; }
        if (Math.abs(den) < 1e-6) { first = true; continue; }
        const py = toY((w + V) / den);
        if (first) { ctx.moveTo(i, py); first = false; }
        else ctx.lineTo(i, py);
    }
    ctx.stroke();

    // w = +-1 stay put under every boost: the two fixed points of the map
    ctx.fillStyle = "rgba(226, 232, 240, 0.8)";
    for (const f of [1, -1]) {
        ctx.beginPath();
        ctx.arc(toX(f), toY(f), 3.5, 0, Math.PI * 2);
        ctx.fill();
    }

    const v = wOf(etaState);
    const vp = wOf(etaState + etaS);

    // Either coordinate going negative means the state has no image on the
    // lemniscate in that frame, since w = x^2 >= 0 there.
    const mappable = v >= 0 && vp >= 0;
    ctx.fillStyle = mappable ? "#10b981" : "#ef4444";
    ctx.beginPath();
    ctx.arc(toX(v), toY(vp), 6, 0, Math.PI * 2);
    ctx.fill();

    // Shade the quadrants that cannot be mapped
    ctx.fillStyle = "rgba(239, 68, 68, 0.07)";
    ctx.fillRect(0, 0, cx, height);
    ctx.fillRect(0, cy, width, height - cy);
}

// --- Dragging the centre of oscillation ---------------------------------
// The amber ring on either hyperbola marks the centre eta_r, shifted by the
// boost in the boosted canvas. Dragging it is the same input as the eta_r
// slider. The pointer is carried onto the hyperbola along the ray through the
// origin, the projection the lemniscate overlay already uses, so the ring
// follows the pointer even when it drifts off the curve.

const DRAG_RADIUS = 14;

function canvasPoint(canvas, ev) {
    const rect = canvas.getBoundingClientRect();
    return [
        (ev.clientX - rect.left) * canvas.width / rect.width,
        (ev.clientY - rect.top) * canvas.height / rect.height
    ];
}

function frameShift(boosted) { return boosted ? -etaS : 0; }

// Screen positions of the grabbable rings on one canvas: the centre on the
// hyperbola, and its image on the lemniscate when that overlay is on.
function centreHandles(canvas, boosted) {
    const P = makePlotter(null, canvas.width, canvas.height);
    const cEta = etaR - frameShift(boosted);
    const out = [P.pt(Math.sinh(cEta), Math.cosh(cEta))];
    if (lemniscate && cEta + etaAmp > 0) {
        const [X, Y, ok] = lemniPoint(cEta);
        if (ok) out.push(P.pt(X, Y));
    }
    return out;
}

// Screen point -> rapidity, by the ray through the origin. Undoes project(),
// so it is correct in the spacelike sector too. Returns null for a point that
// no ray of the sector reaches.
function rapidityAt(canvas, x, y) {
    const [cx, cy] = viewOrigin(canvas.width, canvas.height);
    const scale = viewScale(canvas.width, canvas.height);
    const h = (x - cx) / scale, v = (cy - y) / scale;
    const [s, u] = spacelike ? [v, h] : [h, v];
    if (u <= 0) return null;                  // wrong half of the cone
    const w = s / u;
    if (Math.abs(w) >= 1) return null;        // outside the light cone
    return Math.atanh(w);
}

// Single entry point for both the slider and the drag, so they cannot drift
// apart: clamped and quantized to the slider's own range and step.
function setEtaR(value) {
    const lo = etaRSlider ? parseFloat(etaRSlider.min) : -1.5;
    const hi = etaRSlider ? parseFloat(etaRSlider.max) : 1.5;
    const step = etaRSlider ? parseFloat(etaRSlider.step) : 0.01;
    const v = parseFloat(
        (Math.round(Math.min(hi, Math.max(lo, value)) / step) * step).toFixed(6));
    if (v === etaR) return;                   // no reset on sub-step jitter
    etaR = v;
    if (etaRSlider) etaRSlider.value = String(v);
    resetMotion();
    updateUI();
}

function attachCentreDrag(canvas, boosted) {
    if (!canvas) return;
    canvas.style.touchAction = 'none';        // let a touch drag, not scroll
    let dragging = false;

    const onHandle = (x, y) => centreHandles(canvas, boosted)
        .some(([hx, hy]) => Math.hypot(x - hx, y - hy) <= DRAG_RADIUS);

    canvas.addEventListener('pointerdown', (ev) => {
        const [x, y] = canvasPoint(canvas, ev);
        if (!onHandle(x, y)) return;
        dragging = true;
        canvas.setPointerCapture(ev.pointerId);
        canvas.style.cursor = 'grabbing';
        ev.preventDefault();
    });

    canvas.addEventListener('pointermove', (ev) => {
        const [x, y] = canvasPoint(canvas, ev);
        if (!dragging) {
            canvas.style.cursor = onHandle(x, y) ? 'grab' : 'default';
            return;
        }
        const eta = rapidityAt(canvas, x, y);
        if (eta !== null) setEtaR(eta + frameShift(boosted));
        ev.preventDefault();
    });

    const release = (ev) => {
        if (!dragging) return;
        dragging = false;
        if (canvas.hasPointerCapture(ev.pointerId)) canvas.releasePointerCapture(ev.pointerId);
        canvas.style.cursor = 'grab';
    };
    canvas.addEventListener('pointerup', release);
    canvas.addEventListener('pointercancel', release);
}

// --- Rapidity vs time ---------------------------------------------------

function dashedLine(ctx, x0, y, x1, color) {
    ctx.strokeStyle = color;
    ctx.setLineDash([5, 5]);
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(x0, y); ctx.lineTo(x1, y);
    ctx.stroke();
    ctx.setLineDash([]);
}

function drawTimeFrame(ctx) {
    if (!canvasTime) return;
    const width = canvasTime.width;
    const height = canvasTime.height;

    const padL = 46, padR = 14, padT = 18, padB = 26;
    const plotW = width - padL - padR;
    const plotH = height - padT - padB;
    const y0 = padT + plotH / 2;
    const range = 3.0;
    const scaleY = (plotH / 2) / range;

    const toY = (eta) => y0 - eta * scaleY;

    ctx.clearRect(0, 0, width, height);

    ctx.font = "11px Inter";
    ctx.textAlign = "right";
    ctx.textBaseline = "middle";
    for (let v = -3; v <= 3; v += 1) {
        const y = toY(v);
        if (y < padT - 1 || y > padT + plotH + 1) continue;
        ctx.strokeStyle = v === 0 ? "rgba(255,255,255,0.22)" : "rgba(255,255,255,0.07)";
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(padL, y); ctx.lineTo(padL + plotW, y);
        ctx.stroke();
        ctx.fillStyle = "#64748b";
        ctx.fillText(v.toFixed(0), padL - 8, y);
    }

    ctx.strokeStyle = "rgba(255,255,255,0.18)";
    ctx.beginPath();
    ctx.moveTo(padL, padT); ctx.lineTo(padL, padT + plotH);
    ctx.stroke();

    ctx.textAlign = "left";
    ctx.fillStyle = "#94a3b8";
    ctx.font = "12px Inter";
    ctx.fillText("η", 8, padT + 4);
    ctx.textAlign = "right";
    ctx.fillText("t", width - padR, padT + plotH + 18);
    ctx.textAlign = "left";

    dashedLine(ctx, padL, toY(etaR), padL + plotW, "rgba(16,185,129,0.45)");
    dashedLine(ctx, padL, toY(etaR + etaS), padL + plotW, "rgba(168,85,247,0.45)");

    ctx.save();
    ctx.beginPath();
    ctx.rect(padL, padT, plotW, plotH);
    ctx.clip();

    const drawCurve = (shift, color) => {
        ctx.strokeStyle = color;
        ctx.lineWidth = 2;
        ctx.beginPath();
        history.forEach((p, i) => {
            const x = padL + ((p.t - (time - WINDOW_T)) / WINDOW_T) * plotW;
            const y = toY(p.eta - shift);
            if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
        });
        ctx.stroke();

        ctx.fillStyle = color;
        ctx.beginPath();
        ctx.arc(padL + plotW, toY(etaState - shift), 5, 0, Math.PI * 2);
        ctx.fill();
    };

    drawCurve(0, "#10b981");
    drawCurve(-etaS, "#a855f7");

    ctx.restore();

    ctx.font = "12px Inter";
    ctx.textBaseline = "alphabetic";
    ctx.fillStyle = "#10b981";
    ctx.fillText(`η = ${etaState.toFixed(2)}`, padL + 8, padT + 14);
    ctx.fillStyle = "#a855f7";
    ctx.fillText(`η' = η + ηₛ = ${(etaState + etaS).toFixed(2)}`, padL + 90, padT + 14);
}

function animate() {
    const dt = 0.02;
    time += dt;
    stepMotion(dt);
    history.push({ t: time, eta: etaState });
    while (history.length > 1 && history[0].t < time - WINDOW_T) history.shift();

    if (ctxOrig) drawFrame(ctxOrig, canvasOriginal, false);
    if (ctxBoost) drawFrame(ctxBoost, canvasBoosted, true);
    if (ctxProj) drawProjectiveFrame(ctxProj);
    if (ctxTime) drawTimeFrame(ctxTime);

    requestAnimationFrame(animate);
}

// Init
document.addEventListener('DOMContentLoaded', () => {
    etaSSlider = document.getElementById('etas-slider');
    etaSValDisplay = document.getElementById('etas-val');
    lcToggle = document.getElementById('spacelike-toggle');
    canvasOriginal = document.getElementById('canvas-original');
    canvasBoosted = document.getElementById('canvas-boosted');
    canvasProjective = document.getElementById('canvas-projective');
    canvasTime = document.getElementById('canvas-time');
    mathMatrix = document.getElementById('math-matrix');

    ctxOrig = canvasOriginal ? canvasOriginal.getContext('2d') : null;
    ctxBoost = canvasBoosted ? canvasBoosted.getContext('2d') : null;
    ctxProj = canvasProjective ? canvasProjective.getContext('2d') : null;
    ctxTime = canvasTime ? canvasTime.getContext('2d') : null;

    etaS = parseFloat(etaSSlider.value);

    etaRSlider = document.getElementById('etar-slider');
    etaRValDisplay = document.getElementById('etar-val');
    ampSlider = document.getElementById('amp-slider');
    ampValDisplay = document.getElementById('amp-val');

    mmSlider = document.getElementById('mm-slider');
    mmValDisplay = document.getElementById('mm-val');

    if (etaRSlider) etaR = parseFloat(etaRSlider.value);
    if (ampSlider) etaAmp = parseFloat(ampSlider.value);
    if (mmSlider) MM = parseFloat(mmSlider.value);

    etaSSlider.addEventListener('input', (e) => {
        etaS = parseFloat(e.target.value);
        updateUI();
    });

    if (etaRSlider) {
        etaRSlider.addEventListener('input', (e) => {
            setEtaR(parseFloat(e.target.value));
        });
    }

    if (ampSlider) {
        ampSlider.addEventListener('input', (e) => {
            etaAmp = parseFloat(e.target.value);
            resetMotion();
            updateUI();
        });
    }

    if (mmSlider) {
        mmSlider.addEventListener('input', (e) => {
            MM = parseFloat(e.target.value);
            resetMotion();
            updateUI();
        });
    }

    lemniToggle = document.getElementById('lemniscate-toggle');
    if (lemniToggle) {
        lemniscate = lemniToggle.checked;
        lemniToggle.addEventListener('change', (e) => {
            lemniscate = e.target.checked;
            updateUI();
        });
    }

    if (lcToggle) {
        spacelike = lcToggle.checked;
        lcToggle.addEventListener('change', (e) => {
            spacelike = e.target.checked;
            updateUI();
        });
    }

    attachCentreDrag(canvasOriginal, false);
    attachCentreDrag(canvasBoosted, true);

    resetMotion();
    updateUI();
    animate();
});
