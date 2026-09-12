// Two independent solutions of the same one-dimensional lemniscatic oscillator.
//
//   Route A  integrates the Duffing equation for phi and converts to x^2.
//   Route B  integrates the Euler-Lagrange equation for x directly.
//
// They share only alpha_1D^2 and E_1D. No quantity computed on one route is
// ever used on the other.

let outer = false;   // false: the branch r < 1, true: the branch r > 1
let a2 = 1.0;   // alpha_1D^2
let E1 = 0.15;  // E_1D

const XMAX = 60;          // stop route B once |x| passes this
const PERIODS = 3;        // nominal window length, in periods of phi
const STEPS = 8000;       // integration steps over the window
const PLAY_SECONDS = 7;   // real time taken to replay the window
let U_RANGE = [-0.05, 0.85];   // fixed vertical range of the U(x) panel
let V_RANGE = [-1.2, 1.4];     // fixed vertical range of the V(phi) panel
let PHI_MAX = 1.5;             // fixed horizontal half-range of the V(phi) panel

let current = null;       // last integration
let clock0 = null;

// ---------------------------------------------------------------- utilities

function ellipticK(k2) {
    // AGM; k2 is the parameter m = k^2, required in [0,1)
    let a = 1, b = Math.sqrt(Math.max(1 - k2, 0));
    for (let i = 0; i < 40 && Math.abs(a - b) > 1e-15; i++) {
        const an = 0.5 * (a + b);
        b = Math.sqrt(a * b);
        a = an;
    }
    return Math.PI / (2 * a);
}

function derived() {
    const Fp2 = a2 / 4 + E1;
    const Fm2 = a2 / 4 - E1;
    const Fp = Math.sqrt(Fp2);
    const M2 = Fp2 * Fm2;
    const k2 = 1 - Fm2 / Fp2;
    // period of phi: dn has period 2K for k<1, cn has 4K(1/k) for k>1
    let T;
    if (k2 < 1) {
        T = ellipticK(k2) / Fp;
    } else {
        const k = Math.sqrt(k2);
        T = 2 * ellipticK(1 / k2) / (Fp * k);
    }
    return { Fp2, Fm2, Fp, M2, k2, T };
}

function rk4(f, y, dt) {
    const k1 = f(y);
    const k2 = f([y[0] + dt / 2 * k1[0], y[1] + dt / 2 * k1[1]]);
    const k3 = f([y[0] + dt / 2 * k2[0], y[1] + dt / 2 * k2[1]]);
    const k4 = f([y[0] + dt * k3[0], y[1] + dt * k3[1]]);
    return [
        y[0] + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]),
        y[1] + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
    ];
}

// ------------------------------------------------------------ the two routes

function integrate() {
    const { Fp, T, M2 } = derived();

    // Below the barrier phi is dnoidal and stays positive, so x^2 touches zero
    // once per period of phi and x itself changes sign there: the period of x
    // is 2T, not T. Over the barrier phi is cnoidal and the two periods agree.
    // The window must be a whole number of periods of x, otherwise the replay
    // wraps onto a state with the opposite velocity and the point appears to
    // bounce off the node instead of passing to the other petal.
    const Tx = M2 > 0 ? 2 * T : T;                          // period of x
    const nWin = Math.max(1, Math.round(PERIODS * T / Tx)); // whole periods of x
    const tEnd = nWin * Tx;
    const phiPeriods = tEnd / T;                            // still an integer
    const dt = tEnd / STEPS;

    // Route A: phi'' = -8 phi^3 + 2 a2 phi,  phi(0) = F_+, phi'(0) = 0
    const fA = (y) => [y[1], -8 * y[0] ** 3 + 2 * a2 * y[0]];
    // Below the barrier U(x) has two allowed regions, exchanged by x -> 1/x:
    // the inner one x^2 <= x_-^2 and the outer one x^2 >= 1/x_-^2. In the
    // Duffing variable they are the two wells, phi in [F_-, F_+] and its
    // mirror [-F_+, -F_-], so the branch is selected by the sign of phi(0).
    let yA = [outer ? -Fp : Fp, 0];

    // Route B: Euler-Lagrange for L = m(s)(sdot^2 - a2 s^2)/2, m = 1/(1+s^4).
    // The inversion s -> 1/s is an isometry of this Lagrangian, so the same
    // equation holds in the chart y = 1/x. Swapping charts whenever |s| > 1
    // carries the motion through x = infinity, where the orbit passes from one
    // petal of the lemniscate to the other.
    const fB = (y) => {
        const x = y[0], v = y[1];
        const r = -4 * x ** 3 / (1 + x ** 4);   // m'/m
        return [v, -0.5 * r * v * v - 0.5 * a2 * (r * x * x + 2 * x)];
    };
    // The same initial numbers describe both branches: read in the x-chart they
    // start the motion at the origin, read in the chart y = 1/x they start it
    // at x = infinity, which is the outer branch.
    let yB = [0, Math.sqrt(2 * E1)];
    let chart = outer ? 1 : 0;                  // 0: x-chart, 1: y = 1/x chart

    const ts = [], xA = [], xB = [], phi = [], xraw = [], chartOf = [];
    let escaped = null;
    const yB0 = [yB[0], yB[1]], chart0 = chart;
    let yBend = null, chartEnd = null;

    for (let i = 0; i <= STEPS; i++) {
        const t = i * dt;
        const denom = Fp + yA[0];
        const x2A = Math.abs(denom) < 1e-12 ? Infinity : (Fp - yA[0]) / denom;
        const x2B = chart === 0 ? yB[0] * yB[0] : 1 / (yB[0] * yB[0]);

        ts.push(t);
        xA.push(x2A);
        xB.push(x2B);
        phi.push(yA[0]);
        xraw.push(yB[0]);
        chartOf.push(chart);

        if (escaped === null && x2B > XMAX * XMAX) escaped = t;  // x^2 diverges here
        if (i === STEPS) { yBend = [yB[0], yB[1]]; chartEnd = chart; }

        yA = rk4(fA, yA, dt);
        yB = rk4(fB, yB, dt);
        if (Math.abs(yB[0]) > 1) {              // swap charts; |1/s| < 1, no chatter
            yB = [1 / yB[0], -yB[1] / (yB[0] * yB[0])];
            chart = 1 - chart;
        }
    }

    // How far route B is from returning to its own initial state at t = tEnd.
    // A window cut on a whole period of x closes to integrator accuracy; a
    // window cut on half a period closes in position and fails in velocity,
    // which is exactly the bounce-off-the-node artefact. Comparable only when
    // both ends are read in the same chart.
    let closure = NaN;
    if (yBend && chartEnd === chart0 && Number.isFinite(yBend[0])) {
        closure = Math.hypot(yBend[0] - yB0[0], yBend[1] - yB0[1]);
    }

    return { ts, xA, xB, phi, xraw, chartOf, tEnd, escaped, dt,
             phiPeriods, nWin, closure };
}

// ------------------------------------------------------------------ drawing

const PAD = { l: 54, r: 16, t: 14, b: 28 };

function axes(ctx, W, H, xLabel, yTicks, fmt) {
    ctx.clearRect(0, 0, W, H);
    const pw = W - PAD.l - PAD.r, ph = H - PAD.t - PAD.b;

    ctx.font = '11px Inter';
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    for (const { v, y } of yTicks) {
        ctx.strokeStyle = 'rgba(255,255,255,0.07)';
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(PAD.l, y); ctx.lineTo(PAD.l + pw, y);
        ctx.stroke();
        ctx.fillStyle = '#64748b';
        ctx.fillText(fmt(v), PAD.l - 8, y);
    }

    ctx.strokeStyle = 'rgba(255,255,255,0.18)';
    ctx.beginPath();
    ctx.moveTo(PAD.l, PAD.t); ctx.lineTo(PAD.l, PAD.t + ph);
    ctx.lineTo(PAD.l + pw, PAD.t + ph);
    ctx.stroke();

    ctx.textAlign = 'right';
    ctx.textBaseline = 'alphabetic';
    ctx.fillStyle = '#94a3b8';
    ctx.fillText(xLabel, W - PAD.r, H - 8);
}

function drawX2(data) {
    const cv = document.getElementById('canvas-x2');
    const ctx = cv.getContext('2d');
    const W = cv.width, H = cv.height;
    const pw = W - PAD.l - PAD.r, ph = H - PAD.t - PAD.b;

    const finite = data.xA.filter(Number.isFinite)
        .concat(data.xB.filter(Number.isFinite));
    let yMax = Math.min(Math.max(...finite) * 1.08, 12);
    if (!(yMax > 0)) yMax = 1;

    const toX = (t) => PAD.l + (t / data.tEnd) * pw;
    const toY = (v) => PAD.t + ph - (v / yMax) * ph;

    const ticks = [];
    const nT = 5;
    for (let i = 0; i <= nT; i++) {
        const v = (yMax * i) / nT;
        ticks.push({ v, y: toY(v) });
    }
    axes(ctx, W, H, 't', ticks, (v) => v.toFixed(2));

    const curve = (arr, color, dash) => {
        ctx.strokeStyle = color;
        ctx.lineWidth = 2;
        ctx.setLineDash(dash);
        ctx.beginPath();
        let pen = false;
        for (let i = 0; i < arr.length; i++) {
            const v = arr[i];
            if (!Number.isFinite(v) || v > yMax) { pen = false; continue; }
            const x = toX(data.ts[i]), y = toY(v);
            if (!pen) { ctx.moveTo(x, y); pen = true; } else ctx.lineTo(x, y);
        }
        ctx.stroke();
        ctx.setLineDash([]);
    };

    curve(data.xA, '#10b981', []);
    curve(data.xB, '#f59e0b', [6, 4]);
}

// The lift, drawn honestly: phi runs along the horizontal axis and x along the
// vertical one, so the two branches x = +-sqrt((F+ - phi)/(F+ + phi)) are the
// upper and lower halves of one arch. The white dot is the state and its two
// shadows are the dots on the axes: x in amber, phi in green.
//
// The vertical axis carries x through v = (2/pi) arctan x, all but linear near
// the origin and putting x = infinity on the two edges. Without that the outer
// branch, which runs out to infinity twice per period, spends more than half
// its time off the frame. In the angle psi = 2 arctan x the arch is simply
// phi = F+ cos psi, so it is sampled in psi and never touches a square root.
function drawXPhi(data, i) {
    const cv = document.getElementById('canvas-xphi');
    if (!cv) return;
    const ctx = cv.getContext('2d');
    const W = cv.width, H = cv.height;
    const pw = W - PAD.l - PAD.r, ph = H - PAD.t - PAD.b;
    const { Fp, M2 } = derived();

    const xIn = M2 > 0 ? Math.sqrt((a2 / 4 - Math.sqrt(M2)) / E1) : null;
    const pLo = -Fp * 1.06, pHi = Fp * 1.06;
    const toX = (p) => PAD.l + ((p - pLo) / (pHi - pLo)) * pw;
    const toY = (x) => PAD.t + ph / 2 - (2 / Math.PI) * Math.atan(x) * (ph / 2);
    const swept = (x) => (M2 <= 0 ? true
                        : outer ? Math.abs(x) >= 1 / xIn : Math.abs(x) <= xIn);

    ctx.clearRect(0, 0, W, H);
    ctx.font = '11px Inter';

    ctx.strokeStyle = 'rgba(255,255,255,0.10)';
    ctx.lineWidth = 1;
    ctx.strokeRect(PAD.l, PAD.t, pw, ph);
    ctx.strokeStyle = 'rgba(255,255,255,0.22)';
    ctx.beginPath();
    ctx.moveTo(PAD.l, toY(0)); ctx.lineTo(PAD.l + pw, toY(0));
    ctx.moveTo(toX(0), PAD.t); ctx.lineTo(toX(0), PAD.t + ph);
    ctx.stroke();

    ctx.fillStyle = '#64748b';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    ctx.fillText('\u2212F\u208a', toX(-Fp), toY(0) + 6);
    ctx.fillText('F\u208a', toX(Fp), toY(0) + 6);
    ctx.textAlign = 'left';
    ctx.textBaseline = 'middle';
    for (const v of [-5, -2, -1, 1, 2, 5]) ctx.fillText(String(v), toX(0) + 5, toY(v));
    ctx.fillText('\u221e', toX(0) + 5, PAD.t + 8);
    ctx.fillText('\u2212\u221e', toX(0) + 5, PAD.t + ph - 8);
    ctx.fillStyle = '#94a3b8';
    ctx.textAlign = 'right';
    ctx.textBaseline = 'alphabetic';
    ctx.fillText('\u03c6', W - PAD.r, toY(0) - 6);
    ctx.textAlign = 'left';
    ctx.fillText('(2/\u03c0)arctan x', toX(0) + 6, PAD.t + 24);
    ctx.font = '11px Inter';
    ctx.fillStyle = '#64748b';
    ctx.textAlign = 'right';
    ctx.fillText('vertical axis: (2/\u03c0) arctan x', W - PAD.r - 4, PAD.t + 14);

    // the arch, brighter on the part the motion actually covers
    const arch = (test, colour, width) => {
        ctx.strokeStyle = colour;
        ctx.lineWidth = width;
        ctx.beginPath();
        let pen = false;
        for (let j = 0; j <= 700; j++) {
            const psi = -Math.PI + 2 * Math.PI * j / 700;
            const x = Math.tan(psi / 2);
            if (!test(x)) { pen = false; continue; }
            const px = toX(Fp * Math.cos(psi));
            const py = PAD.t + ph / 2 - (psi / Math.PI) * (ph / 2);
            if (!pen) { ctx.moveTo(px, py); pen = true; } else ctx.lineTo(px, py);
        }
        ctx.stroke();
    };
    arch(() => true, 'rgba(148,163,184,0.45)', 1.8);
    arch(swept, '#e2e8f0', 2.6);

    // the state and its two shadows
    const s = data.xraw[i];
    const trueX = data.chartOf[i] === 0 ? s : 1 / s;
    const phi = data.phi[i];
    const py = toY(trueX);

    ctx.setLineDash([3, 4]);
    ctx.strokeStyle = 'rgba(226,232,240,0.30)';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(toX(phi), py); ctx.lineTo(toX(phi), toY(0));
    ctx.moveTo(toX(phi), py); ctx.lineTo(toX(0), py);
    ctx.stroke();
    ctx.setLineDash([]);

    ctx.fillStyle = '#10b981';
    ctx.beginPath(); ctx.arc(toX(phi), toY(0), 5, 0, Math.PI * 2); ctx.fill();
    ctx.fillStyle = '#f59e0b';
    ctx.beginPath(); ctx.arc(toX(0), py, 5, 0, Math.PI * 2); ctx.fill();
    ctx.fillStyle = '#e2e8f0';
    ctx.beginPath(); ctx.arc(toX(phi), py, 5.5, 0, Math.PI * 2); ctx.fill();
}

// --- Dragging the energy level ------------------------------------------
// The dashed violet line is E_1D itself on the U(x) panel and the Duffing
// level -2M^2 on the V(phi) one, where M^2 = (a^2/4)^2 - E^2, so that level is
// 2E^2 - a^4/8 and inverts to E = sqrt((level + a^4/8)/2). Dragging either line
// is the same input as the E_1D slider, and both go through setE1.

function setE1(value) {
    const es = document.getElementById('e-slider');
    const lo = es ? parseFloat(es.min) : 0.005;
    const hi = es ? parseFloat(es.max) : 0.8;
    const step = es ? parseFloat(es.step) : 0.005;
    if (!Number.isFinite(value)) return;
    const v = parseFloat(
        (Math.round(Math.min(hi, Math.max(lo, value)) / step) * step).toFixed(6));
    if (v === E1) return;                  // no re-integration on sub-step jitter
    E1 = v;
    if (es) es.value = String(v);
    update();
}

const LEVEL_TO_E = {
    'canvas-ulem': (level) => level,
    'canvas-uduf': (level) => Math.sqrt(Math.max(0, (level + a2 * a2 / 8) / 2))
};

function attachLevelDrag(id) {
    const cv = document.getElementById(id);
    if (!cv) return;
    cv.style.touchAction = 'none';
    let dragging = false;

    const yAt = (ev) => {
        const r = cv.getBoundingClientRect();
        return (ev.clientY - r.top) * cv.height / r.height;
    };
    const onLine = (py) => cv._level && Math.abs(py - cv._level.y) <= 9;

    cv.addEventListener('pointerdown', (ev) => {
        if (!onLine(yAt(ev))) return;
        dragging = true;
        cv.setPointerCapture(ev.pointerId);
        cv.style.cursor = 'grabbing';
        ev.preventDefault();
    });

    cv.addEventListener('pointermove', (ev) => {
        const py = yAt(ev);
        if (!dragging) {
            cv.style.cursor = onLine(py) ? 'ns-resize' : 'default';
            return;
        }
        setE1(LEVEL_TO_E[id](cv._level.toValue(py)));
        ev.preventDefault();
    });

    const release = (ev) => {
        if (!dragging) return;
        dragging = false;
        if (cv.hasPointerCapture(ev.pointerId)) cv.releasePointerCapture(ev.pointerId);
        cv.style.cursor = 'ns-resize';
    };
    cv.addEventListener('pointerup', release);
    cv.addEventListener('pointercancel', release);
}

function drawErr(data) {
    const cv = document.getElementById('canvas-err');
    const ctx = cv.getContext('2d');
    const W = cv.width, H = cv.height;
    const pw = W - PAD.l - PAD.r, ph = H - PAD.t - PAD.b;

    const lo = -16, hi = 0;                     // log10 range
    const toX = (t) => PAD.l + (t / data.tEnd) * pw;
    const toY = (L) => PAD.t + ph - ((L - lo) / (hi - lo)) * ph;

    const ticks = [];
    for (let L = lo; L <= hi; L += 4) ticks.push({ v: L, y: toY(L) });
    axes(ctx, W, H, 't', ticks, (L) => '1e' + L);

    ctx.strokeStyle = '#a855f7';
    ctx.lineWidth = 2;
    ctx.beginPath();
    let pen = false;
    for (let i = 0; i < data.ts.length; i++) {
        const a = data.xA[i], b = data.xB[i];
        if (!Number.isFinite(a) || !Number.isFinite(b) || a > 10 || b > 10) { pen = false; continue; }
        const d = Math.abs(a - b);
        const L = Math.max(Math.log10(d + 1e-18), lo);
        const x = toX(data.ts[i]), y = toY(L);
        if (!pen) { ctx.moveTo(x, y); pen = true; } else ctx.lineTo(x, y);
    }
    ctx.stroke();
}

// ------------------------------------------------------------ potential wells

function wellPanel(id, xs, ys, level, dot, xLabel, colour, fixedY, cover) {
    const cv = document.getElementById(id);
    const ctx = cv.getContext('2d');
    const W = cv.width, H = cv.height;
    const P = { l: 46, r: 14, t: 14, b: 26 };
    const pw = W - P.l - P.r, ph = H - P.t - P.b;

    const xLo = xs[0], xHi = xs[xs.length - 1];
    let yLo, yHi;
    if (fixedY) {
        [yLo, yHi] = fixedY;
    } else {
        yLo = Math.min(...ys, level);
        yHi = Math.max(...ys, level);
        const pad = 0.12 * (yHi - yLo || 1);
        yLo -= pad; yHi += pad;
    }

    const toX = (x) => P.l + ((x - xLo) / (xHi - xLo)) * pw;
    const toY = (y) => P.t + ph - ((y - yLo) / (yHi - yLo)) * ph;

    // Where the level sits on screen, and the inverse of toY, so that the
    // pointer handler can drag the line without knowing the panel's scale.
    cv._level = {
        y: toY(level),
        toValue: (py) => yLo + (1 - (py - P.t) / ph) * (yHi - yLo)
    };

    ctx.clearRect(0, 0, W, H);

    // zero line and frame
    ctx.font = '11px Inter';
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    for (let i = 0; i <= 4; i++) {
        const v = yLo + (yHi - yLo) * i / 4;
        const y = toY(v);
        ctx.strokeStyle = 'rgba(255,255,255,0.06)';
        ctx.lineWidth = 1;
        ctx.beginPath(); ctx.moveTo(P.l, y); ctx.lineTo(P.l + pw, y); ctx.stroke();
        ctx.fillStyle = '#64748b';
        ctx.fillText(v.toFixed(2), P.l - 7, y);
    }
    ctx.strokeStyle = 'rgba(255,255,255,0.18)';
    ctx.beginPath();
    ctx.moveTo(P.l, P.t); ctx.lineTo(P.l, P.t + ph); ctx.lineTo(P.l + pw, P.t + ph);
    ctx.stroke();

    // the classically allowed band, U <= level
    ctx.fillStyle = 'rgba(255,255,255,0.05)';
    let run = null;
    for (let i = 0; i < xs.length; i++) {
        const inside = ys[i] <= level;
        if (inside && run === null) run = toX(xs[i]);
        if ((!inside || i === xs.length - 1) && run !== null) {
            ctx.fillRect(run, P.t, toX(xs[i]) - run, ph);
            run = null;
        }
    }

    ctx.save();
    ctx.beginPath();
    ctx.rect(P.l, P.t, pw, ph);
    ctx.clip();

    // potential curve
    ctx.strokeStyle = 'rgba(226,232,240,0.75)';
    ctx.lineWidth = 2;
    ctx.beginPath();
    xs.forEach((x, i) => (i ? ctx.lineTo(toX(x), toY(ys[i])) : ctx.moveTo(toX(x), toY(ys[i]))));
    ctx.stroke();

    // energy level
    ctx.strokeStyle = 'rgba(168,85,247,0.85)';
    ctx.setLineDash([5, 4]);
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.moveTo(P.l, toY(level)); ctx.lineTo(P.l + pw, toY(level));
    ctx.stroke();
    ctx.setLineDash([]);

    // The double cover of the swept interval. phi is even in x, so an interior
    // value of phi has the two preimages +-x, and the two sheets meet exactly
    // where phi = +-F+, that is at x = 0 and at x = infinity. In the angle
    // psi = 2 arctan x the cover is the circle phi = F+ cos psi, drawn here as
    // an ellipse of height 2*OFF pixels: a closed loop over the barrier, an
    // open hairpin below it, where one end is a genuine turning point and the
    // sheets stay apart. Which sheet the dot rides on is the sign of x, the one
    // bit the Duffing equation does not carry.
    const OFF = 7;
    if (cover) {
        const yMid = toY(level);
        const pt = (psi) => [toX(cover.Fp * Math.cos(psi)), yMid - OFF * Math.sin(psi)];
        ctx.strokeStyle = 'rgba(148,163,184,0.85)';
        ctx.lineWidth = 2;
        ctx.beginPath();
        const n = 240;
        for (let j = 0; j <= n; j++) {
            const psi = cover.psi0 + cover.half * (2 * j / n - 1);
            const [px, py] = pt(psi);
            j ? ctx.lineTo(px, py) : ctx.moveTo(px, py);
        }
        ctx.stroke();

        // open ends: a turning point of the motion, where the sheets do not meet
        if (cover.half < Math.PI - 1e-9) {
            ctx.fillStyle = 'rgba(148,163,184,0.9)';
            for (const s of [-1, 1]) {
                const [px, py] = pt(cover.psi0 + s * cover.half);
                ctx.beginPath(); ctx.arc(px, py, 2.6, 0, Math.PI * 2); ctx.fill();
            }
        }

        // sheet labels, only where the two sheets are far enough apart to read
        const wide = cover.psi0 + Math.max(-cover.half,
                     Math.min(cover.half, Math.PI / 2));
        if (Math.abs(Math.sin(wide)) > 0.55) {
            ctx.font = '10px Inter';
            ctx.textAlign = 'center';
            ctx.fillStyle = '#64748b';
            const [lx] = pt(wide);
            ctx.textBaseline = 'alphabetic';
            ctx.fillText('x > 0', lx, yMid - OFF - 5);
            ctx.textBaseline = 'top';
            ctx.fillText('x < 0', lx, yMid + OFF + 4);
        }
    }

    // the particle, riding on the level
    if (Number.isFinite(dot) && dot >= xLo && dot <= xHi) {
        let dy = 0;
        if (cover && cover.sign) {
            const q = Math.min(1, Math.abs(dot) / cover.Fp);   // |sin psi|
            dy = -OFF * cover.sign * Math.sqrt(1 - q * q);
        }
        ctx.fillStyle = colour;
        ctx.beginPath();
        ctx.arc(toX(dot), toY(level) + dy, 5.5, 0, Math.PI * 2);
        ctx.fill();
    }

    ctx.restore();

    ctx.textAlign = 'right';
    ctx.textBaseline = 'alphabetic';
    ctx.fillStyle = '#94a3b8';
    ctx.fillText(xLabel, W - P.r, H - 7);
}

function drawLemniscate(data, i) {
    const cv = document.getElementById('canvas-lemni');
    if (!cv) return;
    const ctx = cv.getContext('2d');
    const W = cv.width, H = cv.height;
    const R = 0.66;                                   // the curve fits in |X|,|Y| < 0.58
    const sc = Math.min(W, H) / (2 * R) * 0.92;
    const cx = W / 2, cy = H / 2;

    const pt = (x) => {
        const d = 1 + x ** 4;
        return [cx + (x ** 3 / d) * sc, cy - (x / d) * sc];
    };

    ctx.clearRect(0, 0, W, H);

    // axes and the symmetry axis X = Y
    ctx.strokeStyle = 'rgba(255,255,255,0.10)';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(cx - R * sc, cy); ctx.lineTo(cx + R * sc, cy);
    ctx.moveTo(cx, cy - R * sc); ctx.lineTo(cx, cy + R * sc);
    ctx.stroke();
    ctx.strokeStyle = 'rgba(234,179,8,0.28)';
    ctx.setLineDash([4, 4]);
    ctx.beginPath();
    ctx.moveTo(cx - R * sc, cy + R * sc); ctx.lineTo(cx + R * sc, cy - R * sc);
    ctx.stroke();
    ctx.setLineDash([]);

    // An arc of the curve, swept by the parameter of the given chart. Sampled
    // uniformly in u = atan(t), not in t: the map t -> point crowds the whole
    // outside of a petal into large |t|, so a uniform grid in t spends its
    // points near the node and leaves the petal itself as a handful of long
    // chords. In u the points spread evenly, and u reaches the node from the
    // other side as t -> +-infinity, which closes the loops.
    const arc = (from, to, colour, width, ch = 0) => {
        const n = 600;
        const uFrom = Math.atan(from), uTo = Math.atan(to);
        ctx.strokeStyle = colour;
        ctx.lineWidth = width;
        ctx.beginPath();
        for (let j = 0; j <= n; j++) {
            const t = Math.tan(uFrom + (uTo - uFrom) * j / n);
            const d = 1 + t ** 4;
            const px = ch === 0 ? cx + (t ** 3 / d) * sc : cx + (t / d) * sc;
            const py = ch === 0 ? cy - (t / d) * sc : cy - (t ** 3 / d) * sc;
            j ? ctx.lineTo(px, py) : ctx.moveTo(px, py);
        }
        ctx.stroke();
    };

    // whole curve, then the arc actually swept
    arc(-1e4, 1e4, 'rgba(148,163,184,0.5)', 1.5);

    const { M2 } = derived();
    const bound = M2 > 0;
    const xMaxLabel = bound ? Math.sqrt((a2 / 4 - Math.sqrt(M2)) / E1) : NaN;
    if (bound) {
        const xIn = Math.sqrt((a2 / 4 - Math.sqrt(M2)) / E1);   // inner turning point
        // The outer branch is the image of the inner one under x -> 1/x, so in
        // the chart y = 1/x it is the same interval.
        arc(-xIn, xIn, '#f59e0b', 3.5, outer ? 1 : 0);
    } else {
        arc(-1e4, 1e4, '#f59e0b', 3.5, 0);   // the whole figure eight
    }

    // tips of the petals: x = +-1, the barrier r = 1
    ctx.fillStyle = 'rgba(234,179,8,0.9)';
    for (const x of [-1, 1]) {
        const [px, py] = pt(x);
        ctx.beginPath(); ctx.arc(px, py, 3, 0, Math.PI * 2); ctx.fill();
    }

    // The particle, drawn from whichever chart is active. In the chart
    // y = 1/x the same point of the curve is (X,Y) = (y, y^3)/(1+y^4), the
    // reflection of the x-chart formula, so the passage through x = infinity
    // is continuous and the orbit moves onto the opposite petal.
    const pointAt = (idx) => {
        const s = data.xraw[idx];
        if (!Number.isFinite(s)) return null;
        const d = 1 + s ** 4;
        const px = data.chartOf[idx] === 0 ? cx + (s ** 3 / d) * sc : cx + (s / d) * sc;
        const py = data.chartOf[idx] === 0 ? cy - (s / d) * sc : cy - (s ** 3 / d) * sc;
        return Number.isFinite(px) && Number.isFinite(py) ? [px, py] : null;
    };

    // A short trail behind the particle. Both petals meet the node along the
    // same vertical tangent, so a single dot crossing x = 0 at full speed is
    // indistinguishable by eye from a dot bouncing off the node; the trail
    // shows which of the two actually happened.
    const span = data.xraw.length - 1;
    const TRAIL = Math.max(8, Math.round(span * 0.04));
    const STRIDE = Math.max(1, Math.round(TRAIL / 80));
    const periodic = Number.isFinite(data.closure) && data.closure < 1e-6;
    ctx.lineWidth = 3;
    ctx.lineCap = 'round';
    for (let k = TRAIL; k > 0; k -= STRIDE) {
        let i0 = i - k, i1 = Math.min(i - k + STRIDE, i);
        if (i0 < 0) {
            if (!periodic) continue;            // window does not close: no wrap
            i0 += span; i1 += span;
        }
        const p0 = pointAt(i0), p1 = pointAt(i1);
        if (!p0 || !p1) continue;
        ctx.strokeStyle = `rgba(16,185,129,${(0.5 * (1 - k / TRAIL)).toFixed(3)})`;
        ctx.beginPath();
        ctx.moveTo(p0[0], p0[1]); ctx.lineTo(p1[0], p1[1]);
        ctx.stroke();
    }

    const here = pointAt(i);
    if (here) {
        ctx.fillStyle = '#10b981';
        ctx.beginPath(); ctx.arc(here[0], here[1], 5.5, 0, Math.PI * 2); ctx.fill();
    }

    ctx.font = '11px Inter';
    ctx.fillStyle = '#64748b';
    ctx.textAlign = 'left';
    ctx.fillText(bound
        ? (outer ? 'x\u00b2 \u2265 ' + (1 / (xMaxLabel * xMaxLabel)).toFixed(3)
                 : 'x\u00b2 \u2264 ' + (xMaxLabel * xMaxLabel).toFixed(3))
        : 'both petals', 10, H - 10);
}

function drawWells(data, i) {
    const { Fp, M2 } = derived();

    drawXPhi(data, i);   // the lift panel shares the frame clock with the wells

    // lemniscatic oscillator: barrier of height a2/4 at x = +-1
    const xs = [], us = [];
    for (let x = -5; x <= 5; x += 0.01) {
        xs.push(x);
        us.push(0.5 * a2 * x * x / (1 + x ** 4));
    }
    // The dot must sit at the true x, not at the coordinate of whichever chart
    // is active: in the chart y = 1/x a large x appears as a small y, and the
    // point would seem to oscillate about the centre while it is in fact past
    // the barrier. Over the barrier it leaves the frame on one side and
    // re-enters on the other, having passed through x = infinity.
    const s = data.xraw[i];
    const trueX = data.chartOf[i] === 0 ? s : 1 / s;

    // Fixed vertical range, so that raising alpha_1D raises the barrier
    // against a stationary frame instead of rescaling the picture.
    wellPanel('canvas-ulem', xs, us, E1, trueX, 'x', '#f59e0b', U_RANGE);

    // Duffing: level -2 M^2
    const ps = [], vs = [];
    for (let p = -PHI_MAX; p <= PHI_MAX; p += PHI_MAX / 300) {
        ps.push(p);
        vs.push(2 * p ** 4 - a2 * p * p);
    }
    // Cover of the swept phi-interval, in the angle psi = 2 arctan x. Below the
    // barrier the sweep is the arc |psi - psi0| <= 2 arctan x_in about psi0 = 0
    // for the inner branch and about psi0 = pi for the outer one, since the two
    // branches are exchanged by x -> 1/x, that is psi -> pi - psi. Over the
    // barrier x runs over the whole projective line and the arc closes.
    const xIn = M2 > 0 ? Math.sqrt((a2 / 4 - Math.sqrt(M2)) / E1) : null;
    const cover = {
        Fp,
        psi0: outer ? Math.PI : 0,
        half: M2 > 0 ? 2 * Math.atan(xIn) : Math.PI,
        sign: Math.sign(trueX)
    };
    wellPanel('canvas-uduf', ps, vs, -2 * M2, data.phi[i], '\u03c6', '#10b981',
              V_RANGE, cover);
    drawLemniscate(data, i);
}

// -------------------------------------------------------------------- update

function update() {
    const d = derived();
    document.getElementById('a2-val').textContent = a2.toFixed(2);
    document.getElementById('e-val').textContent = E1.toFixed(3);

    const data = integrate();
    current = data;
    clock0 = null;

    const regime = d.M2 > 0 ? 'sub-barrier (dnoidal)'
        : d.M2 < 0 ? 'over-barrier (cnoidal)' : 'separatrix';
    const closure = Number.isFinite(data.closure)
        ? data.closure.toExponential(1) : '&mdash;';
    document.getElementById('readout').innerHTML = `
        <tr><td>F<sub>+</sub><sup>2</sup></td><td>${d.Fp2.toFixed(4)}</td></tr>
        <tr><td>F<sub>&minus;</sub><sup>2</sup></td><td>${d.Fm2.toFixed(4)}</td></tr>
        <tr><td>&#119924;<sup>2</sup></td><td>${d.M2.toFixed(4)}</td></tr>
        <tr><td>k<sup>2</sup></td><td>${d.k2.toFixed(4)}</td></tr>
        <tr><td>period of &phi;</td><td>${d.T.toFixed(4)}</td></tr>
        <tr><td>period of x</td><td>${(d.M2 > 0 ? 2 * d.T : d.T).toFixed(4)}</td></tr>
        <tr><td>window</td><td>${data.nWin} &times; period of x</td></tr>
        <tr><td>loop closure</td><td>${closure}</td></tr>
        <tr><td>regime</td><td>${regime}</td></tr>`;

    drawX2(data);
    drawErr(data);

    // Compare only where both curves are bounded: near a finite-time escape
    // x^2 diverges on both routes and an absolute difference says nothing.
    const CAP = 10;
    let worst = 0, counted = 0;
    for (let i = 0; i < data.ts.length; i++) {
        const a = data.xA[i], b = data.xB[i];
        if (!Number.isFinite(a) || !Number.isFinite(b)) continue;
        if (a > CAP || b > CAP) continue;
        worst = Math.max(worst, Math.abs(a - b));
        counted++;
    }

    const v = document.getElementById('verdict');
    if (data.escaped !== null) {
        v.innerHTML = `With &#119924;<sup>2</sup> &lt; 0 the orbit crosses the
            barrier: x<sup>2</sup> first diverges at t &asymp;
            ${data.escaped.toFixed(3)}, where the point passes through the node
            onto the opposite petal. Where x<sup>2</sup> stays below 10 the two
            routes agree to <strong>${worst.toExponential(2)}</strong>.`;
    } else {
        v.innerHTML = `Largest discrepancy over ${data.phiPeriods} periods of
            &phi;, that is ${data.nWin} full ${data.nWin === 1 ? 'period' : 'periods'}
            of x:
            <strong>${worst.toExponential(2)}</strong>.`;
    }
}

function animate(now) {
    if (current) {
        if (clock0 === null) clock0 = now;
        const span = Math.max(current.ts.length - 1, 1);
        const frac = ((now - clock0) / (PLAY_SECONDS * 1000)) % 1;
        drawWells(current, Math.floor(frac * span));
    }
    requestAnimationFrame(animate);
}

document.addEventListener('DOMContentLoaded', () => {
    const a2s = document.getElementById('a2-slider');
    const es = document.getElementById('e-slider');
    a2 = parseFloat(a2s.value);
    E1 = parseFloat(es.value);

    // The U(x) panel keeps a fixed vertical scale, so that moving
    // alpha_1D raises the barrier against a stationary frame. The range must
    // hold the tallest barrier and the highest energy the sliders allow.
    const aMax = parseFloat(a2s.max), aMin = parseFloat(a2s.min);
    const eMax = parseFloat(es.max);

    const top = 1.06 * Math.max(aMax / 4, eMax);
    U_RANGE = [-0.05 * top, top];

    // The V(phi) panel is fixed too: deepest well is -aMax^2/8, the highest
    // level is -2 M^2 at the smallest alpha and the largest energy, and the
    // turning point never exceeds F_+ = sqrt(aMax/4 + eMax).
    PHI_MAX = 1.15 * Math.sqrt(aMax / 4 + eMax);
    const vLo = -1.06 * aMax * aMax / 8;
    const vHi = 1.06 * 2 * (eMax * eMax - (aMin / 4) ** 2);
    V_RANGE = [vLo, Math.max(vHi, 0.15)];
    const ot = document.getElementById('outer-toggle');
    if (ot) {
        outer = ot.checked;
        ot.addEventListener('change', (e) => { outer = e.target.checked; update(); });
    }

    a2s.addEventListener('input', (e) => { a2 = parseFloat(e.target.value); update(); });
    es.addEventListener('input', (e) => setE1(parseFloat(e.target.value)));
    attachLevelDrag('canvas-ulem');
    attachLevelDrag('canvas-uduf');
    update();
    requestAnimationFrame(animate);
});
