// P_4(W) transported along a boost orbit.
//
// The first three sliders fix a two-dimensional oscillator,
//     alpha = sqrt(4 L^2) + dA,      E = alpha L + dE,      L^2 given,
// so that E always clears the minimum of U_eff and the motion exists. From it
// we read the two invariants of the orbit,
//     alpha_1D^2 = alpha^2 - 4 L^2,  E_1D = sqrt(E^2 - alpha^2 L^2),
// and transport the system with
//     E(eta)   = E_1D cosh 4eta + (alpha_1D^2/4) sinh 4eta,
//     alpha^2(eta) = 4[(alpha_1D^2/8) cosh 4eta + (E_1D/2) sinh 4eta + alpha_1D^2/8],
//     L^2(eta)     =  (alpha_1D^2/8) cosh 4eta + (E_1D/2) sinh 4eta - alpha_1D^2/8.
// At eta = 0 this is the representative with L = 0, where P_4 degenerates to
// the cubic P_1D.

let LL0 = 1.0, dA = 0.5, dE = 0.3, frac = 1.0, wMax = 6;

function invariants() {
    const al = Math.sqrt(4 * LL0) + dA;
    const a2 = al * al;
    const E0 = al * Math.sqrt(LL0) + dE;
    const a1 = a2 - 4 * LL0;                       // alpha_1D^2
    const E1 = Math.sqrt(Math.max(E0 * E0 - a2 * LL0, 0));
    return { al, a2, E0, a1, E1 };
}

function transport(eta, a1, E1) {
    const c = Math.cosh(4 * eta), s = Math.sinh(4 * eta);
    return {
        E: E1 * c + (a1 / 4) * s,
        aa: 4 * ((a1 / 8) * c + (E1 / 2) * s + a1 / 8),
        LL: (a1 / 8) * c + (E1 / 2) * s - a1 / 8
    };
}

// the rapidity at which the transported system is the one chosen by the sliders
function etaStar(a1, E1) {
    if (LL0 <= 0) return 0;
    const f = (e) => transport(e, a1, E1).LL - LL0;
    let lo = 0, hi = 0.01;
    while (f(hi) < 0 && hi < 40) hi *= 2;
    for (let i = 0; i < 200; i++) {
        const mid = 0.5 * (lo + hi);
        (f(mid) < 0 ? lo = mid : hi = mid);
    }
    return 0.5 * (lo + hi);
}

function P4(W, LL, E, aa) {
    return -4 * LL * W ** 4 + 8 * E * W ** 3
        - 4 * (aa + 2 * LL) * W ** 2 + 8 * E * W - 4 * LL;
}

// real roots in [wLo, wHi], by sign change and bisection
function roots(LL, E, aa, wLo, wHi) {
    const out = [];
    const N = 4000;
    let prevW = wLo, prevV = P4(wLo, LL, E, aa);
    for (let i = 1; i <= N; i++) {
        const w = wLo + (wHi - wLo) * i / N;
        const v = P4(w, LL, E, aa);
        if (prevV === 0) out.push(prevW);
        else if (prevV * v < 0) {
            let a = prevW, b = w;
            for (let j = 0; j < 80; j++) {
                const m = 0.5 * (a + b);
                (P4(a, LL, E, aa) * P4(m, LL, E, aa) <= 0 ? b = m : a = m);
            }
            out.push(0.5 * (a + b));
        }
        prevW = w; prevV = v;
    }
    return out;
}

// ------------------------------------------------------------------ drawing

function draw() {
    const cv = document.getElementById('cv');
    const ctx = cv.getContext('2d');
    const W = cv.width, H = cv.height;
    const P = { l: 62, r: 18, t: 18, b: 34 };
    const pw = W - P.l - P.r, ph = H - P.t - P.b;

    const { a1, E1 } = invariants();
    const eStar = etaStar(a1, E1);
    const eta = frac * eStar;
    const cur = transport(eta, a1, E1);
    const one = transport(0, a1, E1);

    const wLo = -0.5 * wMax / 6, wHi = wMax;

    // vertical range from the two curves in view
    let yAbs = 1e-9;
    for (let i = 0; i <= 600; i++) {
        const w = wLo + (wHi - wLo) * i / 600;
        yAbs = Math.max(yAbs,
            Math.abs(P4(w, cur.LL, cur.E, cur.aa)),
            Math.abs(P4(w, one.LL, one.E, one.aa)));
    }
    yAbs = Math.min(yAbs, 1e6);

    const toX = (w) => P.l + ((w - wLo) / (wHi - wLo)) * pw;
    const toY = (v) => P.t + ph / 2 - (v / yAbs) * (ph / 2) * 0.92;

    ctx.clearRect(0, 0, W, H);

    // grid
    ctx.font = '11px Inter';
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    for (let i = -2; i <= 2; i++) {
        const v = yAbs * i / 2, y = toY(v);
        ctx.strokeStyle = i === 0 ? 'rgba(255,255,255,0.22)' : 'rgba(255,255,255,0.07)';
        ctx.lineWidth = 1;
        ctx.beginPath(); ctx.moveTo(P.l, y); ctx.lineTo(P.l + pw, y); ctx.stroke();
        ctx.fillStyle = '#64748b';
        ctx.fillText(v.toPrecision(2), P.l - 8, y);
    }
    ctx.strokeStyle = 'rgba(255,255,255,0.18)';
    ctx.beginPath(); ctx.moveTo(toX(0), P.t); ctx.lineTo(toX(0), P.t + ph); ctx.stroke();

    ctx.save();
    ctx.beginPath(); ctx.rect(P.l, P.t, pw, ph); ctx.clip();

    // allowed region of the current system
    const rs = roots(cur.LL, cur.E, cur.aa, Math.max(wLo, 1e-9), wHi);
    ctx.fillStyle = 'rgba(16,185,129,0.10)';
    const bounds = [Math.max(wLo, 0), ...rs, wHi];
    for (let i = 0; i < bounds.length - 1; i++) {
        const mid = 0.5 * (bounds[i] + bounds[i + 1]);
        if (P4(mid, cur.LL, cur.E, cur.aa) > 0)
            ctx.fillRect(toX(bounds[i]), P.t, toX(bounds[i + 1]) - toX(bounds[i]), ph);
    }

    const curve = (st, colour, width, dash) => {
        ctx.strokeStyle = colour; ctx.lineWidth = width; ctx.setLineDash(dash);
        ctx.beginPath();
        for (let i = 0; i <= 1200; i++) {
            const w = wLo + (wHi - wLo) * i / 1200;
            const y = toY(P4(w, st.LL, st.E, st.aa));
            i ? ctx.lineTo(toX(w), y) : ctx.moveTo(toX(w), y);
        }
        ctx.stroke(); ctx.setLineDash([]);
    };

    curve(one, '#64748b', 1.5, [6, 4]);      // the L = 0 profile, for reference
    curve(cur, '#10b981', 2.5, []);

    ctx.fillStyle = '#10b981';
    for (const r of rs) {
        ctx.beginPath(); ctx.arc(toX(r), toY(0), 4.5, 0, Math.PI * 2); ctx.fill();
    }
    ctx.restore();

    ctx.textAlign = 'right'; ctx.textBaseline = 'alphabetic';
    ctx.fillStyle = '#94a3b8';
    ctx.fillText('W', W - P.r, H - 10);

    return { eta, eStar, cur, one, rs, a1, E1 };
}

// ---------------------------------------------------------- potential panels

// U_eff(r) = L^2/(2 s) + (alpha^2/2) s ,  s = m(r) r^2 = r^2/(1+r^4).
// For L = 0 this is the one-dimensional potential (alpha_1D^2/2) s. The turning
// points are the positive roots of P_4, so they are taken from the same routine
// that marks the roots on the main plot.
function Ueff(r, LL, aa) {
    const s = r * r / (1 + r ** 4);
    if (s <= 0) return LL > 0 ? Infinity : 0;
    return LL / (2 * s) + (aa / 2) * s;
}

function potentialPanel(id, LL, E, aa, rMax, colour, xLabel) {
    const cv = document.getElementById(id);
    if (!cv) return;
    const ctx = cv.getContext('2d');
    const W = cv.width, H = cv.height;
    const P = { l: 52, r: 14, t: 14, b: 30 };
    const pw = W - P.l - P.r, ph = H - P.t - P.b;

    const C0 = LL + aa / 4;                       // the barrier at r = 1
    const yHi = Math.max(1.35 * Math.max(E, C0), 1e-6);

    const toX = (r) => P.l + (r / rMax) * pw;
    const toY = (u) => P.t + ph - (u / yHi) * ph;

    ctx.clearRect(0, 0, W, H);

    ctx.font = '11px Inter';
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    for (let i = 0; i <= 4; i++) {
        const u = yHi * i / 4, y = toY(u);
        ctx.strokeStyle = 'rgba(255,255,255,0.06)';
        ctx.lineWidth = 1;
        ctx.beginPath(); ctx.moveTo(P.l, y); ctx.lineTo(P.l + pw, y); ctx.stroke();
        ctx.fillStyle = '#64748b';
        ctx.fillText(u.toPrecision(2), P.l - 7, y);
    }
    ctx.strokeStyle = 'rgba(255,255,255,0.18)';
    ctx.beginPath();
    ctx.moveTo(P.l, P.t); ctx.lineTo(P.l, P.t + ph); ctx.lineTo(P.l + pw, P.t + ph);
    ctx.stroke();

    ctx.save();
    ctx.beginPath(); ctx.rect(P.l, P.t, pw, ph); ctx.clip();

    // classically allowed band, U <= E
    ctx.fillStyle = 'rgba(16,185,129,0.10)';
    let run = null;
    for (let i = 0; i <= 900; i++) {
        const r = rMax * i / 900;
        const inside = Ueff(r, LL, aa) <= E;
        if (inside && run === null) run = toX(r);
        if ((!inside || i === 900) && run !== null) {
            ctx.fillRect(run, P.t, toX(r) - run, ph); run = null;
        }
    }

    // the barrier at r = 1
    ctx.strokeStyle = 'rgba(234,179,8,0.35)';
    ctx.setLineDash([3, 4]); ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(toX(1), P.t); ctx.lineTo(toX(1), P.t + ph); ctx.stroke();
    ctx.setLineDash([]);

    // the potential
    ctx.strokeStyle = colour; ctx.lineWidth = 2;
    ctx.beginPath();
    let pen = false;
    for (let i = 0; i <= 1400; i++) {
        const r = rMax * i / 1400;
        const u = Ueff(r, LL, aa);
        if (!Number.isFinite(u) || u > yHi * 1.6) { pen = false; continue; }
        const x = toX(r), y = toY(u);
        if (!pen) { ctx.moveTo(x, y); pen = true; } else ctx.lineTo(x, y);
    }
    ctx.stroke();

    // energy level
    ctx.strokeStyle = 'rgba(168,85,247,0.9)';
    ctx.setLineDash([5, 4]); ctx.lineWidth = 1.5;
    ctx.beginPath(); ctx.moveTo(P.l, toY(E)); ctx.lineTo(P.l + pw, toY(E)); ctx.stroke();
    ctx.setLineDash([]);

    // Turning points: the non-negative roots of P_4, as radii. W = 0 is a root
    // exactly when L = 0, and then it is the inner turning point, the orbit
    // reaching the origin.
    const rs = roots(LL, E, aa, 0, rMax * rMax);
    ctx.fillStyle = '#a855f7';
    for (const w of rs) {
        const r = Math.sqrt(w);
        if (r > rMax) continue;
        ctx.beginPath(); ctx.arc(toX(r), toY(E), 4.5, 0, Math.PI * 2); ctx.fill();
    }
    ctx.restore();

    ctx.textAlign = 'right'; ctx.textBaseline = 'alphabetic';
    ctx.fillStyle = '#94a3b8';
    ctx.fillText(xLabel, W - P.r, H - 9);
    return rs.length;
}

// ------------------------------------------------------------------- update

function update() {
    document.getElementById('ll-val').textContent = LL0.toFixed(2);
    document.getElementById('da-val').textContent = dA.toFixed(2);
    document.getElementById('de-val').textContent = dE.toFixed(2);
    document.getElementById('s-val').textContent = frac.toFixed(3);
    document.getElementById('wm-val').textContent = wMax;

    const d = draw();
    const rMax = Math.sqrt(wMax);
    const n2d = potentialPanel('cv2d', d.cur.LL, d.cur.E, d.cur.aa, rMax, '#10b981', 'r');
    const n1d = potentialPanel('cv1d', d.one.LL, d.one.E, d.one.aa, rMax, '#64748b', 'x');
    const C0 = d.cur.LL + d.cur.aa / 4;
    const M2 = C0 * C0 - d.cur.E * d.cur.E;

    document.getElementById('state').innerHTML = `
        <tr><td>&eta;<sub>s</sub></td><td>${d.eta.toFixed(4)}</td></tr>
        <tr><td>&eta;<sub>s</sub><sup>&#9733;</sup></td><td>${d.eStar.toFixed(4)}</td></tr>
        <tr><td>E</td><td>${d.cur.E.toFixed(4)}</td></tr>
        <tr><td>&alpha;<sup>2</sup></td><td>${d.cur.aa.toFixed(4)}</td></tr>
        <tr><td>L<sup>2</sup></td><td>${d.cur.LL.toExponential(3)}</td></tr>`;

    document.getElementById('inv').innerHTML = `
        <tr><td>&alpha;<sub>1D</sub><sup>2</sup></td><td>${d.a1.toFixed(4)}</td></tr>
        <tr><td>E<sub>1D</sub></td><td>${d.E1.toFixed(4)}</td></tr>
        <tr><td>&#119924;<sup>2</sup></td><td>${M2.toFixed(4)}</td></tr>
        <tr><td>real roots</td><td>${d.rs.length} in view</td></tr>
        <tr><td>turning pts 2D / 1D</td><td>${n2d} / ${n1d}</td></tr>`;

    const sector = M2 > 0 ? 'timelike, four real roots'
        : M2 < 0 ? 'spacelike, two real roots' : 'light-like';
    document.getElementById('verdict').innerHTML =
        `Sector: <strong>${sector}</strong>. As &eta;<sub>s</sub> &rarr; 0 the outer
         pair of roots runs to 0 and &infin;, the quartic loses its leading
         coefficient &minus;4L<sup>2</sup> and degenerates to the cubic
         P<sub>1D</sub>, drawn dashed.`;
}

document.addEventListener('DOMContentLoaded', () => {
    const bind = (id, set) => {
        const el = document.getElementById(id);
        set(parseFloat(el.value));
        el.addEventListener('input', (e) => { set(parseFloat(e.target.value)); update(); });
    };
    bind('ll', (v) => LL0 = v);
    bind('da', (v) => dA = v);
    bind('de', (v) => dE = v);
    bind('s', (v) => frac = v);
    bind('wm', (v) => wMax = v);
    update();
});
