// src/core.mjs — atmospheric single-scattering transmittance, pure computation.
//
// This module is deliberately free of three.js / DOM / window / navigator /
// requestAnimationFrame and of any top-level side effect. Everything here is
// plain arithmetic on numbers, so `node -e "import('./src/core.mjs')"` works and
// the browser renderer can consume exactly the same results (see main.js).
//
// Conventions (the whole module is dimensionless):
//   * Lengths are expressed in units of the atmospheric scale height H.
//   * Density is normalised to rho0 = 1 at the reference radius R.
//   * Therefore the returned optical depth is tau / (rho0 * H) and the
//     transmittance is exp(-tau) with no hidden unit conversion.
// Only two geometric ratios matter: R/H and r0/H. There is no absolute length
// constant anywhere in this file, which is what makes the k-scaling invariance
// in bench/verify.mjs a structural property rather than a numerical accident.

// ---------------------------------------------------------------------------
// Gauss-Legendre quadrature on [-1, 1]
// ---------------------------------------------------------------------------

const _glCache = new Map();

/** Nodes/weights of the n-point Gauss-Legendre rule on [-1,1]. */
export function gaussLegendreNodes(n) {
  const hit = _glCache.get(n);
  if (hit) return hit;
  const x = new Float64Array(n);
  const w = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    // Chebyshev initial guess (Abramowitz & Stegun 25.4.30)
    let t = Math.cos((Math.PI * (i + 0.75)) / (n + 0.5));
    let dp = 1;
    for (let it = 0; it < 100; it++) {
      let p0 = 1;
      let p1 = t;
      for (let k = 2; k <= n; k++) {
        const p2 = ((2 * k - 1) * t * p1 - (k - 1) * p0) / k;
        p0 = p1;
        p1 = p2;
      }
      const Pn = n === 1 ? t : p1;
      const Pn1 = n === 1 ? 1 : p0;
      dp = (n * (t * Pn - Pn1)) / (t * t - 1);
      const dt = Pn / dp;
      t -= dt;
      if (Math.abs(dt) < 1e-15) break;
    }
    x[i] = t;
    w[i] = 2 / ((1 - t * t) * dp * dp);
  }
  const out = { x, w };
  _glCache.set(n, out);
  return out;
}

/**
 * Composite Gauss-Legendre integral of f over [a, b] with `panels` panels of
 * `nodes` points each.
 */
export function compositeGL(f, a, b, panels, nodes) {
  const { x, w } = gaussLegendreNodes(nodes);
  const h = (b - a) / panels;
  let s = 0;
  for (let p = 0; p < panels; p++) {
    const c = a + h * (p + 0.5);
    const r = h / 2;
    for (let i = 0; i < nodes; i++) s += w[i] * f(c + r * x[i]);
  }
  return s * (h / 2);
}

// ---------------------------------------------------------------------------
// Ray geometry in spherical exponential atmosphere
// ---------------------------------------------------------------------------

/**
 * Radius of the point o + s*d, where |o| = r0 and o·d = r0*mu.
 * @param {number} r0 start radius (H units)
 * @param {number} mu cos(zenith angle at the start point), in [0, 1]
 * @param {number} s  distance travelled along the ray (H units)
 */
export function rayRadius(r0, mu, s) {
  return Math.sqrt(r0 * r0 + 2 * r0 * mu * s + s * s);
}

/** Local density of an exponential atmosphere: rho(h) = exp(-h), h in H units. */
export function densityAtRadius(r, R) {
  return Math.exp(-(r - R));
}

/** d/ds of the ray height above R. Non-negative for mu >= 0 and s >= 0. */
export function rayHeightDerivative(r0, mu, s) {
  return (r0 * mu + s) / rayRadius(r0, mu, s);
}

/**
 * Rigorous truncation plan. For mu >= 0 the height h(s) = |o+sd| - R has a
 * non-decreasing derivative, hence
 *   h(s) - h(sMax) >= h'(sMax) * (s - sMax)   for s >= sMax,
 * so the discarded tail obeys
 *   ∫_{sMax}^∞ exp(-h(s)) ds <= exp(-h(sMax)) / h'(sMax) =: tailBound.
 * sMax is doubled until tailBound <= tol. This is a bound, not a fitted
 * constant: the bench prints the bound next to the measured truncation effect.
 */
export function truncationPlan(r0, mu, R, tol = 1e-16) {
  let sMax = 1;
  for (let i = 0; i < 400; i++) {
    const f = Math.exp(-(rayRadius(r0, mu, sMax) - R));
    const d = Math.max(rayHeightDerivative(r0, mu, sMax), Number.MIN_VALUE);
    const tailBound = f / d;
    if (tailBound <= tol) return { sMax, tailBound, iterations: i };
    sMax *= 2;
  }
  return { sMax, tailBound: Math.exp(-(rayRadius(r0, mu, sMax) - R)), iterations: 400 };
}

/**
 * Optical depth of an exponential atmosphere along a ray that never descends
 * into the planet: r0 >= R and mu >= 0.
 *
 *   tau = ∫_0^∞ exp(-(|o + s d| - R)) ds        (units of rho0 * H)
 *
 * The ray is truncated at an analytically bounded sMax (see truncationPlan)
 * and integrated with composite Gauss-Legendre; the integrand is analytic on
 * the finite interval, so convergence is geometric.
 *
 * @param {number} r0 radius where the ray starts (H units), r0 >= R
 * @param {number} mu cos(zenith angle at the start), 0 <= mu <= 1
 * @param {number} R  reference (ground) radius (H units)
 */
export function opticalDepth(r0, mu, R, opts = {}) {
  const panels = opts.panels ?? 12;
  const nodes = opts.nodes ?? 16;
  const tol = opts.tol ?? 1e-16;
  if (!(r0 >= R - 1e-12)) throw new Error("opticalDepth: ray starts below the reference radius");
  if (mu < -1e-12 || mu > 1 + 1e-12) throw new Error("opticalDepth: mu out of range");
  const { sMax } = truncationPlan(r0, mu, R, tol);
  return compositeGL((s) => densityAtRadius(rayRadius(r0, mu, s), R), 0, sMax, panels, nodes);
}

/** Transmittance T = exp(-tau). */
export function transmittance(r0, mu, R, opts) {
  return Math.exp(-opticalDepth(r0, mu, R, opts));
}

/**
 * Physical entry point: lengths in the same unit (e.g. km), cos(zenith) taken
 * from degrees. This is the only place where a unit conversion happens, so the
 * k-scaling invariance checked in bench/verify.mjs exercises real rounding of
 * r0/H and R/H rather than algebraically identical arguments.
 */
export function opticalDepthPhysical({ R, H, r0, zenithDeg, height = 0, opts }) {
  const R_h = R / H;
  const r0_h = (r0 ?? R) / H + height / H;
  const mu = Math.cos((zenithDeg * Math.PI) / 180);
  return opticalDepth(r0_h, mu, R_h, opts);
}

/**
 * Plane-parallel (flat-atmosphere) optical depth: tau = rho(h0) * sec(chi).
 * This IS elementary closed form and is the approximation that every
 * cheap sky model uses; the point of this repo is to quantify where it breaks.
 */
export function opticalDepthPlaneParallel(r0, mu, R) {
  if (mu <= 0) return Infinity;
  return Math.exp(-(r0 - R)) / mu;
}

/**
 * Horizontal optical depth of an exponential atmosphere, for a ray that starts
 * at radius R and leaves the planet.
 *
 * Exact integral is R * K_1(R) * exp(R); with h(s) = sqrt(R^2+s^2) - R the
 * integrand becomes exp(-s^2/(2R)) and the ray covers only s >= 0, so the
 * closed form is the HALF Gaussian integral:
 *   tau_h = sqrt(pi * R / 2)      (R in H units, i.e. sqrt(pi*R_km*H/2) physical)
 * The exact value is larger by a relative O(H/R) term which the bench measures.
 */
export function horizontalLimitHalfGaussian(R) {
  return Math.sqrt((Math.PI * R) / 2);
}

// ---------------------------------------------------------------------------
// Scattering phase functions and their exact normalisations
// ---------------------------------------------------------------------------

/** Rayleigh phase function, normalised so that the sphere integral is 1. */
export function rayleighPhase(mu) {
  return (3 / (16 * Math.PI)) * (1 + mu * mu);
}

/** Exact integral of the Rayleigh phase over the sphere (analytic). */
export function rayleighSphereIntegral() {
  return 1;
}

/** Henyey-Greenstein phase function with asymmetry g in [0,1). */
export function hgPhase(mu, g) {
  const d = 1 + g * g - 2 * g * mu;
  return ((1 - g * g) / (4 * Math.PI)) / (d * Math.sqrt(d));
}

/**
 * Exact sphere integral of the HG phase function, evaluated from the
 * closed-form primitive rather than hard coded:
 *   ∫ p(μ) dμ = (1-g²)/(4π) * [ d(μ)^{-1/2} / g ]_{-1}^{1},  d = 1+g²-2gμ
 * which simplifies to 1/(2π); multiplied by 2π this is exactly 1 for every
 * g in [0,1). Returning the primitive-built value (not the literal 1) lets the
 * bench show the cancellation is exact.
 */
export function hgSphereIntegral(g) {
  if (g === 0) return 1; // isotropic limit, primitive is singular in g
  const invSqrtD = (mu) => 1 / Math.sqrt(1 + g * g - 2 * g * mu);
  const integralOverMu = ((1 - g * g) / (4 * Math.PI)) * ((invSqrtD(1) - invSqrtD(-1)) / g);
  return 2 * Math.PI * integralOverMu;
}

/** Numerically integrate a phase function over the sphere with n-point GL. */
export function phaseSphereIntegralNumeric(phase, n) {
  const { x, w } = gaussLegendreNodes(n);
  let s = 0;
  for (let i = 0; i < n; i++) s += w[i] * phase(x[i]);
  return 2 * Math.PI * s;
}

// ---------------------------------------------------------------------------
// Rayleigh scattering cross-section (to turn geometric tau into colour)
// ---------------------------------------------------------------------------

/**
 * Rayleigh scattering cross-section
 *   sigma(lambda) = 8*pi^3*(n^2 - 1)^2 / (3 * N^2 * lambda^4)
 * with n the refractive index of air and N the molecular number density.
 * Only the ratio between wavelengths is used downstream, so the "sky is blue"
 * statement is checked as a ratio, not as an absolute calibration.
 */
export function rayleighCrossSection(lambdaNm, { n = 1.000277, N = 2.504e25 } = {}) {
  const lambdaM = lambdaNm * 1e-9;
  const num = 8 * Math.PI ** 3 * (n * n - 1) ** 2;
  const den = 3 * N * N * lambdaM ** 4;
  return num / den;
}

export const SKY_WAVELENGTHS_NM = [440, 550, 680];

// ---------------------------------------------------------------------------
// Convenience wrappers used by tools/gen-sample.mjs and by the renderer
// ---------------------------------------------------------------------------

/** Build a table of tau over zenith angles (deg) and start heights (H units). */
export function transmittanceTable({ R, zenithDeg, heights, opts }) {
  return zenithDeg.map((deg) => {
    const mu = Math.cos((deg * Math.PI) / 180);
    return heights.map((h) => opticalDepth(R + h, mu, R, opts));
  });
}

/** One-stop parameter bundle so node and the browser share the same numbers. */
export function defaultAtmosphere() {
  return { R: 6371, H: 8, name: "us-standard-ish R=6371 km, H=8 km" };
}
