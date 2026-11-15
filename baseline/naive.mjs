// baseline/naive.mjs — independent reference for src/core.mjs.
//
// This file imports NOTHING. It re-implements the ray geometry from scratch and
// evaluates the optical depth with a dense uniform composite trapezoid rule,
// then *predicts its own discretisation error* from the Euler-Maclaurin
// expansion. No quadrature node, tolerance or magic constant is shared with
// src/core.mjs, so agreement between the two is evidence rather than tautology.

/** Independent ray radius: |o + s d| with |o| = r0 and o.d = r0*mu. */
export function radiusN(r0, mu, s) {
  return Math.sqrt(r0 * r0 + 2 * r0 * mu * s + s * s);
}

/** Independent truncation: double until the sampled value is below 1e-18. */
export function cutoffN(r0, mu, R) {
  let s = 1;
  while (Math.exp(-(radiusN(r0, mu, s) - R)) > 1e-18 && s < 1e7) s *= 2;
  return s;
}

/** Dense composite trapezoid rule for tau, uniform step h = sMax / steps. */
export function trapezoidTau(r0, mu, R, { steps = 1 << 14, sMax = null } = {}) {
  const b = sMax ?? cutoffN(r0, mu, R);
  const h = b / steps;
  const f = (s) => Math.exp(-(radiusN(r0, mu, s) - R));
  let sum = 0.5 * (f(0) + f(b));
  for (let i = 1; i < steps; i++) sum += f(i * h);
  return sum * h;
}

/**
 * Euler-Maclaurin estimate of the trapezoid rule's own error:
 *   T - ∫ = (h^2/12)(f'(b) - f'(a)) - (h^4/720)(f'''(b) - f'''(a)) + ...
 * f(s) = exp(-(r(s)-R)),  f'(s) = -f(s) * (r0*mu + s)/r(s).
 */
export function trapezoidErrorEstimate(r0, mu, R, { steps = 1 << 14, sMax = null } = {}) {
  const b = sMax ?? cutoffN(r0, mu, R);
  const h = b / steps;
  const f = (s) => Math.exp(-(radiusN(r0, mu, s) - R));
  const df = (s) => -f(s) * ((r0 * mu + s) / radiusN(r0, mu, s));
  const d3f = (s) => {
    const e = 1e-3 * Math.max(1, b * 1e-3);
    return (df(s + 2 * e) - 2 * df(s + e) + 2 * df(s - e) - df(s - 2 * e)) / (2 * e * e * e);
  };
  const em2 = (h * h) / 12 * (df(b) - df(0));
  const em4 = -(h ** 4) / 720 * (d3f(b) - d3f(0));
  return em2 + em4;
}

/** Trapezoid rule with the Euler-Maclaurin correction applied. */
export function correctedTau(r0, mu, R, opts) {
  return trapezoidTau(r0, mu, R, opts) - trapezoidErrorEstimate(r0, mu, R, opts);
}

/** Trapezoid rule over ONLY [sMax, factor*sMax]; no cancellation with the head. */
export function trapezoidTail(r0, mu, R, { sMax, factor = 8, steps = 1 << 16 } = {}) {
  const a = sMax;
  const b = sMax * factor;
  const h = (b - a) / steps;
  const f = (s) => Math.exp(-(radiusN(r0, mu, s) - R));
  let sum = 0.5 * (f(a) + f(b));
  for (let i = 1; i < steps; i++) sum += f(a + i * h);
  return sum * h;
}
