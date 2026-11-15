// baseline/simpson.mjs — a SECOND independent oracle, with a different
// quadrature family (composite Simpson + Richardson extrapolation) and its own
// geometry code. It imports nothing, and shares no node, step count or
// tolerance with src/core.mjs or baseline/naive.mjs, so it can disagree for
// real reasons.

/** Independent ray radius. */
export function radiusS(r0, mu, s) {
  return Math.sqrt(r0 * r0 + 2 * r0 * mu * s + s * s);
}

/** Composite Simpson rule with an even number n of intervals on [0, b]. */
export function simpsonTau(r0, mu, R, { n = 4096, b = 4000 } = {}) {
  const h = b / n;
  const f = (s) => Math.exp(-(radiusS(r0, mu, s) - R));
  let acc = f(0) + f(b);
  for (let i = 1; i < n; i++) acc += f(i * h) * (i % 2 ? 4 : 2);
  return (acc * h) / 3;
}

/**
 * Richardson extrapolation of two Simpson results at n and 2n intervals.
 * Simpson's error is O(h^4), so halving h divides it by 16 and the correct
 * combination is (16*S(2n) - S(n)) / 15.
 */
export function richardsonTau(r0, mu, R, { n = 4096, b = 4000 } = {}) {
  const s1 = simpsonTau(r0, mu, R, { n, b });
  const s2 = simpsonTau(r0, mu, R, { n: 2 * n, b });
  return (16 * s2 - s1) / 15;
}
