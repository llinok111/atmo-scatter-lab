// bench/bench.mjs — timing of the core quadrature vs the dense trapezoid
// baseline. Writes bench/result.json with cold start, median, p95, stddev
// over >=5 repetitions.

import { writeFileSync, mkdirSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";
import { opticalDepth, opticalDepthPlaneParallel, horizontalLimitHalfGaussian } from "../src/core.mjs";
import { trapezoidTau } from "../baseline/naive.mjs";

const here = dirname(fileURLToPath(import.meta.url));
const root = join(here, "..");
const R_h = 6371 / 8;
const OPTS = { panels: 12, nodes: 16 };
const STEPS = 1 << 14;

const CASES = [
  { r0: R_h, mu: 1 },
  { r0: R_h, mu: Math.cos((60 * Math.PI) / 180) },
  { r0: R_h, mu: Math.cos((85 * Math.PI) / 180) },
  { r0: R_h, mu: 0 },
  { r0: R_h + 5, mu: Math.cos((45 * Math.PI) / 180) },
];

function runCore() {
  let acc = 0;
  for (const c of CASES) acc += opticalDepth(c.r0, c.mu, R_h, OPTS);
  return acc;
}
function runNaive() {
  let acc = 0;
  for (const c of CASES) acc += trapezoidTau(c.r0, c.mu, R_h, { steps: STEPS });
  return acc;
}

function stats(ms) {
  const s = [...ms].sort((a, b) => a - b);
  const n = s.length;
  const mean = s.reduce((a, b) => a + b, 0) / n;
  const variance = s.reduce((a, b) => a + (b - mean) ** 2, 0) / n;
  const p95 = s[Math.min(n - 1, Math.ceil(0.95 * n) - 1)];
  return { min: s[0], median: (s[Math.floor((n - 1) / 2)] + s[Math.ceil((n - 1) / 2)]) / 2, p95, max: s[n - 1], mean, stddev: Math.sqrt(variance), relStddev: Math.sqrt(variance) / mean };
}

const REPS = 7;
const coreMs = [];
const naiveMs = [];
let coldCore = 0;
let coldNaive = 0;

// warmup 1 pass
runCore();
runNaive();

for (let i = 0; i < REPS; i++) {
  let t0 = process.hrtime.bigint();
  runCore();
  let t1 = process.hrtime.bigint();
  coreMs.push(Number(t1 - t0) / 1e6);

  t0 = process.hrtime.bigint();
  runNaive();
  t1 = process.hrtime.bigint();
  naiveMs.push(Number(t1 - t0) / 1e6);

  if (i === 0) {
    coldCore = coreMs[0];
    coldNaive = naiveMs[0];
  }
}

const coreStats = stats(coreMs);
const naiveStats = stats(naiveMs);

const horizontal = opticalDepth(R_h, 0, R_h, OPTS);
const result = {
  node: process.version,
  generatedAt: new Date().toISOString(),
  cases: CASES.length,
  coreQuadrature: { panels: OPTS.panels, nodes: OPTS.nodes, evaluationsPerCase: OPTS.panels * OPTS.nodes },
  baseline: { rule: "uniform composite trapezoid", stepsPerCase: STEPS },
  evalRatio: STEPS / (OPTS.panels * OPTS.nodes),
  coldStartMs: { core: coldCore, baseline: coldNaive },
  samples: REPS,
  coreMs,
  baselineMs: naiveMs,
  core: coreStats,
  baseline: naiveStats,
  medianSpeedup: naiveStats.median / coreStats.median,
  p95Speedup: naiveStats.p95 / coreStats.p95,
  numbers: {
    verticalTau: opticalDepth(R_h, 1, R_h, OPTS),
    horizontalTau: horizontal,
    horizontalHalfGaussianLimit: horizontalLimitHalfGaussian(R_h),
    horizontalRelDiffToLimit: (horizontal - horizontalLimitHalfGaussian(R_h)) / horizontalLimitHalfGaussian(R_h),
    planeParallelAt89Deg: opticalDepthPlaneParallel(R_h, Math.cos((89 * Math.PI) / 180), R_h),
    sphericalAt89Deg: opticalDepth(R_h, Math.cos((89 * Math.PI) / 180), R_h, OPTS),
  },
};

mkdirSync(here, { recursive: true });
writeFileSync(join(here, "result.json"), JSON.stringify(result, null, 2) + "\n");

console.log(`core   median ${coreStats.median.toFixed(4)} ms  (min ${coreStats.min.toFixed(4)}, p95 ${coreStats.p95.toFixed(4)}, relSD ${(coreStats.relStddev * 100).toFixed(2)}%)`);
console.log(`base   median ${naiveStats.median.toFixed(4)} ms  (min ${naiveStats.min.toFixed(4)}, p95 ${naiveStats.p95.toFixed(4)}, relSD ${(naiveStats.relStddev * 100).toFixed(2)}%)`);
console.log(`cold   core ${coldCore.toFixed(4)} ms / baseline ${coldNaive.toFixed(4)} ms`);
console.log(`speedup median ${result.medianSpeedup.toFixed(2)}x  p95-vs-p95 ${result.p95Speedup.toFixed(2)}x  (${result.evalRatio}x fewer integrand evaluations)`);
console.log(`wrote ${join(here, "result.json")}`);
