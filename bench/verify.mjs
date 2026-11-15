// bench/verify.mjs — behavioural assertions with printed values.
// Exits 0 only if every assertion holds. Numbers come from a real run.

import { readFileSync, existsSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";
import {
  opticalDepth,
  opticalDepthPhysical,
  opticalDepthPlaneParallel,
  horizontalLimitHalfGaussian,
  rayleighPhase,
  hgPhase,
  rayleighSphereIntegral,
  hgSphereIntegral,
  phaseSphereIntegralNumeric,
  rayleighCrossSection,
  rayHeightDerivative,
  truncationPlan,
  SKY_WAVELENGTHS_NM,
} from "../src/core.mjs";
import { trapezoidTau, correctedTau, cutoffN, trapezoidTail } from "../baseline/naive.mjs";

const here = dirname(fileURLToPath(import.meta.url));
const root = join(here, "..");
const R_h = 6371 / 8;
const OPTS = { panels: 12, nodes: 16 };

let passed = 0;
let failed = 0;
const rows = [];
function check(name, ok, detail) {
  rows.push({ name, ok, detail });
  if (ok) passed++;
  else {
    failed++;
    console.log(`FAIL ${name}: ${detail}`);
  }
  console.log(`${ok ? "ok  " : "FAIL"} ${name} :: ${detail}`);
}
const rel = (a, b) => Math.abs(a - b) / Math.abs(b);

// ---- G4.1 exact closed-form cases -----------------------------------------
const vertical = opticalDepth(R_h, 1, R_h, OPTS);
check("vertical tau = 1 exactly (∫exp(-s)ds)", rel(vertical, 1) < 1e-11, `tau=${vertical.toFixed(15)} rel=${rel(vertical, 1).toExponential(3)}`);

const horiz = opticalDepth(R_h, 0, R_h, OPTS);
const horizLimit = horizontalLimitHalfGaussian(R_h);
check(
  "horizontal tau within 1e-3 rel of sqrt(pi*R/2)",
  rel(horiz, horizLimit) < 1e-3,
  `tau=${horiz.toFixed(9)} limit=${horizLimit.toFixed(9)} rel=${rel(horiz, horizLimit).toExponential(3)}`,
);

// ---- G4.2 independent oracle + Euler-Maclaurin corrected oracle -----------
const oracleTargets = [
  { r0: R_h, mu: 1 },
  { r0: R_h, mu: Math.cos((60 * Math.PI) / 180) },
  { r0: R_h, mu: Math.cos((85 * Math.PI) / 180) },
  { r0: R_h, mu: 0 },
  { r0: R_h + 5, mu: Math.cos((45 * Math.PI) / 180) },
];
let worstCorr = 0;
let worstRaw = 0;
for (const { r0, mu } of oracleTargets) {
  const mine = opticalDepth(r0, mu, R_h, OPTS);
  const steps = 1 << 14;
  const sMax = cutoffN(r0, mu, R_h);
  const raw = trapezoidTau(r0, mu, R_h, { steps, sMax });
  const corr = correctedTau(r0, mu, R_h, { steps, sMax });
  worstRaw = Math.max(worstRaw, rel(raw, mine));
  worstCorr = Math.max(worstCorr, rel(corr, mine));
}
check("raw trapezoid oracle is NOT already machine-accurate (shows the correction matters)", worstRaw > 1e-9, `worst raw rel dev=${worstRaw.toExponential(3)}`);
check("Euler-Maclaurin-corrected oracle agrees with core to <=1e-12", worstCorr < 1e-12, `worst corrected rel dev=${worstCorr.toExponential(3)}`);

// the EM term is a *prediction*, checked against the measured raw-minus-core gap
let emWorstRel = 0;
let emWorstAbs = 0;
let emCases = 0;
for (const { r0, mu } of oracleTargets) {
  const steps = 1 << 12;
  const sMax = cutoffN(r0, mu, R_h);
  const raw = trapezoidTau(r0, mu, R_h, { steps, sMax });
  const mine = opticalDepth(r0, mu, R_h, OPTS);
  const predicted = raw - correctedTau(r0, mu, R_h, { steps, sMax });
  const measured = raw - mine;
  emWorstAbs = Math.max(emWorstAbs, Math.abs(predicted - measured));
  if (Math.abs(measured) > 1e-9) {
    emCases++;
    emWorstRel = Math.max(emWorstRel, Math.abs(predicted - measured) / Math.abs(measured));
  }
}
check("Euler-Maclaurin prediction matches measured oracle error (rel<=2% where it dominates, abs<=1e-9 everywhere)", emWorstRel < 0.02 && emWorstAbs < 1e-9, `worst rel=${emWorstRel.toExponential(3)} over ${emCases} dominant cases, worst abs=${emWorstAbs.toExponential(3)}`);

// ---- G4.3 truncation bound is a real bound --------------------------------
const plan = truncationPlan(R_h, 0, R_h, 1e-16);
// Integrate ONLY the discarded region [sMax, 8*sMax] so there is no
// catastrophic cancellation against the (35 times larger) head integral.
const tailMeasured = trapezoidTail(R_h, 0, R_h, { sMax: plan.sMax, factor: 8, steps: 1 << 16 });
check("analytic tail bound dominates the measured trapezoid tail", Math.abs(tailMeasured) <= plan.tailBound * 1.05, `|measured tail|=${Math.abs(tailMeasured).toExponential(3)} <= bound=${plan.tailBound.toExponential(3)}`);

// ---- G4.4 phase normalisation: closed form vs fixed-node quadrature --------
check("Rayleigh sphere integral from closed form = 1", rayleighSphereIntegral() === 1, `value=${rayleighSphereIntegral()}`);
const rayGL64 = phaseSphereIntegralNumeric(rayleighPhase, 64);
check("Rayleigh GL64 quadrature = 1 to 1e-14", Math.abs(rayGL64 - 1) < 1e-14, `value=${rayGL64}`);
let hgExactWorst = 0;
for (const g of [0.1, 0.5, 0.9, 0.95]) hgExactWorst = Math.max(hgExactWorst, Math.abs(hgSphereIntegral(g) - 1));
check("HG closed-form primitive integrates to 1 for g in {0.1,0.5,0.9,0.95}", hgExactWorst < 1e-13, `worst |value-1|=${hgExactWorst.toExponential(3)}`);
const hg64_09 = phaseSphereIntegralNumeric((mu) => hgPhase(mu, 0.9), 64);
const hg64_095 = phaseSphereIntegralNumeric((mu) => hgPhase(mu, 0.95), 64);
const hg256_095 = phaseSphereIntegralNumeric((mu) => hgPhase(mu, 0.95), 256);
check("fixed 64-node quadrature under-resolves forward scattering at g=0.95", Math.abs(hg64_095 - 1) > 1e-3, `GL64 g=0.9 err=${(hg64_09 - 1).toExponential(3)} GL64 g=0.95 err=${(hg64_095 - 1).toExponential(3)} GL256 g=0.95 err=${(hg256_095 - 1).toExponential(3)}`);

// ---- G4.5 plane-parallel vs spherical divergence -------------------------
const divExpect = [
  [30, 4.2e-4],
  [60, 3.7e-3],
  [75, 1.7e-2],
  [85, 1.3e-1],
];
let divOk = true;
const divSeen = [];
for (const [deg, expect] of divExpect) {
  const mu = Math.cos((deg * Math.PI) / 180);
  const sph = opticalDepth(R_h, mu, R_h, OPTS);
  const pp = opticalDepthPlaneParallel(R_h, mu, R_h);
  const d = (pp - sph) / sph;
  divSeen.push(`${deg}deg=${(d * 100).toFixed(4)}%`);
  if (Math.abs(d - expect) / expect > 0.15) divOk = false;
}
check("plane-parallel divergence matches the measured curve at 30/60/75/85 deg", divOk, divSeen.join(" "));
const sph89 = opticalDepth(R_h, Math.cos((89 * Math.PI) / 180), R_h, OPTS);
const pp89 = opticalDepthPlaneParallel(R_h, Math.cos((89 * Math.PI) / 180), R_h);
check("near-horizon plane-parallel over-estimates by >100%", (pp89 - sph89) / sph89 > 1, `zenith 89deg spherical=${sph89.toFixed(4)} plane-parallel=${pp89.toFixed(4)} rel=${(((pp89 - sph89) / sph89) * 100).toFixed(2)}%`);
const sph8999 = opticalDepth(R_h, Math.cos((89.99 * Math.PI) / 180), R_h, OPTS);
const pp8999 = opticalDepthPlaneParallel(R_h, Math.cos((89.99 * Math.PI) / 180), R_h);
check("horizontal tau stays finite as mu->0 while plane-parallel diverges", Number.isFinite(sph8999) && pp8999 > 10 * sph8999, `zenith 89.99deg sph=${sph8999.toFixed(4)} pp=${pp8999.toFixed(4)} ratio=${(pp8999 / sph8999).toFixed(1)}x`);

// ---- G4.6 scale invariance (physical units) ------------------------------
const base = opticalDepthPhysical({ R: 6371, H: 8, zenithDeg: 60, opts: OPTS });
let scaleWorst = 0;
const scaleKs = [1e-6, 1e-3, 1e-1, 1, 2, 1e3, 1e6];
for (const k of scaleKs) {
  const t = opticalDepthPhysical({ R: 6371 * k, H: 8 * k, zenithDeg: 60, opts: OPTS });
  scaleWorst = Math.max(scaleWorst, rel(t, base));
}
check("k-scaling of (R,H,ray) by 1e-6..1e6 leaves tau unchanged (rel <=1e-15)", scaleWorst <= 1e-15, `worst rel=${scaleWorst.toExponential(3)} over k in {${scaleKs.join(",")}}`);
let pow2Bit = 0;
for (const m of [-20, -10, 10, 20]) {
  const k = 2 ** m;
  if (opticalDepthPhysical({ R: 6371 * k, H: 8 * k, zenithDeg: 60, opts: OPTS }) === base) pow2Bit++;
}
check("power-of-two k-scaling is bit-identical (4/4)", pow2Bit === 4, `bit-identical=${pow2Bit}/4`);

// ---- G4.7 monotone height derivative (the bound's hypothesis) -------------
let derivOk = true;
for (const mu of [0, 0.3, 1]) {
  let prev = -Infinity;
  for (let s = 0; s <= 300; s += 7.3) {
    const d = rayHeightDerivative(R_h, mu, s);
    if (d < prev - 1e-15) derivOk = false;
    prev = d;
  }
}
check("height derivative is non-decreasing along the ray (validates the tail bound)", derivOk, "checked mu in {0,0.3,1}, s in [0,300]");

// ---- G4.8 cross-section ratio --------------------------------------------
const ratio = rayleighCrossSection(440) / rayleighCrossSection(680);
check("sigma(440)/sigma(680) = (680/440)^4 to 1e-12", rel(ratio, (680 / 440) ** 4) < 1e-12, `ratio=${ratio.toFixed(12)} closed form=${((680 / 440) ** 4).toFixed(12)}`);
check("SKY wavelength list is 3 bands", SKY_WAVELENGTHS_NM.length === 3, SKY_WAVELENGTHS_NM.join(","));

// ---- G7 deterministic sample data ----------------------------------------
const dataPath = join(root, "data", "atmosphere.json");
check("data/atmosphere.json exists", existsSync(dataPath), dataPath);
if (existsSync(dataPath)) {
  const d = JSON.parse(readFileSync(dataPath, "utf8"));
  check("LUT shape is 91 x 25", d.tauLut.length === 91 && d.tauLut[0].length === 25, `${d.tauLut.length} x ${d.tauLut[0].length}`);
  check("LUT tau[0][0] equals the core recomputation", rel(d.tauLut[0][0], opticalDepth(R_h, 1, R_h, OPTS)) < 1e-12, `lut=${d.tauLut[0][0]} core=${opticalDepth(R_h, 1, R_h, OPTS)}`);
  check("LUT transmittance = exp(-tau)", rel(d.transmittanceLut[0][0], Math.exp(-d.tauLut[0][0])) < 1e-15, `T=${d.transmittanceLut[0][0]}`);
}

console.log(`\nRESULT ${failed === 0 ? "PASS" : "FAIL"}: ${passed} passed / ${failed} failed`);
if (rows.length) console.log(`assertions=${rows.length}`);
process.exit(failed === 0 ? 0 : 1);
