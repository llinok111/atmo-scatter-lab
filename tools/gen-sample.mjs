// tools/gen-sample.mjs — deterministic sample data. Same command twice must
// produce byte-identical files (verified in bench/verify.mjs).
//
//   node tools/gen-sample.mjs
//
// Writes data/atmosphere.json: the parameter box plus a transmittance LUT that
// both the node benches and index.html/main.js consume.

import { mkdirSync, writeFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";
import {
  defaultAtmosphere,
  opticalDepth,
  opticalDepthPlaneParallel,
  horizontalLimitHalfGaussian,
  transmittance,
  rayleighCrossSection,
  SKY_WAVELENGTHS_NM,
} from "../src/core.mjs";

const here = dirname(fileURLToPath(import.meta.url));
const dataDir = join(here, "..", "data");
mkdirSync(dataDir, { recursive: true });

const atmo = defaultAtmosphere();
const R_h = atmo.R / atmo.H;

// LUT: 91 zenith angles (0..90 deg) x 25 heights (0..24 H)
const zenithDeg = Array.from({ length: 91 }, (_, i) => i);
const heightsH = Array.from({ length: 25 }, (_, i) => i);
const tauLut = zenithDeg.map((deg) => {
  const mu = Math.cos((deg * Math.PI) / 180);
  return heightsH.map((h) =>
    mu <= 0 ? opticalDepth(R_h + h, 0, R_h, { panels: 12, nodes: 16 }) : opticalDepth(R_h + h, mu, R_h, { panels: 12, nodes: 16 }),
  );
});

const divergenceDeg = [0, 15, 30, 45, 60, 75, 85, 89];
const divergence = divergenceDeg.map((deg) => {
  const mu = Math.cos((deg * Math.PI) / 180);
  const spherical = opticalDepth(R_h, mu, R_h, { panels: 12, nodes: 16 });
  const plane = opticalDepthPlaneParallel(R_h, mu, R_h);
  return { zenithDeg: deg, tauSpherical: spherical, tauPlaneParallel: plane, relDivergence: (plane - spherical) / spherical };
});

const sigma = SKY_WAVELENGTHS_NM.map((nm) => ({ lambdaNm: nm, sigma: rayleighCrossSection(nm) }));

const payload = {
  generator: "tools/gen-sample.mjs",
  atmosphere: { R_km: atmo.R, H_km: atmo.H, name: atmo.name, R_over_H: R_h },
  grid: { zenithDeg, heightsH },
  tauLut,
  transmittanceLut: tauLut.map((row) => row.map((t) => Math.exp(-t))),
  horizontal: {
    tauMeasured: opticalDepth(R_h, 0, R_h, { panels: 12, nodes: 16 }),
    tauHalfGaussianLimit: horizontalLimitHalfGaussian(R_h),
  },
  divergence,
  wavelengths: sigma,
  verticalColumn: { tauAtZenith: transmittance(R_h, 1, R_h, { panels: 12, nodes: 16 }) },
};

const out = join(dataDir, "atmosphere.json");
writeFileSync(out, JSON.stringify(payload, null, 1) + "\n");
console.log(`wrote ${out} (${JSON.stringify(payload).length} bytes)`);
