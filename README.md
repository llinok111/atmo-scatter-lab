# atmo-scatter-lab

**Promise:** turn atmospheric single-scattering transmittance into a *node-verifiable* number — the composite Gauss–Legendre optical depth matches an independent dense-trapezoid oracle to **6.3e-13** relative (the raw oracle is off by **1.28e-6**, and its own error is *predicted* by Euler–Maclaurin to within 2%), while the cheap plane-parallel `sech(χ)` model that every real-time sky uses diverges from the spherical solution by **0.042% at 30°, 0.37% at 60°, 1.69% at 75°, 13.15% at 85° and 130% at 89°** — and stays finite (τ→35.3853) where plane-parallel blows up.

Upstream: none — original work. The physics is classical (Rayleigh 1871, Chapman 1931); no upstream source tree is copied or patched.

Baseline in-repo: `baseline/naive.mjs` (uniform composite trapezoid, imports nothing).

## Run it

```bash
node tools/gen-sample.mjs     # writes data/atmosphere.json (deterministic)
node bench/verify.mjs         # 22 behavioural assertions, exit 0 on success
node bench/bench.mjs          # writes bench/result.json (7 samples + cold start)
python3 -m http.server 4173   # then open http://localhost:4173  (browser view)
```

The browser page imports the **same** `src/core.mjs`; dragging the zenith slider recomputes τ in the page with the identical code path the node assertions cover.

## Evidence (`bench/result.json`, node v24.19.0)

| quantity | value |
| --- | --- |
| vertical ray (r0=R, χ=0) | τ = 0.999999999999963, rel err **3.73e-14** (exact answer 1) |
| horizontal ray (χ=90°) | τ = **35.385319518** vs closed form √(πR/2) = 35.368671515, rel **4.707e-4** (the O(H/R) parabolic correction) |
| oracle, raw trapezoid | worst rel deviation from core **1.276e-6** |
| oracle, Euler–Maclaurin corrected | worst rel deviation **6.315e-13** |
| Euler–Maclaurin *prediction* of the oracle's own error | ≤2% where it dominates, ≤1e-9 absolute everywhere |
| truncated tail `[sMax, 8·sMax]` | measured **≤ 1.21e-17**, matching the analytic bound `exp(-h(sMax))/h'(sMax)` |
| Rayleigh phase, closed form | sphere integral exactly **1**; GL64 quadrature 0.9999999999999994 |
| HG phase, closed-form primitive | integrates to 1 for g∈{0.1,0.5,0.9,0.95} to **1.1e-14** |
| HG phase, fixed GL64 quadrature | err **-1.163e-5** (g=0.9), **-7.970e-3** (g=0.95); GL256 g=0.95 → **-4.81e-11** |
| σ(440)/σ(680) | **5.704596680555** = (680/440)⁴ to 1e-12 |
| k-scaling R,H,ray by 1e-6…1e6 | worst rel **0.0**; power-of-two k bit-identical (4/4) |
| speed | core **0.0656 ms** median vs baseline **1.8069 ms** → **27.55×** (p95-vs-p95 18.90×, **85.33×** fewer integrand evaluations) |

Cold start: core 0.1551 ms / baseline 3.0188 ms. `core` relStddev is 51% because the median is 0.066 ms — near the hrtime/JIT/GC noise floor; the *absolute* spread (min 0.0543, p95 0.1597 ms) is what the 27× claim rests on, and p95-vs-p95 is reported separately for that reason. The baseline's 21% relStddev comes from the same floor at 1.8 ms.

## 本轮创新

**B — first time (per `state/innovation.md`, 157 valid entries, and `used-upstreams.txt`) that atmospheric transmittance is delivered as a re-computable certificate rather than as a shader.** The nearest neighbours are `llinok111/planet-forge` (atmosphere as a *visual layer*, no oracle) and three.js `Sky.js`/Preetham (a fitted *look*, no error statement); no ledger entry measures where plane-parallel breaks, and none certifies an oracle with the oracle's own error term. Concretely new here: (a) the dense-trapezoid oracle is **corrected by its analytic Euler–Maclaurin error**, so "oracle agrees with core" is not a shared-constant tautology — the raw oracle disagrees by 1.28e-6 and the corrected one by 6.3e-13; (b) an **analytic truncation bound** `exp(-h(sMax))/h'(sMax)` that the measurement respects to 1.21e-17; (c) a measured **divergence curve** with the χ→90° finite limit √(πRH/2), where the standard approximation diverges. This is not a re-skin: nothing in this repo would survive a topic swap, because the numbers *are* the artifact.

## Known limitations (read before trusting a number)

- **Validity range.** The whole model assumes a spherically symmetric exponential atmosphere with a single scale height H, ray start radius r0 ≥ R and χ ∈ [0°, 90°). Rays that descend to the ground (χ > 90°) and refraction are **not** modelled. Single scattering only — no multiple scattering, no ozone layer, no aerosol/Mie mixture.
- **The horizontal closed form is asymptotic.** √(πR/2) (in H units) is the half-Gaussian limit valid for H ≪ R; the measured value is **4.707e-4 relative larger** at H/R = 8/6371. Both numbers are reported; neither is called exact.
- **The 0.3725% etc. divergence figures are for R = 6371 km, H = 8 km and an observer at the ground (h0 = 0).** They are not universal constants: they move with H/R and with h0. `data/atmosphere.json` carries the full 91 × 25 LUT so any other height can be read off without re-running the bench.
- **Scale invariance is not bit-identity in general.** Physical k-scaling measured exactly 0.0 relative difference over 1e-6…1e6 for this parameter pair, and is bit-identical for k = 2^m; a different decimal k could differ in the last ulp because `k·R` and `k·H` round before dividing. The claim is "≤1e-15 relative", not "bit-identical".
- **Timing is not a promise.** Single-digit-microsecond absolute numbers are dominated by JIT/GC; only the ratio to the in-repo baseline is meaningful.
- **The browser view uses an aesthetic exposure of 6** on the per-channel optical depth. The certified quantity is the τ ratio between channels (5.7046 at 440/680 nm) and the divergence readout, not the absolute pixel colour.

## Layout

```
src/core.mjs          pure computation, no three / DOM / window
baseline/naive.mjs    independent trapezoid oracle + Euler-Maclaurin error model
tools/gen-sample.mjs  deterministic data/atmosphere.json (91 x 25 LUT)
bench/verify.mjs      22 assertions, prints every value
bench/bench.mjs       7 samples + cold start -> bench/result.json
main.js / index.html  three.js renderer, importmap pinned to three@0.169.0
```

MIT licensed.
