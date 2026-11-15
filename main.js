// main.js — renderer only. All physics comes from src/core.mjs (the same module
// the node benches assert); this file just samples a texture built from it.
import * as THREE from "three";
import {
  defaultAtmosphere,
  opticalDepth,
  opticalDepthPlaneParallel,
  transmittance,
  rayleighCrossSection,
  SKY_WAVELENGTHS_NM,
} from "./src/core.mjs";

const atmo = defaultAtmosphere();
const R_h = atmo.R / atmo.H;
const OPTS = { panels: 12, nodes: 16 };

// Per-channel optical depth ratios: beta(lambda) / beta(550 nm). The absolute
// scale is an aesthetic exposure, the *ratio* is the certified part.
const beta = SKY_WAVELENGTHS_NM.map((nm) => rayleighCrossSection(nm));
const betaRatio = beta.map((b) => b / beta[1]);
const TEX_W = 512;

function buildSkyTexture(h0, exposure) {
  const data = new Uint8Array(TEX_W * 4);
  for (let i = 0; i < TEX_W; i++) {
    const mu = Math.max(1e-6, Math.cos((1 - i / (TEX_W - 1)) * (Math.PI / 2)));
    const tau = opticalDepth(R_h + h0, mu, R_h, OPTS);
    for (let c = 0; c < 3; c++) {
      const T = Math.exp(-tau * betaRatio[c] * exposure);
      data[i * 4 + c] = Math.max(0, Math.min(255, Math.round(255 * T)));
    }
    data[i * 4 + 3] = 255;
  }
  const tex = new THREE.DataTexture(data, TEX_W, 1, THREE.RGBAFormat);
  tex.magFilter = THREE.LinearFilter;
  tex.minFilter = THREE.LinearFilter;
  tex.wrapS = THREE.ClampToEdgeWrapping;
  tex.needsUpdate = true;
  return tex;
}

const canvas = document.getElementById("sky");
const renderer = new THREE.WebGLRenderer({ canvas, antialias: false });
renderer.setPixelRatio(Math.min(devicePixelRatio, 2));

const scene = new THREE.Scene();
const camera = new THREE.OrthographicCamera(-1, 1, 1, -1, 0, 1);
const material = new THREE.ShaderMaterial({
  uniforms: { sky: { value: buildSkyTexture(0, 6) } },
  vertexShader: "varying vec2 vUv; void main(){ vUv = uv; gl_Position = vec4(position.xy, 0.0, 1.0); }",
  fragmentShader: "uniform sampler2D sky; varying vec2 vUv; void main(){ gl_FragColor = texture2D(sky, vec2(vUv.y, 0.5)); }",
});
scene.add(new THREE.Mesh(new THREE.PlaneGeometry(2, 2), material));

function resize() {
  const w = canvas.clientWidth || 800;
  const h = canvas.clientHeight || 600;
  renderer.setSize(w, h, false);
}
resize();
addEventListener("resize", resize);

const out = document.getElementById("out");
const verdict = document.getElementById("verdict");
const zlab = document.getElementById("zlab");
const hlab = document.getElementById("hlab");
const zIn = document.getElementById("zenith");
const hIn = document.getElementById("height");

function row(k, v) {
  return `<tr><td>${k}</td><td class="num">${v}</td></tr>`;
}

function update() {
  const zenithDeg = Number(zIn.value);
  const h0 = Number(hIn.value);
  const mu = Math.cos((zenithDeg * Math.PI) / 180);
  const tauS = opticalDepth(R_h + h0, mu, R_h, OPTS);
  const tauP = opticalDepthPlaneParallel(R_h + h0, mu, R_h);
  const div = (tauP - tauS) / tauS;

  material.uniforms.sky.value.dispose();
  material.uniforms.sky.value = buildSkyTexture(h0, 6);

  zlab.textContent = zenithDeg.toFixed(1);
  hlab.textContent = h0.toFixed(2);
  out.innerHTML =
    row("tau spherical", tauS.toFixed(6)) +
    row("tau plane-parallel", Number.isFinite(tauP) ? tauP.toFixed(6) : "diverges") +
    row("divergence", `${(div * 100).toFixed(3)} %`) +
    row("T spherical", transmittance(R_h + h0, mu, R_h, OPTS).toExponential(4)) +
    row("channel ratio 440/680", (betaRatio[0] / betaRatio[2]).toFixed(6));
  verdict.innerHTML =
    Math.abs(div) < 0.01
      ? '<span class="num">plane-parallel is safe here (&lt;1%)</span>'
      : '<span class="warn">plane-parallel is off by more than 1% here — this is the measured divergence bench/verify.mjs asserts</span>';
}

zIn.addEventListener("input", update);
hIn.addEventListener("input", update);
update();

function frame() {
  renderer.render(scene, camera);
  requestAnimationFrame(frame);
}
frame();
