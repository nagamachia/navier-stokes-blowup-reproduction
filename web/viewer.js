import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

const root = document.getElementById('viewport');
const status = document.getElementById('status');
const datasetSelect = document.getElementById('dataset');
const modeSelect = document.getElementById('mode');
const strideInput = document.getElementById('stride');
const scaleInput = document.getElementById('scale');
const thresholdInput = document.getElementById('threshold');
const strideValue = document.getElementById('strideValue');
const scaleValue = document.getElementById('scaleValue');
const thresholdValue = document.getElementById('thresholdValue');

const scene = new THREE.Scene();
scene.background = new THREE.Color(0x0b1020);
const camera = new THREE.PerspectiveCamera(45, 1, 0.01, 1000);
camera.position.set(8, 8, 8);
const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
root.appendChild(renderer.domElement);
const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
controls.dampingFactor = 0.08;
scene.add(new THREE.AxesHelper(1.5));

let fieldObject = null;
let parsed = null;
let diagnostics = null;
let needsFrame = true;

function resize() {
  const w = Math.max(root.clientWidth, 1);
  const h = Math.max(root.clientHeight, 1);
  renderer.setSize(w, h, false);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
}
new ResizeObserver(resize).observe(root);

function parseLegacyStructuredPoints(text) {
  const lines = text.split(/\r?\n/);
  let dims, origin = [0, 0, 0], spacing = [1, 1, 1], vectorStart = -1, count = 0;
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i].trim();
    if (line.startsWith('DIMENSIONS ')) dims = line.split(/\s+/).slice(1).map(Number);
    else if (line.startsWith('ORIGIN ')) origin = line.split(/\s+/).slice(1).map(Number);
    else if (line.startsWith('SPACING ')) spacing = line.split(/\s+/).slice(1).map(Number);
    else if (line.startsWith('POINT_DATA ')) count = Number(line.split(/\s+/)[1]);
    else if (line.startsWith('VECTORS ')) { vectorStart = i + 1; break; }
  }
  if (!dims || vectorStart < 0 || !count) throw new Error('対応している VTK STRUCTURED_POINTS / VECTORS を検出できませんでした');
  const nums = lines.slice(vectorStart).join(' ').trim().split(/\s+/).map(Number);
  if (nums.length < count * 3) throw new Error('VTK velocity vector data is incomplete');
  return { dims, origin, spacing, count, vectors: Float64Array.from(nums.slice(0, count * 3)) };
}

function idx(i, j, k, nx, ny) { return ((k * ny + j) * nx + i); }
function wrap(i, n) { return (i + n) % n; }

function computeDiagnostics(data) {
  const { dims: [nx, ny, nz], spacing: [dx, dy, dz], vectors } = data;
  const n = nx * ny * nz;
  const speed = new Float64Array(n);
  const vorticity = new Float64Array(n);
  const qcriterion = new Float64Array(n);
  let maxSpeed = 0, maxVorticity = 0, maxQ = 0;

  const comp = (i, j, k, c) => vectors[idx(wrap(i, nx), wrap(j, ny), wrap(k, nz), nx, ny) * 3 + c];
  for (let k = 0; k < nz; k++) for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
    const p = idx(i, j, k, nx, ny);
    const u = vectors[p * 3], v = vectors[p * 3 + 1], w = vectors[p * 3 + 2];
    speed[p] = Math.hypot(u, v, w);
    maxSpeed = Math.max(maxSpeed, speed[p]);

    const duDx = (comp(i + 1, j, k, 0) - comp(i - 1, j, k, 0)) / (2 * dx);
    const duDy = (comp(i, j + 1, k, 0) - comp(i, j - 1, k, 0)) / (2 * dy);
    const duDz = (comp(i, j, k + 1, 0) - comp(i, j, k - 1, 0)) / (2 * dz);
    const dvDx = (comp(i + 1, j, k, 1) - comp(i - 1, j, k, 1)) / (2 * dx);
    const dvDy = (comp(i, j + 1, k, 1) - comp(i, j - 1, k, 1)) / (2 * dy);
    const dvDz = (comp(i, j, k + 1, 1) - comp(i, j, k - 1, 1)) / (2 * dz);
    const dwDx = (comp(i + 1, j, k, 2) - comp(i - 1, j, k, 2)) / (2 * dx);
    const dwDy = (comp(i, j + 1, k, 2) - comp(i, j - 1, k, 2)) / (2 * dy);
    const dwDz = (comp(i, j, k + 1, 2) - comp(i, j, k - 1, 2)) / (2 * dz);

    const wx = dwDy - dvDz, wy = duDz - dwDx, wz = dvDx - duDy;
    vorticity[p] = Math.hypot(wx, wy, wz);
    maxVorticity = Math.max(maxVorticity, vorticity[p]);

    const s11 = duDx, s22 = dvDy, s33 = dwDz;
    const s12 = 0.5 * (duDy + dvDx), s13 = 0.5 * (duDz + dwDx), s23 = 0.5 * (dvDz + dwDy);
    const o12 = 0.5 * (duDy - dvDx), o13 = 0.5 * (duDz - dwDx), o23 = 0.5 * (dvDz - dwDy);
    const s2 = s11*s11 + s22*s22 + s33*s33 + 2*(s12*s12 + s13*s13 + s23*s23);
    const o2 = 2*(o12*o12 + o13*o13 + o23*o23);
    qcriterion[p] = 0.5 * (o2 - s2);
    maxQ = Math.max(maxQ, qcriterion[p]);
  }
  return { speed, vorticity, qcriterion, maxSpeed, maxVorticity, maxQ };
}

function colorForValue(value, maxValue) {
  const c = new THREE.Color();
  const t = maxValue > 0 ? Math.min(Math.max(value / maxValue, 0), 1) : 0;
  c.setHSL((1 - t) * 0.66, 0.9, 0.58);
  return c;
}

function disposeField() {
  if (!fieldObject) return;
  scene.remove(fieldObject);
  fieldObject.geometry.dispose();
  fieldObject.material.dispose();
  fieldObject = null;
}

function frameGeometry(geom) {
  geom.computeBoundingBox();
  const box = geom.boundingBox;
  if (!box || box.isEmpty()) return;
  const center = new THREE.Vector3(); box.getCenter(center);
  const size = new THREE.Vector3(); box.getSize(size);
  const radius = Math.max(size.x, size.y, size.z, 1);
  controls.target.copy(center);
  camera.near = Math.max(radius / 1000, 0.001);
  camera.far = radius * 50;
  camera.updateProjectionMatrix();
  camera.position.copy(center).add(new THREE.Vector3(radius * 1.3, radius * 1.3, radius * 1.3));
  controls.update();
}

function percentileThreshold(values, keepPercent, positiveOnly = false) {
  const samples = [];
  for (let i = 0; i < values.length; i++) if (!positiveOnly || values[i] > 0) samples.push(values[i]);
  if (!samples.length) return Infinity;
  samples.sort((a, b) => a - b);
  const q = Math.max(0, Math.min(1, 1 - keepPercent / 100));
  return samples[Math.floor(q * (samples.length - 1))];
}

function rebuild() {
  if (!parsed || !diagnostics) return;
  disposeField();
  const stride = Number(strideInput.value);
  const scale = Number(scaleInput.value);
  const keepPercent = Number(thresholdInput.value);
  strideValue.textContent = String(stride);
  scaleValue.textContent = scale.toFixed(2);
  thresholdValue.textContent = String(keepPercent);
  const { dims: [nx, ny, nz], origin, spacing, vectors } = parsed;
  const mode = modeSelect.value;
  let geom;

  if (mode === 'velocity') {
    const pos = [], cols = [];
    let shown = 0;
    for (let k = 0; k < nz; k += stride) for (let j = 0; j < ny; j += stride) for (let i = 0; i < nx; i += stride) {
      const p0 = idx(i, j, k, nx, ny), p = p0 * 3;
      const x = origin[0] + i * spacing[0], y = origin[1] + j * spacing[1], z = origin[2] + k * spacing[2];
      const u = vectors[p], v = vectors[p + 1], w = vectors[p + 2];
      const c = colorForValue(diagnostics.speed[p0], diagnostics.maxSpeed);
      pos.push(x, y, z, x + u * scale, y + v * scale, z + w * scale);
      cols.push(c.r, c.g, c.b, c.r, c.g, c.b);
      shown++;
    }
    geom = new THREE.BufferGeometry();
    geom.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
    geom.setAttribute('color', new THREE.Float32BufferAttribute(cols, 3));
    fieldObject = new THREE.LineSegments(geom, new THREE.LineBasicMaterial({ vertexColors: true }));
    status.textContent = `${nx}×${ny}×${nz} / ${shown.toLocaleString()} vectors / max |u|=${diagnostics.maxSpeed.toFixed(4)} / max |ω|=${diagnostics.maxVorticity.toFixed(4)} / max Q=${diagnostics.maxQ.toFixed(4)}`;
  } else {
    const values = mode === 'vorticity' ? diagnostics.vorticity : diagnostics.qcriterion;
    const maxValue = mode === 'vorticity' ? diagnostics.maxVorticity : diagnostics.maxQ;
    const positiveOnly = mode === 'qcriterion';
    const threshold = percentileThreshold(values, keepPercent, positiveOnly);
    const pos = [], cols = [];
    let shown = 0;
    for (let k = 0; k < nz; k += stride) for (let j = 0; j < ny; j += stride) for (let i = 0; i < nx; i += stride) {
      const p = idx(i, j, k, nx, ny), value = values[p];
      if (value < threshold || (positiveOnly && value <= 0)) continue;
      const x = origin[0] + i * spacing[0], y = origin[1] + j * spacing[1], z = origin[2] + k * spacing[2];
      const c = colorForValue(value, maxValue);
      pos.push(x, y, z); cols.push(c.r, c.g, c.b); shown++;
    }
    geom = new THREE.BufferGeometry();
    geom.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
    geom.setAttribute('color', new THREE.Float32BufferAttribute(cols, 3));
    const pointSize = Math.max(...spacing) * 2.2;
    fieldObject = new THREE.Points(geom, new THREE.PointsMaterial({ vertexColors: true, size: pointSize, sizeAttenuation: true }));
    const label = mode === 'vorticity' ? '|ω|' : 'Q';
    status.textContent = `${nx}×${ny}×${nz} / ${shown.toLocaleString()} points / strongest ${keepPercent}% / threshold ${label}=${Number.isFinite(threshold) ? threshold.toFixed(4) : 'n/a'} / max ${label}=${maxValue.toFixed(4)}`;
  }

  scene.add(fieldObject);
  if (needsFrame) { frameGeometry(geom); needsFrame = false; }
}

async function loadDataset() {
  status.textContent = 'VTK を読み込み中…';
  try {
    const res = await fetch(datasetSelect.value, { cache: 'no-store' });
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    parsed = parseLegacyStructuredPoints(await res.text());
    status.textContent = '速度勾配・渦度・Q-criterion を計算中…';
    diagnostics = computeDiagnostics(parsed);
    needsFrame = true;
    rebuild();
  } catch (e) {
    console.error(e);
    status.textContent = `読み込み失敗: ${e.message}`;
  }
}

datasetSelect.addEventListener('change', loadDataset);
modeSelect.addEventListener('change', rebuild);
strideInput.addEventListener('input', rebuild);
scaleInput.addEventListener('input', rebuild);
thresholdInput.addEventListener('input', rebuild);

function animate() { requestAnimationFrame(animate); controls.update(); renderer.render(scene, camera); }
resize(); loadDataset(); animate();
