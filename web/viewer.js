import * as THREE from 'https://cdn.jsdelivr.net/npm/three@0.180.0/build/three.module.js';
import { OrbitControls } from 'https://cdn.jsdelivr.net/npm/three@0.180.0/examples/jsm/controls/OrbitControls.js';

const root = document.getElementById('viewport');
const status = document.getElementById('status');
const datasetSelect = document.getElementById('dataset');
const strideInput = document.getElementById('stride');
const scaleInput = document.getElementById('scale');
const strideValue = document.getElementById('strideValue');
const scaleValue = document.getElementById('scaleValue');

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
  let dims, origin = [0,0,0], spacing = [1,1,1], vectorStart = -1, count = 0;
  for (let i=0; i<lines.length; i++) {
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
  return { dims, origin, spacing, count, vectors: nums.slice(0, count*3) };
}

function colorForMagnitude(m, maxM) {
  const c = new THREE.Color();
  const t = maxM > 0 ? Math.min(m / maxM, 1) : 0;
  c.setHSL((1 - t) * 0.66, 0.9, 0.58);
  return c;
}

function rebuild() {
  if (!parsed) return;
  if (fieldObject) { scene.remove(fieldObject); fieldObject.geometry.dispose(); fieldObject.material.dispose(); }
  const stride = Number(strideInput.value);
  const scale = Number(scaleInput.value);
  strideValue.textContent = String(stride);
  scaleValue.textContent = scale.toFixed(2);
  const {dims:[nx,ny,nz], origin, spacing, vectors} = parsed;
  let maxM = 0;
  for (let i=0;i<vectors.length;i+=3) maxM = Math.max(maxM, Math.hypot(vectors[i],vectors[i+1],vectors[i+2]));
  const pos = [], cols = [];
  let shown = 0;
  for (let k=0;k<nz;k+=stride) for (let j=0;j<ny;j+=stride) for (let i=0;i<nx;i+=stride) {
    const p = ((k*ny + j)*nx + i)*3;
    const x = origin[0] + i*spacing[0], y = origin[1] + j*spacing[1], z = origin[2] + k*spacing[2];
    const u=vectors[p], v=vectors[p+1], w=vectors[p+2];
    const m=Math.hypot(u,v,w), c=colorForMagnitude(m,maxM);
    pos.push(x,y,z, x+u*scale,y+v*scale,z+w*scale);
    cols.push(c.r,c.g,c.b,c.r,c.g,c.b);
    shown++;
  }
  const geom = new THREE.BufferGeometry();
  geom.setAttribute('position', new THREE.Float32BufferAttribute(pos,3));
  geom.setAttribute('color', new THREE.Float32BufferAttribute(cols,3));
  fieldObject = new THREE.LineSegments(geom, new THREE.LineBasicMaterial({vertexColors:true}));
  scene.add(fieldObject);
  geom.computeBoundingBox();
  const box=geom.boundingBox, center=new THREE.Vector3(); box.getCenter(center);
  controls.target.copy(center);
  const size=new THREE.Vector3(); box.getSize(size);
  const radius=Math.max(size.x,size.y,size.z,1);
  camera.near=Math.max(radius/1000,0.001); camera.far=radius*50; camera.updateProjectionMatrix();
  camera.position.copy(center).add(new THREE.Vector3(radius*1.3,radius*1.3,radius*1.3));
  controls.update();
  status.textContent=`${nx}×${ny}×${nz} grid / ${shown.toLocaleString()} vectors shown / max |u|=${maxM.toFixed(4)}`;
}

async function loadDataset() {
  status.textContent='VTK を読み込み中…';
  try {
    const res = await fetch(datasetSelect.value, {cache:'no-store'});
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    parsed = parseLegacyStructuredPoints(await res.text());
    rebuild();
  } catch (e) {
    console.error(e);
    status.textContent=`読み込み失敗: ${e.message}`;
  }
}

datasetSelect.addEventListener('change', loadDataset);
strideInput.addEventListener('input', rebuild);
scaleInput.addEventListener('input', rebuild);

function animate(){ requestAnimationFrame(animate); controls.update(); renderer.render(scene,camera); }
resize(); loadDataset(); animate();
