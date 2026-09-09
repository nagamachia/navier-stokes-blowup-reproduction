import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

const root=document.getElementById('viewport'),status=document.getElementById('status');
const datasetSelect=document.getElementById('dataset'),modeSelect=document.getElementById('mode');
const strideInput=document.getElementById('stride'),scaleInput=document.getElementById('scale'),thresholdInput=document.getElementById('threshold');
const strideValue=document.getElementById('strideValue'),scaleValue=document.getElementById('scaleValue'),thresholdValue=document.getElementById('thresholdValue');
const timeInput=document.getElementById('time'),timeValue=document.getElementById('timeValue'),playButton=document.getElementById('play');

const scene=new THREE.Scene(); scene.background=new THREE.Color(0x0b1020);
const camera=new THREE.PerspectiveCamera(45,1,.01,1000); camera.position.set(8,8,8);
const renderer=new THREE.WebGLRenderer({antialias:true}); renderer.setPixelRatio(Math.min(window.devicePixelRatio||1,2)); root.appendChild(renderer.domElement);
const controls=new OrbitControls(camera,renderer.domElement); controls.enableDamping=true; controls.dampingFactor=.08; scene.add(new THREE.AxesHelper(1.5));
let fieldObject=null,parsed=null,diagnostics=null,needsFrame=true,series=null,seriesBase='',playing=false,playTimer=null,requestId=0;
const frameCache=new Map();

function resize(){const w=Math.max(root.clientWidth,1),h=Math.max(root.clientHeight,1);renderer.setSize(w,h,false);camera.aspect=w/h;camera.updateProjectionMatrix();}
new ResizeObserver(resize).observe(root);

function parseLegacyStructuredPoints(text){
  const lines=text.split(/\r?\n/);let dims,origin=[0,0,0],spacing=[1,1,1],vectorStart=-1,count=0;
  for(let i=0;i<lines.length;i++){const line=lines[i].trim();if(line.startsWith('DIMENSIONS '))dims=line.split(/\s+/).slice(1).map(Number);else if(line.startsWith('ORIGIN '))origin=line.split(/\s+/).slice(1).map(Number);else if(line.startsWith('SPACING '))spacing=line.split(/\s+/).slice(1).map(Number);else if(line.startsWith('POINT_DATA '))count=Number(line.split(/\s+/)[1]);else if(line.startsWith('VECTORS ')){vectorStart=i+1;break;}}
  if(!dims||vectorStart<0||!count)throw new Error('対応している VTK STRUCTURED_POINTS / VECTORS を検出できませんでした');
  const tokens=lines.slice(vectorStart).join(' ').trim().split(/\s+/);const vectors=new Float64Array(count*3);for(let i=0;i<vectors.length;i++){vectors[i]=Number(tokens[i]);if(!Number.isFinite(vectors[i]))throw new Error('VTK velocity vector data is incomplete');}
  return {dims,origin,spacing,count,vectors};
}
const idx=(i,j,k,nx,ny)=>(k*ny+j)*nx+i,wrap=(i,n)=>(i+n)%n;

function computeDiagnostics(data){
  const {dims:[nx,ny,nz],spacing:[dx,dy,dz],vectors}=data,n=nx*ny*nz;const speed=new Float64Array(n),vorticity=new Float64Array(n),qcriterion=new Float64Array(n);let maxSpeed=0,maxVorticity=0,maxQ=0;
  const comp=(i,j,k,c)=>vectors[idx(wrap(i,nx),wrap(j,ny),wrap(k,nz),nx,ny)*3+c];
  for(let k=0;k<nz;k++)for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){
    const p=idx(i,j,k,nx,ny),u=vectors[p*3],v=vectors[p*3+1],w=vectors[p*3+2];speed[p]=Math.hypot(u,v,w);maxSpeed=Math.max(maxSpeed,speed[p]);
    const duDx=(comp(i+1,j,k,0)-comp(i-1,j,k,0))/(2*dx),duDy=(comp(i,j+1,k,0)-comp(i,j-1,k,0))/(2*dy),duDz=(comp(i,j,k+1,0)-comp(i,j,k-1,0))/(2*dz);
    const dvDx=(comp(i+1,j,k,1)-comp(i-1,j,k,1))/(2*dx),dvDy=(comp(i,j+1,k,1)-comp(i,j-1,k,1))/(2*dy),dvDz=(comp(i,j,k+1,1)-comp(i,j,k-1,1))/(2*dz);
    const dwDx=(comp(i+1,j,k,2)-comp(i-1,j,k,2))/(2*dx),dwDy=(comp(i,j+1,k,2)-comp(i,j-1,k,2))/(2*dy),dwDz=(comp(i,j,k+1,2)-comp(i,j,k-1,2))/(2*dz);
    const wx=dwDy-dvDz,wy=duDz-dwDx,wz=dvDx-duDy;vorticity[p]=Math.hypot(wx,wy,wz);maxVorticity=Math.max(maxVorticity,vorticity[p]);
    const s11=duDx,s22=dvDy,s33=dwDz,s12=.5*(duDy+dvDx),s13=.5*(duDz+dwDx),s23=.5*(dvDz+dwDy),o12=.5*(duDy-dvDx),o13=.5*(duDz-dwDx),o23=.5*(dvDz-dwDy);
    const s2=s11*s11+s22*s22+s33*s33+2*(s12*s12+s13*s13+s23*s23),o2=2*(o12*o12+o13*o13+o23*o23);qcriterion[p]=.5*(o2-s2);maxQ=Math.max(maxQ,qcriterion[p]);
  }
  return {speed,vorticity,qcriterion,maxSpeed,maxVorticity,maxQ};
}
function colorForValue(value,maxValue){const c=new THREE.Color(),t=maxValue>0?Math.min(Math.max(value/maxValue,0),1):0;c.setHSL((1-t)*.66,.9,.58);return c;}
function disposeField(){if(!fieldObject)return;scene.remove(fieldObject);fieldObject.geometry.dispose();fieldObject.material.dispose();fieldObject=null;}
function frameGeometry(geom){geom.computeBoundingBox();const box=geom.boundingBox;if(!box||box.isEmpty())return;const center=new THREE.Vector3(),size=new THREE.Vector3();box.getCenter(center);box.getSize(size);const radius=Math.max(size.x,size.y,size.z,1);controls.target.copy(center);camera.near=Math.max(radius/1000,.001);camera.far=radius*50;camera.updateProjectionMatrix();camera.position.copy(center).add(new THREE.Vector3(radius*1.3,radius*1.3,radius*1.3));controls.update();}
function percentileThreshold(values,keepPercent,positiveOnly=false){const samples=[];for(let i=0;i<values.length;i++)if(!positiveOnly||values[i]>0)samples.push(values[i]);if(!samples.length)return Infinity;samples.sort((a,b)=>a-b);return samples[Math.floor(Math.max(0,Math.min(1,1-keepPercent/100))*(samples.length-1))];}

function rebuild(){
  if(!parsed||!diagnostics)return;disposeField();const stride=Number(strideInput.value),scale=Number(scaleInput.value),keepPercent=Number(thresholdInput.value);strideValue.textContent=String(stride);scaleValue.textContent=scale.toFixed(2);thresholdValue.textContent=String(keepPercent);
  const {dims:[nx,ny,nz],origin,spacing,vectors}=parsed,mode=modeSelect.value;let geom;
  if(mode==='velocity'){
    const pos=[],cols=[];let shown=0;for(let k=0;k<nz;k+=stride)for(let j=0;j<ny;j+=stride)for(let i=0;i<nx;i+=stride){const p0=idx(i,j,k,nx,ny),p=p0*3,x=origin[0]+i*spacing[0],y=origin[1]+j*spacing[1],z=origin[2]+k*spacing[2],u=vectors[p],v=vectors[p+1],w=vectors[p+2],c=colorForValue(diagnostics.speed[p0],diagnostics.maxSpeed);pos.push(x,y,z,x+u*scale,y+v*scale,z+w*scale);cols.push(c.r,c.g,c.b,c.r,c.g,c.b);shown++;}
    geom=new THREE.BufferGeometry();geom.setAttribute('position',new THREE.Float32BufferAttribute(pos,3));geom.setAttribute('color',new THREE.Float32BufferAttribute(cols,3));fieldObject=new THREE.LineSegments(geom,new THREE.LineBasicMaterial({vertexColors:true}));status.textContent=`${nx}×${ny}×${nz} / ${shown.toLocaleString()} vectors / max |u|=${diagnostics.maxSpeed.toFixed(4)} / max |ω|=${diagnostics.maxVorticity.toFixed(4)} / max Q=${diagnostics.maxQ.toFixed(4)}`;
  }else{
    const values=mode==='vorticity'?diagnostics.vorticity:diagnostics.qcriterion,maxValue=mode==='vorticity'?diagnostics.maxVorticity:diagnostics.maxQ,positiveOnly=mode==='qcriterion',threshold=percentileThreshold(values,keepPercent,positiveOnly),pos=[],cols=[];let shown=0;
    for(let k=0;k<nz;k+=stride)for(let j=0;j<ny;j+=stride)for(let i=0;i<nx;i+=stride){const p=idx(i,j,k,nx,ny),value=values[p];if(value<threshold||(positiveOnly&&value<=0))continue;const x=origin[0]+i*spacing[0],y=origin[1]+j*spacing[1],z=origin[2]+k*spacing[2],c=colorForValue(value,maxValue);pos.push(x,y,z);cols.push(c.r,c.g,c.b);shown++;}
    geom=new THREE.BufferGeometry();geom.setAttribute('position',new THREE.Float32BufferAttribute(pos,3));geom.setAttribute('color',new THREE.Float32BufferAttribute(cols,3));fieldObject=new THREE.Points(geom,new THREE.PointsMaterial({vertexColors:true,size:Math.max(...spacing)*2.2,sizeAttenuation:true}));const label=mode==='vorticity'?'|ω|':'Q';status.textContent=`${nx}×${ny}×${nz} / ${shown.toLocaleString()} points / strongest ${keepPercent}% / threshold ${label}=${Number.isFinite(threshold)?threshold.toFixed(4):'n/a'} / max ${label}=${maxValue.toFixed(4)}`;
  }
  scene.add(fieldObject);if(needsFrame){frameGeometry(geom);needsFrame=false;}
}

async function fetchVtk(url){if(frameCache.has(url))return frameCache.get(url);const res=await fetch(url,{cache:'no-store'});if(!res.ok)throw new Error(`HTTP ${res.status}`);const data=parseLegacyStructuredPoints(await res.text());if(frameCache.size>=3)frameCache.delete(frameCache.keys().next().value);frameCache.set(url,data);return data;}
async function loadVtk(url,frameCamera=false){const mine=++requestId;status.textContent='VTK を読み込み中…';const data=await fetchVtk(url);if(mine!==requestId)return;parsed=data;status.textContent='速度勾配・渦度・Q-criterion を計算中…';diagnostics=computeDiagnostics(parsed);if(frameCamera)needsFrame=true;rebuild();}
function stopPlayback(){playing=false;if(playTimer){clearTimeout(playTimer);playTimer=null;}playButton.textContent='▶ 再生';}
async function loadSeriesFrame(index,frameCamera=false){if(!series)return;const frame=series.frames[index];timeValue.textContent=`t = ${Number(frame.time).toFixed(3)} (${index+1}/${series.frames.length})`;timeInput.value=String(index);try{await loadVtk(seriesBase+frame.file,frameCamera);}catch(e){console.error(e);status.textContent=`読み込み失敗: ${e.message}`;stopPlayback();}}
function scheduleNext(){if(!playing||!series)return;playTimer=setTimeout(async()=>{let next=(Number(timeInput.value)+1)%series.frames.length;await loadSeriesFrame(next,false);scheduleNext();},650);}

async function loadDataset(){stopPlayback();series=null;timeInput.disabled=true;playButton.disabled=true;timeValue.textContent='—';frameCache.clear();needsFrame=true;const value=datasetSelect.value;
  try{
    if(value.endsWith('.json')){status.textContent='3D時系列 manifest を読み込み中…';const res=await fetch(value,{cache:'no-store'});if(!res.ok)throw new Error(`HTTP ${res.status}`);series=await res.json();if(!Array.isArray(series.frames)||!series.frames.length)throw new Error('時系列フレームがありません');seriesBase=value.slice(0,value.lastIndexOf('/')+1);timeInput.min='0';timeInput.max=String(series.frames.length-1);timeInput.value='0';timeInput.disabled=false;playButton.disabled=false;await loadSeriesFrame(0,true);}
    else await loadVtk(value,true);
  }catch(e){console.error(e);status.textContent=`読み込み失敗: ${e.message}`;}
}

datasetSelect.addEventListener('change',loadDataset);modeSelect.addEventListener('change',rebuild);strideInput.addEventListener('input',rebuild);scaleInput.addEventListener('input',rebuild);thresholdInput.addEventListener('input',rebuild);
timeInput.addEventListener('input',()=>{stopPlayback();loadSeriesFrame(Number(timeInput.value),false);});
playButton.addEventListener('click',()=>{if(playing){stopPlayback();return;}if(!series)return;playing=true;playButton.textContent='⏸ 停止';scheduleNext();});
function animate(){requestAnimationFrame(animate);controls.update();renderer.render(scene,camera);}resize();loadDataset();animate();
