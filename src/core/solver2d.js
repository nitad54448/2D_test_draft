(function(TE){
TE.signal2D=(s,t,f)=>{if(typeof s==='number')return TE.finite(s,'Boundary');TE.assert(s&&typeof s==='object','Invalid boundary signal.');return TE.finite((s.bias??0)+(s.amplitude??0)*Math.cos(2*Math.PI*f*t+(s.phase??0)*Math.PI/180),'Boundary');};
class Solver2D {
 constructor(mesh,materials,boundaries,electrical){this.mesh=mesh;this.materials=materials;this.boundaries=boundaries;this.electrical=electrical;this.frequency=0;
  mesh.map.forEach(id=>TE.assert(Number.isInteger(id)&&materials[id],'Unknown material in map.'));
  TE.assert(['voltage','current','open_circuit'].includes(electrical.kind),'Unknown electrical mode.');
  this.source=mesh.electrode(electrical.sourceSide,electrical.sourceRange);this.sink=mesh.electrode(electrical.sinkSide,electrical.sinkRange);
  const sourceSet=new Set(this.source);TE.assert(!this.sink.some(i=>sourceSet.has(i)),'Electrodes overlap at a node.');
  this.sourceSet=sourceSet;
  for(const side of ['left','right','bottom','top']){const b=boundaries[side];TE.assert(b&&['temperature','flux','convection'].includes(b.kind),'Define all four thermal boundaries.');TE.assert(Number.isFinite(b.h??0)&&(b.h??0)>=0,'Convection coefficient must be nonnegative.');}
 }
 properties(T){
  const m=this.mesh,cap=Array(m.n).fill(0);
  m.cells.forEach(c=>c.nodes.forEach(i=>{const a=this.materials[c.m];cap[i]+=a.density(T[i])*a.heatCapacity(T[i])*c.volume/4;}));
  const p=m.links.map(e=>{const a=this.materials[e.m],Ta=T[e.a],Tb=T[e.b],rho=(a.electricalResistivity(Ta)+a.electricalResistivity(Tb))/2;
   return {g:e.A/(e.L*rho),k:e.A/e.L*(a.thermalConductivity(Ta)+a.thermalConductivity(Tb))/2,alpha:(a.seebeck(Ta)+a.seebeck(Tb))/2};});return {p,cap};
 }
 thermal(t){
  const m=this.mesh,fixed=new Map(),diag=Array(m.n).fill(0),rhs=Array(m.n).fill(0),flux=Array(m.n).fill(0);
  for(const [side,faces] of Object.entries(m.sides)){const b=this.boundaries[side],v=TE.signal2D(b.value,t,this.frequency);
   for(const {node,A} of faces){if(b.kind==='temperature'){TE.assert(v>0,'Temperature must be >0 K.');if(fixed.has(node))TE.assert(Math.abs(fixed.get(node)-v)<1e-8,'Conflicting temperatures at a corner. Use compatible boundary temperatures.');fixed.set(node,v);}
    if(b.kind==='flux'){rhs[node]-=v*A;flux[node]+=v*A;}
    if(b.kind==='convection'){diag[node]+=b.h*A;rhs[node]+=b.h*A*v;}
   }
  }return {fixed,diag,rhs,flux};
 }
 electric(T,p,t){
  const m=this.mesh,e=this.electrical,g=p.map(v=>v.g),zero=Array(m.n).fill(0),rhs=[...zero],fixed=new Map();
  this.source.forEach(i=>fixed.set(i,0));this.sink.forEach(i=>fixed.set(i,0));
  m.links.forEach((l,k)=>{const s=p[k].g*p[k].alpha*(T[l.b]-T[l.a]);rhs[l.a]+=s;rhs[l.b]-=s;});
  const current=(V,seebeck=true)=>m.links.map((l,k)=>g[k]*(V[l.a]-V[l.b]-(seebeck?p[k].alpha*(T[l.b]-T[l.a]):0)));
  const terminal=I=>m.links.reduce((s,l,k)=>s+I[k]*((this.sourceSet.has(l.a)?1:0)-(this.sourceSet.has(l.b)?1:0)),0);
  let V,terminalVoltage;
  if(e.kind==='voltage'){
   terminalVoltage=TE.signal2D(e.value,t,this.frequency);this.sink.forEach(i=>fixed.set(i,terminalVoltage));V=TE.graphSolve(m,g,zero,rhs,fixed);
  }else{
   const base=TE.graphSolve(m,g,zero,rhs,fixed);this.sink.forEach(i=>fixed.set(i,1));
   const unit=TE.graphSolve(m,g,zero,zero,fixed),Ibase=terminal(current(base)),Iunit=terminal(current(unit,false));
   TE.assert(Iunit<0,'Electrodes have no conducting connection.');
   const target=e.kind==='open_circuit'?0:TE.signal2D(e.value,t,this.frequency);
   terminalVoltage=(target-Ibase)/Iunit;V=base.map((v,i)=>v+terminalVoltage*unit[i]);
  }
  const I=current(V);return {V,I,terminalVoltage,current:terminal(I)};
 }
 balance(T,t){
  const m=this.mesh,{p,cap}=this.properties(T),elect=this.electric(T,p,t),source=Array(m.n).fill(0),q=[],work=[];
  m.links.forEach((l,k)=>{const I=elect.I[k],P=I*(elect.V[l.a]-elect.V[l.b]),pel=p[k].alpha*(T[l.a]+T[l.b])/2*I;
   source[l.a]+=-pel+P/2;source[l.b]+=pel+P/2;q.push(pel-p[k].k*(T[l.b]-T[l.a]));work.push(P);
  });return {p,cap,...elect,source,q,work};
 }
 initial(T0){
  let T=T0?[...T0]:Array(this.mesh.n).fill(300);TE.assert(T.length===this.mesh.n,'T0 length mismatch.');
  for(const [i,v] of this.thermal(0).fixed)T[i]=v;this.properties(T);return T;
 }
 implicit(t,target,gammaDt,{initial,tolerance=2e-9,maxIterations=100,relaxation=.85}={}){
  const m=this.mesh,bc=this.thermal(t);let T=[...initial],error=Infinity;
  for(const [i,v] of bc.fixed)T[i]=v;
  for(let iteration=0;iteration<maxIterations;iteration++){
   const b=this.balance(T,t),diag=bc.diag.map((v,i)=>v+(gammaDt?b.cap[i]/gammaDt:0)),rhs=b.source.map((v,i)=>v+bc.rhs[i]+(gammaDt?b.cap[i]/gammaDt*target[i]:0));
   const candidate=TE.graphSolve(m,b.p.map(p=>p.k),diag,rhs,bc.fixed,{initial:T});
   error=Math.max(...candidate.map((v,i)=>Math.abs(v-T[i])));
   TE.assert(candidate.every(v=>Number.isFinite(v)&&v>0),'Nonphysical temperature. Check excitation and material laws.');
   if(error<=tolerance)return candidate;
   T=candidate.map((v,i)=>bc.fixed.has(i)?v:T[i]+relaxation*(v-T[i]));
  }
  throw new Error(`Nonlinear thermal iteration failed (update ${error.toExponential(2)} K). Reduce excitation or refine time steps.`);
 }
 snapshot(T,t=0,steady=false){
  const m=this.mesh,b=this.balance(T,t),Jx=[],Jy=[],qx=[],qy=[];
  m.cells.forEach(c=>{const [a,bb,cc,d]=c.links;Jx.push((b.I[a]+b.I[bb])/(m.dy*m.depth));Jy.push((b.I[cc]+b.I[d])/(m.dx*m.depth));qx.push((b.q[a]+b.q[bb])/(m.dy*m.depth));qy.push((b.q[cc]+b.q[d])/(m.dx*m.depth));});
  const bc=this.thermal(t),res=b.source.map((v,i)=>v+bc.rhs[i]-bc.diag[i]*T[i]);
  m.links.forEach((e,k)=>{const v=b.p[k].k*(T[e.a]-T[e.b]);res[e.a]-=v;res[e.b]+=v;});
  let heatOut=0;for(const [side,faces] of Object.entries(m.sides)){const c=this.boundaries[side],v=TE.signal2D(c.value,t,this.frequency);if(c.kind!=='temperature')for(const {node,A} of faces)heatOut+=A*(c.kind==='flux'?v:c.h*(T[node]-v));}
  for(const [i] of bc.fixed)heatOut+=res[i];
  const electricalPower=-b.current*b.terminalVoltage;
  return {temperature:[...T],voltage:b.V,Jx,Jy,qx,qy,current:b.current,terminalVoltage:b.terminalVoltage,electricalPower,
   energyResidual:steady?heatOut-electricalPower:null,freeResidualWatts:Math.max(0,...res.filter((_,i)=>!bc.fixed.has(i)).map(Math.abs))};
 }
 solveSteady(options={}){
  this.frequency=0;TE.assert(Object.values(this.boundaries).some(b=>b.kind==='temperature'||(b.kind==='convection'&&b.h>0)),'Steady state requires a thermal anchor.');
  const T=this.implicit(0,null,0,{...options,initial:this.initial(options.T0)});
  return {...this.snapshot(T,0,true),method:'steady',converged:true};
 }
 solvePeriodic(frequency,{samples=128,maxPeriods=100,minPeriods=3,periodicAtol=2e-7,periodicRtol=1e-10,tolerance=2e-9,T0,onProgress=()=>{}}={}){
  TE.assert(Number.isFinite(frequency)&&frequency>0,'Frequency must be positive.');TE.assert(Number.isInteger(samples)&&samples>=32&&samples<=2048,'Use 32–2048 integer steps per period.');
  TE.assert(Number.isInteger(maxPeriods)&&maxPeriods>=minPeriods&&minPeriods>=2,'Invalid cycle limit.');
  this.frequency=frequency;let T=this.initial(T0),older=null,previous=null,history,cycle,error=Infinity;const dt=1/(frequency*samples);
  for(cycle=1;cycle<=maxPeriods;cycle++){
   history=[];
   for(let j=0;j<samples;j++){
    history.push([...T]);const target=T.map((v,i)=>older?4*v/3-older[i]/3:v),gamma=older?2/3:1;
    const next=this.implicit(((cycle-1)*samples+j+1)*dt,target,gamma*dt,{initial:T,tolerance});older=T;T=next;
   }
   if(previous){error=0;history.forEach((row,j)=>row.forEach((v,i)=>error=Math.max(error,Math.abs(v-previous[j][i])/(periodicAtol+periodicRtol*Math.max(Math.abs(v),Math.abs(previous[j][i]))))));}
   onProgress({cycle,maxPeriods,error:Number.isFinite(error)?error:null});if(cycle>=minPeriods&&error<=1)break;previous=history;
  }
  TE.assert(cycle<=maxPeriods,`Periodic state not reached in ${maxPeriods} cycles (error ${error.toExponential(2)}).`);
  const traces={},records=history.map((T,i)=>this.snapshot(T,((cycle-1)*samples+i)*dt));
  for(const key of ['temperature','voltage','Jx','Jy','qx','qy','current','terminalVoltage'])traces[key]=records.map(r=>r[key]);
  const harmonics={};for(const [key,value]of Object.entries(traces))harmonics[key]=TE.extractHarmonics(value,3);
  return {...traces,harmonics,time:Array.from({length:samples},(_,i)=>i*dt),frequency,samples,periods:cycle,periodicError:error,finalTemperature:T,method:'BDF2 / nonlinear Picard / matrix-free CG',converged:true};
 }
}
TE.Solver2D=Solver2D;
})(globalThis.TE);
