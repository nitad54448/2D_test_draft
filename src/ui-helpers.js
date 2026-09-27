/* Presentation helpers: displayed rounding never replaces an unchanged raw value. */
(function(TE){
 TE.formatInputNumber=value=>{
  if(!Number.isFinite(value))return String(value);
  return String(Number(value.toPrecision(10)));
 };
 TE.setNumberInput=(element,value)=>{
  const display=TE.formatInputNumber(value);
  element.value=display;element.dataset.rawNumber=String(value);element.dataset.displayNumber=display;
 };
 TE.readNumberInput=element=>{
  TE.assert(element.value.trim()!=='','A numerical value is required.');
  const raw=element.dataset.rawNumber;
  return TE.finite(Number(raw!==undefined&&element.value===element.dataset.displayNumber?raw:element.value),'Numerical value');
 };
 TE.historyForPlot=(result,values)=>{
  const time=result.time.map(v=>v*1000),y=[...values];
  // Only a converged periodic history may be closed back to its first sample.
  if(result.converged){time.push(1000/result.frequency);y.push(values[0]);}
  return {time,values:y};
 };
 TE.remeshConfig=(config,{nx,ny,lx,ly,depth})=>{
  // Validate before allocating or mutating the displayed model.
  TE.assert(Number.isInteger(nx)&&Number.isInteger(ny)&&nx>=2&&ny>=1,'Element counts must be integers: Nx ≥ 2 and Ny ≥ 1.');
  TE.assert((nx+1)*(ny+1)<=1600,`Requested mesh: ${nx*ny} elements, ${(nx+1)*(ny+1)} nodes. The browser limit is 1600 nodes; reduce Nx or Ny.`);
  TE.assert([lx,ly,depth].every(v=>Number.isFinite(v)&&v>0),'Dimensions and depth must be positive.');
  const next={...config,nx,ny,lx,ly,depth};
  next.materialMap=Array.from({length:nx*ny},(_,k)=>{
   const i=k%nx,j=Math.floor(k/nx),oldI=Math.min(config.nx-1,Math.floor((i+.5)/nx*config.nx)),oldJ=Math.min(config.ny-1,Math.floor((j+.5)/ny*config.ny));
   return config.materialMap[oldJ*config.nx+oldI];
  });
  return next;
 };
})(globalThis.TE);
