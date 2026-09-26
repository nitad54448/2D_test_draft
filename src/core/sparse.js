(function(TE){
// Matrix-free symmetric graph Laplacian + positive diagonal. Dirichlet eliminated.
TE.graphSolve=(mesh,conductance,diagonal,rhs,fixed,{rtol=2e-12,maxIter=4000,initial}={})=>{
 const n=mesh.n,idx=new Int32Array(n).fill(-1),free=[];
 for(let i=0;i<n;i++)if(!fixed.has(i)){idx[i]=free.length;free.push(i);}
 const nf=free.length,out=new Float64Array(n);for(const [i,v] of fixed)out[i]=v;
 if(!nf)return Array.from(out);
 const D=new Float64Array(nf),b=new Float64Array(nf),x=new Float64Array(nf);
 for(let k=0;k<nf;k++){const i=free[k];D[k]=diagonal[i]??0;b[k]=rhs[i];x[k]=initial?initial[i]:0;}
 const edges=[];
 mesh.links.forEach((e,l)=>{const g=conductance[l],a=idx[e.a],bb=idx[e.b];TE.assert(Number.isFinite(g)&&g>0,'Invalid face conductance.');
  if(a>=0){D[a]+=g;if(bb<0)b[a]+=g*fixed.get(e.b);}
  if(bb>=0){D[bb]+=g;if(a<0)b[bb]+=g*fixed.get(e.a);}
  if(a>=0&&bb>=0)edges.push([a,bb,g]);
 });
 for(const d of D)TE.assert(d>0&&Number.isFinite(d),'Unanchored or invalid system.');
 const multiply=v=>{const z=Float64Array.from(v,(q,i)=>D[i]*q);for(const [a,bb,g] of edges){z[a]-=g*v[bb];z[bb]-=g*v[a];}return z;};
 const dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0);
 let Ax=multiply(x),r=Float64Array.from(b,(v,i)=>v-Ax[i]),z=Float64Array.from(r,(v,i)=>v/D[i]),p=Float64Array.from(z),rz=dot(r,z);
 const target=rtol*Math.max(Math.sqrt(dot(b,b)),1e-20);
 let iteration=0;
 for(;Math.sqrt(dot(r,r))>target&&iteration<maxIter;iteration++){
  const Ap=multiply(p),den=dot(p,Ap);TE.assert(den>0&&Number.isFinite(den),'CG failed: matrix not positive definite.');
  const a=rz/den;for(let i=0;i<nf;i++){x[i]+=a*p[i];r[i]-=a*Ap[i];}
  z=Float64Array.from(r,(v,i)=>v/D[i]);const next=dot(r,z),beta=next/rz;for(let i=0;i<nf;i++)p[i]=z[i]+beta*p[i];rz=next;
 }
 TE.assert(iteration<maxIter,'Sparse linear solver did not converge.');
 free.forEach((node,i)=>out[node]=x[i]);return Array.from(out);
};
})(globalThis.TE);
