const fs=require('node:fs'),path=require('node:path'),read=p=>fs.readFileSync(path.join(__dirname,p),'utf8');
const core=['materials','fourier','mesh2d','sparse','solver2d','validation','config2d'].map(n=>read('src/core/'+n+'.js')).join('\n');
let html=read('src/index.template.html');for(const [marker,value]of Object.entries({'/*STYLE*/':read('src/styles.css'),'/*CORE*/':core+'\n'+read('src/ui-helpers.js'),'/*WORKER*/':core+'\n'+read('src/worker.js'),'/*APP*/':read('src/app.js')}))html=html.replace(marker,()=>value);
fs.writeFileSync(path.join(__dirname,'index.html'),html);console.log('Built standalone 2D app:',Buffer.byteLength(html),'bytes');
