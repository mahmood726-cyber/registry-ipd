/* Generate compact, internally-consistent bundled examples for the HTML shell.
 * Writes examples/examples.js (window.RIPD_EXAMPLES). Run: node examples/gen_examples.js */
const fs = require('fs'), path = require('path');
function mulberry32(s){let a=s>>>0;return function(){a|=0;a=(a+0x6D2B79F5)|0;let t=Math.imul(a^(a>>>15),1|a);t=(t+Math.imul(t^(t>>>7),61|t))^t;return((t^(t>>>14))>>>0)/4294967296;};}
function km(ipd){const rows=ipd.slice().sort((a,b)=>a.time-b.time);const times=[...new Set(rows.map(r=>r.time))].sort((a,b)=>a-b);let n=rows.length,S=1;const st=[];for(const t of times){const at=rows.filter(r=>Math.abs(r.time-t)<1e-12);const d=at.filter(r=>r.status===1).length,c=at.filter(r=>r.status===0).length;if(d>0)S*=(1-d/n);st.push({t,S,n});n-=(d+c);}return st;}
function evalKM(st,t){let S=1;for(const s of st){if(s.t<=t+1e-9)S=s.S;else break;}return S;}
function nar(ipd,t){return ipd.filter(r=>r.time>=t-1e-9).length;}
function expArm(seed,N,med,cut){const rng=mulberry32(seed);const lam=Math.log(2)/med;const ipd=[];for(let i=0;i<N;i++){const u=Math.max(1e-12,rng());const T=-Math.log(u)/lam;ipd.push({time:Math.min(T,cut),status:T<=cut?1:0});}return ipd;}
function round(x,d){const f=Math.pow(10,d);return Math.round(x*f)/f;}

function tierAarm(seed,N,med,role,label,id){
  const ipd=expArm(seed,N,med,30);
  const st=km(ipd);
  const kmTs=[0,4,8,12,16,20,24];
  const km_points=kmTs.map(t=>({t,S:round(evalKM(st,t),4)}));
  const nar_points=[0,12,24].map(t=>({t,n:nar(ipd,t)}));
  const total_events=ipd.filter(r=>r.status===1).length;
  return {arm_id:id,label,role,N,total_events,follow_up_max:24,km_points,nar_points};
}

const EX={};

// Tier A — rich registry data (Guyot inverse-KM)
(function(){
  const ctl=tierAarm(11,200,11,'comparator','Placebo','OG001');
  const exp=tierAarm(22,200,17,'experimental','Drug A','OG000');
  const hr=round((Math.log(2)/17)/(Math.log(2)/11),3);
  EX.tierA={label:'Tier A — rich (KM points + at-risk)',trial:{
    nct_id:'EXAMPLE-A',source_url:'synthetic example',time_unit:'months',
    arms:[exp,ctl],
    hr:{value:hr,ci_low:0.5,ci_high:0.86,method:'Cox',favors_arm_id:'OG000'}
  }};
})();

// Tier B — median + HR only (parametric + bootstrap envelope)
EX.tierB={label:'Tier B — medium (median + HR, parametric)',trial:{
  nct_id:'EXAMPLE-B',source_url:'synthetic example',time_unit:'months',
  arms:[
    {arm_id:'OG000',label:'Drug B',role:'experimental',N:300,total_events:165,follow_up_max:30,median:{value:16,ci_low:13,ci_high:20.5},km_points:[],nar_points:[]},
    {arm_id:'OG001',label:'Placebo',role:'comparator',N:300,total_events:200,follow_up_max:30,median:{value:10,ci_low:8.5,ci_high:11.8},km_points:[],nar_points:[]}
  ],
  hr:{value:0.625,ci_low:0.50,ci_high:0.78,method:'Cox',favors_arm_id:'OG000'}
}};

// Tier C — HR only (fail closed)
EX.tierC={label:'Tier C — sparse (HR only, refused)',trial:{
  nct_id:'EXAMPLE-C',source_url:'synthetic example',time_unit:'months',
  arms:[
    {arm_id:'OG000',label:'Drug C',role:'experimental',N:null,total_events:null,km_points:[],nar_points:[],median:null},
    {arm_id:'OG001',label:'Placebo',role:'comparator',N:null,total_events:null,km_points:[],nar_points:[],median:null}
  ],
  hr:{value:0.82,ci_low:0.66,ci_high:1.02,method:'Cox',favors_arm_id:'OG000'}
}};

// Real ct.gov trials harvested from AACT (if present) — provenance demos, not synthetic
function addReal(key, file, label) {
  try { EX[key] = { label, trial: JSON.parse(fs.readFileSync(path.join(__dirname, '..', file), 'utf8')) }; }
  catch (e) { /* file not present; skip */ }
}
addReal('realOnc', 'NCT01524783.json', 'REAL ct.gov — NCT01524783 (RADIANT-4, NET; high-event OS)');
addReal('real', 'NCT00725985.json', 'REAL ct.gov — NCT00725985 (cladribine MS; low-event)');

const out='window.RIPD_EXAMPLES = '+JSON.stringify(EX,null,2)+';\n';
fs.writeFileSync(path.join(__dirname,'examples.js'),out);
console.log('wrote examples/examples.js with',Object.keys(EX).length,'examples');
