/* Semantics checks for round-based observations and fixed-budget discovery rates. */
const assert = require('node:assert/strict');
const charts = require('./static/offline_charts.js');
const config = {methods:['a','b'],repeats:4,n_init:4,total_budget:9,batch_size:3,target:{column:'CRR[%]',direction:'max'}};
const trial = (method,values,status='completed') => ({method,status,curve:values.map((value,round)=>({round,evaluations:[4,7,9][round],best_value:value}))});
const results = {trials:[trial('a',[10,30,50]),trial('a',[20,40,60]),trial('a',[999,999],'running'),trial('b',[-10,-5,0])],
  curves:[{method:'a',evaluations:4,hit_rate:0},{method:'a',evaluations:7,hit_rate:.25},{method:'a',evaluations:9,hit_rate:.5}]};
const data=charts.series(config,results);
assert.equal(data.maxRound,2); // Tail batch of 2 still represents one round.
assert.deepEqual(data.best.filter(p=>p.method==='a').map(p=>[p.round,p.value,p.count]),[[0,15,2],[1,35,2],[2,55,2]]);
assert.equal(data.best.find(p=>p.method==='a'&&p.round===2).q25,52.5);
assert.equal(data.best.find(p=>p.method==='b'&&p.round===2).value,0);
assert.equal(data.discovery.at(-1).value,.5); // Keep planned-repeat denominator from API.
const minimum={...config,target:{column:'loss',direction:'min'}};
const minResults={trials:[trial('a',[0,-3,-8]),trial('a',[-2,-5,-10])],curves:[]};
assert.deepEqual(charts.series(minimum,minResults).best.map(p=>p.value),[-1,-4,-9]);
for(const [kind,c,r,reference] of [['best',config,results,70],['best',minimum,minResults,-10],['discovery',config,results],['best',config,{trials:[],curves:[]},0]]) {
  const svg=charts.render(kind,c,r,{colors:['blue','orange'],reference});
  assert.ok(!/NaN|Infinity/.test(svg));
  assert.ok(svg.includes('<svg'));
}
assert.ok(charts.render('best',config,{trials:[],curves:[]}).includes('等待已完成的重复'));
assert.ok(charts.render('best',{...config,target:{column:'<unsafe>',direction:'max'}},results).includes('&lt;unsafe&gt;'));
console.log('Offline charts: round, tail batch, repeats, minimization, zero/negative values, empty state and discovery checks passed.');
