/* Offline visualization reads recorded observations; no model or replay behavior changes. */
const OfflineCharts = (() => {
  const escapeText = value => String(value ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const finite = value => typeof value === 'number' && Number.isFinite(value);
  const display = value => Number(value).toLocaleString('zh-CN', {maximumFractionDigits: 3});
  const dashes = ['', '7 3', '2 3', '9 3 2 3', '4 3', '12 4', '6 2 2 2'];

  function quantile(values, fraction) {
    const sorted = [...values].sort((a, b) => a-b);
    const position = (sorted.length-1)*fraction, low = Math.floor(position);
    return sorted[low] + (sorted[Math.ceil(position)]-sorted[low])*(position-low);
  }

  function series(config, results) {
    const maxRound = Math.ceil((config.total_budget-config.n_init)/config.batch_size);
    const best = [];
    for (const method of config.methods) {
      const complete = (results.trials || []).filter(t => t.method === method && t.status === 'completed');
      for (let round = 0; round <= maxRound; round++) {
        const values = complete.map(t => (t.curve || []).find(p => p.round === round)?.best_value).filter(finite);
        if (values.length) best.push({method, round, value:quantile(values,.5),
          q25:quantile(values,.25), q75:quantile(values,.75), count:values.length});
      }
    }
    const discovery = (results.curves || []).filter(p => finite(p.hit_rate)).map(p => ({
      method:p.method, evaluations:p.evaluations, value:p.hit_rate, count:config.repeats}));
    return {best, discovery, maxRound};
  }

  function render(kind, config, results, {colors, names, reference} = {}) {
    const data = series(config, results), performance = kind === 'best';
    const points = performance ? data.best : data.discovery;
    const W=600, H=300, L=72, R=22, T=34, B=54;
    const start = performance ? 0 : config.n_init, end = performance ? data.maxRound : config.total_budget;
    const domain = points.flatMap(p => performance ? [p.q25,p.q75] : [p.value]);
    if (performance && finite(reference)) domain.push(reference);
    let low=0, high=1;
    if (performance && domain.length) {
      low=Math.min(...domain); high=Math.max(...domain);
      const padding=Math.max((high-low)*.12, Math.max(Math.abs(low),Math.abs(high))*.015,1e-9);
      low-=padding; high+=padding;
    }
    const x=value => L+(value-start)/Math.max(1,end-start)*(W-L-R);
    const y=value => H-B-(value-low)/(high-low)*(H-T-B);
    const path=(rows,key) => rows.map((p,i)=>`${i?'H':'M'}${x(p[performance?'round':'evaluations'])}${i?'V':','}${y(p[key])}`).join('');
    let svg='';
    for(let i=0;i<=4;i++) {
      const value=low+(high-low)*i/4;
      svg+=`<line x1="${L}" x2="${W-R}" y1="${y(value)}" y2="${y(value)}" stroke="#e9eef5"/><text x="${L-10}" y="${y(value)+4}" text-anchor="end" fill="#8190a3">${performance?display(value):Math.round(value*100)+'%'}</text>`;
    }
    const ticks=[...new Set(Array.from({length:5},(_,i)=>Math.round(start+(end-start)*i/4)))];
    for(const tick of ticks) svg+=`<text x="${x(tick)}" y="${H-B+22}" text-anchor="middle" fill="#8190a3">${tick}</text>`;
    svg+=`<text x="${L}" y="17" fill="#61748f">${escapeText(performance?config.target.column:'发现率')}</text>`;
    svg+=`<text x="${(W+L-R)/2}" y="${H-7}" text-anchor="middle" fill="#61748f">${performance?'优化轮次（0 = 初始化）':'累计实验配方数（含初始）'}</text>`;
    if(performance && finite(reference)) svg+=`<line x1="${L}" x2="${W-R}" y1="${y(reference)}" y2="${y(reference)}" stroke="#a4afbd" stroke-dasharray="4 4"/><text x="${W-R}" y="${y(reference)-7}" text-anchor="end" fill="#7f8da0">数据集最优 ${display(reference)}</text>`;
    config.methods.forEach((method,index)=>{
      const rows=points.filter(p=>p.method===method); if(!rows.length)return;
      const color=colors?.[index%colors.length] || '#276be2', dash=dashes[index%dashes.length];
      if(performance) {
        const reverse=[...rows].reverse();
        const lower=`L${x(reverse[0].round)},${y(reverse[0].q25)}`+reverse.slice(1).map(p=>`V${y(p.q25)}H${x(p.round)}`).join('');
        svg+=`<path d="${path(rows,'q75')}${lower}Z" fill="${color}" opacity=".08"/>`;
      }
      svg+=`<path d="${path(rows,'value')}" fill="none" stroke="${color}" stroke-width="2.4" stroke-dasharray="${dash}"><title>${escapeText(names?.[method]||method)}</title></path>`;
      for(const p of rows) svg+=`<circle cx="${x(p[performance?'round':'evaluations'])}" cy="${y(p.value)}" r="2.7" fill="${color}" stroke="white" stroke-width=".7"><title>${escapeText(names?.[method]||method)} · ${performance?'第 '+p.round+' 轮':'累计 '+p.evaluations+' 个配方'} · ${performance?display(p.value):display(p.value*100)+'%'}${performance?' · 已完成 '+p.count+' 次重复':''}</title></circle>`;
    });
    if(!points.length) svg+=`<text x="${(W+L-R)/2}" y="${(T+H-B)/2}" text-anchor="middle" fill="#91a0b2">等待已完成的重复</text>`;
    return `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="${performance?'轮次与已发现的最佳目标性能':'累计实验数与最优配方发现率'}">${svg}</svg>`;
  }
  return {series, render};
})();
if (typeof module !== 'undefined' && module.exports) module.exports = OfflineCharts;
