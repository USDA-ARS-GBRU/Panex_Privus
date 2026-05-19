# ruff: noqa: E501
"""Render self-contained HTML for ``privy interactive --pangenome``."""

from __future__ import annotations

import json
from html import escape
from typing import Any


def render_pangenome_html(data: dict[str, Any]) -> str:
    """Render a self-contained pangenome dashboard HTML document."""
    summary = data["summary"]
    title = escape(str(summary["title"]))
    payload = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>
:root {{
  --ink:#17231f; --muted:#5d6963; --line:#d4ddd7; --paper:#fff; --soft:#f5f7f4;
  --green:#2d765f; --blue:#2d6f9f; --gold:#b8842f; --red:#b75658; --violet:#665da0;
  --shadow:0 14px 40px rgba(25,41,35,.08);
}}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--soft); color:var(--ink); font-family:Inter,ui-sans-serif,system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; line-height:1.45; }}
header {{ padding:26px clamp(18px,3vw,42px) 18px; background:var(--paper); border-bottom:1px solid var(--line); }}
h1 {{ margin:0; font-size:clamp(1.45rem,2.2vw,2.3rem); letter-spacing:0; }}
.subtitle {{ margin:8px 0 0; max-width:980px; color:var(--muted); }}
main {{ padding:20px clamp(14px,2.4vw,34px) 42px; }}
.toolbar {{ display:grid; grid-template-columns:minmax(170px,240px) minmax(170px,240px) minmax(240px,1fr); gap:12px; align-items:end; margin-bottom:16px; }}
label {{ display:block; margin:0 0 6px; color:var(--muted); font-size:.88rem; font-weight:650; }}
select,input[type="search"] {{ width:100%; min-height:42px; border:1px solid #c5d0c9; background:var(--paper); border-radius:7px; padding:8px 10px; color:var(--ink); font:inherit; }}
.metrics {{ display:grid; grid-template-columns:repeat(5,minmax(140px,1fr)); gap:10px; margin-bottom:16px; }}
.metric {{ background:var(--paper); border:1px solid var(--line); border-radius:8px; padding:12px; box-shadow:var(--shadow); }}
.metric span {{ display:block; color:var(--muted); font-size:.82rem; }}
.metric strong {{ display:block; margin-top:3px; font-size:1.35rem; }}
.grid {{ display:grid; grid-template-columns:repeat(2,minmax(300px,1fr)); gap:14px; align-items:start; }}
.panel {{ background:var(--paper); border:1px solid var(--line); border-radius:8px; padding:14px; box-shadow:var(--shadow); }}
.panel h2 {{ margin:0 0 10px; font-size:1rem; }}
.chart {{ width:100%; height:220px; display:block; }}
.tables {{ display:grid; grid-template-columns:minmax(380px,2fr) minmax(320px,1fr); gap:14px; margin-top:14px; align-items:start; }}
.table-wrap {{ overflow:auto; max-height:520px; border:1px solid var(--line); border-radius:8px; }}
table {{ width:100%; border-collapse:collapse; min-width:840px; }}
th,td {{ padding:8px 9px; border-bottom:1px solid #e7ece8; text-align:left; vertical-align:top; font-size:.86rem; }}
th {{ position:sticky; top:0; background:#eef3ef; z-index:1; font-weight:750; }}
tr:hover td {{ background:#f8faf7; }}
.side table {{ min-width:520px; }}
.pill {{ display:inline-block; padding:2px 7px; border-radius:999px; background:#e9f2ee; color:#1d5f4f; font-size:.78rem; font-weight:700; }}
.detail {{ margin-top:10px; color:var(--muted); min-height:40px; }}
.provenance {{ margin-top:14px; }}
.provenance code {{ word-break:break-all; }}
@media (max-width:980px) {{ .toolbar,.metrics,.grid,.tables {{ grid-template-columns:1fr; }} table {{ min-width:760px; }} }}
</style>
</head>
<body>
<header>
  <h1 id="title"></h1>
  <p class="subtitle" id="subtitle"></p>
</header>
<main>
  <section class="toolbar" aria-label="Dashboard controls">
    <div><label for="sourceSelect">Source</label><select id="sourceSelect"></select></div>
    <div><label for="tableMode">Feature Table</label><select id="tableMode"><option value="target_private">Target-private</option><option value="all">All embedded features</option></select></div>
    <div><label for="searchBox">Search features, contigs, types, or categories</label><input id="searchBox" type="search" placeholder="Search embedded features"></div>
  </section>
  <section class="metrics" aria-label="Summary metrics">
    <div class="metric"><span>Sources</span><strong id="metricSources"></strong></div>
    <div class="metric"><span>Total features</span><strong id="metricFeatures"></strong></div>
    <div class="metric"><span>Target-private features</span><strong id="metricTargetPrivate"></strong></div>
    <div class="metric"><span>Off-target-private features</span><strong id="metricOffPrivate"></strong></div>
    <div class="metric"><span>Embedded rows</span><strong id="metricEmbedded"></strong></div>
  </section>
  <section class="grid" aria-label="Interactive pangenome charts">
    <article class="panel"><h2>Composition</h2><canvas class="chart" id="compositionChart"></canvas></article>
    <article class="panel"><h2>Coverage Histogram</h2><canvas class="chart" id="coverageChart"></canvas></article>
    <article class="panel"><h2>Growth Curves</h2><canvas class="chart" id="growthChart"></canvas></article>
    <article class="panel"><h2>Feature Types</h2><canvas class="chart" id="featureTypeChart"></canvas></article>
  </section>
  <section class="tables">
    <article class="panel">
      <h2>Feature Inventory</h2>
      <div class="table-wrap">
        <table>
          <thead><tr><th>Source</th><th>Feature</th><th>Position</th><th>Type</th><th>Length</th><th>Target</th><th>Off-target</th><th>Categories</th></tr></thead>
          <tbody id="featureRows"></tbody>
        </table>
      </div>
      <div class="detail" id="featureDetail">Select a feature row for details.</div>
    </article>
    <div class="side">
      <article class="panel">
        <h2>Top Contigs</h2>
        <div class="table-wrap">
          <table>
            <thead><tr><th>Contig</th><th>Features</th></tr></thead>
            <tbody id="contigRows"></tbody>
          </table>
        </div>
      </article>
      <article class="panel provenance">
        <h2>Provenance</h2>
        <div id="provenance"></div>
      </article>
    </div>
  </section>
</main>
<script>
const DATA = {payload};
const S = DATA.summary;
const state = {{ source: "all", tableMode: "target_private", query: "" }};
const colors = ["#2d765f", "#2d6f9f", "#b8842f", "#b75658", "#665da0", "#4f7f86", "#8b6f47", "#657a56"];
const categoryColors = {{ absent:"#b75658", private:"#b8842f", accessory:"#2d6f9f", core:"#2d765f" }};
function esc(value) {{ return String(value ?? "").replace(/[&<>"']/g, ch => ({{"&":"&amp;","<":"&lt;",">":"&gt;","\\"":"&quot;","'":"&#39;"}}[ch])); }}
function fmt(value) {{ return Number(value || 0).toLocaleString(); }}
function pct(value) {{ return Number(value || 0).toFixed(3); }}
function sourcePayload() {{ return state.source === "all" ? null : DATA.sources.find(src => src.key === state.source) || null; }}
function rowsFor(name) {{ const src = sourcePayload(); return src ? (src[name] || []) : (DATA.aggregate[name] || []); }}
function setupSources() {{
  const select = document.getElementById("sourceSelect");
  select.innerHTML = '<option value="all">All sources</option>';
  for (const src of DATA.sources) {{
    const opt = document.createElement("option"); opt.value = src.key; opt.textContent = src.label; select.appendChild(opt);
  }}
}}
function sizeCanvas(canvas) {{
  const rect = canvas.getBoundingClientRect();
  const dpr = window.devicePixelRatio || 1;
  canvas.width = Math.max(1, Math.floor(rect.width * dpr));
  canvas.height = Math.max(1, Math.floor(rect.height * dpr));
  const ctx = canvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return {{ctx, width: rect.width, height: rect.height}};
}}
function drawGroupedBars(canvasId, rows, valueKey = "n_features") {{
  const canvas = document.getElementById(canvasId);
  const {{ctx, width, height}} = sizeCanvas(canvas);
  ctx.clearRect(0,0,width,height);
  const pad = {{left:44,right:12,top:14,bottom:48}};
  const plotW = width - pad.left - pad.right, plotH = height - pad.top - pad.bottom;
  const labels = [...new Set(rows.map(row => `${{row.source_label || ""}} ${{row.group || row.label || ""}}`.trim()))].slice(0, 18);
  const max = Math.max(1, ...rows.map(row => Number(row[valueKey] || row.count || 0)));
  ctx.strokeStyle = "#d4ddd7"; ctx.beginPath(); ctx.moveTo(pad.left,pad.top); ctx.lineTo(pad.left,pad.top+plotH); ctx.lineTo(pad.left+plotW,pad.top+plotH); ctx.stroke();
  const barW = labels.length ? Math.max(5, plotW / labels.length * .72) : plotW;
  labels.forEach((label, idx) => {{
    const matching = rows.filter(row => `${{row.source_label || ""}} ${{row.group || row.label || ""}}`.trim() === label);
    const total = matching.reduce((a,row)=>a+Number(row[valueKey] || row.count || 0),0);
    let y = pad.top + plotH;
    const x = pad.left + idx * (plotW / Math.max(1, labels.length)) + (plotW / Math.max(1, labels.length) - barW) / 2;
    for (const row of matching) {{
      const value = Number(row[valueKey] || row.count || 0);
      const h = value / max * plotH;
      y -= h;
      ctx.fillStyle = categoryColors[row.category] || colors[idx % colors.length];
      ctx.fillRect(x, y, barW, h);
    }}
    if (!matching.length || matching.length === 1) {{
      ctx.fillStyle = colors[idx % colors.length];
      ctx.fillRect(x, pad.top + plotH - total / max * plotH, barW, total / max * plotH);
    }}
    ctx.save(); ctx.translate(x + barW / 2, pad.top + plotH + 8); ctx.rotate(-Math.PI / 5);
    ctx.fillStyle = "#5d6963"; ctx.font = "11px system-ui, sans-serif"; ctx.textAlign = "center"; ctx.fillText(label.slice(0,22), 0, 0); ctx.restore();
  }});
  ctx.fillStyle = "#5d6963"; ctx.font = "12px system-ui, sans-serif"; ctx.textAlign = "right"; ctx.fillText(fmt(max), pad.left - 6, pad.top + 6);
}}
function drawGrowth() {{
  const canvas = document.getElementById("growthChart");
  const {{ctx, width, height}} = sizeCanvas(canvas);
  ctx.clearRect(0,0,width,height);
  const rows = rowsFor("growth");
  const pad = {{left:44,right:14,top:14,bottom:34}};
  const plotW = width - pad.left - pad.right, plotH = height - pad.top - pad.bottom;
  if (!rows.length) {{ ctx.fillStyle = "#5d6963"; ctx.fillText("No growth_curves.tsv data available.", 18, 32); return; }}
  const seriesKeys = [...new Set(rows.map(row => `${{row.source_label || ""}}:${{row.group}}`))];
  const maxN = Math.max(1, ...rows.map(row => Number(row.n || 0)));
  const maxY = Math.max(1, ...rows.map(row => Number(row.mean_features || 0)));
  ctx.strokeStyle = "#d4ddd7"; ctx.beginPath(); ctx.moveTo(pad.left,pad.top); ctx.lineTo(pad.left,pad.top+plotH); ctx.lineTo(pad.left+plotW,pad.top+plotH); ctx.stroke();
  seriesKeys.forEach((key, idx) => {{
    const items = rows.filter(row => `${{row.source_label || ""}}:${{row.group}}` === key).sort((a,b)=>a.n-b.n);
    ctx.strokeStyle = colors[idx % colors.length]; ctx.lineWidth = 2; ctx.beginPath();
    items.forEach((row, j) => {{
      const x = pad.left + Number(row.n || 0) / maxN * plotW;
      const y = pad.top + plotH - Number(row.mean_features || 0) / maxY * plotH;
      if (j === 0) ctx.moveTo(x,y); else ctx.lineTo(x,y);
    }});
    ctx.stroke();
  }});
  ctx.fillStyle = "#5d6963"; ctx.font = "12px system-ui, sans-serif"; ctx.textAlign = "right"; ctx.fillText(fmt(maxY), pad.left - 6, pad.top + 6);
}}
function renderMetrics() {{
  const src = sourcePayload();
  document.getElementById("metricSources").textContent = src ? src.label : fmt(S.source_count);
  document.getElementById("metricFeatures").textContent = fmt(src ? src.summary.feature_count : S.feature_count);
  document.getElementById("metricTargetPrivate").textContent = fmt(src ? src.summary.target_private_count : S.target_private_count);
  document.getElementById("metricOffPrivate").textContent = fmt(src ? src.summary.offtarget_private_count : S.offtarget_private_count);
  document.getElementById("metricEmbedded").textContent = fmt(src ? src.summary.embedded_features : S.embedded_features);
}}
function activeFeatures() {{
  const base = state.tableMode === "target_private" ? DATA.target_private_features : DATA.features;
  const q = state.query.trim().toLowerCase();
  return base.filter(row => {{
    if (state.source !== "all" && row.source !== state.source) return false;
    if (!q) return true;
    return [row.feature_id,row.source_label,row.feature_type,row.contig,row.full_category,row.target_category,row.offtarget_category].join(" ").toLowerCase().includes(q);
  }});
}}
function renderFeatures() {{
  const tbody = document.getElementById("featureRows");
  const rows = activeFeatures().slice(0, 700);
  tbody.innerHTML = "";
  if (!rows.length) {{ tbody.innerHTML = '<tr><td colspan="8">No embedded features match the current filters.</td></tr>'; return; }}
  for (const row of rows) {{
    const tr = document.createElement("tr");
    const position = row.contig ? `${{row.contig}}:${{fmt(row.start)}}-${{fmt(row.end)}}` : "NA";
    tr.innerHTML = `<td>${{esc(row.source_label)}}</td><td>${{esc(row.feature_id)}}</td><td>${{esc(position)}}</td><td><span class="pill">${{esc(row.feature_type)}}</span></td><td>${{fmt(row.length)}}</td><td>${{fmt(row.target_present_n)}}/${{fmt(row.target_total_n)}}</td><td>${{fmt(row.offtarget_present_n)}}/${{fmt(row.offtarget_total_n)}}</td><td>full=${{esc(row.full_category)}}; target=${{esc(row.target_category)}}; off=${{esc(row.offtarget_category)}}</td>`;
    tr.addEventListener("click", () => {{
      document.getElementById("featureDetail").textContent = `${{row.feature_id}}: target_private=${{row.target_private}}, offtarget_private=${{row.offtarget_private}}, source_type=${{row.source_type}}.`;
    }});
    tbody.appendChild(tr);
  }}
}}
function renderContigs() {{
  const tbody = document.getElementById("contigRows");
  const rows = rowsFor("contig_counts").slice(0, 40);
  tbody.innerHTML = "";
  if (!rows.length) {{ tbody.innerHTML = '<tr><td colspan="2">No coordinate-backed features.</td></tr>'; return; }}
  for (const row of rows) tbody.innerHTML += `<tr><td>${{esc(row.label)}}</td><td>${{fmt(row.count)}}</td></tr>`;
}}
function renderProvenance() {{
  const src = sourcePayload();
  const sources = src ? [src] : DATA.sources;
  let html = "";
  for (const item of sources) {{
    const input = S.inputs.find(x => x.label === item.label);
    html += `<p><strong>${{esc(item.label)}}</strong><br><code>${{esc(input?.path || "")}}</code></p>`;
    html += `<p>source_type=${{esc(item.summary.source_type || "unknown")}}; samples=${{fmt(item.summary.n_samples)}}; targets=${{fmt(item.summary.n_target_samples)}}; off-targets=${{fmt(item.summary.n_offtarget_samples)}}</p>`;
  }}
  document.getElementById("provenance").innerHTML = html;
}}
function updateCharts() {{
  drawGroupedBars("compositionChart", rowsFor("composition"), "n_features");
  drawGroupedBars("coverageChart", rowsFor("coverage"), "n_features");
  drawGrowth();
  drawGroupedBars("featureTypeChart", rowsFor("feature_type_counts"), "count");
}}
function update() {{
  renderMetrics();
  updateCharts();
  renderFeatures();
  renderContigs();
  renderProvenance();
}}
document.getElementById("title").textContent = S.title;
document.getElementById("subtitle").textContent = S.subtitle;
setupSources();
document.getElementById("sourceSelect").addEventListener("change", e => {{ state.source = e.target.value; update(); }});
document.getElementById("tableMode").addEventListener("change", e => {{ state.tableMode = e.target.value; renderFeatures(); }});
document.getElementById("searchBox").addEventListener("input", e => {{ state.query = e.target.value; renderFeatures(); }});
window.addEventListener("resize", updateCharts);
update();
</script>
</body>
</html>
"""
