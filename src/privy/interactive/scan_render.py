# ruff: noqa: E501
"""Render self-contained HTML for ``privy interactive --scan``."""

from __future__ import annotations

import json
from html import escape
from typing import Any

from privy.interactive.branding import FOOTER_HTML, TOP_CREDIT_HTML


def render_scan_html(data: dict[str, Any]) -> str:
    """Render a self-contained scan dashboard HTML document."""
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
  --ink: #17231f;
  --muted: #5a6761;
  --line: #d5ded8;
  --soft: #f5f7f4;
  --paper: #ffffff;
  --green: #28745f;
  --blue: #2e6f9e;
  --gold: #b7832f;
  --red: #b64b52;
  --violet: #6a5d9f;
  --shadow: 0 14px 40px rgba(25, 41, 35, .08);
}}
* {{ box-sizing: border-box; }}
body {{
  margin: 0;
  background: var(--soft);
  color: var(--ink);
  font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
  line-height: 1.45;
}}
header {{
  padding: 26px clamp(18px, 3vw, 42px) 18px;
  background: var(--paper);
  border-bottom: 1px solid var(--line);
}}
h1 {{ margin: 0; font-size: clamp(1.45rem, 2.2vw, 2.3rem); letter-spacing: 0; }}
.subtitle {{ margin: 8px 0 0; max-width: 980px; color: var(--muted); }}
.tool-credit {{ margin: 6px 0 0; color: var(--muted); }}
.tool-credit a, .site-footer a {{ color: #0b5c8e; }}
main {{ padding: 20px clamp(14px, 2.4vw, 34px) 42px; }}
.toolbar {{
  display: grid;
  grid-template-columns: minmax(180px, 260px) minmax(220px, 1fr) minmax(180px, 260px);
  gap: 12px;
  align-items: end;
  margin-bottom: 16px;
}}
label {{ display: block; margin: 0 0 6px; color: var(--muted); font-size: .88rem; font-weight: 650; }}
select, input[type="search"], input[type="range"] {{
  width: 100%;
  min-height: 42px;
  border: 1px solid #c5d0c9;
  background: var(--paper);
  border-radius: 7px;
  padding: 8px 10px;
  color: var(--ink);
  font: inherit;
}}
.metrics {{
  display: grid;
  grid-template-columns: repeat(5, minmax(140px, 1fr));
  gap: 10px;
  margin-bottom: 16px;
}}
.metric {{
  background: var(--paper);
  border: 1px solid var(--line);
  border-radius: 8px;
  padding: 12px;
  box-shadow: var(--shadow);
}}
.metric span {{ display: block; color: var(--muted); font-size: .82rem; }}
.metric strong {{ display: block; margin-top: 3px; font-size: 1.35rem; }}
.grid {{
  display: grid;
  grid-template-columns: repeat(2, minmax(260px, 1fr));
  gap: 14px;
}}
.panel {{
  background: var(--paper);
  border: 1px solid var(--line);
  border-radius: 8px;
  padding: 14px;
  box-shadow: var(--shadow);
}}
.panel h2 {{ margin: 0 0 10px; font-size: 1rem; }}
.chart {{ width: 100%; height: 190px; display: block; }}
.filter-groups {{
  display: grid;
  grid-template-columns: repeat(2, minmax(220px, 1fr));
  gap: 14px;
  margin: 16px 0;
}}
.checkrow {{
  display: flex;
  flex-wrap: wrap;
  gap: 8px;
}}
.checkrow label {{
  display: inline-flex;
  align-items: center;
  gap: 6px;
  margin: 0;
  border: 1px solid var(--line);
  background: #fbfcfa;
  padding: 6px 9px;
  border-radius: 999px;
  color: var(--ink);
  font-size: .86rem;
}}
.tables {{
  display: grid;
  grid-template-columns: minmax(360px, 2fr) minmax(300px, 1fr);
  gap: 14px;
  align-items: start;
  margin-top: 14px;
}}
.table-wrap {{ overflow: auto; max-height: 520px; border: 1px solid var(--line); border-radius: 8px; }}
table {{ width: 100%; border-collapse: collapse; min-width: 820px; }}
th, td {{ padding: 8px 9px; border-bottom: 1px solid #e7ece8; text-align: left; vertical-align: top; font-size: .86rem; }}
th {{ position: sticky; top: 0; background: #eef3ef; z-index: 1; font-weight: 750; }}
tr:hover td {{ background: #f8faf7; }}
.side table {{ min-width: 520px; }}
.pill {{
  display: inline-block;
  padding: 2px 7px;
  border-radius: 999px;
  background: #e9f2ee;
  color: #1d5f4f;
  font-size: .78rem;
  font-weight: 700;
}}
.detail {{
  margin-top: 12px;
  color: var(--muted);
  min-height: 42px;
}}
.muted {{ color: var(--muted); }}
.provenance {{ margin-top: 14px; }}
.provenance code {{ word-break: break-all; }}
.site-footer {{
  padding: 0 clamp(14px, 2.4vw, 34px) 32px;
  color: var(--muted);
  font-size: .92rem;
}}
.site-footer div {{ border-top: 1px solid var(--line); padding-top: 14px; }}
@media (max-width: 980px) {{
  .toolbar, .metrics, .grid, .filter-groups, .tables {{ grid-template-columns: 1fr; }}
  table {{ min-width: 760px; }}
}}
</style>
</head>
<body>
<header>
  <h1 id="title"></h1>
  <p class="tool-credit">{TOP_CREDIT_HTML}</p>
  <p class="subtitle" id="subtitle"></p>
</header>
<main>
  <section class="toolbar" aria-label="Dashboard controls">
    <div>
      <label for="sourceSelect">Source</label>
      <select id="sourceSelect"></select>
    </div>
    <div>
      <label for="searchBox">Search loci, contigs, alleles, or classes</label>
      <input id="searchBox" type="search" placeholder="Search embedded top hits">
    </div>
    <div>
      <label for="scoreSlider">Minimum final score: <span id="scoreLabel">0.00</span></label>
      <input id="scoreSlider" type="range" min="0" max="2" step="0.01" value="0">
    </div>
  </section>
  <section class="metrics" aria-label="Summary metrics">
    <div class="metric"><span>Sources</span><strong id="metricSources"></strong></div>
    <div class="metric"><span>Total hits</span><strong id="metricHits"></strong></div>
    <div class="metric"><span>Total regions</span><strong id="metricRegions"></strong></div>
    <div class="metric"><span>Embedded hit rows</span><strong id="metricEmbeddedHits"></strong></div>
    <div class="metric"><span>Embedded region rows</span><strong id="metricEmbeddedRegions"></strong></div>
  </section>
  <section class="grid" aria-label="Interactive charts">
    <article class="panel"><h2>Strictness Classes</h2><canvas class="chart" id="strictnessChart"></canvas></article>
    <article class="panel"><h2>Variant Types</h2><canvas class="chart" id="variantChart"></canvas></article>
    <article class="panel"><h2>Final Score Distribution</h2><canvas class="chart" id="scoreChart"></canvas></article>
    <article class="panel"><h2>Top Contigs By Hit Count</h2><canvas class="chart" id="contigChart"></canvas></article>
  </section>
  <section class="filter-groups" aria-label="Hit filters">
    <article class="panel"><h2>Strictness Filter</h2><div class="checkrow" id="strictnessFilters"></div></article>
    <article class="panel"><h2>Variant Type Filter</h2><div class="checkrow" id="variantFilters"></div></article>
  </section>
  <section class="tables">
    <article class="panel">
      <h2>Top Hits</h2>
      <div class="table-wrap">
        <table>
          <thead><tr><th>Rank</th><th>Source</th><th>Locus</th><th>Position</th><th>Type</th><th>Strictness</th><th>Target</th><th>Off-target</th><th>Score</th></tr></thead>
          <tbody id="hitRows"></tbody>
        </table>
      </div>
      <div class="detail" id="hitDetail">Select a hit row for details.</div>
    </article>
    <div class="side">
      <article class="panel">
        <h2>Candidate Regions</h2>
        <div class="table-wrap">
          <table>
            <thead><tr><th>Rank</th><th>Source</th><th>Region</th><th>Position</th><th>Loci</th><th>Score</th></tr></thead>
            <tbody id="regionRows"></tbody>
          </table>
        </div>
      </article>
      <article class="panel provenance">
        <h2>QC And Provenance</h2>
        <div id="qcRows"></div>
      </article>
    </div>
  </section>
</main>
<footer class="site-footer"><div>{FOOTER_HTML}</div></footer>
<script>
const DATA = {payload};
const S = DATA.summary;
const state = {{ source: "all", query: "", minScore: 0, strictness: new Set(), variants: new Set() }};
const colors = ["#28745f", "#2e6f9e", "#b7832f", "#b64b52", "#6a5d9f", "#4f7f86", "#8b6f47", "#657a56"];
function esc(value) {{
  return String(value ?? "").replace(/[&<>"']/g, ch => ({{"&":"&amp;","<":"&lt;",">":"&gt;","\\"":"&quot;","'":"&#39;"}}[ch]));
}}
function fmt(value) {{ return Number(value || 0).toLocaleString(); }}
function score(value) {{ return Number(value || 0).toFixed(3); }}
function sourcePayload() {{
  if (state.source === "all") return null;
  return DATA.sources.find(src => src.key === state.source) || null;
}}
function rowsFor(kind) {{
  const src = sourcePayload();
  if (!src) return DATA.aggregate[kind] || [];
  return src[kind] || [];
}}
function allSourceOptions() {{
  const select = document.getElementById("sourceSelect");
  select.innerHTML = '<option value="all">All sources</option>';
  for (const src of DATA.sources) {{
    const opt = document.createElement("option");
    opt.value = src.key;
    opt.textContent = src.label;
    select.appendChild(opt);
  }}
}}
function setupChecks(id, rows, setName) {{
  const box = document.getElementById(id);
  box.innerHTML = "";
  for (const row of rows) {{
    const label = document.createElement("label");
    const input = document.createElement("input");
    input.type = "checkbox";
    input.checked = true;
    input.value = row.label;
    state[setName].add(row.label);
    input.addEventListener("change", () => {{
      if (input.checked) state[setName].add(input.value);
      else state[setName].delete(input.value);
      update();
    }});
    label.appendChild(input);
    label.appendChild(document.createTextNode(`${{row.label}} (${{fmt(row.count)}})`));
    box.appendChild(label);
  }}
}}
function setupFilters() {{
  setupChecks("strictnessFilters", DATA.aggregate.strictness_counts, "strictness");
  setupChecks("variantFilters", DATA.aggregate.variant_type_counts, "variants");
}}
function filteredHits() {{
  const q = state.query.trim().toLowerCase();
  return DATA.hits.filter(row => {{
    if (state.source !== "all" && row.source !== state.source) return false;
    if (Number(row.final_score || 0) < state.minScore) return false;
    if (!state.strictness.has(row.strictness_class)) return false;
    if (!state.variants.has(row.variant_type)) return false;
    if (!q) return true;
    return [row.locus_id, row.contig, row.variant_type, row.strictness_class, row.allele_key, row.source_label].join(" ").toLowerCase().includes(q);
  }});
}}
function filteredRegions() {{
  return DATA.regions.filter(row => state.source === "all" || row.source === state.source);
}}
function renderMetrics() {{
  const src = sourcePayload();
  document.getElementById("metricSources").textContent = src ? src.label : fmt(S.source_count);
  document.getElementById("metricHits").textContent = fmt(src ? src.summary.hit_count : S.hit_count);
  document.getElementById("metricRegions").textContent = fmt(src ? src.summary.region_count : S.region_count);
  document.getElementById("metricEmbeddedHits").textContent = fmt(src ? src.summary.embedded_hit_rows : S.embedded_hit_rows);
  document.getElementById("metricEmbeddedRegions").textContent = fmt(src ? src.summary.embedded_region_rows : S.embedded_region_rows);
}}
function drawBars(canvasId, rows, colorOffset = 0) {{
  const canvas = document.getElementById(canvasId);
  const rect = canvas.getBoundingClientRect();
  const dpr = window.devicePixelRatio || 1;
  canvas.width = Math.max(1, Math.floor(rect.width * dpr));
  canvas.height = Math.max(1, Math.floor(rect.height * dpr));
  const ctx = canvas.getContext("2d");
  ctx.scale(dpr, dpr);
  ctx.clearRect(0, 0, rect.width, rect.height);
  const pad = {{left: 42, right: 12, top: 12, bottom: 48}};
  const width = rect.width - pad.left - pad.right;
  const height = rect.height - pad.top - pad.bottom;
  const max = Math.max(1, ...rows.map(row => Number(row.count || 0)));
  ctx.strokeStyle = "#d5ded8";
  ctx.beginPath();
  ctx.moveTo(pad.left, pad.top);
  ctx.lineTo(pad.left, pad.top + height);
  ctx.lineTo(pad.left + width, pad.top + height);
  ctx.stroke();
  const shown = rows.slice(0, 18);
  const gap = 5;
  const barW = shown.length ? Math.max(5, (width - gap * (shown.length - 1)) / shown.length) : width;
  ctx.textAlign = "center";
  ctx.textBaseline = "top";
  ctx.font = "11px system-ui, sans-serif";
  shown.forEach((row, idx) => {{
    const h = Number(row.count || 0) / max * height;
    const x = pad.left + idx * (barW + gap);
    const y = pad.top + height - h;
    ctx.fillStyle = colors[(idx + colorOffset) % colors.length];
    ctx.fillRect(x, y, barW, h);
    ctx.save();
    ctx.translate(x + barW / 2, pad.top + height + 8);
    ctx.rotate(-Math.PI / 5);
    ctx.fillStyle = "#5a6761";
    ctx.fillText(String(row.label).slice(0, 18), 0, 0);
    ctx.restore();
  }});
  ctx.textAlign = "right";
  ctx.textBaseline = "middle";
  ctx.fillStyle = "#5a6761";
  ctx.fillText(fmt(max), pad.left - 6, pad.top + 4);
}}
function renderHits() {{
  const tbody = document.getElementById("hitRows");
  const rows = filteredHits().slice(0, 500);
  tbody.innerHTML = "";
  if (!rows.length) {{
    tbody.innerHTML = '<tr><td colspan="9" class="muted">No embedded hits match the current filters.</td></tr>';
    return;
  }}
  for (const row of rows) {{
    const tr = document.createElement("tr");
    tr.innerHTML = `<td>${{fmt(row.rank)}}</td><td>${{esc(row.source_label)}}</td><td>${{esc(row.locus_id)}}</td><td>${{esc(row.contig)}}:${{fmt(row.start)}}-${{fmt(row.end)}}</td><td><span class="pill">${{esc(row.variant_type)}}</span></td><td>${{esc(row.strictness_class)}}</td><td>${{fmt(row.target_support_n)}}/${{fmt(row.target_total_n)}} miss ${{fmt(row.target_missing_n)}}</td><td>${{fmt(row.offtarget_support_n)}}/${{fmt(row.offtarget_total_n)}} miss ${{fmt(row.offtarget_missing_n)}}</td><td>${{score(row.final_score)}}</td>`;
    tr.addEventListener("click", () => {{
      document.getElementById("hitDetail").innerHTML = `<strong>${{esc(row.locus_id)}}</strong> ${{esc(row.allele_key)}}; discovery=${{score(row.discovery_score)}}, support=${{score(row.support_score)}}, penalty=${{score(row.penalty_score)}}.`;
    }});
    tbody.appendChild(tr);
  }}
}}
function renderRegions() {{
  const tbody = document.getElementById("regionRows");
  const rows = filteredRegions().slice(0, 250);
  tbody.innerHTML = "";
  if (!rows.length) {{
    tbody.innerHTML = '<tr><td colspan="6" class="muted">No embedded regions for the selected source.</td></tr>';
    return;
  }}
  for (const row of rows) {{
    const tr = document.createElement("tr");
    tr.innerHTML = `<td>${{fmt(row.rank)}}</td><td>${{esc(row.source_label)}}</td><td>${{esc(row.region_id)}}</td><td>${{esc(row.contig)}}:${{fmt(row.start)}}-${{fmt(row.end)}}</td><td>${{fmt(row.n_loci)}}</td><td>${{score(row.final_score)}}</td>`;
    tbody.appendChild(tr);
  }}
}}
function renderQc() {{
  const box = document.getElementById("qcRows");
  const src = sourcePayload();
  const sources = src ? [src] : DATA.sources;
  let html = "";
  for (const item of sources) {{
    html += `<p><strong>${{esc(item.label)}}</strong><br><code>${{esc(S.inputs.find(x => x.label === item.label)?.path || "")}}</code></p>`;
    const rows = item.qc.slice(0, 10);
    if (rows.length) {{
      html += '<table><thead><tr><th>Metric</th><th>Value</th></tr></thead><tbody>';
      for (const row of rows) html += `<tr><td>${{esc(row.metric)}}</td><td>${{esc(row.value)}}</td></tr>`;
      html += '</tbody></table>';
    }}
  }}
  if (DATA.compare.match_counts.length) {{
    html += '<p><strong>Compare match classes</strong></p><table><tbody>';
    for (const row of DATA.compare.match_counts) html += `<tr><td>${{esc(row.label)}}</td><td>${{fmt(row.count)}}</td></tr>`;
    html += '</tbody></table>';
  }}
  box.innerHTML = html;
}}
function updateCharts() {{
  drawBars("strictnessChart", rowsFor("strictness_counts"), 0);
  drawBars("variantChart", rowsFor("variant_type_counts"), 1);
  drawBars("scoreChart", rowsFor("score_bins"), 2);
  drawBars("contigChart", rowsFor("contig_counts"), 3);
}}
function update() {{
  renderMetrics();
  updateCharts();
  renderHits();
  renderRegions();
  renderQc();
}}
document.getElementById("title").textContent = S.title;
document.getElementById("subtitle").textContent = S.subtitle;
allSourceOptions();
setupFilters();
const maxScore = Math.max(1, ...DATA.sources.map(src => Number(src.summary.max_score || 0)));
const slider = document.getElementById("scoreSlider");
slider.max = String(Math.ceil(maxScore * 100) / 100);
slider.addEventListener("input", () => {{
  state.minScore = Number(slider.value || 0);
  document.getElementById("scoreLabel").textContent = state.minScore.toFixed(2);
  update();
}});
document.getElementById("sourceSelect").addEventListener("change", event => {{
  state.source = event.target.value;
  update();
}});
document.getElementById("searchBox").addEventListener("input", event => {{
  state.query = event.target.value;
  renderHits();
}});
window.addEventListener("resize", updateCharts);
update();
</script>
</body>
</html>
"""
