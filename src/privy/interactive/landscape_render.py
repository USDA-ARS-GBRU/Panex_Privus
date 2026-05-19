# ruff: noqa: E501
"""Render self-contained HTML for ``privy interactive --landscape``."""

from __future__ import annotations

import json
from html import escape
from typing import Any


def render_landscape_html(data: dict[str, Any]) -> str:
    """Render a self-contained landscape dashboard HTML document."""
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
.toolbar {{ display:grid; grid-template-columns:repeat(4,minmax(170px,1fr)); gap:12px; align-items:end; margin-bottom:16px; }}
label {{ display:block; margin:0 0 6px; color:var(--muted); font-size:.88rem; font-weight:650; }}
select {{ width:100%; min-height:42px; border:1px solid #c5d0c9; background:var(--paper); border-radius:7px; padding:8px 10px; color:var(--ink); font:inherit; }}
.metrics {{ display:grid; grid-template-columns:repeat(5,minmax(140px,1fr)); gap:10px; margin-bottom:16px; }}
.metric {{ background:var(--paper); border:1px solid var(--line); border-radius:8px; padding:12px; box-shadow:var(--shadow); }}
.metric span {{ display:block; color:var(--muted); font-size:.82rem; }}
.metric strong {{ display:block; margin-top:3px; font-size:1.35rem; }}
.grid {{ display:grid; grid-template-columns:minmax(420px,2fr) minmax(300px,1fr); gap:14px; align-items:start; }}
.panel {{ background:var(--paper); border:1px solid var(--line); border-radius:8px; padding:14px; box-shadow:var(--shadow); }}
.panel h2 {{ margin:0 0 10px; font-size:1rem; }}
.chart {{ width:100%; display:block; }}
#heatmap {{ height:310px; }}
#profile {{ height:220px; }}
#contigChart {{ height:180px; }}
.tables {{ display:grid; grid-template-columns:minmax(360px,2fr) minmax(300px,1fr); gap:14px; margin-top:14px; align-items:start; }}
.table-wrap {{ overflow:auto; max-height:500px; border:1px solid var(--line); border-radius:8px; }}
table {{ width:100%; border-collapse:collapse; min-width:760px; }}
th,td {{ padding:8px 9px; border-bottom:1px solid #e7ece8; text-align:left; vertical-align:top; font-size:.86rem; }}
th {{ position:sticky; top:0; background:#eef3ef; z-index:1; font-weight:750; }}
tr:hover td {{ background:#f8faf7; }}
.side table {{ min-width:480px; }}
.detail {{ margin-top:10px; color:var(--muted); min-height:40px; }}
.pill {{ display:inline-block; padding:2px 7px; border-radius:999px; background:#e9f2ee; color:#1d5f4f; font-size:.78rem; font-weight:700; }}
.provenance {{ margin-top:14px; }}
.provenance code {{ word-break:break-all; }}
.legend {{ display:flex; flex-wrap:wrap; gap:8px; color:var(--muted); font-size:.84rem; margin-top:8px; }}
.swatch {{ width:16px; height:10px; display:inline-block; border-radius:2px; margin-right:4px; vertical-align:middle; }}
@media (max-width:980px) {{ .toolbar,.metrics,.grid,.tables {{ grid-template-columns:1fr; }} table {{ min-width:720px; }} }}
</style>
</head>
<body>
<header>
  <h1 id="title"></h1>
  <p class="subtitle" id="subtitle"></p>
</header>
<main>
  <section class="toolbar" aria-label="Landscape controls">
    <div><label for="contigSelect">Contig</label><select id="contigSelect"></select></div>
    <div><label for="sampleSelect">Sample</label><select id="sampleSelect"></select></div>
    <div><label for="metricSelect">Sample Heatmap Metric</label><select id="metricSelect"></select></div>
    <div><label for="windowMetricSelect">Window Profile Metric</label><select id="windowMetricSelect"></select></div>
  </section>
  <section class="metrics" aria-label="Summary metrics">
    <div class="metric"><span>Windows</span><strong id="metricWindows"></strong></div>
    <div class="metric"><span>Sample windows</span><strong id="metricSampleWindows"></strong></div>
    <div class="metric"><span>Candidate blocks</span><strong id="metricBlocks"></strong></div>
    <div class="metric"><span>Contigs</span><strong id="metricContigs"></strong></div>
    <div class="metric"><span>Samples</span><strong id="metricSamples"></strong></div>
  </section>
  <section class="grid">
    <article class="panel">
      <h2>Sample-By-Window Heatmap</h2>
      <canvas class="chart" id="heatmap"></canvas>
      <div class="legend"><span><span class="swatch" style="background:#eef3ef"></span>low</span><span><span class="swatch" style="background:#2d765f"></span>high</span></div>
      <div class="detail" id="heatmapDetail">Hover over the heatmap for window details.</div>
    </article>
    <article class="panel">
      <h2>Window Profile</h2>
      <canvas class="chart" id="profile"></canvas>
      <div class="detail" id="profileDetail">Profile shows the selected window-level metric across the selected contig.</div>
    </article>
  </section>
  <section class="tables">
    <article class="panel">
      <h2>Candidate Introgression Blocks</h2>
      <div class="table-wrap">
        <table>
          <thead><tr><th>Block</th><th>Sample</th><th>Position</th><th>Donor</th><th>Similarity</th><th>Delta</th><th>Evidence</th></tr></thead>
          <tbody id="blockRows"></tbody>
        </table>
      </div>
      <div class="detail" id="blockDetail">Select a block row for interpretation.</div>
    </article>
    <div class="side">
      <article class="panel">
        <h2>Contig Summary</h2>
        <canvas class="chart" id="contigChart"></canvas>
      </article>
      <article class="panel provenance">
        <h2>QC And Provenance</h2>
        <div id="provenance"></div>
      </article>
    </div>
  </section>
</main>
<script>
const DATA = {payload};
const S = DATA.summary;
const sampleMetrics = [
  ["private_alt_rate", "Private ALT rate"],
  ["nearest_similarity", "Nearest-background similarity"],
  ["missing_rate", "Missing rate"],
  ["nonref_rate", "Non-reference rate"]
];
const windowMetrics = [
  ["target_private_alt_rate", "Target-private ALT rate"],
  ["offtarget_private_alt_rate", "Off-target-private ALT rate"],
  ["density_variants_per_kb", "Variant density per kb"],
  ["target_mean_missing_rate", "Target mean missingness"],
  ["offtarget_mean_missing_rate", "Off-target mean missingness"]
];
const state = {{ contig: "", sample: "all", metric: "private_alt_rate", windowMetric: "target_private_alt_rate" }};
function esc(value) {{ return String(value ?? "").replace(/[&<>"']/g, ch => ({{"&":"&amp;","<":"&lt;",">":"&gt;","\\"":"&quot;","'":"&#39;"}}[ch])); }}
function fmt(value) {{ return Number(value || 0).toLocaleString(); }}
function pct(value) {{ return Number(value || 0).toFixed(3); }}
function setupSelect(id, entries, selected) {{
  const select = document.getElementById(id);
  select.innerHTML = "";
  for (const [value, label] of entries) {{
    const opt = document.createElement("option");
    opt.value = value;
    opt.textContent = label;
    select.appendChild(opt);
  }}
  select.value = selected;
}}
function windowsForContig() {{ return DATA.windows.filter(row => row.contig === state.contig); }}
function samplesForContig() {{
  const rows = DATA.sample_windows.filter(row => row.contig === state.contig);
  return [...new Set(rows.map(row => row.sample))].sort();
}}
function sampleRows() {{
  return DATA.sample_windows.filter(row => row.contig === state.contig && (state.sample === "all" || row.sample === state.sample));
}}
function color(value, max) {{
  const t = Math.max(0, Math.min(1, max ? value / max : 0));
  const lo = [238,243,239], hi = [45,118,95];
  const rgb = lo.map((x, i) => Math.round(x + (hi[i] - x) * t));
  return `rgb(${{rgb[0]}},${{rgb[1]}},${{rgb[2]}})`;
}}
function sizeCanvas(canvas) {{
  const rect = canvas.getBoundingClientRect();
  const dpr = window.devicePixelRatio || 1;
  canvas.width = Math.max(1, Math.floor(rect.width * dpr));
  canvas.height = Math.max(1, Math.floor(rect.height * dpr));
  const ctx = canvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return {{ ctx, width: rect.width, height: rect.height }};
}}
function drawHeatmap() {{
  const canvas = document.getElementById("heatmap");
  const {{ctx, width, height}} = sizeCanvas(canvas);
  ctx.clearRect(0,0,width,height);
  const windows = windowsForContig();
  const rows = sampleRows();
  const samples = state.sample === "all" ? samplesForContig() : [state.sample];
  const pad = {{left: 118, right: 16, top: 12, bottom: 38}};
  const plotW = width - pad.left - pad.right;
  const plotH = height - pad.top - pad.bottom;
  if (!windows.length || !rows.length || !samples.length) {{
    ctx.fillStyle = "#5d6963"; ctx.fillText("No embedded sample-window rows for this selection.", 18, 32); return;
  }}
  const bySampleWindow = new Map(rows.map(row => [`${{row.sample}}|${{row.window_id}}`, row]));
  const max = Math.max(1e-9, ...rows.map(row => Number(row[state.metric] || 0)));
  const cellW = Math.max(1, plotW / windows.length);
  const cellH = Math.max(14, plotH / samples.length);
  ctx.font = "12px system-ui, sans-serif";
  ctx.textAlign = "right"; ctx.textBaseline = "middle";
  samples.forEach((sample, yidx) => {{
    const y = pad.top + yidx * cellH;
    ctx.fillStyle = "#5d6963";
    ctx.fillText(sample, pad.left - 8, y + cellH / 2);
    windows.forEach((win, xidx) => {{
      const row = bySampleWindow.get(`${{sample}}|${{win.window_id}}`);
      const value = row ? Number(row[state.metric] || 0) : 0;
      ctx.fillStyle = color(value, max);
      ctx.fillRect(pad.left + xidx * cellW, y, Math.ceil(cellW), Math.max(1, cellH - 1));
    }});
  }});
  ctx.strokeStyle = "#d4ddd7"; ctx.strokeRect(pad.left, pad.top, plotW, samples.length * cellH);
  ctx.textAlign = "center"; ctx.textBaseline = "top"; ctx.fillStyle = "#5d6963";
  ctx.fillText(`${{state.contig}} windows (${{windows.length}} embedded)`, pad.left + plotW / 2, pad.top + samples.length * cellH + 10);
}}
function drawProfile() {{
  const canvas = document.getElementById("profile");
  const {{ctx, width, height}} = sizeCanvas(canvas);
  ctx.clearRect(0,0,width,height);
  const rows = windowsForContig();
  const pad = {{left: 44, right: 14, top: 14, bottom: 34}};
  const plotW = width - pad.left - pad.right;
  const plotH = height - pad.top - pad.bottom;
  if (!rows.length) return;
  const values = rows.map(row => Number(row[state.windowMetric] || 0));
  const max = Math.max(1e-9, ...values);
  ctx.strokeStyle = "#d4ddd7"; ctx.beginPath(); ctx.moveTo(pad.left, pad.top); ctx.lineTo(pad.left, pad.top + plotH); ctx.lineTo(pad.left + plotW, pad.top + plotH); ctx.stroke();
  ctx.strokeStyle = "#2d6f9f"; ctx.lineWidth = 2; ctx.beginPath();
  values.forEach((value, idx) => {{
    const x = pad.left + (rows.length === 1 ? 0 : idx / (rows.length - 1) * plotW);
    const y = pad.top + plotH - value / max * plotH;
    if (idx === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
  }});
  ctx.stroke();
  ctx.fillStyle = "#5d6963"; ctx.font = "12px system-ui, sans-serif"; ctx.textAlign = "right"; ctx.fillText(max.toFixed(3), pad.left - 6, pad.top + 4);
  ctx.textAlign = "center"; ctx.fillText(state.contig, pad.left + plotW / 2, pad.top + plotH + 14);
}}
function drawContigChart() {{
  const canvas = document.getElementById("contigChart");
  const {{ctx, width, height}} = sizeCanvas(canvas);
  ctx.clearRect(0,0,width,height);
  const rows = DATA.contig_summary.slice(0, 24);
  const pad = {{left: 42, right: 12, top: 10, bottom: 48}};
  const plotW = width - pad.left - pad.right;
  const plotH = height - pad.top - pad.bottom;
  const max = Math.max(1, ...rows.map(row => Number(row.n_windows || 0)));
  const gap = 4, barW = rows.length ? Math.max(4, (plotW - gap * (rows.length - 1)) / rows.length) : plotW;
  rows.forEach((row, idx) => {{
    const h = Number(row.n_windows || 0) / max * plotH;
    const x = pad.left + idx * (barW + gap), y = pad.top + plotH - h;
    ctx.fillStyle = row.contig === state.contig ? "#b8842f" : "#2d765f";
    ctx.fillRect(x, y, barW, h);
    ctx.save(); ctx.translate(x + barW / 2, pad.top + plotH + 8); ctx.rotate(-Math.PI / 5);
    ctx.fillStyle = "#5d6963"; ctx.font = "11px system-ui, sans-serif"; ctx.textAlign = "center"; ctx.fillText(row.contig, 0, 0); ctx.restore();
  }});
}}
function renderBlocks() {{
  const tbody = document.getElementById("blockRows");
  const rows = DATA.candidate_blocks.filter(row => row.contig === state.contig && (state.sample === "all" || row.sample === state.sample)).slice(0, 300);
  tbody.innerHTML = "";
  if (!rows.length) {{ tbody.innerHTML = '<tr><td colspan="7">No embedded candidate blocks for this selection.</td></tr>'; return; }}
  for (const row of rows) {{
    const tr = document.createElement("tr");
    tr.innerHTML = `<td>${{esc(row.block_id)}}</td><td>${{esc(row.sample)}}</td><td>${{esc(row.contig)}}:${{fmt(row.start)}}-${{fmt(row.end)}}</td><td>${{esc(row.candidate_donor)}} <span class="pill">${{esc(row.candidate_donor_role)}}</span></td><td>${{pct(row.mean_donor_similarity)}}</td><td>${{pct(row.mean_similarity_delta)}}</td><td>${{esc(row.evidence_class)}}</td>`;
    tr.addEventListener("click", () => {{
      document.getElementById("blockDetail").textContent = `${{row.block_id}}: ${{row.interpretation || "Candidate block from landscape local-background evidence."}}`;
    }});
    tbody.appendChild(tr);
  }}
}}
function renderProvenance() {{
  let html = `<p><code>${{esc(S.landscape_dir)}}</code></p>`;
  html += '<table><tbody>';
  for (const row of DATA.filter_summary.slice(0, 14)) html += `<tr><td>${{esc(row.metric)}}</td><td>${{esc(row.value)}}</td></tr>`;
  html += '</tbody></table>';
  if (DATA.similarity_pairs.length) {{
    html += '<p><strong>Top pairwise mean similarities</strong></p><table><tbody>';
    for (const row of DATA.similarity_pairs.slice(0, 8)) html += `<tr><td>${{esc(row.pair)}}</td><td>${{pct(row.mean_similarity)}}</td></tr>`;
    html += '</tbody></table>';
  }}
  document.getElementById("provenance").innerHTML = html;
}}
function renderMetrics() {{
  document.getElementById("metricWindows").textContent = fmt(S.window_count);
  document.getElementById("metricSampleWindows").textContent = fmt(S.sample_window_count);
  document.getElementById("metricBlocks").textContent = fmt(S.candidate_block_count);
  document.getElementById("metricContigs").textContent = fmt(S.contigs.length);
  document.getElementById("metricSamples").textContent = fmt(S.samples.length);
}}
function updateSampleOptions() {{
  const current = state.sample;
  const samples = samplesForContig();
  setupSelect("sampleSelect", [["all", "All samples"], ...samples.map(sample => [sample, sample])], samples.includes(current) ? current : "all");
  state.sample = document.getElementById("sampleSelect").value;
}}
function update() {{
  updateSampleOptions();
  drawHeatmap();
  drawProfile();
  drawContigChart();
  renderBlocks();
  renderProvenance();
}}
document.getElementById("title").textContent = S.title;
document.getElementById("subtitle").textContent = S.subtitle;
state.contig = S.contigs[0] || "";
setupSelect("contigSelect", S.contigs.map(contig => [contig, contig]), state.contig);
setupSelect("metricSelect", sampleMetrics, state.metric);
setupSelect("windowMetricSelect", windowMetrics, state.windowMetric);
renderMetrics();
document.getElementById("contigSelect").addEventListener("change", e => {{ state.contig = e.target.value; update(); }});
document.getElementById("sampleSelect").addEventListener("change", e => {{ state.sample = e.target.value; drawHeatmap(); renderBlocks(); }});
document.getElementById("metricSelect").addEventListener("change", e => {{ state.metric = e.target.value; drawHeatmap(); }});
document.getElementById("windowMetricSelect").addEventListener("change", e => {{ state.windowMetric = e.target.value; drawProfile(); }});
document.getElementById("heatmap").addEventListener("mousemove", e => {{
  const rows = sampleRows(); const wins = windowsForContig(); if (!rows.length || !wins.length) return;
  document.getElementById("heatmapDetail").textContent = `${{state.metric}} across ${{wins.length}} embedded windows for ${{state.sample === "all" ? "all samples" : state.sample}}.`;
}});
window.addEventListener("resize", () => {{ drawHeatmap(); drawProfile(); drawContigChart(); }});
update();
</script>
</body>
</html>
"""
