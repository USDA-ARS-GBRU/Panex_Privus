// Panex Privus interactive synteny dashboard (Model C).
// Reads a JSON data blob injected by Privy into #privy-data and renders linked
// riparian + dotplot views with target-private highlighting. Preact + htm + d3-scale.

import { render } from "preact";
import { useState, useMemo } from "preact/hooks";
import { html } from "htm/preact";
import { scaleLinear } from "d3";

const BLOCK_COLOURS = {
  collinear: "#4393c3",
  inversion: "#d6604d",
  translocation: "#9970ab",
  duplication: "#1b7837",
  unaligned: "#bdbdbd",
};
const PRIVATE_COLOUR = "#e7298a";

function loadData() {
  const el = document.getElementById("privy-data");
  try {
    const d = JSON.parse(el.textContent);
    if (d && d.__privy_placeholder__) return null;
    return d;
  } catch (_e) {
    return null;
  }
}

function Ribbon({ block, x, qy, ry, selected, onHover }) {
  const colour = BLOCK_COLOURS[block.block_type] || "#999";
  const mid = (qy + ry) / 2;
  const p =
    `M ${x(block.query_start)},${qy} L ${x(block.query_end)},${qy} ` +
    `C ${x(block.query_end)},${mid} ${x(block.ref_end)},${mid} ${x(block.ref_end)},${ry} ` +
    `L ${x(block.ref_start)},${ry} ` +
    `C ${x(block.ref_start)},${mid} ${x(block.query_start)},${mid} ${x(block.query_start)},${qy} Z`;
  return html`<path
    d=${p}
    fill=${colour}
    fill-opacity=${selected ? 0.85 : 0.42}
    stroke=${selected ? "#111" : "none"}
    stroke-width=${selected ? 1.5 : 0}
    style="cursor:pointer"
    onmouseenter=${() => onHover(block.block_id)}
  />`;
}

function Riparian({ data, x, selected, onHover, privateOnly }) {
  const genomes = data.genomes;
  const ref = data.reference;
  const queries = genomes.filter((g) => g !== ref);
  const rowH = 64;
  const top = 24;
  const trackY = {};
  queries.forEach((g, i) => (trackY[g] = top + i * rowH));
  const refY = top + queries.length * rowH;
  const height = refY + 50;

  const privateSpans = data.regions.filter((r) => r.target_private === true || r.target_private === "True");
  const blocks = privateOnly
    ? data.blocks.filter((b) =>
        privateSpans.some((r) => b.ref_start < r.ref_end && r.ref_start < b.ref_end))
    : data.blocks;

  return html`<svg id="rip-svg" width="100%" viewBox=${`0 0 920 ${height}`} role="img">
    ${privateSpans.map(
      (r) => html`<rect x=${x(r.ref_start)} y=${refY - 4} width=${Math.max(2, x(r.ref_end) - x(r.ref_start))}
        height=${28} fill=${PRIVATE_COLOUR} fill-opacity="0.18" />`
    )}
    ${blocks.map((b) =>
      html`<${Ribbon} block=${b} x=${x} qy=${trackY[b.query_genome] + 14} ry=${refY}
        selected=${b.block_id === selected} onHover=${onHover} />`
    )}
    ${queries.map(
      (g) => html`<g>
        <rect x="40" y=${trackY[g]} width="840" height="14" fill="#333" rx="3" />
        <text x="40" y=${trackY[g] - 4} font-size="11" fill="#333">${g}</text>
      </g>`
    )}
    <rect x="40" y=${refY} width="840" height="14" fill="#000" rx="3" />
    <text x="40" y=${refY + 30} font-size="11" fill="#000">${ref} (reference)</text>
  </svg>`;
}

function Dotplot({ data, selected, onHover }) {
  const size = 340;
  const pad = 36;
  const xs = data.blocks.flatMap((b) => [b.query_start, b.query_end]);
  const ys = data.blocks.flatMap((b) => [b.ref_start, b.ref_end]);
  const sx = scaleLinear().domain([Math.min(0, ...xs), Math.max(1, ...xs)]).range([pad, size - 8]);
  const sy = scaleLinear().domain([Math.min(0, ...ys), Math.max(1, ...ys)]).range([size - pad, 8]);
  return html`<svg width=${size} height=${size} role="img" style="background:#fafafa;border:1px solid #eee">
    ${data.blocks.map((b) => {
      const colour = BLOCK_COLOURS[b.block_type] || "#999";
      const y1 = b.block_type === "inversion" ? b.ref_end : b.ref_start;
      const y2 = b.block_type === "inversion" ? b.ref_start : b.ref_end;
      const on = b.block_id === selected;
      return html`<line x1=${sx(b.query_start)} y1=${sy(y1)} x2=${sx(b.query_end)} y2=${sy(y2)}
        stroke=${colour} stroke-width=${on ? 4 : 2} stroke-opacity=${on ? 1 : 0.8}
        style="cursor:pointer" onmouseenter=${() => onHover(b.block_id)} />`;
    })}
    <text x=${size / 2} y=${size - 6} font-size="10" text-anchor="middle" fill="#666">query →</text>
    <text x="10" y=${size / 2} font-size="10" text-anchor="middle" fill="#666"
      transform=${`rotate(-90 10 ${size / 2})`}>reference →</text>
  </svg>`;
}

function Legend() {
  return html`<div style="display:flex;gap:14px;flex-wrap:wrap;font-size:12px;margin:6px 0">
    ${Object.entries(BLOCK_COLOURS).map(
      ([k, v]) => html`<span><span style=${`display:inline-block;width:12px;height:12px;background:${v};margin-right:4px;border-radius:2px`}></span>${k}</span>`
    )}
    <span><span style=${`display:inline-block;width:12px;height:12px;background:${PRIVATE_COLOUR};opacity:.4;margin-right:4px;border-radius:2px`}></span>target-private region</span>
  </div>`;
}

function Detail({ data, selected }) {
  const b = data.blocks.find((x) => x.block_id === selected);
  if (!b) return html`<div style="color:#888;font-size:12px">Hover a block to inspect it.</div>`;
  return html`<table style="font-size:12px;border-collapse:collapse">
    ${[
      ["block", b.block_id], ["type", b.block_type], ["strand", b.strand],
      ["query", `${b.query_genome}:${b.query_start}-${b.query_end}`],
      ["reference", `${b.ref_genome} ${b.ref_contig}:${b.ref_start}-${b.ref_end}`],
    ].map(([k, v]) => html`<tr><td style="color:#888;padding-right:10px">${k}</td><td>${v}</td></tr>`)}
  </table>`;
}

function exportSvg(id, filename) {
  const svg = document.getElementById(id);
  if (!svg) return;
  const clone = svg.cloneNode(true);
  clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  const text = new XMLSerializer().serializeToString(clone);
  const blob = new Blob([text], { type: "image/svg+xml" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  a.click();
  URL.revokeObjectURL(url);
}

function MicrohapPanel({ microhaplotypes }) {
  if (!microhaplotypes || !microhaplotypes.length) return null;
  const priv = microhaplotypes.filter((m) => m.target_private);
  return html`<div style="margin-top:14px">
    <h3 style="margin:0 0 4px;font-size:14px">Microhaplotypes
      <span style="font-size:11px;color:#666">(${microhaplotypes.length} loci · ${priv.length} target-private)</span>
    </h3>
    <table style="font-size:11px;border-collapse:collapse">
      <tr style="color:#888;text-align:left"><th>locus</th><th>pos</th><th>alleles</th><th>AAF</th><th>private</th></tr>
      ${microhaplotypes.slice(0, 200).map(
        (m) => html`<tr style=${m.target_private ? `background:${PRIVATE_COLOUR}22` : ""}>
          <td style="padding-right:10px">${m.locus_id}</td>
          <td style="padding-right:10px">${m.contig}:${m.start}-${m.end}</td>
          <td style="padding-right:10px">${m.n_alleles}</td>
          <td style="padding-right:10px">${(m.aaf ?? 0).toFixed(3)}</td>
          <td style=${m.target_private ? `color:${PRIVATE_COLOUR};font-weight:600` : "color:#999"}>${m.target_private ? "yes" : "—"}</td>
        </tr>`
      )}
    </table>
  </div>`;
}

function App({ data }) {
  const [selected, setSelected] = useState(null);
  const [privateOnly, setPrivateOnly] = useState(false);
  const x = useMemo(() => {
    const vals = data.blocks.flatMap((b) => [b.query_start, b.query_end, b.ref_start, b.ref_end]);
    return scaleLinear().domain([Math.min(0, ...vals), Math.max(1, ...vals)]).range([40, 880]);
  }, [data]);
  const m = data.meta || {};
  return html`<div style="font-family:DejaVu Sans,Arial,sans-serif;max-width:960px;margin:16px auto;color:#222">
    <h2 style="margin:0 0 2px">Panex Privus — Synteny</h2>
    <div style="font-size:12px;color:#666;margin-bottom:8px">
      reference <b>${data.reference}</b> · ${m.n_blocks ?? data.blocks.length} blocks ·
      ${m.n_regions ?? data.regions.length} regions ·
      <span style=${`color:${PRIVATE_COLOUR}`}>${m.n_target_private ?? 0} target-private</span>
    </div>
    <label style="font-size:12px"><input type="checkbox" checked=${privateOnly}
      onchange=${(e) => setPrivateOnly(e.target.checked)} /> show only blocks in private regions</label>
    <button style="font-size:11px;margin-left:12px;cursor:pointer"
      onclick=${() => exportSvg("rip-svg", "privy_riparian.svg")}>Export SVG</button>
    <${Legend} />
    <${Riparian} data=${data} x=${x} selected=${selected} onHover=${setSelected} privateOnly=${privateOnly} />
    <div style="display:flex;gap:24px;align-items:flex-start;margin-top:8px">
      <${Dotplot} data=${data} selected=${selected} onHover=${setSelected} />
      <${Detail} data=${data} selected=${selected} />
    </div>
    <${MicrohapPanel} microhaplotypes=${data.microhaplotypes} />
    <div style="font-size:10px;color:#aaa;margin-top:14px">Generated by Panex Privus · self-contained, offline</div>
  </div>`;
}

const data = loadData();
const root = document.getElementById("app");
if (!data) {
  render(html`<div style="font-family:sans-serif;margin:24px;color:#888">No data injected.</div>`, root);
} else {
  render(html`<${App} data=${data} />`, root);
}
