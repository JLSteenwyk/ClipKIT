import json
from html import escape
from typing import Any

import numpy as np


def _to_list(values: Any) -> list:
    if values is None:
        return []
    if isinstance(values, np.ndarray):
        return values.tolist()
    return list(values)


def _site_class_values(site_classes: Any) -> list[str]:
    values = []
    for item in _to_list(site_classes):
        values.append(getattr(item, "value", str(item)))
    return values


def _infer_sequence_type(msa) -> str:
    nt_chars = set("ACGTUN-?*X")
    observed = set()
    for row in msa.seq_records.tolist():
        for char in row:
            if not char:
                continue
            observed.add(char.upper())
            if len(observed) > 40:
                break
        if len(observed) > 40:
            break

    if observed and observed.issubset(nt_chars):
        return "nt"
    return "aa"


def _preview_payload(msa, max_sequences: int, max_columns: int) -> dict:
    n_sequences, n_columns = msa.seq_records.shape
    seq_count = min(n_sequences, max_sequences)
    col_count = min(n_columns, max_columns)

    col_indices = list(range(col_count))
    trim_set = {int(i) for i in _to_list(msa._site_positions_to_trim)}

    rows = []
    for idx in range(seq_count):
        info = msa.header_info[idx]
        label = (
            info.get("id")
            or info.get("name")
            or info.get("description")
            or f"seq_{idx + 1}"
        )
        seq = "".join(msa.seq_records[idx, :col_count].tolist())
        rows.append({"label": str(label), "sequence": seq})

    return {
        "shown_sequences": seq_count,
        "total_sequences": n_sequences,
        "shown_columns": col_count,
        "total_columns": n_columns,
        "trimmed_preview_columns": [i for i in col_indices if i in trim_set],
        "rows": rows,
    }


def write_trim_plot_report(
    path: str,
    msa,
    mode: str,
    gaps: float,
    sequence_type: str | None = None,
) -> None:
    trimmed_positions = [int(x) for x in _to_list(msa._site_positions_to_trim)]
    alignment_length = int(msa.original_length)
    trimmed_count = len(trimmed_positions)
    kept_count = max(0, alignment_length - trimmed_count)
    trimmed_percent = (
        round((trimmed_count / alignment_length) * 100.0, 3) if alignment_length else 0.0
    )
    resolved_sequence_type = sequence_type or _infer_sequence_type(msa)

    payload = {
        "mode": mode,
        "gaps": gaps,
        "sequence_type": resolved_sequence_type,
        "alignment_length": alignment_length,
        "trimmed_count": trimmed_count,
        "kept_count": kept_count,
        "trimmed_percent": trimmed_percent,
        "site_gappyness": [float(x) for x in _to_list(msa.site_gappyness)],
        "site_entropy": [float(x) for x in _to_list(msa.site_entropy)],
        "site_classification": _site_class_values(msa.site_classification_types),
        "trimmed_positions": trimmed_positions,
        "preview": _preview_payload(msa, max_sequences=25, max_columns=240),
    }
    payload_json = json.dumps(payload)

    html = f"""<!doctype html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\">
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">
  <title>ClipKIT trim report</title>
  <style>
    :root {{
      --ink: #1b1f24;
      --ink-soft: #4a5562;
      --bg: #f4f5f8;
      --card: #ffffff;
      --line: #d8dde6;
      --accent: #005f9e;
      --warm: #c76a00;
      --trim: #b42318;
      --trim-bg: rgba(180, 35, 24, 0.22);
      --radius: 12px;
      --shadow: 0 10px 32px rgba(24, 39, 75, 0.08);
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      color: var(--ink);
      background:
        radial-gradient(1200px 500px at 10% -10%, #dfeefe 0%, rgba(223, 238, 254, 0) 70%),
        radial-gradient(900px 420px at 100% 0%, #ffefe1 0%, rgba(255, 239, 225, 0) 70%),
        var(--bg);
      font-family: \"Avenir Next\", \"Trebuchet MS\", \"Segoe UI\", sans-serif;
    }}
    .page {{
      max-width: 1200px;
      margin: 22px auto 34px;
      padding: 0 16px;
    }}
    .hero {{
      background: linear-gradient(135deg, #0f3552 0%, #174d73 58%, #276089 100%);
      color: #fff;
      border-radius: var(--radius);
      box-shadow: var(--shadow);
      padding: 20px 22px;
      margin-bottom: 16px;
    }}
    .hero h1 {{
      margin: 0 0 8px;
      font-family: \"Palatino Linotype\", Palatino, \"Book Antiqua\", serif;
      font-size: 1.8rem;
      letter-spacing: 0.3px;
    }}
    .hero p {{ margin: 0; line-height: 1.45; color: #e5f2ff; }}
    .cards {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(170px, 1fr));
      gap: 10px;
      margin: 14px 0 16px;
    }}
    .card {{
      background: var(--card);
      border: 1px solid var(--line);
      border-radius: 10px;
      padding: 10px 12px;
      box-shadow: 0 2px 10px rgba(24, 39, 75, 0.04);
    }}
    .card .label {{
      font-size: 0.75rem;
      text-transform: uppercase;
      letter-spacing: 0.05em;
      color: var(--ink-soft);
      margin-bottom: 4px;
    }}
    .card .value {{
      font-size: 1.25rem;
      font-weight: 700;
    }}
    .panel {{
      background: var(--card);
      border: 1px solid var(--line);
      border-radius: var(--radius);
      padding: 14px;
      box-shadow: var(--shadow);
      margin-bottom: 12px;
    }}
    .panel h2 {{
      margin: 0 0 8px;
      font-size: 1.1rem;
      font-family: \"Palatino Linotype\", Palatino, \"Book Antiqua\", serif;
    }}
    .legend {{
      display: flex;
      gap: 14px;
      flex-wrap: wrap;
      color: var(--ink-soft);
      font-size: 0.88rem;
      margin-bottom: 8px;
    }}
    .dot {{
      display: inline-block;
      width: 10px;
      height: 10px;
      margin-right: 6px;
      border-radius: 999px;
      vertical-align: middle;
    }}
    .legend-note {{
      color: var(--ink-soft);
      margin: 0 0 8px;
      font-size: 0.86rem;
    }}
    .controls {{
      display: flex;
      gap: 8px;
      flex-wrap: wrap;
      margin-bottom: 8px;
    }}
    .btn {{
      border: 1px solid #b9c4d3;
      background: #f8fbff;
      color: #163a59;
      border-radius: 8px;
      padding: 6px 10px;
      font-size: 0.82rem;
      font-weight: 600;
      cursor: pointer;
    }}
    .btn:hover {{
      background: #ecf5ff;
    }}
    .track-wrap {{
      overflow-x: auto;
      border: 1px solid var(--line);
      padding: 8px;
      border-radius: 10px;
      background: #fff;
    }}
    canvas {{ display: block; width: 100%; max-width: 1160px; height: 240px; }}
    .tooltip {{
      margin-top: 8px;
      font-size: 0.9rem;
      color: var(--ink-soft);
      min-height: 1.2em;
    }}
    .preview-wrap {{
      overflow-x: auto;
      border: 1px solid var(--line);
      border-radius: 10px;
      padding: 8px;
      margin-top: 8px;
      background: #fff;
    }}
    table {{
      border-collapse: collapse;
      font-family: Menlo, Monaco, Consolas, \"Liberation Mono\", monospace;
      font-size: 12px;
      border-spacing: 0;
    }}
    th, td {{ border: 1px solid #e6ebf1; padding: 2px 4px; white-space: nowrap; }}
    .pos-head {{ text-align: center; min-width: 16px; color: #5a6472; }}
    .seq-id {{
      position: sticky;
      left: 0;
      background: #fff;
      z-index: 1;
      max-width: 360px;
      text-align: left;
      overflow: hidden;
      text-overflow: ellipsis;
    }}
    .trimmed {{ box-shadow: inset 0 0 0 999px var(--trim-bg); }}
    .trimmed-head {{ background: #fde8e8; }}
    .muted {{ color: var(--ink-soft); margin: 0; }}
    footer {{
      margin-top: 12px;
      color: var(--ink-soft);
      font-size: 0.82rem;
      text-align: center;
      padding: 10px 12px 0;
      border-top: 1px solid var(--line);
    }}
    footer a {{
      color: #0b5ea8;
      text-decoration: none;
      border-bottom: 1px solid rgba(11, 94, 168, 0.35);
    }}
    footer a:hover {{
      text-decoration: underline;
    }}
    @media (max-width: 700px) {{
      .page {{ padding: 0 10px; }}
      .hero {{ padding: 16px; }}
      .card .value {{ font-size: 1.05rem; }}
      canvas {{ height: 190px; }}
    }}
  </style>
</head>
<body>
  <div class=\"page\">
    <section class=\"hero\">
      <h1>ClipKIT Trim Report</h1>
      <p>Per-site diagnostics over the original alignment with explicit highlighting of trimmed columns.</p>
    </section>

    <section class=\"cards\">
      <article class=\"card\"><div class=\"label\">Mode</div><div class=\"value\">{escape(mode)}</div></article>
      <article class=\"card\"><div class=\"label\">Sequence Type</div><div class=\"value\">{escape(resolved_sequence_type)}</div></article>
      <article class=\"card\"><div class=\"label\">Threshold</div><div class=\"value\">{gaps}</div></article>
      <article class=\"card\"><div class=\"label\">Alignment Length</div><div class=\"value\">{alignment_length}</div></article>
      <article class=\"card\"><div class=\"label\">Sites Kept</div><div class=\"value\">{kept_count}</div></article>
      <article class=\"card\"><div class=\"label\">Sites Trimmed</div><div class=\"value\">{trimmed_count} ({trimmed_percent}%)</div></article>
    </section>

    <section class=\"panel\">
      <h2>Per-site Tracks</h2>
      <div class=\"legend\">
        <span><span class=\"dot\" style=\"background: var(--accent);\"></span>gappyness (binned bars)</span>
        <span><span class=\"dot\" style=\"background: var(--warm);\"></span>entropy (line)</span>
        <span><span class=\"dot\" style=\"background: var(--trim);\"></span>trimmed columns (red vertical bands)</span>
      </div>
      <div class=\"controls\">
        <button class=\"btn\" id=\"export-track-png\" type=\"button\">Export Per-site Tracks (PNG)</button>
      </div>
      <p class=\"legend-note\">Trimmed columns are highlighted in red here and in the alignment preview below.</p>
      <div class=\"track-wrap\">
        <canvas id=\"track\" width=\"1200\" height=\"240\"></canvas>
      </div>
      <div class=\"tooltip\" id=\"tooltip\">Hover to inspect a site.</div>
    </section>

    <section class=\"panel\">
      <h2>Alignment Preview</h2>
      <div class=\"controls\">
        <button class=\"btn\" id=\"export-preview-png\" type=\"button\">Export Alignment Preview (PNG)</button>
      </div>
      <p class=\"muted\" id=\"preview-meta\"></p>
      <div class=\"preview-wrap\">
        <table id=\"preview-table\"></table>
      </div>
    </section>

    <footer>
      Citation: Steenwyk et al. 2020, PLOS Biology.
      <a href="https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001007" target="_blank" rel="noopener noreferrer">
        doi:10.1371/journal.pbio.3001007
      </a>
    </footer>
  </div>

  <script id=\"payload\" type=\"application/json\">{payload_json}</script>
  <script>
    const data = JSON.parse(document.getElementById("payload").textContent);
    const n = data.site_gappyness.length;
    const trimmed = new Set(data.trimmed_positions);
    const canvas = document.getElementById("track");
    const ctx = canvas.getContext("2d");
    const tooltip = document.getElementById("tooltip");
    const pad = 26;
    const ratio = window.devicePixelRatio || 1;

    const NT_COLORS = {{
      A: "#a5d6a7", C: "#90caf9", G: "#f8bbd0", T: "#ffe082", U: "#ffe082", N: "#e0e0e0"
    }};
    const AA_COLORS = {{
      A: "#c8e6c9", V: "#c8e6c9", L: "#c8e6c9", I: "#c8e6c9", M: "#c8e6c9", P: "#c8e6c9", W: "#c8e6c9", F: "#c8e6c9",
      G: "#fff9c4", S: "#fff9c4", T: "#fff9c4", Y: "#fff9c4", C: "#fff9c4", N: "#fff9c4", Q: "#fff9c4",
      D: "#ffcdd2", E: "#ffcdd2",
      K: "#bbdefb", R: "#bbdefb", H: "#bbdefb"
    }};

    function residueColor(char, sequenceType) {{
      const c = (char || "").toUpperCase();
      if (!c) return "#ffffff";
      if ("-?*X".includes(c)) return "#f0f0f0";
      if (sequenceType === "nt") return NT_COLORS[c] || "#f7f7f7";
      return AA_COLORS[c] || "#f3f4f6";
    }}

    function downloadCanvas(canvasEl, filename) {{
      const link = document.createElement("a");
      link.href = canvasEl.toDataURL("image/png");
      link.download = filename;
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
    }}

    function setupCanvas() {{
      const cssWidth = canvas.clientWidth || 1200;
      const cssHeight = parseInt(getComputedStyle(canvas).height, 10) || 240;
      canvas.width = Math.floor(cssWidth * ratio);
      canvas.height = Math.floor(cssHeight * ratio);
      ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
      return {{ width: cssWidth, height: cssHeight }};
    }}

    function drawTrack() {{
      const dims = setupCanvas();
      const drawW = dims.width - pad * 2;
      const drawH = dims.height - pad * 2;
      ctx.clearRect(0, 0, dims.width, dims.height);
      ctx.fillStyle = "#fbfcff";
      ctx.fillRect(pad, pad, drawW, drawH);
      ctx.strokeStyle = "#d8dde6";
      ctx.strokeRect(pad, pad, drawW, drawH);

      const bins = [];
      for (let x = 0; x < drawW; x++) {{
        const start = Math.floor(x * n / drawW);
        const end = Math.max(start + 1, Math.floor((x + 1) * n / drawW));
        let gMax = 0;
        let eSum = 0;
        let count = 0;
        let anyTrim = false;
        for (let i = start; i < end; i++) {{
          if (i >= n) break;
          gMax = Math.max(gMax, data.site_gappyness[i]);
          eSum += data.site_entropy[i];
          count += 1;
          if (trimmed.has(i)) anyTrim = true;
        }}
        const eAvg = count ? (eSum / count) : 0;
        bins.push({{ gMax, eAvg, anyTrim }});
      }}

      for (let x = 0; x < drawW; x++) {{
        const bin = bins[x];
        if (bin.anyTrim) {{
          ctx.fillStyle = "rgba(180, 35, 24, 0.18)";
          ctx.fillRect(pad + x, pad, 1, drawH);
        }}
        const gy = pad + drawH - Math.round(bin.gMax * drawH);
        ctx.fillStyle = "#005f9e";
        ctx.fillRect(pad + x, gy, 1, Math.max(1, Math.round(bin.gMax * drawH)));
      }}

      ctx.beginPath();
      ctx.strokeStyle = "#c76a00";
      ctx.lineWidth = 1.6;
      for (let x = 0; x < drawW; x++) {{
        const y = pad + drawH - Math.round(bins[x].eAvg * drawH);
        if (x === 0) ctx.moveTo(pad + x, y);
        else ctx.lineTo(pad + x, y);
      }}
      ctx.stroke();

      ctx.fillStyle = "#5a6472";
      ctx.font = "12px Menlo, Monaco, Consolas, monospace";
      ctx.fillText("1.0", 4, pad + 3);
      ctx.fillText("0.5", 4, pad + Math.round(drawH / 2) + 3);
      ctx.fillText("0.0", 4, pad + drawH + 3);
    }}

    canvas.addEventListener("mousemove", (event) => {{
      const rect = canvas.getBoundingClientRect();
      const drawW = (canvas.clientWidth || 1200) - pad * 2;
      const x = Math.max(pad, Math.min(pad + drawW - 1, Math.round(event.clientX - rect.left)));
      const site = Math.min(n - 1, Math.max(0, Math.floor((x - pad) * n / drawW)));
      tooltip.textContent = `Site ${{site + 1}} | gappyness=${{data.site_gappyness[site].toFixed(4)}} | entropy=${{data.site_entropy[site].toFixed(4)}} | class=${{data.site_classification[site]}} | status=${{trimmed.has(site) ? "trimmed" : "kept"}}`;
    }});

    canvas.addEventListener("mouseleave", () => {{
      tooltip.textContent = "Hover to inspect a site.";
    }});

    function buildPreview() {{
      const preview = data.preview;
      const trimmedPreviewCols = new Set(preview.trimmed_preview_columns);
      document.getElementById("preview-meta").textContent =
        `Showing first ${{preview.shown_sequences}}/${{preview.total_sequences}} sequences and first ${{preview.shown_columns}}/${{preview.total_columns}} columns.`;

      const table = document.getElementById("preview-table");
      table.innerHTML = "";
      const headRow = document.createElement("tr");
      const headLabel = document.createElement("th");
      headLabel.className = "seq-id";
      headLabel.textContent = "seq/site";
      headRow.appendChild(headLabel);
      for (let i = 0; i < preview.shown_columns; i++) {{
        const th = document.createElement("th");
        th.className = "pos-head";
        if ((i + 1) % 10 === 0) th.textContent = String(i + 1);
        if (trimmedPreviewCols.has(i)) th.classList.add("trimmed-head");
        headRow.appendChild(th);
      }}
      table.appendChild(headRow);

      for (const row of preview.rows) {{
        const tr = document.createElement("tr");
        const label = document.createElement("td");
        label.className = "seq-id";
        label.textContent = row.label;
        tr.appendChild(label);
        for (let i = 0; i < row.sequence.length; i++) {{
          const td = document.createElement("td");
          const residue = row.sequence[i];
          td.textContent = residue;
          td.style.backgroundColor = residueColor(residue, data.sequence_type);
          if (trimmedPreviewCols.has(i)) td.classList.add("trimmed");
          tr.appendChild(td);
        }}
        table.appendChild(tr);
      }}
    }}

    function renderPreviewToCanvas() {{
      const preview = data.preview;
      const trimmedPreviewCols = new Set(preview.trimmed_preview_columns);
      const cellW = 14;
      const cellH = 16;
      const labelW = 230;
      const headerH = 18;
      const width = labelW + (preview.shown_columns * cellW);
      const height = headerH + (preview.rows.length * cellH);

      const exportCanvas = document.createElement("canvas");
      exportCanvas.width = width;
      exportCanvas.height = height;
      const ectx = exportCanvas.getContext("2d");

      ectx.fillStyle = "#ffffff";
      ectx.fillRect(0, 0, width, height);
      ectx.font = "11px Menlo, Monaco, Consolas, monospace";
      ectx.textBaseline = "middle";

      ectx.fillStyle = "#f8fafc";
      ectx.fillRect(0, 0, labelW, headerH);
      ectx.strokeStyle = "#d8dde6";
      ectx.strokeRect(0, 0, width, height);
      ectx.fillStyle = "#334155";
      ectx.fillText("seq/site", 6, Math.floor(headerH / 2));

      for (let i = 0; i < preview.shown_columns; i++) {{
        const x = labelW + (i * cellW);
        if (trimmedPreviewCols.has(i)) {{
          ectx.fillStyle = "#fde8e8";
          ectx.fillRect(x, 0, cellW, headerH);
        }}
        ectx.strokeStyle = "#e6ebf1";
        ectx.strokeRect(x, 0, cellW, headerH);
        if ((i + 1) % 10 === 0) {{
          ectx.fillStyle = "#5a6472";
          ectx.fillText(String(i + 1), x + 1, Math.floor(headerH / 2));
        }}
      }}

      for (let r = 0; r < preview.rows.length; r++) {{
        const row = preview.rows[r];
        const y = headerH + (r * cellH);
        ectx.fillStyle = "#ffffff";
        ectx.fillRect(0, y, labelW, cellH);
        ectx.strokeStyle = "#e6ebf1";
        ectx.strokeRect(0, y, labelW, cellH);
        ectx.fillStyle = "#334155";
        ectx.fillText(String(row.label).slice(0, 36), 6, y + Math.floor(cellH / 2));

        for (let c = 0; c < row.sequence.length; c++) {{
          const x = labelW + (c * cellW);
          const residue = row.sequence[c];
          ectx.fillStyle = residueColor(residue, data.sequence_type);
          ectx.fillRect(x, y, cellW, cellH);
          if (trimmedPreviewCols.has(c)) {{
            ectx.fillStyle = "rgba(180, 35, 24, 0.20)";
            ectx.fillRect(x, y, cellW, cellH);
          }}
          ectx.strokeStyle = "#e6ebf1";
          ectx.strokeRect(x, y, cellW, cellH);
          ectx.fillStyle = "#111827";
          ectx.fillText(residue, x + 3, y + Math.floor(cellH / 2));
        }}
      }}

      return exportCanvas;
    }}

    drawTrack();
    buildPreview();
    document.getElementById("export-track-png").addEventListener("click", () => {{
      downloadCanvas(canvas, "clipkit_per_site_tracks.png");
    }});
    document.getElementById("export-preview-png").addEventListener("click", () => {{
      const previewCanvas = renderPreviewToCanvas();
      downloadCanvas(previewCanvas, "clipkit_alignment_preview.png");
    }});
    window.addEventListener("resize", drawTrack);
  </script>
</body>
</html>
"""

    with open(path, "w", encoding="utf-8") as handle:
        handle.write(html)
