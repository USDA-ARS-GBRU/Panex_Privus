"""HTML report generator for privy report.

Converts the Markdown report produced by :func:`~privy.report.markdown.render_markdown_report`
into a self-contained HTML file using the ``Markdown`` library (required dep).
The HTML includes minimal inline CSS for readability without external resources.
"""

from __future__ import annotations

from pathlib import Path


def render_html_report(markdown_path: Path, outdir: Path) -> Path:
    """Convert a Markdown report to a self-contained HTML file.

    Args:
        markdown_path: Path to ``report.md`` produced by
            :func:`~privy.report.markdown.render_markdown_report`.
        outdir: Output directory for ``report.html``.

    Returns:
        Path to the written ``report.html``.
    """
    import markdown  # type: ignore[import-untyped]  # lazy import keeps tests light

    text = markdown_path.read_text(encoding="utf-8")
    body = markdown.markdown(text, extensions=["tables", "fenced_code"])

    title = _extract_title(text)
    out_path = outdir / "report.html"
    out_path.write_text(_html_template(body, title), encoding="utf-8")
    return out_path


def _extract_title(text: str) -> str:
    """Return the text of the first H1 heading, or a generic fallback."""
    for line in text.splitlines():
        if line.startswith("# "):
            return line[2:].strip()
    return "Panex Privus Report"


def _html_template(body: str, title: str) -> str:
    """Wrap an HTML body fragment in a minimal self-contained document."""
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{title}</title>
  <style>
    body {{
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
      max-width: 1200px;
      margin: 2em auto;
      padding: 0 1.5em;
      line-height: 1.6;
      color: #24292e;
    }}
    h1, h2, h3 {{
      border-bottom: 1px solid #eaecef;
      padding-bottom: 0.3em;
    }}
    table {{
      border-collapse: collapse;
      width: 100%;
      margin: 1em 0;
      font-size: 0.9em;
    }}
    th {{
      background: #f6f8fa;
      border: 1px solid #dfe2e5;
      padding: 6px 13px;
      text-align: left;
      font-weight: 600;
    }}
    td {{
      border: 1px solid #dfe2e5;
      padding: 6px 13px;
    }}
    tr:nth-child(even) td {{
      background: #f6f8fa;
    }}
    code {{
      background: #f3f4f6;
      padding: 0.15em 0.4em;
      border-radius: 3px;
      font-size: 0.88em;
      font-family: "SFMono-Regular", Consolas, "Liberation Mono", Menlo, monospace;
    }}
    hr {{
      border: none;
      border-top: 1px solid #eaecef;
      margin: 2em 0;
    }}
    em {{ color: #586069; }}
    ul li {{ margin: 0.3em 0; }}
  </style>
</head>
<body>
{body}
</body>
</html>
"""
