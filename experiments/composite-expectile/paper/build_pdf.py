"""Print-typeset variant of the manuscript, for Chromium's PDF renderer.

Reuses build_paper.build() and swaps the screen stylesheet for a print one:
white ground, single column at a fixed measure, running head and folio via
@page, and break rules that keep floats, tables and propositions whole.
"""
import re, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_paper

PRINT_CSS = """
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&family=Literata:opsz,wght@7..72,400;7..72,600&display=swap">
<style>
@page{size:210mm 297mm;margin:24mm 30mm 22mm 30mm;
  @bottom-center{content:counter(page);font-family:serif;font-size:9pt;color:#444}}
:root{--ink:#111;--ink-2:#3d3d3d;--ink-3:#6a6a6a;--rule:#c9c9c9;--rule-strong:#8a8a8a;
  --accent:#1c5d7d;--accent-2:#9c4a2f;--accent-3:#2d6a4a;
  --sans:"IBM Plex Sans",Arial,sans-serif;--serif:"Literata",Georgia,serif;
  --mono:"IBM Plex Mono",monospace}
*{box-sizing:border-box}
html,body{background:#fff}
body{color:var(--ink);font-family:var(--serif);font-size:10.5pt;line-height:1.55;
  text-align:justify;hyphens:auto;-webkit-hyphens:auto}
.page{max-width:none;margin:0;padding:0}
.col{max-width:none}

header.title{padding:0 0 12pt;border-bottom:.8pt solid var(--ink);margin-bottom:16pt;
  text-align:left}
.kicker{font-family:var(--sans);font-size:7.6pt;font-weight:600;letter-spacing:.14em;
  text-transform:uppercase;color:var(--ink-3);margin-bottom:10pt}
h1{font-family:var(--serif);font-weight:600;font-size:19pt;line-height:1.22;letter-spacing:-.005em;
  margin:0 0 9pt;max-width:none;text-align:left}
.sub{font-size:10.6pt;line-height:1.45;color:var(--ink-2);margin:0;text-align:left;font-style:italic}
.abstract{background:none;border:none;border-top:.6pt solid var(--rule);
  border-bottom:.6pt solid var(--rule);padding:11pt 0;margin:16pt 0 9pt;max-width:none;
  font-size:9.4pt;line-height:1.48}
.abstract h2{font-family:var(--sans);font-size:7.6pt;font-weight:600;letter-spacing:.14em;
  text-transform:uppercase;color:var(--ink-3);margin:0 0 7pt}
.keywords{font-family:var(--serif);font-size:9pt;color:var(--ink-2);margin:0 0 4pt;max-width:none}
.keywords b{color:var(--ink);font-weight:600}

h2{font-family:var(--sans);font-weight:600;font-size:11.6pt;line-height:1.25;
  margin:18pt 0 3pt;display:block;break-after:avoid;text-align:left}
h2 .n{font-family:var(--sans);color:var(--ink);margin-right:8pt}
h3{font-family:var(--sans);font-weight:600;font-size:10pt;margin:12pt 0 2pt;
  break-after:avoid;text-align:left}
h3 .n{color:var(--ink-3);margin-right:6pt;font-family:var(--sans)}
p{margin:0 0 7pt;orphans:3;widows:3}
h2+p,h3+p{margin-top:5pt}
a{color:var(--ink);text-decoration:none}

.eq{display:flex;align-items:center;gap:12pt;margin:10pt 0;max-width:none;break-inside:avoid}
.eq .body{flex:1;font-family:var(--mono);font-size:9pt;line-height:1.65;background:none;
  border-left:none;padding:0;text-align:center;white-space:pre;overflow-x:visible}
.eq .tag{font-family:var(--serif);font-size:9.4pt;color:var(--ink);flex:none}
.eq .hl{color:var(--accent-2)}

.prop{max-width:none;background:none;border:none;border-left:1.6pt solid var(--ink);
  padding:2pt 0 2pt 12pt;margin:11pt 0;break-inside:avoid;font-size:9.7pt}
.prop .lab{font-family:var(--sans);font-size:8pt;font-weight:600;letter-spacing:.07em;
  text-transform:uppercase;color:var(--ink);display:block;margin-bottom:5pt}
.note{max-width:none;border-left:.8pt solid var(--rule-strong);padding:2pt 0 2pt 12pt;
  margin:10pt 0;color:var(--ink-2);font-size:9.3pt;break-inside:avoid}

.tablewrap{overflow-x:visible;margin:12pt 0;max-width:none;break-inside:avoid}
table{border-collapse:collapse;font-family:var(--sans);font-size:8.2pt;
  font-variant-numeric:tabular-nums;width:100%}
caption{text-align:left;font-family:var(--serif);font-size:8.8pt;color:var(--ink);
  padding-bottom:6pt;max-width:none;line-height:1.4}
th,td{padding:3.2pt 6pt;text-align:right;white-space:nowrap}
th:first-child,td:first-child{text-align:left;padding-left:0}
thead th{font-weight:600;font-size:7.8pt;color:var(--ink);border-bottom:.8pt solid var(--ink)}
tbody td{border-bottom:.4pt solid var(--rule)}
tbody tr:last-child td{border-bottom:.8pt solid var(--ink)}
tr.ref td{color:var(--ink-2);font-style:italic}
tr.sub td{color:var(--ink-3);font-size:.94em}
td.win{color:var(--accent-3);font-weight:600}
td.lose{color:var(--accent-2)}

figure{margin:13pt 0;max-width:none;break-inside:avoid}
figure img{width:100%;display:block;border:.4pt solid var(--rule);background:#fff}
figcaption{font-family:var(--serif);font-size:8.8pt;line-height:1.42;color:var(--ink);
  margin-top:6pt;max-width:none;text-align:left}
.fignum{font-weight:600}

.refs{max-width:none;font-size:9pt;line-height:1.42;padding-left:0;list-style:none;counter-reset:r}
.refs li{margin-bottom:5pt;color:var(--ink);padding-left:16pt;text-indent:-16pt;
  break-inside:avoid;text-align:left}
.refs li b{font-weight:600}
footer{margin-top:16pt;padding-top:9pt;border-top:.6pt solid var(--rule);
  font-family:var(--sans);font-size:8.2pt;color:var(--ink-3);max-width:none;text-align:left}
ul,ol{max-width:none}
code{font-family:var(--mono);font-size:.88em;background:none;padding:0}
</style>
"""

if __name__ == "__main__":
    html = build_paper.build()
    # swap the screen head (title + fonts + style) for the print one, keep <title>
    head_end = html.index("<div class=\"page\">")
    html = '<title>Composite M-quantile estimation for extreme-value models</title>\n' \
           + PRINT_CSS + html[head_end:]
    out = Path(__file__).resolve().parent / "manuscript-print.html"
    out.write_text(html)
    print("wrote paper/manuscript-print.html  (%.2f MB)" % (len(html) / 1e6))
