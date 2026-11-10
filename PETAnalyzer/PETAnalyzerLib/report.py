"""Portable reporting without external JavaScript, fonts, or plotting dependencies."""

import csv
import html
import json
from pathlib import Path

from .analysis import BANDS, DIAMETERS

SERIES = (
    ("rc_nominal_mean", "Nominal mean", "#007383"),
    ("rc_a50", "A50 mean", "#2255a4"),
    ("rc_max", "Maximum", "#aa4a19"),
    ("rc_peak", "1 mL peak", "#7a3e86"),
)


def plot(result, profile, *, quality=False):
    rows = result["spheres"]
    series = SERIES
    title = "Activity concentration recovery"
    if quality:
        series = (("contrast_recovery_pct", "Nominal contrast", "#007383"),)
        if "planar_offsets_mm" in result:
            series += (("planar_hot_contrast_pct", "Planar contrast", "#2255a4"),
                       ("planar_background_variability_pct", "Planar variability", "#aa4a19"))
        title, profile = "Contrast recovery and variability (%)", None
    values = [row[key] for row in rows for key, _, _ in series if row[key] is not None]
    low, high = min(0, *values), max(110 if quality else 1.5, *values) * 1.04
    x = lambda diameter: 70 + (diameter - 10) * 540 / 27
    y = lambda value: 305 - (value - low) * 250 / (high - low)
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="700" height="410" viewBox="0 0 700 410">',
             '<rect width="700" height="410" fill="white"/>',
             '<g font-family="sans-serif" font-size="13" fill="#202d3b">',
             f'<text x="70" y="25" font-size="19">{title}</text>']
    for index in range(7):
        value = low + index * (high - low) / 6
        parts.extend((f'<path d="M70 {y(value):.2f}H610" stroke="#d5dde5"/>',
                      f'<text x="16" y="{y(value)+4:.2f}">{value:.2f}</text>'))
    if profile:
        for index, (_, _, color) in enumerate(SERIES[1:len(BANDS[profile][10]) + 1]):
            upper = [(x(d), y(BANDS[profile][d][index][1])) for d in DIAMETERS]
            lower = [(x(d), y(BANDS[profile][d][index][0])) for d in reversed(DIAMETERS)]
            points = " ".join(f"{px:.2f},{py:.2f}" for px, py in upper + lower)
            parts.append(f'<polygon points="{points}" fill="{color}" fill-opacity=".17"/>')
        parts.append(f'<text x="70" y="46">Shaded intervals: {profile} FDG published comparison bands</text>')
    for d in DIAMETERS:
        parts.append(f'<text x="{x(d):.2f}" y="326" text-anchor="middle">{d}</text>')
    parts.append('<text x="340" y="350" text-anchor="middle">Nominal sphere diameter (mm)</text>')
    for index, (key, title, color) in enumerate(series):
        previous = None
        for row in rows:
            value = row[key]
            if value is None:
                previous = None
                continue
            point = (x(row["diameter_mm"]), y(value))
            if previous:
                parts.append(f'<path d="M{previous[0]:.2f} {previous[1]:.2f}L{point[0]:.2f} {point[1]:.2f}" stroke="{color}" stroke-width="2"/>')
            parts.append(f'<circle cx="{point[0]:.2f}" cy="{point[1]:.2f}" r="4" fill="{color}"/>')
            previous = point
        parts.append(f'<text x="{70 + 155 * index}" y="385" fill="{color}">{title}</text>')
    return "".join(parts) + "</g></svg>"


def annotate_bands(result, profile):
    for row in result["spheres"]:
        for key, (low, high) in zip(("rc_a50", "rc_max", "rc_peak"), BANDS[profile][row["diameter_mm"]]):
            value = row[key]
            row[key + "_reference_low"] = low
            row[key + "_reference_high"] = high
            row[key + "_reference_comparison"] = (
                "Unavailable" if value is None else "Within" if low <= value <= high else "Outside"
            )


def render(result, provenance, profile):
    svg = plot(result, profile)
    quality_svg = plot(result, None, quality=True)
    rows = result["spheres"]
    fields = (
        ("diameter_mm", "Diameter (mm)"), ("rc_nominal_mean", "Nominal RC"),
        ("rc_a50", "A50 RC"), ("rc_max", "Maximum RC"), ("rc_peak", "Peak RC"),
        ("a50_ml", "A50 (mL)"), ("contrast_recovery_pct", "Contrast (%)"),
    )
    if "planar_offsets_mm" in result:
        fields += (("planar_hot_contrast_pct", "Planar contrast (%)"),
                   ("planar_background_variability_pct", "Planar BV (%)"))
    if profile:
        fields += (("rc_a50_reference_comparison", "A50 band"), ("rc_max_reference_comparison", "Max band"))
        if profile == "EARL2":
            fields += (("rc_peak_reference_comparison", "Peak band"),)
    def formatted(value):
        if value is None:
            return "Unavailable"
        return f"{value:.3f}" if isinstance(value, float) else html.escape(str(value))

    body = "".join("<tr>" + "".join(f"<td>{formatted(row[key])}</td>" for key, _ in fields) + "</tr>"
                   for row in rows)
    headings = "".join(f"<th>{label}</th>" for _, label in fields)
    reference = (f'<p>{profile} FDG comparison intervals: <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC8994003/">'
                 'historical published Table 1</a>. The EARL2 peak limits were under revision in that source. '
                 'The analysis uses local VOIs and a constrained peak search. Accreditation requires the current official '
                 'protocol and independent assessment.</p>') if profile else ""
    lung = ""
    if "lung_residual_pct" in result:
        lung = f"<li>Lung insert residual: {result['lung_residual_pct']:.2f}%</li>"
    metadata = html.escape(json.dumps(provenance, indent=2, allow_nan=False))
    return f"""<!doctype html><html lang="en"><meta charset="utf-8"><title>PET Analyzer report</title>
<style>body{{font:16px system-ui,sans-serif;color:#202d3b;max-width:950px;margin:2em auto;padding:0 1em}}
table{{border-collapse:collapse;width:100%}}td,th{{padding:.55em;border-bottom:1px solid #ccd8df;text-align:left}}
svg{{width:100%;height:auto}}pre{{white-space:pre-wrap}}@media print{{figure,table{{break-inside:avoid}}}}</style>
<h1>IEC hot-sphere phantom analysis</h1><p>Research report · PET Analyzer 0.1.0</p>
<figure>{svg}</figure><figure>{quality_svg}</figure><h2>Sphere measurements</h2>
<table><thead><tr>{headings}</tr></thead><tbody>{body}</tbody></table>
<h2>Background and lung insert</h2><ul>
<li>Measured background: {result['background_mean_kbq_ml']:.3f} kBq/mL</li>
<li>Background bias: {result['background_bias_pct']:.2f}%</li>
<li>Within-ROI background coefficient of variation: {result['background_cv_pct']:.2f}%</li>{lung}</ul>
<p>Contrast uses the nominal 3D sphere mean and a single background VOI. Peak centers are restricted to the nominal sphere.
Optional planar results use 12 locations on each of five native slices nearest to the nominal offsets. Full NEMA NU 2 compliance also depends on the acquisition, ROI placement and other protocol requirements.</p>{reference}
<h2>Preparation and measurement provenance</h2><pre>{metadata}</pre></html>"""


def export(folder, result, provenance, profile):
    """Export to a new directory to keep previous reports intact."""
    folder = Path(folder)
    folder.mkdir(exist_ok=False)
    payload = {"schema_version": 1, "results": result, "provenance": provenance}
    (folder / "report.html").write_text(render(result, provenance, profile), encoding="utf-8")
    (folder / "contrast.svg").write_text(plot(result, None, quality=True), encoding="utf-8")
    (folder / "recovery.svg").write_text(plot(result, profile), encoding="utf-8")
    (folder / "analysis.json").write_text(json.dumps(payload, indent=2, allow_nan=False), encoding="utf-8")
    with (folder / "spheres.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=result["spheres"][0].keys())
        writer.writeheader()
        writer.writerows(result["spheres"])
    return folder / "report.html"
