#!/usr/bin/env python3
"""Current B28 TOC and graphical abstract, from identity-bound corrected tables.

Only b28_TOC / b28_graphical_abstract assets are written. No fits, table export,
manuscript compilation or method rebinding. All drawings are original vector
schematics; the fold cartoons are not actual sample assignments. Editable SVG
text requires a local licensed Arial; the separate display SVG is outlined.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
import xml.etree.ElementTree as ET

import fitz

# Use a disposable writable cache, not a user's global Matplotlib directory.
if "MPLCONFIGDIR" not in os.environ:
    os.environ["MPLCONFIGDIR"] = tempfile.mkdtemp(prefix="b28-ga-mpl-")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Circle, FancyArrowPatch, Polygon, Rectangle
from matplotlib.text import Text


MAN = Path(__file__).resolve().parents[1]
CAP = MAN.parent
DATA = MAN / "figures" / "plotdata"
OUT = MAN / "figures"
NS = "{http://www.w3.org/2000/svg}"
INK, MUTED, RULE = "#1E2935", "#455868", "#CED7DD"
BLUE, HOLD, PALE = "#0072B2", "#B34D00", "#ECF3F7"
PROTOCOLS = ("random_5fold", "patient_gkf", "logo_xgb")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def corrected_data() -> tuple[dict, dict[str, str]]:
    """Refuse stale plot-data; read only corrected, hash-bound cells."""
    contract_path = MAN / "integration-data.json"
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    if contract.get("schema") != "b28-original-manuscript-integration-1":
        raise RuntimeError("Graphical abstract requires the corrected B28 integration contract")
    if contract.get("scientific_ready_attested") is not False:
        raise RuntimeError("Graphical abstract must not inherit a scientific-ready attestation")
    names = ("b28_pair_counts.csv", "b28_protocol_summary.csv", "b28_ladder.csv",
             "b28_riskcoverage.csv", "b28_riskcoverage_cohort.csv")
    hashes = {"manuscript/integration-data.json": sha(contract_path)}
    tables = {}
    for name in names:
        path = DATA / name
        relative = path.relative_to(CAP).as_posix()
        current = sha(path)
        if contract["generated_sha256"].get(relative) != current:
            raise RuntimeError(f"Corrected table identity mismatch: {relative}")
        hashes[relative] = current
        tables[name] = read_csv(path)
    pairs = tables["b28_pair_counts.csv"]
    totals = {key: sum(int(row[key]) for row in pairs)
              for key in ("n_samples", "n_patients", "n_tumor", "n_normal")}
    if len(pairs) != 6 or totals != dict(n_samples=1191, n_patients=719, n_tumor=689, n_normal=502):
        raise RuntimeError(f"Unexpected corrected population: {totals}")
    if (contract["samples"], contract["patients"], contract["primary_accessions"]) != (1191, 719, 8906):
        raise RuntimeError("Population/accession metadata disagrees with corrected tables")
    abundance = [row for row in tables["b28_ladder.csv"]
                 if row["arm"] == "abundance (8906 accessions)"]
    if len(abundance) != 1:
        raise RuntimeError("No unique accession-aligned abundance result")
    protocol_rows = tables["b28_protocol_summary.csv"]
    pooled = {}
    for key in PROTOCOLS:
        rows = [row for row in protocol_rows if row["protocol"] == key and row["fold_id"] == "pooled"]
        folds = [row for row in protocol_rows if row["protocol"] == key and row["fold_id"] != "pooled"]
        if len(rows) != 1 or len(folds) != (6 if key == "logo_xgb" else 5):
            raise RuntimeError(f"Unexpected corrected protocol inventory: {key}")
        row = rows[0]
        if int(row["n"]) != totals["n_samples"] or sum(int(r["n"]) for r in folds) != totals["n_samples"]:
            raise RuntimeError(f"Protocol population mismatch: {key}")
        pooled[key] = dict(auc=float(row["pooled_auc"]), ece=float(row["pooled_ece"]))
        if not all(math.isfinite(value) and 0 <= value <= 1 for value in pooled[key].values()):
            raise RuntimeError(f"Invalid corrected metric: {key}")
    if not math.isclose(pooled["logo_xgb"]["auc"], float(abundance[0]["auc"]), abs_tol=1e-12):
        raise RuntimeError("LOGO pooled AUROC differs between corrected tables")
    risk = tables["b28_riskcoverage.csv"]
    low = [row for row in risk if math.isclose(float(row["coverage"]), .1)]
    full = [row for row in risk if math.isclose(float(row["coverage"]), 1)]
    mix = [row for row in tables["b28_riskcoverage_cohort.csv"]
           if math.isclose(float(row["coverage"]), .1)]
    if len(low) != 1 or len(full) != 1 or len(mix) != 1:
        raise RuntimeError("Coverage endpoints are not uniquely defined")
    retained = int(mix[0]["n_retained"])
    kept_correct = float(low[0]["retained_acc_conf"]) * retained
    easy_counts = [float(mix[0][key]) * retained for key in ("Luad", "Lscc")]
    if any(not math.isclose(value, round(value), abs_tol=1e-10) for value in [kept_correct, *easy_counts]):
        raise RuntimeError("Retained-set fractions do not recover stored integer counts")
    if retained != math.floor(totals["n_samples"] * .1) or kept_correct != retained:
        raise RuntimeError("Unexpected corrected 10%-coverage endpoint")
    if not math.isclose(float(full[0]["retained_acc_conf"]), float(full[0]["full_acc"]), abs_tol=1e-12):
        raise RuntimeError("Full-coverage accuracy cells disagree")
    return dict(samples=totals["n_samples"], patients=totals["n_patients"],
                accessions=contract["primary_accessions"], cohorts=len(pairs), protocols=pooled,
                logo_macro_auc=float(abundance[0]["macro_auc"]),
                full_accuracy=float(full[0]["full_acc"]),
                retained_accuracy=float(low[0]["retained_acc_conf"]), retained=retained,
                retained_correct=round(kept_correct), easy_retained=sum(round(v) for v in easy_counts)), hashes


def configure_fonts() -> dict:
    faces = {}
    for weight in ("normal", "bold"):
        properties = font_manager.FontProperties(family="Arial", style="normal", weight=weight)
        path = Path(font_manager.findfont(properties, fallback_to_default=False))
        font = font_manager.get_font(str(path))
        if font.family_name != "Arial":
            raise RuntimeError(f"Arial fallback prohibited: {path} resolves to {font.family_name}")
        if weight == "bold" and "bold" not in font.style_name.lower():
            raise RuntimeError("Arial panel-label bold face is missing")
        faces[weight] = dict(family=font.family_name, style=font.style_name,
                             local_path=str(path), sha256=sha(path))
    plt.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial"],
                         "font.size": 8.4, "font.weight": "normal", "text.usetex": False,
                         "pdf.fonttype": 42, "pdf.use14corefonts": False,
                         "svg.fonttype": "none", "savefig.facecolor": "white"})
    return faces


class Plate:
    """Point-coordinate drawing with explicit positions and a text-bounds gate."""
    def __init__(self, width: float, height: float):
        self.width, self.height = width * 72, height * 72
        self.fig = plt.figure(figsize=(width, height), dpi=144)
        self.ax = self.fig.add_axes([0, 0, 1, 1])
        self.ax.set(xlim=(0, self.width), ylim=(self.height, 0))
        self.ax.axis("off")
        self.labels = []

    def text(self, x, y, value, *, size=8.4, color=INK, align="left", tag=None, maxwidth=None):
        artist = self.ax.text(x, y, value, fontsize=size, color=color, ha=align, va="baseline",
                              fontfamily="Arial", fontweight="bold" if tag else "normal",
                              parse_math=False, clip_on=False)
        artist.set_gid(f"panel-label-{tag}" if tag else f"label-{len(self.labels):03d}")
        self.labels.append((artist, maxwidth))
        return artist

    def line(self, x1, y1, x2, y2, *, color=RULE, width=.7, style="-"):
        self.ax.plot([x1, x2], [y1, y2], color=color, linewidth=width, linestyle=style,
                     solid_capstyle="round", clip_on=False)

    def arrow(self, x1, y1, x2, y2, *, color=MUTED):
        self.ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                         mutation_scale=7, linewidth=.8, color=color,
                                         shrinkA=0, shrinkB=0))

    def rect(self, x, y, w, h, *, fill="white", edge=RULE, width=.7):
        self.ax.add_patch(Rectangle((x, y), w, h, facecolor=fill, edgecolor=edge, linewidth=width))

    def check_text(self) -> dict:
        self.fig.canvas.draw()
        renderer = self.fig.canvas.get_renderer()
        boxes, sizes = [], []
        for label, maximum in self.labels:
            font_path = font_manager.findfont(label.get_fontproperties(), fallback_to_default=False)
            font = font_manager.get_font(font_path)
            if font.family_name != "Arial":
                raise RuntimeError(f"Unexpected font for {label.get_text()!r}")
            missing = sorted({ord(c) for c in label.get_text() if not c.isspace()
                              and ord(c) not in font.get_charmap()})
            if missing:
                raise RuntimeError(f"Arial missing glyphs: {missing}: {label.get_text()!r}")
            if label.get_fontweight() != "normal" and not label.get_gid().startswith("panel-label-"):
                raise RuntimeError(f"Only panel labels may be bold: {label.get_text()!r}")
            box = label.get_window_extent(renderer)
            bounds = [value * 72 / self.fig.dpi for value in (box.x0, box.y0, box.x1, box.y1)]
            if min(bounds[:2]) < 3 or bounds[2] > self.width - 3 or bounds[3] > self.height - 3:
                raise RuntimeError(f"Text outside safe canvas: {label.get_text()!r}: {bounds}")
            if maximum is not None and bounds[2] - bounds[0] > maximum:
                raise RuntimeError(f"Text exceeds allocated width: {label.get_text()!r}")
            boxes.append((label.get_text(), bounds))
            sizes.append(label.get_fontsize())
        for index, (name, box) in enumerate(boxes):
            for other_name, other in boxes[index + 1:]:
                if min(box[2], other[2]) - max(box[0], other[0]) > .25 and min(box[3], other[3]) - max(box[1], other[1]) > .25:
                    raise RuntimeError(f"Text overlap: {name!r} / {other_name!r}")
        return dict(text_labels=len(boxes), minimum_pt=min(sizes), maximum_pt=max(sizes),
                    font_fallback=False, missing_glyphs=False, text_overlap=False,
                    safe_canvas_margin_pt=3, only_panel_labels_bold=True)


def patient_icon(p: Plate, x, y):
    p.ax.add_patch(Circle((x + 8, y + 5), 3.5, fill=False, edgecolor=BLUE, linewidth=1))
    p.line(x + 8, y + 9, x + 8, y + 17, color=BLUE, width=1)
    p.line(x + 3, y + 12, x + 13, y + 12, color=BLUE, width=1)
    for xx in (x, x + 12):
        p.rect(xx, y + 22, 7, 5, fill=PALE, edge=BLUE)
        p.line(x + 8, y + 17, xx + 3.5, y + 21, color=BLUE)


def matrix_icon(p: Plate, x, y):
    for row in range(3):
        for col in range(4):
            p.rect(x + col * 5.1, y + row * 5.1, 4, 4, fill=PALE, edge=BLUE, width=.5)


def accession_icon(p: Plate, x, y):
    for row in range(3):
        p.rect(x, y + row * 6, 7, 4.1, fill=PALE, edge=BLUE, width=.6)
        p.line(x + 7, y + row * 6 + 2, x + 12, y + row * 6 + 2, color=BLUE)
        p.rect(x + 12, y + row * 6, 7, 4.1, fill=PALE, edge=BLUE, width=.6)


def split_icon(p: Plate, x: float, y: float, mode: int, scale=1.0):
    """Three cartoon cohorts, two cartoon patients, two cartoon samples each."""
    w, h, dx, dy = 4.0 * scale, 4.8 * scale, 5.2 * scale, 6.6 * scale
    gap = 4.1 * scale
    held = ({(0, 0, 0), (1, 1, 1), (2, 0, 1)} if mode == 0 else
            {(0, 1, 0), (0, 1, 1), (2, 0, 0), (2, 0, 1)} if mode == 1 else
            {(2, col, row) for col in (0, 1) for row in (0, 1)})
    for cohort in range(3):
        left = x + cohort * (2 * dx + gap)
        for patient in range(2):
            for sample in range(2):
                fill = HOLD if (cohort, patient, sample) in held else BLUE
                p.rect(left + patient * dx, y + sample * dy, w, h, fill=fill, edge="white", width=.25)
            if mode == 1:
                p.rect(left + patient * dx - .6 * scale, y - 1.2 * scale,
                       w + 1.2 * scale, dy + h + 2.4 * scale, fill="none", edge=INK, width=.55)
        if mode == 2:
            p.rect(left - 1.2 * scale, y - 1.6 * scale, dx + w + 2.4 * scale,
                   dy + h + 3.2 * scale, fill="none", edge=INK, width=.55)


def full_abstract(data: dict) -> Plate:
    p = Plate(6.5, 4.1)
    p.text(12, 19, "Changing the held-out unit changes the evaluation", size=12.6)
    p.text(12, 35, "Tumor vs source-defined normal · six CPTAC cohorts · corrected labels, patient keys and protein identity", size=8.0)
    p.line(12, 45, 456, 45)
    for x, tag, title in ((12, "A", "Align the measured objects"),
                          (166, "B", "Change the split unit"), (310, "C", "Compare the readouts")):
        p.text(x, 62, tag, size=11.2, tag=tag)
        p.text(x + 15, 62, title, size=8.7, maxwidth=133)
    p.line(154, 74, 154, 177)
    p.line(296, 74, 296, 177)
    matrix_icon(p, 15, 81)
    p.text(48, 88, f"{data['samples']:,} samples", size=11.4, color=BLUE)
    p.text(48, 103, ".C cored-normal corrected", size=7.8, maxwidth=102)
    patient_icon(p, 16, 115)
    p.text(48, 122, f"{data['patients']:,} patients", size=11.4, color=BLUE)
    p.text(48, 137, "True patient keys", size=7.8)
    accession_icon(p, 16, 153)
    p.text(48, 160, f"{data['accessions']:,} stable ENSP", size=10.8, color=BLUE, maxwidth=103)
    p.text(48, 174, "Aligned protein abundance", size=7.7, maxwidth=102)

    titles = ("Sample CV", "Patient CV", "LOGO")
    details = ("Sample-stratified 5-fold", "Patient-grouped 5-fold", "Leave one cohort out")
    for index, (title, detail) in enumerate(zip(titles, details)):
        y = 87 + index * 34
        p.text(166, y, title, size=9.4)
        p.text(166, y + 13, detail, size=7.5, maxwidth=81)
        split_icon(p, 249, y - 8, index, scale=.9)
        p.arrow(285, y + 3, 303, y + 3)
    p.text(343, 79, "Pooled AUROC", size=7.8, align="center")
    p.text(429, 79, "Pooled ECE", size=7.8, align="center")
    for index, protocol in enumerate(PROTOCOLS):
        row = data["protocols"][protocol]
        y = 99 + index * 34
        p.text(343, y, f"{row['auc']:.4f}", size=11.4, align="center", color=BLUE)
        p.text(429, y, f"{row['ece']:.4f}", size=11.4, align="center",
               color=HOLD if index == 2 else INK)
        if index < 2:
            p.line(310, y + 10, 455, y + 10, width=.45)
    p.text(310, 184, f"LOGO macro AUROC {data['logo_macro_auc']:.4f}", size=7.7)
    p.text(310, 195, "ECE: expected calibration error", size=7.5)
    p.rect(12, 185, 5, 5, fill=BLUE, edge=BLUE)
    p.text(21, 190, "train", size=7.5)
    p.rect(49, 185, 5, 5, fill=HOLD, edge=HOLD)
    p.text(58, 190, "test; grids are split cartoons, not actual folds", size=7.5)
    p.line(12, 199, 456, 199)

    p.text(12, 214, "Confidence-ranked abstention changes who remains", size=9.7)
    # Original vector funnel: a selection operation, not assay imagery.
    p.ax.add_patch(Polygon([(14, 224), (37, 224), (28, 234), (28, 241), (23, 241), (23, 234)],
                           closed=True, facecolor=PALE, edgecolor=BLUE, linewidth=.8))
    p.text(49, 229, "100% coverage", size=8.3)
    p.text(49, 242, f"Accuracy {data['full_accuracy']:.4f}", size=8.3)
    p.arrow(158, 233, 205, 233)
    p.text(218, 229, "10% requested coverage", size=8.3)
    p.text(218, 242, f"Accuracy {data['retained_accuracy']:.4f} ({data['retained_correct']}/{data['retained']})", size=8.3)
    p.text(12, 257, f"{data['easy_retained']}/{data['retained']} retained samples are LUAD + LSCC; the retained cohort mix changes.", size=8.0)
    p.line(12, 263, 456, 263, color=HOLD, width=.8)
    p.text(12, 276, "Cancer type × study confounded; all-six feature intersection uses target availability.",
           size=7.9, color=HOLD)
    p.text(12, 286, "Protocol contrast, not target-blind deployment, a biological partition or a coverage guarantee.",
           size=7.5, color=MUTED)
    return p


def short_toc(data: dict) -> Plate:
    p = Plate(3.25, 1.70)
    p.text(117, 13, "Calibration differs across split units", size=9.3, align="center")
    p.text(117, 27, f"{data['samples']:,} samples · {data['patients']:,} patients · {data['accessions']:,} stable ENSP",
           size=7.5, align="center")
    p.line(13, 35, 221, 35)
    centres = (40, 117, 194)
    for index, (x, title, protocol) in enumerate(zip(centres, ("Sample CV", "Patient CV", "LOGO"), PROTOCOLS)):
        p.arrow(x, 35, x, 41)
        p.text(x, 52, title, size=7.8, align="center")
        split_icon(p, x - 21.9, 58, index, scale=1.1)
        p.text(x, 86, f"ECE {data['protocols'][protocol]['ece']:.4f}", size=8.2, align="center",
               color=HOLD if index == 2 else BLUE)
    p.text(117, 97, "Pooled expected calibration error (ECE); cartoon folds", size=7.5, align="center")
    p.line(10, 103, 224, 103, color=HOLD, width=.75)
    p.text(117, 113, "Cancer type × study confounded; target features known", size=7.5,
           align="center", color=HOLD)
    return p


def pdf_font_report(path: Path) -> dict:
    raw = subprocess.run(["pdffonts", str(path)], check=True, capture_output=True, text=True).stdout
    rows = [line.split() for line in raw.splitlines()[2:] if line.strip()]
    if not rows or any(row[-5] != "yes" or "Type 3" in " ".join(row)
                       or "Arial" not in row[0] for row in rows):
        raise RuntimeError(f"Expected only embedded non-Type3 Arial in {path}: {raw}")
    return dict(all_embedded=True, type3=False,
                font_names=sorted({row[0].split("+")[-1] for row in rows}), raw_pdffonts=raw)


def export(p: Plate, stem: str, data: dict) -> dict:
    bounds = p.check_text()
    targets = {kind: OUT / f"{stem}{suffix}" for kind, suffix in
               (("pdf", ".pdf"), ("png", ".png"), ("editable_svg", ".svg"),
                ("outline_svg", "-display-outline.svg"))}
    svg_metadata = {"Date": None, "Title": "B28 corrected evaluation-protocol graphical abstract",
                    "Description": "Original vector schematics, not actual fold assignments or assay images. "
                                   "Stored corrected data; cancer type and study are confounded. "
                                   "Known target feature availability; no deployment or biological-partition claim."}
    pdf_metadata = {"CreationDate": None, "ModDate": None, "Title": svg_metadata["Title"],
                    "Subject": svg_metadata["Description"], "Author": "B28 manuscript figure workflow"}
    p.fig.savefig(targets["pdf"], metadata=pdf_metadata, bbox_inches=None)
    with plt.rc_context({"svg.fonttype": "none", "svg.hashsalt": stem}):
        p.fig.savefig(targets["editable_svg"], metadata=svg_metadata, bbox_inches=None)
    with plt.rc_context({"svg.fonttype": "path", "svg.hashsalt": stem}):
        p.fig.savefig(targets["outline_svg"], metadata=svg_metadata, bbox_inches=None)
    for key in ("editable_svg", "outline_svg"):
        root = ET.parse(targets[key]).getroot()
        texts, images = root.findall(f".//{NS}text"), root.findall(f".//{NS}image")
        if images or (key == "outline_svg" and texts) or (key == "editable_svg" and not texts):
            raise RuntimeError(f"Wrong text/raster content in {targets[key]}")
        if key == "editable_svg" and len(texts) != bounds["text_labels"]:
            raise RuntimeError("Native SVG lost editable labels")
    fonts = pdf_font_report(targets["pdf"])
    with fitz.open(targets["pdf"]) as document:
        if len(document) != 1 or abs(document[0].rect.width - p.width) > .01 or abs(document[0].rect.height - p.height) > .01:
            raise RuntimeError("PDF dimensions changed")
        page = document[0]
        if page.get_images():
            raise RuntimeError("Graphical abstract contains unexpected raster data")
        text = page.get_text()
        required = [f"{data[key]:,}" for key in ("samples", "patients", "accessions")]
        required += [f"{data['protocols'][key]['ece']:.4f}" for key in PROTOCOLS]
        if stem == "b28_graphical_abstract":
            required += [f"{data['protocols'][key]['auc']:.4f}" for key in PROTOCOLS]
            required += [f"{data['logo_macro_auc']:.4f}", f"{data['full_accuracy']:.4f}",
                         f"{data['retained_correct']}/{data['retained']}", f"{data['easy_retained']}/{data['retained']}"]
        if any(value not in text for value in required):
            raise RuntimeError(f"Bound numbers missing in {targets['pdf']}: {required}")
        page.get_pixmap(matrix=fitz.Matrix(300 / 72, 300 / 72), alpha=False).save(targets["png"])
        for block in page.get_text("dict")["blocks"]:
            for line in block.get("lines", []):
                for span in line["spans"]:
                    if "Bold" in span["font"] and span["text"] not in ("A", "B", "C"):
                        raise RuntimeError(f"Unexpected non-panel bold text: {span}")
    plt.close(p.fig)
    return dict(stem=stem, width_in=p.width / 72, height_in=p.height / 72,
                width_pt=p.width, height_pt=p.height, png_dpi=300,
                typography_and_geometry=bounds, pdf_fonts=fonts,
                outputs_sha256={path.relative_to(CAP).as_posix(): sha(path) for path in targets.values()},
                editable_svg_text_elements=bounds["text_labels"], display_svg_text_elements=0,
                raster_images_in_pdf_and_svg=0, figures_are_original_vector_schematics=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, help="Optional report outside the scientific contract")
    arguments = parser.parse_args()
    data, sources = corrected_data()
    fonts = configure_fonts()
    OUT.mkdir(parents=True, exist_ok=True)
    results = [export(short_toc(data), "b28_TOC", data),
               export(full_abstract(data), "b28_graphical_abstract", data)]
    sources[Path(__file__).resolve().relative_to(CAP).as_posix()] = sha(Path(__file__).resolve())
    sources["manuscript/scripts/R/make_b28_figures.R"] = sha(MAN / "scripts/R/make_b28_figures.R")
    if corrected_data()[1] != {key: value for key, value in sources.items()
                               if not key.startswith("manuscript/scripts/")}:
        raise RuntimeError("Corrected input changed during graphical-abstract generation")
    record = dict(schema="b28-current-graphical-abstract-v1", status="AUTOMATED_CHECKS_PASS_VISUAL_REVIEW_REQUIRED",
                  sources_sha256=sources, bound_numbers=data, local_font_resolution=fonts, assets=results,
                  scientific_ready_attested=False, scientific_nonclaims=[
                      "No target-blind or unknown-cohort deployment validation.",
                      "Cancer type and study effects cannot be separated in this design.",
                      "AUROC increments are not biological signal fractions.",
                      "Confidence ranking is not calibrated confidence or a coverage guarantee.",
                      "Macro and pooled AUROC are different estimands; no confidence interval is invented.",
                      "Split drawings are cartoons, not the actual five-fold/six-cohort partitions."],
                  commercial_font_files_distributed=False, source_tables_modified=False,
                  other_figures_modified=False, manuscript_compiled=False, method_contract_rebound=False)
    if arguments.report:
        arguments.report.parent.mkdir(parents=True, exist_ok=True)
        arguments.report.write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(record, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
