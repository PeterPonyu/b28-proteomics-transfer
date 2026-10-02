#!/usr/bin/env python3
"""Geometric PDF checks and complete raster inventory, not a visual-review substitute."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import fitz
from PIL import Image

CAP = Path(__file__).resolve().parents[2]
MAN = CAP/"manuscript"
OUT = CAP/"output/verification"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--render",action="store_true",help="render all pages and figures with Poppler")
    args=parser.parse_args()
    paths=[MAN/"manuscript.pdf",MAN/"SI.pdf",*sorted((MAN/"figures").glob("b28_*.pdf"))]
    render=OUT/"rendered-poppler"
    render.mkdir(parents=True, exist_ok=True)
    records=[]
    errors=[]
    for path in paths:
        doc=fitz.open(path)
        pages=[]
        for number,page in enumerate(doc,1):
            spans=[s for b in page.get_text("dict")["blocks"] for line in b.get("lines",[])
                   for s in line["spans"] if s["text"].strip()]
            outside=[]
            for span in spans:
                r=fitz.Rect(span["bbox"])
                if r.x0 < -1 or r.y0 < -1 or r.x1 > page.rect.width+1 or r.y1 > page.rect.height+1:
                    outside.append(span["text"])
            overlaps=[]
            if path.parent.name=="figures":
                for i,a in enumerate(spans):
                    rect_a=fitz.Rect(a["bbox"])
                    for b in spans[i+1:]:
                        rect_b=fitz.Rect(b["bbox"])
                        area=min(rect_a.get_area(),rect_b.get_area())
                        if area>0 and (rect_a & rect_b).get_area() > .2*area:
                            overlaps.append([a["text"],b["text"]])
            tiny=[s["text"] for s in spans if s["size"]<6 and len(s["text"].strip())>2]
            if outside or overlaps or tiny:
                errors.append(dict(file=str(path.relative_to(CAP)),page=number,outside=outside,
                                   potential_text_overlaps=overlaps,body_text_under_6pt=tiny))
            pages.append(dict(page=number,n_text_spans=len(spans),outside_page=len(outside),
                              potential_text_overlaps=len(overlaps),text_under_6pt=len(tiny)))
        if args.render:
            subprocess.run(["pdftoppm","-r","110","-png",str(path),str(render/path.stem)],check=True,
                           stdout=subprocess.DEVNULL)
        rasters=sorted(render.glob(path.stem+"-*.png"))
        if args.render:
            assert len(rasters)==len(doc),f"Incomplete page rendering: {path}"
        raster_checks = []
        for raster in rasters:
            with Image.open(raster) as image:
                ink = image.convert("L").point(lambda value: 255 if value < 245 else 0)
                box = ink.getbbox()
                touches_edge = box is not None and (box[0] <= 1 or box[1] <= 1 or
                               box[2] >= image.width-1 or box[3] >= image.height-1)
                raster_checks.append(dict(path=str(raster.relative_to(CAP)), ink_bbox=box,
                                          width=image.width, height=image.height,
                                          touches_page_edge=touches_edge))
                if box is None or touches_edge:
                    errors.append(dict(file=str(path.relative_to(CAP)), raster=raster.name,
                                       blank=box is None, ink_touches_page_edge=touches_edge))
        records.append(dict(path=str(path.relative_to(CAP)),sha256=sha(path),pages=pages,
                            rendered_files=[dict(path=str(p.relative_to(CAP)),sha256=sha(p)) for p in rasters],
                            raster_clipping_screen=raster_checks))
    report=dict(status="automated geometry and complete rendering passed" if not errors else "layout defects detected",
                rendered_in_this_run=args.render, visual_review_completed=False,
                scientific_ready_attested=False, records=records, errors=errors,
                limitations="Text boxes and raster edge-ink checks cannot detect every internal clip, collision, poor visual balance, missing stroke, or misleading graphical encoding")
    (OUT/"layout-QA.json").write_text(json.dumps(report,indent=2)+"\n")
    print(f"Layout: {sum(len(r['pages']) for r in records)} PDF pages across {len(records)} files, {len(errors)} geometric/text defects")
    if errors:
        raise SystemExit(json.dumps(errors,indent=2))


if __name__ == "__main__":
    main()
