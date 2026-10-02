#!/usr/bin/env python3
"""Render complete Article/SI contact sheets; no model or manuscript changes."""
from pathlib import Path
import fitz
from PIL import Image, ImageDraw

CAP = Path(__file__).resolve().parents[2]
OUT = CAP / "output/verification/rendered"
OUT.mkdir(parents=True, exist_ok=True)
for name in ("manuscript", "SI"):
    doc = fitz.open(CAP/f"manuscript/{name}.pdf")
    thumbs=[]
    for i,page in enumerate(doc):
        path=OUT/f"{name}-{i+1:02}.png"
        page.get_pixmap(matrix=fitz.Matrix(1.5,1.5),alpha=False).save(path)
        im=Image.open(path).convert("RGB")
        im.thumbnail((306,396))
        tile=Image.new("RGB",(326,426),"#e5e5e5")
        tile.paste(im,((326-im.width)//2,20))
        ImageDraw.Draw(tile).text((10,5),f"{name} page {i+1}",fill="black")
        thumbs.append(tile)
    for batch in range(0,len(thumbs),6):
        chosen=thumbs[batch:batch+6]
        sheet=Image.new("RGB",(3*326,((len(chosen)+2)//3)*426),"white")
        for j,im in enumerate(chosen):
            sheet.paste(im,((j%3)*326,(j//3)*426))
        sheet.save(OUT/f"{name}-contact-{batch//6+1}.png")
print(OUT)
