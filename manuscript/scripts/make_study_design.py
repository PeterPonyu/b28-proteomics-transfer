#!/usr/bin/env python3
"""Study-design schematic for the corrected six-cohort CPTAC composition.

Panel integers are read from current accession-aligned plot-data.
Panel C is a split-unit cartoon, not a real fold assignment and not tissue art.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import pandas as pd
from matplotlib.patches import Rectangle, Circle, FancyArrowPatch

ROOT = Path(__file__).resolve().parents[1]
CSV = ROOT / "figures" / "plotdata" / "b28_pair_counts.csv"
OUT = ROOT / "figures"

TRAIN = "#9CC5DD"
HOLD = "#D55E00"
COLS = ["#0072B2", "#E69F00", "#999999"]
INSIDE_MIN = 15


def _read_counts() -> pd.DataFrame:
    data = pd.read_csv(CSV)
    data = data[data.cohort != "TOTAL"].copy()
    # First CSV row (CCRCC) at the top of the bars.
    return data.iloc[::-1].reset_index(drop=True)


def _label_segments(ax, starts, values, y_pos, fill_index: int) -> None:
    """Print every positive integer. Small segments sit just right of the bar."""
    for y, start, value in zip(y_pos, starts, values):
        n = int(value)
        if n <= 0:
            continue
        if n >= INSIDE_MIN:
            ink = "white" if fill_index == 0 else "black"
            ax.text(
                start + value / 2.0,
                y,
                str(n),
                ha="center",
                va="center",
                fontsize=8.3,
                color=ink,
                clip_on=False,
            )
        else:
            ax.text(
                start + value + 2.8,
                y,
                str(n),
                ha="left",
                va="center",
                fontsize=8.3,
                color="black",
                clip_on=False,
            )


def _stacked(ax, data: pd.DataFrame, keys: list[str], labs: list[str], xlabel: str, xmax: float) -> None:
    y_pos = list(range(len(data)))
    left = [0.0] * len(data)
    for k, (key, lab) in enumerate(zip(keys, labs)):
        values = [float(v) for v in data[key]]
        ax.barh(y_pos, values, left=left, color=COLS[k], label=lab, height=0.64)
        _label_segments(ax, left, values, y_pos, k)
        left = [a + b for a, b in zip(left, values)]
    ax.set_yticks(y_pos)
    ax.set_yticklabels(data.cohort.str.upper())
    ax.set_xlabel(xlabel)
    ax.set_xlim(0, xmax)
    ax.set_ylim(-0.55, len(data) - 0.45)
    ncol = 2 if len(labs) == 2 else 3
    ax.legend(
        frameon=False,
        ncol=ncol,
        loc="upper center",
        bbox_to_anchor=(0.44, 1.22),
        fontsize=8.3,
        handlelength=1.1,
        columnspacing=0.8,
    )


def _cell(ax, x, y, w, h, held: bool) -> None:
    ax.add_patch(
        Rectangle(
            (x, y),
            w,
            h,
            facecolor=HOLD if held else TRAIN,
            edgecolor="white",
            lw=0.55,
        )
    )


def _panel_c(ax) -> None:
    """Random vs patient-grouped vs LOGO, distinguished by grouping unit."""
    ax.set_xlim(0, 3)
    ax.set_ylim(0, 1)
    ax.axis("off")

    titles = (
        "Random sample folds",
        "Patient-grouped folds",
        "Leave-one-cohort-out",
    )
    units = (
        "Grouping unit: sample",
        "Grouping unit: patient",
        "Grouping unit: cohort",
    )
    lines = (
        "Samples split independently\nPatient pairs may cross folds",
        "Each patient stays in one fold\nCohorts appear in train and test",
        "One complete cohort held out\nCancer and acquisition shift together",
    )
    # (cohort, patient, row) with row 0 = tumour (top), row 1 = normal.
    random_hold = {(0, 0, 0), (1, 1, 1), (2, 0, 1)}
    patient_hold = {(0, 1), (2, 0)}

    cell_w, cell_h = 0.095, 0.105
    y_top = 0.68
    for j in range(3):
        x0 = j + 0.06
        ax.text(x0 + 0.44, 1.01, titles[j], ha="center", fontsize=9)
        ax.text(x0 + 0.44, 0.90, units[j], ha="center", fontsize=8.3, color="#333333")
        for c in range(3):
            for p in range(2):
                x = x0 + 0.08 + c * 0.28 + p * 0.112
                for r in range(2):
                    y = y_top - r * (cell_h + 0.018)
                    if j == 0:
                        held = (c, p, r) in random_hold
                    elif j == 1:
                        held = (c, p) in patient_hold
                    else:
                        held = c == 2
                    _cell(ax, x, y, cell_w, cell_h, held)
                if j == 1:
                    ax.add_patch(
                        Rectangle(
                            (x - 0.008, y_top - cell_h - 0.026),
                            cell_w + 0.016,
                            2 * cell_h + 0.044,
                            fill=False,
                            edgecolor="#444444",
                            lw=0.9,
                        )
                    )
            if j == 2:
                x_left = x0 + 0.08 + c * 0.28 - 0.012
                ax.add_patch(
                    Rectangle(
                        (x_left, y_top - cell_h - 0.032),
                        2 * cell_w + 0.112 - 0.095 + 0.024,
                        2 * cell_h + 0.056,
                        fill=False,
                        edgecolor="#444444",
                        lw=0.9,
                    )
                )
        ax.text(x0 + 0.44, 0.36, lines[j], ha="center", va="top", fontsize=8.4, linespacing=1.45)

    ax.text(
        1.5,
        0.04,
        "Blue: training    Orange: held out    Boxes: grouping unit (schematic, not a real fold)",
        ha="center",
        fontsize=8.3,
    )


def _physical_sampling(fig):
    ax=fig.add_axes([.04,.765,.935,.18]);ax.set_xlim(0,10);ax.set_ylim(0,2);ax.axis('off')
    # Small vector marks identify the three measured objects in the workflow:
    # CPTAC sample, mass-spectrometry measurement, and cohort-level transfer.
    for x, label, color in [(1.48,'CPTAC','#0072B2'),(4.30,'MS','#D55E00'),(9.05,'LOGO','#009E73')]:
        ax.text(x,1.94,label,ha='center',va='center',fontsize=8.3,color=color,
                bbox=dict(boxstyle='round,pad=.28',fc='white',ec=color,lw=.8))
    ax.text(1.4,2.24,'Tissue samples',ha='center',fontsize=9)
    for y,color,label in [(1.25,'#0072B2','Tumour'),(.40,'#E69F00','Source-defined\nnormal')]:
        ax.add_patch(Rectangle((.05,y-.24),.75,.48,facecolor=color,alpha=.12,edgecolor=color,lw=.8))
        for dx,dy in [(.15,-.10),(.35,.08),(.55,-.08)]:ax.add_patch(Circle((dx,y+dy),.08,facecolor=color,alpha=.75,edgecolor='white',lw=.6))
        ax.text(.95,y,label,fontsize=8.3,va='center')
    for x1,x2 in [(2.8,3.35),(5.25,5.75),(7.70,8.15)]:ax.add_patch(FancyArrowPatch((x1,.90),(x2,.90),arrowstyle='->',mutation_scale=11,lw=1,color='#666666'))
    ax.text(4.3,2.24,'Protein measurement',ha='center',fontsize=9)
    ax.text(4.3,1.42,'Extract proteins',ha='center',fontsize=8.3)
    ax.text(4.3,1.12,'Digest to peptides',ha='center',fontsize=8.3)
    ax.text(4.3,0.82,'Mass spectrometry',ha='center',fontsize=8.3)
    # Spectrum sits below the labels so peaks do not cross the words.
    ax.plot([3.70,4.90],[.18,.18],color='#D55E00',lw=1.0)
    ax.plot([3.95,3.95],[.18,.42],color='#D55E00',lw=1.0)
    ax.plot([4.30,4.30],[.18,.58],color='#D55E00',lw=1.0)
    ax.plot([4.65,4.65],[.18,.34],color='#D55E00',lw=1.0)
    ax.text(6.77,2.24,'Quantified abundance',ha='center',fontsize=9)
    for i in range(3):
        for j in range(5):ax.add_patch(Rectangle((5.88+j*.32,.38+i*.33),.30,.30,facecolor='#E6EFF5',edgecolor='#ACBDC8',lw=.5))
    ax.text(6.65,1.58,'Proteins',ha='center',fontsize=8.3)
    ax.text(6.72,.03,'1,191 samples × 8,906 accessions',ha='center',fontsize=8.3)
    ax.text(9.05,2.24,'Transfer test',ha='center',fontsize=9)
    ax.text(9.05,1.2,'New cohort',ha='center',fontsize=8.3)
    ax.text(9.05,.78,'Tissue biology +',ha='center',fontsize=8.3)
    ax.text(9.05,.42,'acquisition conditions',ha='center',fontsize=8.3)
    fig.text(.51,.731,'Sampling schematic; no tissue image, spectrum or molecular mechanism is inferred',ha='center',fontsize=8.3,color='#444444')
    fig.text(.51,.697,'Fixed 8,906-feature intersection uses all six cohorts before LOGO; target feature availability is known',
             ha='center',fontsize=8.3,color='#374151')


def main() -> None:
    data = _read_counts()
    expected = {
        "Ccrcc": (110, 84, 84, 26, 0),
        "Luad": (111, 102, 102, 9, 0),
        "Ucec": (103, 49, 30, 73, 19),
        "Hnscc": (110, 72, 68, 42, 2),
        "Lscc": (110, 102, 102, 8, 0),
        "Pdac": (145, 93, 84, 61, 9),
    }
    for row in data.itertuples(index=False):
        got = (
            int(row.n_tumor),
            int(row.n_normal),
            int(row.n_paired_patients),
            int(row.n_tumor_only_patients),
            int(row.n_normal_only_patients),
        )
        if got != expected[row.cohort]:
            raise SystemExit(f"pair-count mismatch for {row.cohort}: {got}")

    # Fail closed for both used faces instead of Matplotlib's silent DejaVu fallback.
    for weight in ("normal", "bold"):
        path = font_manager.findfont(
            font_manager.FontProperties(family="Arial", weight=weight),
            fallback_to_default=False,
        )
        font = font_manager.get_font(path)
        if font.family_name != "Arial" or ("bold" in font.style_name.lower()) != (weight == "bold"):
            raise SystemExit(f"B28 requires exact Arial {weight}; no substitution")
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial"],
            "font.weight": "normal",
            "axes.titleweight": "normal",
            "axes.labelweight": "normal",
            "font.size": 9,
            "pdf.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    fig = plt.figure(figsize=(6.5, 6.05))
    _physical_sampling(fig)
    grid = fig.add_gridspec(2, 2, height_ratios=[1.18, 0.88], wspace=0.4, hspace=0.55)

    ax_a = fig.add_subplot(grid[0, 0])
    _stacked(ax_a, data, ["n_tumor", "n_normal"], ["Tumour", "Normal"], "Samples", 260)

    ax_b = fig.add_subplot(grid[0, 1])
    _stacked(
        ax_b,
        data,
        ["n_paired_patients", "n_tumor_only_patients", "n_normal_only_patients"],
        ["Both classes", "Tumour only", "Normal only"],
        "Patients",
        178,
    )

    _panel_c(fig.add_subplot(grid[1, :]))
    fig.subplots_adjust(left=0.10, right=0.985, top=0.60, bottom=0.045)
    # Figure-relative letters share a row baseline, independent of both legends.
    for letter, x, y in (("A", .04, .979), ("B", .04, .680),
                         ("C", .55, .680), ("D", .04, .278)):
        fig.text(x, y, letter, ha="left", va="top", weight="bold", fontsize=11)
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"b28_F0_design.{ext}", dpi=350)
    plt.close(fig)


if __name__ == "__main__":
    main()
