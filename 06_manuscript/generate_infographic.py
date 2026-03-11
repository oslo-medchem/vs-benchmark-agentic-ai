#!/usr/bin/env python3
"""Generate a compact, high-resolution infographic for journal submission.

Portrait layout with 4 sections:
  1. Concept: AI Agent + Skill File
  2. Autonomous Pipeline (6 stages)
  3. Protocol Comparison (Naive vs Skill table)
  4. Results (ROC curves + key metrics)

Outputs: infographic.{png,eps,pdf,svg} at 600 DPI
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import json
from pathlib import Path

# ── Paths ────────────────────────────────────────────────────
EVAL_DIR = Path(__file__).resolve().parent.parent / "05_evaluation"
OUT_DIR = Path(__file__).resolve().parent / "figures"
OUT_DIR.mkdir(exist_ok=True)

# ── Colours ──────────────────────────────────────────────────
BG_DARK   = "#0F172A"
C_BLUE    = "#2563EB"
C_BLUE_L  = "#3B82F6"
C_PURPLE  = "#7C3AED"
C_PURPLE_L= "#8B5CF6"
C_SKY     = "#0EA5E9"
C_SKY_D   = "#0284C7"
C_GREEN   = "#10B981"
C_GREEN_D = "#059669"
C_RED     = "#EF4444"
C_SLATE   = "#64748B"
C_SLATE_D = "#475569"
C_TEXT    = "#1E293B"
C_MUTED   = "#94A3B8"
C_BORDER  = "#E2E8F0"
C_BG_ALT  = "#F1F5F9"
C_BG_CARD = "#F8FAFC"


def rbox(ax, x, y, w, h, fc, text="", fs=7, fw="normal", tc="white",
         ec="none", lw=0, rad=0.008, zorder=2, linespacing=1.2):
    """Rounded rectangle with centred text."""
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0,rounding_size={rad}",
        facecolor=fc, edgecolor=ec, linewidth=lw, zorder=zorder,
    )
    ax.add_patch(patch)
    if text:
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
                fontsize=fs, fontweight=fw, color=tc, zorder=zorder + 1,
                linespacing=linespacing)


def section_bar(ax, y, label, color, number):
    """Draw a numbered section header bar."""
    rbox(ax, 0.03, y, 0.94, 0.022, color, rad=0.005, zorder=3)
    ax.text(0.06, y + 0.011, f"{number:02d}", ha="center", va="center",
            fontsize=8, fontweight="bold", color="white",
            bbox=dict(boxstyle="round,pad=0.15", fc=color, ec="white", lw=0.8),
            zorder=4)
    ax.text(0.10, y + 0.011, label, ha="left", va="center",
            fontsize=8.5, fontweight="bold", color="white", zorder=4)


def make_infographic():
    """Build the infographic."""

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "figure.facecolor": "white",
    })

    roc_naive = json.loads((EVAL_DIR / "results" / "roc_naive.json").read_text())
    roc_skill = json.loads((EVAL_DIR / "results" / "roc_skill.json").read_text())

    fig, ax = plt.subplots(figsize=(7.5, 13.0), dpi=300)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("auto")
    ax.axis("off")

    # Outer border
    rbox(ax, 0.005, 0.005, 0.99, 0.99, "white", ec=C_BORDER, lw=1.5,
         rad=0.008, zorder=0)

    # ════════════════════════════════════════════════════════════
    # HEADER — y: 0.94–1.00
    # ════════════════════════════════════════════════════════════
    rbox(ax, 0.01, 0.94, 0.98, 0.05, BG_DARK, rad=0.008, zorder=2)
    ax.text(0.50, 0.974, "Agentic AI with Structured Skill Files",
            ha="center", va="center", fontsize=13, fontweight="bold",
            color="white", zorder=5)
    ax.text(0.50, 0.955, "Autonomously Generates GPU-Accelerated Virtual Screening",
            ha="center", va="center", fontsize=9, color=C_MUTED, zorder=5)
    ax.text(0.50, 0.944, "FPR2 (CHEMBL4227)  ·  PDB 7T6S  ·  Claude Code Opus 4.6",
            ha="center", va="center", fontsize=6.5, color="#64748B", zorder=5)

    # ════════════════════════════════════════════════════════════
    # SECTION 1: CONCEPT — y: 0.86–0.93
    # ════════════════════════════════════════════════════════════
    section_bar(ax, 0.915, "THE CONCEPT", C_BLUE, 1)

    # Agent + Skill → Pipeline
    rbox(ax, 0.06, 0.87, 0.24, 0.038, C_BLUE_L,
         "Claude Code\n(LLM Agent)", fs=8, fw="bold",
         ec=C_BLUE, lw=1.2, rad=0.006)

    ax.text(0.33, 0.889, "+", ha="center", va="center",
            fontsize=14, fontweight="bold", color=C_TEXT, zorder=5)

    rbox(ax, 0.36, 0.87, 0.24, 0.038, C_PURPLE_L,
         "Skill File\n(212 lines)", fs=8, fw="bold",
         ec=C_PURPLE, lw=1.2, rad=0.006)

    ax.annotate("", xy=(0.72, 0.889), xytext=(0.64, 0.889),
                arrowprops=dict(arrowstyle="-|>", color=C_GREEN_D, lw=2), zorder=5)

    rbox(ax, 0.73, 0.87, 0.22, 0.038, C_GREEN,
         "Rigorous VS\nPipeline", fs=7.5, fw="bold",
         ec=C_GREEN_D, lw=1.2, rad=0.006)

    ax.text(0.50, 0.858, "Domain knowledge injection via plain-text — "
            "no fine-tuning, no RAG",
            ha="center", va="center", fontsize=6, color=C_SLATE,
            fontstyle="italic", zorder=5)

    # ════════════════════════════════════════════════════════════
    # SECTION 2: PIPELINE — y: 0.66–0.85
    # ════════════════════════════════════════════════════════════
    section_bar(ax, 0.835, "AUTONOMOUS PIPELINE", C_SKY, 2)

    stages = [
        ("ChEMBL API Query",      "382 FPR2 actives (pChEMBL ≥ 5)"),
        ("Butina Clustering",     "100 diverse actives (TC 0.4, ECFP4)"),
        ("Decoy Generation",      "4,780 property-matched decoys (MW/LogP/HBD/HBA/RotB)"),
        ("Structural Filtering",  "PAINS A/B/C + Brenk filters (762 compounds removed)"),
        ("Ligand Preparation",    "RDKit ETKDGv3 + Meeko 0.7.1 → PDBQT"),
        ("GPU Docking",           "Uni-Dock 1.1.3  ·  2× RTX 4500 Ada  ·  ~5,000 ligands"),
    ]

    row_h = 0.022
    y_start = 0.808
    for i, (label, detail) in enumerate(stages):
        y = y_start - i * (row_h + 0.004)
        bg = C_BG_ALT if i % 2 == 0 else "white"
        rbox(ax, 0.06, y, 0.88, row_h, bg, ec=C_BORDER, lw=0.4,
             rad=0.004, zorder=1)
        # Number circle
        circle = plt.Circle((0.09, y + row_h / 2), 0.009,
                             facecolor=C_SKY, edgecolor=C_SKY_D,
                             linewidth=0.6, zorder=3)
        ax.add_patch(circle)
        ax.text(0.09, y + row_h / 2, str(i + 1), ha="center", va="center",
                fontsize=6.5, fontweight="bold", color="white", zorder=4)
        # Label
        ax.text(0.115, y + row_h / 2, label, ha="left", va="center",
                fontsize=7, fontweight="bold", color=C_TEXT, zorder=3)
        # Detail
        ax.text(0.35, y + row_h / 2, detail, ha="left", va="center",
                fontsize=6.5, color=C_SLATE, zorder=3)

    # Stats bar
    sb_y = y_start - 6 * (row_h + 0.004) - 0.002
    rbox(ax, 0.06, sb_y, 0.88, 0.018, C_BG_ALT, ec=C_SKY, lw=0.8,
         rad=0.004, zorder=1)
    for text, sx in [("18 scripts", 0.21), ("~3,100 lines", 0.46),
                     ("zero human editing", 0.73)]:
        ax.text(sx, sb_y + 0.009, text, ha="center", va="center",
                fontsize=6.5, fontweight="bold", color=C_SKY_D, zorder=3)
    for sx in [0.33, 0.60]:
        ax.plot([sx, sx], [sb_y + 0.003, sb_y + 0.015],
                color=C_MUTED, lw=0.5, zorder=2)

    # ════════════════════════════════════════════════════════════
    # SECTION 3: PROTOCOL COMPARISON — y: 0.48–0.65
    # ════════════════════════════════════════════════════════════
    s3_top = sb_y - 0.018
    section_bar(ax, s3_top, "PROTOCOL COMPARISON", C_PURPLE, 3)

    rows = [
        ("Ligand 3D",       "OpenBabel gen3d",    "RDKit ETKDGv3"),
        ("PDBQT Prep",      "OpenBabel",          "Meeko 0.7.1"),
        ("Receptor Prep",   "OpenBabel",          "Meeko mk_prepare_receptor"),
        ("Struct. Filters",  "None",               "PAINS A/B/C + Brenk"),
        ("Box Size",        "20 Å",               "25 Å"),
        ("Exhaustiveness",  "8",                  "32"),
    ]

    th_y = s3_top - 0.004 - 0.020
    rbox(ax, 0.06, th_y, 0.88, 0.020, C_PURPLE, rad=0.004, zorder=2)
    for text, x in [("Parameter", 0.17), ("Naive", 0.46), ("Skill-Guided", 0.77)]:
        ax.text(x, th_y + 0.010, text, ha="center", va="center",
                fontsize=6.5, fontweight="bold", color="white", zorder=3)

    for i, (param, naive, skill) in enumerate(rows):
        y = th_y - (i + 1) * 0.020
        bg = C_BG_ALT if i % 2 == 0 else "white"
        rbox(ax, 0.06, y, 0.88, 0.020, bg, ec=C_BORDER, lw=0.3,
             rad=0.003, zorder=1)
        ax.text(0.17, y + 0.010, param, ha="center", va="center",
                fontsize=6.5, fontweight="bold", color=C_TEXT, zorder=3)
        ax.text(0.46, y + 0.010, naive, ha="center", va="center",
                fontsize=6.5, color=C_SLATE, zorder=3)
        ax.text(0.77, y + 0.010, skill, ha="center", va="center",
                fontsize=6.5, fontweight="bold", color=C_GREEN_D, zorder=3)

    # ════════════════════════════════════════════════════════════
    # SECTION 4: RESULTS — y: 0.04–0.46
    # ════════════════════════════════════════════════════════════
    s4_top = th_y - (len(rows) + 1.5) * 0.020
    section_bar(ax, s4_top, "RESULTS", C_GREEN, 4)

    # ROC plot — left half
    roc_bottom_abs = 0.12
    roc_top_abs = s4_top - 0.01
    roc_height = roc_top_abs - roc_bottom_abs
    ax_roc = fig.add_axes([0.09, roc_bottom_abs, 0.42, roc_height])

    fpr_n = np.array(roc_naive["fpr"])
    tpr_n = np.array(roc_naive["tpr"])
    fpr_s = np.array(roc_skill["fpr"])
    tpr_s = np.array(roc_skill["tpr"])

    ax_roc.plot([0, 1], [0, 1], "--", color="#CBD5E1", lw=0.6, zorder=1)
    ax_roc.fill_between(fpr_n, tpr_n, alpha=0.06, color=C_SLATE, zorder=2)
    ax_roc.plot(fpr_n, tpr_n, color=C_SLATE, lw=1.8,
                label="Naive (AUC = 0.658)", zorder=3)
    ax_roc.fill_between(fpr_s, tpr_s, alpha=0.08, color=C_RED, zorder=2)
    ax_roc.plot(fpr_s, tpr_s, color=C_RED, lw=1.8,
                label="Skill (AUC = 0.710)", zorder=3)

    ax_roc.set_xlabel("False Positive Rate", fontsize=7, labelpad=2)
    ax_roc.set_ylabel("True Positive Rate", fontsize=7, labelpad=2)
    ax_roc.set_xlim(0, 1)
    ax_roc.set_ylim(0, 1)
    ax_roc.tick_params(labelsize=5.5, length=2, pad=1)
    ax_roc.legend(fontsize=6, loc="lower right", frameon=True,
                  fancybox=True, framealpha=0.95, edgecolor=C_BORDER)
    ax_roc.set_title("ROC Curves — FPR2 Virtual Screening",
                     fontsize=7.5, fontweight="bold", pad=3, color=C_TEXT)
    for spine in ax_roc.spines.values():
        spine.set_color(C_BORDER)
        spine.set_linewidth(0.6)

    # Right half: metric cards — placed in main ax coordinates
    cx = 0.76  # centre x of cards
    cw = 0.34  # card width
    c_left = cx - cw / 2

    # Card 1: ΔAUC
    c1_y = s4_top - 0.035
    c1_h = 0.065
    rbox(ax, c_left, c1_y, cw, c1_h, C_BG_CARD,
         ec=C_GREEN, lw=2, rad=0.006, zorder=2)
    ax.text(cx, c1_y + c1_h * 0.80, "ROC AUC Improvement",
            ha="center", va="center", fontsize=7, fontweight="bold",
            color=C_TEXT, zorder=3)
    ax.text(cx, c1_y + c1_h * 0.42, "+0.048",
            ha="center", va="center", fontsize=18, fontweight="bold",
            color=C_GREEN, zorder=3)
    ax.text(cx, c1_y + c1_h * 0.12, "DeLong p = 0.004",
            ha="center", va="center", fontsize=6.5, fontweight="bold",
            color=C_SLATE_D, zorder=3)

    # Card 2: AUC values (two mini cards side by side)
    c2_y = c1_y - 0.050
    c2_h = 0.038
    half = (cw - 0.01) / 2

    rbox(ax, c_left, c2_y, half, c2_h, C_BG_CARD,
         ec=C_SLATE, lw=1, rad=0.005, zorder=2)
    ax.text(c_left + half / 2, c2_y + c2_h * 0.72, "Naive",
            ha="center", va="center", fontsize=6, color=C_SLATE, zorder=3)
    ax.text(c_left + half / 2, c2_y + c2_h * 0.30, "0.658",
            ha="center", va="center", fontsize=12, fontweight="bold",
            color=C_SLATE_D, zorder=3)

    rbox(ax, c_left + half + 0.01, c2_y, half, c2_h, C_BG_CARD,
         ec=C_RED, lw=1, rad=0.005, zorder=2)
    ax.text(c_left + half + 0.01 + half / 2, c2_y + c2_h * 0.72, "Skill",
            ha="center", va="center", fontsize=6, color=C_RED, zorder=3)
    ax.text(c_left + half + 0.01 + half / 2, c2_y + c2_h * 0.30, "0.710",
            ha="center", va="center", fontsize=12, fontweight="bold",
            color=C_RED, zorder=3)

    # Card 3: Additional metrics
    c3_y = c2_y - 0.042
    c3_h = 0.032
    rbox(ax, c_left, c3_y, cw, c3_h, C_BG_CARD,
         ec=C_BORDER, lw=0.8, rad=0.005, zorder=2)
    ax.text(cx, c3_y + c3_h * 0.65, "EF₁₀%: 2.2 (skill) vs 2.1 (naive)",
            ha="center", va="center", fontsize=6, color=C_TEXT, zorder=3)
    ax.text(cx, c3_y + c3_h * 0.28, "Mann–Whitney U: p < 10⁻¹⁰ (both)",
            ha="center", va="center", fontsize=6, color=C_TEXT, zorder=3)

    # Card 4: Library composition
    c4_y = c3_y - 0.042
    c4_h = 0.032
    rbox(ax, c_left, c4_y, cw, c4_h, C_BG_CARD,
         ec=C_BORDER, lw=0.8, rad=0.005, zorder=2)
    ax.text(cx, c4_y + c4_h * 0.65,
            "Naive: 4,632 docked (89 actives, 4,543 decoys)",
            ha="center", va="center", fontsize=6, color=C_TEXT, zorder=3)
    ax.text(cx, c4_y + c4_h * 0.28,
            "Skill: 4,082 docked (93 actives, 3,989 decoys)",
            ha="center", va="center", fontsize=6, color=C_TEXT, zorder=3)

    # Card 5: Confidence intervals
    c5_y = c4_y - 0.042
    c5_h = 0.032
    rbox(ax, c_left, c5_y, cw, c5_h, C_BG_CARD,
         ec=C_BORDER, lw=0.8, rad=0.005, zorder=2)
    ax.text(cx, c5_y + c5_h * 0.65,
            "95% CI: [0.591, 0.713] naive",
            ha="center", va="center", fontsize=6, color=C_TEXT, zorder=3)
    ax.text(cx, c5_y + c5_h * 0.28,
            "95% CI: [0.650, 0.748] skill",
            ha="center", va="center", fontsize=6, color=C_TEXT, zorder=3)

    # ════════════════════════════════════════════════════════════
    # FOOTER
    # ════════════════════════════════════════════════════════════
    from matplotlib.patches import Rectangle as Rect
    footer = Rect((0.01, 0.005), 0.98, 0.04,
                  transform=fig.transFigure, facecolor=BG_DARK,
                  edgecolor="none", zorder=10, clip_on=False)
    fig.patches.append(footer)
    fig.text(0.50, 0.025,
             "Lightweight skill files democratise rigorous large-scale "
             "virtual screening",
             ha="center", va="center", fontsize=8.5, fontweight="bold",
             color="white", zorder=11)

    return fig


if __name__ == "__main__":
    print("Generating infographic...")
    fig = make_infographic()

    stem = "infographic"
    for fmt in ("png", "eps", "pdf", "svg"):
        out = OUT_DIR / f"{stem}.{fmt}"
        dpi = 600 if fmt == "png" else None
        fig.savefig(out, dpi=dpi, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)

    # Copy to evaluation figures
    eval_fig_dir = EVAL_DIR / "figures"
    eval_fig_dir.mkdir(exist_ok=True)
    import shutil
    for fmt in ("png", "eps", "pdf", "svg"):
        shutil.copy2(OUT_DIR / f"{stem}.{fmt}", eval_fig_dir / f"{stem}.{fmt}")

    for f in sorted(OUT_DIR.glob(f"{stem}.*")):
        print(f"  {f.name:40s} ({f.stat().st_size:,} bytes)")
    print("Done.")
