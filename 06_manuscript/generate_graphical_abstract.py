#!/usr/bin/env python3
"""Generate compact, high-quality graphical abstract.

Single-canvas landscape, tight Nature-quality layout:
  Left: Agent + Skill | Centre: Pipeline | Right: ΔAUC banner + ROC plot

Outputs: graphical_abstract.{png,eps,pdf,svg} at 600 DPI
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
import numpy as np
import json
from pathlib import Path

EVAL_DIR = Path(__file__).resolve().parent.parent / "05_evaluation"
OUT_DIR = Path(__file__).resolve().parent / "figures"
OUT_DIR.mkdir(exist_ok=True)

# Palette
C_BLUE   = "#2563EB"
C_PURPLE = "#7C3AED"
C_SKY    = "#0EA5E9"
C_SKY_D  = "#0369A1"
C_GREEN  = "#059669"
C_RED    = "#DC2626"
C_SLATE  = "#64748B"
C_TEXT   = "#0F172A"
C_MUTED  = "#94A3B8"
C_BORDER = "#CBD5E1"


def pill(ax, x, y, w, h, fc, text="", fs=7.5, fw="bold", tc="white",
         ec="none", lw=0, rad=0.012):
    p = FancyBboxPatch((x, y), w, h,
                       boxstyle=f"round,pad=0,rounding_size={rad}",
                       fc=fc, ec=ec, lw=lw, zorder=2)
    ax.add_patch(p)
    if text:
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
                fontsize=fs, fontweight=fw, color=tc, zorder=3,
                linespacing=1.1)


def make_graphical_abstract():
    roc_n = json.loads((EVAL_DIR / "results" / "roc_naive.json").read_text())
    roc_s = json.loads((EVAL_DIR / "results" / "roc_skill.json").read_text())

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica Neue", "DejaVu Sans"],
        "font.size": 8,
    })

    fig, ax = plt.subplots(figsize=(9.5, 4.2), dpi=300, facecolor="white")
    ax.set_xlim(0, 9.5)
    ax.set_ylim(0, 4.2)
    ax.set_aspect("equal")
    ax.axis("off")

    # ── Title ────────────────────────────────────────────────
    ax.text(4.75, 4.02,
            "Agentic AI with Structured Skill Files for Virtual Screening",
            ha="center", va="center", fontsize=12, fontweight="bold",
            color=C_TEXT)
    ax.text(4.75, 3.78,
            "FPR2 · PDB 7T6S · Claude Code Opus 4.6 · "
            "18 scripts · ~3,100 LOC · zero human editing",
            ha="center", va="center", fontsize=6, color=C_MUTED)
    ax.plot([0.2, 9.3], [3.65, 3.65], color=C_BORDER, lw=0.7, zorder=1)

    # ════════════════════════════════════════════════════════════
    # LEFT: Agent + Skill (x: 0.15–2.30)
    # ════════════════════════════════════════════════════════════
    lx, lw_ = 0.15, 2.00

    pill(ax, lx, 2.15, lw_, 0.80, C_BLUE,
         "Claude Code\n(LLM Agent)", fs=9.5, ec="#1D4ED8", lw=1.2, rad=0.06)

    ax.text(lx + lw_ / 2, 1.85, "+", ha="center", va="center",
            fontsize=15, fontweight="bold", color=C_TEXT, zorder=5)

    pill(ax, lx, 0.85, lw_, 0.80, C_PURPLE,
         "Skill File\n(212 lines)", fs=9.5, ec="#6D28D9", lw=1.2, rad=0.06)

    ax.text(lx + lw_ / 2, 0.58, "Domain knowledge injection",
            ha="center", va="center", fontsize=6.5, color=C_MUTED,
            fontstyle="italic")

    # Arrow left → centre
    ax.annotate("", xy=(2.55, 1.85), xytext=(2.28, 1.85),
                arrowprops=dict(arrowstyle="-|>", color=C_BLUE,
                                lw=2.2, mutation_scale=14), zorder=5)

    # ════════════════════════════════════════════════════════════
    # CENTRE: Pipeline (x: 2.60–5.55)
    # ════════════════════════════════════════════════════════════
    stages = [
        ("ChEMBL Query",  "100 actives"),
        ("Decoy Gen.",    "4,780 matched"),
        ("PAINS/Brenk",   "762 removed"),
        ("GPU Docking",   "2× RTX 4500"),
    ]

    px, pw_, ph = 2.65, 2.70, 0.52
    gap = 0.13
    y0 = 3.00

    ax.text(px + pw_ / 2, 3.45, "Autonomous Pipeline",
            ha="center", va="center", fontsize=7.5, fontweight="bold",
            color=C_SKY_D)

    for i, (label, note) in enumerate(stages):
        y = y0 - i * (ph + gap)
        pill(ax, px, y, pw_, ph, C_SKY, "", ec=C_SKY_D, lw=0.8, rad=0.04)
        ax.text(px + 0.15, y + ph / 2, label, ha="left", va="center",
                fontsize=8.5, fontweight="bold", color="white", zorder=3)
        ax.text(px + pw_ - 0.15, y + ph / 2, note, ha="right", va="center",
                fontsize=6.5, color="#E0F2FE", zorder=3)
        if i < len(stages) - 1:
            ax.annotate("", xy=(px + pw_ / 2, y - 0.01),
                        xytext=(px + pw_ / 2, y + 0.01),
                        arrowprops=dict(arrowstyle="-|>", color=C_BORDER,
                                        lw=1.2, mutation_scale=10), zorder=2)

    # Arrow centre → right
    ax.annotate("", xy=(5.75, 1.85), xytext=(5.48, 1.85),
                arrowprops=dict(arrowstyle="-|>", color=C_SKY,
                                lw=2.2, mutation_scale=14), zorder=5)

    # ════════════════════════════════════════════════════════════
    # RIGHT: ΔAUC callout (top) + ROC (bottom)  x: 5.80–9.30
    # ════════════════════════════════════════════════════════════

    # ── ΔAUC callout banner ──────────────────────────────────
    bx, bw, bh = 5.90, 3.25, 1.10
    by = 2.45
    pill(ax, bx, by, bw, bh, "#F0FDF4", "",
         ec=C_GREEN, lw=2, rad=0.08)

    bcx = bx + bw / 2
    ax.text(bcx - 0.70, by + bh / 2, "+0.048",
            ha="center", va="center", fontsize=22, fontweight="bold",
            color=C_GREEN, zorder=5)

    ax.text(bcx + 0.65, by + bh * 0.70, "ROC AUC",
            ha="center", va="center", fontsize=8, fontweight="bold",
            color=C_TEXT, zorder=5)
    ax.text(bcx + 0.65, by + bh * 0.45, "0.658 → 0.710",
            ha="center", va="center", fontsize=7.5, fontweight="bold",
            color=C_SLATE, zorder=5)
    ax.text(bcx + 0.65, by + bh * 0.20, "DeLong p = 0.004",
            ha="center", va="center", fontsize=6.5, color=C_SLATE,
            zorder=5)

    # ── ROC plot ─────────────────────────────────────────────
    pos = ax.get_position()
    roc_l = pos.x0 + pos.width * 0.63
    roc_b = pos.y0 + pos.height * 0.16
    roc_w = pos.width * 0.32
    roc_h = pos.height * 0.40

    ax_roc = fig.add_axes([roc_l, roc_b, roc_w, roc_h])

    fpr_n, tpr_n = np.array(roc_n["fpr"]), np.array(roc_n["tpr"])
    fpr_s, tpr_s = np.array(roc_s["fpr"]), np.array(roc_s["tpr"])

    ax_roc.plot([0, 1], [0, 1], "--", color="#E2E8F0", lw=0.5, zorder=1)
    ax_roc.fill_between(fpr_n, tpr_n, alpha=0.05, color=C_SLATE)
    ax_roc.plot(fpr_n, tpr_n, color=C_SLATE, lw=1.5,
                label="Naive (0.658)", zorder=3)
    ax_roc.fill_between(fpr_s, tpr_s, alpha=0.07, color=C_RED)
    ax_roc.plot(fpr_s, tpr_s, color=C_RED, lw=1.5,
                label="Skill (0.710)", zorder=3)

    ax_roc.set_xlim(0, 1); ax_roc.set_ylim(0, 1)
    ax_roc.set_xlabel("FPR", fontsize=6, labelpad=1)
    ax_roc.set_ylabel("TPR", fontsize=6, labelpad=1)
    ax_roc.tick_params(labelsize=5, length=1.5, pad=1)
    ax_roc.legend(fontsize=5.5, loc="lower right", frameon=True,
                  framealpha=0.95, edgecolor=C_BORDER, handlelength=1,
                  borderpad=0.3, handletextpad=0.3)
    for sp in ax_roc.spines.values():
        sp.set_color(C_BORDER); sp.set_linewidth(0.5)

    # ── Footer ───────────────────────────────────────────────
    bar = Rectangle((0.02, 0.005), 0.96, 0.06, transform=fig.transFigure,
                    fc=C_TEXT, ec="none", alpha=0.92, zorder=10,
                    clip_on=False)
    fig.patches.append(bar)
    fig.text(0.50, 0.035,
             "Lightweight skill files democratise rigorous "
             "large-scale virtual screening",
             ha="center", va="center", fontsize=8.5, fontweight="bold",
             color="white", zorder=11)

    return fig


if __name__ == "__main__":
    print("Generating graphical abstract...")
    fig = make_graphical_abstract()

    stem = "graphical_abstract"
    for fmt in ("png", "eps", "pdf", "svg"):
        out = OUT_DIR / f"{stem}.{fmt}"
        dpi = 600 if fmt == "png" else None
        fig.savefig(out, dpi=dpi, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)

    eval_fig_dir = EVAL_DIR / "figures"
    eval_fig_dir.mkdir(exist_ok=True)
    import shutil
    for fmt in ("png", "eps", "pdf", "svg"):
        shutil.copy2(OUT_DIR / f"{stem}.{fmt}", eval_fig_dir / f"{stem}.{fmt}")

    for f in sorted(OUT_DIR.glob(f"{stem}.*")):
        print(f"  {f.name:40s} ({f.stat().st_size:,} bytes)")
    print("Done.")
