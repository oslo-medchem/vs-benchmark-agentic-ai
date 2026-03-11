#!/usr/bin/env python3
"""Generate graphical abstract for Bioinformatics Application Note.

Layout (3 columns, left to right):
  1. AI Agent + Skill File
  2. Autonomous Pipeline (4 stages)
  3. Results: ROC curves + ΔAUC callout

Outputs: graphical_abstract.{png,eps,pdf,svg}
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.gridspec import GridSpec
import numpy as np
import json
from pathlib import Path

# ── Paths ────────────────────────────────────────────────────
EVAL_DIR = Path(__file__).resolve().parent.parent / "05_evaluation"
OUT_DIR = Path(__file__).resolve().parent / "figures"
OUT_DIR.mkdir(exist_ok=True)

# ── Colours ──────────────────────────────────────────────────
C_AGENT = "#3B82F6"       # blue
C_SKILL = "#8B5CF6"       # purple
C_PIPE = "#0EA5E9"        # sky blue
C_NAIVE = "#64748B"       # slate
C_SKILLR = "#EF4444"      # red
C_BG = "#F8FAFC"          # background
C_ACCENT = "#10B981"      # emerald
C_TEXT = "#1E293B"         # dark text
C_MUTED = "#94A3B8"       # muted text


def rbox(ax, x, y, w, h, color, text, fs=8, fw="bold", tc="white",
         ec="none", lw=0, rad=0.03, zorder=2):
    """Rounded box with centred text."""
    box = FancyBboxPatch((x, y), w, h,
                         boxstyle=f"round,pad=0,rounding_size={rad}",
                         facecolor=color, edgecolor=ec, linewidth=lw,
                         zorder=zorder)
    ax.add_patch(box)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=fs, fontweight=fw, color=tc, zorder=zorder + 1)


def draw_arrow(ax, x1, y1, x2, y2, color=C_MUTED, lw=1.5):
    """Simple arrow."""
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw),
                zorder=3)


def make_graphical_abstract():
    """Create the graphical abstract."""

    roc_naive = json.loads((EVAL_DIR / "results" / "roc_naive.json").read_text())
    roc_skill = json.loads((EVAL_DIR / "results" / "roc_skill.json").read_text())

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 9,
    })

    # Use GridSpec for clean 3-column layout
    fig = plt.figure(figsize=(11, 4.0), dpi=300, facecolor="white")
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1.0, 1.0, 1.3],
                  left=0.02, right=0.98, bottom=0.12, top=0.88,
                  wspace=0.08)

    # ── Column 1: Agent + Skill ──────────────────────────────
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_xlim(0, 4)
    ax1.set_ylim(0, 5)
    ax1.set_aspect("equal")
    ax1.axis("off")

    # Background card
    rbox(ax1, 0.1, 0.1, 3.8, 4.8, C_BG, "", ec="#E2E8F0", lw=1, rad=0.12, zorder=0)

    # Section label
    ax1.text(2.0, 4.55, "Input", ha="center", va="center",
             fontsize=8, fontweight="bold", color=C_MUTED, zorder=5)

    # Agent box
    rbox(ax1, 0.4, 2.8, 3.2, 1.2, C_AGENT, "Claude Code\n(LLM Agent)",
         fs=10, fw="bold", ec="#2563EB", lw=1.5, rad=0.08)

    # "+" symbol
    ax1.text(2.0, 2.35, "+", ha="center", va="center", fontsize=18,
             fontweight="bold", color=C_TEXT, zorder=5)

    # Skill file box
    rbox(ax1, 0.4, 1.0, 3.2, 1.2, C_SKILL, "Skill File\n(212 lines)",
         fs=10, fw="bold", ec="#7C3AED", lw=1.5, rad=0.08)

    # Label
    ax1.text(2.0, 0.45, "Domain Knowledge Injection",
             ha="center", va="center", fontsize=7, color=C_MUTED,
             fontstyle="italic", fontweight="bold", zorder=5)

    # ── Column 2: Pipeline stages ────────────────────────────
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.set_xlim(0, 4)
    ax2.set_ylim(0, 5)
    ax2.set_aspect("equal")
    ax2.axis("off")

    # Background card
    rbox(ax2, 0.1, 0.1, 3.8, 4.8, C_BG, "", ec="#E2E8F0", lw=1, rad=0.12, zorder=0)

    # Section label
    ax2.text(2.0, 4.55, "Autonomous Pipeline", ha="center", va="center",
             fontsize=8, fontweight="bold", color=C_PIPE, zorder=5)
    ax2.text(2.0, 4.20, "18 scripts  |  ~3,100 lines  |  zero human editing",
             ha="center", va="center", fontsize=5.5, color=C_MUTED, zorder=5)

    stages = [
        "ChEMBL\nData Retrieval",
        "Active Curation &\nDecoy Generation",
        "PAINS / Brenk\nFiltering",
        "GPU Docking\n(Uni-Dock)",
    ]

    box_w = 3.0
    box_h = 0.68
    box_x = 0.5
    y_top = 3.40
    y_gap = 0.20  # gap between boxes

    for i, label in enumerate(stages):
        y = y_top - i * (box_h + y_gap)
        rbox(ax2, box_x, y, box_w, box_h, C_PIPE, label,
             fs=8, fw="bold", ec="#0284C7", lw=1, rad=0.06)
        # Down arrow
        if i < len(stages) - 1:
            ay1 = y
            ay2 = y - y_gap
            draw_arrow(ax2, box_x + box_w / 2, ay1, box_x + box_w / 2, ay2,
                       color=C_MUTED, lw=1.5)

    # GPU note
    ax2.text(box_x + box_w + 0.05, y_top - 3 * (box_h + y_gap) + box_h / 2,
             "2× RTX\n4500 Ada", ha="left", va="center", fontsize=6,
             color=C_MUTED, fontweight="bold", zorder=5)

    # ── Column 3: Results ────────────────────────────────────
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.set_xlim(0, 5.2)
    ax3.set_ylim(0, 5)
    ax3.set_aspect("equal")
    ax3.axis("off")

    # Background card
    rbox(ax3, 0.1, 0.1, 5.0, 4.8, C_BG, "", ec="#E2E8F0", lw=1, rad=0.12, zorder=0)

    # Section label
    ax3.text(2.6, 4.55, "Results", ha="center", va="center",
             fontsize=8, fontweight="bold", color=C_ACCENT, zorder=5)

    # ROC curves as inset axes within ax3
    # Get ax3 position in figure coords
    pos3 = ax3.get_position()
    roc_left = pos3.x0 + 0.02
    roc_bottom = pos3.y0 + 0.08
    roc_w = pos3.width * 0.55
    roc_h = pos3.height * 0.72

    ax_roc = fig.add_axes([roc_left, roc_bottom, roc_w, roc_h])

    fpr_n = np.array(roc_naive["fpr"])
    tpr_n = np.array(roc_naive["tpr"])
    fpr_s = np.array(roc_skill["fpr"])
    tpr_s = np.array(roc_skill["tpr"])

    ax_roc.plot([0, 1], [0, 1], "--", color="#CBD5E1", lw=0.8, zorder=1)
    ax_roc.fill_between(fpr_n, tpr_n, alpha=0.08, color=C_NAIVE, zorder=2)
    ax_roc.plot(fpr_n, tpr_n, color=C_NAIVE, lw=2.0,
                label="Naive (AUC = 0.658)", zorder=3)
    ax_roc.fill_between(fpr_s, tpr_s, alpha=0.10, color=C_SKILLR, zorder=2)
    ax_roc.plot(fpr_s, tpr_s, color=C_SKILLR, lw=2.0,
                label="Skill (AUC = 0.710)", zorder=3)

    ax_roc.set_xlabel("False Positive Rate", fontsize=7, labelpad=2)
    ax_roc.set_ylabel("True Positive Rate", fontsize=7, labelpad=2)
    ax_roc.set_xlim(0, 1)
    ax_roc.set_ylim(0, 1)
    ax_roc.tick_params(labelsize=6, length=2, pad=1)
    ax_roc.legend(fontsize=6.5, loc="lower right", frameon=True,
                  fancybox=True, framealpha=0.95, edgecolor="#E2E8F0")
    for spine in ax_roc.spines.values():
        spine.set_color("#CBD5E1")
        spine.set_linewidth(0.8)

    # ΔAUC callout — right of ROC, inside the results card
    cx = 4.25
    rbox(ax3, 3.35, 0.8, 1.75, 3.0, "white", "",
         ec=C_ACCENT, lw=2, rad=0.1, zorder=4)

    ax3.text(cx, 3.45, "ROC AUC", ha="center", va="center",
             fontsize=8.5, fontweight="bold", color=C_TEXT, zorder=5)

    ax3.text(cx, 2.65, "+0.048", ha="center", va="center",
             fontsize=17, fontweight="bold", color=C_ACCENT, zorder=5)

    ax3.text(cx, 2.05, "0.658 → 0.710", ha="center", va="center",
             fontsize=7, fontweight="bold", color=C_MUTED, zorder=5)

    ax3.text(cx, 1.40, "DeLong\np = 0.004", ha="center", va="center",
             fontsize=7, fontweight="bold", color=C_TEXT, zorder=5)

    # ── Title banner ─────────────────────────────────────────
    fig.text(0.50, 0.95,
             "Agentic AI with Structured Skill Files for Virtual Screening",
             ha="center", va="center", fontsize=12, fontweight="bold",
             color=C_TEXT)

    # ── Bottom takeaway bar ──────────────────────────────────
    # Use figure-level rectangle
    from matplotlib.patches import Rectangle
    bar = Rectangle((0.02, 0.01), 0.96, 0.055, transform=fig.transFigure,
                    facecolor=C_TEXT, edgecolor="none", alpha=0.9,
                    zorder=10, clip_on=False)
    fig.patches.append(bar)
    fig.text(0.50, 0.035,
             "Lightweight skill files democratise rigorous large-scale virtual screening",
             ha="center", va="center", fontsize=8.5, fontweight="bold",
             color="white", zorder=11)

    # ── Connecting arrows between columns ────────────────────
    # Use figure-level annotation for cross-axes arrows
    # Arrow from col1 to col2
    fig.patches.append(
        FancyArrowPatch(
            (0.345, 0.50), (0.365, 0.50),
            transform=fig.transFigure,
            arrowstyle="-|>", color=C_AGENT, lw=2.5,
            mutation_scale=15, zorder=10, clip_on=False
        )
    )
    # Arrow from col2 to col3
    fig.patches.append(
        FancyArrowPatch(
            (0.665, 0.50), (0.685, 0.50),
            transform=fig.transFigure,
            arrowstyle="-|>", color=C_PIPE, lw=2.5,
            mutation_scale=15, zorder=10, clip_on=False
        )
    )

    return fig


if __name__ == "__main__":
    print("Generating graphical abstract...")
    fig = make_graphical_abstract()

    stem = "graphical_abstract"
    for fmt in ("png", "eps", "pdf", "svg"):
        out = OUT_DIR / f"{stem}.{fmt}"
        dpi = 600 if fmt == "png" else None
        fig.savefig(out, dpi=dpi, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)

    # Also copy to evaluation figures
    eval_fig_dir = EVAL_DIR / "figures"
    eval_fig_dir.mkdir(exist_ok=True)
    import shutil
    for fmt in ("png", "eps", "pdf", "svg"):
        shutil.copy2(OUT_DIR / f"{stem}.{fmt}", eval_fig_dir / f"{stem}.{fmt}")

    for f in sorted(OUT_DIR.glob(f"{stem}.*")):
        print(f"  {f.name:40s} ({f.stat().st_size:,} bytes)")
    print("Done.")
