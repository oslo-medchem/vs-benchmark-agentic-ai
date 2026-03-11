#!/usr/bin/env python3
"""
Replace placeholder macros in manuscript.tex with actual metric values
after compute_metrics.py and statistical_analysis.py have run.
"""
import json
import re
from pathlib import Path

BASE = Path(__file__).parent
RESULTS = BASE / "05_evaluation" / "results"
MANUSCRIPT = BASE / "06_manuscript" / "manuscript.tex"


def fmt(v, decimals=4):
    if v is None:
        return "N/A"
    return f"{float(v):.{decimals}f}"


def fmt2(v):
    return fmt(v, 2)


def fmt4(v):
    return fmt(v, 4)


def fmt_p(v):
    if v is None:
        return "N/A"
    v = float(v)
    if v < 0.001:
        return "$<$0.001"
    return f"{v:.3f}"


def load(fname):
    with open(RESULTS / fname) as f:
        return json.load(f)


def main():
    m_n = load("metrics_naive.json")
    m_s = load("metrics_skill.json")
    comp = load("comparison_stats.json")
    boot = load("bootstrap_stats.json")

    ci_n_roc = m_n.get("bootstrap_ci", {}).get("roc_auc", {})
    ci_s_roc = m_s.get("bootstrap_ci", {}).get("roc_auc", {})
    ci_s_bed = m_s.get("bootstrap_ci", {}).get("bedroc_20", {})

    delong = comp.get("delong_auc", {})
    perm_bed = comp.get("permutation_bedroc_20", {})
    mw_s = comp.get("mannwhitney_skill", {})

    replacements = {
        r"\\RAUCNAIVE\b": fmt4(m_n.get("roc_auc")),
        r"\\RAUCSKILL\b": fmt4(m_s.get("roc_auc")),
        r"\\RAUCNAIVECI\b": f"[{fmt4(ci_n_roc.get('ci_low'))}, {fmt4(ci_n_roc.get('ci_high'))}]",
        r"\\RAUCSKILLCI\b": f"[{fmt4(ci_s_roc.get('ci_low'))}, {fmt4(ci_s_roc.get('ci_high'))}]",
        r"\\DELTAAUC\b": f"{float(delong.get('delta_auc', 0)):+.4f}",
        r"\\DELONGP\b": fmt_p(delong.get("p_value")),
        r"\\BEDROC20NAIVE\b": fmt4(m_n.get("bedroc_20")),
        r"\\BEDROC20SKILL\b": fmt4(m_s.get("bedroc_20")),
        r"\\DELTABEDROC\b": f"{float(perm_bed.get('observed_delta', 0)):+.4f}",
        r"\\PBEDROC\b": fmt_p(perm_bed.get("p_value_permutation")),
        r"\\EFONEPERCENTNAIVE\b": fmt2(m_n.get("ef_1pct")),
        r"\\EFONEPERCENTSKILL\b": fmt2(m_s.get("ef_1pct")),
        r"\\DELTAEFONEM\b": f"{float(comp.get('permutation_ef_1pct', {}).get('observed_delta', 0)):+.2f}",
        r"\\AUPRNAIVE\b": fmt4(m_n.get("aupr")),
        r"\\AUPRSKILL\b": fmt4(m_s.get("aupr")),
        r"\\MANNWP\b": fmt_p(mw_s.get("p_value")),
        r"\\MANNWR\b": fmt4(mw_s.get("rank_biserial_r")),
    }

    content = MANUSCRIPT.read_text()
    n_replaced = 0
    for pattern, replacement in replacements.items():
        new_content, n = re.subn(pattern, replacement.replace("\\", "\\\\"), content)
        if n > 0:
            print(f"  {pattern[2:].rstrip(r'\b')} → {replacement}")
        n_replaced += n
        content = new_content

    MANUSCRIPT.write_text(content)
    print(f"\nReplaced {n_replaced} placeholders in {MANUSCRIPT}")

    # Also update SI Table S1
    _update_si_table(m_n, m_s, comp)


def _update_si_table(m_n, m_s, comp):
    """Update placeholder rows in SI.tex."""
    si_path = BASE / "06_manuscript" / "supporting_information" / "SI.tex"
    if not si_path.exists():
        return

    def row(name, key_n, key_s, ci_key=None):
        vn = fmt4(m_n.get(key_n))
        vs = fmt4(m_s.get(key_s))
        ci_n = m_n.get("bootstrap_ci", {}).get(ci_key, {}) if ci_key else {}
        ci_s = m_s.get("bootstrap_ci", {}).get(ci_key, {}) if ci_key else {}
        cin = f"[{fmt4(ci_n.get('ci_low'))}, {fmt4(ci_n.get('ci_high'))}]" if ci_n else "---"
        cis = f"[{fmt4(ci_s.get('ci_low'))}, {fmt4(ci_s.get('ci_high'))}]" if ci_s else "---"
        p_key = f"permutation_{ci_key}" if ci_key else "delong_auc"
        pv = comp.get(p_key, {}).get("p_value_permutation") or comp.get("delong_auc", {}).get("p_value")
        sig = "* " if pv and float(pv) < 0.05 else ""
        return f"  {name} & {vn} & {cin} & {vs} & {cis} & {sig}\\\\"

    rows = [
        row("ROC AUC", "roc_auc", "roc_auc", "roc_auc"),
        row("BEDROC($\\alpha$=20)", "bedroc_20", "bedroc_20", "bedroc_20"),
        row("BEDROC($\\alpha$=80)", "bedroc_80", "bedroc_80"),
        row("EF$_{1\\%}$", "ef_1pct", "ef_1pct", "ef_1pct"),
        row("EF$_{5\\%}$", "ef_5pct", "ef_5pct", "ef_5pct"),
        row("EF$_{10\\%}$", "ef_10pct", "ef_10pct"),
        row("LogAUC", "log_auc", "log_auc", "log_auc"),
        row("pROC AUC (1\\%)", "proc_auc_1pct", "proc_auc_1pct"),
        row("pROC AUC (5\\%)", "proc_auc_5pct", "proc_auc_5pct"),
        row("AUPR", "aupr", "aupr", "aupr"),
        row("RIE($\\alpha$=20)", "rie_20", "rie_20"),
    ]

    si_content = si_path.read_text()
    table_content = "\n".join(rows)
    # Replace placeholder rows (lines with [PLACEHOLDER])
    new_content = re.sub(
        r"(\\midrule\n)(.*?)(\\bottomrule)",
        r"\1" + table_content + r"\n\3",
        si_content,
        flags=re.DOTALL,
        count=1,
    )
    si_path.write_text(new_content)
    print(f"  Updated SI table: {si_path}")


if __name__ == "__main__":
    main()
