#!/usr/bin/env python3
"""
QC: Verify property distributions of actives vs decoys match.
Generates 6-panel histograms and KS-test for each property.
"""

import csv
import numpy as np
from pathlib import Path
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ACTIVES_CSV = Path(__file__).parent.parent / "01_active_curation" / "actives_curated.csv"
DECOYS_CSV = Path(__file__).parent / "decoys_curated.csv"
FIG_OUT = Path(__file__).parent / "qc_property_distributions.pdf"
REPORT = Path(__file__).parent / "qc_property_report.md"

PROPERTIES = ["MW", "LogP", "HBD", "HBA", "RotBonds", "NetCharge"]

def load_props(csv_path, smiles_col="std_smiles"):
    """Load properties from CSV."""
    data = {p: [] for p in PROPERTIES}
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            for p in PROPERTIES:
                data[p].append(float(row[p]))
    return data

def main():
    actives = load_props(ACTIVES_CSV)
    decoys = load_props(DECOYS_CSV, smiles_col="smiles")

    print(f"Actives: {len(actives['MW'])} compounds")
    print(f"Decoys: {len(decoys['MW'])} compounds")

    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    axes = axes.flatten()

    report_lines = ["# Property Match QC Report\n"]
    report_lines.append("| Property | KS Statistic | p-value | Match? |")
    report_lines.append("|----------|-------------|---------|--------|")

    all_pass = True
    for i, prop in enumerate(PROPERTIES):
        ax = axes[i]
        a_vals = np.array(actives[prop])
        d_vals = np.array(decoys[prop])

        # KS test
        ks_stat, ks_p = stats.ks_2samp(a_vals, d_vals)
        match = "PASS" if ks_p > 0.05 else "FAIL"
        if ks_p <= 0.05:
            all_pass = False
        report_lines.append(f"| {prop} | {ks_stat:.4f} | {ks_p:.4f} | {match} |")

        # Histogram
        bins = np.histogram_bin_edges(np.concatenate([a_vals, d_vals]), bins="auto")
        ax.hist(a_vals, bins=bins, alpha=0.6, color="blue", label=f"Actives (n={len(a_vals)})", density=True)
        ax.hist(d_vals, bins=bins, alpha=0.4, color="gray", label=f"Decoys (n={len(d_vals)})", density=True)
        ax.set_xlabel(prop, fontsize=10)
        ax.set_ylabel("Density", fontsize=9)
        ax.set_title(f"{prop} (KS p={ks_p:.3f})", fontsize=10)
        ax.legend(fontsize=7)

    plt.tight_layout()
    plt.savefig(FIG_OUT, dpi=300, bbox_inches="tight")
    print(f"Saved figure to {FIG_OUT}")
    plt.close()

    report_lines.append(f"\nOverall: **{'ALL PASS' if all_pass else 'SOME FAIL'}**\n")
    with open(REPORT, "w") as f:
        f.write("\n".join(report_lines))
    print(f"Report saved to {REPORT}")

    # Print results
    for line in report_lines:
        print(line)

if __name__ == "__main__":
    main()
