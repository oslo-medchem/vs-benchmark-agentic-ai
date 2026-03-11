#!/usr/bin/env python3
"""
Merge actives + decoys into a single shuffled library with ground truth labels.
"""

import csv
import random
from pathlib import Path

ACTIVES_SMI = Path(__file__).parent.parent / "01_active_curation" / "actives_100.smi"
DECOYS_SMI = Path(__file__).parent.parent / "02_decoy_generation" / "decoys_local_5000.smi"
LABELS_CSV = Path(__file__).parent / "library_labels.csv"
LIBRARY_SMI = Path(__file__).parent / "library_combined.smi"

SEED = 42

def main():
    random.seed(SEED)

    # Load actives
    actives = []
    with open(ACTIVES_SMI) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                actives.append({
                    "smiles": parts[0],
                    "compound_id": parts[1],
                    "label": "active",
                    "source_id": parts[2] if len(parts) > 2 else "",
                })
    print(f"Loaded {len(actives)} actives")

    # Load decoys
    decoys = []
    with open(DECOYS_SMI) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                decoys.append({
                    "smiles": parts[0],
                    "compound_id": parts[1],
                    "label": "decoy",
                    "source_id": parts[2] if len(parts) > 2 else "",
                })
    print(f"Loaded {len(decoys)} decoys")

    # Combine and shuffle
    library = actives + decoys
    random.shuffle(library)
    print(f"Combined library: {len(library)} compounds (shuffled, seed={SEED})")

    # Save labels CSV
    with open(LABELS_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["compound_id", "label", "smiles", "source_id"])
        writer.writeheader()
        for entry in library:
            writer.writerow(entry)
    print(f"Labels saved to {LABELS_CSV}")

    # Save combined SMILES
    with open(LIBRARY_SMI, "w") as f:
        for entry in library:
            f.write(f"{entry['smiles']}\t{entry['compound_id']}\n")
    print(f"Library SMILES saved to {LIBRARY_SMI}")

    # Summary
    n_act = sum(1 for e in library if e["label"] == "active")
    n_dec = sum(1 for e in library if e["label"] == "decoy")
    print(f"\nSummary: {n_act} actives + {n_dec} decoys = {len(library)} total")
    print(f"Active ratio: {n_act/len(library)*100:.1f}%")

if __name__ == "__main__":
    main()
