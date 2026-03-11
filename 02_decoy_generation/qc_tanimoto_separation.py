#!/usr/bin/env python3
"""
QC: Verify Tanimoto separation between decoys and actives.
All decoy-active pairs should have Tc < 0.35.
"""

import csv
import numpy as np
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.DataStructs import BulkTanimotoSimilarity

RDLogger.logger().setLevel(RDLogger.ERROR)

ACTIVES_CSV = Path(__file__).parent.parent / "01_active_curation" / "actives_curated.csv"
DECOYS_CSV = Path(__file__).parent / "decoys_curated.csv"

MAX_TC = 0.35

def main():
    # Load actives
    fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    active_fps = []
    with open(ACTIVES_CSV) as f:
        for row in csv.DictReader(f):
            mol = Chem.MolFromSmiles(row["std_smiles"])
            if mol:
                active_fps.append(fpgen.GetFingerprint(mol))

    # Load decoys
    decoy_fps = []
    decoy_ids = []
    with open(DECOYS_CSV) as f:
        for row in csv.DictReader(f):
            mol = Chem.MolFromSmiles(row["smiles"])
            if mol:
                decoy_fps.append(fpgen.GetFingerprint(mol))
                decoy_ids.append(row["compound_id"])

    print(f"Checking {len(decoy_fps)} decoys against {len(active_fps)} actives...")

    max_tcs = []
    violations = []
    for i, dfp in enumerate(decoy_fps):
        sims = BulkTanimotoSimilarity(dfp, active_fps)
        max_tc = max(sims)
        max_tcs.append(max_tc)
        if max_tc >= MAX_TC:
            violations.append((decoy_ids[i], max_tc))

    max_tcs = np.array(max_tcs)
    print(f"\nMax Tc distribution:")
    print(f"  Min:    {max_tcs.min():.4f}")
    print(f"  Median: {np.median(max_tcs):.4f}")
    print(f"  Mean:   {max_tcs.mean():.4f}")
    print(f"  Max:    {max_tcs.max():.4f}")
    print(f"  Violations (Tc >= {MAX_TC}): {len(violations)}")

    if violations:
        print(f"\nFirst 10 violations:")
        for did, tc in violations[:10]:
            print(f"  {did}: Tc = {tc:.4f}")

    # Result
    if len(violations) == 0:
        print(f"\nPASS: All decoys have Tc < {MAX_TC} vs actives")
    else:
        print(f"\nWARNING: {len(violations)} decoys have Tc >= {MAX_TC}")
        print(f"({len(violations)/len(decoy_fps)*100:.1f}% of decoys)")

if __name__ == "__main__":
    main()
