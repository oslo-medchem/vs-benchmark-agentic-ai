#!/usr/bin/env python3
"""
Fetch FPR2 (CHEMBL4227) active compounds from ChEMBL.
Queries bioactivity data with pChEMBL >= 5.0, saves raw CSV.
"""

import csv
import sys
from pathlib import Path
from chembl_webresource_client.new_client import new_client

TARGET_CHEMBL_ID = "CHEMBL4227"
MIN_PCHEMBL = 5.0
OUTPUT = Path(__file__).parent / "actives_raw.csv"

def main():
    activity = new_client.activity

    print(f"Querying ChEMBL for {TARGET_CHEMBL_ID} activities (pChEMBL >= {MIN_PCHEMBL})...")

    # Query all activities for FPR2 with pChEMBL >= 5
    results = activity.filter(
        target_chembl_id=TARGET_CHEMBL_ID,
        pchembl_value__gte=MIN_PCHEMBL,
        standard_relation="=",
        standard_type__in=["EC50", "IC50", "Ki", "Kd"],
    )

    rows = []
    for r in results:
        if r.get("canonical_smiles") and r.get("pchembl_value"):
            rows.append({
                "molecule_chembl_id": r["molecule_chembl_id"],
                "canonical_smiles": r["canonical_smiles"],
                "standard_type": r["standard_type"],
                "standard_value": r["standard_value"],
                "standard_units": r["standard_units"],
                "pchembl_value": float(r["pchembl_value"]),
                "assay_chembl_id": r["assay_chembl_id"],
                "assay_type": r.get("assay_type", ""),
                "assay_description": r.get("assay_description", ""),
                "document_chembl_id": r.get("document_chembl_id", ""),
                "document_year": r.get("document_year", ""),
            })

    print(f"Retrieved {len(rows)} activity records")

    # Save to CSV
    if not rows:
        print("ERROR: No activities found!", file=sys.stderr)
        sys.exit(1)

    fieldnames = list(rows[0].keys())
    with open(OUTPUT, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Saved to {OUTPUT}")

    # Summary statistics
    unique_mols = len(set(r["molecule_chembl_id"] for r in rows))
    pchembl_vals = [r["pchembl_value"] for r in rows]
    types = {}
    for r in rows:
        types[r["standard_type"]] = types.get(r["standard_type"], 0) + 1

    print(f"\nSummary:")
    print(f"  Total records: {len(rows)}")
    print(f"  Unique molecules: {unique_mols}")
    print(f"  pChEMBL range: {min(pchembl_vals):.2f} - {max(pchembl_vals):.2f}")
    print(f"  Activity types: {types}")

if __name__ == "__main__":
    main()
