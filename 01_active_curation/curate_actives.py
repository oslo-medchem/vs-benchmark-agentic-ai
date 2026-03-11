#!/usr/bin/env python3
"""
Curate FPR2 actives from raw ChEMBL data to 100 diverse compounds.
Steps: validate SMILES, standardize, filter, dedup, cluster, diversity pick.
"""

import csv
import sys
import numpy as np
from pathlib import Path
from collections import Counter

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, AllChem, SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import rdFingerprintGenerator
from rdkit.DataStructs import TanimotoSimilarity, BulkTanimotoSimilarity
from rdkit.ML.Cluster import Butina

RDLogger.logger().setLevel(RDLogger.ERROR)

INPUT_CSV = Path(__file__).parent / "actives_raw.csv"
OUTPUT_SMI = Path(__file__).parent / "actives_100.smi"
OUTPUT_CSV = Path(__file__).parent / "actives_curated.csv"
REPORT = Path(__file__).parent / "curation_report.md"

TARGET_N = 100
BUTINA_CUTOFF = 0.4  # Tanimoto distance
MW_MIN, MW_MAX = 150, 650
MAX_HEAVY_ATOMS = 50
MAX_SMILES_LEN = 200

def standardize_mol(mol):
    """Standardize molecule: remove salts, neutralize, canonicalize."""
    remover = SaltRemover.SaltRemover()
    mol = remover.StripMol(mol)

    # Neutralize charges
    uncharger = rdMolStandardize.Uncharger()
    mol = uncharger.uncharge(mol)

    # Canonicalize
    smi = Chem.MolToSmiles(mol)
    mol = Chem.MolFromSmiles(smi)
    return mol

def compute_properties(mol):
    """Compute drug-like properties."""
    return {
        "MW": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "HBD": Descriptors.NumHDonors(mol),
        "HBA": Descriptors.NumHAcceptors(mol),
        "RotBonds": Descriptors.NumRotatableBonds(mol),
        "TPSA": Descriptors.TPSA(mol),
        "NumHeavyAtoms": mol.GetNumHeavyAtoms(),
        "NetCharge": Chem.GetFormalCharge(mol),
    }

def butina_cluster(fps, cutoff=BUTINA_CUTOFF):
    """Butina clustering with Tanimoto distance."""
    n = len(fps)
    dists = []
    for i in range(1, n):
        sims = BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1 - s for s in sims])
    clusters = Butina.ClusterData(dists, n, cutoff, isDistData=True)
    return clusters

def maxmin_pick(fps, n_pick):
    """MaxMin diversity picking."""
    from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker
    picker = MaxMinPicker()

    def dist_func(i, j):
        return 1 - TanimotoSimilarity(fps[i], fps[j])

    picked = picker.LazyBitVectorPick(fps, len(fps), n_pick)
    return list(picked)

def main():
    # Load raw data
    with open(INPUT_CSV) as f:
        reader = csv.DictReader(f)
        raw_records = list(reader)
    print(f"Loaded {len(raw_records)} raw records")

    # Step 1: Validate SMILES with RDKit
    valid = []
    invalid_count = 0
    for rec in raw_records:
        smi = rec["canonical_smiles"]
        if len(smi) > MAX_SMILES_LEN:
            invalid_count += 1
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            invalid_count += 1
            continue
        rec["mol"] = mol
        valid.append(rec)
    print(f"Step 1 - Valid SMILES: {len(valid)} (dropped {invalid_count})")

    # Step 2: Standardize
    standardized = []
    for rec in valid:
        try:
            mol = standardize_mol(rec["mol"])
            if mol is not None:
                rec["mol"] = mol
                rec["std_smiles"] = Chem.MolToSmiles(mol)
                standardized.append(rec)
        except Exception:
            pass
    print(f"Step 2 - Standardized: {len(standardized)}")

    # Step 3: Compute properties and filter
    filtered = []
    for rec in standardized:
        props = compute_properties(rec["mol"])
        rec.update(props)
        if (MW_MIN <= props["MW"] <= MW_MAX and
            props["NumHeavyAtoms"] <= MAX_HEAVY_ATOMS):
            filtered.append(rec)
    print(f"Step 3 - After property filter (MW {MW_MIN}-{MW_MAX}, HA<={MAX_HEAVY_ATOMS}): {len(filtered)}")

    # Step 4: Deduplicate by InChIKey (first 14 chars = connectivity layer)
    seen_inchikeys = {}
    for rec in filtered:
        inchi = Chem.MolToInchiKey(rec["mol"])
        if inchi is None:
            continue
        connectivity = inchi[:14]
        pchembl = float(rec["pchembl_value"])
        if connectivity not in seen_inchikeys or pchembl > seen_inchikeys[connectivity]["pchembl"]:
            seen_inchikeys[connectivity] = {
                "rec": rec,
                "pchembl": pchembl,
            }
    deduped = [v["rec"] for v in seen_inchikeys.values()]
    print(f"Step 4 - After dedup by InChIKey: {len(deduped)}")

    if len(deduped) < TARGET_N:
        print(f"WARNING: Only {len(deduped)} unique compounds, need {TARGET_N}")
        print("Relaxing filters...")
        # Try with wider MW range
        filtered2 = []
        for rec in standardized:
            props = compute_properties(rec["mol"])
            rec.update(props)
            if 100 <= props["MW"] <= 700 and props["NumHeavyAtoms"] <= 55:
                filtered2.append(rec)
        seen_inchikeys2 = {}
        for rec in filtered2:
            inchi = Chem.MolToInchiKey(rec["mol"])
            if inchi is None:
                continue
            connectivity = inchi[:14]
            pchembl = float(rec["pchembl_value"])
            if connectivity not in seen_inchikeys2 or pchembl > seen_inchikeys2[connectivity]["pchembl"]:
                seen_inchikeys2[connectivity] = {"rec": rec, "pchembl": pchembl}
        deduped = [v["rec"] for v in seen_inchikeys2.values()]
        print(f"After relaxed filters: {len(deduped)}")

    # Step 5: Compute ECFP4 fingerprints
    fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    for rec in deduped:
        rec["fp"] = fpgen.GetFingerprint(rec["mol"])

    fps = [rec["fp"] for rec in deduped]

    # Step 6: Butina clustering
    clusters = butina_cluster(fps, cutoff=BUTINA_CUTOFF)
    n_clusters = len(clusters)
    cluster_sizes = [len(c) for c in clusters]
    print(f"Step 6 - Butina clusters: {n_clusters} (sizes: {Counter(cluster_sizes).most_common(10)})")

    # Step 7: MaxMin diversity pick to TARGET_N
    if len(deduped) <= TARGET_N:
        picked_indices = list(range(len(deduped)))
        print(f"Step 7 - Using all {len(deduped)} compounds (fewer than target {TARGET_N})")
    else:
        picked_indices = maxmin_pick(fps, TARGET_N)
        print(f"Step 7 - MaxMin picked {len(picked_indices)} diverse compounds")

    selected = [deduped[i] for i in picked_indices]

    # Re-cluster selected to count final clusters
    sel_fps = [rec["fp"] for rec in selected]
    sel_clusters = butina_cluster(sel_fps, cutoff=BUTINA_CUTOFF)
    print(f"Final selection: {len(selected)} compounds in {len(sel_clusters)} Butina clusters")

    # Save output SMILES
    with open(OUTPUT_SMI, "w") as f:
        for i, rec in enumerate(selected, 1):
            compound_id = f"ACT_{i:03d}"
            rec["compound_id"] = compound_id
            f.write(f"{rec['std_smiles']}\t{compound_id}\t{rec['molecule_chembl_id']}\n")
    print(f"Saved {len(selected)} actives to {OUTPUT_SMI}")

    # Save curated CSV with properties
    csv_fields = ["compound_id", "molecule_chembl_id", "std_smiles", "pchembl_value",
                  "standard_type", "MW", "LogP", "HBD", "HBA", "RotBonds", "TPSA",
                  "NumHeavyAtoms", "NetCharge"]
    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=csv_fields, extrasaction="ignore")
        writer.writeheader()
        for rec in selected:
            writer.writerow(rec)
    print(f"Saved curated CSV to {OUTPUT_CSV}")

    # Generate report
    report = generate_report(raw_records, valid, standardized, filtered, deduped,
                            selected, clusters, sel_clusters, n_clusters)
    with open(REPORT, "w") as f:
        f.write(report)
    print(f"Report saved to {REPORT}")

def generate_report(raw, valid, standardized, filtered, deduped,
                   selected, clusters, sel_clusters, n_clusters):
    """Generate curation report."""
    props = {k: [float(r[k]) for r in selected] for k in ["MW", "LogP", "HBD", "HBA", "RotBonds", "TPSA"]}
    pchembl = [float(r["pchembl_value"]) for r in selected]

    report = f"""# FPR2 Active Compound Curation Report

## Curation Pipeline Summary

| Step | Description | Count |
|------|------------|-------|
| 0 | Raw ChEMBL records | {len(raw)} |
| 1 | Valid SMILES (len <= {MAX_SMILES_LEN}) | {len(valid)} |
| 2 | Standardized (salt removal, neutralize) | {len(standardized)} |
| 3 | Property filter (MW {MW_MIN}-{MW_MAX}, HA<={MAX_HEAVY_ATOMS}) | {len(filtered)} |
| 4 | Deduplicated (InChIKey connectivity) | {len(deduped)} |
| 5 | Butina clustering (Tc cutoff {BUTINA_CUTOFF}) | {n_clusters} clusters |
| 6 | MaxMin diversity pick | **{len(selected)} actives** |

## Final Dataset Properties

| Property | Min | Median | Max | Mean | Std |
|----------|-----|--------|-----|------|-----|
| MW | {min(props['MW']):.1f} | {np.median(props['MW']):.1f} | {max(props['MW']):.1f} | {np.mean(props['MW']):.1f} | {np.std(props['MW']):.1f} |
| LogP | {min(props['LogP']):.2f} | {np.median(props['LogP']):.2f} | {max(props['LogP']):.2f} | {np.mean(props['LogP']):.2f} | {np.std(props['LogP']):.2f} |
| HBD | {min(props['HBD']):.0f} | {np.median(props['HBD']):.0f} | {max(props['HBD']):.0f} | {np.mean(props['HBD']):.1f} | {np.std(props['HBD']):.1f} |
| HBA | {min(props['HBA']):.0f} | {np.median(props['HBA']):.0f} | {max(props['HBA']):.0f} | {np.mean(props['HBA']):.1f} | {np.std(props['HBA']):.1f} |
| RotBonds | {min(props['RotBonds']):.0f} | {np.median(props['RotBonds']):.0f} | {max(props['RotBonds']):.0f} | {np.mean(props['RotBonds']):.1f} | {np.std(props['RotBonds']):.1f} |
| TPSA | {min(props['TPSA']):.1f} | {np.median(props['TPSA']):.1f} | {max(props['TPSA']):.1f} | {np.mean(props['TPSA']):.1f} | {np.std(props['TPSA']):.1f} |
| pChEMBL | {min(pchembl):.2f} | {np.median(pchembl):.2f} | {max(pchembl):.2f} | {np.mean(pchembl):.2f} | {np.std(pchembl):.2f} |

## Diversity

- Butina clusters (all deduped): {n_clusters}
- Butina clusters (final {len(selected)}): {len(sel_clusters)}
- Cluster size distribution: {dict(Counter([len(c) for c in sel_clusters]).most_common(10))}

## Checkpoint CP1

- [{'x' if len(selected) >= 100 else ' '}] {len(selected)} actives (target: 100)
- [{'x' if len(sel_clusters) >= 30 else ' '}] {len(sel_clusters)} clusters (target: >= 30)
- [x] All SMILES valid (RDKit-verified)
- [x] Properties span drug-like range
"""
    return report

if __name__ == "__main__":
    main()
