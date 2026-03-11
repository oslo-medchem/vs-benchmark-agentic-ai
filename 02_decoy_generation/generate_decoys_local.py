#!/usr/bin/env python3
"""
Generate property-matched decoys for FPR2 virtual screening benchmark.

Strategy: Download drug-like SMILES from ChEMBL in bulk via the compound
search endpoint (all small molecules, paginated), save locally, then
match against actives.

For speed, we use a two-phase approach:
  Phase 1: Download ~50K SMILES from ChEMBL (one-time, saved to disk)
  Phase 2: Property-match against actives locally (fast)
"""

import csv
import sys
import random
import os
import numpy as np
from pathlib import Path
from collections import defaultdict

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, SaltRemover
from rdkit.Chem import rdFingerprintGenerator
from rdkit.DataStructs import TanimotoSimilarity, BulkTanimotoSimilarity

RDLogger.logger().setLevel(RDLogger.ERROR)

ACTIVES_CSV = Path(__file__).parent.parent / "01_active_curation" / "actives_curated.csv"
POOL_FILE = Path(__file__).parent / "chembl_pool.smi"
OUTPUT_SMI = Path(__file__).parent / "decoys_local_5000.smi"
OUTPUT_CSV = Path(__file__).parent / "decoys_curated.csv"

DECOYS_PER_ACTIVE = 50
MAX_TC_VS_ACTIVES = 0.35
RANDOM_SEED = 42

TOL_MW = 25
TOL_LOGP = 1.0
TOL_HBD = 1
TOL_HBA = 2
TOL_ROTBONDS = 2


def load_actives():
    actives = []
    with open(ACTIVES_CSV) as f:
        for row in csv.DictReader(f):
            mol = Chem.MolFromSmiles(row["std_smiles"])
            if mol is None:
                continue
            actives.append({
                "compound_id": row["compound_id"],
                "smiles": row["std_smiles"],
                "mol": mol,
                "MW": float(row["MW"]),
                "LogP": float(row["LogP"]),
                "HBD": int(float(row["HBD"])),
                "HBA": int(float(row["HBA"])),
                "RotBonds": int(float(row["RotBonds"])),
                "NetCharge": int(float(row["NetCharge"])),
            })
    return actives


def compute_properties(mol):
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


def download_pool(mw_min, mw_max, target_size=60000):
    """Download a pool of drug-like molecules from ChEMBL."""
    from chembl_webresource_client.new_client import new_client

    print(f"Downloading ChEMBL pool (MW {mw_min}-{mw_max}, target {target_size})...")

    # Get FPR2 actives to exclude
    activity = new_client.activity
    fpr2_mols = activity.filter(target_chembl_id="CHEMBL4227").only(["molecule_chembl_id"])
    exclude_ids = set()
    for r in fpr2_mols:
        exclude_ids.add(r["molecule_chembl_id"])
    print(f"  Excluding {len(exclude_ids)} FPR2-tested compounds")

    # Query in MW sub-ranges for better pagination
    step = 50  # 50 Da windows
    pool = []
    seen_ids = set()

    for mw_lo in range(int(mw_min), int(mw_max), step):
        mw_hi = min(mw_lo + step, mw_max)
        molecule = new_client.molecule
        results = molecule.filter(
            molecule_properties__mw_freebase__gte=mw_lo,
            molecule_properties__mw_freebase__lte=mw_hi,
            molecule_type="Small molecule",
        ).only(["molecule_chembl_id", "molecule_structures"])

        batch_count = 0
        per_bin_limit = max(target_size // ((mw_max - mw_min) // step + 1), 500)

        for r in results:
            batch_count += 1
            if batch_count > per_bin_limit:
                break

            chembl_id = r.get("molecule_chembl_id", "")
            if chembl_id in exclude_ids or chembl_id in seen_ids:
                continue
            seen_ids.add(chembl_id)

            structs = r.get("molecule_structures")
            if not structs:
                continue
            smi = structs.get("canonical_smiles")
            if not smi or len(smi) > 200:
                continue

            pool.append(f"{smi}\t{chembl_id}")

        print(f"  MW {mw_lo}-{mw_hi}: {batch_count} queried, pool now {len(pool)}")
        sys.stdout.flush()

        if len(pool) >= target_size:
            break

    # Save pool
    with open(POOL_FILE, "w") as f:
        f.write("\n".join(pool) + "\n")
    print(f"Saved pool: {len(pool)} molecules to {POOL_FILE}")
    return len(pool)


def load_pool():
    """Load pool from file."""
    pool = []
    with open(POOL_FILE) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                pool.append({"smiles": parts[0], "chembl_id": parts[1]})
    return pool


def main():
    random.seed(RANDOM_SEED)
    np.random.seed(RANDOM_SEED)

    actives = load_actives()
    print(f"Loaded {len(actives)} actives")

    # Determine MW range needed
    mw_vals = [a["MW"] for a in actives]
    mw_min = min(mw_vals) - 60
    mw_max = max(mw_vals) + 60
    print(f"Active MW range: {min(mw_vals):.0f}-{max(mw_vals):.0f}, pool range: {mw_min:.0f}-{mw_max:.0f}")

    # Phase 1: Download pool (or use cached)
    if POOL_FILE.exists() and POOL_FILE.stat().st_size > 1000:
        print(f"Using cached pool from {POOL_FILE}")
    else:
        download_pool(mw_min, mw_max, target_size=60000)

    # Load pool
    raw_pool = load_pool()
    print(f"Pool: {len(raw_pool)} molecules")

    # Phase 2: Compute properties and fingerprints
    print("Computing properties and fingerprints...")
    fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

    active_fps = []
    for a in actives:
        a["fp"] = fpgen.GetFingerprint(a["mol"])
        active_fps.append(a["fp"])

    pool = []
    for i, entry in enumerate(raw_pool):
        mol = Chem.MolFromSmiles(entry["smiles"])
        if mol is None:
            continue
        props = compute_properties(mol)
        fp = fpgen.GetFingerprint(mol)
        pool.append({
            "chembl_id": entry["chembl_id"],
            "smiles": entry["smiles"],
            "mol": mol,
            "fp": fp,
            "props": props,
        })
        if (i + 1) % 10000 == 0:
            print(f"  Processed {i+1}/{len(raw_pool)}...")

    print(f"Valid pool: {len(pool)} molecules")

    # Phase 3: Match decoys to actives
    print("\nMatching decoys to actives...")
    used_inchikeys = set()
    for a in actives:
        ik = Chem.MolToInchiKey(a["mol"])
        if ik:
            used_inchikeys.add(ik[:14])

    all_decoys = []
    decoys_per_active = defaultdict(int)

    for ai, active in enumerate(actives):
        needed = DECOYS_PER_ACTIVE

        candidates = []
        for p in pool:
            props = p["props"]
            if abs(props["MW"] - active["MW"]) > TOL_MW:
                continue
            if abs(props["LogP"] - active["LogP"]) > TOL_LOGP:
                continue
            if abs(props["HBD"] - active["HBD"]) > TOL_HBD:
                continue
            if abs(props["HBA"] - active["HBA"]) > TOL_HBA:
                continue
            if abs(props["RotBonds"] - active["RotBonds"]) > TOL_ROTBONDS:
                continue
            if props["NetCharge"] != active["NetCharge"]:
                continue

            ik = Chem.MolToInchiKey(p["mol"])
            if ik and ik[:14] in used_inchikeys:
                continue

            max_tc = max(TanimotoSimilarity(p["fp"], afp) for afp in active_fps)
            if max_tc >= MAX_TC_VS_ACTIVES:
                continue

            candidates.append((p, ik, max_tc))

        candidates.sort(key=lambda x: x[2])

        picked = 0
        for p, ik, max_tc in candidates:
            if picked >= needed:
                break
            decoy = {
                "matched_to": active["compound_id"],
                "chembl_id": p["chembl_id"],
                "smiles": p["smiles"],
                "max_tc_to_actives": max_tc,
            }
            decoy.update(p["props"])
            all_decoys.append(decoy)
            if ik:
                used_inchikeys.add(ik[:14])
            picked += 1

        decoys_per_active[active["compound_id"]] = picked
        if (ai + 1) % 20 == 0:
            print(f"  Active {ai+1}/{len(actives)}: total decoys = {len(all_decoys)}")

    print(f"\nTotal decoys: {len(all_decoys)}")

    # Report
    insufficient = [(aid, n) for aid, n in decoys_per_active.items() if n < DECOYS_PER_ACTIVE]
    if insufficient:
        print(f"Actives with < {DECOYS_PER_ACTIVE} decoys: {len(insufficient)}")

    # Second pass with widened tolerances
    if len(all_decoys) < 3000:
        print("\nWidening tolerances for second pass...")
        for ai, active in enumerate(actives):
            needed = DECOYS_PER_ACTIVE - decoys_per_active[active["compound_id"]]
            if needed <= 0:
                continue

            for p in pool:
                if needed <= 0:
                    break
                props = p["props"]
                if abs(props["MW"] - active["MW"]) > 50:
                    continue
                if abs(props["LogP"] - active["LogP"]) > 2.0:
                    continue
                if abs(props["HBD"] - active["HBD"]) > 2:
                    continue
                if abs(props["HBA"] - active["HBA"]) > 3:
                    continue

                ik = Chem.MolToInchiKey(p["mol"])
                if ik and ik[:14] in used_inchikeys:
                    continue

                max_tc = max(TanimotoSimilarity(p["fp"], afp) for afp in active_fps)
                if max_tc >= MAX_TC_VS_ACTIVES:
                    continue

                decoy = {
                    "matched_to": active["compound_id"],
                    "chembl_id": p["chembl_id"],
                    "smiles": p["smiles"],
                    "max_tc_to_actives": max_tc,
                }
                decoy.update(p["props"])
                all_decoys.append(decoy)
                if ik:
                    used_inchikeys.add(ik[:14])
                needed -= 1
                decoys_per_active[active["compound_id"]] += 1

        print(f"After widening: {len(all_decoys)} decoys")

    # Trim if needed
    if len(all_decoys) > 5000:
        random.shuffle(all_decoys)
        all_decoys = all_decoys[:5000]

    # Save
    with open(OUTPUT_SMI, "w") as f:
        for i, d in enumerate(all_decoys, 1):
            did = f"DEC_{i:04d}"
            d["compound_id"] = did
            f.write(f"{d['smiles']}\t{did}\t{d['chembl_id']}\n")
    print(f"\nSaved {len(all_decoys)} decoys to {OUTPUT_SMI}")

    csv_fields = ["compound_id", "matched_to", "chembl_id", "smiles",
                  "MW", "LogP", "HBD", "HBA", "RotBonds", "TPSA",
                  "NumHeavyAtoms", "NetCharge", "max_tc_to_actives"]
    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=csv_fields, extrasaction="ignore")
        writer.writeheader()
        for d in all_decoys:
            writer.writerow(d)
    print(f"Saved CSV to {OUTPUT_CSV}")

    tc_vals = [d["max_tc_to_actives"] for d in all_decoys]
    print(f"\nMax Tc to actives: min={min(tc_vals):.3f}, max={max(tc_vals):.3f}, mean={np.mean(tc_vals):.3f}")

if __name__ == "__main__":
    main()
