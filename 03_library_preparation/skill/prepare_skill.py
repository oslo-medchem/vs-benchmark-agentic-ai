#!/usr/bin/env python3
"""
Skill-guided protocol: RDKit ETKDGv3 + Meeko 0.7.1 ligand preparation.
Applies PAINS + Brenk filtering. Writes PDBQT to pdbqt/ directory.
"""
import csv
import os
import sys
import subprocess
import tempfile
import multiprocessing as mp
from pathlib import Path
from functools import partial

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

# Paths
BASE = Path(__file__).parent
LIBRARY_CSV = BASE.parent / "library_labels.csv"
OUTPUT_DIR = BASE / "pdbqt"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_FILE = BASE / "preparation_log.csv"
FILTER_REPORT = BASE / "filter_report.csv"

MAX_WORKERS = max(1, mp.cpu_count() - 2)
TIMEOUT = 90  # seconds

# PAINS + Brenk filter catalog
def build_filter():
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    return FilterCatalog(params)


def prepare_one(row, output_dir, apply_filter=True):
    compound_id = row["compound_id"]
    smiles = row["smiles"]
    output_pdbqt = output_dir / f"{compound_id}.pdbqt"

    # Skip if already prepared and non-zero
    if output_pdbqt.exists() and output_pdbqt.stat().st_size > 100:
        first4 = open(output_pdbqt, "rb").read(4)
        if first4 != b"\x00\x00\x00\x00":
            return compound_id, "skipped", 0, ""

    catalog = build_filter()

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return compound_id, "failed", 0, "invalid_smiles"

    # PAINS/Brenk filter
    if apply_filter:
        entry = catalog.GetFirstMatch(mol)
        if entry is not None:
            filter_name = entry.GetDescription()
            return compound_id, "filtered", 0, filter_name

    # Add explicit hydrogens
    mol_h = Chem.AddHs(mol, addCoords=False)

    # Generate 3D conformer (ETKDGv3)
    params_etkdg = AllChem.ETKDGv3()
    params_etkdg.randomSeed = 42
    params_etkdg.numThreads = 1
    result = AllChem.EmbedMolecule(mol_h, params_etkdg)
    if result == -1:
        # Fallback: ETKDG
        params_etkdg2 = AllChem.ETKDGv3()
        params_etkdg2.randomSeed = 0
        result = AllChem.EmbedMolecule(mol_h, params_etkdg2)
        if result == -1:
            return compound_id, "failed", 0, "embed_failed"

    # Optimize geometry
    AllChem.MMFFOptimizeMolecule(mol_h, maxIters=2000)

    # Write 3D SDF
    with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False, mode="w") as sdf_f:
        sdf_path = sdf_f.name
    writer = Chem.SDWriter(sdf_path)
    writer.write(mol_h)
    writer.close()

    # Meeko: SDF → PDBQT
    try:
        result2 = subprocess.run(
            ["mk_prepare_ligand.py", "-i", sdf_path,
             "-o", str(output_pdbqt)],
            capture_output=True, text=True, timeout=TIMEOUT
        )
        os.unlink(sdf_path)
        if result2.returncode != 0 or not output_pdbqt.exists() or output_pdbqt.stat().st_size < 50:
            # Fallback: obabel
            result3 = subprocess.run(
                ["obabel", sdf_path, "-O", str(output_pdbqt),
                 "--partialcharge", "gasteiger"],
                capture_output=True, text=True, timeout=TIMEOUT
            )
            if output_pdbqt.exists() and output_pdbqt.stat().st_size > 50:
                return compound_id, "ok_obabel_fallback", output_pdbqt.stat().st_size, ""
            return compound_id, "failed", 0, f"meeko: {result2.stderr[:100]}"
    except Exception as e:
        try:
            os.unlink(sdf_path)
        except Exception:
            pass
        return compound_id, "failed", 0, str(e)[:100]

    size = output_pdbqt.stat().st_size if output_pdbqt.exists() else 0
    return compound_id, "ok", size, ""


def main():
    print(f"Workers: {MAX_WORKERS}")

    # Load library
    rows = list(csv.DictReader(open(LIBRARY_CSV)))
    print(f"Total compounds: {len(rows)}")

    # Process
    prep_func = partial(prepare_one, output_dir=OUTPUT_DIR, apply_filter=True)

    log_rows = []
    filter_rows = []
    ok = failed = filtered = skipped = 0

    with mp.Pool(MAX_WORKERS) as pool:
        for i, result in enumerate(pool.imap(prep_func, rows, chunksize=20)):
            compound_id, status, size, msg = result
            log_rows.append({"compound_id": compound_id, "status": status,
                              "pdbqt_size": size, "note": msg})
            if status in ("ok", "ok_obabel_fallback", "skipped"):
                ok += 1
            elif status == "filtered":
                filtered += 1
                filter_rows.append({"compound_id": compound_id, "filter_hit": msg})
            else:
                failed += 1

            if (i + 1) % 200 == 0:
                print(f"  {i+1}/{len(rows)} ok={ok} filtered={filtered} failed={failed}")

    print(f"\nDone: {ok} ok, {filtered} filtered, {failed} failed")

    # Write log
    with open(LOG_FILE, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["compound_id", "status", "pdbqt_size", "note"])
        w.writeheader()
        w.writerows(log_rows)

    # Write filter report
    with open(FILTER_REPORT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["compound_id", "filter_hit"])
        w.writeheader()
        w.writerows(filter_rows)

    # Summary
    pdbqt_count = len(list(OUTPUT_DIR.glob("*.pdbqt")))
    print(f"PDBQT files in output dir: {pdbqt_count}")


if __name__ == "__main__":
    main()
