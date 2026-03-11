#!/usr/bin/env python3
"""
Naive protocol: Prepare library compounds using OpenBabel defaults.
OpenBabel gen3d → SDF, then SDF → PDBQT (-h --partialcharge gasteiger).
No filtering, no salt removal, no standardization beyond OpenBabel defaults.
"""

import csv
import subprocess
import sys
import tempfile
from pathlib import Path

LIBRARY_SMI = Path(__file__).parent.parent / "library_combined.smi"
OUTPUT_DIR = Path(__file__).parent / "pdbqt"
LOG_FILE = Path(__file__).parent / "preparation_log.csv"
TIMEOUT = 60  # seconds per compound

def prepare_naive(smiles, compound_id, output_dir):
    """Prepare a single compound using OpenBabel naive protocol."""
    output_pdbqt = output_dir / f"{compound_id}.pdbqt"

    try:
        with tempfile.NamedTemporaryFile(suffix=".smi", mode="w", delete=False) as smi_f:
            smi_f.write(f"{smiles}\n")
            smi_path = smi_f.name

        with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as sdf_f:
            sdf_path = sdf_f.name

        # Step 1: SMILES → 3D SDF via gen3d
        result = subprocess.run(
            ["obabel", smi_path, "-O", sdf_path, "--gen3d", "-h"],
            capture_output=True, text=True, timeout=TIMEOUT,
        )
        if result.returncode != 0 or not Path(sdf_path).exists():
            return False, f"gen3d failed: {result.stderr[:200]}"

        # Check SDF has content
        if Path(sdf_path).stat().st_size < 50:
            return False, "Empty SDF"

        # Step 2: SDF → PDBQT with Gasteiger charges
        result = subprocess.run(
            ["obabel", sdf_path, "-O", str(output_pdbqt),
             "-h", "--partialcharge", "gasteiger"],
            capture_output=True, text=True, timeout=TIMEOUT,
        )
        if result.returncode != 0 or not output_pdbqt.exists():
            return False, f"PDBQT conversion failed: {result.stderr[:200]}"

        if output_pdbqt.stat().st_size < 50:
            return False, "Empty PDBQT"

        return True, "OK"

    except subprocess.TimeoutExpired:
        return False, "Timeout"
    except Exception as e:
        return False, str(e)[:200]
    finally:
        for p in [smi_path, sdf_path]:
            try:
                Path(p).unlink(missing_ok=True)
            except Exception:
                pass

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load library
    compounds = []
    with open(LIBRARY_SMI) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                compounds.append((parts[0], parts[1]))
    print(f"Preparing {len(compounds)} compounds (naive protocol)")

    # Process
    results = []
    success = 0
    fail = 0
    for i, (smiles, cid) in enumerate(compounds):
        ok, msg = prepare_naive(smiles, cid, OUTPUT_DIR)
        results.append({"compound_id": cid, "success": ok, "message": msg})
        if ok:
            success += 1
        else:
            fail += 1

        if (i + 1) % 100 == 0:
            print(f"  [{i+1}/{len(compounds)}] Success: {success}, Failed: {fail}")

    # Save log
    with open(LOG_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["compound_id", "success", "message"])
        writer.writeheader()
        writer.writerows(results)

    rate = success / len(compounds) * 100 if compounds else 0
    print(f"\nDone. Success: {success}/{len(compounds)} ({rate:.1f}%)")
    print(f"Failed: {fail}")
    print(f"PDBQT files in: {OUTPUT_DIR}")
    print(f"Log: {LOG_FILE}")

if __name__ == "__main__":
    main()
