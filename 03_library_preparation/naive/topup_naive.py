#!/usr/bin/env python3
"""
Top-up: Re-prepare the ~1128 null-byte naive PDBQT files using OpenBabel.
Reads library_labels.csv for SMILES, skips files that already have real content.
"""
import csv
import os
import subprocess
import sys
import tempfile
import multiprocessing as mp
from pathlib import Path
from functools import partial

BASE = Path(__file__).parent
LIBRARY_CSV = BASE.parent / "library_labels.csv"
OUTPUT_DIR = BASE / "pdbqt"
TIMEOUT = 60

MAX_WORKERS = max(1, mp.cpu_count() - 2)


def is_null_or_empty(path):
    if not path.exists() or path.stat().st_size == 0:
        return True
    first4 = open(path, "rb").read(4)
    return first4 == b"\x00\x00\x00\x00"


def prepare_one(row, output_dir):
    compound_id = row["compound_id"]
    smiles = row["smiles"]
    output_pdbqt = output_dir / f"{compound_id}.pdbqt"

    if not is_null_or_empty(output_pdbqt):
        return compound_id, "skipped"

    try:
        with tempfile.NamedTemporaryFile(suffix=".smi", mode="w", delete=False) as smi_f:
            smi_f.write(f"{smiles}\n")
            smi_path = smi_f.name

        with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as sdf_f:
            sdf_path = sdf_f.name

        # SMILES → 3D SDF
        r1 = subprocess.run(
            ["obabel", smi_path, "-O", sdf_path, "--gen3d", "-h"],
            capture_output=True, text=True, timeout=TIMEOUT
        )
        os.unlink(smi_path)

        if r1.returncode != 0 or not Path(sdf_path).exists() or Path(sdf_path).stat().st_size < 50:
            try:
                os.unlink(sdf_path)
            except Exception:
                pass
            return compound_id, "failed_gen3d"

        # SDF → PDBQT
        r2 = subprocess.run(
            ["obabel", sdf_path, "-O", str(output_pdbqt),
             "-h", "--partialcharge", "gasteiger"],
            capture_output=True, text=True, timeout=TIMEOUT
        )
        os.unlink(sdf_path)

        if r2.returncode != 0 or not output_pdbqt.exists() or output_pdbqt.stat().st_size < 50:
            return compound_id, "failed_pdbqt"

        return compound_id, "ok"

    except Exception as e:
        return compound_id, f"error: {str(e)[:80]}"


def main():
    rows = list(csv.DictReader(open(LIBRARY_CSV)))
    print(f"Total compounds: {len(rows)}")

    # Find which need top-up
    need_topup = [r for r in rows
                  if is_null_or_empty(OUTPUT_DIR / f"{r['compound_id']}.pdbqt")]
    print(f"Need preparation: {len(need_topup)}")

    if not need_topup:
        print("Nothing to do.")
        return

    prep_func = partial(prepare_one, output_dir=OUTPUT_DIR)
    ok = failed = skipped = 0

    with mp.Pool(MAX_WORKERS) as pool:
        for i, (cid, status) in enumerate(
            pool.imap(prep_func, need_topup, chunksize=20)
        ):
            if status == "ok":
                ok += 1
            elif status == "skipped":
                skipped += 1
            else:
                failed += 1
            if (i + 1) % 200 == 0:
                print(f"  {i+1}/{len(need_topup)} ok={ok} skip={skipped} fail={failed}")

    print(f"\nDone: ok={ok}, skipped={skipped}, failed={failed}")
    total_real = sum(
        1 for r in rows
        if not is_null_or_empty(OUTPUT_DIR / f"{r['compound_id']}.pdbqt")
    )
    print(f"Total real PDBQT files now: {total_real}/{len(rows)}")


if __name__ == "__main__":
    main()
