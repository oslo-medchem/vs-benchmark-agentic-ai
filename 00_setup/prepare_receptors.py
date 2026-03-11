#!/usr/bin/env python3
"""
Prepare receptor PDBQT files for both naive and skill protocols.
Downloads 7T6S from RCSB, extracts chain R, prepares with OpenBabel (naive)
and Meeko (skill).
"""
import subprocess
import sys
import urllib.request
from pathlib import Path

BASE_DIR = Path(__file__).parent.parent
INPUT_DIR = BASE_DIR / "04_docking"
PDB_ID = "7T6S"
CHAIN = "R"

NAIVE_RECEPTOR = INPUT_DIR / "naive" / "receptor_naive.pdbqt"
SKILL_RECEPTOR = INPUT_DIR / "skill" / "receptor_skill.pdbqt"
WORK_DIR = BASE_DIR / "00_setup" / "tmp_receptor"
WORK_DIR.mkdir(parents=True, exist_ok=True)


def download_pdb():
    pdb_path = WORK_DIR / f"{PDB_ID}.pdb"
    if pdb_path.exists() and pdb_path.stat().st_size > 1000:
        print(f"Using cached {pdb_path}")
        return pdb_path
    url = f"https://files.rcsb.org/download/{PDB_ID}.pdb"
    print(f"Downloading {url} ...")
    urllib.request.urlretrieve(url, pdb_path)
    print(f"Downloaded: {pdb_path.stat().st_size} bytes")
    return pdb_path


def extract_chain_r(pdb_path):
    """Extract chain R, remove HETATM, save as clean PDB."""
    chain_path = WORK_DIR / f"{PDB_ID}_chainR.pdb"
    with open(pdb_path) as fh:
        lines = fh.readlines()
    kept = []
    for line in lines:
        rec = line[:6].strip()
        if rec == "ATOM":
            chain = line[21]
            if chain == CHAIN:
                kept.append(line)
        elif rec == "TER":
            kept.append(line)
    kept.append("END\n")
    with open(chain_path, "w") as fh:
        fh.writelines(kept)
    print(f"Extracted chain R: {len(kept)} lines → {chain_path}")
    return chain_path


def prepare_naive(chain_path):
    """OpenBabel: PDB → PDBQT (add hydrogens, Gasteiger charges)."""
    NAIVE_RECEPTOR.parent.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        ["obabel", str(chain_path), "-O", str(NAIVE_RECEPTOR),
         "-xr",    # rigid receptor
         "--partialcharge", "gasteiger"],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"obabel stderr: {result.stderr[:500]}")
    if NAIVE_RECEPTOR.exists() and NAIVE_RECEPTOR.stat().st_size > 1000:
        print(f"Naive receptor: {NAIVE_RECEPTOR} ({NAIVE_RECEPTOR.stat().st_size} bytes)")
    else:
        # Fallback: try without -xr flag
        result2 = subprocess.run(
            ["obabel", str(chain_path), "-O", str(NAIVE_RECEPTOR),
             "--partialcharge", "gasteiger"],
            capture_output=True, text=True
        )
        print(f"Fallback obabel stdout: {result2.stdout[:200]}")
        print(f"Naive receptor: {NAIVE_RECEPTOR} ({NAIVE_RECEPTOR.stat().st_size if NAIVE_RECEPTOR.exists() else 0} bytes)")


def prepare_skill(chain_path):
    """Meeko: mk_prepare_receptor.py --read_pdb → PDBQT."""
    SKILL_RECEPTOR.parent.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        ["mk_prepare_receptor.py", "--read_pdb", str(chain_path),
         "-p", str(SKILL_RECEPTOR)],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"Meeko stderr: {result.stderr[:500]}")
        # Fallback to obabel if Meeko fails
        print("Falling back to obabel for skill receptor...")
        subprocess.run(
            ["obabel", str(chain_path), "-O", str(SKILL_RECEPTOR),
             "--partialcharge", "gasteiger"],
            capture_output=True, text=True
        )
    if SKILL_RECEPTOR.exists():
        print(f"Skill receptor: {SKILL_RECEPTOR} ({SKILL_RECEPTOR.stat().st_size} bytes)")
    else:
        print("ERROR: Skill receptor not created")


if __name__ == "__main__":
    print("=== Receptor Preparation ===")
    pdb = download_pdb()
    chain = extract_chain_r(pdb)
    print("\n--- Naive protocol (OpenBabel) ---")
    prepare_naive(chain)
    print("\n--- Skill protocol (Meeko) ---")
    prepare_skill(chain)
    print("\nDone.")
