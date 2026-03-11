# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Purpose

Full virtual screening benchmark comparing **naive defaults** vs **Claude Code skill-guided protocols** for FPR2 (CHEMBL4227, PDB 7T6S chain R). ~100 actives + ~5000 property-matched decoys, docked with Uni-Dock on GPU. Primary metrics: ROC AUC, BEDROC, enrichment factors.

## Environment

```bash
eval "$(conda shell.bash hook)" && conda activate vina_dock
# Required packages: vina, meeko, rdkit, openbabel, unidock, scipy, numpy,
#                   seaborn, scikit-learn, chembl_webresource_client, matplotlib, pandas
```

Hardware: 2x NVIDIA RTX 4500 Ada (24 GB VRAM), CUDA 12.9, Uni-Dock 1.1.3

## Pipeline (in order)

```bash
# 1. Curate actives from ChEMBL (pChEMBL >= 5, then Butina cluster + diversity pick to 100)
cd 01_active_curation
python fetch_fpr2_actives.py     # → actives_raw.csv
python curate_actives.py         # → actives_100.smi, actives_curated.csv

# 2. Generate property-matched decoys (MW/LogP/HBD/HBA/RotBond tolerances, TC < 0.35 vs actives)
cd ../02_decoy_generation
python generate_decoys_local.py  # → decoys_local_5000.smi (uses chembl_pool.smi as source)
python qc_property_match.py      # → qc_property_distributions.pdf, qc_property_report.md
python qc_tanimoto_separation.py # Verifies TC separation

# 3. Merge library, then prepare PDBQT files for both protocols
cd ../03_library_preparation
python prepare_combined_library.py   # → library_combined.smi, library_labels.csv (ground truth)
cd naive && python prepare_naive.py && cd ..   # OpenBabel gen3d → PDBQT (no filters)
cd skill && python prepare_skill.py && cd ..   # RDKit ETKDGv3 + Meeko + PAINS/Brenk filter

# 4. Dock on GPU (run both protocols in parallel)
cd ../04_docking
bash naive/run_unidock.sh &     # receptor_naive.pdbqt, box 20A, exhaustiveness=8
bash skill/run_unidock.sh &     # receptor_skill.pdbqt, box 25A, exhaustiveness=32
wait
python parse_scores.py          # → scores_naive.csv, scores_skill.csv

# 5. Evaluate
cd ../05_evaluation
python compute_metrics.py       # → results/metrics_*.json, roc_*.json
python statistical_analysis.py  # → results/bootstrap_stats.json, comparison_stats.json
python generate_figures.py      # → figures/fig1-7_*.pdf
```

## Data Flow

```
ChEMBL API
  └─ 01_active_curation/actives_100.smi
       └─ 03_library_preparation/library_combined.smi  ←── 02_decoy_generation/decoys_local_5000.smi
            ├─ naive/pdbqt/*.pdbqt
            └─ skill/pdbqt/*.pdbqt (PAINS/Brenk filtered compounds removed)

04_docking/
  ├─ naive/scores_naive.csv  (compound_id → best_score)
  └─ skill/scores_skill.csv

05_evaluation/
  ├─ results/*.json           (ROC AUC, BEDROC, EF1/5/10)
  └─ figures/fig1-7_*.pdf
```

Ground truth labels are in `03_library_preparation/library_labels.csv` (compound_id, label: active/decoy).

## Protocol Differences

| Parameter | Naive | Skill-Guided |
|-----------|-------|-------------|
| Receptor prep | OpenBabel | Meeko (`mk_prepare_receptor.py --read_pdb`) |
| Ligand 3D | OpenBabel `gen3d` | RDKit ETKDGv3, then Meeko |
| Filters | None | RDKit FilterCatalog (PAINS + Brenk) |
| Box size | 20 Å | 25 Å |
| Exhaustiveness | 8 | 32 |

## Key Notes

- `library_labels.csv` is the single source of ground truth — do not regenerate it without re-running the full pipeline from step 3
- Skill-prepared library will have fewer compounds than naive due to PAINS/Brenk filtering; evaluation scripts must handle this size mismatch
- Uni-Dock output score files go to `04_docking/{naive,skill}/results/` before `parse_scores.py` collects them
- `chembl_pool.smi` in `02_decoy_generation/` is a pre-downloaded ~50K compound pool; re-downloading is slow (paginated ChEMBL API)
- See parent CLAUDE.md files for Meeko preparation gotchas and Uni-Dock/Vina CLI notes
