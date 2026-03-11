# Reproduction Guide

Detailed instructions for reproducing the FPR2 virtual screening benchmark from scratch, or evaluating from pre-computed scores.

## Quick Start: Evaluate from Pre-Computed Scores

If you only want to reproduce the figures, metrics, and statistical analysis (no GPU required):

```bash
conda env create -f environment.yml
conda activate vina_dock

cd 05_evaluation
python compute_metrics.py        # ~30 s
python statistical_analysis.py   # ~2 min (bootstrap)
python generate_figures.py       # ~20 s

# Outputs:
#   results/metrics_{naive,skill}.json
#   results/bootstrap_stats.json, comparison_stats.json
#   figures/fig_main.pdf, figS1-S5.pdf
```

The pre-computed docking scores in `04_docking/{naive,skill}/scores_{naive,skill}.csv` and ground truth labels in `03_library_preparation/library_labels.csv` are included in the repository.

## Full Reproduction (from scratch)

### Prerequisites

**Software:**
- Conda (Miniforge or Miniconda recommended)
- NVIDIA GPU with CUDA >= 12.0 (for Uni-Dock)
- Uni-Dock 1.1.3 (GPU-accelerated AutoDock Vina fork)
- pdfLaTeX + BibTeX (for manuscript compilation)

**Hardware (tested configuration):**
- 2x NVIDIA RTX 4500 Ada (24 GB VRAM each)
- CUDA 12.9
- 64 GB RAM
- The two GPUs run naive and skill docking protocols in parallel

### Environment Setup

```bash
# Create conda environment
conda env create -f environment.yml
conda activate vina_dock

# Verify key packages
python -c "import rdkit; print(f'RDKit {rdkit.__version__}')"
python -c "import meeko; print(f'Meeko {meeko.__version__}')"
which unidock  # Should resolve to conda env bin
```

**Meeko installation note:** Meeko 0.7.1 requires `gemmi` for receptor preparation. Both are installed via pip in `environment.yml`. If you see import errors for gemmi, run `pip install gemmi==0.7.5`.

### Stage 0: Receptor Preparation (~2 min)

```bash
cd 00_setup
python prepare_receptors.py
# Outputs:
#   04_docking/naive/receptor_naive.pdbqt  (OpenBabel)
#   04_docking/skill/receptor_skill.pdbqt  (Meeko)
```

Downloads PDB 7T6S from RCSB, extracts chain R, and prepares receptor PDBQT files using both protocols.

### Stage 1: Active Curation (~5 min)

```bash
cd 01_active_curation
python fetch_fpr2_actives.py     # Queries ChEMBL API → actives_raw.csv
python curate_actives.py         # → actives_100.smi, actives_curated.csv
```

Fetches FPR2 (CHEMBL4227) bioactivity data with pChEMBL >= 5.0, then curates to 100 diverse compounds via SMILES validation, standardisation, deduplication, Butina clustering, and MaxMin diversity picking.

### Stage 2: Decoy Generation (~15 min)

```bash
cd 02_decoy_generation
python generate_decoys_local.py  # → decoys_local_5000.smi
python qc_property_match.py      # → qc_property_distributions.pdf
python qc_tanimoto_separation.py # Verifies TC < 0.35
```

Generates 5000 property-matched decoys from `chembl_pool.smi` (pre-downloaded ~50K ChEMBL compounds). Matching tolerances: MW +/-25, LogP +/-1.0, HBD +/-1, HBA +/-1, RotBond +/-1. Tanimoto cutoff: TC < 0.35 vs all actives.

**Note:** If `chembl_pool.smi` is missing, `generate_decoys_local.py` will download it from ChEMBL (~30 min due to paginated API).

### Stage 3: Library Preparation (~30 min)

```bash
cd 03_library_preparation
python prepare_combined_library.py   # → library_combined.smi, library_labels.csv

# Naive protocol (OpenBabel gen3d → PDBQT)
cd naive && python prepare_naive.py && cd ..

# Skill protocol (RDKit ETKDGv3 + Meeko + PAINS/Brenk filter)
cd skill && python prepare_skill.py && cd ..
```

`library_labels.csv` is the single source of ground truth for evaluation. The skill protocol applies PAINS and Brenk filters, so its PDBQT set will be smaller than naive (~4100 vs ~4900 ligands).

### Stage 4: Docking (~45 min with 2 GPUs)

```bash
cd 04_docking

# Run both protocols in parallel on separate GPUs
CUDA_VISIBLE_DEVICES=0 bash naive/run_unidock.sh &
CUDA_VISIBLE_DEVICES=1 bash skill/run_unidock.sh &
wait

# Parse scores
python parse_scores.py
# → naive/scores_naive.csv, skill/scores_skill.csv
```

**Single GPU:** Run sequentially instead (~90 min total).

**Uni-Dock parameters:**
| Parameter | Naive | Skill |
|-----------|-------|-------|
| Box size | 20 A | 25 A |
| Exhaustiveness | 8 | 32 |
| Scoring | Vina | Vina |
| Num modes | 9 | 9 |

Box centre: (81.893, 115.414, 96.832) — from co-crystallised FUI ligand in 7T6S.

### Stage 5: Evaluation (~3 min)

```bash
cd 05_evaluation
python compute_metrics.py        # ROC AUC, BEDROC, EF, bootstrap CIs
python statistical_analysis.py   # DeLong, permutation, Mann-Whitney, BH FDR
python generate_figures.py       # PDF + SVG figures
```

### Stage 6: Manuscript Compilation (optional)

```bash
cd ..
python update_manuscript_values.py   # Inject metric values into LaTeX

cd 06_manuscript
pdflatex -interaction=nonstopmode manuscript.tex
bibtex manuscript
pdflatex -interaction=nonstopmode manuscript.tex
pdflatex -interaction=nonstopmode manuscript.tex
# → manuscript.pdf

cd supporting_information
pdflatex -interaction=nonstopmode SI.tex
pdflatex -interaction=nonstopmode SI.tex
# → SI.pdf
```

**Or run everything at once:**

```bash
bash run_evaluation_pipeline.sh   # Stages 4 (parse only) through 6
```

## Troubleshooting

### Meeko receptor preparation fails
`mk_prepare_receptor.py` uses `--read_pdb` (not `--pdb`). It writes PDBQT via the `-p` flag. Ensure `gemmi` is installed (`pip install gemmi`).

### Meeko ligand preparation fails
Input SDF must have explicit hydrogens. The skill preparation script adds them with `rdkit.Chem.AddHs(mol, addCoords=True)` before calling Meeko.

### NumPy 2.x compatibility
`np.trapz` was removed in NumPy 2.0. The evaluation scripts use `getattr(np, 'trapezoid', np.trapz)` for backward compatibility.

### Uni-Dock GPU errors
- Ensure CUDA toolkit version matches your driver: `nvidia-smi` shows driver CUDA, `nvcc --version` shows toolkit CUDA
- Uni-Dock 1.1.3 requires CUDA >= 12.0
- If out-of-memory, reduce batch size or run with `--num_modes 5`

### Docking score of 0.0
Compounds that fail to dock (no output PDBQT) are assigned score 0.0 by `parse_scores.py`. This places them at the worst rank, which is the conservative default for benchmark evaluation.
