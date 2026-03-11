# FPR2 Virtual Screening Benchmark: AI Skill-Guided vs Naive Protocols

Controlled benchmark comparing **naive default** protocols against **Claude Code AI skill-guided** protocols for structure-based virtual screening of the formyl peptide receptor 2 (FPR2, CHEMBL4227, PDB [7T6S](https://www.rcsb.org/structure/7T6S) chain R).

100 curated ChEMBL actives + 5000 property-matched decoys, docked with [Uni-Dock](https://github.com/dptech-corp/Uni-Dock) (GPU-accelerated AutoDock Vina). Primary metric: **ROC AUC** (0.658 naive vs 0.710 skill, DeLong p = 0.004).

## Directory Structure

```
vs_benchmark/
├── 00_setup/                      # Receptor preparation (PDB download, chain extraction)
│   ├── prepare_receptors.py       #   Downloads 7T6S, extracts chain R, makes both PDBQT files
│   └── tmp_receptor/              #   Raw and processed PDB files
│
├── 01_active_curation/            # ChEMBL active compound curation
│   ├── fetch_fpr2_actives.py      #   Query ChEMBL for FPR2 bioactivity (pChEMBL >= 5)
│   ├── curate_actives.py          #   Standardise, cluster, diversity-pick → 100 actives
│   ├── actives_raw.csv            #   Raw ChEMBL query results
│   ├── actives_curated.csv        #   Curated compounds with properties
│   └── actives_100.smi            #   Final 100 diverse actives (SMILES)
│
├── 02_decoy_generation/           # Property-matched decoy generation
│   ├── generate_decoys_local.py   #   Match MW/LogP/HBD/HBA/RotBond, TC < 0.35
│   ├── qc_property_match.py       #   QC: property distribution comparison (KS test)
│   ├── qc_tanimoto_separation.py  #   QC: verify Tanimoto separation
│   ├── chembl_pool.smi            #   Pre-downloaded ~50K ChEMBL compound pool (4.4 MB)
│   ├── decoys_local_5000.smi      #   Generated decoys (SMILES)
│   └── qc_property_distributions.pdf  # QC figure
│
├── 03_library_preparation/        # Ligand preparation (both protocols)
│   ├── prepare_combined_library.py    # Merge actives + decoys, generate ground truth
│   ├── library_combined.smi       #   Shuffled combined library
│   ├── library_labels.csv         #   Ground truth (compound_id → active/decoy)
│   ├── naive/                     #   OpenBabel gen3d → PDBQT (no filters)
│   │   └── prepare_naive.py
│   └── skill/                     #   RDKit ETKDGv3 + Meeko + PAINS/Brenk filter
│       └── prepare_skill.py
│
├── 04_docking/                    # GPU docking with Uni-Dock
│   ├── parse_scores.py            #   Extract best score per compound from output PDBQT
│   ├── naive/
│   │   ├── run_unidock.sh         #   Box 20A, exhaustiveness 8
│   │   ├── receptor_naive.pdbqt
│   │   └── scores_naive.csv       #   Pre-computed docking scores
│   └── skill/
│       ├── run_unidock.sh         #   Box 25A, exhaustiveness 32
│       ├── receptor_skill.pdbqt
│       └── scores_skill.csv       #   Pre-computed docking scores
│
├── 05_evaluation/                 # Metrics, statistics, and figures
│   ├── compute_metrics.py         #   ROC AUC, BEDROC, EF, AUPR + bootstrap CIs
│   ├── statistical_analysis.py    #   DeLong, permutation, Mann-Whitney, BH FDR
│   ├── generate_figures.py        #   Publication-quality PDF + SVG figures
│   ├── results/                   #   JSON metrics, bootstrap stats, ROC data, .npy arrays
│   └── figures/                   #   fig_main.pdf, figS1-S5 (PDF + SVG)
│
├── 06_manuscript/                 # LaTeX manuscript and supporting information
│   ├── manuscript.tex
│   ├── manuscript.pdf
│   ├── references.bib
│   ├── cover_letter.tex
│   └── supporting_information/
│       ├── SI.tex
│       └── SI.pdf
│
├── SKILL.md                       # Claude Code skill file used in the experiment
├── CLAUDE.md                      # Claude Code project instructions
├── run_evaluation_pipeline.sh     # One-command pipeline (parse → metrics → figures → manuscript)
├── update_manuscript_values.py    # Inject computed metrics into LaTeX macros
├── environment.yml                # Conda environment specification
├── REPRODUCTION.md                # Detailed reproduction guide with troubleshooting
├── LICENSE                        # MIT licence
└── .zenodo.json                   # Zenodo deposit metadata
```

## Quick Start: Evaluate from Pre-Computed Scores

No GPU required. Pre-computed docking scores are included in the repository.

```bash
conda env create -f environment.yml
conda activate vina_dock

cd 05_evaluation
python compute_metrics.py         # ~30 s → results/*.json
python statistical_analysis.py    # ~2 min → bootstrap_stats.json
python generate_figures.py        # ~20 s → figures/*.pdf
```

## Full Reproduction

See [REPRODUCTION.md](REPRODUCTION.md) for complete instructions, hardware requirements, per-stage runtimes, and troubleshooting.

**Summary of stages:**

| Stage | Script(s) | Time | GPU |
|-------|-----------|------|-----|
| 0. Receptor prep | `00_setup/prepare_receptors.py` | 2 min | No |
| 1. Active curation | `01_active_curation/fetch_fpr2_actives.py`, `curate_actives.py` | 5 min | No |
| 2. Decoy generation | `02_decoy_generation/generate_decoys_local.py`, QC scripts | 15 min | No |
| 3. Library preparation | `03_library_preparation/prepare_combined_library.py`, naive + skill | 30 min | No |
| 4. Docking | `04_docking/{naive,skill}/run_unidock.sh` | 45 min | Yes |
| 5. Evaluation | `05_evaluation/compute_metrics.py`, stats, figures | 3 min | No |
| 6. Manuscript | `update_manuscript_values.py`, pdflatex | 1 min | No |

## Protocol Comparison

| Parameter | Naive | Skill-Guided |
|-----------|-------|--------------|
| Receptor preparation | OpenBabel | Meeko (`mk_prepare_receptor.py`) |
| Ligand 3D generation | OpenBabel `gen3d` | RDKit ETKDGv3 |
| Ligand PDBQT conversion | OpenBabel | Meeko (`mk_prepare_ligand.py`) |
| Chemical filters | None | PAINS + Brenk (RDKit FilterCatalog) |
| Docking box | 20 A | 25 A |
| Exhaustiveness | 8 | 32 |

## What the Skill File Does

The `SKILL.md` file is a Claude Code skill — a structured reference document that an AI coding agent consults when executing virtual screening tasks. It encodes expert best practices for:

- **Target preparation**: protonation, H-bond optimization, grid definition
- **Library preparation**: 3D conformer generation, PAINS/Brenk filtering, tautomer enumeration
- **Docking parameters**: box sizing (10-15 A beyond ligand), higher exhaustiveness
- **Post-processing**: enrichment metric thresholds, hit triage cascades

The "skill-guided" protocol in this benchmark follows the skill's recommendations, while the "naive" protocol uses tool defaults without consulting the skill.

## Key Results

| Metric | Naive | Skill | p-value (BH-adjusted) |
|--------|-------|-------|-----------------------|
| **ROC AUC** | 0.658 | 0.710 | **0.004** |
| BEDROC (alpha=20) | 0.135 | 0.109 | 0.309 |
| EF 1% | 2.99 | 0.00 | 0.067 |
| AUPR | 0.042 | 0.039 | n.s. |

The skill-guided protocol achieves significantly higher ROC AUC (+0.052, DeLong test). Lower early enrichment metrics (BEDROC, EF1%) for the skill protocol are an artefact of PAINS/Brenk filtering: 7 actives removed by filtering receive worst-rank scores, suppressing early enrichment. See the manuscript for full discussion.

## Citation

```bibtex
@software{gani2026fpr2_vs_benchmark,
  author    = {Gani, Osman A.B.S.M.},
  title     = {{FPR2} Virtual Screening Benchmark: Agentic {AI} Skill-Guided vs Naive Protocols},
  year      = {2026},
  publisher = {GitHub},
  url       = {https://github.com/PLACEHOLDER/fpr2-vs-benchmark},
  license   = {MIT}
}
```

## Licence

[MIT](LICENSE)
