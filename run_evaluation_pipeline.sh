#!/bin/bash
# Post-docking evaluation pipeline
# Runs: parse_scores → compute_metrics → statistical_analysis → generate_figures
# Then updates manuscript placeholders with real values

set -e
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
eval "$(conda shell.bash hook)"
conda activate vina_dock

cd "$BASE"

echo "=== [$(date +%H:%M)] Step 1: Parse docking scores ==="
python3 04_docking/parse_scores.py
echo "Scores written:"
wc -l 04_docking/naive/scores_naive.csv 04_docking/skill/scores_skill.csv

echo ""
echo "=== [$(date +%H:%M)] Step 2: Compute metrics (bootstrap CIs) ==="
cd 05_evaluation
python3 compute_metrics.py 2>&1 | tee /tmp/metrics.log

echo ""
echo "=== [$(date +%H:%M)] Step 3: Statistical analysis ==="
python3 statistical_analysis.py 2>&1 | tee /tmp/stats.log

echo ""
echo "=== [$(date +%H:%M)] Step 4: Generate figures (PDF + SVG) ==="
python3 generate_figures.py 2>&1 | tee /tmp/figures.log

echo ""
echo "=== [$(date +%H:%M)] Step 5: Update manuscript placeholders ==="
cd "$BASE"
python3 update_manuscript_values.py 2>&1 | tee /tmp/update_ms.log

echo ""
echo "=== [$(date +%H:%M)] Step 6: Compile manuscript ==="
cd 06_manuscript
pdflatex -interaction=nonstopmode manuscript.tex > /tmp/pdflatex1.log 2>&1
bibtex manuscript >> /tmp/pdflatex1.log 2>&1 || true
pdflatex -interaction=nonstopmode manuscript.tex >> /tmp/pdflatex1.log 2>&1
pdflatex -interaction=nonstopmode manuscript.tex >> /tmp/pdflatex1.log 2>&1

if [ -f manuscript.pdf ]; then
    echo "Manuscript compiled: $(ls -lh manuscript.pdf | awk '{print $5}')"
else
    echo "WARNING: manuscript.pdf not generated, check /tmp/pdflatex1.log"
fi

echo ""
echo "=== DONE: $(date) ==="
echo "Figures: $BASE/05_evaluation/figures/"
echo "Metrics: $BASE/05_evaluation/results/"
echo "Manuscript: $BASE/06_manuscript/manuscript.pdf"
