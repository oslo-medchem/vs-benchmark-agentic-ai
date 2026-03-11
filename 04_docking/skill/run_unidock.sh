#!/usr/bin/env bash
# Skill-guided protocol: Uni-Dock docking
# Box 25Å, exhaustiveness=32, Vina scoring
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VS_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"

RECEPTOR="$SCRIPT_DIR/receptor_skill.pdbqt"
PDBQT_DIR="$VS_DIR/03_library_preparation/skill/pdbqt"
RESULTS_DIR="$SCRIPT_DIR/results"
LIGAND_LIST="$SCRIPT_DIR/ligand_list_skill.txt"

mkdir -p "$RESULTS_DIR"

# Box center from FUI crystal pose in 7T6S chain R
CENTER_X=81.893
CENTER_Y=115.414
CENTER_Z=96.832
SIZE=25
EXHAUSTIVENESS=32
NUM_MODES=9

# Build ligand index file
echo "Building ligand list..."
find "$PDBQT_DIR" -name "*.pdbqt" -size +50c | sort > "$LIGAND_LIST"
N=$(wc -l < "$LIGAND_LIST")
echo "Ligands to dock: $N"

echo "Starting Uni-Dock (skill protocol)..."
START_TIME=$(date +%s)

unidock \
  --receptor "$RECEPTOR" \
  --ligand_index "$LIGAND_LIST" \
  --center_x $CENTER_X \
  --center_y $CENTER_Y \
  --center_z $CENTER_Z \
  --size_x $SIZE \
  --size_y $SIZE \
  --size_z $SIZE \
  --exhaustiveness $EXHAUSTIVENESS \
  --num_modes $NUM_MODES \
  --dir "$RESULTS_DIR" \
  --scoring vina \
  2>&1 | tee "$SCRIPT_DIR/unidock_skill.log"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
echo "Elapsed: ${ELAPSED}s" | tee -a "$SCRIPT_DIR/unidock_skill.log"
echo "Skill docking complete."
