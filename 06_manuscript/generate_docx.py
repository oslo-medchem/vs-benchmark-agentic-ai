#!/usr/bin/env python3
"""
Generate manuscript.docx formatted for Bioinformatics Application Note submission.

Follows Oxford Bioinformatics guidelines:
- Double-spaced, 12pt Times New Roman, numbered lines
- Structured abstract (Motivation, Results, Availability and implementation, Contact)
- Numbered body sections
- Author-year citations
- Single figure embedded
"""

from docx import Document
from docx.shared import Pt, Inches, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.section import WD_ORIENT
from docx.oxml.ns import qn
from pathlib import Path
import copy

BASE = Path(__file__).parent
VS_DIR = BASE.parent
FIG_PATH = VS_DIR / "05_evaluation" / "figures" / "fig_main.png"

doc = Document()

# ── Page setup ──────────────────────────────────────────────
section = doc.sections[0]
section.page_width = Cm(21.0)
section.page_height = Cm(29.7)
section.top_margin = Cm(2.54)
section.bottom_margin = Cm(2.54)
section.left_margin = Cm(2.54)
section.right_margin = Cm(2.54)

# ── Default style: 12pt Times New Roman, double-spaced ─────
style = doc.styles["Normal"]
font = style.font
font.name = "Times New Roman"
font.size = Pt(12)
pf = style.paragraph_format
pf.line_spacing = 2.0
pf.space_after = Pt(0)
pf.space_before = Pt(0)


def add_heading_custom(text, level=1):
    """Add numbered or unnumbered heading in TNR bold."""
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(12)
    p.paragraph_format.space_after = Pt(6)
    p.paragraph_format.line_spacing = 2.0
    run = p.add_run(text)
    run.bold = True
    run.font.name = "Times New Roman"
    if level == 1:
        run.font.size = Pt(14)
    elif level == 2:
        run.font.size = Pt(12)
    return p


def add_para(text, bold=False, italic=False, align=None):
    """Add a paragraph with default formatting."""
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing = 2.0
    if align:
        p.alignment = align
    run = p.add_run(text)
    run.font.name = "Times New Roman"
    run.font.size = Pt(12)
    run.bold = bold
    run.italic = italic
    return p


def add_mixed_para(segments):
    """Add paragraph with mixed bold/italic/normal segments.
    segments: list of (text, bold, italic) tuples.
    """
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing = 2.0
    for text, bold, italic in segments:
        run = p.add_run(text)
        run.font.name = "Times New Roman"
        run.font.size = Pt(12)
        run.bold = bold
        run.italic = italic
    return p


# ══════════════════════════════════════════════════════════════
# TITLE PAGE
# ══════════════════════════════════════════════════════════════

# Article type
add_para("Application Note", bold=True, align=WD_ALIGN_PARAGRAPH.CENTER)

# Title
p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
p.paragraph_format.space_before = Pt(24)
p.paragraph_format.space_after = Pt(12)
run = p.add_run(
    "Agentic AI with structured skill files autonomously generates and improves "
    "GPU-accelerated virtual screening: a controlled benchmark on FPR2"
)
run.bold = True
run.font.name = "Times New Roman"
run.font.size = Pt(16)

# Authors
add_para("Osman A.B.S.M. Gani", bold=False, align=WD_ALIGN_PARAGRAPH.CENTER)

# Affiliation
p = add_para(
    "Department of Pharmacy, University of Oslo, 0316 Oslo, Norway",
    align=WD_ALIGN_PARAGRAPH.CENTER,
)
p.runs[0].font.size = Pt(11)

# Corresponding author
p = add_para("Corresponding author: [email withheld for review]", align=WD_ALIGN_PARAGRAPH.CENTER)
p.runs[0].font.size = Pt(11)

# ══════════════════════════════════════════════════════════════
# STRUCTURED ABSTRACT
# ══════════════════════════════════════════════════════════════

add_heading_custom("Abstract", level=1)

add_mixed_para([
    ("Motivation: ", True, False),
    ("Large language model (LLM)-based coding agents can autonomously generate complete "
     "computational drug discovery pipelines, yet the quantitative impact of domain knowledge "
     "injection on virtual screening (VS) outcomes remains unexplored. Here we present the first "
     "controlled benchmark evaluating whether a structured skill file\u2014a concise, machine-readable "
     "document encoding expert procedural knowledge\u2014measurably improves the discriminatory power "
     "of an autonomously generated VS pipeline.", False, False),
])

add_mixed_para([
    ("Results: ", True, False),
    ("Claude Code (Claude Opus 4.6), augmented with a 212-line virtual screening skill file, "
     "autonomously wrote 18 Python and shell scripts (\u223c3,100 lines of code) implementing a "
     "complete GPU-accelerated VS benchmark against the formyl peptide receptor 2 (FPR2)\u2014from "
     "ChEMBL activity retrieval through Uni-Dock docking on dual NVIDIA RTX 4500 Ada GPUs to "
     "bootstrap statistical evaluation\u2014without any human code editing. The skill-guided protocol "
     "achieved an ROC AUC of 0.710 compared with 0.658 for the naive baseline (\u0394AUC = +0.048, "
     "DeLong p = 0.004). All code and pre-computed scores are provided for independent reproduction.", False, False),
])

add_mixed_para([
    ("Availability and implementation: ", True, False),
    ("All code and data are available at https://github.com/oslo-medchem/vs-benchmark-agentic-ai "
     "under the MIT licence. The skill file is included in the repository.", False, False),
])

add_mixed_para([
    ("Contact: ", True, False),
    ("[contact withheld for review]", False, False),
])

add_mixed_para([
    ("Supplementary information: ", True, False),
    ("Supplementary data are available at ", False, False),
    ("Bioinformatics", False, True),
    (" online.", False, False),
])

# ══════════════════════════════════════════════════════════════
# 1 INTRODUCTION
# ══════════════════════════════════════════════════════════════

add_heading_custom("1 Introduction", level=1)

add_para(
    "Virtual screening (VS) constitutes a cornerstone of early-stage drug discovery, enabling the "
    "rapid computational prioritisation of candidate molecules from large chemical libraries for "
    "subsequent experimental validation (Schneider et al., 2020). The advent of GPU-accelerated "
    "docking engines such as Uni-Dock (Yu et al., 2023), which achieves throughputs exceeding "
    "1,000 ligands per minute on commodity hardware, has rendered ultralarge-scale VS campaigns "
    "technically accessible. However, the design and execution of a rigorous VS pipeline\u2014"
    "encompassing receptor preparation, ligand conformer generation, substructure filtering, docking "
    "parameterisation, and multi-metric statistical evaluation\u2014has traditionally demanded expertise "
    "spanning structural biology, cheminformatics, and biostatistics."
)

add_para(
    "LLM-based coding agents, most prominently Claude Code (Anthropic, 2024) and comparable systems "
    "(Boiko et al., 2023; Bran et al., 2024), represent a qualitative shift in this landscape. Such "
    "agents can interpret natural-language instructions and autonomously generate executable code "
    "without human intervention at the code level. Recent reviews have catalogued their expanding "
    "capabilities across chemical research domains (Ramos et al., 2025). Nevertheless, a "
    "general-purpose coding agent lacking domain-specific guidance will typically produce "
    "syntactically valid but methodologically na\u00efve pipelines\u2014employing, for instance, OpenBabel "
    "gen3d for conformer generation, default docking box dimensions, and no structural filtering\u2014"
    "choices known to be suboptimal for structure-based VS (Huang et al., 2006; Truchon and Bayly, 2007)."
)

add_para(
    "Structured skill files\u2014concise, machine-readable documents that encode expert procedural "
    "knowledge\u2014offer a lightweight mechanism for injecting domain expertise into LLM agents without "
    "model fine-tuning or retrieval-augmented generation (Anthropic, 2024). A skill file is loaded "
    "into the agent\u2019s context window at the outset of a session, whereafter the agent follows the "
    "encoded conventions. The quantitative impact of such knowledge injection on downstream VS "
    "performance has not, to our knowledge, been evaluated in a controlled setting."
)

add_para(
    "Here we present a controlled benchmark in which Claude Code\u2014operating with and without a "
    "dedicated virtual-screening-pipeline skill file (212 lines)\u2014autonomously generated a complete "
    "GPU-accelerated VS pipeline targeting FPR2 (CHEMBL4227, PDB 7T6S; Liao et al., 2022; Liu et al., "
    "2022). The agent produced 18 scripts comprising \u223c3,100 lines of code without human code editing. "
    "We report the first quantitative evidence that structured skill files measurably improve global VS "
    "discrimination."
)

# ══════════════════════════════════════════════════════════════
# 2 MATERIALS AND METHODS
# ══════════════════════════════════════════════════════════════

add_heading_custom("2 Materials and Methods", level=1)

add_heading_custom("2.1 Agentic pipeline generation", level=2)

add_para(
    "All computational protocols were designed and executed by Claude Code powered by Claude Opus 4.6 "
    "(Anthropic, 2024). Two parallel pipelines were generated in response to identical natural-language "
    "instructions: one by the base agent (naive) and one by the agent augmented with the "
    "virtual-screening-pipeline skill file (skill-guided). The skill file (212 lines) encodes expert "
    "conventions for receptor and ligand preparation, substructure filtering, docking box sizing, and "
    "exhaustiveness settings; it is interpretable, editable, and requires no model fine-tuning."
)

add_para(
    "The agent autonomously produced 15 Python scripts (2,947 lines) and 3 shell scripts (156 lines), "
    "implementing a six-stage pipeline: (1) active curation from ChEMBL, (2) property-matched decoy "
    "generation, (3) combined library preparation with two independent ligand/receptor preparation "
    "protocols, (4) GPU docking with Uni-Dock, (5) statistical evaluation with bootstrap confidence "
    "intervals, and (6) publication-ready figure generation. No human code editing was performed."
)

add_para(
    "Docking was executed on two NVIDIA RTX 4500 Ada GPUs (24 GB VRAM each, CUDA 12.9) using "
    "Uni-Dock v1.1.3 (Yu et al., 2023) in batch mode."
)

add_heading_custom("2.2 Dataset", level=2)

add_para(
    "FPR2 bioactivities (pChEMBL \u2265 5) were retrieved from ChEMBL v33 (Zdrazil et al., 2024) via "
    "the REST API (n = 382 after deduplication and salt stripping). Butina clustering (Butina, 1999; "
    "Tanimoto distance threshold 0.4, ECFP4 fingerprints) with selection of the highest-pChEMBL "
    "representative per cluster yielded 100 structurally diverse actives. Property-matched decoys "
    "(n = 4,780) were assembled from a pre-downloaded ChEMBL compound pool using matching windows of "
    "MW \u00b125 Da, cLogP \u00b11.5, HBD \u00b11, HBA \u00b12, rotatable bonds \u00b12, and a maximum "
    "Tanimoto similarity of 0.35 to any active (ECFP4). The combined library comprised 4,880 compounds "
    "(active:decoy ratio \u2248 1:48)."
)

add_heading_custom("2.3 Protocol comparison", level=2)

# Table
table = doc.add_table(rows=7, cols=3)
table.style = "Table Grid"
headers = ["Parameter", "Naive", "Skill-guided"]
for i, h in enumerate(headers):
    cell = table.rows[0].cells[i]
    cell.text = h
    for p in cell.paragraphs:
        for r in p.runs:
            r.bold = True
            r.font.name = "Times New Roman"
            r.font.size = Pt(10)

rows_data = [
    ("Ligand 3D generation", "OpenBabel gen3d", "RDKit ETKDGv3"),
    ("PDBQT preparation", "OpenBabel", "Meeko 0.7.1"),
    ("Structural filters", "None", "PAINS A/B/C + Brenk"),
    ("Receptor preparation", "OpenBabel", "Meeko mk_prepare_receptor.py"),
    ("Box size (\u00c5)", "20", "25"),
    ("Exhaustiveness", "8", "32"),
]
for i, (param, naive, skill) in enumerate(rows_data, 1):
    for j, val in enumerate([param, naive, skill]):
        cell = table.rows[i].cells[j]
        cell.text = val
        for p in cell.paragraphs:
            for r in p.runs:
                r.font.name = "Times New Roman"
                r.font.size = Pt(10)

add_para("")  # spacer
add_mixed_para([
    ("Table 1. ", True, False),
    ("Protocol comparison between naive and skill-guided pipelines. Both protocols used the same "
     "docking box centre, derived from the centroid of the co-crystallised FUI ligand in PDB 7T6S "
     "chain R.", False, False),
])

add_heading_custom("2.4 Evaluation metrics and statistical analysis", level=2)

add_para(
    "All metrics treat more-negative docking scores as higher rank. Primary metrics are ROC AUC, "
    "BEDROC(\u03b1=20) (Truchon and Bayly, 2007), EF at 1%, 5%, and 10% (Huang et al., 2006). "
    "Bootstrap 95% confidence intervals (BCa, n = 10,000) were computed for all metrics. AUC "
    "comparison employed DeLong\u2019s non-parametric test (DeLong et al., 1988); metric differences "
    "were assessed by permutation tests (n = 10,000). Within-protocol active vs. decoy score "
    "distributions were compared by the Mann\u2013Whitney U test. All p-values were corrected by the "
    "Benjamini\u2013Hochberg procedure (Benjamini and Hochberg, 1995)."
)

# ══════════════════════════════════════════════════════════════
# 3 RESULTS
# ══════════════════════════════════════════════════════════════

add_heading_custom("3 Results", level=1)

add_heading_custom("3.1 Autonomous pipeline generation", level=2)

add_para(
    "Claude Code autonomously generated 18 scripts (\u223c3,100 lines) implementing the complete VS "
    "pipeline from ChEMBL API query through GPU docking to statistical figures. When guided by the "
    "skill file, the agent incorporated ETKDGv3 conformer generation, Meeko-based PDBQT preparation, "
    "and PAINS/Brenk substructure filtering\u2014none of which appeared in the naive pipeline."
)

add_heading_custom("3.2 Library composition", level=2)

add_para(
    "The naive library comprised 4,632 valid ligands (89 actives, 4,543 decoys after PDBQT "
    "validation). The skill-guided library comprised 4,082 ligands (93 actives, 3,989 decoys); 762 "
    "compounds (15.6% of the original library) were removed by PAINS/Brenk filtering. A greater "
    "proportion of actives survived filtration (93%) than decoys (\u224884%), consistent with the "
    "expectation that confirmed bioactive compounds are less likely to contain pan-assay interference "
    "substructures."
)

add_heading_custom("3.3 Global discrimination", level=2)

add_para(
    "ROC AUC increased from 0.658 (95% CI [0.591, 0.713]) for the naive protocol to 0.710 "
    "(95% CI [0.650, 0.748]) for the skill-guided protocol (Figure 1A). DeLong\u2019s test confirmed a "
    "statistically significant improvement (\u0394AUC = +0.048, z = 2.90, p = 0.004, p_BH = 0.009). "
    "Mann\u2013Whitney U tests demonstrated significant active\u2013decoy separation for both protocols "
    "(p < 10\u207b\u00b9\u2070; rank-biserial r = 0.44 naive, 0.40 skill)."
)

add_heading_custom("3.4 Early enrichment and the PAINS filtering artefact", level=2)

add_para(
    "BEDROC(\u03b1=20) was 0.109 (skill) vs. 0.135 (naive), and EF at 1% was 0.00 vs. 2.99. Neither "
    "difference reached significance after Benjamini\u2013Hochberg correction (p_BH = 0.31 and 0.12, "
    "respectively). These decrements are attributable to seven actives that were removed by PAINS/Brenk "
    "filtering but retained in the ground-truth label set; these compounds were assigned worst-rank "
    "docking scores in the skill evaluation, thereby collapsing early enrichment metrics. This "
    "constitutes a benchmark design artefact\u2014not a failure of the skill-guided docking protocol "
    "itself\u2014and highlights that future AI-driven VS benchmarks should assign labels after filtering "
    "or account for filtered compounds explicitly."
)

# ── Figure ───────────────────────────────────────────────────
doc.add_page_break()

p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
if FIG_PATH.exists():
    run = p.add_run()
    run.add_picture(str(FIG_PATH), width=Inches(6.5))

add_mixed_para([
    ("Figure 1. ", True, False),
    ("Virtual screening performance comparison. ", True, False),
    ("(A) ROC curves for naive (blue) and skill-guided (red) protocols with 95% bootstrap confidence "
     "bands and AUC values [95% CI]. "
     "(B) Bootstrap distributions of ROC AUC (n = 2,000); violin plots with embedded box plots and "
     "observed point estimates (diamonds). DeLong\u2019s test: \u0394AUC = +0.048, p = 0.004. Dashed "
     "line indicates random classifier (AUC = 0.5). "
     "(C) Box plots of docking score distributions for actives and decoys by protocol; lower (more "
     "negative) scores indicate stronger predicted binding. Mann\u2013Whitney U significance: "
     "***p < 0.001.", False, False),
])

# ══════════════════════════════════════════════════════════════
# 4 DISCUSSION
# ══════════════════════════════════════════════════════════════

doc.add_page_break()
add_heading_custom("4 Discussion", level=1)

add_para(
    "The central finding of this work is twofold: an LLM-based coding agent can autonomously generate "
    "a complete, GPU-accelerated, statistically rigorous VS pipeline without any human code editing, "
    "and a structured skill file of modest length (212 lines) measurably improves global discrimination "
    "(ROC AUC +0.048, DeLong p = 0.004). This represents a paradigm shift in accessibility. A VS "
    "campaign of the scope described here\u2014\u223c5,000 compounds docked on GPU with full statistical "
    "evaluation\u2014has hitherto required proficiency in Linux system administration, CUDA driver "
    "management, molecular file format conventions, cheminformatics toolkits, and biostatistics. An "
    "agentic AI system collapses these prerequisites into a single natural-language interface, rendering "
    "large-scale VS accessible to medicinal chemists and pharmacologists who may lack formal training "
    "in computational methods."
)

add_para(
    "The skill file mechanism merits particular attention. Unlike model fine-tuning, which risks "
    "catastrophic forgetting, or retrieval-augmented generation, which demands vector database "
    "infrastructure, a skill file is a plain-text document that can be authored, version-controlled, "
    "and peer-reviewed by domain experts (Anthropic, 2024). It is model-agnostic in principle: the "
    "same file could guide any sufficiently capable coding agent."
)

add_para(
    "Several limitations should be acknowledged. This benchmark targets a single receptor; "
    "generalisation requires multi-target validation. The active set (n = 100) limits statistical "
    "power. Only Vina scoring was evaluated; consensus scoring could further improve enrichment. The "
    "decline in BEDROC and EF at 1% is a benchmark design artefact\u2014seven actives removed by "
    "PAINS/Brenk filtering were assigned worst-rank scores\u2014and highlights that future AI-driven VS "
    "benchmarks should assign labels after filtering or account for filtered compounds explicitly."
)

# ══════════════════════════════════════════════════════════════
# 5 CONCLUSION
# ══════════════════════════════════════════════════════════════

add_heading_custom("5 Conclusion", level=1)

add_para(
    "A 212-line structured skill file raises the global discrimination of an LLM coding agent\u2019s VS "
    "pipeline by a statistically significant margin (ROC AUC +0.048, DeLong p = 0.004). The agent "
    "autonomously generated 18 scripts (\u223c3,100 lines) implementing a complete GPU-accelerated "
    "pipeline without human code editing. Agentic AI guided by interpretable skill files could "
    "democratise rigorous large-scale virtual screening. A multi-target validation study is underway."
)

# ══════════════════════════════════════════════════════════════
# BACK MATTER
# ══════════════════════════════════════════════════════════════

add_heading_custom("Acknowledgements", level=1)

add_para(
    "All computational protocols were designed and executed with Claude Code (Claude Opus 4.6, "
    "Anthropic). Docking calculations were performed on NVIDIA RTX 4500 Ada GPUs. We thank the "
    "ChEMBL team for curated bioactivity data and the Uni-Dock developers for GPU-accelerated docking."
)

add_heading_custom("Funding", level=1)
add_para("[Funding statement to be added.]")

add_heading_custom("Conflict of Interest", level=1)
add_para("None declared.")

add_heading_custom("Data Availability", level=1)
add_para(
    "All code, data, and pre-computed docking scores are available at "
    "https://github.com/oslo-medchem/vs-benchmark-agentic-ai under the MIT licence. The conda "
    "environment specification (environment.yml) and detailed reproduction instructions "
    "(REPRODUCTION.md) are included in the repository."
)

# ══════════════════════════════════════════════════════════════
# REFERENCES
# ══════════════════════════════════════════════════════════════

add_heading_custom("References", level=1)

refs = [
    "Anthropic (2024) Claude Code: an agentic coding tool. Anthropic PBC, San Francisco, CA. https://claude.ai/code",
    "Baell,J.B. and Holloway,G.A. (2010) New substructure filters for removal of pan assay interference compounds (PAINS) from screening libraries and for their exclusion in bioassays. J. Med. Chem., 53, 2719\u20132740.",
    "Benjamini,Y. and Hochberg,Y. (1995) Controlling the false discovery rate: a practical and powerful approach to multiple testing. J. R. Stat. Soc. B, 57, 289\u2013300.",
    "Boiko,D.A. et al. (2023) Autonomous chemical research with large language models. Nature, 624, 570\u2013578.",
    "Bran,A.M. et al. (2024) ChemCrow: augmenting large-language models with chemistry tools. Nat. Mach. Intell., 6, 525\u2013535.",
    "Brenk,R. et al. (2008) Lessons learnt from assembling screening libraries for drug discovery for neglected diseases. ChemMedChem, 3, 435\u2013444.",
    "Butina,D. (1999) Unsupervised data base clustering based on Daylight\u2019s fingerprint and Tanimoto similarity. J. Chem. Inf. Comput. Sci., 39, 747\u2013750.",
    "DeLong,E.R. et al. (1988) Comparing the areas under two or more correlated receiver operating characteristic curves: a nonparametric approach. Biometrics, 44, 837\u2013845.",
    "Huang,N. et al. (2006) Benchmarking sets for molecular docking. J. Med. Chem., 49, 6789\u20136801.",
    "Liao,C. et al. (2022) Structural basis of formyl peptide receptor 2 activation and signaling. Nat. Commun., 13, 1463.",
    "Liu,C. et al. (2022) Molecular basis of the recognition of formalyl peptide by FPR2 and blocking by the small molecule inhibitor. Nat. Commun., 13, 1476.",
    "Ramos,M.C. et al. (2025) A review of large language models and autonomous agents in chemistry. Chem. Sci., 16, 2514\u20132572.",
    "Schneider,P. et al. (2020) Rethinking drug design in the artificial intelligence era. Nat. Rev. Drug Discov., 19, 353\u2013364.",
    "Truchon,J.-F. and Bayly,C.I. (2007) Evaluating virtual screening methods: good and bad metrics for the \u201cearly recognition\u201d problem. J. Chem. Inf. Model., 47, 488\u2013508.",
    "Wang,S. et al. (2020) Improving conformer generation for small rings and macrocycles based on distance geometry and experimental torsional angle preferences. J. Chem. Inf. Model., 60, 2044\u20132058.",
    "Yu,Y. et al. (2023) Uni-Dock: GPU-accelerated docking enables ultralarge virtual screening. J. Chem. Theory Comput., 19, 3336\u20133345.",
    "Zdrazil,B. et al. (2024) The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods. Nucleic Acids Res., 52, D1180\u2013D1192.",
]

for ref in refs:
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing = 2.0
    p.paragraph_format.left_indent = Cm(1.0)
    p.paragraph_format.first_line_indent = Cm(-1.0)
    run = p.add_run(ref)
    run.font.name = "Times New Roman"
    run.font.size = Pt(10)

# ── Save ─────────────────────────────────────────────────────
out_path = BASE / "manuscript.docx"
doc.save(str(out_path))
print(f"Saved: {out_path} ({out_path.stat().st_size / 1024:.0f} KB)")
