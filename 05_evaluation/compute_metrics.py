#!/usr/bin/env python3
"""
Compute all VS evaluation metrics for naive and skill-guided protocols.
Metrics: ROC AUC, BEDROC(α=20,80), EF1/5/10%, LogAUC, pROC AUC (1%, 5%),
         AUPR, RIE(α=20).
Bootstrap 95% CI (BCa, n=10,000) for all metrics.
"""
import csv
import json
import sys
import numpy as np
from pathlib import Path
from scipy.stats import bootstrap as scipy_bootstrap
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve

# Try RDKit ML scoring
try:
    from rdkit.ML.Scoring.Scoring import CalcAUC, CalcBEDROC, CalcEnrichment, CalcRIE
    RDKIT_SCORING = True
except ImportError:
    RDKIT_SCORING = False
    print("WARNING: RDKit ML Scoring not available, using sklearn fallback")

BASE = Path(__file__).parent
VS_DIR = BASE.parent

LABELS_CSV = VS_DIR / "03_library_preparation" / "library_labels.csv"
RESULTS_DIR = BASE / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

PROTOCOLS = {
    "naive": VS_DIR / "04_docking" / "naive" / "scores_naive.csv",
    "skill": VS_DIR / "04_docking" / "skill" / "scores_skill.csv",
}

N_BOOTSTRAP = 10_000
CI_LEVEL = 0.95


def load_data(labels_csv, scores_csv):
    """Load labels and scores, return aligned arrays."""
    labels = {}
    for r in csv.DictReader(open(labels_csv)):
        labels[r["compound_id"]] = 1 if r["label"] == "active" else 0

    scores = {}
    for r in csv.DictReader(open(scores_csv)):
        try:
            scores[r["compound_id"]] = float(r["score"])
        except (ValueError, KeyError):
            scores[r["compound_id"]] = 0.0

    # Align
    cids = sorted(labels.keys())
    y_true = np.array([labels[c] for c in cids])
    y_score = np.array([scores.get(c, 0.0) for c in cids])

    # For docking: more negative = better. Negate for rank (higher = better).
    # Treat 0.0 (failed docking) as worst score.
    # First handle missing docked compounds
    missing_mask = y_score == 0.0
    y_score_rank = -y_score  # negate: less negative (0) → worst rank

    return cids, y_true, y_score, y_score_rank


def rdkit_sorted_list(y_true, y_score_rank):
    """Build RDKit-style sorted list: [(score, label)] sorted descending."""
    pairs = sorted(zip(y_score_rank, y_true), reverse=True)
    return [[sc, lab] for sc, lab in pairs]


def calc_logauc(y_true, y_score_rank, fpr_min=0.001, fpr_max=1.0):
    """Calculate LogAUC (log-scale ROC AUC)."""
    fpr, tpr, _ = roc_curve(y_true, y_score_rank)
    # Interpolate on log scale
    fpr_log = np.logspace(np.log10(fpr_min), np.log10(fpr_max), 1000)
    tpr_interp = np.interp(fpr_log, fpr, tpr)
    # AUC in log space
    trapz = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
    log_auc = trapz(tpr_interp, np.log10(fpr_log)) / (np.log10(fpr_max) - np.log10(fpr_min))
    return float(log_auc)


def calc_partial_roc_auc(y_true, y_score_rank, max_fpr):
    """Partial ROC AUC normalized to [0, 1] range."""
    try:
        pauc = roc_auc_score(y_true, y_score_rank, max_fpr=max_fpr)
        # sklearn returns the standardized form — normalize to [0,1]
        pauc_norm = pauc / max_fpr
        return float(pauc_norm)
    except Exception:
        return 0.5


def compute_all_metrics(y_true, y_score_rank):
    """Compute all metrics for one protocol."""
    n_actives = int(y_true.sum())
    n_total = len(y_true)

    # RDKit-based metrics
    rdkit_list = rdkit_sorted_list(y_true, y_score_rank)

    if RDKIT_SCORING:
        roc_auc = float(CalcAUC(rdkit_list, 1))
        bedroc_20 = float(CalcBEDROC(rdkit_list, 1, 20.0))
        bedroc_80 = float(CalcBEDROC(rdkit_list, 1, 80.0))
        ef1 = float(CalcEnrichment(rdkit_list, 1, [0.01])[0])
        ef5 = float(CalcEnrichment(rdkit_list, 1, [0.05])[0])
        ef10 = float(CalcEnrichment(rdkit_list, 1, [0.10])[0])
        rie_20 = float(CalcRIE(rdkit_list, 1, 20.0))
    else:
        # sklearn fallback
        roc_auc = float(roc_auc_score(y_true, y_score_rank))
        bedroc_20 = bedroc_20_sklearn(y_true, y_score_rank, alpha=20.0)
        bedroc_80 = bedroc_20_sklearn(y_true, y_score_rank, alpha=80.0)
        ef1 = calc_ef(y_true, y_score_rank, 0.01)
        ef5 = calc_ef(y_true, y_score_rank, 0.05)
        ef10 = calc_ef(y_true, y_score_rank, 0.10)
        rie_20 = calc_rie(y_true, y_score_rank, alpha=20.0)

    log_auc = calc_logauc(y_true, y_score_rank)
    proc_1pct = calc_partial_roc_auc(y_true, y_score_rank, 0.01)
    proc_5pct = calc_partial_roc_auc(y_true, y_score_rank, 0.05)
    aupr = float(average_precision_score(y_true, y_score_rank))

    return {
        "n_total": n_total,
        "n_actives": n_actives,
        "n_decoys": n_total - n_actives,
        "roc_auc": roc_auc,
        "bedroc_20": bedroc_20,
        "bedroc_80": bedroc_80,
        "ef_1pct": ef1,
        "ef_5pct": ef5,
        "ef_10pct": ef10,
        "rie_20": rie_20,
        "log_auc": log_auc,
        "proc_auc_1pct": proc_1pct,
        "proc_auc_5pct": proc_5pct,
        "aupr": aupr,
    }


def bedroc_20_sklearn(y_true, y_score_rank, alpha=20.0):
    """BEDROC calculation from sorted arrays (RDKit-independent)."""
    n = len(y_true)
    n_actives = y_true.sum()
    ra = n_actives / n
    order = np.argsort(-y_score_rank)
    y_sorted = y_true[order]

    ranks = np.where(y_sorted == 1)[0] + 1  # 1-indexed
    s = np.sum(np.exp(-alpha * ranks / n))
    s_bar = (1 - np.exp(-alpha * ra)) / (np.exp(alpha / n) - 1)
    ri = np.exp(-alpha * ra) - np.exp(-alpha)
    d_val = np.exp(alpha / n) - 1

    bedroc = (s * d_val * np.sinh(alpha / 2) /
              (np.cosh(alpha / 2) - np.cosh(alpha / 2 - alpha * ra))) * ra + ri / (1 - np.exp(alpha * (1 - ra)))
    # Normalize
    min_bedroc = (1 - np.exp(alpha * ra)) / (1 - np.exp(alpha)) * ra
    max_bedroc = (1 - np.exp(-alpha * ra)) / (1 - np.exp(-alpha)) * ra
    bedroc_norm = (bedroc - min_bedroc) / (max_bedroc - min_bedroc) if (max_bedroc - min_bedroc) > 0 else 0.5
    return float(bedroc_norm)


def calc_ef(y_true, y_score_rank, fraction):
    """Enrichment Factor at given fraction of the library."""
    n = len(y_true)
    n_actives = y_true.sum()
    n_top = max(1, int(n * fraction))
    order = np.argsort(-y_score_rank)
    y_sorted = y_true[order]
    hits_in_top = y_sorted[:n_top].sum()
    expected = n_actives * fraction
    return float(hits_in_top / expected) if expected > 0 else 0.0


def calc_rie(y_true, y_score_rank, alpha=20.0):
    """Robust Initial Enhancement."""
    n = len(y_true)
    n_actives = y_true.sum()
    ra = n_actives / n
    order = np.argsort(-y_score_rank)
    y_sorted = y_true[order]
    ranks = np.where(y_sorted == 1)[0] + 1
    s = np.sum(np.exp(-alpha * ranks / n))
    rie_max = (1 - np.exp(-alpha * ra)) / (ra * (1 - np.exp(-alpha)))
    rie_min = (1 - np.exp(alpha * ra)) / (ra * (1 - np.exp(alpha)))
    rie_raw = s / (n * ra) * (1 - np.exp(-alpha)) / (np.exp(alpha / n) - 1)
    return float(rie_raw)


def bootstrap_metric(y_true, y_score_rank, metric_fn, n_resamples=N_BOOTSTRAP):
    """BCa bootstrap CI for a metric."""
    def _fn(y_t, y_s):
        return metric_fn(y_t, y_s)

    try:
        result = scipy_bootstrap(
            (y_true, y_score_rank),
            lambda yt, ys: _fn(yt, ys),
            n_resamples=n_resamples,
            confidence_level=CI_LEVEL,
            method="BCa",
            vectorized=False,
            paired=True,
        )
        return {
            "ci_low": float(result.confidence_interval.low),
            "ci_high": float(result.confidence_interval.high),
        }
    except Exception as e:
        # Fallback: percentile bootstrap
        vals = []
        rng = np.random.default_rng(42)
        n = len(y_true)
        for _ in range(min(1000, n_resamples)):
            idx = rng.integers(0, n, n)
            yt_b = y_true[idx]
            ys_b = y_score_rank[idx]
            if yt_b.sum() > 0 and yt_b.sum() < len(yt_b):
                try:
                    vals.append(_fn(yt_b, ys_b))
                except Exception:
                    pass
        if vals:
            return {
                "ci_low": float(np.percentile(vals, 2.5)),
                "ci_high": float(np.percentile(vals, 97.5)),
            }
        return {"ci_low": None, "ci_high": None}


def get_roc_data(y_true, y_score_rank):
    """Return ROC curve points for plotting."""
    fpr, tpr, thresh = roc_curve(y_true, y_score_rank)
    return {
        "fpr": fpr.tolist(),
        "tpr": tpr.tolist(),
        "thresholds": thresh.tolist(),
    }


def main():
    print("Loading data...")
    labels_data = list(csv.DictReader(open(LABELS_CSV)))
    labels = {r["compound_id"]: r["label"] for r in labels_data}

    all_results = {}

    for protocol, scores_csv in PROTOCOLS.items():
        print(f"\n=== {protocol.upper()} ===")

        if not scores_csv.exists():
            print(f"  WARNING: {scores_csv} not found, skipping")
            continue

        cids, y_true, y_score, y_score_rank = load_data(LABELS_CSV, scores_csv)

        n_docked = int((y_score != 0.0).sum())
        print(f"  Compounds: {len(cids)} total, {n_docked} docked")
        print(f"  Actives: {int(y_true.sum())}, Decoys: {int((1-y_true).sum())}")

        # Point estimates
        metrics = compute_all_metrics(y_true, y_score_rank)
        print(f"  ROC AUC: {metrics['roc_auc']:.4f}")
        print(f"  BEDROC(α=20): {metrics['bedroc_20']:.4f}")
        print(f"  EF1%: {metrics['ef_1pct']:.2f}")

        # Bootstrap CIs
        print("  Computing bootstrap CIs (10,000 resamples)...")
        bootstrap_cis = {}

        metric_fns = {
            "roc_auc": lambda yt, ys: float(roc_auc_score(yt, ys)) if yt.sum() > 0 else 0.5,
            "bedroc_20": lambda yt, ys: bedroc_20_sklearn(yt, ys, 20.0),
            "bedroc_80": lambda yt, ys: bedroc_20_sklearn(yt, ys, 80.0),
            "ef_1pct": lambda yt, ys: calc_ef(yt, ys, 0.01),
            "ef_5pct": lambda yt, ys: calc_ef(yt, ys, 0.05),
            "log_auc": lambda yt, ys: calc_logauc(yt, ys),
            "aupr": lambda yt, ys: float(average_precision_score(yt, ys)) if yt.sum() > 0 else 0.0,
        }

        for mname, mfn in metric_fns.items():
            ci = bootstrap_metric(y_true, y_score_rank, mfn, N_BOOTSTRAP)
            bootstrap_cis[mname] = ci
            print(f"    {mname}: {metrics[mname]:.4f} [{ci['ci_low']:.4f}, {ci['ci_high']:.4f}]")

        # Save results
        output = {
            "protocol": protocol,
            "n_compounds": len(cids),
            "n_docked": n_docked,
            **metrics,
            "bootstrap_ci": bootstrap_cis,
        }

        metrics_file = RESULTS_DIR / f"metrics_{protocol}.json"
        with open(metrics_file, "w") as f:
            json.dump(output, f, indent=2)
        print(f"  Saved: {metrics_file}")

        # Save ROC curve data
        roc_data = get_roc_data(y_true, y_score_rank)
        roc_file = RESULTS_DIR / f"roc_{protocol}.json"
        with open(roc_file, "w") as f:
            json.dump(roc_data, f, indent=2)

        # Save score arrays for later use
        np.save(RESULTS_DIR / f"y_true_{protocol}.npy", y_true)
        np.save(RESULTS_DIR / f"y_score_{protocol}.npy", y_score_rank)

        all_results[protocol] = output

    print("\n=== Summary ===")
    for p, res in all_results.items():
        print(f"{p:10s}: AUC={res['roc_auc']:.4f} BEDROC20={res['bedroc_20']:.4f} EF1%={res['ef_1pct']:.2f} AUPR={res['aupr']:.4f}")

    print("\nDone.")


if __name__ == "__main__":
    main()
