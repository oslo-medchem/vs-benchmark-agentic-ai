#!/usr/bin/env python3
"""
Statistical analysis comparing naive vs skill-guided VS protocols.
Tests: DeLong's AUC comparison, bootstrap permutation, Mann-Whitney U,
       Benjamini-Hochberg FDR, Net Reclassification Improvement.
Inputs: results/{y_true,y_score}_{naive,skill}.npy, metrics_{naive,skill}.json
Outputs: results/comparison_stats.json, results/bootstrap_stats.json
"""
import json
import numpy as np
from pathlib import Path
from scipy import stats
from scipy.stats import mannwhitneyu, bootstrap as scipy_bootstrap
from sklearn.metrics import roc_auc_score

BASE = Path(__file__).parent
RESULTS = BASE / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

N_BOOTSTRAP = 10_000
N_PERM = 10_000
RNG = np.random.default_rng(42)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_data(protocol):
    yt = np.load(RESULTS / f"y_true_{protocol}.npy")
    ys = np.load(RESULTS / f"y_score_{protocol}.npy")
    return yt, ys


def load_metrics(protocol):
    with open(RESULTS / f"metrics_{protocol}.json") as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# DeLong's test (compare two correlated ROC curves on the same sample)
# Reference: DeLong et al. Biometrika 1988
# ---------------------------------------------------------------------------

def _delong_placement(y_true, y_score):
    """Compute placement values for DeLong's test."""
    n_pos = int(y_true.sum())
    n_neg = int((1 - y_true).sum())
    pos_scores = y_score[y_true == 1]
    neg_scores = y_score[y_true == 0]
    # V10: for each positive, fraction of negatives it outranks
    v10 = np.array([np.mean(ps > neg_scores) + 0.5 * np.mean(ps == neg_scores) for ps in pos_scores])
    # V01: for each negative, fraction of positives it is outranked by
    v01 = np.array([np.mean(ns < pos_scores) + 0.5 * np.mean(ns == pos_scores) for ns in neg_scores])
    return v10, v01, n_pos, n_neg


def delong_test(y_true, y_score_a, y_score_b):
    """
    DeLong's test for comparing two correlated ROC AUC estimates.
    Returns: z_stat, p_value, delta_auc, ci_low, ci_high (95% CI on delta AUC)
    """
    v10_a, v01_a, n_pos, n_neg = _delong_placement(y_true, y_score_a)
    v10_b, v01_b, _, _ = _delong_placement(y_true, y_score_b)

    auc_a = v10_a.mean()
    auc_b = v10_b.mean()
    delta = auc_b - auc_a

    # Covariance matrix of the two AUC estimators
    s_10_aa = np.var(v10_a, ddof=1) / n_pos
    s_10_bb = np.var(v10_b, ddof=1) / n_pos
    s_10_ab = np.cov(v10_a, v10_b, ddof=1)[0, 1] / n_pos

    s_01_aa = np.var(v01_a, ddof=1) / n_neg
    s_01_bb = np.var(v01_b, ddof=1) / n_neg
    s_01_ab = np.cov(v01_a, v01_b, ddof=1)[0, 1] / n_neg

    var_delta = s_10_aa + s_10_bb - 2 * s_10_ab + s_01_aa + s_01_bb - 2 * s_01_ab
    if var_delta <= 0:
        return 0.0, 1.0, float(delta), None, None

    se = np.sqrt(var_delta)
    z = delta / se
    p = 2 * (1 - stats.norm.cdf(abs(z)))

    ci_low = delta - 1.96 * se
    ci_high = delta + 1.96 * se

    return float(z), float(p), float(delta), float(ci_low), float(ci_high)


# ---------------------------------------------------------------------------
# Bootstrap permutation test for a generic metric difference
# ---------------------------------------------------------------------------

def _bedroc_sklearn(y_true, y_score_rank, alpha=20.0):
    n = len(y_true)
    n_actives = y_true.sum()
    if n_actives == 0 or n_actives == n:
        return 0.5
    ra = n_actives / n
    order = np.argsort(-y_score_rank)
    y_sorted = y_true[order]
    ranks = np.where(y_sorted == 1)[0] + 1
    s = np.sum(np.exp(-alpha * ranks / n))
    d_val = np.exp(alpha / n) - 1
    ri = np.exp(-alpha * ra) - np.exp(-alpha)
    bedroc = (s * d_val * np.sinh(alpha / 2) /
              (np.cosh(alpha / 2) - np.cosh(alpha / 2 - alpha * ra))) * ra + ri / (1 - np.exp(alpha * (1 - ra)))
    min_b = (1 - np.exp(alpha * ra)) / (1 - np.exp(alpha)) * ra
    max_b = (1 - np.exp(-alpha * ra)) / (1 - np.exp(-alpha)) * ra
    return float((bedroc - min_b) / (max_b - min_b)) if (max_b - min_b) > 0 else 0.5


def _calc_ef(y_true, y_score_rank, fraction):
    n = len(y_true)
    n_actives = y_true.sum()
    n_top = max(1, int(n * fraction))
    order = np.argsort(-y_score_rank)
    hits = y_true[order[:n_top]].sum()
    expected = n_actives * fraction
    return float(hits / expected) if expected > 0 else 0.0


METRIC_FNS = {
    "roc_auc": lambda yt, ys: float(roc_auc_score(yt, ys)) if 0 < yt.sum() < len(yt) else 0.5,
    "bedroc_20": lambda yt, ys: _bedroc_sklearn(yt, ys, 20.0),
    "ef_1pct": lambda yt, ys: _calc_ef(yt, ys, 0.01),
    "ef_5pct": lambda yt, ys: _calc_ef(yt, ys, 0.05),
}


def permutation_test_delta(y_true_a, y_score_a, y_true_b, y_score_b, metric_fn, n_perm=N_PERM):
    """
    Two-sided permutation test for H0: metric(A) == metric(B).
    Permutes the protocol labels (not compound labels).
    Returns: observed_delta, p_value
    """
    obs_a = metric_fn(y_true_a, y_score_a)
    obs_b = metric_fn(y_true_b, y_score_b)
    observed_delta = obs_b - obs_a

    # Permutation: randomly swap which score vector goes to "a" vs "b"
    all_scores_a = np.concatenate([y_score_a, y_score_b])
    n = len(y_score_a)
    null_deltas = []
    for _ in range(n_perm):
        perm = RNG.permutation(len(all_scores_a))
        s_a = all_scores_a[perm[:n]]
        s_b = all_scores_a[perm[n:]]
        # y_true stays the same for each (same compounds)
        d = metric_fn(y_true_b, s_b) - metric_fn(y_true_a, s_a)
        null_deltas.append(d)

    null_deltas = np.array(null_deltas)
    p_val = float(np.mean(np.abs(null_deltas) >= np.abs(observed_delta)))

    return float(observed_delta), p_val


# ---------------------------------------------------------------------------
# Mann-Whitney U test (actives vs decoys score distributions)
# ---------------------------------------------------------------------------

def mannwhitney_test(y_true, y_score_raw):
    """
    Compare raw docking scores of actives vs decoys.
    Lower (more negative) score = better binding.
    H1: actives have lower (better) docking scores than decoys.
    Returns: U, p_value, r (rank-biserial correlation effect size)
    """
    actives_scores = y_score_raw[y_true == 1]
    decoys_scores = y_score_raw[y_true == 0]

    # Exclude failed dockings (score == 0.0)
    actives_scores = actives_scores[actives_scores != 0.0]
    decoys_scores = decoys_scores[decoys_scores != 0.0]

    if len(actives_scores) == 0 or len(decoys_scores) == 0:
        return None, None, None

    # alternative='less': test that actives scores < decoys scores
    result = mannwhitneyu(actives_scores, decoys_scores, alternative='less')
    U = float(result.statistic)
    p = float(result.pvalue)

    # Rank-biserial correlation as effect size: r = 1 - 2U/(n1*n2)
    n1, n2 = len(actives_scores), len(decoys_scores)
    r = 1 - 2 * U / (n1 * n2)

    return U, p, float(r)


# ---------------------------------------------------------------------------
# Benjamini-Hochberg FDR correction
# ---------------------------------------------------------------------------

def bh_correction(p_values):
    """Benjamini-Hochberg FDR correction. Returns corrected p-values."""
    n = len(p_values)
    sorted_idx = np.argsort(p_values)
    sorted_p = np.array(p_values)[sorted_idx]
    corrected = np.zeros(n)
    for i in range(n - 1, -1, -1):
        if i == n - 1:
            corrected[i] = sorted_p[i]
        else:
            corrected[i] = min(corrected[i + 1], sorted_p[i] * n / (i + 1))
    corrected = np.minimum(corrected, 1.0)
    # Map back to original order
    result = np.zeros(n)
    result[sorted_idx] = corrected
    return result.tolist()


# ---------------------------------------------------------------------------
# Net Reclassification Improvement at EF1% threshold
# ---------------------------------------------------------------------------

def calc_nri_at_ef_threshold(y_true_a, y_score_a, y_true_b, y_score_b, fraction=0.01):
    """
    NRI at EF1% threshold.
    Compounds in top 1% of skill that were NOT in top 1% of naive (up-events for actives).
    Returns: NRI, NRI_actives, NRI_decoys, 95% CI via bootstrap.
    """
    n = len(y_true_a)
    n_top = max(1, int(n * fraction))

    order_a = np.argsort(-y_score_a)
    order_b = np.argsort(-y_score_b)

    in_top_a = np.zeros(n, dtype=bool)
    in_top_b = np.zeros(n, dtype=bool)
    in_top_a[order_a[:n_top]] = True
    in_top_b[order_b[:n_top]] = True

    # Reclassification among actives: moved up - moved down
    act = y_true_a == 1
    up_act = np.sum(in_top_b[act] & ~in_top_a[act])
    down_act = np.sum(~in_top_b[act] & in_top_a[act])
    n_act = act.sum()

    # Reclassification among decoys (up = wrong, down = correct)
    dec = y_true_a == 0
    up_dec = np.sum(in_top_b[dec] & ~in_top_a[dec])
    down_dec = np.sum(~in_top_b[dec] & in_top_a[dec])
    n_dec = dec.sum()

    nri_act = float((up_act - down_act) / n_act) if n_act > 0 else 0.0
    nri_dec = float((down_dec - up_dec) / n_dec) if n_dec > 0 else 0.0
    nri = nri_act + nri_dec

    return {
        "nri": float(nri),
        "nri_actives": float(nri_act),
        "nri_decoys": float(nri_dec),
        "n_actives_up": int(up_act),
        "n_actives_down": int(down_act),
    }


# ---------------------------------------------------------------------------
# Bootstrap CI for delta metric
# ---------------------------------------------------------------------------

def bootstrap_delta_ci(y_true_a, y_score_a, y_true_b, y_score_b, metric_fn,
                        n_resamples=N_BOOTSTRAP):
    """Bootstrap 95% CI for (metric_b - metric_a)."""
    n = len(y_true_a)
    deltas = []
    for _ in range(n_resamples):
        idx = RNG.integers(0, n, n)
        yt_a = y_true_a[idx]
        ys_a = y_score_a[idx]
        yt_b = y_true_b[idx]
        ys_b = y_score_b[idx]
        if yt_a.sum() == 0 or yt_b.sum() == 0:
            continue
        try:
            d = metric_fn(yt_b, ys_b) - metric_fn(yt_a, ys_a)
            deltas.append(d)
        except Exception:
            pass

    if len(deltas) < 100:
        return {"ci_low": None, "ci_high": None, "n_valid": len(deltas)}

    deltas = np.array(deltas)
    return {
        "ci_low": float(np.percentile(deltas, 2.5)),
        "ci_high": float(np.percentile(deltas, 97.5)),
        "n_valid": len(deltas),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("Loading data...")
    yt_n, ys_n = load_data("naive")
    yt_s, ys_s = load_data("skill")
    m_n = load_metrics("naive")
    m_s = load_metrics("skill")

    # Raw scores (before negation) — needed for Mann-Whitney
    # y_score saved in compute_metrics is already negated (y_score_rank = -y_score_raw)
    # So raw score = -y_score_rank (more negative raw = better)
    ys_n_raw = -ys_n
    ys_s_raw = -ys_s

    comparison = {}
    bootstrap_stats = {}

    # -----------------------------------------------------------------------
    # 1. DeLong's test — AUC comparison
    # -----------------------------------------------------------------------
    print("\n=== DeLong's Test (AUC naive vs skill) ===")
    z, p_delong, delta_auc, ci_low, ci_high = delong_test(yt_n, ys_n, ys_s)
    print(f"  AUC naive:  {m_n['roc_auc']:.4f}")
    print(f"  AUC skill:  {m_s['roc_auc']:.4f}")
    print(f"  ΔAUC:       {delta_auc:.4f} [{ci_low:.4f}, {ci_high:.4f}]")
    print(f"  Z={z:.3f}, p={p_delong:.4f}")
    comparison["delong_auc"] = {
        "auc_naive": m_n["roc_auc"],
        "auc_skill": m_s["roc_auc"],
        "delta_auc": delta_auc,
        "ci_95_low": ci_low,
        "ci_95_high": ci_high,
        "z_stat": z,
        "p_value": p_delong,
    }

    # -----------------------------------------------------------------------
    # 2. Permutation tests for BEDROC and EF differences
    # -----------------------------------------------------------------------
    print("\n=== Permutation Tests (10,000 permutations) ===")
    perm_p_values = {}
    for mname, fn in METRIC_FNS.items():
        print(f"  {mname}...", end=" ", flush=True)
        obs_delta, p_perm = permutation_test_delta(yt_n, ys_n, yt_s, ys_s, fn)
        print(f"Δ={obs_delta:.4f}, p={p_perm:.4f}")
        perm_p_values[mname] = p_perm
        comparison[f"permutation_{mname}"] = {
            "observed_delta": obs_delta,
            "p_value_permutation": p_perm,
            "metric_naive": m_n.get(mname),
            "metric_skill": m_s.get(mname),
        }

    # -----------------------------------------------------------------------
    # 3. Mann-Whitney U test (actives vs decoys within each protocol)
    # -----------------------------------------------------------------------
    print("\n=== Mann-Whitney U (actives vs decoys score distributions) ===")
    for protocol, yt, ys_raw in [("naive", yt_n, ys_n_raw), ("skill", yt_s, ys_s_raw)]:
        U, p_mw, r = mannwhitney_test(yt, ys_raw)
        if U is not None:
            print(f"  {protocol}: U={U:.0f}, p={p_mw:.2e}, r={r:.4f}")
        else:
            print(f"  {protocol}: insufficient data")
        comparison[f"mannwhitney_{protocol}"] = {
            "U_statistic": U,
            "p_value": p_mw,
            "rank_biserial_r": r,
            "interpretation": "actives ranked better than decoys (p < 0.05 = significant)" if p_mw and p_mw < 0.05 else "not significant",
        }

    # -----------------------------------------------------------------------
    # 4. Bootstrap CIs for metric deltas
    # -----------------------------------------------------------------------
    print("\n=== Bootstrap CIs for Metric Deltas ===")
    for mname, fn in METRIC_FNS.items():
        print(f"  {mname}...", end=" ", flush=True)
        ci = bootstrap_delta_ci(yt_n, ys_n, yt_s, ys_s, fn, N_BOOTSTRAP)
        print(f"Δ CI=[{ci['ci_low']:.4f}, {ci['ci_high']:.4f}]" if ci["ci_low"] is not None else "failed")
        bootstrap_stats[mname] = ci

    # -----------------------------------------------------------------------
    # 5. Benjamini-Hochberg FDR correction
    # -----------------------------------------------------------------------
    print("\n=== Benjamini-Hochberg FDR Correction ===")
    all_p_raw = {}
    all_p_raw["delong_auc"] = comparison["delong_auc"]["p_value"]
    for mname in METRIC_FNS:
        all_p_raw[f"perm_{mname}"] = comparison[f"permutation_{mname}"]["p_value_permutation"]
    all_p_raw["mw_naive"] = comparison["mannwhitney_naive"]["p_value"] or 1.0
    all_p_raw["mw_skill"] = comparison["mannwhitney_skill"]["p_value"] or 1.0

    keys = list(all_p_raw.keys())
    raw_vals = [all_p_raw[k] for k in keys]
    corrected_vals = bh_correction(raw_vals)

    bh_results = {}
    for k, p_raw, p_corr in zip(keys, raw_vals, corrected_vals):
        bh_results[k] = {"p_raw": p_raw, "p_bh": p_corr, "significant_bh": p_corr < 0.05}
        sig = "* " if p_corr < 0.05 else "  "
        print(f"  {sig}{k}: p_raw={p_raw:.4f} → p_BH={p_corr:.4f}")

    comparison["bh_correction"] = bh_results

    # -----------------------------------------------------------------------
    # 6. Net Reclassification Improvement at EF1%
    # -----------------------------------------------------------------------
    print("\n=== NRI at EF1% threshold ===")
    nri_result = calc_nri_at_ef_threshold(yt_n, ys_n, yt_s, ys_s, fraction=0.01)
    print(f"  NRI total: {nri_result['nri']:.4f}")
    print(f"  NRI actives: {nri_result['nri_actives']:.4f}  (actives up: {nri_result['n_actives_up']}, down: {nri_result['n_actives_down']})")
    print(f"  NRI decoys: {nri_result['nri_decoys']:.4f}")
    comparison["nri_ef1pct"] = nri_result

    # -----------------------------------------------------------------------
    # Save outputs
    # -----------------------------------------------------------------------
    with open(RESULTS / "comparison_stats.json", "w") as f:
        json.dump(comparison, f, indent=2, default=lambda x: None if np.isnan(x) else float(x))

    with open(RESULTS / "bootstrap_stats.json", "w") as f:
        json.dump(bootstrap_stats, f, indent=2)

    print(f"\nSaved: {RESULTS}/comparison_stats.json")
    print(f"Saved: {RESULTS}/bootstrap_stats.json")
    print("\nDone.")


if __name__ == "__main__":
    main()
