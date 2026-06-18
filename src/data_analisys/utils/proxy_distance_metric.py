"""
Distance Evaluation Metrics
============================
Implements the compendium-wide batch correction evaluation metrics described in:

    sim(L1, L2)     — multi-axial metadata similarity score
    G_d_inter       — global inter-study weighted Euclidean distance
    G_d_intra       — global intra-study weighted Euclidean distance
    Ratio_global    — G_d_intra / G_d_inter (integration headline metric)

These metrics operate in a PCA-reduced expression space and use a weighted
label-agreement similarity to up-weight distances between biologically
comparable samples.
"""

from __future__ import annotations

import warnings
from typing import Dict, Optional

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances

from src.data_analisys.utils.cluster_exploration_utils_final import get_gsm_id  # noqa: E402
from src.constants import RNA_USED
# ---------------------------------------------------------------------------
# Label axis weights  (edit to reflect biological importance)
# ---------------------------------------------------------------------------
# DEFAULT_AXIS_WEIGHTS: Dict[str, float] = {
#     "tissue":              3.0,
#     "developmental_stage": 2.0,
#     "treatment":           1.5,
#     "treatment_intensity": 1.0,
#     "ecotype":             1.0,
#     "medium":              0.5,
#     "modification":        0.5,
# }

DEFAULT_AXIS_WEIGHTS = {
    "tissue":              4.0,
    "developmental_stage": 3.0,
    "modification":        2.5,
    "ecotype":             2.0,
    "treatment":           2.0,
    "treatment_intensity": 1.0,
    "medium":              0.5,
}

# Values treated as "unknown / not annotated" — excluded from similarity
_UNKNOWN_VALUES = {"unspecified", "unknown", "na", "nan", "none", ""}


# ---------------------------------------------------------------------------
# Core similarity function
# ---------------------------------------------------------------------------

def compute_sim(
    labels_s: Dict[str, Optional[str]],
    labels_i: Dict[str, Optional[str]],
    weights: Dict[str, float],
) -> float:
    """
    Compute the multi-axial metadata similarity score sim(L_s, L_i).

    For each label axis i with weight w_i, both samples must have a known
    (non-unspecified) annotation for the axis to contribute.  The score is
    the fraction of valid shared axes on which the two labels agree,
    weighted by w_i.

        sim = sum_i  w_i * 1[ L1^i == L2^i  AND  v1^i  AND  v2^i ]
              --------------------------------------------------------
              sum_i  w_i * 1[                      v1^i  AND  v2^i ]

    Returns 0.0 when no axis has valid annotations for both samples.
    """
    numerator = 0.0
    denominator = 0.0

    for axis, w in weights.items():
        v1 = labels_s.get(axis)
        v2 = labels_i.get(axis)

        # validity flags: annotation exists and is not an "unknown" placeholder
        valid1 = (v1 is not None) and (str(v1).strip().lower() not in _UNKNOWN_VALUES)
        valid2 = (v2 is not None) and (str(v2).strip().lower() not in _UNKNOWN_VALUES)

        if valid1 and valid2:
            denominator += w
            if str(v1).strip() == str(v2).strip():
                numerator += w

    return numerator / denominator if denominator > 0.0 else 0.0


# ---------------------------------------------------------------------------
# Build fast look-up structures from the labels_map dict
# ---------------------------------------------------------------------------

def _build_label_lookup(
    labels_map: Dict[str, Dict[str, str]],
    samples: list[str],
) -> Dict[str, Dict[str, Optional[str]]]:
    """
    Returns {sample_id: {axis: value, ...}} for every sample in `samples`.
    Missing entries are stored as None.
    """
    lookup: Dict[str, Dict[str, Optional[str]]] = {}
    for s in samples:
        lookup[s] = {axis: labels_map[axis].get(s) for axis in labels_map}
    return lookup


# ---------------------------------------------------------------------------
# PCA distance helper
# ---------------------------------------------------------------------------

def _pca_distances(
    expr_df: pd.DataFrame,
    n_components: int = 50,
    use_precomputed_pca: Optional[np.ndarray] = None,
) -> tuple[np.ndarray, list[str]]:
    """
    Project `expr_df` (genes × samples) into PCA space and return:
        - distance matrix  (n_samples × n_samples), Euclidean
        - ordered list of sample IDs matching the matrix rows/cols
    """
    samples = list(expr_df.columns)
    X = expr_df.values.T  # → (n_samples, n_genes)

    if use_precomputed_pca is not None:
        coords = use_precomputed_pca
    else:
        n_comp = min(n_components, X.shape[0] - 1, X.shape[1])
        pca = PCA(n_components=n_comp, random_state=42)
        coords = pca.fit_transform(X)  # (n_samples, n_comp)

    dist_matrix = pairwise_distances(coords, metric="euclidean")
    return dist_matrix, samples

def compute_mean_pairwise_distance(dist_matrix: np.ndarray) -> float:
    """
    Compute the global mean pairwise distance Dist_bar(S).
 
        Dist_bar(S) = 2 / (|S| * (|S| - 1))  *  sum_{j < k} Dist(x_j, x_k)
 
    This is simply the mean of the upper triangle of the full distance matrix,
    which is equivalent to the closed-form expression above.
 
    Parameters
    ----------
    dist_matrix : np.ndarray
        Square (n x n) pairwise distance matrix.
 
    Returns
    -------
    float — mean pairwise distance across all distinct sample pairs.
    """
    n = dist_matrix.shape[0]
    if n < 2:
        warnings.warn("Cannot compute mean pairwise distance with fewer than 2 samples.")
        return float("nan")
 
    upper_tri = dist_matrix[np.triu_indices(n, k=1)]   # j < k pairs only
    return float(upper_tri.mean())


# ---------------------------------------------------------------------------
# Main metric function
# ---------------------------------------------------------------------------
def compute_global_distance_metrics(
    expr_df: pd.DataFrame,
    labels_map: Dict[str, Dict[str, str]],
    study_map: pd.DataFrame,
    axis_weights: Optional[Dict[str, float]] = None,
    n_pca_components: int = 50,
    precomputed_pca: Optional[np.ndarray] = None,
    verbose: bool = True,
) -> Dict[str, float]:

    if axis_weights is None:
        axis_weights = DEFAULT_AXIS_WEIGHTS

    # samples = list(expr_df.columns)

    # PCA distance matrix
    dist_matrix, ordered_samples = _pca_distances(
        expr_df,
        n_pca_components,
        precomputed_pca,
    )

    label_lookup = _build_label_lookup(
        labels_map,
        ordered_samples,
    )

    dist_bar = compute_mean_pairwise_distance(dist_matrix)

    similar_sum = 0.0
    dissimilar_sum = 0.0

    n_pairs = 0

    for i in range(len(ordered_samples)):
        sample_i = ordered_samples[i]
        labels_i = label_lookup[sample_i]

        for j in range(i + 1, len(ordered_samples)):
            sample_j = ordered_samples[j]
            labels_j = label_lookup[sample_j]

            similarity = compute_sim(
                labels_i,
                labels_j,
                axis_weights,
            )

            distance = dist_matrix[i, j]

            similar_sum += distance * similarity
            dissimilar_sum += distance * (1.0 - similarity)

            n_pairs += 1

    separation_score = (
        similar_sum / dissimilar_sum
        if dissimilar_sum > 0
        else float("nan")
    )

    if verbose:
        print(
            f"[DistMetrics] Dist_bar(S) = {dist_bar:.6f}"
        )
        print(
            f"[DistMetrics] BiologicalSeparation = "
            f"{separation_score:.6f}"
        )

    return {
        "Dist_bar": dist_bar,
        "BiologicalSeparation": separation_score,
    }

def compute_global_distance_metrics_ratio(
    expr_df: pd.DataFrame,
    labels_map: Dict[str, Dict[str, str]],
    study_map: pd.DataFrame,
    axis_weights: Optional[Dict[str, float]] = None,
    n_pca_components: int = 50,
    precomputed_pca: Optional[np.ndarray] = None,
    sim_threshold: float = 0.0,
    verbose: bool = True,
) -> Dict[str, float]:
    """
    Compute Dist_bar(S), G_d_inter, G_d_intra, and Ratio_global.
 
    Both G_d_inter and G_d_intra are normalised by Dist_bar(S) — the mean
    pairwise distance across all distinct sample pairs — making them
    scale-invariant and comparable across different expression matrices.
 
    Parameters
    ----------
    expr_df : pd.DataFrame
        Gene expression matrix — index = genes, columns = sample IDs.
    labels_map : dict
        Nested dict  {axis_name: {sample_id: label_value, ...}, ...}
        as produced by make_df_from_labels / load_labels_study.
    study_map : pd.DataFrame
        DataFrame with index = sample_id, column "StudyID" = study identifier.
        Matches SAMPLE_STUDY_MAP from constants.py.
    axis_weights : dict, optional
        Per-axis weights for sim().  Defaults to DEFAULT_AXIS_WEIGHTS.
    n_pca_components : int
        Number of PCA dimensions for the distance space.
    precomputed_pca : np.ndarray, optional
        Pre-computed PCA coordinates (n_samples × k).  If supplied,
        the PCA step is skipped.
    sim_threshold : float
        Only include (s, x_i) pairs whose sim() >= sim_threshold in the
        averages.  Default 0.0 keeps all pairs (matches the paper formula).
    verbose : bool
        Print progress info.
 
    Returns
    -------
    dict with keys:
        "Dist_bar"       — global mean pairwise distance (normalisation factor)
        "G_d_inter"      — global normalised inter-study weighted distance
        "G_d_intra"      — global normalised intra-study weighted distance
        "Ratio_global"   — G_d_intra / G_d_inter
    """
    if axis_weights is None:
        axis_weights = DEFAULT_AXIS_WEIGHTS
 
    samples = list(expr_df.columns)
    n = len(samples)
 
    if verbose:
        print(f"[DistMetrics] {n} samples | PCA dims: {n_pca_components}")
 
    # --- 1. Distance matrix in PCA space -----------------------------------
    dist_matrix, ordered_samples = _pca_distances(
        expr_df, n_pca_components, precomputed_pca
    )
    idx = {s: i for i, s in enumerate(ordered_samples)}
 
    # --- 2. Global normalisation factor  Dist_bar(S) -----------------------
    dist_bar = compute_mean_pairwise_distance(dist_matrix)
    if verbose:
        print(f"  Dist_bar(S)  = {dist_bar:.6f}  (mean of {n*(n-1)//2} pairwise distances)")
 
    if not np.isfinite(dist_bar) or dist_bar == 0.0:
        warnings.warn("Dist_bar(S) is zero or undefined — normalisation skipped.")
        norm = 1.0
    else:
        norm = dist_bar
 
    # --- 3. Label and study look-ups ----------------------------------------
    label_lookup = _build_label_lookup(labels_map, ordered_samples)
 
    def get_study(sample: str) -> str:
        if sample in study_map.index:
            return str(study_map.at[sample, "StudyID"])
        return "Unknown_Study"
 
    study_lookup = {s: get_study(s) for s in ordered_samples}
 
    # --- 4. Per-sample accumulation -----------------------------------------
    inter_contributions = []
    intra_contributions = []
 
    for s in ordered_samples:
        s_idx = idx[s]
        s_study = study_lookup[s]
        s_labels = label_lookup[s]
 
        inter_weighted_dist = 0.0
        inter_count = 0
        intra_weighted_dist = 0.0
        intra_count = 0
 
        for xi in ordered_samples:
            if xi == s:
                continue
 
            xi_study = study_lookup[xi]
            xi_labels = label_lookup[xi]
 
            similarity = compute_sim(s_labels, xi_labels, axis_weights)
 
            if similarity < sim_threshold:
                continue
 
            distance = dist_matrix[s_idx, idx[xi]]
            weighted = distance * similarity
 
            if xi_study != s_study:
                # Inter-study pool  X_inter(s)
                inter_weighted_dist += weighted
                inter_count += 1
            else:
                # Intra-study pool  X_intra(s)
                intra_weighted_dist += weighted
                intra_count += 1
 
        if inter_count > 0:
            inter_contributions.append(inter_weighted_dist / inter_count)
 
        if intra_count > 0:
            intra_contributions.append(intra_weighted_dist / intra_count)
 
    # --- 5. Global means, divided by Dist_bar(S) ----------------------------
    #
    #   The formula divides the outer sum by  |S| * Dist_bar(S).
    #   np.mean() handles the 1/|S| part; we then divide by norm = Dist_bar(S).
    #
    if not inter_contributions:
        warnings.warn("No inter-study pairs found — G_d_inter is undefined.")
        G_inter = float("nan")
    else:
        G_inter = float(np.mean(inter_contributions)) / norm
 
    if not intra_contributions:
        warnings.warn("No intra-study pairs found — G_d_intra is undefined.")
        G_intra = float("nan")
    else:
        G_intra = float(np.mean(intra_contributions)) / norm
 
    # ratio = G_intra / G_inter if (np.isfinite(G_inter) and G_inter > 0) else float("nan")

    separation_score = (
        (G_inter - G_intra) / (G_inter + G_intra)
        if (
            np.isfinite(G_inter)
            and np.isfinite(G_intra)
            and (G_inter + G_intra) > 0
        )
        else float("nan")
    )
 
    if verbose:
        print(f"  G_d_inter    = {G_inter:.6f}  (n_samples contributing: {len(inter_contributions)})")
        print(f"  G_d_intra    = {G_intra:.6f}  (n_samples contributing: {len(intra_contributions)})")
        print(f"  SeparationScore  = {separation_score:.4f}") 
    return {
        "Dist_bar":    dist_bar,
        "G_d_inter":   G_inter,
        "G_d_intra":   G_intra,
        "SeparationScore":  separation_score,
    }




# ---------------------------------------------------------------------------
# Convenience wrapper that plugs directly into the existing pipeline loop
# ---------------------------------------------------------------------------

def run_distance_evaluation(
    data_df: pd.DataFrame,
    labels_dict: Dict[str, Dict[str, str]],
    sample_study_map: pd.DataFrame,
    experiment_name: str = "",
    n_pca_components: int = 50,
    axis_weights: Optional[Dict[str, float]] = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Drop-in wrapper for the stage loop in __main__.

    Returns a single-row pd.DataFrame with columns
    ['experiment', 'G_d_inter', 'G_d_intra', 'Ratio_global'].

    Usage inside the existing loop
    --------------------------------
    from src.data_analisys.distance_evaluation_metrics import run_distance_evaluation

    dist_metrics = run_distance_evaluation(
        data_df=df,
        labels_dict=labels_map,
        sample_study_map=SAMPLE_STUDY_MAP,
        experiment_name=file,
    )
    all_dist_metrics[file] = dist_metrics
    """
    if verbose:
        print(f"\n[DistMetrics] Running for stage: '{experiment_name}'")
    if RNA_USED:
        data_df.columns = [get_gsm_id(col.split('_')[1]) for col in data_df.columns]
    if RNA_USED:
        try:
            sample_study_map.index = [get_gsm_id(ind.split('_')[1]) for ind in sample_study_map.index]
        except IndexError:
            pass
    results = compute_global_distance_metrics(
        expr_df=data_df,
        labels_map=labels_dict,
        study_map=sample_study_map,
        axis_weights=axis_weights,
        n_pca_components=n_pca_components,
        verbose=verbose,
    )

    row = {"experiment": experiment_name, **results}
    return pd.DataFrame([row])