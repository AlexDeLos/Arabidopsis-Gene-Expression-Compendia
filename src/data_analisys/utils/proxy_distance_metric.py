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
from scipy.stats import spearmanr, pearsonr
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
def _parse_treatments(value: str) -> set[str]:
    return {
        t.strip().lower()
        for t in str(value).split("+")
        if t.strip()
    }

def _treatment_similarity(
    t1: str,
    t2: str,
) -> float:

    s1 = _parse_treatments(t1)
    s2 = _parse_treatments(t2)

    if not s1 or not s2:
        return 0.0

    return len(s1 & s2) / len(s1 | s2)

def _intensity_similarity(
    treatment1: str,
    treatment2: str,
    intensity1: str,
    intensity2: str,
) -> float:

    overlap = _treatment_similarity(
        treatment1,
        treatment2,
    )

    if overlap == 0:
        return 0.0

    try:
        i1 = float(intensity1)
        i2 = float(intensity2)

        max_i = max(i1, i2)

        if max_i == 0:
            return overlap

        intensity_sim = 1.0 - abs(i1 - i2) / max_i

        return overlap * intensity_sim

    except Exception:
        return overlap
    
def compute_sim(
    labels_s: Dict[str, Optional[str]],
    labels_i: Dict[str, Optional[str]],
    weights: Dict[str, float],
) -> float:
    """
    Compute the weighted metadata similarity score sim(L_s, L_i).

    Each metadata axis contributes according to its biological importance
    (weight). An axis only contributes when both samples have valid
    (non-unknown) annotations.

    For most axes, similarity is binary:

        similarity_axis =
            1.0  if L_s^i == L_i^i
            0.0  otherwise

    Two special cases are handled:

    1. Treatment
    Treatment labels may contain multiple treatments separated by "+".
    Similarity is computed using the Jaccard overlap between treatment sets:

        sim_treatment =
            |T_s ∩ T_i| / |T_s ∪ T_i|

    Examples:
        Heat vs Heat               -> 1.0
        Heat vs Salt               -> 0.0
        Heat+Salt vs Heat          -> 0.5
        Heat+Salt vs Heat+Salt     -> 1.0

    2. Treatment Intensity
    Treatment intensity is only considered similar when the corresponding
    treatment labels share at least one treatment. Intensity similarity is
    computed from the relative difference between intensity values and
    contributes a fractional score between 0 and 1.

    The final similarity score is the weighted average across all valid axes:

                    Σ_i w_i · sim_i
        sim(L_s,L_i)= ----------------
                        Σ_i w_i

    where sim_i is either a binary match (standard labels) or a fractional
    similarity (treatment-related labels).

    Returns
    -------
    float
        Similarity score in the range [0, 1].

        1.0 indicates identical metadata across all valid axes.
        0.0 indicates no metadata similarity across valid axes.

        Returns 0.0 when no axis has valid annotations for both samples.
    """
    numerator = 0.0
    denominator = 0.0

    for axis, w in weights.items():

        v1 = labels_s.get(axis)
        v2 = labels_i.get(axis)

        valid1 = (
            v1 is not None
            and str(v1).strip().lower() not in _UNKNOWN_VALUES
        )

        valid2 = (
            v2 is not None
            and str(v2).strip().lower() not in _UNKNOWN_VALUES
        )

        if not (valid1 and valid2):
            continue

        denominator += w

        # -----------------------------
        # treatment
        # -----------------------------

        if axis == "treatment":

            numerator += (
                w *
                _treatment_similarity(
                    str(v1),
                    str(v2),
                )
            )

        # -----------------------------
        # treatment intensity
        # -----------------------------

        elif axis == "treatment_intensity":

            numerator += (
                w *
                _intensity_similarity(
                    labels_s.get("treatment", ""),
                    labels_i.get("treatment", ""),
                    str(v1),
                    str(v2),
                )
            )

        # -----------------------------
        # standard labels
        # -----------------------------

        elif str(v1).strip() == str(v2).strip():

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
    Project expr_df (samples × genes) into PCA space and return:

        - distance matrix (n_samples × n_samples)
        - ordered list of sample IDs matching matrix rows/cols
    """

    samples = list(expr_df.index)

    # Samples × Genes
    X = expr_df.values

    if use_precomputed_pca is not None:

        coords = use_precomputed_pca

        if coords.shape[0] != len(samples):
            raise ValueError(
                "Precomputed PCA rows do not match number of samples. "
                f"Got {coords.shape[0]} rows for {len(samples)} samples."
            )

    else:

        n_comp = min(
            n_components,
            X.shape[0] - 1,
            X.shape[1],
        )

        pca = PCA(n_components=n_comp)

        # returns (n_samples, n_components)
        coords = pca.fit_transform(X)

    dist_matrix = pairwise_distances(
        coords,
        metric="euclidean",
    )

    return dist_matrix, samples

def compute_mean_pairwise_distance(
    dist_matrix: np.ndarray,
) -> float:
    """
    Mean distance across all unique sample pairs.
    """

    n = dist_matrix.shape[0]

    if n < 2:
        warnings.warn(
            "Cannot compute mean pairwise distance "
            "with fewer than 2 samples."
        )
        return float("nan")

    total_distance = dist_matrix.sum() / 2.0

    n_pairs = n * (n - 1) / 2.0

    return float(total_distance / n_pairs)

# ---------------------------------------------------------------------------
# Main metric function
# ---------------------------------------------------------------------------



def compute_global_distance_metrics(
    expr_df: pd.DataFrame,
    labels_map: Dict[str, Dict[str, str]],
    study_map: pd.DataFrame,
    axis_weights: Optional[Dict[str, float]] = None,
    n_pca_components: int | None = None,
    precomputed_pca: Optional[np.ndarray] = None,
    verbose: bool = True,
) -> Dict[str, object]:
    """
    Parameters
    ----------
    expr_df
        Samples × Genes expression matrix.

        Rows:
            Samples

        Columns:
            Genes/features
    """

    if axis_weights is None:
        axis_weights = DEFAULT_AXIS_WEIGHTS

    # --------------------------------------------------
    # PCA distance matrix
    # --------------------------------------------------
    if n_pca_components:
        dist_matrix, ordered_samples = _pca_distances(
            expr_df,
            n_pca_components,
            precomputed_pca,
        )
    else:
        # Samples × Genes convention:
        # rows = samples
        # columns = genes
        ordered_samples = list(expr_df.index)

        dist_matrix = pairwise_distances(
            expr_df,
            metric="euclidean",
        )
    print("distance matrix built")
    # --------------------------------------------------
    # Sanity checks
    # --------------------------------------------------
    if dist_matrix.shape != (
        len(ordered_samples),
        len(ordered_samples),
    ):
        raise ValueError(
            "Distance matrix shape does not match number "
            f"of samples. Got dist_matrix.shape="
            f"{dist_matrix.shape} and "
            f"{len(ordered_samples)} samples."
        )
    label_lookup = _build_label_lookup(
        labels_map,
        ordered_samples,
    )

    study_lookup = {
        sample: str(study_map.at[sample, "StudyID"])
        if sample in study_map.index
        else "Unknown_Study"
        for sample in ordered_samples
    }
    print("labels set up")
    
    dist_bar = compute_mean_pairwise_distance(
        dist_matrix
    )
    print("mean distance calculated")

    # --------------------------------------------------
    # Biological separation metric
    # --------------------------------------------------

    similar_sum = 0.0
    dissimilar_sum = 0.0

    # --------------------------------------------------
    # Inter / intra study metrics
    # --------------------------------------------------

    inter_weighted_sum = 0.0
    intra_weighted_sum = 0.0

    inter_pairs = 0
    intra_pairs = 0

    # --------------------------------------------------
    # Similarity-distance correlation data
    # --------------------------------------------------

    pairwise_records = []

    # --------------------------------------------------
    # Pairwise loop
    # --------------------------------------------------

    for i in range(len(ordered_samples)):

        sample_i = ordered_samples[i]
        labels_i = label_lookup[sample_i]
        study_i = study_lookup[sample_i]

        for j in range(i + 1, len(ordered_samples)):

            sample_j = ordered_samples[j]
            labels_j = label_lookup[sample_j]
            study_j = study_lookup[sample_j]

            similarity = compute_sim(
                labels_i,
                labels_j,
                axis_weights,
            )

            distance = float(dist_matrix[i, j])

            # ------------------------------------------
            # Store pair for plotting/correlation
            # ------------------------------------------

            pairwise_records.append(
                {
                    "Sample_A": sample_i,
                    "Sample_B": sample_j,
                    "Study_A": study_i,
                    "Study_B": study_j,
                    "Distance": distance,
                    "Similarity": float(similarity),
                    "SameStudy": study_i == study_j,
                }
            )

            # ------------------------------------------
            # Biological separation
            # ------------------------------------------

            similar_sum += distance * similarity

            dissimilar_sum += (
                distance
                * (1.0 - similarity)
            )

            # ------------------------------------------
            # Inter / intra study metrics
            # ------------------------------------------

            weighted_distance = (
                distance * similarity
            )

            if study_i == study_j:

                intra_weighted_sum += (
                    weighted_distance
                )
                intra_pairs += 1

            else:

                inter_weighted_sum += (
                    weighted_distance
                )
                inter_pairs += 1

    # --------------------------------------------------
    # Pairwise dataframe
    # --------------------------------------------------

    pairwise_df = pd.DataFrame(
        pairwise_records
    )

    # --------------------------------------------------
    # Similarity-distance correlations
    # --------------------------------------------------

    if len(pairwise_df) > 1:

        try:
            spearman_corr, spearman_p = spearmanr(
                pairwise_df["Distance"],
                pairwise_df["Similarity"],
            )
            print(f"got a spearman_p value of {spearman_p}")
        except Exception:
            spearman_corr = np.nan
            spearman_p = np.nan

        try:
            pearson_corr, pearson_p = pearsonr(
                pairwise_df["Distance"],
                pairwise_df["Similarity"],
            )
            print(f"got a pearson_p value of {pearson_p}")
        except Exception:
            pearson_corr = np.nan
            pearson_p = np.nan

    else:

        spearman_corr = np.nan
        spearman_p = np.nan

        pearson_corr = np.nan
        pearson_p = np.nan

    # --------------------------------------------------
    # Biological separation
    # --------------------------------------------------

    biological_separation = (
        similar_sum / dissimilar_sum
        if dissimilar_sum > 0
        else float("nan")
    )

    # --------------------------------------------------
    # G_d_inter
    # --------------------------------------------------

    G_d_inter = (
        (inter_weighted_sum / inter_pairs)
        / dist_bar
        if inter_pairs > 0 and dist_bar > 0
        else float("nan")
    )

    # --------------------------------------------------
    # G_d_intra
    # --------------------------------------------------

    G_d_intra = (
        (intra_weighted_sum / intra_pairs)
        / dist_bar
        if intra_pairs > 0 and dist_bar > 0
        else float("nan")
    )

    # --------------------------------------------------
    # Logging
    # --------------------------------------------------

    if verbose:

        print(
            f"[DistMetrics] Dist_bar(S) = "
            f"{dist_bar:.6f}"
        )

        print(
            f"[DistMetrics] G_d_inter = "
            f"{G_d_inter:.6f}"
        )

        print(
            f"[DistMetrics] G_d_intra = "
            f"{G_d_intra:.6f}"
        )

        print(
            f"[DistMetrics] BiologicalSeparation = "
            f"{biological_separation:.6f}"
        )

        print(
            f"[DistMetrics] Spearman(sim,dist) = "
            f"{spearman_corr:.6f}, p of {spearman_p}"
        )

        print(
            f"[DistMetrics] Pearson(sim,dist) = "
            f"{pearson_corr:.6f} p of {pearson_p:.6f}"
        )

    return {
        "Dist_bar": dist_bar,
        "G_d_inter": G_d_inter,
        "G_d_intra": G_d_intra,
        "BiologicalSeparation": biological_separation,
        "SimilarityDistanceSpearman": float(
            spearman_corr
        ),
        "SimilarityDistanceSpearmanP": float(
            spearman_p
        ),
        "SimilarityDistancePearson": float(
            pearson_corr
        ),
        "SimilarityDistancePearsonP": float(
            pearson_p
        ),
        "PairwiseSimilarityDistanceDF": pairwise_df,
    }



# ---------------------------------------------------------------------------
# Convenience wrapper that plugs directly into the existing pipeline loop
# ---------------------------------------------------------------------------

def run_distance_evaluation(
    data_df: pd.DataFrame,
    labels_dict: Dict[str, Dict[str, str]],
    sample_study_map: pd.DataFrame,
    experiment_name: str = "",
    n_pca_components: int | None = None,
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
        data_df.index = [get_gsm_id(col.split('_')[1]) for col in data_df.index]
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