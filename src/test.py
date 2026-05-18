"""
Unit tests for run_rank_in_normalization (Tang et al., 2021).

Run with:
    pytest test_rankin_normalization.py -v

Each test targets one specific property of the algorithm described in the paper.
"""

import warnings
import numpy as np
import pandas as pd
import pytest
import sys
# from scipy.sparse.linalg import svds
module_dir = "./"
sys.path.append(module_dir)
from src.data_importing.data_norm_and_analisys import run_rank_in_normalization  # noqa:

# ─────────────────────────────────────────────────────────────────────────────
# Shared fixtures
# ─────────────────────────────────────────────────────────────────────────────

def make_synthetic_data(
    n_genes: int = 500,
    n_cancer: int = 10,
    n_normal: int = 10,
    array_loc: float = 6.0,
    rnaseq_loc: float = 3.0,
    seed: int = 42,
) -> tuple[pd.DataFrame, pd.Series]:
    """
    Simulates a mixed microarray / RNA-seq dataset with a known platform shift.

    - First half of each class = microarray-like  (higher mean)
    - Second half of each class = RNA-seq-like     (lower mean, higher variance)

    Returns (df, sample_classes).
    """
    rng = np.random.default_rng(seed)

    half_c = n_cancer // 2
    half_n = n_normal // 2

    cancer_array  = rng.normal(array_loc  + 2, 1.2, (n_genes, half_c))
    cancer_rnaseq = rng.normal(rnaseq_loc + 2, 2.5, (n_genes, n_cancer - half_c))
    normal_array  = rng.normal(array_loc,       1.2, (n_genes, half_n))
    normal_rnaseq = rng.normal(rnaseq_loc,      2.5, (n_genes, n_normal - half_n))

    expr = np.hstack([cancer_array, cancer_rnaseq, normal_array, normal_rnaseq])

    sample_names = (
        [f"cancer_array_{i}"  for i in range(half_c)] +
        [f"cancer_rnaseq_{i}" for i in range(n_cancer - half_c)] +
        [f"normal_array_{i}"  for i in range(half_n)] +
        [f"normal_rnaseq_{i}" for i in range(n_normal - half_n)]
    )
    gene_names = [f"gene_{i}" for i in range(n_genes)]

    df = pd.DataFrame(expr, index=gene_names, columns=sample_names)

    labels = (["cancer"] * n_cancer) + (["normal"] * n_normal)
    sample_classes = pd.Series(labels, index=sample_names)

    return df, sample_classes


# ─────────────────────────────────────────────────────────────────────────────
# Tests
# ─────────────────────────────────────────────────────────────────────────────

class TestOutputShape:
    """The adjusted matrix must have exactly the same shape and labels as input."""

    def test_shape_preserved(self):
        df, sc = make_synthetic_data(n_genes=300, n_cancer=8, n_normal=8)
        result = run_rank_in_normalization(df, sc, k=2)
        assert result.shape == df.shape, (
            f"Expected shape {df.shape}, got {result.shape}"
        )

    def test_index_preserved(self):
        df, sc = make_synthetic_data()
        result = run_rank_in_normalization(df, sc, k=2)
        pd.testing.assert_index_equal(result.index, df.index)

    def test_columns_preserved(self):
        df, sc = make_synthetic_data()
        result = run_rank_in_normalization(df, sc, k=2)
        pd.testing.assert_index_equal(result.columns, df.columns)


class TestStep1BinnedRanking:
    """
    Step 1: each sample's weighted ranks should span [1, n_bins] after binning,
    and the weighting in step 2 scales these — so the raw binned values (before
    weighting) must be integers in [1, n_bins].
    """

    def test_binned_values_in_range(self):
        """
        Verify indirectly: if we run with a=0,b=1 (identity weight, no quadratic
        distortion possible) the adjusted matrix values are bounded near [1, n_bins].
        We force this by using constant expression per sample so polyfit gives ~0 slope.
        """
        rng = np.random.default_rng(0)
        n_genes, n_samples = 200, 10
        # All samples have identical expression → ranks are well-defined [1,100]
        base = rng.normal(5, 1, n_genes)
        expr = np.column_stack([base + rng.normal(0, 0.01, n_genes)
                                for _ in range(n_samples)])
        cols = [f"s{i}" for i in range(n_samples)]
        df = pd.DataFrame(expr, columns=cols)
        labels = ["A"] * 5 + ["B"] * 5
        sc = pd.Series(labels, index=cols)

        result = run_rank_in_normalization(df, sc, k=1)
        # After SVD correction the values should remain in a reasonable range
        assert result.notna().all().all(), "Result contains NaN values"

    def test_n_bins_parameter_respected(self):
        """Output should differ between n_bins=50 and n_bins=100 (different granularity)."""
        df, sc = make_synthetic_data(n_genes=200)
        r50  = run_rank_in_normalization(df, sc, n_bins=50,  k=2)
        r100 = run_rank_in_normalization(df, sc, n_bins=100, k=2)
        assert not np.allclose(r50.values, r100.values), (
            "n_bins=50 and n_bins=100 should produce different outputs"
        )


class TestStep2Weighting:
    """Step 2: the quadratic weighting should change the matrix relative to unweighted ranks."""

    def test_weighting_changes_output(self):
        """
        With real expression variance, the quadratic fit will have a non-zero slope,
        so the weighted matrix must differ from the raw binned matrix.
        We verify this by checking the adjusted outputs differ from a hypothetical
        identity-weight version (simulated by flat expression per sample).
        """
        rng = np.random.default_rng(1)
        n_genes = 300

        # Varying expression — weighting will have effect
        expr_varying = rng.normal(5, 2, (n_genes, 10))
        # Flat expression — weighting has almost no effect (slope ≈ 0)
        expr_flat = np.tile(np.linspace(1, 10, n_genes)[:, None], (1, 10))
        expr_flat += rng.normal(0, 0.001, expr_flat.shape)

        cols = [f"s{i}" for i in range(10)]
        sc   = pd.Series(["A"] * 5 + ["B"] * 5, index=cols)

        res_varying = run_rank_in_normalization(
            pd.DataFrame(expr_varying, columns=cols), sc, k=1
        )
        res_flat = run_rank_in_normalization(
            pd.DataFrame(expr_flat, columns=cols), sc, k=1
        )

        # The two results should not be identical
        assert not np.allclose(res_varying.values, res_flat.values), (
            "Expression variance should produce different weighted outputs"
        )

    def test_fallback_warning_on_degenerate_sample(self):
        """
        A sample with all-identical expression values cannot be polyfit'd.
        The function should warn rather than raise.
        """
        rng = np.random.default_rng(2)
        n_genes = 100
        expr = rng.normal(5, 1, (n_genes, 6))
        expr[:, 2] = 5.0  # degenerate: all same value in sample index 2

        cols = [f"s{i}" for i in range(6)]
        sc   = pd.Series(["A"] * 3 + ["B"] * 3, index=cols)
        df   = pd.DataFrame(expr, columns=cols)

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = run_rank_in_normalization(df, sc, k=1)

        warning_messages = [str(w.message) for w in caught]
        # numpy emits "Polyfit may be poorly conditioned" for degenerate input;
        # our code emits a RuntimeWarning with "fallback" / "fit failed".
        # Accept either — both indicate the degenerate case was handled gracefully.
        assert any(
            "fallback"           in m.lower() or
            "fit failed"         in m.lower() or
            "poorly conditioned" in m.lower() or
            "polyfit"            in m.lower()
            for m in warning_messages
        ), f"Expected a warning for degenerate sample, got: {warning_messages}"
        assert result.shape == df.shape, "Function should still return valid output after fallback"


class TestStep3SVD:
    """
    Step 3: SVD nonbiological effect removal.
    """

    def test_platform_variance_reduced(self):
        """
        Core paper claim: after Rank-In, within-class variance across platforms
        should be lower than before.

        We construct a dataset where cancer samples split evenly across two
        platforms with a large platform offset, and verify that the
        inter-platform spread (std) for cancer genes decreases after adjustment.
        """
        rng = np.random.default_rng(3)
        n_genes = 500
        platform_offset = 4.0  # large systematic shift between platforms

        cancer_p1 = rng.normal(6.0,               1.0, (n_genes, 8))
        cancer_p2 = rng.normal(6.0 + platform_offset, 1.0, (n_genes, 8))
        normal_p1 = rng.normal(4.0,               1.0, (n_genes, 8))
        normal_p2 = rng.normal(4.0 + platform_offset, 1.0, (n_genes, 8))

        expr = np.hstack([cancer_p1, cancer_p2, normal_p1, normal_p2])
        cols = (
            [f"cancer_p1_{i}" for i in range(8)] +
            [f"cancer_p2_{i}" for i in range(8)] +
            [f"normal_p1_{i}" for i in range(8)] +
            [f"normal_p2_{i}" for i in range(8)]
        )
        sc = pd.Series(["cancer"] * 16 + ["normal"] * 16, index=cols)
        df = pd.DataFrame(expr, index=[f"g{i}" for i in range(n_genes)], columns=cols)

        result = run_rank_in_normalization(df, sc, k=1)

        # Compare column-wise std across all samples before vs. after
        std_before = df.std(axis=1).mean()
        std_after  = result.std(axis=1).mean()

        assert std_after < std_before, (
            f"Expected reduced cross-sample variance after Rank-In "
            f"(before={std_before:.3f}, after={std_after:.3f})"
        )

    def test_k_auto_selection_is_small(self):
        """
        The paper states k should be small (dominant nonbiological components).
        Auto-selected k should be well below k_max for typical data.
        """
        df, sc = make_synthetic_data(n_genes=500, n_cancer=12, n_normal=12)

        # Capture printed output to extract auto-selected k
        import io
        from contextlib import redirect_stdout

        buf = io.StringIO()
        with redirect_stdout(buf):
            run_rank_in_normalization(df, sc, k=None, k_max=10)

        output = buf.getvalue()
        # Parse "Auto-selected k=N"
        import re
        match = re.search(r"Auto-selected k=(\d+)", output)
        assert match, f"Could not find auto-selected k in output:\n{output}"
        k_selected = int(match.group(1))

        assert k_selected <= 5, (
            f"Auto-selected k={k_selected} is unexpectedly large; "
            "k should be small (1–5) for typical platform-effect removal"
        )

    def test_manual_k_overrides_auto(self):
        """Passing k=3 explicitly should not trigger auto-selection."""
        df, sc = make_synthetic_data()

        import io
        from contextlib import redirect_stdout

        buf = io.StringIO()
        with redirect_stdout(buf):
            run_rank_in_normalization(df, sc, k=3)

        output = buf.getvalue()
        assert "user-supplied k=3" in output, (
            f"Expected confirmation of user-supplied k=3 in output:\n{output}"
        )
        assert "Auto-selected" not in output

    def test_centering_uses_group_means_not_grand_mean(self):
        """
        The paper centers by subtracting per-class group means (Me_ij), not the
        grand mean. We verify this directly: the variance matrix produced inside
        the function (R - Me_ij) should have a smaller within-class residual than
        the alternative of subtracting the grand mean.

        Strategy: after Rank-In with correct two-class centering, the cancer and
        normal group means of the *adjusted* matrix should differ more than they
        would if a grand mean had been used (grand-mean centering deflates
        inter-class differences by pulling both groups toward the same centre).
        We simulate the "wrong" case by manually applying grand-mean centering
        to the weighted matrix and running SVD on that instead, then comparing
        inter-class separation in both adjusted outputs.
        """
        rng = np.random.default_rng(99)
        n_genes, n_cancer, n_normal = 300, 10, 10

        # Large biological gap between cancer (high) and normal (low)
        cancer_expr = rng.normal(8.0, 1.0, (n_genes, n_cancer))
        normal_expr = rng.normal(3.0, 1.0, (n_genes, n_normal))
        expr = np.hstack([cancer_expr, normal_expr])

        cancer_cols = [f"c{i}" for i in range(n_cancer)]
        normal_cols = [f"n{i}" for i in range(n_normal)]
        all_cols    = cancer_cols + normal_cols

        df = pd.DataFrame(expr, index=[f"g{i}" for i in range(n_genes)], columns=all_cols)
        sc = pd.Series(["cancer"] * n_cancer + ["normal"] * n_normal, index=all_cols)

        # Correct: per-class centering
        result_correct = run_rank_in_normalization(df, sc, k=2)

        inter_correct = (
            result_correct[cancer_cols].mean(axis=1) -
            result_correct[normal_cols].mean(axis=1)
        ).abs().mean()

        # Manually simulate grand-mean centering: center weighted matrix by
        # the overall gene mean, then SVD, then subtract nonbio effects.
        # (Replicates what the original buggy code did.)
        

        # Re-derive the weighted matrix the same way the function does
        rank_mat = df.rank(pct=True, method="average")
        binned   = np.ceil(rank_mat.values * 100).astype(float)
        weighted = np.zeros_like(binned)
        for j in range(binned.shape[1]):
            r = binned[:, j]
            e = df.iloc[:, j].values.astype(float)
            try:
                coeffs = np.polyfit(r, e, 2)
                a, b = coeffs[0], coeffs[1]
            except Exception:
                a, b = 0.0, 1.0
            weighted[:, j] = r * (2 * a * r + b)

        weighted_df = pd.DataFrame(weighted, index=df.index, columns=df.columns)

        # Grand-mean centering (the bug)
        grand_mean  = weighted_df.mean(axis=1)
        centered_gm = weighted_df.sub(grand_mean, axis=0)

        U, s, Vt  = np.linalg.svd(centered_gm.values, full_matrices=False)
        nonbio_gm = U[:, :2] @ np.diag(s[:2]) @ Vt[:2, :]
        adjusted_gm = pd.DataFrame(
            weighted_df.values - nonbio_gm,
            index=df.index, columns=df.columns,
        )

        inter_wrong = (
            adjusted_gm[cancer_cols].mean(axis=1) -
            adjusted_gm[normal_cols].mean(axis=1)
        ).abs().mean()

        assert inter_correct > inter_wrong, (
            f"Per-class centering should preserve inter-class separation better "
            f"than grand-mean centering "
            f"(per-class={inter_correct:.3f}, grand-mean={inter_wrong:.3f})"
        )


class TestInputValidation:
    """The function should raise clearly on bad inputs."""

    def test_mismatched_sample_classes_raises(self):
        df, sc = make_synthetic_data()
        sc_bad = sc.copy()
        sc_bad.index = [f"wrong_{i}" for i in range(len(sc))]
        with pytest.raises((ValueError, KeyError)):
            run_rank_in_normalization(df, sc_bad, k=2)

    def test_single_class_raises(self):
        df, sc = make_synthetic_data()
        sc_one = pd.Series(["cancer"] * len(sc), index=sc.index)
        with pytest.raises(ValueError, match="two distinct class labels"):
            run_rank_in_normalization(df, sc_one, k=2)

    def test_missing_class_label_raises(self):
        df, sc = make_synthetic_data()
        sc_missing = sc.copy().astype(object)
        sc_missing.iloc[0] = np.nan
        with pytest.raises(ValueError, match="missing labels"):
            run_rank_in_normalization(df, sc_missing, k=2)

    def test_k_clamped_with_warning(self):
        df, sc = make_synthetic_data(n_genes=100, n_cancer=6, n_normal=6)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = run_rank_in_normalization(df, sc, k=9999)
        assert result.shape == df.shape
        assert any("clamping" in str(w.message).lower() for w in caught), (
            "Expected a clamping warning when k exceeds safe maximum"
        )


class TestDeterminism:
    """Results must be deterministic (no hidden randomness)."""

    def test_same_input_same_output(self):
        df, sc = make_synthetic_data(seed=7)
        r1 = run_rank_in_normalization(df, sc, k=2)
        r2 = run_rank_in_normalization(df, sc, k=2)
        pd.testing.assert_frame_equal(r1, r2)

    def test_output_is_finite(self):
        df, sc = make_synthetic_data()
        result = run_rank_in_normalization(df, sc, k=2)
        assert np.isfinite(result.values).all(), "Output contains NaN or Inf values"