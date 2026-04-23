"""
diff_expr.py
------------
Runs limma-based differential expression analysis for each treatment vs. Control,
optionally restricted to a single tissue and/or a fixed set of samples.

Entry point: diff_exp_combine_tissues(...)
"""

import os
import sys
import re
import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import Formula, pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

module_dir = "./"
sys.path.append(module_dir)
from src.constants import CLUSTER_RUN  # noqa: E402
from src.constants_labeling import TreatmentEnum  # noqa: E402


def diff_exp_combine_tissues(
    treatments: list[str],
    save_dir: str,
    data_type: str,
    design: pd.DataFrame,
    out_dir: str,
    samples: list[str] | None = None,
    pure: bool = False,
    tissue: str | None = None,
    filter_low_combination: int = 10,
) -> None:
    """
    Runs limma differential expression for every treatment in `treatments`.

    For each treatment a CSV of differentially expressed genes is written to
    `out_dir` with columns [ID, t, logFC, adj.P.Val], sorted by adjusted p-value.

    Parameters
    ----------
    treatments : list of str
        Treatment labels to test (each tested vs. Control).
    save_dir : str
        Directory that contains `{data_type}.csv` (genes × samples).
    data_type : str
        Name of the expression matrix file (without .csv extension).
    design : pd.DataFrame
        Labels DataFrame with at least columns [sample_id, treatment, tissue].
        Index must match sample IDs used as column headers in the expression matrix.
    out_dir : str
        Directory where output files are written; created if it does not exist.
    samples : list of str or None
        If given, restrict analysis to these sample IDs (sanity-check mode).
    pure : bool
        If True, only keep samples whose treatment exactly matches the target
        treatment (no compound treatments).
    tissue : str or None
        If given, restrict to samples whose tissue label contains this string.
    filter_low_combination : int
        Minimum number of samples required per (treatment × tissue) group.
        Groups below this threshold are dropped before model fitting.
    """
    os.makedirs(out_dir, exist_ok=True)

    for treatment in treatments:
        try:
            tissue_label = tissue if tissue is not None else "All-Tissues"
            output_filename = f"{tissue_label}_{treatment}"
            file_to_output = f"{out_dir}{output_filename}_genes.csv"

            print(f"--- Starting analysis for treatment: {treatment} ---")
            if os.path.isfile(file_to_output):
                print("    File already exists, skipping.")
                continue

            # ------------------------------------------------------------------
            # 1. Load expression data
            # ------------------------------------------------------------------
            data = pd.read_csv(f"{save_dir}/{data_type}.csv", index_col=0)

            # ------------------------------------------------------------------
            # 2. Build sample masks
            # ------------------------------------------------------------------
            is_tissue = (
                design["tissue"].str.contains(tissue, na=False)
                if tissue
                else design["tissue"].str.contains("", na=False)
            )
            is_treatment = design["treatment"].str.contains(treatment, na=False)
            is_only_treatment = design["treatment"].apply(
                lambda x: len(x) == len(treatment) and treatment in x  # noqa: B023
            )
            if pure:
                is_treatment = is_only_treatment
            is_control = design["treatment"].str.contains(TreatmentEnum.CONTROL, na=False)

            if samples is not None:
                is_study = design["Sample_ID"].apply(lambda x: x in samples)
                design_filtered = design[(is_treatment | is_control) & is_study & is_tissue].copy()
            else:
                design_filtered = design[(is_treatment | is_control) & is_tissue].copy()

            # ------------------------------------------------------------------
            # 3. Synchronise samples between expression data and design
            # ------------------------------------------------------------------
            common_samples = list(set(design_filtered.index) & set(data.columns))
            design_filtered = design_filtered[design_filtered.index.isin(common_samples)]
            data_filtered = data[common_samples]

            # 1. Drop duplicate columns from data
            data_filtered = data_filtered.T.drop_duplicates().T
            
            # 2. IMPORTANT: Update design_filtered to ONLY include samples 
            # that still exist in data_filtered columns
            design_filtered = design_filtered[design_filtered.index.isin(data_filtered.columns)]

            # 3. Now drop duplicate rows from design (if any)
            design_filtered = design_filtered.drop_duplicates()

            # 4. Filter by count combination
            design_filtered = (
                design_filtered.sort_values(by="Sample_ID")
                .reset_index(drop=False)
                .groupby(["treatment", "tissue"])
                .filter(lambda x: len(x) >= filter_low_combination)
            )
            
            # 5. This will now work because both are perfectly synced
            data_filtered = data_filtered[design_filtered["Sample_ID"]]

            print(f"    {len(design_filtered)} samples aligned for '{treatment}' vs. {TreatmentEnum.CONTROL}.")

            # ------------------------------------------------------------------
            # 4. Build metadata for the model
            # ------------------------------------------------------------------
            metadata = design_filtered[["Sample_ID", "treatment", "tissue"]].copy()
            del design_filtered
            def r_make_names(s):
                # 1. Replace any non-alphanumeric character (like spaces) with a dot
                s = re.sub(r'[^a-zA-Z0-9_]', '.', str(s))
                # 2. If it starts with a digit or a dot, prepend 'X' (R requirement)
                if s[0].isdigit() or s[0] == '.':
                    s = 'X' + s
                return s
            is_control_mask = metadata["treatment"].str.contains("Control", na=False)
            metadata["Target"] = "treatment"
            metadata.loc[is_control_mask, "Target"] = "Control"
            del metadata["treatment"]

            metadata["Target"] = metadata["Target"].apply(r_make_names)
            
            if "tissue" in metadata.columns:
                metadata["tissue"] = metadata["tissue"].apply(r_make_names)

            single_tissue = len(set(metadata["tissue"])) == 1
            if single_tissue:
                del metadata["tissue"]

            if len(set(metadata["Target"])) == 1:
                print(f"    Not enough unique Target levels for '{treatment}', skipping.")
                continue

            # ------------------------------------------------------------------
            # 5. Import R libraries
            # ------------------------------------------------------------------
            base = importr("base")
            stats = importr("stats")
            if CLUSTER_RUN:
                limma = importr("limma")
                writexl = importr("writexl")
            else:
                limma = importr("limma", lib_loc="/home/alex/R/x86_64-pc-linux-gnu-library/4.5/")
                writexl = importr("writexl", lib_loc="/home/alex/R/x86_64-pc-linux-gnu-library/4.5/")

            # ------------------------------------------------------------------
            # 6. Convert Python objects to R
            # ------------------------------------------------------------------
            with localconverter(ro.default_converter + pandas2ri.converter):
                r_data = ro.conversion.py2rpy(data_filtered)
                r_metadata = ro.conversion.py2rpy(metadata)
                genes = ro.StrVector(data_filtered.index.tolist())

            # Clean Target column to valid R names
            target_col = r_metadata.rx2("Target")
            clean_target_col = base.make_names(target_col)
            target_idx = r_metadata.names.index("Target")
            r_metadata[target_idx] = clean_target_col

            # ------------------------------------------------------------------
            # 7. Fit linear model
            # ------------------------------------------------------------------
            r_formula = Formula("~0 + Target") if single_tissue else Formula("~0 + Target + tissue")
            r_design = stats.model_matrix(r_formula, data=r_metadata)

            print("    Fitting linear model...")
            fit = limma.lmFit(r_data, r_design)

            # ------------------------------------------------------------------
            # 8. Define and fit contrast
            # ------------------------------------------------------------------
            target_levels = base.levels(base.factor(r_metadata.rx2("Target")))
            if len(target_levels) < 2:
                print(f"    Not enough target levels for '{treatment}', skipping.")
                continue

            contrast_str = f"Target{target_levels[-1]}-Target{target_levels[0]}"
            print(f"    Contrast: {contrast_str}")
            contrast_matrix = limma.makeContrasts(contrast_str, levels=r_design)

            fit2 = limma.contrasts_fit(fit, contrast_matrix)
            fit2 = limma.eBayes(fit2)

            # ------------------------------------------------------------------
            # 9. Extract and save results
            # ------------------------------------------------------------------
            print("    Extracting results...")
            r_output = limma.topTreat(fit2, coef=1, genelist=genes, number=np.inf)
            writexl.write_xlsx(r_output, f"{out_dir}{output_filename}.xlsx")

            diff_exp_results = pd.read_excel(f"{out_dir}{output_filename}.xlsx")
            genes_to_save = (
                diff_exp_results.sort_values(by=["adj.P.Val"])[["ID", "t", "logFC", "adj.P.Val"]]
            )
            genes_to_save.to_csv(file_to_output, index=False)
            print(f"    Saved → {file_to_output}")

        except Exception as e:
            print(f"ERROR processing '{treatment}': {e}")
