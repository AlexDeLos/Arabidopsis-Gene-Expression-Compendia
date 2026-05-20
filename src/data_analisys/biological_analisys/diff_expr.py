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
from src.constants import CLUSTER_RUN, RNA_USED,DEBUG  # noqa: E402
from src.constants_labeling import TreatmentEnum  # noqa: E402
from src.data_analisys.utils.cluster_exploration_utils_final import get_gsm_id  # noqa: E402

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

    # ------------------------------------------------------------------
    # [0] One-time design diagnostics — printed before the treatment loop
    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("DIAGNOSTIC [0]: Label DataFrame (design)")
    print("=" * 60)
    print(f"  Total samples in design:        {len(design)}")
    print(f"  Design index name:              {design.index.name!r}")
    print(f"  Design index dtype:             {design.index.dtype}")
    print(f"  Design index sample (first 3):  {design.index[:3].tolist()}")
    print(f"  Columns present:                {design.columns.tolist()}")
    print(f"  Unique tissues ({len(design['tissue'].unique())}):  "
          f"{sorted(design['tissue'].dropna().unique())}")
    print("  Treatment value counts (top 15):")
    print(design["treatment"].value_counts().head(15).to_string())
    print("  Null counts per column:")
    print(design.isnull().sum().to_string())
    print("=" * 60 + "\n")

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
            data = pd.read_csv(os.path.join(save_dir, f"{data_type}.csv"), index_col=0)
            if DEBUG:
                # Keep only the samples we need + a random gene subset for speed
                debug_genes = data.index[:500]   # first 500 genes — deterministic, fast
                data = data.loc[debug_genes]
                print(f"  [DEBUG] Expression matrix capped to {data.shape}")
            if RNA_USED:
                data.columns = [get_gsm_id(col.split('_')[1]) for col in data.columns]
            print(f"  [1] Expression matrix shape:          {data.shape}")
            print(f"  [1] Expression col dtype:             {data.columns.dtype}")
            print(f"  [1] Expression col sample (first 3):  {list(data.columns[:3])}")
            print(f"  [1] Design index sample (first 3):    {design.index[:3].tolist()}")
            overlap = len(set(design.index) & set(data.columns))
            print(f"  [1] Index/col overlap:                {overlap}")
            # Show examples of mismatches — critical for catching format divergence
            in_design_not_data = set(design.index) - set(data.columns)
            in_data_not_design = set(data.columns) - set(design.index)
            print(f"  [1] In design but NOT in data:        "
                  f"{len(in_design_not_data)}  e.g. {sorted(in_design_not_data)[:3]}")
            print(f"  [1] In data but NOT in design:        "
                  f"{len(in_data_not_design)}  e.g. sorted({sorted(in_data_not_design)[:3]})")

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

            print(f"  [2] is_tissue True count:             {is_tissue.sum()}")
            print(f"  [2] is_treatment True count:          {is_treatment.sum()}")
            print(f"  [2] is_control True count:            {is_control.sum()}")
            if samples is not None:
                print(f"  [2] is_study True count:              {is_study.sum()}")
            print(f"  [2] design_filtered rows:             {len(design_filtered)}")
            print(f"  [2] Treatments in filtered:")
            print(design_filtered["treatment"].value_counts().to_string())
            # Catch samples that are labelled as BOTH treatment and control
            treatment_ids = set(design[is_treatment].index)
            control_ids   = set(design[is_control].index)
            print(f"  [2] Treatment-only samples:           {len(treatment_ids - control_ids)}")
            print(f"  [2] Control-only samples:             {len(control_ids - treatment_ids)}")
            print(f"  [2] Samples labelled BOTH:            {len(treatment_ids & control_ids)}")

            # ------------------------------------------------------------------
            # 3. Synchronise samples between expression data and design
            # ------------------------------------------------------------------
            common_samples = list(set(design_filtered.index) & set(data.columns))
            design_filtered = design_filtered[design_filtered.index.isin(common_samples)]
            data_filtered = data[common_samples]

            print(f"  [3] common_samples count:             {len(common_samples)}")
            # Show examples of why samples are lost — format mismatch shows up here
            not_in_data   = set(design_filtered.index) - set(data.columns)
            not_in_design = set(data.columns) - set(design_filtered.index)
            print(f"  [3] design_filtered samples not in data:   "
                  f"{len(not_in_data)}  e.g. {sorted(not_in_data)[:3]}")
            print(f"  [3] data columns not in design_filtered:   "
                  f"{len(not_in_design)}  e.g. {sorted(not_in_design)[:3]}")

            # Drop duplicate columns from data
            before_dedup = data_filtered.shape[1]
            # data_filtered = data_filtered.T.drop_duplicates().T
            data_filtered = data_filtered.loc[:, ~data_filtered.columns.duplicated(keep='first')]

            # Keep only design rows whose sample is still in data after dedup
            design_filtered = design_filtered[design_filtered.index.isin(data_filtered.columns)]

            # Drop duplicate rows from design
            design_filtered = design_filtered[~design_filtered.index.duplicated(keep='first')]

            print(f"  [3] Duplicate cols dropped from data: {before_dedup - data_filtered.shape[1]}")
            print(f"  [3] After dedup — data cols:          {data_filtered.shape[1]}")
            print(f"  [3] After dedup — design rows:        {len(design_filtered)}")

            # ------------------------------------------------------------------
            # 4. Filter by minimum group size
            # ------------------------------------------------------------------
            # Show ALL groups and their sizes BEFORE filtering so we can see
            # what is about to be dropped
            all_groups = design_filtered.groupby(["treatment", "tissue"]).size()
            print(f"  [4] All (treatment × tissue) groups BEFORE filter "
                  f"(threshold={filter_low_combination}):")
            print(all_groups.to_string())

            dropped_groups = all_groups[all_groups < filter_low_combination]
            if not dropped_groups.empty:
                print(f"  [4] Groups DROPPED (< {filter_low_combination} samples):")
                print(dropped_groups.to_string())
            else:
                print(f"  [4] No groups dropped — all meet the threshold.")

            design_filtered = (
                design_filtered.sort_values(by="Sample_ID")
                .reset_index(drop=False)          # index → "Sample_ID" column
                .groupby(["treatment", "tissue"])
                .filter(lambda x: len(x) >= filter_low_combination)
            )

            # Select expression columns using the Sample_ID column (not the index)
            data_filtered = data_filtered[design_filtered["Sample_ID"].values]

            print(f"  [4] After group filter — rows remaining: {len(design_filtered)}")
            if len(design_filtered) > 0:
                print("  [4] Groups surviving filter:")
                print(design_filtered.groupby(["treatment", "tissue"]).size().to_string())

            print(f"    {len(design_filtered)} samples aligned for "
                  f"'{treatment}' vs. {TreatmentEnum.CONTROL}.")

            # ------------------------------------------------------------------
            # 4b. Build metadata for the model
            # ------------------------------------------------------------------
            metadata = design_filtered[["treatment", "tissue"]].copy()
            metadata.index = design_filtered["Sample_ID"].values

            def r_make_names(s):
                s = re.sub(r'[^a-zA-Z0-9_]', '.', str(s))
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
                print(f"    Not enough unique Target levels for '{treatment}' "
                      f"(n={len(metadata)}), skipping.")
                print(f"    Target values present: {metadata['Target'].value_counts().to_dict()}")
                continue

            print(f"  [5] Final sample count:               {len(metadata)}")
            print(f"  [5] Unique Target values:             {sorted(metadata['Target'].unique())}")
            print(f"  [5] Target value counts:              {metadata['Target'].value_counts().to_dict()}")
            print(f"  [5] Unique tissues in model:          "
                  f"{sorted(metadata['tissue'].unique()) if 'tissue' in metadata.columns else 'single (dropped)'}")
            assert(len(metadata) == len(design_filtered))
            del design_filtered

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