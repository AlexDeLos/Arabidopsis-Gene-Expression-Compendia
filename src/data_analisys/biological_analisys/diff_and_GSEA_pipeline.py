"""
diff_and_GSEA_pipeline.py
-------------------------
Orchestration script for the full differential expression + GSEA pipeline.

Entry point:   run_diff_exp_and_enrichment(...)
Spider plots:  get_spider_plots(...)
"""

import os
import sys

import numpy as np
import pandas as pd
from difflib import SequenceMatcher
from nltk.stem import PorterStemmer

module_dir = "./"
sys.path.append(module_dir)

from src.constants import (  # noqa: E402
    CORE_DATA_DIR,
    EXPERIMENT_NAME,
    FIGURES_DIR,
    LABELS_PATH,
    PROCESSED_DATA_FOLDER,
    DEBUG
)
from src.constants_labeling import ( # noqa: E402
    TreatmentEnum
)

from src.data_analisys.biological_analisys.diff_expr import diff_exp_combine_tissues  # noqa: E402
from src.data_analisys.biological_analisys.plot_enrichment_new import (  # noqa: E402
    create_gsea_spider_plot,
    plot_enrichment_scatter_interactive,
)
from src.data_analisys.biological_analisys.pr_rank_gene_enrich import (  # noqa: E402
    get_go_data,
    perform_gsea_enrichment,
)
from src.data_analisys.utils.cluster_exploration_utils_final import (  # noqa: E402
    load_labels_study,
    make_df_from_labels,
)

# =============================================================================
# Constants
# =============================================================================
ITERATIONS = 1000

GO_OBO_FILE     = f"{CORE_DATA_DIR}go-basic.obo"
ANNOTATION_FILE = f"{CORE_DATA_DIR}tair.gaf.gz"

# ------------------------------------------------------------------
# Pass 1: load obodag unfiltered, just to discover GO root IDs
# ------------------------------------------------------------------
GO_NAME_OVERRIDES: dict[TreatmentEnum, str | None] = {
    # Fuzzy match fails — correct name provided manually
    TreatmentEnum.DEHYDRATION:         "response to water deprivation",
    TreatmentEnum.BIOTIC:              "response to biotic stimulus",
    TreatmentEnum.ABIOTIC:             "response to abiotic stimulus",
    TreatmentEnum.LOW_LIGHT:           "response to light intensity",
    TreatmentEnum.OTHER_LIGHT:         "response to light stimulus",
    TreatmentEnum.CUT:                 "response to wounding",
    TreatmentEnum.NUTRIENT_DEFICIENCY: "response to nutrient levels",
    
    # No meaningful GO root — exclude from GSEA entirely
    TreatmentEnum.OTHER:               None,
    TreatmentEnum.CONTROL:             None,
    TreatmentEnum.UNKNOWN:             None,
    TreatmentEnum.UNSPECIFIED:         None,
}

_stemmer = PorterStemmer()

def _stem_phrase(phrase: str) -> str:
    """Stem every word in a phrase and rejoin."""
    return " ".join(_stemmer.stem(w) for w in phrase.lower().split())

def find_go_root_for_treatment(treatment: TreatmentEnum, obodag) -> str | None:
    # Check for explicit override first
    if treatment in GO_NAME_OVERRIDES:
        override = GO_NAME_OVERRIDES[treatment]
        if override is None:
            return None  # explicitly excluded
        query = override
    else:
        query = f"response to {treatment.value.lower()}"

    stemmed_query = _stem_phrase(query)
    candidates = [
        (go_id, term) for go_id, term in obodag.items()
        if term.name.lower().startswith("response to")
        and term.namespace == "biological_process"
    ]

    if not candidates:
        return None

    def score(item):
        go_id, term = item
        stemmed_name = _stem_phrase(term.name)
        name_similarity = SequenceMatcher(
            None, stemmed_query, stemmed_name
        ).ratio()
        depth_penalty = len(term.parents) * 0.01
        return name_similarity - depth_penalty

    best_go_id, best_term = max(candidates, key=score)
    best_score = score((best_go_id, best_term))

    print(f"  {treatment.value!r} → '{best_term.name}' ({best_go_id}, score={best_score:.3f})")

    if best_score < 0.6:
        print(f"  Warning: low confidence match for '{treatment.value}', consider adding an override.")
        return None

    return best_go_id


unfiltered_obodag, _ = get_go_data(GO_OBO_FILE, ANNOTATION_FILE)

TREATMENT_GO_MAP: dict[TreatmentEnum, str] = {}
for treatment in TreatmentEnum:
    go_id = find_go_root_for_treatment(treatment, unfiltered_obodag)
    if go_id:
        TREATMENT_GO_MAP[treatment] = go_id

# Derived automatically
TREATMENTS:      list[str]      = [t.value for t in TREATMENT_GO_MAP]
print(TREATMENTS)
STRESS_GO_ROOTS: dict[str, str] = {go: t.value for t, go in TREATMENT_GO_MAP.items()}
STRESS_IDS:      set[str]       = set(STRESS_GO_ROOTS.keys())

# ------------------------------------------------------------------
# Pass 2: reload with proper descendant filtering now that STRESS_IDS is known
# ------------------------------------------------------------------
obodag, geneid2gos = get_go_data(
    GO_OBO_FILE,
    ANNOTATION_FILE,
    stress_root_go_ids=STRESS_IDS,
)
#: Sample IDs used for the single-study sanity check (subset of GSE36649).
SANITY_CHECK_SAMPLES = [
    "GSM1027688", "GSM1027720", "GSM1027701", "GSM1027875", "GSM1027745",
    "GSM1027710", "GSM1027860", "GSM1027768", "GSM1027781", "GSM1027855",
    "GSM1027682", "GSM1027801", "GSM1027843", "GSM1027746", "GSM1027854",
    "GSM1027695", "GSM1027690", "GSM1027804", "GSM1027730", "GSM1027719",
    "GSM1027773", "GSM1027871", "GSM1027691", "GSM1027856", "GSM1027738",
    "GSM1027742", "GSM1027800", "GSM1027806", "GSM1027815", "GSM1027863",
    "GSM1027776", "GSM1027760", "GSM1027882", "GSM1027817", "GSM1027876",
    "GSM1027739", "GSM1027721", "GSM1027861", "GSM1027704", "GSM1027826",
    "GSM1027799", "GSM1027808", "GSM1027759", "GSM1027818", "GSM1027830",
    "GSM1027733", "GSM1027751", "GSM1027715", "GSM1027697", "GSM1027728",
    "GSM1027785", "GSM1027824", "GSM1027835", "GSM1027740", "GSM1027767",
    "GSM1027884", "GSM1027716", "GSM1027868", "GSM1027880", "GSM1027849",
    "GSM1027698", "GSM1027807", "GSM1027689", "GSM1027812", "GSM1027680",
    "GSM1027709", "GSM1027862", "GSM1027877", "GSM1027832", "GSM1027726",
    "GSM1027839", "GSM1027753", "GSM1027718", "GSM1027810", "GSM1027797",
    "GSM1027834", "GSM1027761", "GSM1027802", "GSM1027732", "GSM1027821",
    "GSM1027793", "GSM1027684", "GSM1027853", "GSM1027712", "GSM1027703",
    "GSM1027789", "GSM1027692", "GSM1027722", "GSM1027833", "GSM1027827",
    "GSM1027699", "GSM1027873", "GSM1027735", "GSM1027774", "GSM1027841",
    "GSM1027820", "GSM1027743", "GSM1027696", "GSM1027837", "GSM1027846",
    "GSM1027744", "GSM1027694", "GSM1027850", "GSM1027828", "GSM1027805",
    "GSM1027847", "GSM1027794", "GSM1027829", "GSM1027811", "GSM1027845",
    "GSM1027842", "GSM1027754", "GSM1027678", "GSM1027859", "GSM1027787",
    "GSM1027708", "GSM1027852", "GSM1027792", "GSM1027874", "GSM1027878",
    "GSM1027883", "GSM1027700", "GSM1027848", "GSM1027679", "GSM1027879",
    "GSM1027851", "GSM1027756", "GSM1027724", "GSM1027764", "GSM1027723",
    "GSM1027750", "GSM1027798", "GSM1027775", "GSM1027791", "GSM1027795",
    "GSM1027831", "GSM1027731", "GSM1027685", "GSM1027870", "GSM1027736",
    "GSM1027836", "GSM1027872", "GSM1027734", "GSM1027866", "GSM1027725",
    "GSM1027749", "GSM1027757", "GSM1027687", "GSM1027788", "GSM1027681",
    "GSM1027858", "GSM1027881", "GSM1027840", "GSM1027786", "GSM1027763",
    "GSM1027705", "GSM1027778", "GSM1027765", "GSM1027867", "GSM1027825",
    "GSM1027752", "GSM1027803", "GSM1027790", "GSM1027823", "GSM1027809",
    "GSM1027822", "GSM1027783", "GSM1027771", "GSM1027755", "GSM1027747",
    "GSM1027814", "GSM1027686", "GSM1027844", "GSM1027741", "GSM1027706",
    "GSM1027762", "GSM1027780", "GSM1027702", "GSM1027777", "GSM1027714",
    "GSM1027782", "GSM1027717", "GSM1027865", "GSM1027713", "GSM1027758",
    "GSM1027838", "GSM1027796", "GSM1027769", "GSM1027857", "GSM1027770",
    "GSM1027707", "GSM1027748", "GSM1027813", "GSM1027784", "GSM1027766",
    "GSM1027869", "GSM1027772", "GSM1027779", "GSM1027816", "GSM1027727",
    "GSM1027693", "GSM1027737", "GSM1027864", "GSM1027683", "GSM1027819",
    "GSM1027729", "GSM1027711",
]


# =============================================================================
# Spider plot helper
# =============================================================================

def get_spider_plots(
    path: str,
    results_path: str,
    data_types: list[str],
    Fulls: list[bool],
    tissues: list[str | None],
    pure_val: bool,
    filter_val: int
) -> None:
    """
    Generate radar / spider plots that compare GSEA statistics for each GO
    root term across all valid experiment configurations.

    One SVG per GO root term is written to `path`.

    Parameters
    ----------
    path : str
        Output directory for spider plot SVGs.
    results_path : str
        Root directory that contains the per-experiment GSEA result folders.
    data_types : list of str
        All normalization stages to consider.
    Fulls : list of bool
        Which study-scope modes to consider (True = all studies, False = sanity).
    tissues : list of str or None
        Which tissue filters to consider (None = all tissues).
    pure_val : bool
        Whether the run used pure (single-treatment) samples.
    filter_val : int
        Minimum-group-size filter value used for this run.
    full : bool
        If True, include all `data_types`; otherwise restrict to a curated
        subset (combat_seq, tissue_normalized*, imputed).
    """
    os.makedirs(path, exist_ok=True)

    allowed_types = data_types
    pure_str = "pure" if pure_val else "mixed"

    for term, stress in STRESS_GO_ROOTS.items():
        print(f"  Collecting spider data for {term} ({stress})...")
        all_rows: list[pd.DataFrame] = []

        for data_type in allowed_types:
            for Full in Fulls:
                for tissue in tissues:
                    tissue_str  = tissue if tissue else "All_tissues"
                    full_str    = "full" if Full else "sanity"
                    exp_name    = (
                        f"{EXPERIMENT_NAME}_{data_type}_{tissue_str}_"
                        f"{full_str}_{pure_str}_min_group_{filter_val}"
                    )
                    out_dir = (
                        f"{results_path}GSEA_enrichment_{exp_name}/"
                    )
                    csv_file    = (
                        f"{out_dir}{stress}_gsea_go_enrichment_results_{ITERATIONS}.csv"
                    )

                    if not os.path.isfile(csv_file):
                        continue

                    try:
                        gsea_df  = pd.read_csv(csv_file)
                        
                        term_row = gsea_df[gsea_df["go_id"] == term].copy()
                        if term_row.empty:
                            continue

                        tissue_desc = "Single Tissue" if tissue else "Combined Tissues"
                        study_desc  = "Combined Studies" if Full else "Single Study"
                        term_row["Name"] = f"{data_type} {tissue_desc} {study_desc}"
                        all_rows.append(term_row)

                    except Exception as exc:
                        print(f"    Warning: could not process {csv_file}: {exc}")

        if not all_rows:
            print(f"  No data for {term}, skipping spider plot.")
            continue

        df = pd.concat(all_rows, ignore_index=True)
        if df.empty:
            continue

        plot_filename = (
            f"{path}{stress.replace(' ', '_')}_spider_plot.svg"
        )
        try:
            create_gsea_spider_plot(df, save_path=plot_filename, term=stress)
        except Exception as exc:
            print(f"  Error generating spider plot for {term}: {exc}")


# =============================================================================
# Main pipeline
# =============================================================================

def run_diff_exp_and_enrichment(
    save_dir: str = PROCESSED_DATA_FOLDER,
    data_types: list[str] | None = None,
    pures: list[bool] | None = None,
    Fulls: list[bool] | None = None,
    filter_low_combination: list[int] | None = None,
    tissues: list[str | None] | None = None,
    just_plot: bool = False,
) -> None:
    """
    Run differential expression and GSEA for all combinations of the supplied
    experiment parameters, then generate per-combination scatter plots and
    cross-combination spider plots.

    Skip steps that have already been completed (crash-safe / resumable):
    - Differential expression: skipped if ``{out_dir}done.txt`` exists.
    - GSEA: skipped if the results CSV already exists.

    Parameters
    ----------
    save_dir : str
        Directory containing ``{data_type}.csv`` expression matrices.
    data_types : list of str
        Normalisation stages to iterate over.
    pures : list of bool
        Whether to use pure (single-treatment) samples only.
    Fulls : list of bool
        True = all studies, False = sanity-check samples only.
    filter_low_combination : list of int
        Minimum-group-size thresholds to iterate over.
    tissues : list of str or None
        Tissue filters to iterate over (None = all tissues combined).
    just_plot : bool
        If True, skip computation and regenerate plots from existing CSVs only.
    """
    if data_types is None:
        data_types = ["combat_seq_norm", "rankin", "filter"]
    if pures is None:
        pures = [False]
    if Fulls is None:
        Fulls = [True, False]
    if filter_low_combination is None:
        filter_low_combination = [0, 15]
    if tissues is None:
        tissues = [None]

    labels = make_df_from_labels(load_labels_study(LABELS_PATH))

    if DEBUG:
        print("\n  *** DEBUG MODE — using random subset of data ***\n")
        labels = subsample_labels_for_debug(
            labels,
            n_per_group=DEBUG_N_SAMPLES_PER_GROUP,
            treatments=TREATMENTS,
            seed=DEBUG_SEED,
        )

    print("\n" + "="*60)
    print("DIAGNOSTIC: Label DataFrame")
    print("="*60)
    print(f"  Total samples in design:       {len(labels)}")
    print(f"  Design index name:             {labels.index.name}")
    print(f"  Design index sample (first 3): {labels.index[:3].tolist()}")
    print(f"  Columns:                       {labels.columns.tolist()}")
    print(f"  Unique tissues:                {sorted(labels['tissue'].unique())}")
    print("\n  Treatment value counts (top 20):")
    print(labels['treatment'].value_counts().head(20).to_string())
    print("\n  Null counts per column:")
    print(labels.isnull().sum().to_string())
    print("="*60 + "\n")

    for fil in filter_low_combination:
        for pure in pures:
            pure_str = "pure" if pure else "mixed"

            for data_type in data_types:
                for Full in Fulls:
                    for tissue in tissues:
                        tissue_str     = tissue if tissue else "All_tissues"
                        tissue_display = tissue if tissue is not None else "All-Tissues"
                        full_str       = "full" if Full else "sanity"

                        exp_name        = (
                            f"{EXPERIMENT_NAME}_{data_type}_{tissue_str}_"
                            f"{full_str}_{pure_str}_min_group_{fil}"
                        )
                        diff_exp_outdir = f"{FIGURES_DIR}dif_expression_results/{exp_name}/"
                        done_file       = f"{diff_exp_outdir}done.txt"

                        # ------------------------------------------------------
                        # Differential expression
                        # ------------------------------------------------------
                        if not os.path.isfile(done_file):
                            if not just_plot:
                                target_samples = None if Full else SANITY_CHECK_SAMPLES
                                diff_exp_combine_tissues(
                                    TREATMENTS,
                                    save_dir,
                                    data_type,
                                    out_dir=diff_exp_outdir,
                                    design=labels,
                                    samples=target_samples,
                                    pure=pure,
                                    tissue=tissue,
                                    filter_low_combination=fil,
                                )
                            produced = any(
                                f.endswith("_genes.csv")
                                for f in os.listdir(diff_exp_outdir)
                                if os.path.isfile(os.path.join(diff_exp_outdir, f))
                            )
                            if produced:
                                with open(done_file, "w") as fh:
                                    fh.write("completed")
                                print(f"  Done → {diff_exp_outdir}")
                            else:
                                print(f"  Warning: no gene CSVs produced for {diff_exp_outdir}, not marking done.")
                        else:
                            print(f"  Already done: {diff_exp_outdir}")

                        # ------------------------------------------------------
                        # GSEA per stress treatment
                        # ------------------------------------------------------
                        gsea_outdir = (
                            f"{FIGURES_DIR}GSEA_enrichment_results/"
                            f"GSEA_enrichment_{exp_name}/"
                        )

                        for stress in TREATMENTS:
                            try:
                                gsea_csv = (
                                    f"{gsea_outdir}{stress}_gsea_go_enrichment_results"
                                    f"_{ITERATIONS}.csv"
                                )

                                if (not just_plot) and (not os.path.isfile(gsea_csv)):
                                    # Keep only the five root term sub-DAGs
                                    diff_csv = (
                                        f"{diff_exp_outdir}"
                                        f"{tissue_display}_{stress}_genes.csv"
                                    )
                                    diff_results = pd.read_csv(diff_csv)
                                    diff_results["rank"] = (
                                        diff_results["logFC"]
                                        * (-np.log10(diff_results["adj.P.Val"]))
                                    )

                                    gsea_df = perform_gsea_enrichment(
                                        ranked_gene_df=diff_results,
                                        gene_col="ID",
                                        rank_col="rank",
                                        obodag=obodag,
                                        geneid2gos=geneid2gos,
                                        keys=None,
                                        stress=stress,
                                        out_path=gsea_outdir,
                                        permutations=ITERATIONS,
                                    )
                                    gsea_df.sort_values(by="FDR q-val", inplace=True)
                                    gsea_df.to_csv(gsea_csv, index=False)
                                    print(f"    GSEA saved → {gsea_csv}")

                                else:
                                    print(f"    Loading pre-existing GSEA: {gsea_csv}")
                                    gsea_df = pd.read_csv(gsea_csv)

                                # Scatter plot
                                plot_title = (
                                    f"GSEA for {stress} | {tissue_display} | "
                                    f"{full_str} | {pure_str} treatments | "
                                    f"filter >{fil}"
                                )
                                plot_out = (
                                    f"{FIGURES_DIR}plots_enrichment/"
                                    f"{EXPERIMENT_NAME}/{full_str}/{tissue_display}/"
                                    f"{data_type}/{fil}/{pure_str}/{stress}.html"
                                )
                                plot_enrichment_scatter_interactive(
                                    gsea_df,
                                    save_path=plot_out,
                                    title=plot_title,
                                    treatments=TREATMENTS,
                                    normalizations=data_types,
                                )

                            except Exception as exc:
                                print(f"  Error during GSEA for {stress}: {exc}")

            # ------------------------------------------------------------------
            # Spider plots — one set per (filter × pure) combination
            # ------------------------------------------------------------------
            spider_dir = (
                f"{FIGURES_DIR}GSEA_radar_comparison/"
                f"{EXPERIMENT_NAME}_{fil}_{pure_str}/"
            )
            get_spider_plots(
                path=spider_dir,
                results_path=f"{FIGURES_DIR}GSEA_enrichment_results/",
                data_types=data_types,
                Fulls=Fulls,
                tissues=tissues,
                pure_val=pure,
                filter_val=fil
            )

    print("Pipeline complete.")


# =============================================================================
# Debug flag — set to True to run locally on a small random subset
# =============================================================================
DEBUG_SEED  = None  # Set to None for a completely different random subset every run!
DEBUG_N_SAMPLES_PER_GROUP = 10  # samples per (treatment × tissue) group to keep

def subsample_labels_for_debug(
    labels: pd.DataFrame,
    n_per_group: int,
    treatments: list[str],
    seed: int | None = None,
) -> pd.DataFrame:
    """
    Subsamples the labels dataframe for local debugging.
    If seed is None, it uses OS entropy to guarantee a unique set of 
    samples every single time you execute the script.
    """
    # Initialize Generator. None means truly random every time.
    rng = np.random.default_rng(seed)

    # Keep only relevant treatments + Control
    is_active = labels["treatment"].isin(treatments + ["Control"])
    subset = labels[is_active].copy()

    sampled_indices = []
    
    # Iterate through each group and sample safely from their true indices
    for position, grp in subset.groupby(["treatment", "tissue"]):
        n = min(len(grp), n_per_group)
        # Sample directly from the index labels (preserves Sample_IDs)
        chosen_samples = rng.choice(grp.index, size=n, replace=False)
        sampled_indices.extend(chosen_samples)
        
    # Extract the final dataframe matching our randomly chosen samples
    final_subset = labels.loc[sampled_indices].copy()

    print(f"  [DEBUG] Labels subsampled: {len(labels)} → {len(final_subset)} samples")
    print("  [DEBUG] Treatment counts in subset:")
    print(final_subset["treatment"].value_counts().to_string())
    print("  [DEBUG] Tissue counts in subset:")
    print(final_subset["tissue"].value_counts().to_string())
    
    return final_subset

# =============================================================================
# CLI entry point
# =============================================================================

if __name__ == "__main__":
    run_diff_exp_and_enrichment(just_plot=False, data_types=['filter_norm', 'combat_seq_norm'],Fulls=[True], filter_low_combination=[0])
