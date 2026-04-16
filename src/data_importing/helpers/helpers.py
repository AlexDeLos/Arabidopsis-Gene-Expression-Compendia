import json
import logging
import os
import sys
from collections.abc import Callable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
import seaborn as sns

module_dir = "./"
sys.path.append(module_dir)
from src.constants import GEO_DOWNLOAD_DIR, METADATA_OUTPUT_DIR, PROCESSED_DATA_FOLDER  # noqa: E402


def setup_directories():
    """Create output directories if they don't exist."""
    logging.info("Setting up output directories...")
    os.makedirs(GEO_DOWNLOAD_DIR, exist_ok=True)
    os.makedirs(METADATA_OUTPUT_DIR, exist_ok=True)
    os.makedirs(PROCESSED_DATA_FOLDER, exist_ok=True)


def process_metadata(geo_accession, gse, gsm, save_path=METADATA_OUTPUT_DIR):
    """Filters and saves metadata for a given sample to the ./metadata folder."""
    excluded_keys = ["contact", "date", "_id", "taxid", "data_", "status", "type", "contributor", "relation", "geo_accession", "row_count", "organism", "label", "supplementary_file"]
    filtered_study_meta = {k: v for k, v in gse.metadata.items() if not any(ex in k for ex in excluded_keys)}
    filtered_sample_meta = {k: v for k, v in gsm.metadata.items() if not any(ex in k for ex in excluded_keys)}
    final_metadata = {"study_id": geo_accession, "sample_id": gsm.name, "platform": gse.metadata["platform_id"][0], "study_metadata": filtered_study_meta, "sample_metadata": filtered_sample_meta}
    output_path = os.path.join(save_path, f"{geo_accession}_{gsm.name}.json")
    os.makedirs(save_path, exist_ok=True)
    try:
        with open(output_path, "w") as fp:
            json.dump(final_metadata, fp, indent=4)
    except Exception as e:
        logging.error(f"Failed to save metadata for {gsm.name}: {e}")


def process_sample_data(geo_accession, gsm, probe_to_gene_map):
    """
    Processes expression data for a sample and returns it as a DataFrame
    with a single, uniquely named column.
    """
    if gsm.table.empty:
        logging.warning(f"Sample {gsm.name} from study {geo_accession} has no data table. Skipping.")
        return None

    df = gsm.table.copy()
    df["gene_symbol"] = df["ID_REF"].map(probe_to_gene_map)
    df.dropna(subset=["gene_symbol"], inplace=True)
    df["gene_symbol"] = df["gene_symbol"].map(lambda x: x.split("_")[0])
    if df.empty:
        logging.warning(f"No genes could be mapped for sample {gsm.name}. Trying a diferent mapping.")
        # TODO: a lot of the entries here already have proper mappings, so just let them through
        df = gsm.table.copy()
        df["gene_symbol"] = df["ID_REF"].map(lambda x: x.split("_")[0])
        # del df['ID_REF']
        # return None

    # Filter for Arabidopsis thaliana genes (AtXgXXXXX format)
    at_gene_regex = r"^At[1-5MC]g\d{5}$"
    df = df[df["gene_symbol"].str.match(at_gene_regex, case=False)]
    if df.empty:
        logging.error(f"No Arabidopsis thaliana genes found in {gsm.name}. Skipping.")
        return None

    df.set_index("gene_symbol", inplace=True)
    df = df[["VALUE"]]

    # Handle duplicate genes by averaging their expression
    df = df.groupby(df.index).mean()

    # Create a unique name for the sample's data column
    unique_sample_id = f"{geo_accession}_{gsm.name}"
    return df.rename(columns={"VALUE": unique_sample_id})


def handle_duplicates(sample_df):

    if sample_df.index.duplicated().any():
        # print('Duplicates in df:', in_df[in_df.index.duplicated(keep=False)])

        # Create a dictionary to keep track of the count of duplicates
        # duplicate_count = sample_df.index.value_counts().to_dict()

        # Group by the index column and calculate the mean of the values
        sample_df = sample_df.groupby("gene_symbol")[sample_df.columns[0]].mean().reset_index()  # assuming there is only one value
    return sample_df


def get_geo_list(path: str):
    read = pd.read_csv(path)
    read = read.loc[read["depository_source"] == "GEO"]
    read = read.loc[read["species"] == "Arabidopsis thaliana"]
    return list(read["depository_accession"])


def mapping(x):
    if type(x) is str:
        return x.upper()
    return x


def predicate(gene: str, chromosome: str) -> bool:
    return str("AT" + chromosome + "G") in gene


def get_first_indexs(df_index, chromo: list[str]):
    array = []
    for i in chromo:
        gene: str = next(filter(lambda x: predicate(x, str(i)), df_index))
        array.append(df_index.get_loc(gene))
    return array


def plot_platform_intensities(df: pl.DataFrame, save_dir, key, lengend=False):
    """
    Function 2: For each platform (by id), plots an intensity plot of the data.
    Uses a different color for each sample with transparency.

    Parameters:
    - df: A Polars DataFrame containing the data.
    - bucket_size: The bin width for the histogram (intensity bucket).
    """
    os.makedirs(save_dir, exist_ok=True)
    # Get unique platforms
    platforms = df["id"].unique().to_list()

    for platform_id in platforms:
        # platform_id ='GPL198'
        # 1. Filter data for this specific platform using Polars syntax
        platform_subset = df.filter(pl.col("id") == platform_id)

        # 2. Explode the 'data' column using Polars syntax.
        # This expands the list[float] column into individual rows.
        exploded_data = platform_subset.explode("data")

        # 3. Plotting
        # Create a dictionary of numpy arrays for Seaborn
        # This avoids converting the entire object to a Pandas DataFrame
        data_array = np.nan_to_num(exploded_data["data"].to_numpy())
        plot_data = {"data": data_array, "sample_id": exploded_data["sample_id"].to_numpy(), "study_id": exploded_data["study_id"].to_numpy()}
        bucket_size = data_array.max() / 300

        plt.figure(figsize=(12, 7))

        sns.histplot(
            data=plot_data,
            x="data",
            hue=key,  # Different color for each sample
            binwidth=bucket_size,  # The requested bucket size
            element="step",  # 'step' style is cleaner for overlapping distributions
            fill=True,  # Fill the area under the curve/step
            alpha=0.1,  # Transparency to see overlaps
            # kde=False              # Disable KDE
            legend=lengend,
        )
        df.remove(pl.col("id") == platform_id)
        plt.title(f"Intensity Distribution: {platform_id}")
        plt.xlabel("Intensity Value")
        plt.ylabel("Frequency")
        plt.grid(True, linestyle=":", alpha=0.6)
        plt.tight_layout()
        plt.savefig(save_dir + f"platform_intensity_{platform_id}.svg")
        plt.close()
        # return


class tracker:
    def __init__(self, storage) -> None:
        self.schema_dic: dict = {"id": str, "sample_id": str, "study_id": str, "data": list[float]}
        self.platform_counter: pl.DataFrame = pl.DataFrame(schema=self.schema_dic)  # .set_index('id')
        os.makedirs(storage, exist_ok=True)
        self.storage_loc: str = storage + "tracker.json"

    def update_counter(self, sample_id, sample_info, study_id, data) -> None:
        # new_data = np.array(data)
        new_record: pl.DataFrame = pl.DataFrame([{"id": sample_info.name, "sample_id": sample_id, "study_id": study_id, "data": data}], schema=self.schema_dic)  # .set_index('id')
        new_df = pl.concat([self.platform_counter, new_record])
        self.platform_counter = new_df

    def store(self):
        self.platform_counter.write_json(self.storage_loc)

    def load(self):
        self.platform_counter = pl.read_json(self.storage_loc).match_to_schema(self.schema_dic)

    def log_transform_high_intensity_studies(self):
        df = self.platform_counter
        df.explode("data")["data"].to_numpy().max()
        self.platform_counter = df.with_columns(
            pl.when(pl.col("data").list.max().max().over("study_id") > 300)
            .then(
                # Apply log2 to every element in the list using list.eval
                pl.col("data").list.eval((pl.element() + 1).log(2))
            )
            .otherwise(pl.col("data"))
        )

    def plot(self):
        df = self.platform_counter.clone()

        # plot_platform_usage(df,'test_plots/')
        # plot_platform_intensities(df,'test_plots/intensity_plots_by_sample/','sample_id')
        plot_platform_intensities(df, "test_plots/intensity_plots_by_study/", "study_id", False)


def plot_tracker_results(json_file="tracker_results.json", output_dir="."):
    """
    Reads tracker JSON and saves SVG plots:
    1. Pie Chart: Total Studies (Used vs Skipped)
    2. Pie Chart: Total Samples (Used vs Skipped)
    3. Bar Charts: Platforms (Filtered by >0 samples used)
    """
    # 1. Load Data
    try:
        with open(json_file) as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f"Error: {json_file} not found.")
        return

    totals = data["totals"]
    platform_data = data["platform_counts"]

    # --- PLOT 1: PIE CHART (STUDIES) ---
    studies_used = totals["total_studies_used"]
    studies_seen = totals["total_studies_seen"]
    studies_skipped = studies_seen - studies_used

    plt.figure(figsize=(6, 6))
    plt.pie([studies_used, studies_skipped], labels=[f"Used ({studies_used})", f"Skipped ({studies_skipped})"], autopct="%1.1f%%", colors=sns.color_palette("pastel")[0:2], startangle=140)
    plt.title(f"Total Studies Seen: {studies_seen}")
    plt.savefig(f"{output_dir}/tracker_pie_studies.svg", format="svg")
    plt.close()
    print(f"Saved {output_dir}/tracker_pie_studies.svg")

    # --- PLOT 2: PIE CHART (SAMPLES) ---
    samples_used = totals["total_samples_used"]
    samples_seen = totals["total_sample_seen"]
    samples_skipped = samples_seen - samples_used

    plt.figure(figsize=(6, 6))
    plt.pie([samples_used, samples_skipped], labels=[f"Used ({samples_used})", f"Skipped ({samples_skipped})"], autopct="%1.1f%%", colors=sns.color_palette("pastel")[2:4], startangle=140)
    plt.title(f"Total Samples Seen: {samples_seen}")
    plt.savefig(f"{output_dir}/tracker_pie_samples.svg", format="svg")
    plt.close()
    print(f"Saved {output_dir}/tracker_pie_samples.svg")

    # --- PLOT 3: BAR CHARTS (PLATFORMS) ---
    # Filter and Reshape Data
    records = []

    for platform, stats in platform_data.items():
        # CRITICAL FILTER: Only include platforms where we actually used samples
        if stats["samples_used"] > 0:
            # Add Seen Data
            records.append({"Platform": platform, "Count": stats["studies_seen"], "Metric": "Studies", "Status": "Seen"})
            records.append({"Platform": platform, "Count": stats["samples_seen"], "Metric": "Samples", "Status": "Seen"})

            # Add Used Data
            records.append({"Platform": platform, "Count": stats["studies_used"], "Metric": "Studies", "Status": "Used"})
            records.append({"Platform": platform, "Count": stats["samples_used"], "Metric": "Samples", "Status": "Used"})

    if not records:
        print("No platforms found with used samples. Skipping bar plots.")
        return

    df = pd.DataFrame(records)

    # Sort by Used Samples to make the plot readable
    # We find the order of platforms based on 'Samples' + 'Used' count
    sorter = df[(df["Metric"] == "Samples") & (df["Status"] == "Used")].sort_values("Count", ascending=False)
    order = sorter["Platform"].tolist()

    # Create Subplots (1 Row, 2 Columns)
    sns.set_theme(style="whitegrid")
    _fig, axes = plt.subplots(1, 2, figsize=(16, 8), sharey=True)

    # Subplot A: Studies
    studies_df = df[df["Metric"] == "Studies"]
    sns.barplot(data=studies_df, x="Count", y="Platform", hue="Status", order=order, ax=axes[0], palette="muted")
    axes[0].set_title("Studies per Platform (Filtered)")
    axes[0].set_xlabel("Number of Studies")

    # Subplot B: Samples
    samples_df = df[df["Metric"] == "Samples"]
    sns.barplot(data=samples_df, x="Count", y="Platform", hue="Status", order=order, ax=axes[1], palette="muted")
    axes[1].set_title("Samples per Platform (Filtered)")
    axes[1].set_xlabel("Number of Samples")

    plt.tight_layout()
    plt.savefig(f"{output_dir}/tracker_platforms.svg", format="svg")
    plt.close()
    print(f"Saved {output_dir}/tracker_platforms.svg")


def plot_tracker_results_RNA(json_file="rnaseq_tracker_results.json", output_dir="."):
    """
    Reads tracker JSON and saves SVG plots:
    1. Pie Chart: Total Studies (Used vs Skipped)
    2. Pie Chart: Total Samples (Used vs Skipped)
    3. Bar Charts: Platforms (Filtered by >0 samples used)
    4. Histogram: Distribution of Samples per Study
    """
    # 1. Load Data
    try:
        with open(json_file) as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f"Error: {json_file} not found.")
        return

    totals = data.get("totals", {})
    platform_data = data.get("platform_counts", {})
    # sample_dist = data.get('sample_distribution', {}) # Format: {"12": 5, "24": 1}
    dist_per_platform = data.get("sample_distribution_per_platform", {})
    # Use a style suitable for publication
    sns.set_theme(style="whitegrid")

    # --- PLOT 1: PIE CHART (STUDIES) ---
    studies_used = totals.get("total_studies_used", 0)
    studies_seen = totals.get("total_studies_seen", 0)
    studies_skipped = studies_seen - studies_used

    if studies_seen > 0:
        plt.figure(figsize=(6, 6))
        plt.pie([studies_used, studies_skipped], labels=[f"Used ({studies_used})", f"Skipped ({studies_skipped})"], autopct="%1.1f%%", colors=sns.color_palette("pastel")[0:2], startangle=140)
        plt.title(f"Total Studies Seen: {studies_seen}")
        plt.savefig(f"{output_dir}/tracker_pie_studies.svg", format="svg")
        plt.close()
        print(f"Saved {output_dir}/tracker_pie_studies.svg")

    # --- PLOT 2: PIE CHART (SAMPLES) ---
    samples_used = totals.get("total_samples_used", 0)
    samples_seen = totals.get("total_sample_seen", 0)
    samples_skipped = samples_seen - samples_used

    if samples_seen > 0:
        plt.figure(figsize=(6, 6))
        plt.pie([samples_used, samples_skipped], labels=[f"Used ({samples_used})", f"Skipped ({samples_skipped})"], autopct="%1.1f%%", colors=sns.color_palette("pastel")[2:4], startangle=140)
        plt.title(f"Total Samples Seen: {samples_seen}")
        plt.savefig(f"{output_dir}/tracker_pie_samples.svg", format="svg")
        plt.close()
        print(f"Saved {output_dir}/tracker_pie_samples.svg")

    # --- PLOT 3: BAR CHARTS (PLATFORMS) ---
    records = []
    for platform, stats in platform_data.items():
        if stats["samples_with_raw"] > 0:
            records.append({"Platform": platform, "Count": stats["studies_seen"], "Metric": "Studies", "Status": "Seen"})
            records.append({"Platform": platform, "Count": stats["samples_seen"], "Metric": "Samples", "Status": "Seen"})
            records.append({"Platform": platform, "Count": stats["studies_with_raw"], "Metric": "Studies", "Status": "Used"})
            records.append({"Platform": platform, "Count": stats["samples_with_raw"], "Metric": "Samples", "Status": "Used"})

    if records:
        df = pd.DataFrame(records)
        sorter = df[(df["Metric"] == "Samples") & (df["Status"] == "Used")].sort_values("Count", ascending=False)
        order = sorter["Platform"].tolist()

        _fig, axes = plt.subplots(1, 2, figsize=(16, 8), sharey=True)  # Share Y axis (Platforms)

        sns.barplot(data=df[df["Metric"] == "Studies"], x="Count", y="Platform", hue="Status", order=order, ax=axes[0], palette="muted")
        axes[0].set_title("Studies per Platform")

        sns.barplot(data=df[df["Metric"] == "Samples"], x="Count", y="Platform", hue="Status", order=order, ax=axes[1], palette="muted")
        axes[1].set_title("Samples per Platform")

        plt.tight_layout()
        plt.savefig(f"{output_dir}/tracker_platforms.svg", format="svg")
        plt.close()
        print(f"Saved {output_dir}/tracker_platforms.svg")
    else:
        print("No platform data available to plot.")
    if dist_per_platform:
        # 1. Identify Top 6 Platforms by Total Studies Seen
        # Sorting format: (PlatformName, {stats...})
        sorted_platforms = sorted(platform_data.items(), key=lambda item: item[1]["studies_seen"], reverse=True)
        top_6_platforms = [p[0] for p in sorted_platforms[:6]]

        # 2. Reconstruct Data for DataFrame
        plot_records = []
        for platform in top_6_platforms:
            # Get the frequency dict for this platform: {"12": 5, "24": 1}
            freq_dict = dist_per_platform.get(platform, {})

            for size_str, count in freq_dict.items():
                size = int(size_str)
                # Add 'count' rows for this size
                plot_records.extend([{"Platform": platform, "SampleSize": size}] * count)

        if plot_records:
            df_hist = pd.DataFrame(plot_records)

            plt.figure(figsize=(12, 7))

            # 3. Plot Stacked Histogram
            sns.histplot(
                data=df_hist,
                x="SampleSize",
                hue="Platform",
                multiple="stack",  # STACKS the bars instead of overlapping them
                binwidth=1,  # Forces 1 bar per integer (1, 2, 3 samples...)
                palette="viridis",  # High contrast palette
                edgecolor="white",  # Adds a thin white line between stack segments
                linewidth=0.5,
            )

            plt.title("Distribution of Samples per Study (Top 6 Platforms)")
            plt.xlabel("Number of Samples in Study")
            plt.ylabel("Count of Studies")

            # Limit X to a reasonable range (e.g. 0-60) to see the "1, 2, 3" bars clearly
            # Most RNA-seq studies are small; long tails will be cut off but the detail is preserved
            plt.xlim(0, 60)

            # Optional: Force integer ticks on X axis for clarity
            plt.xticks(range(0, 61, 5))

            plt.tight_layout()
            plt.savefig(f"{output_dir}/tracker_histogram_top6_stacked.svg", format="svg")
            plt.close()
            print(f"Saved {output_dir}/tracker_histogram_top6_stacked.svg")
        else:
            print("Top platforms found, but no sample distribution data available.")
    # --- PLOT 4: HISTOGRAM (SAMPLE SIZES) ---
    # Reconstruct the raw list of sample sizes from the frequency dictionary
    # dict: {"6": 10, "12": 5} -> list: [6, 6, ..., 12, 12, ...]
    # sample_sizes = []
    # if sample_dist:
    #     for size_str, count in sample_dist.items():
    #         size = int(size_str)
    #         sample_sizes.extend([size] * count)
    #     plt.figure(figsize=(10, 6))
    #     sns.histplot(sample_sizes, kde=False, bins=30, color='skyblue', log_scale=(False, False)) # Can toggle log_scale=(True, False) for X axis

    #     plt.title("Distribution of Samples per Study")
    #     plt.xlabel("Number of Samples in Study")
    #     plt.ylabel("Count of Studies")

    #     # Add a vertical line for the mean/median
    #     median_val = float(np.median(sample_sizes))
    #     plt.axvline(median_val, color='r', linestyle='--', label=f'Median: {int(median_val)}')
    #     plt.legend()
    #     plt.savefig(f"{output_dir}/tracker_histogram.svg", format='svg')
    #     plt.close()
    #     print(f"Saved {output_dir}/tracker_histogram.svg")
    # else:
    #     print("No sample distribution data available to plot.")

    if dist_per_platform:
        # 1. Identify Top 6 Platforms by Total Studies Seen
        # Sorting format: (PlatformName, {stats...})
        sorted_platforms = sorted(platform_data.items(), key=lambda item: item[1]["studies_seen"], reverse=True)
        top_6_platforms = [p[0] for p in sorted_platforms[:6]]

        # 2. Reconstruct Data for DataFrame
        plot_records = []
        for platform in top_6_platforms:
            # Get the frequency dict for this platform: {"12": 5, "24": 1}
            freq_dict = dist_per_platform.get(platform, {})
            for size_str, count in freq_dict.items():
                size = int(size_str)
                # Add 'count' rows for this size
                plot_records.extend([{"Platform": platform, "SampleSize": size}] * count)

        if plot_records:
            df_hist = pd.DataFrame(plot_records)

            plt.figure(figsize=(12, 7))

            # 3. Plot Overlapping Histogram
            # element="step" creates the outline look which handles overlap better than bars
            sns.histplot(
                data=df_hist,
                x="SampleSize",
                hue="Platform",
                element="step",
                bins=30,
                common_norm=False,  # Normalize each platform independently? False = raw counts
                log_scale=(False, False),  # Set (True, False) for log X axis if needed
                palette="bright",
                alpha=0.3,
            )

            plt.title("Distribution of Samples per Study (Top 6 Platforms)")
            plt.xlabel("Number of Samples in Study")
            plt.ylabel("Count of Studies")
            plt.xlim(0, 100)  # Optional: limit X axis if there are massive outliers

            plt.savefig(f"{output_dir}/tracker_histogram_top6.svg", format="svg")
            plt.close()
            print(f"Saved {output_dir}/tracker_histogram_top6.svg")
        else:
            print("Top platforms found, but no sample distribution data available.")
    else:
        print("No per-platform sample distribution data found (Running with old Tracker?).")


def combine_files_microarray(folder: str, new_file_name: str, new_file_location: str, combine_genes: bool = False, combination_method: str | Callable = "mean") -> tuple[pd.DataFrame, dict[str, str]]:
    """
    Combines individual GSE RMA CSV files into a single large expression matrix
    and generates a mapping of Samples to Study IDs.

    Args:
        folder (str): Root directory containing GSE subfolders.
        new_file_name (str): Output filename (e.g. 'combined_data.csv').
        new_file_location (str): Output directory.
        combine_genes (bool): If True, merges 'Gene.1', 'Gene.2' into 'Gene'.
        combination_method (str or callable): Method to merge overlapping values.

    Returns:
        Tuple[pd.DataFrame, Dict]: (The combined expression dataframe, Dictionary {SampleID: StudyID})
    """

    # 1. Setup Output Directory
    if not os.path.exists(new_file_location):
        try:
            os.makedirs(new_file_location)
        except OSError as e:
            print(f"Error creating directory {new_file_location}: {e}")
            raise e

    dataframes = []
    sample_to_study_map = {}  # Dictionary to store Sample -> Study relationship

    print(f"Scanning '{folder}' for processed files...")

    # 2. Iterate over the folder structure
    for root, dirs, _files in os.walk(folder):
        for d in dirs:
            # The folder name 'd' is assumed to be the Study ID (e.g., GSE12345)
            study_id = d

            expected_csv_name = f"{d}_RMA_LocusID.csv"
            file_path = os.path.join(root, d, expected_csv_name)

            if os.path.exists(file_path):
                print(f"  - Found: {expected_csv_name} (Study: {study_id})")
                try:
                    df = pd.read_csv(file_path, index_col=0)

                    # Requirement: Set all genes to full capitalization
                    df.index = df.index.str.upper()

                    # Local deduplication
                    if df.index.duplicated().any():
                        df = df.groupby(df.index).mean()

                    # --- NEW LOGIC: Map Samples to Study ---
                    # df.columns contains the GSM (Sample) IDs
                    for sample_id in df.columns:
                        sample_to_study_map[sample_id] = study_id
                    # ---------------------------------------

                    dataframes.append(df)

                except Exception as e:
                    print(f"    ! Error reading {expected_csv_name}: {e}")

    if not dataframes:
        msg = "No valid CSV files found to combine."
        raise ValueError(msg)

    # 3. Combine Dataframes
    print(f"\nMerging {len(dataframes)} datasets... (This may take a moment)")

    combined_df = pd.concat(dataframes, axis=1, join="outer", sort=True)

    # 4. Optional: Combine Gene Variants (Gene.1 -> Gene)
    if combine_genes:
        print(f"  - Consolidating gene variants (e.g., 'Gene.1' -> 'Gene') using method: {combination_method}...")
        cleaned_index = combined_df.index.str.split(".").str[0]
        combined_df.index = cleaned_index
        combined_df = combined_df.groupby(combined_df.index).agg(combination_method)

    elif combined_df.index.duplicated().any():
        print("  - resolving duplicate indices in final merge (using mean)...")
        combined_df = combined_df.groupby(combined_df.index).mean()

    # 5. Save Data Matrix
    output_path = os.path.join(new_file_location, new_file_name)
    print(f"Saving combined matrix to: {output_path}")

    try:
        combined_df.to_csv(output_path)
        print("SUCCESS: Data matrix saved.")
        print(f"Final Dimensions: {combined_df.shape[0]} Genes x {combined_df.shape[1]} Samples")
    except Exception as e:
        print(f"Error saving data file: {e}")

    # 6. Save Sample Map
    # We generate a filename based on the input name: "data.csv" -> "data_sample_map.csv"
    base_name = os.path.splitext(new_file_name)[0]
    map_filename = f"{base_name}_sample_map.csv"
    map_output_path = os.path.join(new_file_location, map_filename)

    print(f"Saving sample map to: {map_output_path}")
    try:
        # Convert dict to DataFrame for easy saving
        map_df = pd.DataFrame.from_dict(sample_to_study_map, orient="index", columns=["StudyID"])
        map_df.index.name = "SampleID"
        map_df.to_csv(map_output_path)
        print("SUCCESS: Sample map saved.")
    except Exception as e:
        print(f"Error saving map file: {e}")

    return combined_df, sample_to_study_map


# Force matplotlib to not use any Xwindow/GUI backend.
# This must be called BEFORE importing pyplot.
