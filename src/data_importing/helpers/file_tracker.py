import json
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

module_dir = "./"
sys.path.append(module_dir)

from src.constants import STATUS_DOWNLOADED, STATUS_ERROR, STATUS_IGNORE, STATUS_LOCKED, STATUS_PROCESSED,STORAGE_DIR  # noqa: E402


# Ensure these match your constants file
class FileTracker:
    def __init__(self, root_dir):
        """
        Initializes the file-based tracker.
        :param root_dir: Path to the 'tracker_status' directory.
        """
        self.tracker_dir = root_dir
        os.makedirs(self.tracker_dir, exist_ok=True)

    def _get_file_path(self, gse_id, ext=".txt"):
        return os.path.join(self.tracker_dir, f"{gse_id}{ext}")

    # --- CORE STATUS METHODS ---
    def get_status(self, gse_id):
        path = self._get_file_path(gse_id)
        if not os.path.exists(path):
            return -1
        try:
            with open(path) as f:
                return int(f.read().strip())
        except Exception as _e:
            return STATUS_ERROR

    def set_status(self, gse_id, status_code):
        path = self._get_file_path(gse_id)
        with open(path, "w") as f:
            f.write(str(status_code))

    # --- MISSING BOOLEAN CHECKS ---
    def is_locked(self, gse_id):
        return self.get_status(gse_id) == STATUS_LOCKED

    def is_processed(self, gse_id):
        return self.get_status(gse_id) == STATUS_PROCESSED

    def is_downloaded(self, gse_id):
        # Checks if status is Downloaded OR Processed (since processed implies downloaded) Not true if run and delete == true
        s = self.get_status(gse_id)
        return s == STATUS_DOWNLOADED  # or s == STATUS_PROCESSED

    def is_ignored(self, gse_id):
        return self.get_status(gse_id) == STATUS_IGNORE

    def is_error(self, gse_id):
        return self.get_status(gse_id) == STATUS_ERROR

    # --- MARKER METHODS ---
    def mark_processed(self, gse_id):
        print(f"Marking {gse_id} as processed")
        self.set_status(gse_id, STATUS_PROCESSED)

    def mark_downloaded(self, gse_id):
        print(f"Marking {gse_id} as downloaded")
        self.set_status(gse_id, STATUS_DOWNLOADED)

    def mark_ignore(self, gse_id):
        print(f"Marking {gse_id} as ignore")
        self.set_status(gse_id, STATUS_IGNORE)

    def mark_error(self, gse_id):
        print(f"Marking {gse_id} as error")
        self.set_status(gse_id, STATUS_ERROR)

    # --- ECOTYPE TRACKING ---
    def mark_ecotype(self, gse_id, ecotype: str):
        """
        Persists the detected ecotype for a study so GEO metadata never needs
        re-fetching on subsequent runs. Writes to: tracker_dir/GSE123_ecotype.txt
        """
        path = self._get_file_path(gse_id, ext="_ecotype.txt")
        with open(path, "w") as f:
            f.write(ecotype.strip())

    def get_ecotype(self, gse_id) -> str | None:
        """
        Returns the previously stored ecotype, or None if not yet detected.
        None means the caller should fetch from GEO and then call mark_ecotype().
        """
        path = self._get_file_path(gse_id, ext="_ecotype.txt")
        if not os.path.exists(path):
            return None
        try:
            with open(path) as f:
                return f.read().strip()
        except Exception:
            return None

    def has_ecotype(self, gse_id) -> bool:
        """Returns True if an ecotype has already been detected and stored."""
        return os.path.exists(self._get_file_path(gse_id, ext="_ecotype.txt"))

    def save_study_metadata(self, gse_id, platform, num_samples, has_raw):
        """
        New method to persist stats per-study.
        Writes to: tracker_dir/GSE123_meta.json
        """
        data = {"gse_id": gse_id, "platform": platform, "num_samples": num_samples, "has_raw": has_raw}
        path = self._get_file_path(gse_id, ext="_meta.json")
        with open(path, "w") as f:
            json.dump(data, f)

    # --- DUMMY METHODS (For Compatibility) ---
    def save_to_json(self, path=None):
        pass  # Auto-save happens instantly on set_status

    @classmethod
    def load_from_json(cls, path):
        # If your code tries to load from a JSON path, we simply
        # redirect it to open the directory instead.
        # We assume 'path' ends in '.../tracker.json', so we strip the filename
        # to get the directory.
        if os.path.isdir(path):
            return cls(path)
        # Fallback: use the parent directory of the json file
        return cls(os.path.dirname(path) + "/tracker_status")

    def generate_detailed_report(self):
        """
        Scans all metadata files and returns a Pandas DataFrame containing
        the status and metadata for every study.
        """
        data = []
        files = os.listdir(self.tracker_dir)

        for filename in files:
            if filename.endswith("_meta.json"):
                try:
                    filepath = os.path.join(self.tracker_dir, filename)
                    with open(filepath) as f:
                        entry = json.load(f)
                        # Ensure keys exist
                        data.append(
                            {"gse_id": entry.get("gse_id", "Unknown"), "platform": entry.get("platform", "Unknown"), "num_samples": entry.get("num_samples", 0), "has_raw": entry.get("has_raw", False)}
                        )
                except Exception:
                    pass  # Skip corrupted files

        return pd.DataFrame(data)

    # ---------------- PLOTTING METHODS ----------------

    def get_pie_charts(self, save_path=f"{STORAGE_DIR}/outputs/scanner_plots/RNA-seq/tracker_pie_charts.svg"):
        """
        Function 1: Produces pie charts showing SRA (Raw Data) availability
        for both Studies and Samples.
        """
        df = self.generate_detailed_report()
        if df.empty:
            return

        _fig, axes = plt.subplots(1, 2, figsize=(12, 6))

        # --- Chart 1: Studies with SRA ---
        study_counts = df["has_raw"].value_counts()
        labels = [f"Raw Data ({study_counts.get(True, 0)})", f"No Raw ({study_counts.get(False, 0)})"]
        axes[0].pie(study_counts, labels=labels, autopct="%1.1f%%", colors=["#66b3ff", "#ff9999"], startangle=90)
        axes[0].set_title("Studies with SRA Available")

        # --- Chart 2: Samples with SRA ---
        # Group by availability and sum the number of samples
        sample_counts = df.groupby("has_raw")["num_samples"].sum()
        labels_samp = [f"Raw Data ({sample_counts.get(True, 0)})", f"No Raw ({sample_counts.get(False, 0)})"]
        axes[1].pie(sample_counts, labels=labels_samp, autopct="%1.1f%%", colors=["#99ff99", "#ffcc99"], startangle=90)
        axes[1].set_title("Samples with SRA Available")

        plt.tight_layout()
        plt.savefig(save_path)
        plt.close()
        print(f"Pie charts saved to {save_path}")

    def produce_study_dis(self, save_path=f"{STORAGE_DIR}/outputs/scanner_plots/RNA-seq/tracker_histogram_top6_stacked.svg"):
        """
        Function 2: Stacked Histogram of Study Sizes (Samples per Study),
        split by the Top 6 Platforms.
        """
        df = self.generate_detailed_report()
        if df.empty:
            return

        # 1. Identify Top 6 Platforms
        top_platforms = df["platform"].value_counts().nlargest(6).index.tolist()

        # 2. Filter Data (Keep only top 6, label others as 'Other' if you wanted,
        # but usually we just show the top 6 for cleanliness)
        df_filtered = df[df["platform"].isin(top_platforms)]

        # 3. Create Plot
        plt.figure(figsize=(12, 6))
        sns.histplot(
            data=df_filtered,
            x="num_samples",
            hue="platform",
            multiple="stack",
            palette="viridis",
            edgecolor=".3",
            linewidth=0.5,
            binwidth=1,  # Adjust based on your data spread (e.g., 1 for small studies, 5 for larger)
            log_scale=(False, False),  # Set (True, False) if x-axis is too wide
        )

        plt.title("Distribution of Samples per Study (Top 6 Platforms)")
        plt.xlabel("Number of Samples in Study")
        plt.ylabel("Number of Studies")
        plt.xlim(0, 50)  # Limit x-axis to focus on common study sizes (optional)

        plt.savefig(save_path)
        plt.close()
        print(f"Study distribution saved to {save_path}")

    def produce_platform_dis(self, save_path=f"{STORAGE_DIR}/outputs/scanner_plots/RNA-seq/tracker_platforms.svg"):
        """
        Function 3: Bar plot showing Total Samples per Platform.
        """
        df = self.generate_detailed_report()
        if df.empty:
            return

        # 1. Aggregate Data
        platform_stats = df.groupby("platform")["num_samples"].sum().reset_index()

        # 2. Sort and take Top 20 for readability
        platform_stats = platform_stats.sort_values("num_samples", ascending=False).head(20)

        # 3. Create Plot
        plt.figure(figsize=(10, 8))
        sns.barplot(data=platform_stats, y="platform", x="num_samples", palette="magma")

        plt.title("Total Samples per Platform (Top 20)")
        plt.xlabel("Total Samples")
        plt.ylabel("Platform (GPL)")
        plt.grid(axis="x", linestyle="--", alpha=0.6)

        plt.tight_layout()
        plt.savefig(save_path)
        plt.close()
        print(f"Platform distribution saved to {save_path}")
