#!/usr/bin/env python3
"""
rnaseq_processing_status_report.py

Builds a current-state table + plots for the RNA-seq processing pipeline,
combining:
  1) The full desired study list (e.g. study_ids/RNA_seq_ids.txt)
  2) The FileTracker status directory (file_tracker/) -- per-study status
     files, _missing_samples.txt files, and _meta.json files
  3) (optional) The sample map CSV used to build the final expression matrix
     (e.g. Salmon_RNAseq_Combined_TPM_sample_map.csv), to cross-check how
     many samples for each study actually made it into the matrix.

It does NOT import your project's src/ code -- it reads the tracker
directory directly using the same file-naming convention as
FileTracker (file_tracker.py), so it is safe to run standalone anywhere
that can see the tracker folder (e.g. on DAIC).

Usage
-----
python rnaseq_processing_status_report.py \
    --study-list /path/to/RNA_seq_ids.txt \
    --tracker-dir /path/to/rnaseq_data/file_tracker/ \
    --sample-map /path/to/Salmon_RNAseq_Combined_TPM_sample_map.csv \
    --output-dir ./rnaseq_status_report/

Only --study-list and --tracker-dir are required. --sample-map is optional;
if omitted, the matrix cross-check columns are skipped.

Outputs (written to --output-dir)
----------------------------------
- status_table.csv            One row per desired study, full state detail
- status_summary.csv          Counts of studies per status category
- missing_samples_detail.csv  One row per (study, missing sample) for
                               semi-complete studies
- status_overview.png         Bar chart: studies per status
- status_by_platform.png      (only if meta.json platform info exists)
                               Stacked bar of status per platform
"""

import argparse
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")  # safe for headless cluster runs
import matplotlib.pyplot as plt
import pandas as pd

# ----------------------------------------------------------------------------
# Status codes -- must match src/constants.py
# ----------------------------------------------------------------------------
STATUS_LOCKED = 0
STATUS_DOWNLOADED = 1
STATUS_PROCESSED = 2
STATUS_IGNORE = 3
STATUS_ERROR = 4
STATUS_SEMI = 5
STATUS_NOT_STARTED = -1  # no status file found at all

STATUS_LABELS = {
    STATUS_NOT_STARTED: "Not yet processed",
    STATUS_LOCKED: "Locked (in progress)",
    STATUS_DOWNLOADED: "Downloaded (not processed)",
    STATUS_PROCESSED: "Processed",
    STATUS_IGNORE: "Skipped (ignored)",
    STATUS_ERROR: "Error",
    STATUS_SEMI: "Semi-complete (missing samples)",
}

# Fixed display order for plots/tables
STATUS_ORDER = [
    STATUS_PROCESSED,
    STATUS_SEMI,
    STATUS_DOWNLOADED,
    STATUS_LOCKED,
    STATUS_ERROR,
    STATUS_IGNORE,
    STATUS_NOT_STARTED,
]

STATUS_COLORS = {
    STATUS_PROCESSED: "#4CAF50",
    STATUS_SEMI: "#FFC107",
    STATUS_DOWNLOADED: "#2196F3",
    STATUS_LOCKED: "#9C27B0",
    STATUS_ERROR: "#F44336",
    STATUS_IGNORE: "#9E9E9E",
    STATUS_NOT_STARTED: "#E0E0E0",
}


# ----------------------------------------------------------------------------
# Loading helpers
# ----------------------------------------------------------------------------
def load_study_list(path):
    """Loads desired GSE IDs from a file that is either comma-separated
    on one line, or one ID per line (or a mix)."""
    with open(path) as f:
        raw = f.read()
    # Split on both commas and newlines, drop empties/whitespace
    ids = [tok.strip() for chunk in raw.splitlines() for tok in chunk.split(",")]
    ids = [i for i in ids if i]
    # De-duplicate while preserving order
    seen = set()
    out = []
    for i in ids:
        if i not in seen:
            seen.add(i)
            out.append(i)
    return out


def get_status(tracker_dir, gse_id):
    path = os.path.join(tracker_dir, f"{gse_id}.txt")
    if not os.path.exists(path):
        return STATUS_NOT_STARTED
    try:
        with open(path) as f:
            return int(f.read().strip())
    except Exception:
        return STATUS_ERROR


def get_missing_samples(tracker_dir, gse_id):
    """Reads the CURRENT _missing_samples.txt (not timestamped archives)."""
    path = os.path.join(tracker_dir, f"{gse_id}_missing_samples.txt")
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]


def get_meta(tracker_dir, gse_id):
    path = os.path.join(tracker_dir, f"{gse_id}_meta.json")
    if not os.path.exists(path):
        return {}
    try:
        with open(path) as f:
            return json.load(f)
    except Exception:
        return {}


def get_ecotype(tracker_dir, gse_id):
    path = os.path.join(tracker_dir, f"{gse_id}_ecotype.txt")
    if not os.path.exists(path):
        return None
    try:
        with open(path) as f:
            return f.read().strip()
    except Exception:
        return None


def guess_study_column(df):
    """Tries to find which column in the sample map identifies the GSE/study."""
    candidates = ["study", "gse", "gse_id", "GSE", "GSE_ID", "series_id", "study_id", "Study"]
    for c in candidates:
        if c in df.columns:
            return c
    return None


# ----------------------------------------------------------------------------
# Main report builder
# ----------------------------------------------------------------------------
def build_status_table(study_ids, tracker_dir, sample_map_path=None, study_col=None):
    rows = []
    for gse in study_ids:
        status = get_status(tracker_dir, gse)
        missing = get_missing_samples(tracker_dir, gse)
        meta = get_meta(tracker_dir, gse)
        ecotype = get_ecotype(tracker_dir, gse)

        rows.append(
            {
                "gse_id": gse,
                "status_code": status,
                "status_label": STATUS_LABELS.get(status, f"Unknown ({status})"),
                "platform": meta.get("platform"),
                "num_samples_meta": meta.get("num_samples"),
                "has_raw": meta.get("has_raw"),
                "num_missing_samples": len(missing),
                "missing_sample_ids": ";".join(missing) if missing else "",
                "ecotype": ecotype,
            }
        )

    df = pd.DataFrame(rows)

    # --- Optional cross-check against the sample map used for the final matrix ---
    if sample_map_path is not None and os.path.exists(sample_map_path):
        smap = pd.read_csv(sample_map_path, index_col=0)
        if study_col is None:
            study_col = guess_study_column(smap)
        if study_col is None:
            print(
                "WARNING: could not auto-detect a study/GSE column in the sample map "
                f"(columns found: {list(smap.columns)}). Skipping matrix cross-check. "
                "Pass --study-col explicitly to enable it."
            )
        else:
            counts = smap[study_col].value_counts()
            df["num_samples_in_matrix"] = df["gse_id"].map(counts).fillna(0).astype(int)
            # Flag studies marked PROCESSED but with zero samples actually in the matrix
            df["in_matrix_but_marked_incomplete"] = (df["num_samples_in_matrix"] > 0) & (
                df["status_code"] != STATUS_PROCESSED
            )
            df["processed_but_missing_from_matrix"] = (df["status_code"] == STATUS_PROCESSED) & (
                df["num_samples_in_matrix"] == 0
            )

    return df


def write_outputs(df, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    # 1. Full per-study table
    table_path = os.path.join(output_dir, "status_table.csv")
    df.sort_values(["status_code", "gse_id"]).to_csv(table_path, index=False)
    print(f"Wrote {table_path}  ({len(df)} studies)")

    # 2. Summary counts
    summary = (
        df["status_code"]
        .value_counts()
        .rename_axis("status_code")
        .reset_index(name="num_studies")
    )
    summary["status_label"] = summary["status_code"].map(STATUS_LABELS)
    summary = summary.sort_values("status_code")
    summary_path = os.path.join(output_dir, "status_summary.csv")
    summary.to_csv(summary_path, index=False)
    print(f"Wrote {summary_path}")
    print("\n=== Current processing state ===")
    print(summary[["status_label", "num_studies"]].to_string(index=False))
    print(f"Total desired studies: {len(df)}")

    # 3. Missing-sample detail (long format: one row per missing sample)
    semi = df[df["status_code"] == STATUS_SEMI]
    detail_rows = []
    for _, r in semi.iterrows():
        ids = r["missing_sample_ids"].split(";") if r["missing_sample_ids"] else []
        for sid in ids:
            detail_rows.append({"gse_id": r["gse_id"], "missing_sample_id": sid})
    detail_path = os.path.join(output_dir, "missing_samples_detail.csv")
    pd.DataFrame(detail_rows).to_csv(detail_path, index=False)
    print(f"Wrote {detail_path}  ({len(detail_rows)} missing samples across {len(semi)} studies)")

    # 4. Cross-check flags, if present
    if "processed_but_missing_from_matrix" in df.columns:
        n_flag = int(df["processed_but_missing_from_matrix"].sum())
        if n_flag:
            print(
                f"\nWARNING: {n_flag} studies are marked PROCESSED in the tracker but have "
                "0 samples in the final matrix's sample map. See "
                "'processed_but_missing_from_matrix' column in status_table.csv."
            )

    # --- Plot 1: overview bar chart ---
    plot_path = os.path.join(output_dir, "status_overview.png")
    ordered = summary.set_index("status_code").reindex(STATUS_ORDER).dropna(how="all")
    ordered["num_studies"] = ordered["num_studies"].fillna(0).astype(int)
    ordered["status_label"] = ordered.index.map(STATUS_LABELS)

    fig, ax = plt.subplots(figsize=(9, 5.5))
    colors = [STATUS_COLORS[c] for c in ordered.index]
    bars = ax.bar(ordered["status_label"], ordered["num_studies"], color=colors)
    for b, val in zip(bars, ordered["num_studies"]):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height(), str(val), ha="center", va="bottom")
    ax.set_ylabel("Number of studies")
    ax.set_title(f"RNA-seq pipeline status ({len(df)} studies total)")
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=150)
    plt.close(fig)
    print(f"Wrote {plot_path}")

    # --- Plot 2: status by platform (only if we have platform metadata) ---
    if df["platform"].notna().any():
        plat_path = os.path.join(output_dir, "status_by_platform.png")
        pivot = (
            df[df["platform"].notna()]
            .groupby(["platform", "status_code"])
            .size()
            .unstack(fill_value=0)
        )
        # keep only columns that actually occur, in fixed order
        cols = [c for c in STATUS_ORDER if c in pivot.columns]
        pivot = pivot[cols]
        pivot.columns = [STATUS_LABELS[c] for c in cols]
        top_platforms = pivot.sum(axis=1).sort_values(ascending=False).head(15).index
        pivot = pivot.loc[top_platforms]

        fig, ax = plt.subplots(figsize=(10, 7))
        pivot.plot(
            kind="barh",
            stacked=True,
            ax=ax,
            color=[STATUS_COLORS[c] for c in cols],
        )
        ax.set_xlabel("Number of studies")
        ax.set_ylabel("Platform")
        ax.set_title("Processing status by platform (top 15 by study count)")
        ax.legend(loc="lower right", fontsize=8)
        plt.tight_layout()
        plt.savefig(plat_path, dpi=150)
        plt.close(fig)
        print(f"Wrote {plat_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study-list", required=True, help="Path to full desired study ID list (comma- or newline-separated GSE IDs)")
    parser.add_argument("--tracker-dir", required=True, help="Path to the FileTracker status directory")
    parser.add_argument("--sample-map", default=None, help="Optional path to the sample map CSV used for the final matrix")
    parser.add_argument("--study-col", default=None, help="Optional column name in the sample map identifying the GSE/study (auto-detected if omitted)")
    parser.add_argument("--output-dir", default="./rnaseq_status_report/", help="Directory to write the table/plots to")
    args = parser.parse_args()

    if not os.path.isdir(args.tracker_dir):
        print(f"ERROR: tracker dir not found: {args.tracker_dir}", file=sys.stderr)
        sys.exit(1)

    study_ids = load_study_list(args.study_list)
    print(f"Loaded {len(study_ids)} desired studies from {args.study_list}")

    df = build_status_table(
        study_ids,
        args.tracker_dir,
        sample_map_path=args.sample_map,
        study_col=args.study_col,
    )
    write_outputs(df, args.output_dir)


if __name__ == "__main__":
    main()