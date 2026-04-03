#!/usr/bin/env python3
"""
monitor_jobs.py — RNA-seq pipeline SLURM job monitor
=====================================================
Watches one or more SLURM array job log files, parses study/batch state in
real-time, cross-references against the FileTracker on disk, and writes a
self-refreshing HTML dashboard.

Usage:
    python3 monitor_jobs.py 12348004 12351666
    python3 monitor_jobs.py 12348004 --log-dir ~/Dataset_fusion_Microarray/logs_slurm
    python3 monitor_jobs.py 12348004 --tracker-dir /tudelft.net/.../rnaseq_data/file_tracker
    python3 monitor_jobs.py 12348004 --once          # single snapshot, no loop
    python3 monitor_jobs.py 12348004 --interval 30   # refresh every 30s (default 20)
    python3 monitor_jobs.py 12348004 --out dashboard.html

Tracker vs log comparison
--------------------------
The FileTracker (.txt files, one per GSE) is the ground truth written by the
pipeline. The log is the in-progress narrative which can be stale or truncated.
Discrepancies are highlighted as warnings in the dashboard:

  LOG says processed  + TRACKER says downloaded → tracker write may have failed
  LOG says error      + TRACKER says processed  → stale log from a previous run
  LOG says ignore     + TRACKER says not_tried  → mark_ignore may have silently failed
  TRACKER has studies + no log entry at all     → studies processed by a prior job
"""

import argparse
import glob
import os
import re
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone, timedelta
from pathlib import Path

# ──────────────────────────────────────────────────────────────────────────────
# TIME UTILITIES
# ──────────────────────────────────────────────────────────────────────────────

# The log timestamps come from bash `date`, e.g. "Wed Apr  1 14:07:51 CEST 2026"
# strptime doesn't handle timezone abbreviations reliably, so we strip the
# timezone name and parse the rest as local time. Since both start and end come
# from the same node we only need the difference, so the absolute TZ doesn't matter.
_DATE_FMT = "%a %b %d %H:%M:%S %Y"
_RE_TZ_NAME = re.compile(r"\s+[A-Z]{2,5}\s+")  # e.g. " CEST " or " UTC "

def _parse_log_timestamp(ts: str):
    """Parse a bash `date` string into a datetime, or return None."""
    if not ts:
        return None
    try:
        clean = _RE_TZ_NAME.sub(" ", ts).strip()
        return datetime.strptime(clean, _DATE_FMT)
    except ValueError:
        return None

def _format_duration(seconds: float) -> str:
    """Format elapsed seconds as a human-readable string, e.g. '2h 14m' or '37m 05s'."""
    seconds = int(seconds)
    h, rem = divmod(seconds, 3600)
    m, s   = divmod(rem, 60)
    if h > 0:
        return f"{h}h {m:02d}m"
    if m > 0:
        return f"{m}m {s:02d}s"
    return f"{s}s"

# ──────────────────────────────────────────────────────────────────────────────
# TRACKER READING
# Read FileTracker state directly from .txt files — no project imports needed.
# One file per GSE: contains a single integer matching STATUS_* constants.
# ──────────────────────────────────────────────────────────────────────────────

# Mirrors src/constants.py — do not import from the project to keep this script
# fully standalone.
_TRACKER_INT_TO_NAME = {
    0: "not_tried",
    1: "downloaded",
    2: "processed",
    3: "ignore",
    4: "error",
}

# Maps tracker state name → (log state name it aligns with, display label, color)
TRACKER_STATES = {
    "not_tried":  ("#475569", "Not tried"),
    "downloaded": ("#a78bfa", "Downloaded"),
    "processed":  ("#4ade80", "✓ Processed"),
    "ignore":     ("#94a3b8", "Ignored"),
    "error":      ("#f87171", "Error"),
    "unknown":    ("#334155", "Unknown"),
}

# Log state → closest equivalent tracker state (for comparison)
_LOG_TO_TRACKER_EQUIV = {
    "processed":   "processed",
    "samplesheet": "downloaded",   # samplesheet done = still downloaded in tracker
    "downloaded":  "downloaded",
    "processing":  "downloaded",   # mid-flight; tracker might say downloaded already
    "no_fastq":    "ignore",
    "ignore":      "ignore",
    "skip_compat": "ignore",
    "error":       "error",
    "unknown":     "not_tried",
}

# Discrepancy rules: (log_equiv, tracker_state) → warning message
# Only pairs that actually indicate a real problem are listed.
_DISCREPANCY_RULES = {
    ("processed",  "downloaded"): "Log says processed but tracker still shows downloaded — tracker write may have failed",
    ("processed",  "not_tried"):  "Log says processed but tracker shows not_tried — very unexpected",
    ("processed",  "error"):      "Log says processed but tracker shows error — stale log or race condition",
    ("error",      "processed"):  "Log says error but tracker shows processed — log is stale from a prior run",
    ("ignore",     "not_tried"):  "Log says ignored but tracker not updated — mark_ignore may have failed",
    ("ignore",     "error"):      "Log says ignored but tracker shows error — minor inconsistency",
    ("downloaded", "processed"):  "Log says still downloading but tracker already shows processed — log may be from older job",
}


def read_tracker_dir(tracker_dir: str) -> dict:
    """
    Read all GSExxx.txt files from the tracker directory.
    Returns {gse_id: tracker_state_name}.
    Does not raise — returns empty dict if directory doesn't exist.
    """
    states = {}
    if not tracker_dir or not os.path.isdir(tracker_dir):
        return states
    for fname in os.listdir(tracker_dir):
        if not fname.startswith("GSE") or not fname.endswith(".txt"):
            continue
        gse_id = fname[:-4]  # strip .txt
        path = os.path.join(tracker_dir, fname)
        try:
            with open(path, "r") as fh:
                raw = fh.read().strip()
            try:
                code = int(raw)
                states[gse_id] = _TRACKER_INT_TO_NAME.get(code, "unknown")
            except ValueError:
                # Some entries are written as "ignore" or "error" strings
                states[gse_id] = raw if raw in _TRACKER_INT_TO_NAME.values() else "unknown"
        except OSError:
            states[gse_id] = "unknown"
    return states


def compare_log_vs_tracker(log_studies: dict, tracker_states: dict) -> list:
    """
    Cross-reference per-study log states against tracker states.
    Returns a list of discrepancy dicts:
        {gse_id, log_state, tracker_state, message}
    """
    discrepancies = []
    for gse_id, log_state in log_studies.items():
        tracker_state = tracker_states.get(gse_id)
        if tracker_state is None:
            # Study appears in log but has no tracker file yet — not a problem
            # if the job is still running; only flag if job is finished
            continue
        log_equiv = _LOG_TO_TRACKER_EQUIV.get(log_state, "unknown")
        key = (log_equiv, tracker_state)
        if key in _DISCREPANCY_RULES:
            discrepancies.append({
                "gse_id":        gse_id,
                "log_state":     log_state,
                "tracker_state": tracker_state,
                "message":       _DISCREPANCY_RULES[key],
            })
    return discrepancies

# ──────────────────────────────────────────────────────────────────────────────
# LOG PARSING
# ──────────────────────────────────────────────────────────────────────────────

# Patterns matched against each log line
_RE_STUDY       = re.compile(r"=== Processing study: (GSE\d+) ===")
_RE_DOWNLOADED  = re.compile(r"Download completed for (GSE\d+)")
_RE_SAMPLESHEET = re.compile(r"DONE generating sample sheet for: (GSE\d+)")
_RE_NO_FASTQ    = re.compile(r"No valid FASTQ pairs found for (GSE\d+)")
_RE_IGNORE      = re.compile(r"Marking (GSE\d+) as ignore")
_RE_ERROR       = re.compile(r"Marking (GSE\d+) as error")
# _RE_PROCESSED   = re.compile(r"Saved (GSE\d+) counts to ")
_RE_PROCESSED   = re.compile(r"Marking (GSE\d+) as processed")
_RE_SKIP_COMPAT = re.compile(r"\[!\] Skipping (GSE\d+):")
_RE_SUPERSERIES = re.compile(r"SuperSeries detected")
_RE_BATCH_RUN   = re.compile(r"Running nf-core/rnaseq \(Batch Mode\) in .*(batch_\S+)\.\.\.")
_RE_BATCH_OK    = re.compile(r"Demultiplexing batch results")
_RE_BATCH_ERR   = re.compile(r"Nextflow Batch Error")
_RE_BATCH_FATAL = re.compile(r"\[!\] No bad samples identified")
_RE_NF_START    = re.compile(r"N E X T F L O W  ~  version")
_RE_NF_DONE     = re.compile(r"Pipeline completed")
_RE_NF_ERROR    = re.compile(r"Pipeline completed with errors")
_RE_RETRY       = re.compile(r"\[Retry (\d+)/(\d+)\] Re-downloading (\d+) failed SRRs")
_RE_ARRAY_JOB   = re.compile(r"--- ARRAY JOB #(\d+) ---")
_RE_BATCH_IDS   = re.compile(r"IDs: (\[.*?\])")
_RE_JOB_ID      = re.compile(r"Job ID: (\d+)")
_RE_NODE        = re.compile(r"Node: (\S+)")
_RE_DATE        = re.compile(r"^Date: (.+)$")
_RE_JOB_DONE    = re.compile(r"^Job finished")
_RE_TMP_NODEV   = re.compile(r"nodev.*mount option|chmod.*operation not permitted")

STUDY_STATES = {
    "processing":   ("#60a5fa", "Processing"),
    "downloaded":   ("#a78bfa", "Downloaded"),
    "samplesheet":  ("#34d399", "Samplesheet ready"),
    "no_fastq":     ("#f59e0b", "No FASTQs"),
    "ignore":       ("#94a3b8", "Ignored"),
    "error":        ("#f87171", "Error"),
    "processed":    ("#4ade80", "✓ Processed"),
    "skip_compat":  ("#fb923c", "Incompatible"),
    "unknown":      ("#475569", "Unknown"),
}

BATCH_STATES = {
    "running":   ("#60a5fa", "Running"),
    "success":   ("#4ade80", "Success"),
    "error":     ("#f87171", "Error"),
    "unknown":   ("#94a3b8", "Pending"),
}


def parse_log_file(path: str) -> dict:
    """Parse a single log file and return a structured result dict."""
    result = {
        "path":         path,
        "job_id":       None,
        "array_index":  None,
        "node":         None,
        "start_time":   None,
        "end_time":     None,
        "duration_secs": None,   # computed after parsing; None if job still running
        "finished":     False,
        "studies":      {},   # gse_id -> state string
        "study_order":  [],   # ordered list of gse_ids seen
        "batches":      {},   # batch_id -> state
        "retries":      defaultdict(list),  # gse_id -> list of retry counts
        "errors":       [],   # free-text error lines
        "warnings":     [],
        "nf_runs":      0,
        "nodev_error":  False,
        "raw_lines":    0,
    }

    if not os.path.exists(path):
        result["errors"].append(f"Log file not found: {path}")
        return result

    current_study = None
    current_batch = None

    try:
        with open(path, "r", errors="replace") as fh:
            for line in fh:
                result["raw_lines"] += 1
                line = line.rstrip()

                m = _RE_JOB_ID.search(line)
                if m: result["job_id"] = m.group(1)

                m = _RE_ARRAY_JOB.search(line)
                if m: result["array_index"] = m.group(1)

                m = _RE_NODE.search(line)
                if m: result["node"] = m.group(1)

                m = _RE_DATE.match(line)
                if m:
                    ts = m.group(1).strip()
                    if result["start_time"] is None:
                        result["start_time"] = ts
                    else:
                        result["end_time"] = ts

                if _RE_JOB_DONE.match(line):
                    result["finished"] = True

                m = _RE_BATCH_IDS.search(line)
                if m:
                    try:
                        result["studies_in_batch"] = eval(m.group(1))
                    except Exception:
                        pass

                # Study state transitions
                m = _RE_STUDY.search(line)
                if m:
                    current_study = m.group(1)
                    if current_study not in result["studies"]:
                        result["study_order"].append(current_study)
                    result["studies"][current_study] = "processing"

                m = _RE_DOWNLOADED.search(line)
                if m:
                    gse = m.group(1)
                    if result["studies"].get(gse) not in ("samplesheet", "processed"):
                        result["studies"][gse] = "downloaded"

                m = _RE_SAMPLESHEET.search(line)
                if m:
                    result["studies"][m.group(1)] = "samplesheet"

                m = _RE_NO_FASTQ.search(line)
                if m:
                    gse = m.group(1)
                    if result["studies"].get(gse) not in ("processed",):
                        result["studies"][gse] = "no_fastq"

                m = _RE_SKIP_COMPAT.search(line)
                if m:
                    result["studies"][m.group(1)] = "skip_compat"

                m = _RE_IGNORE.search(line)
                if m:
                    gse = m.group(1)
                    if result["studies"].get(gse) not in ("processed",):
                        result["studies"][gse] = "ignore"

                m = _RE_ERROR.search(line)
                if m:
                    result["studies"][m.group(1)] = "error"

                m = _RE_PROCESSED.search(line)
                if m:
                    result["studies"][m.group(1)] = "processed"

                # Batch state
                m = _RE_BATCH_RUN.search(line)
                if m:
                    current_batch = m.group(1)
                    result["batches"][current_batch] = "running"
                    if _RE_NF_START.search(line):
                        result["nf_runs"] += 1

                if _RE_NF_START.search(line):
                    result["nf_runs"] += 1

                if _RE_BATCH_OK.search(line) and current_batch:
                    result["batches"][current_batch] = "success"

                if _RE_BATCH_ERR.search(line) and current_batch:
                    result["batches"][current_batch] = "error"
                    result["errors"].append(f"Batch {current_batch}: Nextflow error")

                if _RE_TMP_NODEV.search(line):
                    result["nodev_error"] = True
                    result["warnings"].append("⚠ /tmp nodev error detected — APPTAINER_TMPDIR not set correctly")

                m = _RE_RETRY.search(line)
                if m and current_study:
                    result["retries"][current_study].append(
                        f"Retry {m.group(1)}/{m.group(2)}: {m.group(3)} SRRs"
                    )

    except Exception as e:
        result["errors"].append(f"Failed to read {path}: {e}")

    # Compute duration from timestamps
    t_start = _parse_log_timestamp(result["start_time"])
    t_end   = _parse_log_timestamp(result["end_time"])
    if t_start and t_end:
        result["duration_secs"] = (t_end - t_start).total_seconds()
    elif t_start and not result["finished"]:
        # Job still running — compute elapsed time against now
        result["duration_secs"] = (datetime.now() - t_start).total_seconds()

    return result


def find_log_files(job_ids: list, log_dir: str) -> dict:
    """Find all log files for the given job IDs. Returns {task_id: path}."""
    found = {}
    for jid in job_ids:
        pattern = os.path.join(log_dir, f"stdout-*-{jid}_*.txt")
        matches = glob.glob(pattern)
        # Also try pattern where job_id is the array parent
        if not matches:
            pattern2 = os.path.join(log_dir, f"stdout-*{jid}*.txt")
            matches = glob.glob(pattern2)
        for path in sorted(matches):
            # Extract task ID from filename
            m = re.search(r"_(\d+)\.txt$", path)
            task_id = m.group(1) if m else os.path.basename(path)
            found[f"{jid}_{task_id}"] = path
    return found


# ──────────────────────────────────────────────────────────────────────────────
# HTML GENERATION
# ──────────────────────────────────────────────────────────────────────────────

def _bar(value, total, color):
    pct = (value / total * 100) if total > 0 else 0
    return f'<div class="bar-fill" style="width:{pct:.1f}%;background:{color}"></div>'


def _state_counts(studies_dict):
    counts = defaultdict(int)
    for state in studies_dict.values():
        counts[state] += 1
    return counts


def generate_html(all_results: list, refresh_s: int, generated_at: str,
                  tracker_states = None) -> str:
    total_studies_all = sum(len(r["studies"]) for r in all_results)
    total_processed   = sum(1 for r in all_results for s in r["studies"].values() if s == "processed")
    total_errors      = sum(1 for r in all_results for s in r["studies"].values() if s == "error")
    total_ignored     = sum(1 for r in all_results for s in r["studies"].values() if s in ("ignore", "skip_compat", "no_fastq"))
    active_jobs       = sum(1 for r in all_results if not r["finished"])

    # ── TRACKER SIDEBAR ──────────────────────────────────────────────────────
    # Aggregate tracker states globally — this is the ground truth view
    tracker_states = tracker_states or {}
    tracker_counts = defaultdict(int)
    for state in tracker_states.values():
        tracker_counts[state] += 1
    tracker_total = len(tracker_states)

    tracker_sidebar_html = ""
    if tracker_states:
        tk_legend = ""
        for state, (color, label) in TRACKER_STATES.items():
            n = tracker_counts.get(state, 0)
            if n == 0:
                continue
            pct = n / tracker_total * 100 if tracker_total else 0
            tk_legend += f'''
            <div class="legend-item">
              <span class="legend-dot" style="background:{color}"></span>
              <span class="legend-label">{label}</span>
              <span class="legend-count">{n}</span>
            </div>
            <div class="tk-bar-row">
              <div class="tk-bar"><div class="tk-bar-fill" style="width:{pct:.1f}%;background:{color}"></div></div>
              <span class="tk-pct">{pct:.0f}%</span>
            </div>'''

        # Count discrepancies across all jobs
        total_discrepancies = sum(
            len(compare_log_vs_tracker(r["studies"], tracker_states))
            for r in all_results
        )
        disc_badge = ""
        if total_discrepancies:
            disc_badge = f'<div class="disc-total">⚠ {total_discrepancies} discrepanc{"y" if total_discrepancies==1 else "ies"} found</div>'

        tracker_sidebar_html = f'''
        <div class="sidebar-card">
          <div class="sidebar-title">Tracker on disk ({tracker_total} GSEs)</div>
          {disc_badge}
          {tk_legend}
        </div>'''

    # Aggregate state counts across all results
    agg_states = defaultdict(int)
    for r in all_results:
        for state in r["studies"].values():
            agg_states[state] += 1

    def donut_segments(counts, total):
        """Return SVG path segments for a donut chart."""
        if total == 0:
            return '<circle cx="60" cy="60" r="45" fill="none" stroke="#1e2433" stroke-width="18"/>'
        order = ["processed", "samplesheet", "downloaded", "processing", "no_fastq", "skip_compat", "ignore", "error", "unknown"]
        segs = []
        angle = -90  # start at top
        r, cx, cy, sw = 45, 60, 60, 18
        circumference = 2 * 3.14159 * r
        for state in order:
            n = counts.get(state, 0)
            if n == 0:
                continue
            color = STUDY_STATES[state][0]
            sweep = (n / total) * 360
            # SVG arc
            start_rad = angle * 3.14159 / 180
            end_rad   = (angle + sweep) * 3.14159 / 180
            x1 = cx + r * (end_rad - start_rad if False else 1) * 0  # simplified
            # Use stroke-dasharray technique on a circle instead
            pct = n / total
            segs.append((pct, color, state, n))
            angle += sweep

        # Build using stroke-dashoffset technique
        result = []
        offset = 0
        for pct, color, state, n in segs:
            dash = pct * circumference
            result.append(
                f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="none" stroke="{color}" '
                f'stroke-width="{sw}" stroke-dasharray="{dash:.2f} {circumference:.2f}" '
                f'stroke-dashoffset="-{offset:.2f}" style="transform:rotate(-90deg);transform-origin:center"/>'
            )
            offset += dash
        return "\n".join(result)

    donut_svg = donut_segments(agg_states, total_studies_all)

    # Legend items
    legend_html = ""
    for state, (color, label) in STUDY_STATES.items():
        n = agg_states.get(state, 0)
        if n == 0:
            continue
        legend_html += f'''
        <div class="legend-item">
          <span class="legend-dot" style="background:{color}"></span>
          <span class="legend-label">{label}</span>
          <span class="legend-count">{n}</span>
        </div>'''

    # Per-job cards
    job_cards_html = ""
    for r in all_results:
        counts   = _state_counts(r["studies"])
        n_total  = len(r["studies"])
        finished = r["finished"]
        status_badge = ('<span class="badge badge-green">Finished</span>' if finished
                        else '<span class="badge badge-blue anim-pulse">Running</span>')

        # Duration display
        dur_secs = r.get("duration_secs")
        if dur_secs is not None:
            dur_str = _format_duration(dur_secs)
            if finished:
                duration_html = f'<span class="duration-chip dur-done">⏱ {dur_str}</span>'
            else:
                duration_html = f'<span class="duration-chip dur-live">⟳ {dur_str} elapsed</span>'
        else:
            duration_html = ""

        nodev_warn = ""
        if r.get("nodev_error"):
            nodev_warn = '<div class="inline-warn">⚠ /tmp nodev error — set APPTAINER_TMPDIR</div>'

        # Study rows — with tracker comparison column
        study_rows = ""
        discrepancies = compare_log_vs_tracker(r["studies"], tracker_states)
        disc_by_gse = {d["gse_id"]: d for d in discrepancies}

        for gse_id in r.get("study_order", r["studies"].keys()):
            log_state = r["studies"].get(gse_id, "unknown")
            log_color, log_label = STUDY_STATES.get(log_state, STUDY_STATES["unknown"])
            retries = r["retries"].get(gse_id, [])
            retry_html = ""
            if retries:
                retry_html = f'<span class="retry-tag">{retries[-1]}</span>'

            # Tracker column
            tk_state = tracker_states.get(gse_id)
            if tk_state:
                tk_color, tk_label = TRACKER_STATES.get(tk_state, TRACKER_STATES["unknown"])
                tracker_cell = f'<span class="state-pill" style="--c:{tk_color}">{tk_label}</span>'
            elif tracker_states:
                # Tracker dir was given but this GSE has no file yet
                tracker_cell = '<span class="muted small">no file</span>'
            else:
                tracker_cell = ""  # no tracker dir given

            # Discrepancy warning
            disc = disc_by_gse.get(gse_id)
            disc_cell = ""
            if disc:
                disc_cell = f'<td class="disc-cell" title="{disc["message"]}">⚠</td>'
            elif tracker_states:
                disc_cell = "<td></td>"

            tracker_col = f"<td>{tracker_cell}</td>{disc_cell}" if tracker_states else ""

            study_rows += f'''
            <tr{"" if not disc else ' class="row-warn"'}>
              <td class="mono">{gse_id}</td>
              <td><span class="state-pill" style="--c:{log_color}">{log_label}</span>{retry_html}</td>
              {tracker_col}
            </tr>'''

        # Studies table header
        extra_headers = "<th>Tracker</th><th></th>" if tracker_states else ""

        # Batch rows
        batch_rows = ""
        for batch_id, state in r.get("batches", {}).items():
            color, label = BATCH_STATES.get(state, BATCH_STATES["unknown"])
            batch_rows += f'''
            <tr>
              <td class="mono small">{batch_id}</td>
              <td><span class="state-pill" style="--c:{color}">{label}</span></td>
            </tr>'''
        if not batch_rows:
            batch_rows = '<tr><td colspan="2" class="muted small">No batches started yet</td></tr>'

        # Mini stacked progress bar
        bar_html = '<div class="mini-bar">'
        bar_order = ["processed", "samplesheet", "downloaded", "processing", "no_fastq", "skip_compat", "ignore", "error"]
        for s in bar_order:
            n = counts.get(s, 0)
            if n and n_total:
                pct = n / n_total * 100
                color = STUDY_STATES[s][0]
                bar_html += f'<div style="width:{pct:.1f}%;background:{color}" title="{STUDY_STATES[s][1]}: {n}"></div>'
        bar_html += "</div>"

        errs_html = ""
        for e in r.get("errors", []):
            errs_html += f'<div class="error-line">{e}</div>'

        path_name = os.path.basename(r["path"])
        job_cards_html += f'''
        <div class="job-card">
          <div class="job-header">
            <div>
              <span class="job-title">Job {r.get("job_id","?")} · Array task {r.get("array_index","?")}</span>
              <span class="job-node">{r.get("node","")}</span>
            </div>
            <div style="display:flex;align-items:center;gap:8px">{duration_html}{status_badge}</div>
          </div>
          <div class="job-meta mono small">
            {path_name}<br>
            Started: {r.get("start_time","—")} &nbsp;·&nbsp; {"Ended: " + r.get("end_time","—") if finished else "still running"}
          </div>
          {nodev_warn}
          {bar_html}
          <div class="counts-row">
            {"".join(f'<span class="count-chip" style="--c:{STUDY_STATES[s][0]}">{counts.get(s,0)} {STUDY_STATES[s][1]}</span>' for s in bar_order if counts.get(s,0))}
          </div>

          <div class="section-grid">
            <div>
              <div class="sub-title">Studies ({n_total})</div>
              <table class="inner-table">
                <thead><tr><th>GSE</th><th>Log</th>{extra_headers}</tr></thead>
                <tbody>{study_rows}</tbody>
              </table>
              {f'<div class="disc-summary">⚠ {len(discrepancies)} discrepanc{"y" if len(discrepancies)==1 else "ies"} — hover ⚠ for details</div>' if discrepancies else ""}
            </div>
            <div>
              <div class="sub-title">Nextflow Batches ({len(r.get("batches",{}))})</div>
              <table class="inner-table">{batch_rows}</table>
            </div>
          </div>
          {errs_html}
        </div>'''

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta http-equiv="refresh" content="{refresh_s}">
<title>RNA-seq Monitor</title>
<style>
@import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=IBM+Plex+Sans:wght@400;500;600;700&display=swap');
:root {{
  --bg:#0a0d12; --s1:#111520; --s2:#181d2a; --border:#1f2638;
  --text:#dde3f0; --muted:#4a5577; --mono:'IBM Plex Mono',monospace;
  --sans:'IBM Plex Sans',sans-serif;
}}
*{{box-sizing:border-box;margin:0;padding:0}}
body{{background:var(--bg);color:var(--text);font-family:var(--sans);font-size:14px;line-height:1.6}}
a{{color:inherit}}

/* ── HEADER ── */
header{{padding:28px 40px 22px;border-bottom:1px solid var(--border);display:flex;align-items:center;justify-content:space-between;flex-wrap:wrap;gap:12px}}
.header-left h1{{font-size:22px;font-weight:700;letter-spacing:-.02em}}
.header-left p{{color:var(--muted);font-size:12px;font-family:var(--mono);margin-top:4px}}
.header-stats{{display:flex;gap:24px}}
.stat{{text-align:center}}
.stat-n{{font-size:28px;font-weight:700;line-height:1}}
.stat-l{{font-size:11px;color:var(--muted);font-family:var(--mono);text-transform:uppercase;letter-spacing:.1em}}
.stat-n.green{{color:#4ade80}} .stat-n.red{{color:#f87171}} .stat-n.blue{{color:#60a5fa}} .stat-n.amber{{color:#f59e0b}}

/* ── LAYOUT ── */
.container{{max-width:1400px;margin:0 auto;padding:32px 40px;display:grid;grid-template-columns:260px 1fr;gap:28px;align-items:start}}
@media(max-width:900px){{.container{{grid-template-columns:1fr;padding:20px}}}}

/* ── SIDEBAR ── */
.sidebar{{position:sticky;top:20px;display:flex;flex-direction:column;gap:16px}}
.sidebar-card{{background:var(--s1);border:1px solid var(--border);border-radius:10px;padding:20px}}
.sidebar-title{{font-size:11px;font-family:var(--mono);text-transform:uppercase;letter-spacing:.15em;color:var(--muted);margin-bottom:14px}}
.donut-wrap{{display:flex;justify-content:center;margin-bottom:14px}}
.legend-item{{display:flex;align-items:center;gap:8px;padding:3px 0}}
.legend-dot{{width:9px;height:9px;border-radius:50%;flex-shrink:0}}
.legend-label{{flex:1;color:#94a3b8;font-size:13px}}
.legend-count{{font-family:var(--mono);font-size:12px;color:var(--text);font-weight:600}}

/* ── JOB CARDS ── */
.jobs-col{{display:flex;flex-direction:column;gap:20px}}
.job-card{{background:var(--s1);border:1px solid var(--border);border-radius:10px;padding:22px;display:flex;flex-direction:column;gap:14px}}
.job-header{{display:flex;justify-content:space-between;align-items:flex-start;gap:8px}}
.job-title{{font-weight:700;font-size:15px}}
.job-node{{font-family:var(--mono);font-size:11px;color:var(--muted);display:block;margin-top:2px}}
.job-meta{{color:var(--muted);line-height:1.6}}

/* ── BADGES ── */
.badge{{font-family:var(--mono);font-size:10px;text-transform:uppercase;letter-spacing:.1em;padding:3px 9px;border-radius:4px;font-weight:600}}
.badge-green{{background:rgba(74,222,128,.15);color:#4ade80;border:1px solid rgba(74,222,128,.3)}}
.badge-blue{{background:rgba(96,165,250,.15);color:#60a5fa;border:1px solid rgba(96,165,250,.3)}}
@keyframes pulse{{0%,100%{{opacity:1}}50%{{opacity:.5}}}}
.anim-pulse{{animation:pulse 2s infinite}}

/* ── PROGRESS BAR ── */
.mini-bar{{display:flex;height:7px;border-radius:4px;overflow:hidden;background:var(--border)}}
.mini-bar>div{{height:100%;transition:width .4s}}
.counts-row{{display:flex;flex-wrap:wrap;gap:6px}}
.count-chip{{font-family:var(--mono);font-size:10px;padding:2px 8px;border-radius:4px;background:color-mix(in srgb, var(--c) 15%, transparent);color:var(--c);border:1px solid color-mix(in srgb, var(--c) 30%, transparent)}}

/* ── INNER TABLE ── */
.section-grid{{display:grid;grid-template-columns:1fr 1fr;gap:16px}}
@media(max-width:700px){{.section-grid{{grid-template-columns:1fr}}}}
.sub-title{{font-size:11px;font-family:var(--mono);text-transform:uppercase;letter-spacing:.12em;color:var(--muted);margin-bottom:8px}}
.inner-table{{width:100%;border-collapse:collapse;font-size:12px}}
.inner-table td{{padding:4px 6px;border-bottom:1px solid var(--border);vertical-align:middle}}
.inner-table tr:last-child td{{border-bottom:none}}

/* ── STATE PILL ── */
.state-pill{{display:inline-block;font-size:11px;padding:1px 7px;border-radius:3px;background:color-mix(in srgb,var(--c) 15%,transparent);color:var(--c);border:1px solid color-mix(in srgb,var(--c) 25%,transparent);white-space:nowrap}}
.retry-tag{{font-family:var(--mono);font-size:10px;color:#f59e0b;margin-left:5px;opacity:.8}}

.mono{{font-family:var(--mono)}} .small{{font-size:11px}} .muted{{color:var(--muted)}}
.error-line{{background:rgba(248,113,113,.1);border:1px solid rgba(248,113,113,.2);color:#f87171;border-radius:5px;padding:7px 10px;font-size:12px;font-family:var(--mono)}}
.inline-warn{{background:rgba(245,158,11,.1);border:1px solid rgba(245,158,11,.2);color:#f59e0b;border-radius:5px;padding:7px 10px;font-size:12px}}
.refresh-note{{font-size:11px;color:var(--muted);font-family:var(--mono);text-align:center;padding-bottom:20px}}
/* ── TRACKER ── */
.tk-bar-row{{display:flex;align-items:center;gap:6px;margin:-2px 0 6px 17px}}
.tk-bar{{flex:1;height:3px;background:var(--border);border-radius:2px;overflow:hidden}}
.tk-bar-fill{{height:100%;border-radius:2px;transition:width .4s}}
.tk-pct{{font-family:var(--mono);font-size:9px;color:var(--muted);min-width:24px;text-align:right}}
.disc-total{{background:rgba(245,158,11,.12);border:1px solid rgba(245,158,11,.25);color:#f59e0b;border-radius:5px;padding:5px 9px;font-size:11px;font-family:var(--mono);margin-bottom:10px}}
.disc-summary{{font-size:11px;color:#f59e0b;font-family:var(--mono);margin-top:6px;padding:4px 6px;background:rgba(245,158,11,.08);border-radius:4px}}
.disc-cell{{color:#f59e0b;font-size:13px;cursor:help;text-align:center;width:20px}}
.row-warn td{{background:rgba(245,158,11,.04)}}
.inner-table thead th{{font-family:var(--mono);font-size:9px;text-transform:uppercase;letter-spacing:.12em;color:var(--muted);padding:4px 6px;border-bottom:1px solid var(--border)}}
/* ── DURATION ── */
.duration-chip{{font-family:var(--mono);font-size:11px;font-weight:600;padding:3px 10px;border-radius:4px;white-space:nowrap}}
.dur-done{{background:rgba(74,222,128,.12);color:#4ade80;border:1px solid rgba(74,222,128,.25)}}
.dur-live{{background:rgba(96,165,250,.12);color:#60a5fa;border:1px solid rgba(96,165,250,.25);animation:pulse 2s infinite}}
</style>
</head>
<body>

<header>
  <div class="header-left">
    <h1>RNA-seq Pipeline Monitor</h1>
    <p>Generated {generated_at} · auto-refresh every {refresh_s}s</p>
  </div>
  <div class="header-stats">
    <div class="stat"><div class="stat-n blue">{len(all_results)}</div><div class="stat-l">Tasks</div></div>
    <div class="stat"><div class="stat-n">{active_jobs}</div><div class="stat-l">Running</div></div>
    <div class="stat"><div class="stat-n green">{total_processed}</div><div class="stat-l">Processed</div></div>
    <div class="stat"><div class="stat-n amber">{total_ignored}</div><div class="stat-l">Ignored</div></div>
    <div class="stat"><div class="stat-n red">{total_errors}</div><div class="stat-l">Errors</div></div>
  </div>
</header>

<div class="container">
  <aside class="sidebar">
    <div class="sidebar-card">
      <div class="sidebar-title">All Studies ({total_studies_all})</div>
      <div class="donut-wrap">
        <svg viewBox="0 0 120 120" width="120" height="120">
          {donut_svg}
          <text x="60" y="56" text-anchor="middle" font-size="18" font-weight="700" fill="#dde3f0" font-family="IBM Plex Sans">{total_studies_all}</text>
          <text x="60" y="70" text-anchor="middle" font-size="9" fill="#4a5577" font-family="IBM Plex Mono">STUDIES</text>
        </svg>
      </div>
      {legend_html}
    </div>
    {tracker_sidebar_html}
  </aside>

  <main class="jobs-col">
    {job_cards_html}
    <div class="refresh-note">Page auto-refreshes every {refresh_s}s — {generated_at}</div>
  </main>
</div>

</body>
</html>"""


# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Monitor RNA-seq SLURM job logs and generate an HTML dashboard."
    )
    parser.add_argument(
        "job_ids", nargs="+",
        help="SLURM job IDs to monitor (e.g. 12348004 12351666)"
    )
    parser.add_argument(
        "--log-dir", default=None,
        help="Directory containing log files. Defaults to "
             "~/Dataset_fusion_Microarray/logs_slurm"
    )
    parser.add_argument(
        "--tracker-dir", default=None,
        help="Path to the FileTracker directory (contains GSExxx.txt files). "
             "Defaults to /tudelft.net/staff-umbrella/GeneExpressionStorage/rnaseq_data/file_tracker"
    )
    parser.add_argument(
        "--out", default="pipeline_monitor.html",
        help="Output HTML file (default: pipeline_monitor.html)"
    )
    parser.add_argument(
        "--interval", type=int, default=20,
        help="Refresh interval in seconds (default: 20)"
    )
    parser.add_argument(
        "--once", action="store_true",
        help="Run once and exit (no loop)"
    )
    args = parser.parse_args()

    log_dir = args.log_dir or os.path.expanduser(
        "~/Dataset_fusion_Microarray/logs_slurm"
    )
    tracker_dir = args.tracker_dir or (
        "/tudelft.net/staff-umbrella/GeneExpressionStorage/rnaseq_data/file_tracker"
    )

    if not os.path.isdir(log_dir):
        print(f"Error: log directory not found: {log_dir}", file=sys.stderr)
        print("Use --log-dir to specify the correct path.", file=sys.stderr)
        sys.exit(1)

    if not os.path.isdir(tracker_dir):
        print(f"Warning: tracker directory not found: {tracker_dir}", file=sys.stderr)
        print("Continuing without tracker comparison. Use --tracker-dir to specify.", file=sys.stderr)
        tracker_dir = None

    print(f"Monitoring job IDs: {args.job_ids}")
    print(f"Log directory:      {log_dir}")
    print(f"Tracker directory:  {tracker_dir or '(not found — skipping comparison)'}")
    print(f"Output:             {args.out}")
    print(f"Refresh interval:   {args.interval}s")
    print()

    while True:
        log_files = find_log_files(args.job_ids, log_dir)
        tracker_states = read_tracker_dir(tracker_dir) if tracker_dir else {}

        if not log_files:
            print(f"[{datetime.now(timezone.utc).strftime('%H:%M:%S')}] "
                  f"No log files found for job IDs {args.job_ids} in {log_dir}")
        else:
            results = []
            for key, path in sorted(log_files.items()):
                r = parse_log_file(path)
                results.append(r)

            now_str = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
            html = generate_html(results, args.interval, now_str,
                                 tracker_states=tracker_states)

            with open(args.out, "w") as fh:
                fh.write(html)

            # Count discrepancies for terminal summary
            total_disc = sum(
                len(compare_log_vs_tracker(r["studies"], tracker_states))
                for r in results
            )

            total   = sum(len(r["studies"]) for r in results)
            done    = sum(1 for r in results for s in r["studies"].values() if s == "processed")
            errors  = sum(1 for r in results for s in r["studies"].values() if s == "error")
            running = sum(1 for r in results if not r["finished"])
            disc_str = f" | ⚠ {total_disc} discrepancies" if total_disc else ""
            print(
                f"[{datetime.now(timezone.utc).strftime('%H:%M:%S')}] "
                f"{len(results)} tasks | {total} studies (log) | "
                f"{len(tracker_states)} studies (tracker) | "
                f"✓ {done} processed | ✗ {errors} errors | "
                f"⟳ {running} running{disc_str} → {args.out}"
            )

        if args.once:
            break
        time.sleep(args.interval)


if __name__ == "__main__":
    main()