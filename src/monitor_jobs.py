#!/usr/bin/env python3
"""
monitor_jobs.py — RNA-seq pipeline SLURM job monitor
=====================================================
Watches one or more SLURM array job log files, parses study/batch state in
real-time, cross-references against the FileTracker on disk, and writes a
self-refreshing HTML dashboard.

Usage:
    python3 monitor_jobs.py 19290
    python3 monitor_jobs.py 19290 --log-dir ~/Dataset_fusion_Microarray/logs_slurm
    python3 monitor_jobs.py 19290 --tracker-dir /tudelft.net/.../rnaseq_data/file_tracker
    python3 monitor_jobs.py 19290 --once          # single snapshot, no loop
    python3 monitor_jobs.py 19290 --interval 30   # refresh every 30s (default 20)
    python3 monitor_jobs.py 19290 --out dashboard.html

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
    """Read all GSExxx.txt files from the tracker directory."""
    states = {}
    if not tracker_dir or not os.path.isdir(tracker_dir):
        return states
    for fname in os.listdir(tracker_dir):
        if not fname.startswith("GSE") or not fname.endswith(".txt"):
            continue
        gse_id = fname[:-4]
        path = os.path.join(tracker_dir, fname)
        try:
            with open(path, "r") as fh:
                raw = fh.read().strip()
            try:
                code = int(raw)
                states[gse_id] = _TRACKER_INT_TO_NAME.get(code, "unknown")
            except ValueError:
                states[gse_id] = raw if raw in _TRACKER_INT_TO_NAME.values() else "unknown"
        except OSError:
            states[gse_id] = "unknown"
    return states


def compare_log_vs_tracker(log_studies: dict, tracker_states: dict) -> list:
    """Cross-reference per-study log states against tracker states."""
    discrepancies = []
    for gse_id, log_state in log_studies.items():
        tracker_state = tracker_states.get(gse_id)
        if tracker_state is None:
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

# Patterns matched against each log line (Updated for multiple-space formatting)
_RE_STUDY       = re.compile(r"=== Processing study: (GSE\d+) ===")
_RE_DOWNLOADED  = re.compile(r"Download completed for (GSE\d+)")
_RE_SAMPLESHEET = re.compile(r"DONE generating sample sheet for: (GSE\d+)")
_RE_NO_FASTQ    = re.compile(r"No valid FASTQ pairs found for (GSE\d+)")
_RE_IGNORE      = re.compile(r"Marking (GSE\d+) as ignore")
_RE_ERROR       = re.compile(r"Marking (GSE\d+) as error")
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
_RE_NF_RETRY    = re.compile(r"Removing bad sample from samplesheet: (\S+)")
_RE_ARRAY_JOB   = re.compile(r"--- ARRAY JOB #(\d+) ---")
_RE_ARRAY_TASK  = re.compile(r"Array Task ID:\s+(\d+)")
_RE_BATCH_IDS   = re.compile(r"IDs:\s+(\[.*?\])")
_RE_JOB_ID      = re.compile(r"Job ID:\s+(\d+)")
_RE_NODE        = re.compile(r"Node:\s+(\S+)")
_RE_DATE        = re.compile(r"^Date:\s+(.+)$")
_RE_JOB_DONE    = re.compile(r"^(?:Job (finished|Completed)|Trimming batch directory)")
_RE_TMP_NODEV   = re.compile(r"nodev.*mount option|chmod.*operation not permitted")
_RE_SKIP_PROCESSED = re.compile(r"Skipped — already processed \(\d+\):\s+(\[.*?\])")
_RE_SKIP_IGNORED   = re.compile(r"Skipped — ignored \(\d+\):\s+(\[.*?\])")
_RE_SKIP_ERROR     = re.compile(r"Skipped — error \(\d+\):\s+(\[.*?\])")

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
        "duration_secs": None,   
        "finished":     False,
        "studies":      {},   
        "study_order":  [],   
        "batches":      {},   
        "retries":      defaultdict(list),  
        "errors":       [],   
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

                m = _RE_ARRAY_TASK.search(line)
                if m: result["array_index"] = m.group(1)

                m = _RE_ARRAY_JOB.search(line)
                if m and not result["array_index"]: result["array_index"] = m.group(1)

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
                    gse = m.group(1)
                    if result["studies"].get(gse) not in ("processed",):
                        result["studies"][gse] = "error"

                m = _RE_PROCESSED.search(line)
                if m:
                    result["studies"][m.group(1)] = "processed"

                m = _RE_RETRY.search(line)
                if m and current_study:
                    result["retries"][current_study].append(m.group(1))

                m = _RE_NF_RETRY.search(line)
                if m:
                    bad_sample = m.group(1)
                    gse = bad_sample.split("_")[0] if "_" in bad_sample else bad_sample
                    result["retries"][gse].append("NF_Retry")

                # Nextflow batch processing state
                m = _RE_BATCH_RUN.search(line)
                if m:
                    current_batch = m.group(1)
                    result["batches"][current_batch] = "running"
                    
                if _RE_BATCH_OK.search(line) and current_batch:
                    result["batches"][current_batch] = "success"

                if _RE_BATCH_ERR.search(line) and current_batch:
                    result["batches"][current_batch] = "error"

                if _RE_NF_START.search(line):
                    result["nf_runs"] += 1

                if _RE_TMP_NODEV.search(line):
                    result["nodev_error"] = True

                # Process skipped studies from top-of-log report
                m = _RE_SKIP_PROCESSED.search(line)
                if m:
                    try:
                        for gse in eval(m.group(1)):
                            if gse not in result["studies"]:
                                result["study_order"].append(gse)
                            result["studies"][gse] = "processed"
                    except Exception: pass
                
                m = _RE_SKIP_IGNORED.search(line)
                if m:
                    try:
                        for gse in eval(m.group(1)):
                            if gse not in result["studies"]:
                                result["study_order"].append(gse)
                            result["studies"][gse] = "ignore"
                    except Exception: pass

                m = _RE_SKIP_ERROR.search(line)
                if m:
                    try:
                        for gse in eval(m.group(1)):
                            if gse not in result["studies"]:
                                result["study_order"].append(gse)
                            result["studies"][gse] = "error"
                    except Exception: pass

    except Exception as e:
        result["errors"].append(f"Failed to read {path}: {e}")

    # Compute duration from timestamps
    t_start = _parse_log_timestamp(result["start_time"])
    t_end   = _parse_log_timestamp(result["end_time"])

    # If job is running, end_time defaults to NOW to calculate elapsed time
    if t_start and not result["finished"]:
        t_end = datetime.now()

    if t_start and t_end:
        result["duration_secs"] = (t_end - t_start).total_seconds()

    return result

# ──────────────────────────────────────────────────────────────────────────────
# DASHBOARD GENERATION
# ──────────────────────────────────────────────────────────────────────────────

def generate_html(all_results: list, refresh_s: int, generated_at: str, tracker_states: dict = None) -> str:
    """Generate the HTML dashboard from the parsed results."""
    if tracker_states is None:
        tracker_states = {}

    cards_html = ""
    for r in all_results:
        job_label = f'{r["job_id"]}_{r["array_index"]}' if r["array_index"] else r["job_id"]
        if not job_label:
            job_label = os.path.basename(r["path"])

        status_class = "card-done" if r["finished"] else "card-live"
        status_dot = "⚫ Finished" if r["finished"] else "🟢 Running"
        
        # Determine extra column header if tracker exists
        extra_headers = "<th>Tracker</th>" if tracker_states else ""

        # Batch rows
        batch_rows = ""
        for bid, bstate in r.get("batches", {}).items():
            b_color, b_label = BATCH_STATES.get(bstate, BATCH_STATES["unknown"])
            b_name = bid.replace("batch_", "")
            batch_rows += f'<tr><td>{b_name}</td><td><span class="state-pill" style="--c:{b_color}">{b_label}</span></td></tr>'

        errs_html = ""
        if r["errors"]:
            errs_list = "".join(f"<li>{e}</li>" for e in r["errors"][:5])
            if len(r["errors"]) > 5: errs_list += "<li>...</li>"
            errs_html = f'<div class="error-box"><ul>{errs_list}</ul></div>'

        dur_secs = r.get("duration_secs")
        if dur_secs is not None:
            dur_str = _format_duration(dur_secs)
            if r["finished"]:
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
                retry_html = f'<span class="retry-tag" title="Retry instances">{" | ".join(retries)}</span>'
                
            # Tracker column
            tk_state = tracker_states.get(gse_id)
            if tk_state:
                tk_color, tk_label = TRACKER_STATES.get(tk_state, TRACKER_STATES["unknown"])
                tracker_cell = f'<span class="state-pill" style="--c:{tk_color}">{tk_label}</span>'
            elif tracker_states:
                tracker_cell = '<span class="muted small">no file</span>'
            else:
                tracker_cell = "" # no tracker dir given

            # Discrepancy warning
            disc = disc_by_gse.get(gse_id)
            warn_icon = f'<span class="warn-icon" title="{disc["message"]}">⚠</span>' if disc else ""

            tracker_td = f"<td>{tracker_cell}</td>" if tracker_states else ""
            study_rows += f'<tr><td>{gse_id} {warn_icon} {retry_html}</td><td><span class="state-pill" style="--c:{log_color}">{log_label}</span></td>{tracker_td}</tr>'

        n_total = len(r.get("studies", []))

        cards_html += f'''
        <div class="card {status_class}">
            <div class="card-header">
                <div class="card-title">Job {job_label} {nodev_warn}</div>
                <div class="card-meta">{status_dot} • {r.get("node", "unknown node")}</div>
                <div style="margin-top:4px">{duration_html}</div>
            </div>
            
            <div class="card-cols">
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

    # Tracker overview sidebar
    tracker_sidebar_html = ""
    if tracker_states:
        tk_counts = defaultdict(int)
        for s in tracker_states.values():
            tk_counts[s] += 1
        tracker_total = len(tracker_states)
        
        tk_legend = ""
        for state, (color, label) in TRACKER_STATES.items():
            if state == "unknown": continue
            n = tk_counts[state]
            if n == 0: continue
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
        <div class="sidebar">
            <h2>FileTracker State</h2>
            <div class="tk-total">{tracker_total} studies tracked</div>
            {disc_badge}
            <div class="tk-legend">
                {tk_legend}
            </div>
            <div class="tk-hint">
                Shows ground truth across the entire project drive, not just monitored logs.
            </div>
        </div>
        '''

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
            --text:#dde3f0; --muted:#4a5577;
            --mono:'IBM Plex Mono',monospace; --sans:'IBM Plex Sans',sans-serif;
        }}
        *{{box-sizing:border-box;margin:0;padding:0}}
        body{{background:var(--bg);color:var(--text);font-family:var(--sans);font-size:14px;line-height:1.6}}
        a{{color:inherit}}
        .header{{background:var(--s1);padding:16px 24px;border-bottom:1px solid var(--border);display:flex;justify-content:space-between;align-items:center}}
        .header h1{{font-size:18px;font-weight:600}}
        .layout{{display:flex;gap:24px;padding:24px;max-width:1600px;margin:0 auto;align-items:flex-start}}
        
        /* Sidebar */
        .sidebar{{background:var(--s1);border:1px solid var(--border);border-radius:8px;padding:20px;width:300px;flex-shrink:0}}
        .sidebar h2{{font-size:14px;text-transform:uppercase;letter-spacing:1px;color:var(--muted);margin-bottom:8px}}
        .tk-total{{font-size:24px;font-weight:600;margin-bottom:16px}}
        .disc-total{{background:rgba(245,158,11,0.1);color:#f59e0b;padding:8px;border-radius:6px;font-size:12px;margin-bottom:16px;border:1px solid rgba(245,158,11,0.2);}}
        .tk-legend{{display:flex;flex-direction:column;gap:12px}}
        .legend-item{{display:flex;align-items:center;font-size:13px}}
        .legend-dot{{width:10px;height:10px;border-radius:50%;margin-right:10px}}
        .legend-label{{flex:1}}
        .legend-count{{font-family:var(--mono);color:var(--muted)}}
        .tk-bar-row{{display:flex;align-items:center;gap:8px}}
        .tk-bar{{flex:1;height:4px;background:var(--s2);border-radius:2px;overflow:hidden}}
        .tk-bar-fill{{height:100%;border-radius:2px}}
        .tk-pct{{font-family:var(--mono);font-size:11px;color:var(--muted);width:24px;text-align:right}}
        .tk-hint{{font-size:11px;color:var(--muted);margin-top:20px;line-height:1.4}}

        /* Main Content */
        .main-content{{flex:1;display:flex;flex-direction:column;gap:16px}}
        .card{{background:var(--s1);border:1px solid var(--border);border-radius:8px;overflow:hidden}}
        .card-done{{border-left:4px solid var(--muted)}}
        .card-live{{border-left:4px solid #3b82f6}}
        .card-header{{padding:16px 20px;border-bottom:1px solid var(--border);background:var(--s2)}}
        .card-title{{font-family:var(--mono);font-weight:600;font-size:15px;display:flex;align-items:center;gap:12px}}
        .card-meta{{font-size:12px;color:var(--muted);margin-top:4px}}
        .card-cols{{display:grid;grid-template-columns:minmax(300px, 1.5fr) 1fr;gap:0;}}
        .card-cols > div{{padding:20px;border-right:1px solid var(--border)}}
        .card-cols > div:last-child{{border-right:none}}
        .sub-title{{font-size:12px;text-transform:uppercase;letter-spacing:1px;color:var(--muted);margin-bottom:12px;font-weight:600}}
        
        /* Tables */
        .inner-table{{width:100%;border-collapse:collapse;font-size:13px}}
        .inner-table th{{text-align:left;padding-bottom:8px;color:var(--muted);font-weight:500;border-bottom:1px solid var(--border)}}
        .inner-table td{{padding:6px 0;border-bottom:1px solid var(--border);font-family:var(--mono)}}
        .inner-table tr:last-child td{{border-bottom:none}}
        
        .state-pill{{display:inline-block;padding:2px 8px;border-radius:12px;background:color-mix(in srgb, var(--c) 15%, transparent);color:var(--c);font-size:11px;font-family:var(--sans);font-weight:600;white-space:nowrap}}
        .retry-tag{{display:inline-block;background:#334155;color:#cbd5e1;padding:1px 6px;border-radius:4px;font-size:10px;margin-left:8px;vertical-align:middle;font-family:var(--sans)}}
        .error-box{{background:rgba(248,113,113,0.1);padding:12px 20px;border-top:1px solid rgba(248,113,113,0.2)}}
        .error-box ul{{margin:0;padding-left:20px;color:#fca5a5;font-family:var(--mono);font-size:12px}}
        .duration-chip{{display:inline-block;padding:2px 8px;border-radius:4px;font-size:11px;font-family:var(--mono)}}
        .dur-live{{background:rgba(59,130,246,0.15);color:#93c5fd}}
        .dur-done{{background:rgba(148,163,184,0.15);color:#cbd5e1}}
        .inline-warn{{display:inline-block;background:rgba(245,158,11,0.15);color:#fcd34d;padding:2px 8px;border-radius:4px;font-size:11px;margin-left:12px;font-family:var(--sans)}}
        .warn-icon{{color:#fbbf24;cursor:help;margin-left:4px;font-family:var(--sans);font-weight:bold;font-size:14px;}}
        .disc-summary{{font-size:12px;color:#fbbf24;margin-top:12px;padding-top:12px;border-top:1px dashed var(--border);}}
        .muted{{color:var(--muted)}}
        .small{{font-size:11px}}
    </style>
</head>
<body>
    <div class="header">
        <h1>RNA-seq Pipeline Monitor</h1>
        <div style="font-size:12px;color:var(--muted)">Auto-refresh: {refresh_s}s</div>
    </div>
    
    <div class="layout">
        {tracker_sidebar_html}
        <main class="main-content">
            {cards_html}
            <div style="text-align:center;color:var(--muted);font-size:12px;margin-top:20px">Last updated: {generated_at}</div>
        </main>
    </div>
</body>
</html>"""

# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────

def find_log_files(job_ids: list, log_dir: str) -> dict:
    """Find all log files for the given job IDs. Returns {task_id: path}."""
    found = {}
    for jid in job_ids:
        pattern = os.path.join(log_dir, f"stdout-*-{jid}_*.txt")
        matches = glob.glob(pattern)
        if not matches:
            pattern2 = os.path.join(log_dir, f"stdout-*{jid}*.txt")
            matches = glob.glob(pattern2)
        for path in sorted(matches):
            m = re.search(r"_(\d+)\.txt$", path)
            task_id = m.group(1) if m else os.path.basename(path)
            found[f"{jid}_{task_id}"] = path
    return found

def main():
    parser = argparse.ArgumentParser(
        description="Monitor RNA-seq SLURM job logs and generate an HTML dashboard."
    )
    parser.add_argument(
        "job_ids", nargs="+",
        help="SLURM job IDs to monitor (e.g. 19290)"
    )
    parser.add_argument(
        "--log-dir", default=os.path.join(os.getcwd(), "logs_slurm"),
        help="Directory containing SLURM log files."
    )
    parser.add_argument(
        "--tracker-dir", default=None,
        help="Path to file_tracker dir to cross-reference ground truth."
    )
    parser.add_argument("--out", default="dashboard.html", help="Output HTML file path.")
    parser.add_argument("--interval", type=int, default=20, help="Refresh interval in seconds.")
    parser.add_argument("--once", action="store_true", help="Run once and exit.")
    args = parser.parse_args()

    tracker_states = {}

    while True:
        log_files = find_log_files(args.job_ids, args.log_dir)
        
        if not log_files:
            print(f"[{datetime.now(timezone.utc).strftime('%H:%M:%S')}] No logs found for {args.job_ids} in {args.log_dir}")
        else:
            if args.tracker_dir:
                tracker_states = read_tracker_dir(args.tracker_dir)

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