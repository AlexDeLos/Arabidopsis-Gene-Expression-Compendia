import os
import json
import time
from datetime import datetime, timedelta
import sys

module_dir = './'
sys.path.append(module_dir)

from src.constants import *

def monitor_tracker(tracker_dir):
    if not os.path.exists(tracker_dir):
        print(f"Error: Tracker directory not found: {tracker_dir}")
        return

    print(f"\n{'='*60}")
    print(f"  CLUSTER MONITORING DASHBOARD")
    print(f"  Source: {tracker_dir}")
    print(f"  Time:   {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'='*60}\n")

    # --- COUNTERS ---
    stats = {
        STATUS_NOT_TRIED: 0,
        STATUS_DOWNLOADED: 0,
        STATUS_PROCESSED: 0,
        STATUS_IGNORE: 0,
        STATUS_ERROR: 0
    }
    
    platform_stats = {} # {"GPL123": {"studies": 0, "samples": 0}}
    errors = []
    
    # --- SPEED METRICS ---
    now = time.time()
    files_last_hour = 0
    files_last_24h = 0
    
    # Get all files
    try:
        all_files = os.listdir(tracker_dir)
    except Exception as e:
        print(f"Error reading directory: {e}")
        return

    print(f"Scanning {len(all_files)} files... please wait...\n")

    # --- SCAN LOOP ---
    for filename in all_files:
        filepath = os.path.join(tracker_dir, filename)
        
        # 1. Process Status Files (.txt)
        if filename.endswith(".txt"):
            try:
                # Get modification time for speed calculation
                mtime = os.path.getmtime(filepath)
                if now - mtime < 3600: # 1 hour
                    files_last_hour += 1
                if now - mtime < 86400: # 24 hours
                    files_last_24h += 1

                with open(filepath, 'r') as f:
                    code = int(f.read().strip())
                    
                if code in stats:
                    stats[code] += 1
                
                if code == STATUS_ERROR:
                    errors.append(filename.replace(".txt", ""))

            except:
                pass # Ignore unreadable files

        # 2. Process Metadata Files (.json) - For Platform stats
        elif filename.endswith("_meta.json"):
            try:
                with open(filepath, 'r') as f:
                    data = json.load(f)
                    plat = data.get("platform", "Unknown")
                    samps = data.get("num_samples", 0)
                    
                    if plat not in platform_stats:
                        platform_stats[plat] = {"studies": 0, "samples": 0}
                    
                    platform_stats[plat]["studies"] += 1
                    platform_stats[plat]["samples"] += samps
            except:
                pass

    # --- REPORTING ---
    total_processed = sum(stats.values())
    
    # 1. STATUS BARS
    print("--- JOB STATUS ---")
    print(f"TOTAL FILES SEEN: {total_processed}")
    
    # Helper for progress bar
    def print_bar(label, count, total, color_code=""):
        percent = (count / total * 100) if total > 0 else 0
        bar_len = 30
        filled = int(bar_len * percent / 100)
        bar = "█" * filled + "-" * (bar_len - filled)
        print(f"{label:<15} |{bar}| {count:>5} ({percent:>5.1f}%)")

    print_bar("SUCCESS (2)", stats[STATUS_PROCESSED], total_processed)
    print_bar("IGNORED (3)", stats[STATUS_IGNORE], total_processed)
    print_bar("ERRORS  (4)", stats[STATUS_ERROR], total_processed)
    print_bar("DOWNLD  (1)", stats[STATUS_DOWNLOADED], total_processed)
    print_bar("PENDING (0)", stats[STATUS_NOT_TRIED], total_processed)

    # 2. SPEED ESTIMATION
    print(f"\n--- CLUSTER SPEED ---")
    print(f"Activity (Last 60 mins): {files_last_hour} files updated/created")
    print(f"Activity (Last 24 hrs):  {files_last_24h} files updated/created")
    if files_last_hour > 0:
        est_daily = files_last_hour * 24
        print(f"Est. Daily Pace:         ~{est_daily} studies/day")
    else:
        print("Status:                  IDLE (No activity in last hour)")

    # 3. PLATFORM STATISTICS
    print(f"\n--- TOP PLATFORMS (By Sample Count) ---")
    print(f"{'PLATFORM':<15} {'STUDIES':<10} {'SAMPLES':<10}")
    print("-" * 40)
    
    # Sort by sample count descending
    sorted_plats = sorted(platform_stats.items(), key=lambda item: item[1]['samples'], reverse=True)
    
    for plat, data in sorted_plats[:5]: # Show top 5
        print(f"{plat:<15} {data['studies']:<10} {data['samples']:<10}")

    # 4. ERROR LOG
    if errors:
        print(f"\n--- FAILED STUDIES ({len(errors)}) ---")
        print(f"First 10 failures: {', '.join(errors[:10])}")
        if len(errors) > 10:
            print(f"...and {len(errors) - 10} more.")
        print("To see all errors, run: grep -r '4' path/to/tracker")
    else:
        print("\n--- NO ERRORS DETECTED ---")

    print(f"\n{'='*60}\n")

if __name__ == "__main__":
    # CONFIGURATION
    # Update this path to match your actual tracker storage location
    ROOT_STORAGE = "./new_storage/"
    TRACKER_DIR = os.path.join(ROOT_STORAGE, "rnaseq_data/file_tracker")
    
    monitor_tracker(TRACKER_DIR)