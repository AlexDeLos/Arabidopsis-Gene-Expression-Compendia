import os
import json
import sys
module_dir = './'
sys.path.append(module_dir)

from src.constants import *

# Ensure these match your constants file
class FileTracker:
    def __init__(self, root_dir):
        """
        Initializes the file-based tracker.
        :param root_dir: Path to the 'tracker_status' directory.
        """
        self.tracker_dir = root_dir
        os.makedirs(self.tracker_dir, exist_ok=True)
        
        # Compatibility: These are technically useless in parallel jobs, 
        # but we keep them so your code doesn't crash if it tries to access them.
        self.totals = {
            'total_studies_seen': 0, 'total_sample_seen': 0,
            'total_samples_used': 0, 'total_studies_used': 0
        }
        # In parallel mode, we cannot maintain 'states' sets in memory.
        # We will generate them dynamically if requested.
        self.states = {'ignore': set(), 'downloaded': set(), 'processed': set()}

    def _get_file_path(self, gse_id, ext=".txt"):
        return os.path.join(self.tracker_dir, f"{gse_id}{ext}")

    # --- CORE STATUS METHODS ---
    def get_status(self, gse_id):
        path = self._get_file_path(gse_id)
        if not os.path.exists(path):
            return STATUS_NOT_TRIED
        try:
            with open(path, 'r') as f:
                return int(f.read().strip())
        except:
            return STATUS_ERROR

    def set_status(self, gse_id, status_code):
        path = self._get_file_path(gse_id)
        with open(path, 'w') as f:
            f.write(str(status_code))

    # --- MISSING BOOLEAN CHECKS ---
    def is_processed(self, gse_id):
        return self.get_status(gse_id) == STATUS_PROCESSED

    def is_downloaded(self, gse_id):
        # Checks if status is Downloaded OR Processed (since processed implies downloaded)
        s = self.get_status(gse_id)
        return s == STATUS_DOWNLOADED or s == STATUS_PROCESSED

    def is_ignored(self, gse_id):
        return self.get_status(gse_id) == STATUS_IGNORE

    # --- MARKER METHODS ---
    def mark_processed(self, gse_id):
        self.set_status(gse_id, STATUS_PROCESSED)

    def mark_downloaded(self, gse_id):
        self.set_status(gse_id, STATUS_DOWNLOADED)

    def mark_ignore(self, gse_id):
        self.set_status(gse_id, STATUS_IGNORE)

    def mark_error(self, gse_id):
        self.set_status(gse_id, STATUS_ERROR)

    # --- ADAPTED PLATFORM LOGIC (The Distributed Fix) ---
    def update_platform(self, platform, num_samples, has_raw):
        """
        COMPATIBILITY METHOD:
        In the old tracker, this updated a global dictionary.
        In the new tracker, this does nothing because memory is not shared.
        
        To actually save this data, use 'save_study_metadata' below.
        """
        pass 

    def save_study_metadata(self, gse_id, platform, num_samples, has_raw):
        """
        New method to persist stats per-study.
        Writes to: tracker_dir/GSE123_meta.json
        """
        data = {
            "gse_id": gse_id,
            "platform": platform,
            "num_samples": num_samples,
            "has_raw": has_raw
        }
        path = self._get_file_path(gse_id, ext="_meta.json")
        with open(path, 'w') as f:
            json.dump(data, f)

    # --- DUMMY METHODS (For Compatibility) ---
    def save_to_json(self, path=None):
        pass # Auto-save happens instantly on set_status

    @classmethod
    def load_from_json(cls, path):
        # If your code tries to load from a JSON path, we simply 
        # redirect it to open the directory instead.
        # We assume 'path' ends in '.../tracker.json', so we strip the filename
        # to get the directory.
        if os.path.isdir(path):
            return cls(path)
        else:
            # Fallback: use the parent directory of the json file
            return cls(os.path.dirname(path) + "/tracker_status")

    # --- REPORTING (Replaces the old 'totals' dict) ---
    def generate_detailed_report(self):
        """
        Scans all files to rebuild the 'platform_counts' and 'totals'
        dictionaries that used to exist in memory.
        """
        stats = {
            "statuses": {0:0, 1:0, 2:0, 3:0, 4:0},
            "platforms": {},
            "total_samples": 0
        }
        
        print("Generating report from distributed files...")
        files = os.listdir(self.tracker_dir)
        
        for f in files:
            full_path = os.path.join(self.tracker_dir, f)
            
            # Count Statuses (.txt files)
            if f.endswith(".txt"):
                try:
                    with open(full_path, 'r') as file:
                        code = int(file.read().strip())
                        stats["statuses"][code] += 1
                except: pass
            
            # Aggregate Metadata (.json files)
            elif f.endswith("_meta.json"):
                try:
                    with open(full_path, 'r') as file:
                        data = json.load(file)
                        plat = data.get("platform", "Unknown")
                        samps = data.get("num_samples", 0)
                        
                        if plat not in stats["platforms"]:
                            stats["platforms"][plat] = {"studies": 0, "samples": 0}
                        
                        stats["platforms"][plat]["studies"] += 1
                        stats["platforms"][plat]["samples"] += samps
                        stats["total_samples"] += samps
                except: pass
                
        return stats