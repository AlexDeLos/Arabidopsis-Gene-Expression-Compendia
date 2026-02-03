
import json
import os

class RNASeq_tracker:
    def __init__(self,location) -> None:
        self.platform_counts: dict = {}
        self.loc = location
        self.totals: dict = {
            'total_studies_seen': 0, 'total_sample_seen': 0,
            'total_samples_used': 0, 'total_studies_used': 0
        }
        self.states: dict = {'ignore': set(), 'downloaded': set(), 'processed': set()}
        
        # NEW STRUCTURE: { "Platform_Name": { "12": 5, "24": 1 } }
        self.sample_distribution_per_platform: dict = {}

    def update_platform(self, platform: str, samples: int, has_raw: bool):
        if platform not in self.platform_counts:
            self.platform_counts[platform] = {
                'studies_seen': 0, 'samples_seen': 0,
                'studies_with_raw': 0, 'samples_with_raw': 0
            }
        
        # --- NEW: Update Per-Platform Histogram ---
        if platform not in self.sample_distribution_per_platform:
            self.sample_distribution_per_platform[platform] = {}
            
        # Use string key for JSON compatibility
        s_str = str(samples)
        if s_str in self.sample_distribution_per_platform[platform]:
            self.sample_distribution_per_platform[platform][s_str] += 1
        else:
            self.sample_distribution_per_platform[platform][s_str] = 1

        # Update Totals
        self.totals['total_studies_seen'] += 1
        self.totals['total_sample_seen'] += samples
        self.platform_counts[platform]['studies_seen'] += 1
        self.platform_counts[platform]['samples_seen'] += samples
        if has_raw:
            self.totals['total_studies_used'] += 1
            self.totals['total_samples_used'] += samples
            self.platform_counts[platform]['studies_with_raw'] += 1
            self.platform_counts[platform]['samples_with_raw'] += samples

    def mark_ignore(self, gse_id):
        self.states['ignore'].add(gse_id); self.states['downloaded'].discard(gse_id); self.states['processed'].discard(gse_id)
        self.save_to_json(self.loc)
    def mark_downloaded(self, gse_id):
        self.states['downloaded'].add(gse_id); self.states['ignore'].discard(gse_id); self.states['processed'].discard(gse_id)
        self.save_to_json(self.loc)
    def mark_processed(self, gse_id):
        self.states['processed'].add(gse_id); self.states['downloaded'].discard(gse_id); self.states['ignore'].discard(gse_id)
        self.save_to_json(self.loc)
    def is_ignored(self, gse_id): return gse_id in self.states['ignore']
    def is_processed(self, gse_id): return gse_id in self.states['processed']
    def is_downloaded(self, gse_id): return gse_id in self.states['downloaded']

    def save_to_json(self, filename):
        serializable_states = {k: list(v) for k, v in self.states.items()}
        data = {
            "totals": self.totals, 
            "platform_counts": self.platform_counts, 
            "states": serializable_states,
            "sample_distribution_per_platform": self.sample_distribution_per_platform
        }
        with open(filename, 'w') as f: json.dump(data, f, indent=4)
    
    @classmethod
    def load_from_json(cls, filename="rnaseq_tracker_results.json"):
        if not os.path.exists(filename): return cls()
        with open(filename, 'r') as f: data = json.load(f)
        tracker = cls()
        tracker.totals = data.get("totals", tracker.totals)
        tracker.platform_counts = data.get("platform_counts", {})
        
        loaded_states = data.get("states", {})
        tracker.states = {k: set(loaded_states.get(k, [])) for k in ['ignore', 'downloaded', 'processed']}
        
        # Load the new per-platform distribution
        tracker.sample_distribution_per_platform = data.get("sample_distribution_per_platform", {})
        
        return tracker