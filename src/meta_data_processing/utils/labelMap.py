import json
from typing import Dict, Optional
import os
from src.constants import LABELS

def load_json(path:str):
    with open(path, 'r') as file:
        object = json.load(file)
    return object

class LabelMap:
    """
    Manages persistent mapping of raw terms to grounded ontology terms.
    Structure: { "raw_term": ["grounded_term", ["StudyID1", "StudyID2"]] }
    """
    def __init__(self, path: Optional[str] = None):
        self.path = path
        
        # If a path is provided but the folder doesn't exist, create it
        if path is not None and not os.path.exists(path):
            os.makedirs(path, exist_ok=True)

        # Generate specific maps dynamically based on LABELS
        for label in LABELS:
            attr_name = f"map_{label}" 
            
            if path is None:
                setattr(self, attr_name, {})
            else:
                file_path = os.path.join(path, f"{attr_name}.json")
                try:
                    # Assuming load_json is defined elsewhere in your file
                    setattr(self, attr_name, load_json(file_path))
                except Exception:
                    print(f"Warning: {file_path} not found or invalid. Generating a new one.")
                    setattr(self, attr_name, {})
                    with open(file_path, 'w') as f:
                        json.dump({}, f)

    def _update_entry(self, map_dict: Dict, raw_label: str, grounded_label: str, study_id: str):
        """Helper to safely update the mapping structure with study tracking."""
        if raw_label not in map_dict:
            # Create new entry: [Value, [StudyID]]
            map_dict[raw_label] = [grounded_label, [study_id]]
        else:
            entry = map_dict[raw_label]
            # Migration handling: Convert old string format to new list format if needed
            if isinstance(entry, str):
                entry = [entry, []]
            
            # entry[0] is the label. We keep the existing label (first winner stays).
            # entry[1] is the list of studies.
            if study_id and study_id not in entry[1]:
                entry[1].append(study_id)
            
            map_dict[raw_label] = entry

    def _get_value(self, category: str, raw_label: str) -> Optional[str]:
        """Helper to retrieve just the string label dynamically."""
        map_dict = getattr(self, f"map_{category}", None)
        if map_dict is not None and raw_label in map_dict:
            entry = map_dict[raw_label]
            if isinstance(entry, list):
                return entry[0]
            return entry # Fallback if it's still a string
        return None

    def add_entry(self, category: str, raw_label: str, grounded_label: str, study_id: str) -> None:
        """Dynamically add an entry to the correct map."""
        map_dict = getattr(self, f"map_{category}", None)
        if map_dict is not None:
            self._update_entry(map_dict, raw_label, grounded_label, study_id)
        else:
            print(f"Warning: Category '{category}' is not a valid label map.")

    def add_mapping_dict(self, mapping: Dict[str, str], category: str, study_id: str):
        """Batch update from LLM results using dynamic fetching."""
        map_dict = getattr(self, f"map_{category}", None)
        if map_dict is not None:
            for raw, grounded in mapping.items():
                self._update_entry(map_dict, raw, grounded, study_id)
        else:
            print(f"Warning: Category '{category}' is not a valid label map.")

    def add_mapping(self, og_sample: Dict, grounded_sample: Dict, study_id: str) -> None:
        """
        Updates the internal maps based on the difference between raw and grounded samples.
        Dynamically loops through all LABELS.
        """
        for category in LABELS:
            raw_val = og_sample.get(category)
            ground_val = grounded_sample.get(category)

            # Fix for `#TODO: broke it by making everything a set`
            if isinstance(raw_val, set): raw_val = list(raw_val)
            if isinstance(ground_val, set): ground_val = list(ground_val)

            # Standardize to lists so we can reliably check lengths
            if isinstance(raw_val, str): raw_val = [raw_val]
            if isinstance(ground_val, str): ground_val = [ground_val]
            
            # Heuristic: Only map if 1-to-1 to avoid list confusion
            if raw_val and ground_val and len(raw_val) == 1 and len(ground_val) == 1:
                r_str = str(raw_val[0])
                g_str = str(ground_val[0])
                self.add_entry(category, r_str, g_str, study_id)

    def save_map(self):
        """Dynamically saves all maps based on LABELS."""
        if self.path:
            for label in LABELS:
                file_path = os.path.join(self.path, f"map_{label}.json")
                map_dict = getattr(self, f"map_{label}", {})
                with open(file_path, 'w') as f:
                    json.dump(map_dict, f, indent=4)