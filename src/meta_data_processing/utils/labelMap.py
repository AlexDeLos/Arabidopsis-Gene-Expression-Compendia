import json
from typing import Dict, Optional


def load_json(path:str):
    with open(path, 'r') as file:
        object = json.load(file)
    return object

class LabelMap:
    """
    Manages persistent mapping of raw terms to grounded ontology terms.
    Structure: { "raw_term": ["grounded_term", ["StudyID1", "StudyID2"]] }
    """
    def __init__(self, path:Optional[str]=None):
        self.path = path
        if path is None:
            self.map_treatment = {}
            self.map_tissue = {}
            self.map = {} 
        else:
            try:
                self.map_treatment = load_json(path+'/map_treatment.json')
                self.map_tissue = load_json(path+'/map_tissue.json')
                self.map = load_json(path+'/map.json')
            except:
                print('Warning: LabelMap paths not found, starting empty.')
                self.map_treatment = {}
                self.map_tissue = {}
                self.map = {}

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

    def _get_value(self, map_dict: Dict, raw_label: str) -> Optional[str]:
        """Helper to retrieve just the string label from the complex structure."""
        if raw_label in map_dict:
            entry = map_dict[raw_label]
            if isinstance(entry, list):
                return entry[0]
            return entry # Fallback if it's still a string
        return None

    def add(self, label:str, id:str, study_id:str)->None:
        self._update_entry(self.map, label, id, study_id)

    def add_treatment(self, label:str, id:str, study_id:str)->None:
        self._update_entry(self.map_treatment, label, id, study_id)

    def add_tissue(self, label:str, id:str, study_id:str)->None:
        self._update_entry(self.map_tissue, label, id, study_id)
    
    def add_mapping_dict(self, mapping: Dict[str, str], category: str, study_id: str):
        """Batch update from LLM results."""
        target_map = None
        if category == 'treatment': target_map = self.map_treatment
        elif category == 'tissue': target_map = self.map_tissue
        elif category == 'medium': target_map = self.map
        
        if target_map is not None:
            for raw, grounded in mapping.items():
                self._update_entry(target_map, raw, grounded, study_id)

    def add_mapping(self, og_sample: Dict, grounded_sample: Dict, study_id: str) -> None:
        """
        Updates the internal maps based on the difference between raw and grounded samples.
        """
        #TODO: broke it by making everything a set
        # 1. Tissue
        raw_tissue = og_sample.get('tissue')
        if isinstance(raw_tissue, list): raw_tissue = raw_tissue[0] if raw_tissue else None
        ground_tissue = grounded_sample.get('tissue')
        
        if raw_tissue and ground_tissue:
            self.add_tissue(str(raw_tissue), str(ground_tissue), study_id)

        # 2. Medium
        raw_medium = og_sample.get('medium')
        if isinstance(raw_medium, list): raw_medium = raw_medium[0] if raw_medium else None
        ground_medium = grounded_sample.get('medium')
        
        if raw_medium and ground_medium:
            self.add(str(raw_medium), str(ground_medium), study_id)

        # 3. Treatments
        raw_treats = og_sample.get('treatment', [])
        if isinstance(raw_treats, str): raw_treats = [raw_treats]
        ground_treats = grounded_sample.get('treatment', [])
        if isinstance(ground_treats, str): ground_treats = [ground_treats]

        # Heuristic: Only map if 1-to-1 to avoid list confusion
        if len(raw_treats) == 1 and len(ground_treats) == 1:
            r_val = str(raw_treats.pop())
            g_val = str(ground_treats.pop())
            self.add_treatment(r_val, g_val, study_id)

    def save_map(self):
        if self.path:
            with open(self.path+'/map_treatment.json', 'w') as f:
                json.dump(self.map_treatment, f, indent=4)
            with open(self.path+'/map_tissue.json', 'w') as f:
                json.dump(self.map_tissue, f, indent=4)
            with open(self.path+'/map.json', 'w') as f:
                json.dump(self.map, f, indent=4)