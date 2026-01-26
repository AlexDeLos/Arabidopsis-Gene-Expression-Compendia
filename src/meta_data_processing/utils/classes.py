import json
import spacy
from typing import List, Dict, Optional, Tuple
from sentence_transformers import SentenceTransformer, util
from spacy.cli import download

# Import constants directly
import sys
module_dir = './'
sys.path.append(module_dir)
from src.constants import VALID_TREATMENTS, VALID_TISSUES, VALID_MEDIUMS, VALID_TREATMENTS_ALT

def load_json(path:str):
    with open(path, 'r') as file:
        object = json.load(file)
    return object

class LabelMap:
    """
    Manages persistent mapping of raw terms to grounded ontology terms.
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

    def add(self, label:str, id)->None:
        self.map[label] = id

    def add_treatment(self, label:str, id)->None:
        self.map_treatment[label] = id

    def add_tissue(self, label:str, id)->None:
        self.map_tissue[label] = id

    def add_mapping_dict(self, mapping: Dict[str, str], category: str):
        """Batch update the map from LLM results"""
        if category == 'treatment':
            self.map_treatment.update(mapping)
        elif category == 'tissue':
            self.map_tissue.update(mapping)
        elif category == 'medium':
            self.map.update(mapping)

    def save_map(self):
        if self.path:
            with open(self.path+'/map_treatment.json', 'w') as f:
                json.dump(self.map_treatment, f, indent=4)
            with open(self.path+'/map_tissue.json', 'w') as f:
                json.dump(self.map_tissue, f, indent=4)
            with open(self.path+'/map.json', 'w') as f:
                json.dump(self.map, f, indent=4)

    def add_mapping(self, og_sample: Dict, grounded_sample: Dict) -> None:
        """
        Updates the internal maps based on the difference between raw and grounded samples.
        """
        
        # --- 1. Handle Tissue (One-to-One) ---
        raw_tissue = og_sample.get('tissue')
        # Normalize: Extractor might return list or string
        if isinstance(raw_tissue, list): 
            raw_tissue = raw_tissue[0] if raw_tissue else None
            
        ground_tissue = grounded_sample.get('tissue')
        
        # Only map if we have valid strings and they aren't identical
        if raw_tissue and ground_tissue:
            raw_str = str(raw_tissue)
            ground_str = str(ground_tissue)
            self.map_tissue[raw_str] = ground_str


        # --- 2. Handle Medium (One-to-One) ---
        raw_medium = og_sample.get('medium')
        if isinstance(raw_medium, list): 
            raw_medium = raw_medium[0] if raw_medium else None
            
        ground_medium = grounded_sample.get('medium')
        
        if raw_medium and ground_medium:
            raw_str = str(raw_medium)
            ground_str = str(ground_medium)
            self.map[raw_str] = ground_str


        # --- 3. Handle Treatments (Many-to-Many) ---
        raw_treats = og_sample.get('treatment', [])
        if isinstance(raw_treats, str): raw_treats = [raw_treats]
        
        ground_treats = grounded_sample.get('treatment', [])
        if isinstance(ground_treats, str): ground_treats = [ground_treats]

        # HEURISTIC: Safety Check
        # We can only safely infer a direct mapping from the final output lists 
        # if there is exactly ONE item.
        # If there are multiple (e.g. Raw=["A", "B"] -> Ground=["Y", "Z"]), 
        # we cannot know if A->Y or A->Z because lists might be re-sorted.
        if len(raw_treats) == 1 and len(ground_treats) == 1:
            r_val = str(raw_treats[0])
            g_val = str(ground_treats[0])
            
            # Don't map specific terms to generic "Other" catch-alls permanently
            if "Other" not in g_val:
                self.map_treatment[r_val] = g_val
        return

class GroundingOptimizer:
    def __init__(self):
        print("Initializing GroundingOptimizer...")
        self.model = SentenceTransformer('FremyCompany/BioLORD-2023')
        
        try:
            self.nlp = spacy.load("en_core_web_sm", disable=['parser', 'ner'])
        except OSError:
            from spacy.cli import download
            download("en_core_web_sm")
            self.nlp = spacy.load("en_core_web_sm", disable=['parser', 'ner'])
        
        self.valid_treatments = VALID_TREATMENTS
        self.valid_treatments_alt = VALID_TREATMENTS_ALT 
        self.valid_tissues = VALID_TISSUES
        self.valid_mediums = VALID_MEDIUMS

        print("Pre-computing ontology vectors...")
        self.treatment_vecs = self._encode_ontology(self.valid_treatments)
        self.treatment_vecs_alt = self._encode_ontology(self.valid_treatments_alt)
        self.tissue_vecs = self._encode_ontology(self.valid_tissues)
        self.medium_vecs = self._encode_ontology(self.valid_mediums)
        
    def _lemmatize(self, text: str) -> str:
        if not text: return ""
        doc = self.nlp(str(text).lower())
        return " ".join([token.lemma_ for token in doc])

    def _encode_ontology(self, terms: list):
        lemmatized_terms = [self._lemmatize(t) for t in terms]
        return self.model.encode(lemmatized_terms)

    # --- NEW: Heuristic Logic ---
    def _check_heuristics(self, text: str, category: str) -> Optional[str]:
        """
        Fast checks to catch obvious items or extract keywords embedded in noise.
        """
        text_lower = text.lower().strip()
        
        # 1. STOPLIST: Handle generic/noise terms immediately
        # "conditions" -> No stress
        IGNORED_TERMS = {
            'conditions', 'condition', 'time', 'date', 'replicate', 'rep', 
            'genotype', 'analysis', 'samples', 'study', 'experiment'
        }
        if text_lower in IGNORED_TERMS:
            if category == 'treatment': return "No stress"
            if category == 'tissue': return "unknown" 
            return None

        # 2. SUBSTRING RESCUE (Reverse Lookup)
        # Check if a known valid label exists inside the messy string.
        # e.g. "transcripts expressed in the phloem tissue" -> contains "phloem"
        
        if category == 'tissue':
            # Iterate through known tissues (sorted by length to match "guard cell" before "cell")
            for tissue in sorted(self.valid_tissues, key=len, reverse=True):
                # We use word boundary check or simple inclusion?
                # Simple inclusion is risky for short words (e.g. 'bud' in 'budget'), 
                # but safe for biological terms if we are careful.
                if tissue.lower() in text_lower:
                    return tissue
                    
        elif category == 'treatment':
            # Check Alternative names (e.g. "Drought")
            # If "Drought" is in text, return "Drought Stress"
            for i, alt_name in enumerate(self.valid_treatments_alt):
                if alt_name.lower() in text_lower:
                    return self.valid_treatments[i] # Return the canonical name
            
            # Special manual overrides
            if 'mock' in text_lower or 'control' in text_lower:
                return "No stress"

        elif category == 'medium':
            for medium in sorted(self.valid_mediums, key=len, reverse=True):
                if medium.lower() in text_lower:
                    return medium
                
        return None

    def find_semantic_match(self, raw_term: str, category: str, threshold: float = 0.88) -> Optional[str]:
        if not raw_term: return None
        
        # 1. Try Heuristics First (Fast & fixes "phloem tissue")
        heuristic_match = self._check_heuristics(raw_term, category)
        if heuristic_match:
            return heuristic_match
        
        # 2. Try Vector Search (Smart but needs similarity)
        label, score = self.get_best_match_with_score(raw_term, category)
        if score > threshold:
            return label
            
        return None

    def get_best_match_with_score(self, text: str, category: str) -> Tuple[Optional[str], float]:
        # (Same implementation as previous step)
        if not text: return None, 0.0

        clean_text = self._lemmatize(text)
        query_vec = self.model.encode([clean_text])[0]

        if category == 'treatment':
            scores = util.cos_sim(query_vec, self.treatment_vecs)[0]
            best_idx = scores.argmax().item()
            best_score = scores[best_idx].item()
            best_label = self.valid_treatments[best_idx]

            scores_alt = util.cos_sim(query_vec, self.treatment_vecs_alt)[0]
            idx_alt = scores_alt.argmax().item()
            if scores_alt[idx_alt].item() > best_score:
                best_score = scores_alt[idx_alt].item()
                best_label = self.valid_treatments[idx_alt]
            return best_label, best_score

        elif category == 'tissue':
            vectors = self.tissue_vecs
            labels = self.valid_tissues
        elif category == 'medium':
            vectors = self.medium_vecs
            labels = self.valid_mediums
        else:
            return None, 0.0

        scores = util.cos_sim(query_vec, vectors)[0]
        best_idx = scores.argmax().item()
        return labels[best_idx], scores[best_idx].item()

    def batch_process_study(self, data: Dict, extracted_samples: List[Dict], llm_func_treat, llm_func_tis, label_map: LabelMap):
        # (Same implementation as previous step, ensuring find_semantic_match is called)
        local_cache = { 'treatment': {}, 'tissue': {}, 'medium': {} }
        unknowns = { 'treatment': set(), 'tissue': set(), 'medium': set() }

        # --- PASS 1: Try Maps, Heuristics, and Vectors ---
        for sample in extracted_samples:
            
            # Treatments
            raw_treats = sample.get('treatment', [])
            if isinstance(raw_treats, str): raw_treats = [raw_treats]
            for t in raw_treats:
                if t in local_cache['treatment']: continue
                if t in label_map.map_treatment:
                    local_cache['treatment'][t] = label_map.map_treatment[t]
                    continue
                
                match = self.find_semantic_match(t, 'treatment', threshold=0.88)
                if match:
                    local_cache['treatment'][t] = match
                    label_map.add_treatment(t, match)
                else:
                    unknowns['treatment'].add(t)

            # Tissue
            raw_tissue = sample.get('tissue', 'unknown')
            if isinstance(raw_tissue, list): raw_tissue = raw_tissue[0] if raw_tissue else "unknown"
            
            if raw_tissue not in local_cache['tissue']:
                if raw_tissue in label_map.map_tissue:
                    local_cache['tissue'][raw_tissue] = label_map.map_tissue[raw_tissue]
                else:
                    match = self.find_semantic_match(raw_tissue, 'tissue', threshold=0.88)
                    if match:
                        local_cache['tissue'][raw_tissue] = match
                        label_map.add_tissue(raw_tissue, match)
                    else:
                        unknowns['tissue'].add(raw_tissue)

            # Medium
            raw_medium = sample.get('medium', 'unspecified')
            if isinstance(raw_medium, list): raw_medium = raw_medium[0] if raw_medium else "unspecified"
            
            if raw_medium not in local_cache['medium']:
                if raw_medium in label_map.map:
                    local_cache['medium'][raw_medium] = label_map.map[raw_medium]
                else:
                    match = self.find_semantic_match(raw_medium, 'medium', threshold=0.80)
                    if match:
                        local_cache['medium'][raw_medium] = match
                        label_map.add(raw_medium, match)
                    else:
                        # Medium usually doesn't go to LLM
                        local_cache['medium'][raw_medium] = "unspecified" 

        # --- PASS 2: Batch LLM Calls ---
        study_context = data.get('study_metadata', {})
        
        if unknowns['treatment']:
            print(f"  > Sending {len(unknowns['treatment'])} treatments to LLM...")
            try:
                results = None#llm_func_treat(list(unknowns['treatment']), study_context)
                if results:
                    local_cache['treatment'].update(results)
                    label_map.add_mapping_dict(results, 'treatment')
            except Exception as e:
                print(f"    LLM Error (Treatment): {e}")

        if unknowns['tissue']:
            print(f"  > Sending {len(unknowns['tissue'])} tissues to LLM...")
            try:
                results = None#llm_func_tis(list(unknowns['tissue']), study_context)
                if results:
                    local_cache['tissue'].update(results)
                    label_map.add_mapping_dict(results, 'tissue')
            except Exception as e:
                print(f"    LLM Error (Tissue): {e}")

        # --- PASS 3: Apply ---
        final_samples = []
        for sample in extracted_samples:
            new_sample = sample.copy()
            
            # Treatment
            raw_treats = sample.get('treatment', [])
            if isinstance(raw_treats, str): raw_treats = [raw_treats]
            final_treats = []
            for t in raw_treats:
                val = local_cache['treatment'].get(t, "Other stress")
                final_treats.append(val)
            new_sample['treatment'] = sorted(list(set(final_treats)))
            
            # Tissue
            raw_tissue = sample.get('tissue')
            if isinstance(raw_tissue, list): raw_tissue = raw_tissue[0] if raw_tissue else "unknown"
            new_sample['tissue'] = local_cache['tissue'].get(raw_tissue, "unknown")
            
            # Medium
            raw_medium = sample.get('medium')
            if isinstance(raw_medium, list): raw_medium = raw_medium[0] if raw_medium else "unspecified"
            new_sample['medium'] = local_cache['medium'].get(raw_medium, "unspecified")

            final_samples.append(new_sample)

        return final_samples

## OLD
class GroundingOptimizer_old:
    def __init__(self):
        print("Initializing GroundingOptimizer...")
        
        # 1. LOAD MODEL
        print("Loading BioLORD model...")
        self.model = SentenceTransformer('FremyCompany/BioLORD-2023')
        
        # 2. LOAD LEMMATIZER
        try:
            self.nlp = spacy.load("en_core_web_sm", disable=['parser', 'ner'])
        except OSError:
            print("Downloading spacy model...")
            download("en_core_web_sm")
            self.nlp = spacy.load("en_core_web_sm", disable=['parser', 'ner'])
        
        # 3. LOAD ONTOLOGIES
        self.valid_treatments = VALID_TREATMENTS
        # New: Alternative treatments list (e.g. "Drought" vs "Drought Stress")
        self.valid_treatments_alt = VALID_TREATMENTS_ALT 
        
        self.valid_tissues = VALID_TISSUES
        self.valid_mediums = VALID_MEDIUMS

        # 4. PRE-COMPUTE VECTORS
        print("Pre-computing ontology vectors...")
        self.treatment_ontology_vectors = self._encode_ontology(self.valid_treatments)
        # New: Encode the alternative list
        self.treatment_ontology_vectors_alt = self._encode_ontology(self.valid_treatments_alt)
        
        self.tissue_ontology_vectors = self._encode_ontology(self.valid_tissues)
        self.medium_ontology_vectors = self._encode_ontology(self.valid_mediums)
        
    def _lemmatize(self, text: str) -> str:
        if not text: return ""
        doc = self.nlp(str(text).lower())
        return " ".join([token.lemma_ for token in doc])

    def _encode_ontology(self, terms: list):
        lemmatized_terms = [self._lemmatize(t) for t in terms]
        return self.model.encode(lemmatized_terms)

    def find_semantic_match(self, raw_term: str, ontology: str, threshold: float = 0.85) -> Optional[str]:
        if not raw_term:
            return None
        label, score = self.get_best_match_with_score(raw_term, ontology)
        if score > threshold:
            return label
        return None

    def get_best_match_with_score(self, text: str, category: str) -> Tuple[Optional[str], float]:
        """
        Calculates similarity of text against the specified ontology category.
        For Treatments, checks both Canonical and Alt lists and returns the best score.
        """
        if not text:
            return None, 0.0

        # CRITICAL: Lemmatize input before embedding
        clean_text = self._lemmatize(text)
        query_vec = self.model.encode([clean_text])[0]

        best_label = None
        best_score = 0.0

        if category == 'treatment':
            # Check 1: Canonical List ("Drought Stress")
            scores_main = util.cos_sim(query_vec, self.treatment_ontology_vectors)[0]
            idx_main = scores_main.argmax().item()
            score_main = scores_main[idx_main].item()

            # Check 2: Alternative List ("Drought")
            scores_alt = util.cos_sim(query_vec, self.treatment_ontology_vectors_alt)[0]
            idx_alt = scores_alt.argmax().item()
            score_alt = scores_alt[idx_alt].item()

            # Logic: Compare scores, but map index back to Canonical list
            # We assume indices are perfectly aligned (i.e., VALID_TREATMENTS[0] corresponds to VALID_TREATMENTS_ALT[0])
            if score_alt > score_main:
                best_score = score_alt
                # Return the canonical name even if the match was found in the alt list
                best_label = self.valid_treatments[idx_alt]
            else:
                best_score = score_main
                best_label = self.valid_treatments[idx_main]

            return best_label, best_score

        elif category == 'tissue':
            vectors = self.tissue_ontology_vectors
            labels = self.valid_tissues
        elif category == 'medium':
            vectors = self.medium_ontology_vectors
            labels = self.valid_mediums
        else:
            return None, 0.0

        # Standard processing for Tissue/Medium
        scores = util.cos_sim(query_vec, vectors)[0]
        best_idx = scores.argmax().item()
        best_score = scores[best_idx].item()
        best_label = labels[best_idx]
        
        return best_label, best_score

    def batch_process_study(self, data: Dict, extracted_samples: List[Dict], llm_func, llm_func_tis, label_map_obj: LabelMap):
        # (Same implementation as provided previously)
        mapping_cache_tissue = {}
        mapping_cache_treatment = {}
        mapping_cache_medium = {}
        
        unknown_treatments_for_llm = set()
        unknown_tissues_for_llm = set()

        def resolve_term(term, specific_map, ontology_type):
            if not term: return None
            if specific_map and term in specific_map:
                return specific_map[term]
            match = self.find_semantic_match(term, ontology_type, threshold=0.85)
            if match:
                return match
            return None

        # --- STEP 1: Process Fields locally ---
        for sample in extracted_samples:
            # A. Process TREATMENT
            treatments = sample.get('treatment', [])
            if isinstance(treatments, str): treatments = [treatments]
            
            for t in treatments:
                if t in mapping_cache_treatment: continue 
                resolved = resolve_term(t, label_map_obj.map_treatment, 'treatment')
                if resolved:
                    mapping_cache_treatment[t] = resolved
                else:
                    unknown_treatments_for_llm.add(t)

            # B. Process TISSUE
            tissue = sample.get('tissue')
            if isinstance(tissue, list): tissue = tissue[0] if tissue else "unknown"
            if tissue and tissue not in mapping_cache_tissue:
                resolved_tis = resolve_term(tissue, label_map_obj.map_tissue, 'tissue')
                if resolved_tis:
                    mapping_cache_tissue[tissue] = resolved_tis
                else:
                    unknown_tissues_for_llm.add(tissue)

            # C. Process MEDIUM
            medium = sample.get('medium')
            if isinstance(medium, list): medium = medium[0] if medium else "unspecified"
            if medium and medium not in mapping_cache_medium:
                resolved_m = resolve_term(medium, label_map_obj.map, 'medium')
                if resolved_m:
                    mapping_cache_medium[medium] = resolved_m
                else:
                    mapping_cache_medium[medium] = medium

        # --- STEP 2: Batch LLM Calls (Fallback) ---
        data_study = data.get('study_metadata', {})
        
        if unknown_treatments_for_llm:
            print(f"Querying LLM for {len(unknown_treatments_for_llm)} unique treatments...")
            try:
                llm_result = llm_func(list(unknown_treatments_for_llm), data_study)
                if llm_result: mapping_cache_treatment.update(llm_result)
            except Exception as e:
                print(f"LLM Batch Error (Treatments): {e}")
        
        if unknown_tissues_for_llm:
            print(f"Querying LLM for {len(unknown_tissues_for_llm)} unique tissues...")
            try:
                llm_result = llm_func_tis(list(unknown_tissues_for_llm), data_study)
                if llm_result: mapping_cache_tissue.update(llm_result)
            except Exception as e:
                print(f"LLM Batch Error (Tissues): {e}")

        # --- STEP 3: Apply mappings and Return ---
        final_samples = []
        for sample in extracted_samples:
            new_sample = sample.copy()
            
            # Treatments
            new_treatments = []
            raw_treats = sample.get('treatment', [])
            if isinstance(raw_treats, str): raw_treats = [raw_treats]
            
            for t in raw_treats:
                val = mapping_cache_treatment.get(t)
                if not val and t in unknown_treatments_for_llm:
                    val = "Other stress"
                new_treatments.append(val if val else "Other stress")
            new_sample['treatment'] = sorted(list(set(new_treatments)))
            
            # Tissue
            raw_tissue = sample.get('tissue')
            if isinstance(raw_tissue, list): raw_tissue = raw_tissue[0] if raw_tissue else "unknown"
            new_sample['tissue'] = mapping_cache_tissue.get(raw_tissue, "unknown")
                
            # Medium
            raw_medium = sample.get('medium')
            if isinstance(raw_medium, list): raw_medium = raw_medium[0] if raw_medium else "unspecified"
            new_sample['medium'] = mapping_cache_medium.get(raw_medium, "unspecified")

            final_samples.append(new_sample)
            
        return final_samples