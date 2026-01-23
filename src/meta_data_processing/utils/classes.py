import json
import spacy
from spacy.cli import download
from typing import List, Dict, Optional, Tuple, Set
from sentence_transformers import SentenceTransformer, util
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np

# Import constants directly to ensure access to Enum-derived lists
import sys
module_dir = './'
sys.path.append(module_dir)
from src.constants import VALID_TREATMENTS, VALID_TISSUES, VALID_MEDIUMS

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
            self.map = {} # Used for generic or Medium
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
    
    def add_mapping(self, og_sample, grounded_sample)->None:
        """
        Updates the internal maps based on the difference between raw and grounded samples.
        """
        # 1. Treatments (List)
        raw_treats = og_sample.get('treatment', [])
        ground_treats = grounded_sample.get('treatment', [])
        
        # We can only map if lists are same length, or logic gets complex.
        # Simple heuristic: If we resolved it, map raw -> resolved.
        # (For accurate tracking, the extractor ideally returns a dict, 
        # but here we assume the grounder maintained order or we map knowns)
        for t in raw_treats:
            # This is a simplification; in batch_process we fill the map explicitly
            pass 
            
    def save_map(self):
        if self.path:
            with open(self.path+'/map_treatment.json', 'w') as f:
                json.dump(self.map_treatment, f, indent=4)
            with open(self.path+'/map_tissue.json', 'w') as f:
                json.dump(self.map_tissue, f, indent=4)
            with open(self.path+'/map.json', 'w') as f:
                json.dump(self.map, f, indent=4)


class GroundingOptimizer:
    def __init__(self):
        print("Initializing GroundingOptimizer...")
        
        # 1. LOAD MODEL (BioLORD is best for biomedical ontologies)
        # It aligns 'Arabidopsis' closer to 'Plant' than 'Gun'
        print("Loading BioLORD model...")
        self.model = SentenceTransformer('FremyCompany/BioLORD-2023')
        
        # 2. LOAD LEMMATIZER (Spacy)
        # Crucial for mapping "shoots" -> "shoot", "leaves" -> "leaf"
        try:
            self.nlp = spacy.load("en_core_web_sm", disable=['parser', 'ner'])
        except OSError:
            print("Downloading spacy model...")
            download("en_core_web_sm")
            self.nlp = spacy.load("en_core_web_sm", disable=['parser', 'ner'])
        
        # 3. LOAD ONTOLOGIES
        self.valid_treatments = VALID_TREATMENTS
        self.valid_tissues = VALID_TISSUES
        self.valid_mediums = VALID_MEDIUMS

        # 4. PRE-COMPUTE VECTORS (With Lemmatization)
        print("Pre-computing ontology vectors...")
        self.treatment_ontology_vectors = self._encode_ontology(self.valid_treatments)
        self.tissue_ontology_vectors = self._encode_ontology(self.valid_tissues)
        self.medium_ontology_vectors = self._encode_ontology(self.valid_mediums)
        
    def _lemmatize(self, text: str) -> str:
        """
        Converts 'shoots' -> 'shoot', 'leaves' -> 'leaf', 'grown' -> 'grow'
        to ensure plural/tense differences don't break similarity.
        """
        if not text: return ""
        doc = self.nlp(str(text).lower())
        return " ".join([token.lemma_ for token in doc])

    def _encode_ontology(self, terms: list):
        """Helper to encode list of terms after lemmatizing"""
        lemmatized_terms = [self._lemmatize(t) for t in terms]
        return self.model.encode(lemmatized_terms)

    def find_semantic_match(self, raw_term: str, ontology: str, threshold: float = 0.85) -> Optional[str]:
        """
        Checks if the raw term is semantically identical to a known label
        without using an LLM. Returns the label if match > threshold.
        """
        if not raw_term:
            return None
            
        label, score = self.get_best_match_with_score(raw_term, ontology)
        
        if score > threshold:
            return label
        return None

    def get_best_match_with_score(self, text: str, category: str) -> Tuple[Optional[str], float]:
        """
        Calculates similarity of text against the specified ontology category.
        Returns: (best_matching_label, similarity_score)
        """
        if not text:
            return None, 0.0

        if category == 'treatment':
            vectors = self.treatment_ontology_vectors
            labels = self.valid_treatments
        elif category == 'tissue':
            vectors = self.tissue_ontology_vectors
            labels = self.valid_tissues
        elif category == 'medium':
            vectors = self.medium_ontology_vectors
            labels = self.valid_mediums
        else:
            return None, 0.0

        # CRITICAL: Lemmatize input before embedding to match ontology
        clean_text = self._lemmatize(text)
        
        # Encode input
        query_vec = self.model.encode([clean_text])[0]
        
        # Calculate Cosine Similarity
        # util.cos_sim is often faster/cleaner than sklearn for SentenceTransformers
        scores = util.cos_sim(query_vec, vectors)[0]
        
        # Find best
        best_idx = scores.argmax().item()
        best_score = scores[best_idx].item()
        best_label = labels[best_idx]
        
        return best_label, best_score

    def batch_process_study(self, data: Dict, extracted_samples: List[Dict], llm_func, llm_func_tis, label_map_obj: LabelMap):
        """
        Process a list of extracted samples, resolving terms via:
        1. Existing LabelMap (Memory)
        2. Vector Similarity (BioLORD)
        3. LLM (Batch Fallback)
        """
        
        # Local cache for this study to prevent redundant lookups
        mapping_cache_tissue = {}
        mapping_cache_treatment = {}
        mapping_cache_medium = {}
        
        # Sets to collect unknowns for LLM batching
        unknown_treatments_for_llm = set()
        unknown_tissues_for_llm = set()

        # --- HELPER: Resolve logic ---
        def resolve_term(term, specific_map, ontology_type):
            if not term: return None
            
            # 1. Check LabelMap
            if specific_map and term in specific_map:
                return specific_map[term]
            
            # 2. Check Vector Match (BioLORD)
            # We use a slightly stricter threshold for auto-mapping (0.85) to avoid bad data
            match = self.find_semantic_match(term, ontology_type, threshold=0.85)
            if match:
                return match
            
            return None

        # --- STEP 1: Process Fields locally ---
        
        for sample in extracted_samples:
            
            # A. Process TREATMENT (List)
            treatments = sample.get('treatment', [])
            # Ensure it's a list
            if isinstance(treatments, str): treatments = [treatments]
            
            for t in treatments:
                if t in mapping_cache_treatment: continue 
                
                resolved = resolve_term(t, label_map_obj.map_treatment, 'treatment')
                if resolved:
                    mapping_cache_treatment[t] = resolved
                else:
                    unknown_treatments_for_llm.add(t)

            # B. Process TISSUE (String)
            tissue = sample.get('tissue')
            if isinstance(tissue, list): 
                tissue = tissue[0] if tissue else "unknown"
                
            if tissue and tissue not in mapping_cache_tissue:
                resolved_tis = resolve_term(tissue, label_map_obj.map_tissue, 'tissue')
                if resolved_tis:
                    mapping_cache_tissue[tissue] = resolved_tis
                else:
                    unknown_tissues_for_llm.add(tissue)

            # C. Process MEDIUM (String)
            medium = sample.get('medium')
            if isinstance(medium, list):
                medium = medium[0] if medium else "unspecified"
                
            if medium and medium not in mapping_cache_medium:
                # Mediums often map to generic map or self-contained logic
                resolved_m = resolve_term(medium, label_map_obj.map, 'medium')
                if resolved_m:
                    mapping_cache_medium[medium] = resolved_m
                else:
                    # If vector fails for medium, we often default to the raw term 
                    # or map to "Unspecified" if it's too weird. 
                    # For now, let's keep raw if no match.
                    mapping_cache_medium[medium] = medium

        # --- STEP 2: Batch LLM Calls (Fallback) ---
        data_study = data.get('study_metadata', {})
        
        # A. Treatments LLM
        if unknown_treatments_for_llm:
            print(f"Querying LLM for {len(unknown_treatments_for_llm)} unique treatments...")
            try:
                llm_result = llm_func(list(unknown_treatments_for_llm), data_study)
                if llm_result:
                    mapping_cache_treatment.update(llm_result)
            except Exception as e:
                print(f"LLM Batch Error (Treatments): {e}")
        
        # B. Tissues LLM
        if unknown_tissues_for_llm:
            print(f"Querying LLM for {len(unknown_tissues_for_llm)} unique tissues...")
            try:
                llm_result = llm_func_tis(list(unknown_tissues_for_llm), data_study)
                if llm_result:
                    mapping_cache_tissue.update(llm_result)
            except Exception as e:
                print(f"LLM Batch Error (Tissues): {e}")

        # --- STEP 3: Apply mappings and Return ---
        final_samples = []
        for sample in extracted_samples:
            new_sample = sample.copy()
            
            # 1. Update Treatments
            new_treatments = []
            raw_treats = sample.get('treatment', [])
            if isinstance(raw_treats, str): raw_treats = [raw_treats]
            
            for t in raw_treats:
                # Check cache (populated by Vector or LLM)
                val = mapping_cache_treatment.get(t)
                
                # If still missing, check if it was in the unknown set (meaning LLM failed to return it)
                if not val and t in unknown_treatments_for_llm:
                    val = "Other stress" # Default fallback
                
                if val:
                    new_treatments.append(val)
                else:
                    # If it wasn't unknown (maybe skipped), keep raw? 
                    # Usually better to map to Other if uncertain.
                    new_treatments.append("Other stress")
                    
            new_sample['treatment'] = sorted(list(set(new_treatments)))
            
            # 2. Update Tissue
            raw_tissue = sample.get('tissue')
            if isinstance(raw_tissue, list): raw_tissue = raw_tissue[0] if raw_tissue else "unknown"
            
            if raw_tissue in mapping_cache_tissue:
                new_sample['tissue'] = mapping_cache_tissue[raw_tissue]
            else:
                new_sample['tissue'] = "unknown"
                
            # 3. Update Medium
            raw_medium = sample.get('medium')
            if isinstance(raw_medium, list): raw_medium = raw_medium[0] if raw_medium else "unspecified"
            
            if raw_medium in mapping_cache_medium:
                new_sample['medium'] = mapping_cache_medium[raw_medium]
            else:
                new_sample['medium'] = "unspecified"

            final_samples.append(new_sample)
            
        return final_samples