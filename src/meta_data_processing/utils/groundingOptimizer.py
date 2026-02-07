import spacy
from typing import List, Dict, Optional, Tuple
from sentence_transformers import SentenceTransformer, util
from spacy.cli import download
import torch
import re

# Import constants directly
import sys
module_dir = './'
sys.path.append(module_dir)
from src.constants import VALID_TREATMENTS, VALID_TISSUES, VALID_MEDIUMS, VALID_TREATMENTS_ALT,EXPLICIT_KEYWORDS
from src.meta_data_processing.utils.labelMap import LabelMap

class GroundingOptimizer:
    def __init__(self):
        print("Initializing GroundingOptimizer...")
        self.model = SentenceTransformer('FremyCompany/BioLORD-2023')
        
        try:
            self.nlp = spacy.load("en_core_web_sm", disable=['parser', 'ner'])
        except OSError:
            download("en_core_web_sm")
            self.nlp = spacy.load("en_core_web_sm", disable=['parser', 'ner'])
        
        self.valid_treatments = EXPLICIT_KEYWORDS['treatment']
        self.valid_tissues = EXPLICIT_KEYWORDS['tissue']
        self.valid_mediums = EXPLICIT_KEYWORDS['medium']

        print("Pre-computing ontology vectors...")
        self.treatment_vecs = self._encode_ontology(self.valid_treatments)
        self.tissue_vecs = self._encode_ontology(self.valid_tissues)
        self.medium_vecs = self._encode_ontology(self.valid_mediums)
        
    def _lemmatize(self, text: str) -> str:
        if not text: return ""
        doc = self.nlp(str(text).lower())
        return " ".join([token.lemma_ for token in doc])

    def _encode_ontology(self, terms: list):
        lemmatized_terms = [self._lemmatize(t) for t in terms]
        return self.model.encode(lemmatized_terms,device='cuda')

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
            
            # Special manual overrides
            if 'mock' in text_lower or 'control' in text_lower:
                return "No stress"

        elif category == 'medium':
            for medium in sorted(self.valid_mediums, key=len, reverse=True):
                if medium.lower() in text_lower:
                    return medium
                
        return None

    def find_semantic_match(self, string: str, words: List[str], threshold: float) -> Tuple[bool, float, Optional[str], Optional[str]]:
        """
        Checks if any n-gram in the input string is semantically similar to any word in the golden list.
        
        Args:
            string (str): The raw metadata string to check.
            words (List[str]): The list of golden keywords.
            threshold (float): Cosine similarity threshold (0.0 to 1.0).

        Returns:
            Tuple[bool, float, str, str]: (is_match, score, found_candidate, matched_keyword)
        """
        if not string or not words:
            return False, 0.0, None, None

        # 1. OPTIMIZATION: Exact string match (O(1))
        # Checks if the word exists literally in the string before doing math.
        string_lower = string.lower()
        for w in words:
            w_lower = w.lower()
            if w_lower in string_lower:
                return True, 1.0, w_lower, w 

        # 2. OPTIMIZATION: Use Cached Vectors
        # We check if the 'words' list passed in is EXACTLY one of our pre-computed lists
        # using the 'is' operator (pointer comparison).
        target_embeddings = None
        
        if words is self.valid_treatments:
            target_embeddings = self.treatment_vecs
        elif words is self.valid_tissues:
            target_embeddings = self.tissue_vecs
        elif words is self.valid_mediums:
            target_embeddings = self.medium_vecs
        else:
            # Fallback: Encode the list on the fly if it's a custom list
            # We use the model's current device (CPU or CUDA) automatically
            target_embeddings = self.model.encode(words, convert_to_tensor=True)

        # 3. Generate Candidates (N-grams)
        # We look for 1-grams, 2-grams, and 3-grams
        tokens = re.findall(r'\b\w+\b', string_lower)
        if not tokens:
            return False, 0.0, None, None

        candidates = set(tokens)
        
        # Add 2-grams (e.g., "osmotic stress")
        if len(tokens) > 1:
            candidates.update([" ".join(tokens[i:i+2]) for i in range(len(tokens)-1)])
        
        # Add 3-grams (e.g., "root apical meristem")
        if len(tokens) > 2:
            candidates.update([" ".join(tokens[i:i+3]) for i in range(len(tokens)-2)])
        
        candidates_list = list(candidates)
        
        # 4. Encode Candidates
        candidate_embeddings = self.model.encode(candidates_list, convert_to_tensor=True)

        # 5. Compute Cosine Similarity Matrix
        # Shape: [num_candidates x num_golden_words]
        # We ensure both tensors are on the same device for the calculation
        if target_embeddings.device != candidate_embeddings.device:
            target_embeddings = target_embeddings.to(candidate_embeddings.device)

        cosine_scores = util.cos_sim(candidate_embeddings, target_embeddings)

        # 6. Find Best Match Indices
        # Get the single highest value in the entire matrix
        max_score = torch.max(cosine_scores).item()
        
        if max_score > threshold:
            # Find the coordinates (row, col) of the max value
            max_idx = torch.argmax(cosine_scores)
            
            # Deconstruct the flat index into row (candidate) and col (golden word)
            num_cols = cosine_scores.shape[1]
            row_idx = (max_idx // num_cols).item()
            col_idx = (max_idx % num_cols).item()
            
            best_candidate = candidates_list[row_idx]
            best_keyword = words[col_idx]
            
            return True, max_score, best_candidate, best_keyword

        return False, max_score, None, None

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
                    # label_map.add_treatment(t, match)
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
                        # label_map.add_tissue(raw_tissue, match)
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
                        # label_map.add(raw_medium, match)
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
                    # label_map.add_mapping_dict(results, 'treatment')
            except Exception as e:
                print(f"    LLM Error (Treatment): {e}")

        if unknowns['tissue']:
            print(f"  > Sending {len(unknowns['tissue'])} tissues to LLM...")
            try:
                results = None#llm_func_tis(list(unknowns['tissue']), study_context)
                if results:
                    local_cache['tissue'].update(results)
                    # label_map.add_mapping_dict(results, 'tissue')
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
            new_sample['treatment'] = sorted(final_treats)
            
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