import spacy
from typing import List, Dict, Optional, Tuple
from sentence_transformers import SentenceTransformer, util
from spacy.cli import download
import torch
import copy
import re

# Import constants directly
import sys
module_dir = './'
sys.path.append(module_dir)
from src.constants import LABELS,EXPLICIT_KEYWORDS,BUCKET_KEYWORDS
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
        
        print("Pre-computing ontology vectors...")
        self.vectors = {'bucket':{},'explicit':{}}
        for label in LABELS:
            self.vectors['bucket'][label] = self._encode_ontology(EXPLICIT_KEYWORDS[label])

            self.vectors['explicit'][label] = self._encode_ontology(EXPLICIT_KEYWORDS[label])
    def _lemmatize(self, text: str) -> str:
        if not text: return ""
        doc = self.nlp(str(text).lower())
        return " ".join([token.lemma_ for token in doc])

    def _encode_ontology(self, terms: list):
        lemmatized_terms = [self._lemmatize(t) for t in terms]
        return torch.tensor(self.model.encode(lemmatized_terms,device='cuda'))

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
            for tissue in sorted(EXPLICIT_KEYWORDS['tissue'], key=len, reverse=True):
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
            for medium in sorted(EXPLICIT_KEYWORDS['medium'], key=len, reverse=True):
                if medium.lower() in text_lower:
                    return medium
                
        return None

    def find_semantic_match(self, string: str, words: List[str],target_embeddings, threshold: float) -> Tuple[bool, float, Optional[str], Optional[str]]:
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
        vectored = False
        for label in LABELS:
            if words is BUCKET_KEYWORDS[label]:
                target_embeddings = self.vectors['bucket'][label]
                vectored= True
        if not vectored:
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
        
        # Find the coordinates (row, col) of the max value
        max_idx = torch.argmax(cosine_scores)
        
        # Deconstruct the flat index into row (candidate) and col (golden word)
        num_cols = cosine_scores.shape[1]
        row_idx = (max_idx // num_cols).item()
        col_idx = (max_idx % num_cols).item()
        
        best_candidate = candidates_list[row_idx]
        best_keyword = words[col_idx]
        
        return max_score> threshold, max_score, best_candidate, best_keyword


    def get_best_match_with_score(self, text: str, category: str) -> Tuple[Optional[str], float]:
        # (Same implementation as previous step)
        if not text: return None, 0.0

        clean_text = self._lemmatize(text)
        query_vec = self.model.encode([clean_text])[0]

        scores = util.cos_sim(query_vec, self.vectors['explicit'][category])[0]
        best_idx = scores.argmax().item()
        best_score = scores[best_idx].item()
        best_label = EXPLICIT_KEYWORDS[category][best_idx]
        return best_label, best_score
    def _get_canonical_term(self, category: str, term: str) -> str:
        """
        Maps a specific term (e.g., 'Dehydration', 'Dark') to its canonical bucket 
        term (e.g., 'Dehydration Stress', 'Low Light Stress').
        """
        # 1. Trivial Check: If it's not treatment, we currently assume 1-to-1 mapping
        if category != 'treatment':
            return term

        # 2. Import Enums locally to avoid top-level dependency issues
        from src.constants import TreatmentEnum, TreatmentEnum_alt

        # 3. Check if term is already Canonical (optimization)
        if term in [t.value for t in TreatmentEnum]:
            return term

        # 4. Map Synonyms (Alt) -> Canonical
        # We match them by their Enum Variable Name (e.g., both have .DROUGHT)
        for alt in TreatmentEnum_alt:
            if alt.value == term:
                # Found the synonym (e.g., "Drought"), get the canonical name
                if hasattr(TreatmentEnum, alt.name):
                    return getattr(TreatmentEnum, alt.name).value

        # 5. Handle manual edge cases (from the ['light', 'dark'] list in constants)
        term_lower = term.lower()
        if term_lower == 'dark': return TreatmentEnum.LOW_LIGHT.value
        if term_lower == 'light': return TreatmentEnum.OTHER_LIGHT.value

        # Default: If no mapping found, return original (or "Other")
        return term

    def batch_process_study(self, data: Dict, extracted_samples: List[Dict], label_map: LabelMap)->List:
        
        local_cache = { 'treatment': {}, 'tissue': {}, 'medium': {} }
        unknowns = { 'treatment': set(), 'tissue': set(), 'medium': set() }
        final_samples = []

        # --- PASS 1: Try Maps, Heuristics, and Vectors ---
        for sample in extracted_samples:
            
            for key,val in sample.items():
                if key == 'sample_id':
                    continue
                if val !=set():
                    val_ = re.sub(r'[^a-zA-Z0-9 ,]','',val.pop())
                    sample[key]= set(val_.split(','))
            for label in LABELS:
                raw_tissues = sample.get(label, 'unspecified')
                for raw_tissue in raw_tissues:
                    best_fit:str = 'unknown'
                    max_score = 0.55
                    if raw_tissue not in local_cache[label]:
                    
                        if raw_tissue in label_map.map_tissue:
                            local_cache[label][raw_tissue] = label_map.map_tissue[raw_tissue]
                        else:
                            _, score, _, best_keyword = self.find_semantic_match(raw_tissue, BUCKET_KEYWORDS[label],self.vectors['bucket'][label], threshold=0.88)
                            _, score_, _, best_keyword_ = self.find_semantic_match(raw_tissue, EXPLICIT_KEYWORDS[label], self.vectors['explicit'][label],threshold=0.88)
                            if score_ > score:
                                score = score_
                                best_keyword = best_keyword_
                            if score>=max_score:
                                best_fit = str(best_keyword)
                                max_score = score
                        if best_fit:
                            canonical_fit = self._get_canonical_term(label, best_fit)
                            local_cache[label][raw_tissue] = canonical_fit
                        else:
                            raise ValueError('this was not planned to be reached')
                            # unknowns[label].add(raw_tissues)

        # --- PASS 3: Apply ---
        final_samples = []
        for sample in extracted_samples:
            new_sample = sample.copy()
            for label in LABELS:
                # Treatment
                raw = sample.get(label, set())
                if raw == set():
                    new_sample[label]=["unspecified"]
                else:
                    final_treats = []
                    for t in raw:
                        val = local_cache[label].get(t, "unknown")
                        final_treats.append(val)
                    new_sample[label] = list(sorted(final_treats))
            final_samples.append(new_sample)

        return final_samples