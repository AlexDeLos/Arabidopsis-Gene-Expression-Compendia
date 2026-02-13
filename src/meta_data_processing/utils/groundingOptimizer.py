import spacy
from typing import List, Dict, Optional, Tuple
from sentence_transformers import SentenceTransformer, util
from spacy.cli import download
import torch
import copy
import re
import sys

module_dir = './'
sys.path.append(module_dir)

# Import the new map alongside the others
from src.constants import LABELS, EXPLICIT_KEYWORDS, BUCKET_KEYWORDS, CANONICAL_MAP
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
        self.vectors = {'bucket':{}, 'explicit':{}}
        for label in LABELS:
            # Vectors are built on the config generated in constants.py
            self.vectors['bucket'][label] = self._encode_ontology(BUCKET_KEYWORDS[label])
            self.vectors['explicit'][label] = self._encode_ontology(EXPLICIT_KEYWORDS[label])

    def _lemmatize(self, text: str) -> str:
        if not text: return ""
        doc = self.nlp(str(text).lower())
        return " ".join([token.lemma_ for token in doc])

    def _encode_ontology(self, terms: list):
        lemmatized_terms = [self._lemmatize(t) for t in terms]
        return torch.tensor(self.model.encode(lemmatized_terms, device='cuda' if torch.cuda.is_available() else 'cpu'))

    # --- UPDATED: Heuristics using the Synonyms ---
    def _check_heuristics(self, text: str, category: str) -> Optional[str]:
        text_lower = text.lower().strip()
        
        IGNORED_TERMS = {
            'conditions', 'condition', 'time', 'date', 'replicate', 'rep', 
            'genotype', 'analysis', 'samples', 'study', 'experiment'
        }
        if text_lower in IGNORED_TERMS:
            if category == 'treatment': return "No stress"
            if category == 'tissue': return "unknown" 
            return None

        # Optimization: Check directly against the generated EXPLICIT list
        # We sort by length to match "Root system" before "Root"
        sorted_candidates = sorted(EXPLICIT_KEYWORDS[category], key=len, reverse=True)
        
        for candidate in sorted_candidates:
            # Word boundary check is safer: matches "leaf" but not "flea"
            # Using simple inclusion for now as per your original logic
            if candidate.lower() in text_lower:
                return candidate
        
        return None

    # --- NEW: Generic Collapsing Method ---
    def _get_canonical_term(self, category: str, term: str) -> str:
        """
        Collapses a synonym (e.g. 'Dark') to its canonical enum value 
        (e.g. 'Low Light Stress') using the generated CANONICAL_MAP.
        """
        if category not in CANONICAL_MAP:
            return term
        
        # 1. Direct Lookup
        if term in CANONICAL_MAP[category]:
            return CANONICAL_MAP[category][term]
            
        # 2. Case-insensitive Lookup (Fallback)
        term_lower = term.lower()
        for syn, canonical in CANONICAL_MAP[category].items():
            if syn.lower() == term_lower:
                return canonical
                
        return term # Return original if no mapping found

    def find_semantic_match(self, string: str, words: List[str], target_embeddings, threshold: float) -> Tuple[bool, float, Optional[str], Optional[str]]:
        # ... (This method remains exactly the same as your provided code) ...
        # [Use previous implementation]
        if not string or not words: return False, 0.0, None, None
        string_lower = string.lower()
        for w in words:
            if w.lower() in string_lower: return True, 1.0, w.lower(), w 

        vectored = False
        for label in LABELS:
            if words is BUCKET_KEYWORDS.get(label) or words is EXPLICIT_KEYWORDS.get(label):
                vectored= True
        
        if not vectored:
            target_embeddings = self.model.encode(words, convert_to_tensor=True)

        tokens = re.findall(r'\b\w+\b', string_lower)
        if not tokens: return False, 0.0, None, None

        candidates = set(tokens)
        if len(tokens) > 1: candidates.update([" ".join(tokens[i:i+2]) for i in range(len(tokens)-1)])
        if len(tokens) > 2: candidates.update([" ".join(tokens[i:i+3]) for i in range(len(tokens)-2)])
        
        candidates_list = list(candidates)
        candidate_embeddings = self.model.encode(candidates_list, convert_to_tensor=True)

        if target_embeddings.device != candidate_embeddings.device:
            target_embeddings = target_embeddings.to(candidate_embeddings.device)

        cosine_scores = util.cos_sim(candidate_embeddings, target_embeddings)
        max_idx = torch.argmax(cosine_scores)
        num_cols = cosine_scores.shape[1]
        row_idx = (max_idx // num_cols).item()
        col_idx = (max_idx % num_cols).item()
        
        best_candidate = candidates_list[row_idx]
        best_keyword = words[col_idx]
        max_score = torch.max(cosine_scores).item()
        
        return max_score > threshold, max_score, best_candidate, best_keyword
    
    
    def get_best_match_with_score(self, text: str, category: str) -> Tuple[Optional[str], float]:
        """
        Finds the best match for a string within the explicit ontology (including synonyms),
        but returns the CANONICAL (bucket) term.
        """
        if not text: 
            return None, 0.0

        # 1. Prepare Query
        # We ensure the query vector is created on the same device as the stored vectors
        clean_text = self._lemmatize(text)
        target_device = self.vectors['explicit'][category].device
        query_vec = self.model.encode(clean_text, convert_to_tensor=True, device=target_device)

        # 2. Compute Cosine Similarity against the EXPLICIT list (which includes synonyms)
        scores = util.cos_sim(query_vec, self.vectors['explicit'][category])[0]
        
        # 3. Find Best Match
        best_idx = torch.argmax(scores).item()
        best_score = scores[best_idx].item()
        
        # 4. Retrieve the matched term (This matches the index in EXPLICIT_KEYWORDS)
        # Note: This might be a synonym (e.g., "Dark")
        matched_raw_term = EXPLICIT_KEYWORDS[category][best_idx]
        
        # 5. Collapse to Canonical (e.g., "Dark" -> "Low Light Stress")
        # This ensures the extractor returns the standardized Enum value
        canonical_term = self._get_canonical_term(category, matched_raw_term)
        
        return canonical_term, best_score
    def batch_process_study(self, data: Dict, extracted_samples: List[Dict], label_map: LabelMap) -> List:
        
        local_cache = {label: {} for label in LABELS}
        final_samples = []

        # --- PASS 1: Grounding ---
        for sample in extracted_samples:
            
            # Clean sets (kept from original)
            for key, val in sample.items():
                if key == 'sample_id': continue
                if val != set():
                    val_ = re.sub(r'[^a-zA-Z0-9 ,]', '', val.pop())
                    sample[key] = set(val_.split(','))

            for label in LABELS:
                raw_tissues = sample.get(label, 'unspecified')
                if isinstance(raw_tissues, str): raw_tissues = {raw_tissues}

                for raw_tissue in raw_tissues:
                    best_fit: str = 'unknown'
                    max_score = 0.55
                    
                    if raw_tissue not in local_cache[label]:
                        if raw_tissue in label_map.map_tissue:
                            local_cache[label][raw_tissue] = label_map.map_tissue[raw_tissue]
                        else:
                            # 1. Bucket Match (Canonical)
                            _, score, _, best_keyword = self.find_semantic_match(raw_tissue, BUCKET_KEYWORDS[label], self.vectors['bucket'][label], threshold=0.88)
                            
                            # 2. Explicit Match (Synonyms)
                            _, score_, _, best_keyword_ = self.find_semantic_match(raw_tissue, EXPLICIT_KEYWORDS[label], self.vectors['explicit'][label], threshold=0.88)
                            
                            if score_ > score:
                                score = score_
                                best_keyword = best_keyword_
                            
                            if score >= max_score:
                                best_fit = str(best_keyword)
                                max_score = score
                        
                        if best_fit and best_fit != 'unknown':
                            # --- COLLAPSE SYNONYMS HERE ---
                            # This now uses the map generated in constants.py
                            canonical_fit = self._get_canonical_term(label, best_fit)
                            local_cache[label][raw_tissue] = canonical_fit
                        else:
                            local_cache[label][raw_tissue] = "unknown"

        # --- PASS 2: Apply ---
        for sample in extracted_samples:
            new_sample = sample.copy()
            for label in LABELS:
                raw = sample.get(label, set())
                if isinstance(raw, str): raw = {raw}
                
                if not raw or raw == {'unspecified'}:
                    new_sample[label] = ["unspecified"]
                else:
                    final_treats = []
                    for t in raw:
                        val = local_cache[label].get(t, "unknown")
                        final_treats.append(val)
                    # Sort and Deduplicate
                    new_sample[label] = list(sorted(set(final_treats)))
            final_samples.append(new_sample)

        return final_samples