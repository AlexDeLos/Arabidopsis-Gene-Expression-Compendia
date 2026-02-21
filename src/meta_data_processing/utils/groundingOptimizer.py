import spacy
from typing import List, Dict, Optional, Tuple
from sentence_transformers import SentenceTransformer, util
from spacy.cli import download
import torch
from src.constants import UNIQUE_LABELS, CONTROL_MAP
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
        self.model = SentenceTransformer('FremyCompany/BioLORD-2023',device='cuda')
        
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
        return torch.tensor(self.model.encode(lemmatized_terms, device='cuda' if torch.cuda.is_available() else 'cpu', convert_to_tensor=True))

    # --- PUBLIC: Canonicalization Helper ---
    def canonicalize_term(self, term: str, category: str) -> str:
        """
        Public wrapper to collapse a synonym (e.g. 'Dark') to its canonical enum value 
        (e.g. 'Low Light Stress') using the generated CANONICAL_MAP.
        """
        if not term: return "unknown"
        if category not in CANONICAL_MAP: return term
        
        # 1. Direct Lookup
        if term in CANONICAL_MAP[category]:
            return CANONICAL_MAP[category][term]
            
        # 2. Case-insensitive Lookup (Fallback)
        term_lower = term.lower()
        for syn, canonical in CANONICAL_MAP[category].items():
            if syn.lower() == term_lower:
                return canonical
                
        # 3. If no map found, return original (or could return 'unknown')
        return term 

    def find_semantic_match(self, string: str, words: List[str], target_embeddings, threshold: float) -> Tuple[bool, float, Optional[str], Optional[str]]:
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
        Finds the best match within the EXPLICIT ontology (Synonyms included),
        but returns the CANONICAL (bucket) term.
        """
        if not text: 
            return None, 0.0

        # 1. Prepare Query
        clean_text = self._lemmatize(text)
        target_device = self.vectors['explicit'][category].device
        query_vec = self.model.encode(clean_text, convert_to_tensor=True, device=target_device)

        # 2. Compute Cosine Similarity against EXPLICIT list
        scores = util.cos_sim(query_vec, self.vectors['explicit'][category])[0]
        
        # 3. Find Best Match
        best_idx = torch.argmax(scores).item()
        best_score = scores[best_idx].item()
        
        # 4. Retrieve matched raw term (e.g. "Dark")
        matched_raw_term = EXPLICIT_KEYWORDS[category][best_idx]
        
        # 5. Collapse to Canonical (e.g. "Low Light Stress")
        canonical_term = self.canonicalize_term(matched_raw_term, category)
        
        return canonical_term, best_score
    def remove_redundant_unknowns(self,samples):
        cleaned_samples = []
        
        for sample in samples:
            cleaned_sample = {}
            for key, value in sample.items():
                # Check if the value is a list (to avoid touching strings like sample_id)
                if isinstance(value, list):
                    # Filter out the 'unknown' elements
                    valid_labels = [item for item in value if item != 'unknown']
                    
                    # If there are valid labels left, use the filtered list
                    if len(valid_labels) > 0:
                        cleaned_sample[key] = valid_labels
                    else:
                        # If the list was ONLY 'unknown', keep it as is
                        cleaned_sample[key] = value
                else:
                    # Keep non-list items (like 'sample_id') exactly the same
                    cleaned_sample[key] = value
                    
            cleaned_samples.append(cleaned_sample)
            
        return cleaned_samples
    
    def batch_process_study(self, data: Dict, extracted_samples: List[Dict], label_map: LabelMap) -> List:

        local_cache = {label: {} for label in LABELS}
        final_samples = []

        # --- PASS 1: Grounding & Scoring ---
        for sample in extracted_samples:
            
            for key, val in sample.items():
                if key == 'sample_id': continue
                if val and val != set():
                    if isinstance(val, set):
                        val_str = list(val)[0]
                    else:
                        val_str = str(val)
                    val_ = re.sub(r'[^a-zA-Z0-9 ,./-]', '', val_str)
                    sample[key] = set(val_.split(','))

            for label in LABELS:
                raw_terms = sample.get(label, 'unspecified')
                if isinstance(raw_terms, str): raw_terms = {raw_terms}

                for raw_term in raw_terms:
                    best_fit: str = 'unknown'
                    max_score = 0.55
                    
                    if raw_term not in local_cache[label]:
                        mapped_val = label_map._get_value(label, raw_term)
                        
                        if mapped_val is not None:
                            # Trust label map completely. Assign a massive score (2.0) so it always wins.
                            canonical_fit = self.canonicalize_term(mapped_val, label)
                            local_cache[label][raw_term] = {"val": canonical_fit, "score": 2.0}
                        else:
                            # 1. Bucket Match (Canonical)
                            _, score, _, best_keyword = self.find_semantic_match(raw_term, BUCKET_KEYWORDS[label], self.vectors['bucket'][label], threshold=0.88)
                            
                            # 2. Explicit Match (Synonyms)
                            _, score_, _, best_keyword_ = self.find_semantic_match(raw_term, EXPLICIT_KEYWORDS[label], self.vectors['explicit'][label], threshold=0.88)
                            
                            if score_ > score:
                                score = score_
                                best_keyword = best_keyword_
                            
                            if score >= max_score:
                                best_fit = str(best_keyword)
                                max_score = score
                        
                            if best_fit and best_fit != 'unknown':
                                canonical_fit = self.canonicalize_term(best_fit, label)
                                local_cache[label][raw_term] = {"val": canonical_fit, "score": max_score}
                            else:
                                local_cache[label][raw_term] = {"val": "unknown", "score": 0.0}

        # --- PASS 2: Apply Rules (Unique & Control) ---
        for sample in extracted_samples:
            new_sample = sample.copy()
            for label in LABELS:
                raw = sample.get(label, set())
                if isinstance(raw, str): raw = {raw}
                
                if not raw or raw == {'unspecified'}:
                    new_sample[label] = ["unspecified"]
                    continue
                
                # Fetch mapped items and their scores
                mapped_items = []
                for t in raw:
                    item = local_cache[label].get(t, {"val": "unknown", "score": 0.0})
                    mapped_items.append(item)

                # RULE 1: Highest-score Unique Label
                if label in UNIQUE_LABELS and len(mapped_items) > 1:
                    # Sort descending by score
                    mapped_items.sort(key=lambda x: x["score"], reverse=True)
                    # Pick the highest scoring item that isn't 'unknown' (or fallback to unknown if that's all there is)
                    best_item = next((item for item in mapped_items if item["val"] != "unknown"), mapped_items[0])
                    final_vals = [best_item["val"]]
                else:
                    final_vals = [item["val"] for item in mapped_items]

                # Clean up into a unique list
                final_vals = list(sorted(set(final_vals)))

                # RULE 2: Control Dominance
                control_val = CONTROL_MAP.get(label)
                if control_val and control_val in final_vals:
                    final_vals = [control_val]
                
                new_sample[label] = final_vals
                
            final_samples.append(new_sample)

        return self.remove_redundant_unknowns(final_samples)