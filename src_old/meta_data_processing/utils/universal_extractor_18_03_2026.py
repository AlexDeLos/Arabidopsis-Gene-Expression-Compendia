import re
from typing import Dict, List, Any, Optional, Set, Tuple
from sentence_transformers import util
from src.constants import TissueEnum, TreatmentEnum, MediumEnum, EXPLICIT_KEYWORDS, AREA_KEYWORDS, LABELS
from src.meta_data_processing.utils.groundingOptimizer import GroundingOptimizer
import sys
import torch
module_dir = './'
sys.path.append(module_dir)

# --- PRE-COMPILED REGEX PATTERNS ---
SPLIT_PATTERN = re.compile(r"[\[\]:;.,']+")
CLEAN_CANDIDATE_PATTERN = re.compile(r'\b(rep\d+|col-0|atgen)\b')
NOISE_FILTER_PATTERN = re.compile(r'[^a-zA-Z0-9\s]')
TOKEN_BOUNDARY_PATTERN = re.compile(r'\b\w+\b')

def complex_split(text: str) -> list[str]:
    if not text: return []
    return [token.strip() for token in SPLIT_PATTERN.split(text) if token.strip()]

def condense_candidates(candidates: Set[str], optimizer: GroundingOptimizer, category: str) -> Set[str]:
    """
    Selects the single best term from candidates and returns its CANONICAL form.
    """
    if not candidates:
        return set()

    cleaned_candidates = set()
    for cand in candidates:
        clean = CLEAN_CANDIDATE_PATTERN.sub('', cand.lower()).strip()
        if len(clean) > 2:
            cleaned_candidates.add(clean)

    if not cleaned_candidates:
        return set()

    best_term = None
    best_score = -1.0
    
    for term in cleaned_candidates:
        # get_best_match_with_score now returns the CANONICAL term
        label, score = optimizer.get_best_match_with_score(term, category)
        
        if score > best_score:
            best_score = score
            best_term = label # This is already canonical (e.g. "Low Light Stress")
        elif abs(score - best_score) < 0.02:
            # If scores tie, prefer the shorter canonical name (heuristic)
            if best_term and len(label) < len(best_term):
                best_term = label

    return set([best_term] if best_term else [])

def extract_valid_candidates(candidate_set: Set[str], optimizer: GroundingOptimizer, target_category: str, study_info:bool, top_k: int = 5) -> List[Tuple[str, float]]:
    if target_category == 'treatment': contrast_category = 'tissue'
    elif target_category == 'tissue': contrast_category = 'treatment' 
    elif target_category == 'medium': contrast_category = 'tissue' 
    else: contrast_category = None

    NOISE_STOPLIST = {
        "wild type", "wild-type", "columbia", "col-0", "mutant", "genotype", 
        "transgenic", "plants", "analysis", "data", "replicate", "study", 
        "experiment", "grown", "growth", "independent", "agb1", "gpa1",
        "control", "mock", "buffer", "treated", "samples", "using", "total", "rna",
        "rep1", "rep2", "rep3", "atgen", "unknown"
    }

    valid_raw_terms = {} 
    golden_keywords_list = EXPLICIT_KEYWORDS.get(target_category, [])

    for term in candidate_set:
        term_lower = term.lower()
        
        if any(bad in term_lower for bad in NOISE_STOPLIST): continue
        if len(term) < 2: continue
        if len(term_lower.split(' ')) > 10: continue

        # Check Golden Keywords (Synonyms included)
        has_golden_keyword, socre, _, _ = has_golden_key_word(term_lower, golden_keywords_list, optimizer, 0.75)
        
        if has_golden_keyword:
            valid_raw_terms[term] = socre 
            continue

        _, target_score = optimizer.get_best_match_with_score(term, category=target_category)
        
        if target_score < 0.75: continue

        if contrast_category:
            _, contrast_score = optimizer.get_best_match_with_score(term, category=contrast_category)
            if contrast_score > target_score: continue

        valid_raw_terms[term] = target_score

    sorted_candidates = sorted(valid_raw_terms.items(), key=lambda x: x[1], reverse=True)
    return sorted_candidates[:top_k]

def has_golden_key_word(string: str, words: list[str], optimized: GroundingOptimizer, th: float) -> Tuple[bool, float, Optional[str], Optional[str]]:
    if not string or not words: return False, 0.0, None, None

    string_lower = string.lower()
    for w in words:
        if w.lower() in string_lower:
            return True, 1.0, w.lower(), w 

    vectored = False
    for label in LABELS:
        if words is EXPLICIT_KEYWORDS[label]:
            target_embeddings = optimized.vectors['explicit'][label]
            vectored= True
    if not vectored:
        target_embeddings = optimized.model.encode(words, convert_to_tensor=True)

    tokens = TOKEN_BOUNDARY_PATTERN.findall(string_lower)
    if not tokens: return False, 0.0, None, None

    candidates = set(tokens)
    if len(tokens) > 1: candidates.update([" ".join(tokens[i:i+2]) for i in range(len(tokens)-1)])
    if len(tokens) > 2: candidates.update([" ".join(tokens[i:i+3]) for i in range(len(tokens)-2)])
    
    candidates_list = list(candidates)
    if not candidates_list: return False, 0.0, None, None

    candidate_embeddings = optimized.model.encode(candidates_list, convert_to_tensor=True)
    
    if target_embeddings.device != candidate_embeddings.device:
        target_embeddings = target_embeddings.to(candidate_embeddings.device)

    cosine_scores = util.cos_sim(candidate_embeddings, target_embeddings)
    max_val, max_idx_flat = torch.max(cosine_scores.flatten(), dim=0)
    max_score = max_val.item()
    
    if max_score > th:
        num_cols = cosine_scores.shape[1]
        row_idx = (max_idx_flat // num_cols).item()
        col_idx = (max_idx_flat % num_cols).item()
        best_candidate = candidates_list[row_idx]
        best_keyword = words[col_idx]
        return True, max_score, best_candidate, best_keyword

    return False, max_score, None, None


import re
import torch
from typing import Dict, List, Set, Any, Tuple
from sentence_transformers import util
from src.constants import LABELS, LABEL_CONFIG, EXPLICIT_KEYWORDS

class UniversalExtractor:
    def __init__(self, optimizer):
        self.optimizer = optimizer
        self.model = optimizer.model
        # Pre-compute trigger vectors for column identification
        self.trigger_vectors = {}
        self.search_col_vectors = {}
        self._column_cache = {}
        for label in LABELS:
            triggers = LABEL_CONFIG[label].get('search_triggers', [])
            # Combine triggers with explicit keywords for better coverage
            combined_triggers = list(set(triggers + EXPLICIT_KEYWORDS.get(label, [])))
            self.trigger_vectors[label] = self.model.encode(combined_triggers, convert_to_tensor=True)
            cols = LABEL_CONFIG[label].get('priority_cols', [])
            self.search_col_vectors[label] = self.model.encode(cols, convert_to_tensor=True)

    def extract(self, sample_metadata: Dict, study_metadata: Dict, study_id: str) -> Dict[str, Set[str]]:
        extracted = {}
        
        for label_type in LABELS:
            # Step 1: Check Explicit Columns (Column names matching the label_type)
            hits = self._extract_from_matched_columns(sample_metadata, label_type)
            
            # Step 2: Fallback to Priority Text (Sample Title/Characteristics)
            if not hits:
                priority_text = self._get_text_blob(sample_metadata, ['title', 'characteristics_ch1', 'source_name_ch1'])
                hits = self._semantic_ngram_search(priority_text, label_type)
            
            # Step 3: Fallback to Study Level Metadata
            if not hits:
                study_text = self._get_text_blob(study_metadata, ['summary', 'overall_design'])
                hits = self._semantic_ngram_search(study_text, label_type)
            
            # Requirement: If empty, set to unspecified
            extracted[label_type] = hits if hits else {"unspecified"}
            
        return extracted

    def _get_text_blob(self, metadata: Dict, keys: List[str]) -> str:
        """Helper to collect text from specific metadata fields."""
        bits = []
        for k in keys:
            val = metadata.get(k, "")
            if isinstance(val, list): bits.extend([str(v) for v in val])
            else: bits.append(str(val))
        return " ".join(bits).lower()

    def _is_relevant_column(self, key: str, label_type: str, threshold: float = 0.82) -> bool:
        cache_key = f"rel_{key}_{label_type}"
        if cache_key in self._column_cache:
            return self._column_cache[cache_key]

        if label_type not in self.trigger_vectors: 
            return False

        clean_key = key.replace('_', ' ').lower()
        key_vec = self.model.encode(clean_key, convert_to_tensor=True)
        scores = util.cos_sim(key_vec, self.trigger_vectors[label_type])
        result = torch.max(scores).item() > threshold
        
        self._column_cache[cache_key] = result
        return result

    def _is_priority_column(self, key: str, label_type: str, threshold: float = 0.82) -> bool:
        cache_key = f"prio_{key}_{label_type}"
        if cache_key in self._column_cache:
            return self._column_cache[cache_key]

        if label_type not in self.search_col_vectors:
            return False

        clean_key = key.replace('_', ' ').lower()
        key_vec = self.model.encode(clean_key, convert_to_tensor=True)
        scores = util.cos_sim(key_vec, self.search_col_vectors[label_type])
        result = torch.max(scores).item() > threshold
        
        self._column_cache[cache_key] = result
        return result
    
    def _is_relevant_element_old(self, text: str, label_type: str, threshold: float = 0.82) -> bool:
        """
        Modified version of _is_relevant_column.
        Evaluates if a specific string (like a sub-key from a list element) is relevant to the label type.
        """
        cache_key = f"rel_elem_{text}_{label_type}"
        if cache_key in self._column_cache:
            return self._column_cache[cache_key]

        if label_type not in self.trigger_vectors: 
            return False

        clean_text = text.replace('_', ' ').lower().strip()
        text_vec = self.model.encode(clean_text, convert_to_tensor=True)
        
        scores = util.cos_sim(text_vec, self.trigger_vectors[label_type])
        result = torch.max(scores).item() > threshold
        
        self._column_cache[cache_key] = result
        return result

    def _is_relevant_element(self, text: str, label_type: str, threshold: float = 0.65) -> bool:
        clean_text = text.replace('_', ' ').lower().strip()
        
        # 1. FAST PATH: Instant Keyword Match (Bypasses vector dilution)
        # If the text explicitly says "stress" or "treatment", we know it's relevant!
        triggers = LABEL_CONFIG[label_type].get('search_triggers', [])
        for trigger in triggers:
            # Matches whole words to prevent accidental substring matches
            if re.search(fr'\b{re.escape(trigger)}\b', clean_text):
                return True

        # 2. VECTOR FALLBACK: Topic Truncation for implied meanings
        words = clean_text.split()
        if len(words) > 4:
            clean_text = " ".join(words[:4]) # Truncate tighter to 4 words

        cache_key = f"rel_elem_{clean_text}_{label_type}"
        if cache_key in self._column_cache:
            return self._column_cache[cache_key]

        if label_type not in self.trigger_vectors: 
            return False

        text_vec = self.model.encode(clean_text, convert_to_tensor=True)
        scores = util.cos_sim(text_vec, self.trigger_vectors[label_type])
        
        # Lowered threshold allows for slight dilution
        result = torch.max(scores).item() > threshold
        
        self._column_cache[cache_key] = result
        return result
    def _extract_from_matched_columns(self, metadata: Dict, label_type: str) -> Set[str]:
        """Scans keys and priority lists for semantic matches to the label type and extracts values."""
        found = set()
        for key, value in metadata.items():
            vals = value if isinstance(value, list) else [value]
            
            # --- CASE 1: The key is a priority column (e.g., 'characteristics_ch1') ---
            if self._is_priority_column(key, label_type):
                for val in vals:
                    val_str = str(val).strip()
                    # Removed replace('-',' ') here to protect terms like 'Col-0' or 'UV-B'
                    val_str_clean = val_str.replace('_', ' ')
                    
                    if ":" in val_str_clean:
                        sub_key, actual_val = val_str_clean.split(":", 1) 
                        
                        if self._is_relevant_element(sub_key.strip(), label_type):
                            actual_val = actual_val.strip()
                            
                            # THE FIX: If the value is long, extract ONLY the keywords
                            if len(actual_val.split()) > 3:
                                ngrams = self._semantic_ngram_search(actual_val, label_type, threshold=0.75)
                                found.update(ngrams)
                            else:
                                found.add(actual_val)
                    else:
                        # NO COLON exists (e.g., pure paragraph in a protocol column)
                        # Check if the start of the string indicates it's relevant to this label_type
                        if self._is_relevant_element(val_str_clean, label_type):
                            if len(val_str_clean.split()) > 3:
                                # This will extract ONLY "Osmotic stress" from your massive paragraph!
                                ngrams = self._semantic_ngram_search(val_str_clean, label_type, threshold=0.75)
                                found.update(ngrams)
                            else:
                                found.add(val_str_clean)
                                
            # --- CASE 2: The key itself is highly relevant (e.g., key='tissue') ---
            elif self._is_relevant_column(key, label_type):
                for val in vals:
                    val_str = str(val).strip()
                    
                    if ":" in val_str:
                        val_str = val_str.split(":", 1)[-1].strip()
                        
                    if len(val_str.split()) > 3:
                        ngrams = self._semantic_ngram_search(val_str, label_type, threshold=0.75)
                        found.update(ngrams)
                    else:
                        found.add(val_str)
                        
        return {f for f in found if len(f) > 1}

    def _extract_from_matched_columns_old(self, metadata: Dict, label_type: str) -> Set[str]:
        """Scans keys and priority lists for semantic matches to the label type and extracts values."""
        found = set()
        for key, value in metadata.items():
            vals = value if isinstance(value, list) else [value]
            
            # --- CASE 1: The key is a priority column (e.g., 'characteristics_ch1') ---
            if self._is_priority_column(key, label_type):
                for val in vals:
                    val_str = str(val).strip()
                    val_str = val_str.replace('_',' ').replace('-',' ')
                    
                    # Elements often come as 'sub_key: value' (e.g., 'tissue: whole plant')
                    if ":" in val_str:
                        # Split only on the first colon
                        sub_key, actual_val = val_str.split(":", 1) 
                        
                        # Check if the sub_key (e.g., 'developmental stage') is relevant
                        if self._is_relevant_element(sub_key.strip(), label_type):
                            actual_val = actual_val.strip()
                            
                            if len(actual_val.split()) > 3:
                                ngrams = self._semantic_ngram_search(actual_val, label_type, threshold=0.75)
                                found.update(ngrams)
                            else:
                                found.add(actual_val)
                    else:
                        # If no colon exists, check if the whole string contains relevant trigger words
                        if self._is_relevant_element(val_str, label_type):
                            if len(val_str.split()) > 3:
                                ngrams = self._semantic_ngram_search(val_str, label_type, threshold=0.75)
                                found.update(ngrams)
                            else:
                                found.add(val_str)
                                
            # --- CASE 2: The key itself is highly relevant (e.g., key='tissue') ---
            elif self._is_relevant_column(key, label_type):
                for val in vals:
                    val_str = str(val).strip()
                    
                    if ":" in val_str:
                        val_str = val_str.split(":", 1)[-1].strip()
                        
                    if len(val_str.split()) > 3:
                        ngrams = self._semantic_ngram_search(val_str, label_type, threshold=0.75)
                        found.update(ngrams)
                    else:
                        found.add(val_str)
                        
        return {f for f in found if len(f) > 1}
    def _semantic_ngram_search(self, text: str, label_type: str, threshold: float = 0.75) -> Set[str]:
        if not text: return set()
        
        # 1. Strip out noise words that poison the vectors
        noise_words = r'\b(genotype|mutant|wild-type|wildtype|age|weeks|days|old|developmental stage|cell type)\b'
        clean_text = re.sub(noise_words, '', text, flags=re.IGNORECASE)
        
        # 2. THE FIX FOR TITLES: Break apart squished string formatting
        # Converts "AtGen_6-2521_Osmoticstress-Roots" -> "AtGen 6 2521 Osmoticstress Roots"
        clean_text = clean_text.replace('_', ' ').replace('-', ' ')
        
        # Tokenize (Now keeping forward slash for things like "1/2 MS")
        words = re.findall(r'\b[\w/]+\b', clean_text)
        if not words: return set()

        # Generate N-grams (1 to 3 words)
        candidates = []
        for n in range(1, 4):
            for i in range(len(words) - n + 1):
                candidates.append(" ".join(words[i:i+n]))
        
        if not candidates: return set()

        unique_candidates = list(set(candidates))
        candidate_vectors = self.model.encode(unique_candidates, convert_to_tensor=True)
        target_vectors = self.optimizer.vectors['explicit'][label_type]
        
        cosine_matrix = util.cos_sim(candidate_vectors, target_vectors)
        matches = torch.where(cosine_matrix > threshold)
        
        results = set()
        for idx in matches[0]:
            results.add(unique_candidates[idx.item()])
            
        return results
