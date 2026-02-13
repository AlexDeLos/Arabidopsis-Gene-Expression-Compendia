import re
from typing import Dict, List, Any, Optional, Set, Tuple
from sentence_transformers import util
from src.constants import TissueEnum, TreatmentEnum, MediumEnum, EXPLICIT_KEYWORDS, AREA_KEYWORDS
from src.meta_data_processing.utils.groundingOptimizer import GroundingOptimizer
import sys
import torch
import numpy as np
from functools import lru_cache

module_dir = './'
sys.path.append(module_dir)
from src.constants import LABELS
# --- PRE-COMPILED REGEX PATTERNS (Performance Optimization) ---
# Compiling these once at module level is much faster than re-compiling per function call
SPLIT_PATTERN = re.compile(r"[\[\]:;.,']+")
CLEAN_CANDIDATE_PATTERN = re.compile(r'\b(rep\d+|col-0|atgen)\b')
NOISE_FILTER_PATTERN = re.compile(r'[^a-zA-Z0-9\s]')
TOKEN_BOUNDARY_PATTERN = re.compile(r'\b\w+\b')

def complex_split(text: str) -> list[str]:
    """
    Splits text by [, ], :, ,, and ' characters.
    Removes empty strings and trims whitespace.
    """
    if not text:
        return []
    
    # Use pre-compiled pattern
    return [token.strip() for token in SPLIT_PATTERN.split(text) if token.strip()]

def condense_candidates(candidates: Set[str], optimizer: GroundingOptimizer, category: str) -> Set[str]:
    """
    Selects the single best term from a list of overlapping candidates 
    based on Semantic Purity.
    """
    if not candidates:
        return set()

    # 1. Cleaning Step
    cleaned_candidates = set()
    for cand in candidates:
        # Use pre-compiled regex for cleaning
        clean = CLEAN_CANDIDATE_PATTERN.sub('', cand.lower()).strip()
        if len(clean) > 2:
            cleaned_candidates.add(clean)

    if not cleaned_candidates:
        return set()

    # 2. Semantic Scoring Competition
    best_term = None
    best_score = -1.0
    
    # Note: If optimizer supported batch processing here, it would be faster.
    # Keeping as-is to preserve exact functionality of 'get_best_match_with_score'.
    for term in cleaned_candidates:
        label, score = optimizer.get_best_match_with_score(term, category)
        
        if score > best_score:
            best_score = score
            best_term = term
        elif abs(score - best_score) < 0.02:
            if len(term) < len(best_term):
                best_term = term

    return set([best_term] if best_term else [])

def extract_valid_candidates(candidate_set: Set[str], optimizer: GroundingOptimizer, target_category: str,study_info:bool, top_k: int = 5) -> List[Tuple[str, float]]:
    """
    Filters candidates using Keyword VIP Pass + Vector Similarity.
    """
    
    # 1. Define Distractors
    if target_category == 'treatment':
        contrast_category = 'tissue'
    elif target_category == 'tissue':
        contrast_category = 'treatment' 
    elif target_category == 'medium':
        contrast_category = 'tissue' 
    else:
        contrast_category = None

    NOISE_STOPLIST = {
        "wild type", "wild-type", "columbia", "col-0", "mutant", "genotype", 
        "transgenic", "plants", "analysis", "data", "replicate", "study", 
        "experiment", "grown", "growth", "independent", "agb1", "gpa1",
        "control", "mock", "buffer", "treated", "samples", "using", "total", "rna",
        "rep1", "rep2", "rep3", "atgen", "unknown"
    }

    valid_raw_terms = {} 
    
    # Get Golden Keywords once per call
    golden_keywords_list = EXPLICIT_KEYWORDS.get(target_category, [])

    for term in candidate_set:
        term_lower = term.lower()
        
        # A. Noise Filter
        if any(bad in term_lower for bad in NOISE_STOPLIST):
            continue
        if len(term) < 2: 
            continue
        if len(term_lower.split(' ')) > 10:
            continue

        # B. VIP Pass (Golden Keywords)
        # Optimized has_golden_key_word usage
        has_golden_keyword, socre, _, _ = has_golden_key_word(term_lower, golden_keywords_list, optimizer, 0.75)
        
        if has_golden_keyword:
            valid_raw_terms[term] = socre 
            continue

        # --- Standard Vector Logic ---
        _, target_score = optimizer.get_best_match_with_score(term, category=target_category)
        
        if target_score < 0.75:
            continue

        if contrast_category:
            _, contrast_score = optimizer.get_best_match_with_score(term, category=contrast_category)
            if contrast_score > target_score:
                continue

        valid_raw_terms[term] = target_score

    # 3. Sort by Score and Return Top K
    sorted_candidates = sorted(valid_raw_terms.items(), key=lambda x: x[1], reverse=True)
    return sorted_candidates[:top_k]


def has_golden_key_word(string: str, words: list[str], optimized: GroundingOptimizer, th: float) -> Tuple[bool, float, Optional[str], Optional[str]]:
    """
    Checks if any n-gram in the input string is semantically similar to any word in the golden list.
    OPTIMIZED: Minimizes CPU-GPU transfers and leverages efficient tensor ops.
    """
    if not string or not words:
        return False, 0.0, None, None

    # 1. Exact string match (O(1) logic)
    string_lower = string.lower()
    for w in words:
        w_lower = w.lower()
        if w_lower in string_lower:
            return True, 1.0, w_lower, w 

    # 2. Use Cached Vectors
    vectored = False
    for label in LABELS:
        if words is EXPLICIT_KEYWORDS[label]:
            target_embeddings = optimized.vectors['explicit'][label]
            vectored= True
    if not vectored:
        # Fallback: Encode the list on the fly if it's a custom list
        # We use the model's current device (CPU or CUDA) automatically
        target_embeddings = optimized.model.encode(words, convert_to_tensor=True)
    # 3. Generate Candidates (N-grams)
    # Use pre-compiled regex for tokenization
    tokens = TOKEN_BOUNDARY_PATTERN.findall(string_lower)
    if not tokens:
         return False, 0.0, None, None

    candidates = set(tokens)
    
    if len(tokens) > 1:
        candidates.update([" ".join(tokens[i:i+2]) for i in range(len(tokens)-1)])
    if len(tokens) > 2:
        candidates.update([" ".join(tokens[i:i+3]) for i in range(len(tokens)-2)])
    
    candidates_list = list(candidates)
    if not candidates_list:
        return False, 0.0, None, None

    # 4. Encode Candidates
    # batch encoding is efficient here
    candidate_embeddings = optimized.model.encode(candidates_list, convert_to_tensor=True)
    
    # 5. Compute Cosine Similarity Matrix
    # Optimization: Ensure both tensors are on the same device.
    # Removed .cpu() call which forces synchronization and data transfer.
    if target_embeddings.device != candidate_embeddings.device:
        target_embeddings = target_embeddings.to(candidate_embeddings.device)

    cosine_scores = util.cos_sim(candidate_embeddings, target_embeddings)

    # 6. Find Best Match Indices
    # Use torch.max to get values on device immediately
    max_val, max_idx_flat = torch.max(cosine_scores.flatten(), dim=0)
    max_score = max_val.item()
    
    if max_score > th:
        # Calculate indices purely with tensors to stay efficient
        num_cols = cosine_scores.shape[1]
        
        # .item() pulls the scalar to CPU only at the very end
        row_idx = (max_idx_flat // num_cols).item()
        col_idx = (max_idx_flat % num_cols).item()
        
        best_candidate = candidates_list[row_idx]
        best_keyword = words[col_idx]
        
        return True, max_score, best_candidate, best_keyword

    return False, max_score, None, None


class UniversalExtractor:
    def __init__(self):
        self.optimizer = GroundingOptimizer()
        self.known_tissues = {t.value.lower() for t in TissueEnum}
        self.known_treatments = {t.value.lower() for t in TreatmentEnum}
        self.known_mediums = {m.value.lower() for m in MediumEnum}

        self.trigger_keywords = {
            'treatment': ['treatment', 'treated', 'stress', 'condition', 'exposed to', 'exposure', 'incubated', 'temperature', 'growth condition']+ EXPLICIT_KEYWORDS.get('treatment', []),
            'tissue': ['tissue', 'organ', 'source', 'derived from', 'cells', 'cell type', 'organism part']+ EXPLICIT_KEYWORDS.get('tissue', []),
            'medium': ['medium', 'growth medium', 'grown on', 'cultured in', 'substrate']+ EXPLICIT_KEYWORDS.get('medium', [])
        }
        
        self.explicit_vectors = {
            cat: self.optimizer.model.encode(terms) 
            for cat, terms in AREA_KEYWORDS.items()
        }
        self.last_study:str = ''
        self.extracted_data_study = {"tissue": [], "treatment": [], "medium": []}

    def extract(self, sample_metadata: Dict, study_metadata: Dict, study_id:str) -> Dict:
        if study_id != self.last_study:
            extracted_data = {"tissue": set(), "treatment": set(), "medium": set()}
            self.extracted_data_study = extracted_data.copy()
        else:
            extracted_data = self.extracted_data_study.copy()

        for category in ['tissue', 'treatment', 'medium']:
            if extracted_data[category] != set():
                continue
            
            explicit_hits = self._check_for_explicit_columns(sample_metadata, category)
            if explicit_hits:
                extracted_data[category] = explicit_hits
                continue 

            priority_text = self._get_priority_text(sample_metadata, category)
            scanned_hits = self._scan_semantic_match(priority_text, category=category)
            if scanned_hits:
                final_hits = condense_candidates(scanned_hits, self.optimizer, category)
                extracted_data[category] = set(final_hits)
                continue

            broad_text_sample = self._get_broad_text(sample_metadata, {})
            scanned_broad_sample = self._scan_semantic_match(broad_text_sample, category=category)
            if scanned_broad_sample:
                final_hits = condense_candidates(scanned_broad_sample, self.optimizer, category)
                extracted_data[category] = set(final_hits)
                continue

            broad_text = self._get_broad_text({}, study_metadata)
            scanned_broad = self._scan_semantic_match(broad_text, category=category)
            final_hits = condense_candidates(scanned_broad, self.optimizer, category)
            extracted_data[category] = set(final_hits)

        self.last_study = study_id
        return extracted_data

    def _get_priority_text(self, metadata: Dict, category: str) -> str:
        text_bits = []
        text_bits.append(str(metadata.get('title', '')))
        text_bits.append(str(metadata.get('source_name_ch1', '')))
        text_bits.append(str(metadata.get('characteristics_ch1', '')))
        text_bits.append(str(metadata.get('growth_protocol_ch1', '')))
        text_bits.append(str(metadata.get('characteristics', '')))
        for key, val in metadata.items():
            if self._is_key_semantic_match(key, category):
                if isinstance(val, list): text_bits.extend([str(v) for v in val])
                else: text_bits.append(str(val))
        return " ".join(text_bits)
    
    # Optimization: Cache the semantic match results.
    # Metadata keys (e.g., "title", "source_name_ch1") repeat thousands of times across samples.
    # This prevents re-encoding the same short string over and over.
    @lru_cache(maxsize=1024)
    def _is_key_semantic_match(self, key_text: str, category: str, threshold: float = 0.85) -> bool:
        clean_key = key_text.replace('_', ' ').replace('-', ' ')
        clean_key = re.sub(r'(?<!^)(?=[A-Z])', ' ', clean_key).lower().strip()
        
        # Note: We access self.optimizer via self (this works fine with lru_cache on instance methods in Python 3.8+)
        # However, to be strictly safe and efficient, we assume self.optimizer doesn't change.
        key_vector = self.optimizer.model.encode(clean_key)
        
        # Use appropriate category vector from cache
        category_vectors = self.explicit_vectors[category]
        
        scores = util.cos_sim(key_vector, category_vectors)[0]
        return bool(scores.max() > threshold)
    
    def _check_for_explicit_columns(self, metadata: Dict, category: str) -> Set[str]:
        """
        Scans metadata columns to find values that are likely the category we want.
        """
        found_values = set()
        generic_container_keys = {
            'source_name_ch1', 'characteristics_ch1', 'characteristics', 
            'growth_protocol_ch1', 'description', 'title'
        }

        for key, val in metadata.items():
            raw_values = val if isinstance(val, list) else [str(val)]
            
            # STRATEGY A: Explicit Column Name Match
            if self._is_key_semantic_match(key, category, threshold=0.85):
                found_values.update([str(v) for v in raw_values])
                continue 

            # STRATEGY B: Generic Column Parsing
            if key in generic_container_keys:
                for item in raw_values:
                    item_str = str(item).strip()
                    
                    if ':' in item_str:
                        parts = item_str.split(':', 1)
                        sub_key = parts[0].strip()
                        sub_val = parts[1].strip().strip('"').strip("'")
                        
                        if self._is_key_semantic_match(sub_key, category, threshold=0.90):
                            if len(sub_val) > 1:
                                found_values.add(sub_val)
                                continue 

        return found_values

    def _get_broad_text(self, sample_meta: Dict, study_meta: Dict) -> str:
        text_bits = []
        text_bits.append(str(sample_meta.get('description', '')))
        text_bits.append(str(sample_meta.get('growth_protocol_ch1', '')))
        text_bits.append(str(sample_meta.get('treatment_protocol_ch1', '')))
        text_bits.append(str(study_meta.get('summary', '')))
        text_bits.append(str(study_meta.get('overall_design', '')))
        return " ".join(text_bits)

    def _scan_for_keyword_contexts(self, text: str, category: str) -> Set[str]:
        candidates = set()
        triggers = self.trigger_keywords.get(category, [])
        delimiters = r",:;\[\]\"'()."
        
        for kw in triggers:
            # Note: Regex construction inside loop is necessary as 'kw' changes.
            pattern = fr'(?:^|[{delimiters}])\s*([^{delimiters}]*?{re.escape(kw)}[^{delimiters}]*)'
            matches = re.findall(pattern, text, re.IGNORECASE)
            
            for m in matches:
                cleaned = m.strip()
                cleaned = cleaned.strip('-_')
                if 2 < len(cleaned) < 120: 
                    candidates.add(cleaned)
                    
        return candidates

    def _scan_semantic_match(self, text: str, category: str, study_text:bool = False) -> Set[str]:
        # --- PHASE A: Keyword-Anchored Context Scan ---
        raw_context_candidates_alt:set = set(complex_split(text))
        
        valid_context_tuples = extract_valid_candidates(
            candidate_set=raw_context_candidates_alt,
            optimizer=self.optimizer,
            target_category=category,
            top_k=5,
            study_info = study_text
        )
        valid_context_strings = [x[0] for x in valid_context_tuples]

        # --- PHASE B: N-Gram Scan (Broad Scan) ---
        # Optimization: Use pre-compiled regex for cleaning
        clean_text = NOISE_FILTER_PATTERN.sub(' ', text.lower())
        tokens = clean_text.split()
        tokens = tokens[:600] 
        
        STOPWORDS = {
            "the", "and", "of", "in", "to", "with", "for", "on", "at", "by", "from", 
            "was", "were", "are", "is", "an", "a", "or", "that", "this", "using", 
            "analysis", "data", "samples", "expression", "profiling", "study", "results",
            "between", "under", "during", "after", "before", "total", "rna"
        }

        raw_ngram_candidates = set()
        max_n = 3 
        
        # This loop is generally fast in Python for < 2000 items, manual optimization often yields diminishing returns vs readability.
        for n in range(1, max_n + 1):
            for i in range(len(tokens) - n + 1):
                chunk_tokens = tokens[i : i + n]
                
                if chunk_tokens[0] in STOPWORDS or chunk_tokens[-1] in STOPWORDS:
                    continue
                if any(len(t) < 2 for t in chunk_tokens):
                    continue

                chunk_str = " ".join(chunk_tokens)
                raw_ngram_candidates.add(chunk_str)

        if not raw_ngram_candidates and not valid_context_strings:
            return set()

        valid_ngram_tuples = extract_valid_candidates(
            candidate_set=raw_ngram_candidates,
            optimizer=self.optimizer,
            target_category=category,
            top_k=5,
            study_info = study_text
        )
        valid_ngram_strings = [x[0] for x in valid_ngram_tuples]


        # --- PHASE C: Merge and Condense ---
        merged_candidates = []
        seen = set()
        
        if valid_context_strings:
            for cand in valid_context_strings:
                if cand not in seen:
                    merged_candidates.append(cand)
                    seen.add(cand)
            return set(valid_context_strings) # Early return prioritizes context
        
        if valid_ngram_strings:
            for cand in valid_ngram_strings:
                if cand not in seen:
                    merged_candidates.append(cand)
                    seen.add(cand)

        return condense_candidates(merged_candidates, self.optimizer, category)