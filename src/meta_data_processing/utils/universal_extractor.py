import re
from typing import Dict, List, Any, Optional, Set, Tuple
from sentence_transformers import util
from src.constants import TissueEnum, TreatmentEnum, MediumEnum, EXPLICIT_KEYWORDS, AREA_KEYWORDS, LABELS
from src.meta_data_processing.utils.groundingOptimizer import GroundingOptimizer
import sys
import torch
from functools import lru_cache

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


class UniversalExtractor:
    def __init__(self):
        self.optimizer = GroundingOptimizer()
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

        for category in LABELS:
            if extracted_data[category] != set():
                continue
            
            # 1. Explicit Columns (High Confidence)
            explicit_hits = self._check_for_explicit_columns(sample_metadata, category)
            if explicit_hits:
                # IMPORTANT: Canonicalize the explicit hits! 
                # If column says "Dark", we save "Low Light Stress"
                canonical_hits = {self.optimizer.canonicalize_term(h, category) for h in explicit_hits}
                extracted_data[category] = canonical_hits
                continue 

            # 2. Priority Text Scan
            priority_text = self._get_priority_text(sample_metadata, category)
            scanned_hits = self._scan_semantic_match(priority_text, category=category)
            if scanned_hits:
                final_hits = condense_candidates(scanned_hits, self.optimizer, category)
                extracted_data[category] = set(final_hits)
                continue

            # 3. Broad Sample Text Scan
            broad_text_sample = self._get_broad_text(sample_metadata, {})
            scanned_broad_sample = self._scan_semantic_match(broad_text_sample, category=category)
            if scanned_broad_sample:
                final_hits = condense_candidates(scanned_broad_sample, self.optimizer, category)
                extracted_data[category] = set(final_hits)
                continue

            # 4. Study Level Fallback
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
    
    @lru_cache(maxsize=1024)
    def _is_key_semantic_match(self, key_text: str, category: str, threshold: float = 0.85) -> bool:
        clean_key = key_text.replace('_', ' ').replace('-', ' ')
        clean_key = re.sub(r'(?<!^)(?=[A-Z])', ' ', clean_key).lower().strip()
        key_vector = self.optimizer.model.encode(clean_key)
        category_vectors = self.explicit_vectors[category]
        scores = util.cos_sim(key_vector, category_vectors)[0]
        return bool(scores.max() > threshold)
    
    def _check_for_explicit_columns(self, metadata: Dict, category: str) -> Set[str]:
        found_values = set()
        generic_container_keys = {
            'source_name_ch1', 'characteristics_ch1', 'characteristics', 
            'growth_protocol_ch1', 'description', 'title'
        }

        for key, val in metadata.items():
            raw_values = val if isinstance(val, list) else [str(val)]
            
            if self._is_key_semantic_match(key, category, threshold=0.85):
                found_values.update([str(v) for v in raw_values])
                continue 

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

    def _scan_semantic_match(self, text: str, category: str, study_text:bool = False) -> Set[str]:
        # PHASE A: Keyword-Anchored
        raw_context_candidates_alt = set(complex_split(text))
        valid_context_tuples = extract_valid_candidates(raw_context_candidates_alt, self.optimizer, category, study_text, top_k=5)
        valid_context_strings = [x[0] for x in valid_context_tuples]

        # PHASE B: N-Gram Scan
        clean_text = NOISE_FILTER_PATTERN.sub(' ', text.lower())
        tokens = clean_text.split()[:600] 
        STOPWORDS = {"the", "and", "of", "in", "to", "with", "for", "on", "at", "by", "from", "was", "were", "are", "is", "an", "a", "or", "using", "analysis", "data", "samples", "study", "results", "between", "under", "total", "rna"}

        raw_ngram_candidates = set()
        max_n = 3 
        
        for n in range(1, max_n + 1):
            for i in range(len(tokens) - n + 1):
                chunk_tokens = tokens[i : i + n]
                if chunk_tokens[0] in STOPWORDS or chunk_tokens[-1] in STOPWORDS: continue
                if any(len(t) < 2 for t in chunk_tokens): continue
                raw_ngram_candidates.add(" ".join(chunk_tokens))

        if not raw_ngram_candidates and not valid_context_strings:
            return set()

        valid_ngram_tuples = extract_valid_candidates(raw_ngram_candidates, self.optimizer, category, study_text, top_k=5)
        valid_ngram_strings = [x[0] for x in valid_ngram_tuples]

        # PHASE C: Merge
        merged_candidates = []
        seen = set()
        
        if valid_context_strings:
            for cand in valid_context_strings:
                if cand not in seen:
                    merged_candidates.append(cand)
                    seen.add(cand)
            return set(valid_context_strings) 
        
        if valid_ngram_strings:
            for cand in valid_ngram_strings:
                if cand not in seen:
                    merged_candidates.append(cand)
                    seen.add(cand)

        return condense_candidates(merged_candidates, self.optimizer, category)