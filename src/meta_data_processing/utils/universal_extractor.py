import re
from typing import Dict, List, Any, Optional, Set, Tuple
from sentence_transformers import util
from src.constants import TissueEnum, TreatmentEnum, MediumEnum, GenotypeEnum, EXPLICIT_KEYWORDS, AREA_KEYWORDS, LABELS
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
    if target_category == 'treatment':  contrast_category = 'tissue'
    elif target_category == 'tissue':   contrast_category = 'treatment'
    elif target_category == 'medium':   contrast_category = 'tissue'
    elif target_category == 'genotype': contrast_category = 'treatment'
    else: contrast_category = None

    # Hard blocklist — whole-word matching only (prevents "rep1" blocking "replicate 1",
    # and prevents dropping "mock inoculation" because it contains "mock").
    # NOTE: genotype-signal terms ("col-0", "mutant", "transgenic", "columbia",
    # "wild type", "wild-type") are intentionally absent — they are valid signal
    # for the genotype label and are handled by extract_valid_candidates' contrast
    # scoring when target_category != 'genotype'.
    NOISE_TOKENS = {
        "plants", "analysis", "data", "replicate", "study", "experiment",
        "grown", "growth", "independent", "agb1", "gpa1",
        "rep1", "rep2", "rep3", "atgen", "unknown",
        "total", "rna", "using", "samples", "buffer",
    }

    # Soft-strip words: present in a phrase but not the meaningful part.
    # Stripped before scoring so "treated with NaCl" -> "NaCl" gets a fair chance,
    # but the *original* term is kept as the candidate key for the grounder.
    STRIP_WORDS = {"treated", "treatment", "with", "under", "by", "the", "a", "an",
                   "condition", "conditions", "exposure"}

    valid_raw_terms: Dict[str, float] = {}
    golden_keywords_list = EXPLICIT_KEYWORDS.get(target_category, [])

    for term in candidate_set:
        term_lower = term.lower()
        term_words = set(re.findall(r'\b\w+\b', term_lower))

        # Whole-word noise block (not substring — avoids over-filtering)
        if term_words & NOISE_TOKENS:
            continue
        if len(term) < 2:
            continue
        if len(term_lower.split()) > 10:
            continue

        # Soft-strip then score; keep original as key for grounder
        stripped_words = [w for w in term_lower.split() if w not in STRIP_WORDS]
        eval_term = " ".join(stripped_words).strip() or term_lower

        # Golden keyword fast-path (exact or semantic synonym match)
        has_golden_keyword, score, _, _ = has_golden_key_word(eval_term, golden_keywords_list, optimizer, 0.75)

        if has_golden_keyword:
            valid_raw_terms[term] = score
            continue

        _, target_score = optimizer.get_best_match_with_score(eval_term, category=target_category)

        # Per-label threshold: medium terms are short and lexically unusual,
        # so a slightly lower bar prevents "MS", "B5", "agar" from being dropped.
        per_label_threshold = 0.70 if target_category == 'medium' else 0.75
        if target_score < per_label_threshold:
            continue

        if contrast_category:
            _, contrast_score = optimizer.get_best_match_with_score(eval_term, category=contrast_category)
            if contrast_score > target_score: continue

        valid_raw_terms[term] = target_score

    sorted_candidates = sorted(valid_raw_terms.items(), key=lambda x: x[1], reverse=True)
    return sorted_candidates[:top_k]

def has_golden_key_word(string: str, words: list[str], optimized: GroundingOptimizer, th: float) -> Tuple[bool, float, Optional[str], Optional[str]]:
    if not string or not words: return False, 0.0, None, None

    string_lower = string.lower()
    # FIX: return the longest (most specific) substring match rather than the first.
    # Mirrors the fix applied to find_semantic_match in groundingOptimizer.py.
    substring_matches = [(w, w.lower()) for w in words if w.lower() in string_lower]
    if substring_matches:
        best_w, best_w_lower = max(substring_matches, key=lambda x: len(x[1]))
        return True, 1.0, best_w_lower, best_w

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
            all_hits: Set[str] = set()

            # Source 1: Explicit column key:value parsing
            col_hits = self._extract_from_matched_columns(sample_metadata, label_type)
            all_hits.update(col_hits)

            # Source 2: Priority text blob from sample metadata
            priority_text = self._get_text_blob(sample_metadata, ['title', 'characteristics_ch1', 'source_name_ch1'])
            if priority_text:
                all_hits.update(self._semantic_ngram_search(priority_text, label_type))

            # Source 3: Study-level text — always searched for medium/genotype since
            # those are often only described in the study summary, not per-sample fields.
            # For treatment/tissue, only supplement when sample level found nothing useful.
            study_text = self._get_text_blob(study_metadata, ['summary', 'overall_design'])
            if study_text:
                if label_type in ('medium', 'genotype') or not all_hits:
                    all_hits.update(self._semantic_ngram_search(study_text, label_type))

            # Post-filter: contrast-category disambiguation + noise removal
            if all_hits:
                valid = extract_valid_candidates(
                    candidate_set=all_hits,
                    optimizer=self.optimizer,
                    target_category=label_type,
                    study_info=bool(study_metadata),
                    top_k=5
                )
                final_hits = {term for term, _ in valid}
            else:
                final_hits = set()

            extracted[label_type] = final_hits if final_hits else {"unspecified"}

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

        # Pre-build a flat set of all trigger words for this label for fast sub-key bypass
        trigger_words = set(w.lower() for w in LABEL_CONFIG[label_type].get('search_triggers', []))

        for key, value in metadata.items():
            vals = value if isinstance(value, list) else [value]
            
            # --- CASE 1: The key is a priority column (e.g., 'characteristics_ch1') ---
            if self._is_priority_column(key, label_type):
                for val in vals:
                    val_str = str(val).strip()
                    val_str_clean = val_str.replace('_', ' ')
                    
                    if ":" in val_str_clean:
                        sub_key, actual_val = val_str_clean.split(":", 1)
                        sub_key_clean = sub_key.strip().lower()
                        actual_val = actual_val.strip()

                        # FIX: fast-path bypass — if ANY trigger word appears in the sub-key,
                        # skip the vector gate entirely. This catches sub-keys like
                        # "growth condition", "treatment protocol", "experimental factor"
                        # that were scoring below threshold and silently dropping valid values.
                        sub_key_has_trigger = any(t in sub_key_clean for t in trigger_words)

                        if sub_key_has_trigger or self._is_relevant_element(sub_key_clean, label_type, threshold=0.60):
                            if len(actual_val.split()) > 3:
                                ngrams = self._semantic_ngram_search(actual_val, label_type, threshold=0.75)
                                found.update(ngrams)
                            else:
                                found.add(actual_val)
                    else:
                        if self._is_relevant_element(val_str_clean, label_type):
                            if len(val_str_clean.split()) > 3:
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
        
        # Noise regex strips words that pollute treatment/tissue/medium searches.
        # For genotype, these words ARE the signal — skip stripping entirely.
        if label_type != 'genotype':
            noise_words = r'\b(genotype|mutant|wild-type|wildtype|age|weeks|days|old|developmental stage|cell type)\b'
            clean_text = re.sub(noise_words, '', text, flags=re.IGNORECASE)
        else:
            clean_text = text
        
        # Break apart squished string formatting
        # e.g. "AtGen_6-2521_Osmoticstress-Roots" -> "AtGen 6 2521 Osmoticstress Roots"
        clean_text = clean_text.replace('_', ' ').replace('-', ' ')
        
        # Tokenize (keeping forward slash for things like "1/2 MS")
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

    def _semantic_ngram_search_old(self, text: str, label_type: str, threshold: float = 0.80) -> Set[str]:
        """
        Chunks text into n-grams and compares them to the valid ontology values
        defined in EXPLICIT_KEYWORDS.
        """
        if not text: return set()
        
        # 1. Strip out highly common, distracting prefixes/words that confuse the vectorizer
        noise_words = r'\b(genotype|mutant|wild-type|wildtype|age|weeks|days|old|developmental stage|cell type)\b' #TODO move this to the constants file
        clean_text = re.sub(noise_words, '', text, flags=re.IGNORECASE)
        
        words = re.findall(r'\b[\w/-]+\b', clean_text)
        if not words: return set()

        # Generate N-grams (1 to 3 words)
        candidates = []
        for n in range(1, 4):
            for i in range(len(words) - n + 1):
                candidates.append(" ".join(words[i:i+n]))
        
        if not candidates: return set()

        # Deduplicate and Encode Candidates in one batch for speed
        unique_candidates = list(set(candidates))
        candidate_vectors = self.model.encode(unique_candidates, convert_to_tensor=True)
        
        # Get target vectors (The valid synonyms/canonical names for this label)
        # These are pre-computed in your GroundingOptimizer
        target_vectors = self.optimizer.vectors['explicit'][label_type]
        
        # Compute Similarity Matrix (Candidates x Ontology)
        cosine_matrix = util.cos_sim(candidate_vectors, target_vectors)
        
        # Find matches above threshold
        matches = torch.where(cosine_matrix > threshold)
        
        results = set()
        for idx in matches[0]:
            # We return the raw string found in text. 
            # The grounding step will later map this to the Enum.
            results.add(unique_candidates[idx.item()])
            
        return results