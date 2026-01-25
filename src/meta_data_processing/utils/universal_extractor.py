import re
from typing import Dict, List, Any, Optional, Set, Tuple
from sentence_transformers import util
from src.constants import TissueEnum, TreatmentEnum, MediumEnum
from src.meta_data_processing.utils.classes import GroundingOptimizer

# In src/universal_extractor.py
def condense_candidates(candidates: List[str], optimizer: GroundingOptimizer, category: str) -> List[str]:
    """
    Selects the single best term from a list of overlapping candidates 
    based on Semantic Purity.
    """
    if not candidates:
        return []

    # 1. Cleaning Step: Remove obvious garbage (pure numbers, short codes)
    # This turns "2521 osmotic stress" -> "osmotic stress"
    cleaned_candidates = set()
    for cand in candidates:
        # Regex: Remove leading/trailing digits, genotypes (Col-0), and noise
        # Leaves the core biological phrase
        clean = re.sub(r'^\d+\s+|\s+\d+$|\b(rep\d+|col-0|atgen)\b', '', cand.lower()).strip()
        if len(clean) > 2:
            cleaned_candidates.add(clean)

    if not cleaned_candidates:
        return []

    # 2. Semantic Scoring Competition
    best_term = None
    best_score = -1.0
    results = []

    for term in cleaned_candidates:
        # Ask BioLORD: "How much does this look like a Treatment?"
        # 'osmotic stress' -> 0.85
        # 'osmotic stress roots' -> 0.72 (Polluted by Tissue term)
        label, score = optimizer.get_best_match_with_score(term, category)
        results.append((label,score))
        
        # Heuristic: Prefer shorter strings if scores are very close (within 0.02)
        # This avoids picking "osmotic stress condition" over "osmotic stress"
        if score > best_score:
            best_score = score
            best_term = term
        elif abs(score - best_score) < 0.02:
            if len(term) < len(best_term):
                best_term = term
    results.sort(key=lambda tup: tup[1],reverse=True)
    # if results[0][1]> 0.75 and results[1][1]> 0.75 and results[0][0] !=results[1][0]:
    #     return [results[0][0],results[1][0]]

    return [best_term] if best_term else []
def extract_valid_candidates(candidate_set: Set[str], optimizer: GroundingOptimizer, target_category: str, top_k: int = 3) -> List[Tuple[str, float]]:
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
        "rep1", "rep2", "rep3", "atgen"
    }

    # 2. Define Golden Keywords (The VIP List)
    # If the text contains these, we trust the text over the vector model.
    golden_keywords = []
    if target_category == 'treatment':
        golden_keywords = ['stress', 'treatment', 'condition', 'exposed', 'drought', 'heat', 'salt', 'cold', 'mock', 'aba', 'nacl']
    elif target_category == 'tissue':
        golden_keywords = ['tissue', 'organ', 'root', 'leaf', 'shoot', 'flower', 'seedling', 'rosette', 'guard cell']
    elif target_category == 'medium':
        golden_keywords = ['medium', 'agar', 'soil', 'plate', 'hydroponic', 'ms']

    valid_raw_terms = {} 

    for term in candidate_set:
        term_lower = term.lower()
        
        # A. Noise Filter
        if any(bad in term_lower for bad in NOISE_STOPLIST):
            continue
        # Length check (relaxed for short codes like "ABA")
        if len(term) < 2: 
            continue

        # B. VIP Pass (Golden Keywords)
        # If the term literally contains "stress" or "root", we KEEP it immediately.
        # We assume explicit text > vector semantics for composite strings like "Root_Stress_12h"
        has_golden_keyword = any(kw in term_lower for kw in golden_keywords)
        
        if has_golden_keyword:
            # Assign a "Perfect" score to ensure it floats to the top
            # We skip the contrast check entirely because "Root_Stress" matches both, 
            # and we don't want the Tissue vector to kill the Treatment candidate.
            valid_raw_terms[term] = 0.99
            continue

        # --- Standard Vector Logic (For terms without keywords) ---
        
        # C. Get Target Score
        _, target_score = optimizer.get_best_match_with_score(term, category=target_category)
        
        # D. Threshold Check (Strict for non-keywords)
        if target_score < 0.82:
            continue

        # E. Contrast Check (Only for ambiguous non-keyword terms)
        if contrast_category:
            _, contrast_score = optimizer.get_best_match_with_score(term, category=contrast_category)
            
            # If explicit keyword is missing, we rely on the tug-of-war
            if contrast_score > target_score:
                continue

        valid_raw_terms[term] = target_score

    # 3. Sort by Score and Return Top K
    sorted_candidates = sorted(valid_raw_terms.items(), key=lambda x: x[1], reverse=True)
    return sorted_candidates[:top_k]

def extract_valid_candidates_old(candidate_set: Set[str], optimizer: GroundingOptimizer, target_category: str, top_k: int = 3) -> List[Tuple[str, float]]:
    """
    Filters candidates using Vector Similarity + Keyword Boosting.
    Prioritizes terms containing category-specific substrings (e.g. 'stress').
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
        "rep1", "rep2", "rep3" # Added specific replicate markers from your example
    }

    # --- NEW: Golden Keywords for Boosting ---
    # If these substrings appear, we boost the score significantly.
    # This catches "osmoticstress" (contains 'stress') even if the vector is weak.
    golden_keywords = []
    if target_category == 'treatment':
        golden_keywords = ['stress', 'treatment', 'condition', 'exposed', 'drought', 'heat', 'salt', 'cold']
    elif target_category == 'tissue':
        golden_keywords = ['tissue', 'organ', 'root', 'leaf', 'shoot', 'flower', 'seedling']
    elif target_category == 'medium':
        golden_keywords = ['medium', 'agar', 'soil', 'plate', 'hydroponic']

    valid_raw_terms = {} 

    for term in candidate_set:
        term_lower = term.lower()
        
        # A. Noise Filter
        if any(bad in term_lower for bad in NOISE_STOPLIST):
            continue
        if len(term) < 3: 
            continue

        # B. Check for Keyword Boost
        # We use a simple substring check so "osmoticstress" triggers on "stress"
        has_golden_keyword = any(kw in term_lower for kw in golden_keywords)

        # C. Get Target Score
        _, target_score = optimizer.get_best_match_with_score(term, category=target_category)
        
        # --- LOGIC CHANGE: Apply Boost ---
        # If it has a keyword, we artificially inflate the score.
        # This helps it survive the threshold AND win the contrast check.
        effective_score = target_score
        if has_golden_keyword:
            effective_score += 0.15 # Significant boost (e.g., 0.70 -> 0.85)

        # D. Threshold Check
        # If we have a keyword, we tolerate lower vector quality (0.70). 
        # Otherwise, we keep the strict 0.82.
        threshold = 0.70 if has_golden_keyword else 0.82
        
        if effective_score < threshold:
            continue

        # E. Contrast Check (The Distractor Trap)
        if contrast_category:
            _, contrast_score = optimizer.get_best_match_with_score(term, category=contrast_category)
            
            # If the term has a golden keyword for the TARGET, we ignore the contrast
            # unless the contrast score is overwhelmingly high (e.g. 0.98).
            # Example: "Osmoticstress-Roots"
            #   Target (Treatment): 0.75 + 0.15 Boost = 0.90
            #   Contrast (Tissue): 0.85 (matches 'roots')
            #   Result: 0.90 > 0.85 -> It stays as a Treatment candidate.
            if not has_golden_keyword and contrast_score > effective_score:
                continue
            elif has_golden_keyword and contrast_score > (effective_score + 0.1):
                # Only discard a boosted term if the contrast is MUCH stronger
                continue

        # Save the EFFECTIVE score so boosted terms float to the top of the list
        valid_raw_terms[term] = effective_score

    # 3. Sort by Score and Return Top K
    sorted_candidates = sorted(valid_raw_terms.items(), key=lambda x: x[1], reverse=True)
    return sorted_candidates[:top_k]

class UniversalExtractor:
    def __init__(self):
        self.optimizer = GroundingOptimizer()
        self.known_tissues = {t.value.lower() for t in TissueEnum}
        self.known_treatments = {t.value.lower() for t in TreatmentEnum}
        self.known_mediums = {m.value.lower() for m in MediumEnum}

        # --- Anchor Vectors (Keep existing) ---
        anchors = {
            'tissue': ['tissue', 'organ', 'source', 'cell type', 'tissue type', 'organism part'],
            'treatment': ['treatment', 'stress', 'condition', 'treatment protocol', 'treated with', 'exposure'],
            'medium': ['medium', 'growth medium', 'substrate', 'growth condition', 'culture medium']
        }
        self.explicit_vectors = {
            cat: self.optimizer.model.encode(terms) 
            for cat, terms in anchors.items()
        }

        # --- NEW: Trigger Keywords for Context Extraction ---
        # If these words appear in a sentence, we extract the whole clause.
        self.trigger_keywords = {
            'treatment': ['treatment', 'treated', 'stress', 'condition', 'exposed to', 'incubated'],
            'tissue': ['tissue', 'organ', 'source', 'derived from', 'cells'],
            'medium': ['medium', 'growth medium', 'grown on', 'cultured in', 'substrate']
        }
        self.last_study:str = ''
        self.extracted_data_study = {"tissue": [], "treatment": [], "medium": []}

    # ... (extract, _is_semantic_key_match, _check_for_explicit_columns, _get_priority_text, _get_broad_text methods remain the same) ...
    def extract(self, sample_metadata: Dict, study_metadata: Dict,study_id:str) -> Dict:
        # (This remains unchanged from previous version)
        if study_id != self.last_study:
            extracted_data = {"tissue": [], "treatment": [], "medium": []}
            self.extracted_data_study = extracted_data.copy()
        else:
            extracted_data = self.extracted_data_study.copy()
        for category in ['tissue', 'treatment', 'medium']:
            if extracted_data[category] !=[]:
                continue
            explicit_hits = self._check_for_explicit_columns(sample_metadata, category)
            if explicit_hits:
                extracted_data[category] = explicit_hits
                continue 

            priority_text = self._get_priority_text(sample_metadata, category)
            scanned_hits = self._scan_semantic_match(priority_text, category=category)
            if scanned_hits:
                final_hits = condense_candidates(scanned_hits,self.optimizer,category)
                extracted_data[category] = list(set(final_hits))
                continue
            broad_text_sample = self._get_broad_text(sample_metadata, {})
            scanned_broad_sample = self._scan_semantic_match(broad_text_sample, category=category)
            if scanned_broad_sample:
                final_hits = condense_candidates(scanned_broad_sample,self.optimizer,category)
                extracted_data[category] = list(set(final_hits))
                continue

            broad_text = self._get_broad_text({}, study_metadata)
            scanned_broad = self._scan_semantic_match(broad_text, category=category)
            final_hits = condense_candidates(scanned_broad,self.optimizer,category)
            extracted_data[category] = list(set(final_hits))
            #here we learned from the study, hence it sahould be constant across samples
            self.extracted_data_study[category] = list(set(final_hits))
        self.last_study = study_id
        return extracted_data

    def _get_priority_text(self, metadata: Dict, category: str) -> str:
        # (Keep existing)
        text_bits = []
        text_bits.append(str(metadata.get('title', '')))
        text_bits.append(str(metadata.get('source_name_ch1', '')))
        for key, val in metadata.items():
            if self._is_semantic_key_match(key, category):
                if isinstance(val, list): text_bits.extend([str(v) for v in val])
                else: text_bits.append(str(val))
        return " ".join(text_bits)
    
    def _is_semantic_key_match(self, key_text: str, category: str, threshold: float = 0.85) -> bool:
        # (Keep existing)
        clean_key = key_text.replace('_', ' ').replace('-', ' ')
        clean_key = re.sub(r'(?<!^)(?=[A-Z])', ' ', clean_key).lower().strip()
        key_vector = self.optimizer.model.encode(clean_key)
        scores = util.cos_sim(key_vector, self.explicit_vectors[category])[0]
        return bool(scores.max() > threshold)
    
    def _check_for_explicit_columns(self, metadata: Dict, category: str) -> List[str]:
        # (Keep existing)
        found_values = []
        for key, val in metadata.items():
            if self._is_semantic_key_match(key, category, threshold=0.85):
                if isinstance(val, list): found_values.extend([str(v) for v in val])
                else: found_values.append(str(val))
        target_list_keys = ['characteristics_ch1', 'characteristics', 'source_name_ch1']
        for list_key in target_list_keys:
            if list_key in metadata:
                val_list = metadata[list_key]
                if not isinstance(val_list, list): val_list = [str(val_list)]
                for item in val_list:
                    item_str = str(item)
                    if ':' in item_str:
                        parts = item_str.split(':', 1)
                        key_part = parts[0].strip()
                        val_part = parts[1].strip()
                        if self._is_semantic_key_match(key_part, category):
                            val_part = val_part.strip('"').strip("'")
                            if len(val_part) > 1:
                                found_values.append(val_part)
        return list(set(found_values))

    def _get_broad_text(self, sample_meta: Dict, study_meta: Dict) -> str:
        # (Keep existing)
        text_bits = []
        text_bits.append(str(sample_meta.get('description', '')))
        text_bits.append(str(sample_meta.get('growth_protocol_ch1', '')))
        text_bits.append(str(sample_meta.get('treatment_protocol_ch1', '')))
        text_bits.append(str(study_meta.get('summary', '')))
        text_bits.append(str(study_meta.get('overall_design', '')))
        return " ".join(text_bits)

    def _scan_for_keyword_contexts(self, text: str, category: str) -> Set[str]:
        """
        Uses Regex to find clauses containing trigger words, even if embedded in 
        identifiers without spaces (e.g. 'Osmoticstress').
        """
        candidates = set()
        triggers = self.trigger_keywords.get(category, [])
        
        # 1. Define Delimiters
        # We treat these chars as Hard Boundaries for a context chunk.
        # We include [ ] ' " to handle Python list string representations often found in metadata.
        # We include . , ; : as standard sentence delimiters.
        # We do NOT include _ or - because we want to capture "Osmotic-stress" as one unit.
        delimiters = r",:;\[\]\"'()."
        
        for kw in triggers:
            # Construct Regex
            # 1. (?:^|[{delimiters}]) -> Start at beginning of string OR after a delimiter
            # 2. \s* -> Optional whitespace
            # 3. ([^{delimiters}]*? ... [^{delimiters}]*) -> Capture group:
            #    Matches content that is NOT a delimiter, containing the keyword {kw}
            #    We REMOVED \b to allow matching "stress" inside "Osmoticstress"
            
            pattern = fr'(?:^|[{delimiters}])\s*([^{delimiters}]*?{re.escape(kw)}[^{delimiters}]*)'
            
            matches = re.findall(pattern, text, re.IGNORECASE)
            
            for m in matches:
                cleaned = m.strip()
                # Clean up leading/trailing connectors often captured at the edge
                cleaned = cleaned.strip('-_')
                
                # Length check:
                # We allow slightly longer chunks (up to 120) to accommodate long file IDs
                if 2 < len(cleaned) < 120: 
                    candidates.add(cleaned)
                    
        return candidates

    def _scan_semantic_match_old(self, text: str, category: str) -> List[str]:
        if not self.optimizer or not text.strip(): return []
        
        # 1. Clean Text lightly (keep punctuation for context scanning)
        # We need the punctuation for the _scan_for_keyword_contexts logic
        # But we create a cleaned version for n-grams
        
        candidates = set()

        # --- PHASE A: Keyword-Anchored Context Scan ---
        # This catches "no treatment" or "abscisic acid treatment"
        context_candidates = self._scan_for_keyword_contexts(text, category)
        candidates.update(context_candidates)

        # --- PHASE B: N-Gram Scan (The original logic) ---
        clean_text = re.sub(r'[^a-zA-Z0-9\s]', ' ', text.lower())
        tokens = clean_text.split()
        tokens = tokens[:600] 
        STOPWORDS = {"the", "and", "of", "in", "to", "with", "for", "on", "at", "by", "from", "was", "were", "are", "is", "an", "a", "or", "that", "this", "using", "analysis", "data", "samples", "expression", "profiling", "study", "results", "between", "under", "during", "after", "before", "total", "rna"}
        
        max_n = 3 
        for n in range(1, max_n + 1):
            for i in range(len(tokens) - n + 1):
                chunk_tokens = tokens[i : i + n]
                if chunk_tokens[0] in STOPWORDS or chunk_tokens[-1] in STOPWORDS: continue
                if any(len(t) < 2 for t in chunk_tokens): continue
                candidates.add(" ".join(chunk_tokens))

        if not candidates: return []
        
        # 3. Filter Candidates using the Optimizer
        # Pass ALL candidates (Context + N-grams) to the semantic filter
        top_candidates = extract_valid_candidates(candidates, self.optimizer, category, top_k=5)
        if top_candidates == []:
            # IF THERE ARE NO CONDIDATE
            text_split = text.replace('[',']')
            text_split = text_split.split(']')
            for te in text_split:
                # TODO??
                pass
        return [item[0] for item in top_candidates]
    
    def _scan_semantic_match(self, text: str, category: str) -> List[str]:
        
        # --- PHASE A: Keyword-Anchored Context Scan ---
        # 1. Extract candidates based on regex triggers (e.g., "treatment: [X]")
        raw_context_candidates = self._scan_for_keyword_contexts(text, category)
        
        # 2. Process/Filter Context candidates immediately
        # We assume these are higher fidelity, so we check them in isolation.
        valid_context_tuples = extract_valid_candidates(
            candidate_set=raw_context_candidates,
            optimizer=self.optimizer,
            target_category=category,
            top_k=5 # Keep top 5 specific context hits
        )
        # Extract just the string strings
        valid_context_strings = [x[0] for x in valid_context_tuples]


        # --- PHASE B: N-Gram Scan (Broad Scan) ---
        # 1. Clean Text (keep punctuation for context scanning, but remove for n-grams)
        clean_text = re.sub(r'[^a-zA-Z0-9\s]', ' ', text.lower())
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
        for n in range(1, max_n + 1):
            for i in range(len(tokens) - n + 1):
                chunk_tokens = tokens[i : i + n]
                
                # Filter A: Stopword boundaries
                if chunk_tokens[0] in STOPWORDS or chunk_tokens[-1] in STOPWORDS:
                    continue
                
                # Filter B: Length check
                if any(len(t) < 2 for t in chunk_tokens):
                    continue

                chunk_str = " ".join(chunk_tokens)
                raw_ngram_candidates.add(chunk_str)

        if not raw_ngram_candidates and not valid_context_strings:
            return []

        # 2. Process/Filter N-Gram candidates separately
        valid_ngram_tuples = extract_valid_candidates(
            candidate_set=raw_ngram_candidates,
            optimizer=self.optimizer,
            target_category=category,
            top_k=5 # Keep top 5 broad hits
        )
        valid_ngram_strings = [x[0] for x in valid_ngram_tuples]


        # --- PHASE C: Merge and Condense ---
        # We prioritize Context hits, then append N-gram hits.
        # We use a set to deduplicate while preserving priority order.
        
        merged_candidates = []
        seen = set()
        
        # Add Context hits first (Priority)
        if valid_context_strings != []:
            for cand in valid_context_strings:
                if cand not in seen:
                    merged_candidates.append(cand)
                    seen.add(cand)
            return valid_context_strings
        
        # Add N-gram hits second
        if valid_ngram_strings!= []:
            for cand in valid_ngram_strings:
                if cand not in seen:
                    merged_candidates.append(cand)
                    seen.add(cand)

        return condense_candidates(merged_candidates, self.optimizer, category)