import re
from typing import Dict, List, Any, Optional, Set, Tuple
from collections import defaultdict

# Import your ontologies
from src.constants import TissueEnum, TreatmentEnum, MediumEnum
from meta_data_processing.utils.groundingOptimizer import GroundingOptimizer
class UniversalExtractor:
    def __init__(self):
        """
        Args:
            grounding_optimizer: An instance of GroundingOptimizer. 
                                 Required for Phase 3b (Semantic Scanning).
        """
        self.optimizer = GroundingOptimizer()

        # 1. Load Ontologies for Keyword Scanning
        self.known_tissues = {t.value.lower() for t in TissueEnum}
        self.known_treatments = {t.value.lower() for t in TreatmentEnum}
        self.known_mediums = {m.value.lower() for m in MediumEnum}

        # 2. Define High-Priority Fields
        self.priority_fields = [
            'characteristics_ch1', 
            'source_name_ch1', 
            'title', 
            'description', 
            'growth_protocol_ch1',
            'treatment_protocol_ch1'
        ]

    def extract(self, sample_metadata: Dict, study_metadata: Dict) -> Dict:
        result = {
            "tissue": "unknown",
            "treatment": [],
            "medium": "unspecified"
        }

        # --- PHASE 1: Flatten Data ---
        searchable_text = []
        def unpack(obj):
            if isinstance(obj, list): return [str(x) for x in obj]
            if isinstance(obj, str): return [obj]
            return []

        for key in self.priority_fields:
            if key in sample_metadata:
                searchable_text.extend(unpack(sample_metadata[key]))

        if study_metadata:
            for key in ['title', 'summary', 'overall_design']:
                if key in study_metadata:
                    searchable_text.extend(unpack(study_metadata[key]))

        # --- PHASE 2: explicit Key-Value Extraction ---
        kv_extraction = self._extract_key_value_pairs(searchable_text)
        
        if kv_extraction.get('tissue'): result['tissue'] = kv_extraction['tissue']
        if kv_extraction.get('treatment'): result['treatment'].extend(kv_extraction['treatment'])
        if kv_extraction.get('medium'): result['medium'] = kv_extraction['medium']

        # --- PHASE 3: Ontology Scan (Text Match + Semantic) ---
        full_text_blob = " ".join(searchable_text).lower()

        # A. Tissue Scan
        if result['tissue'] == "unknown":
            # 1. Try Direct Text Match
            found_tissue = self._scan_text_match(full_text_blob, self.known_tissues)
            # 2. If fail, Try Vector Match
            if not found_tissue:
                found_tissue = self._scan_semantic_match(full_text_blob, category='tissue')
            
            if found_tissue:
                result['tissue'] = found_tissue

        # B. Medium Scan
        if result['medium'] == "unspecified":
            found_medium = self._scan_text_match(full_text_blob, self.known_mediums)
            if not found_medium:
                found_medium = self._scan_semantic_match(full_text_blob, category='medium')
            
            if found_medium:
                result['medium'] = found_medium

        # C. Treatment Scan (Accumulative)
        # 1. Text Match
        found_treatments = self._scan_text_match_multi(full_text_blob, self.known_treatments)
        
        # 2. Semantic Match (Optional: usually text match is enough for known ontology terms, 
        # but you can enable this if you want to catch synonyms like "water deficit" -> "Drought")
        semantic_treatments = self._scan_semantic_match(full_text_blob, category='treatment')
        if semantic_treatments: found_treatments.append(semantic_treatments)

        current_treatments = set(result['treatment'])
        current_treatments.update(found_treatments)
        
        # --- PHASE 4: Final Cleanup ---
        if not current_treatments:
            current_treatments.add("No stress")
        
        if len(current_treatments) > 1 and "No stress" in current_treatments:
            current_treatments.remove("No stress")

        result['treatment'] = list(current_treatments)

        return result

    def _extract_key_value_pairs(self, text_list: List[str]) -> Dict[str, Any]:
        """Parses strings like 'Tissue: Root'."""
        found = defaultdict(list)
        patterns = {
            'tissue': re.compile(r'(?:tissue|organ|source|cell type)\s*:\s*(.+)', re.IGNORECASE),
            'treatment': re.compile(r'(?:treatment|stress|condition)\s*:\s*(.+)', re.IGNORECASE),
            'medium': re.compile(r'(?:medium|growth|substrate)\s*:\s*(.+)', re.IGNORECASE)
        }

        for text in text_list:
            for key, pattern in patterns.items():
                match = pattern.search(text)
                if match:
                    clean_val = match.group(1).split(';')[0].strip()
                    if len(clean_val) < 50:
                        found[key].append(clean_val)
        
        return {
            'tissue': found['tissue'][0] if found['tissue'] else None,
            'medium': found['medium'][0] if found['medium'] else None,
            'treatment': found['treatment']
        }

    def _is_negated(self, text: str, term: str, window: int = 20) -> bool:
        """
        Checks if a term is preceded by 'no', 'not', 'non-', 'without'.
        """
        # Find index of term
        try:
            idx = text.index(term)
        except ValueError:
            return False

        # Look at the window of characters before the term
        start = max(0, idx - window)
        preceding_text = text[start:idx]

        # Check for negation keywords
        negations = ["no ", "not ", "non-", "without ", "lack of "]
        for neg in negations:
            if neg in preceding_text:
                return True
        return False

    def _scan_text_match(self, text: str, ontology_set: set) -> Optional[str]:
        """
        Scans for exact string presence, handling negations.
        Returns the Longest matching term (prioritizing specificity).
        """
        best_match = None
        longest_len = 0
        
        for term in ontology_set:
            if term in text:
                # Check negation
                if self._is_negated(text, term):
                    continue
                
                # We prioritize longest match (e.g. "Low Light Stress" > "Stress")
                if len(term) > longest_len:
                    longest_len = len(term)
                    best_match = term
        return best_match

    def _scan_text_match_multi(self, text: str, ontology_set: set) -> List[str]:
        """Finds ALL matching keywords (for treatments)."""
        found = []
        for term in ontology_set:
            if term in text:
                if not self._is_negated(text, term):
                    found.append(term)
        return found

    def _scan_semantic_match_0(self, text: str, category: str) -> Optional[str]:
        """
        Splits text into N-grams (1-3 words) and checks semantic similarity against the ontology.
        Returns the label with the highest similarity score found in the text.
        """
        if not self.optimizer:
            return None

        # 1. Clean and Tokenize
        # Remove special chars but keep spaces/hyphens
        clean_text = re.sub(r'[^a-zA-Z0-9\s-]', ' ', text.lower())
        tokens = clean_text.split()
        
        # Limit processing to first N tokens to prevent lag on massive descriptions
        # 500 words is usually enough to cover the biological summary
        tokens = tokens[:500] 

        # 2. Define Stopwords (Fast set for filtering)
        # We don't want to waste compute embedding "and the" or "samples were"
        STOPWORDS = {
            "the", "and", "of", "in", "to", "with", "for", "on", "at", "by", "from", 
            "was", "were", "are", "is", "an", "a", "or", "that", "this", "using", 
            "analysis", "data", "samples", "expression", "profiling", "study", "results"
        }

        # 3. Generate N-gram Candidates (Sliding Window of 1, 2, and 3 words)
        candidates = set()
        max_n = 3  # Most ontology terms are 1-3 words long

        for n in range(1, max_n + 1):
            for i in range(len(tokens) - n + 1):
                chunk_tokens = tokens[i : i + n]
                
                # Filter A: Check if start/end is a stopword (e.g. skip "grown in")
                if chunk_tokens[0] in STOPWORDS or chunk_tokens[-1] in STOPWORDS:
                    continue
                
                # Filter B: Skip if any token is very short/garbage (unless it's a specific chemical abbreviation)
                if any(len(t) < 2 for t in chunk_tokens):
                    continue

                # Reconstruct string
                chunk_str = " ".join(chunk_tokens)
                candidates.add(chunk_str)

        if not candidates:
            return None

        # 4. Score Candidates
        best_label = None
        best_score = 0.0
        extract_valid_candidates(candidates,self.optimizer,category)

        return best_label
    
    def _scan_semantic_match(self, text: str, category: str) -> Optional[str]:
        """
        Splits text into N-grams (1-3 words) and applies Contrastive Semantic Filtering.
        Returns the best matching ontology label.
        """
        

        if not self.optimizer:
            return None

        # 1. Clean and Tokenize
        clean_text = re.sub(r'[^a-zA-Z0-9\s-]', ' ', text.lower())
        tokens = clean_text.split()
        
        # Limit processing
        tokens = tokens[:500] 

        # 2. Define Stopwords
        STOPWORDS = {
            "the", "and", "of", "in", "to", "with", "for", "on", "at", "by", "from", 
            "was", "were", "are", "is", "an", "a", "or", "that", "this", "using", 
            "analysis", "data", "samples", "expression", "profiling", "study", "results",
            "between", "under", "during", "after", "before"
        }

        # 3. Generate N-gram Candidates
        candidates = set()
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
                candidates.add(chunk_str)

        if not candidates:
            return None

        # 4. Score and Filter Candidates (The new call)
        # We ask for the Top 1 match
        top_candidates = extract_valid_candidates(
            candidate_set=candidates,
            optimizer=self.optimizer,
            target_category=category,
            top_k=1
        )

        # 5. Return the label of the best match, or None
        if top_candidates:
            # top_candidates is [(Label, Score)]
            return top_candidates[0][0]
            
        return None


def extract_valid_candidates(candidate_set: Set[str], optimizer, target_category: str, top_k: int = 3) -> List[Tuple[str, float]]:
    """
    Generic Semantic Filter.
    1. Scans n-grams against the Target Category.
    2. Checks against a 'Contrast Category' (Distractor) to prevent cross-contamination.
    3. Returns top K unique ontology labels sorted by confidence score.
    """
    
    # 1. Define Distractors (The Semantic Tug-of-War)
    # If we want X, we make sure it's not Y.
    if target_category == 'treatment':
        contrast_category = 'tissue'
    elif target_category == 'tissue':
        contrast_category = 'treatment' 
    elif target_category == 'medium':
        # Mediums (Soil, Agar) are often confused with Tissues (Root, Seedling)
        contrast_category = 'tissue' 
    else:
        contrast_category = None

    # 2. Hardcoded Stoplist for generic metadata noise
    NOISE_STOPLIST = {
        "wild type", "wild-type", "columbia", "col-0", "mutant", "genotype", 
        "transgenic", "plants", "analysis", "data", "replicate", "study", 
        "experiment", "grown", "growth", "independent", "agb1", "gpa1",
        "control", "mock", "buffer" # Added generic procedural terms
    }

    verified_candidates = {} # Maps Standard Label -> Max Score found

    for term in candidate_set:
        term_lower = term.lower()
        
        # A. Fast Noise Filter
        if any(bad in term_lower for bad in NOISE_STOPLIST):
            continue
        if len(term) < 3: 
            continue

        # B. Get Target Score
        # Assumes optimizer.get_best_match_with_score returns (Label, FloatScore)
        target_label, target_score = optimizer.get_best_match_with_score(term, category=target_category)
        
        # C. Threshold Check (Base Validity)
        # We use 0.82 as a "High Confidence" semantic threshold
        if target_score < 0.82:
            continue

        # D. Contrast Check (The Distractor Trap)
        if contrast_category:
            _, contrast_score = optimizer.get_best_match_with_score(term, category=contrast_category)
            
            # Rule: If the term is MORE similar to the distractor than the target, discard it.
            # E.g. Target=Treatment ("Guard Cells") -> Score 0.40
            #      Contrast=Tissue ("Guard Cells")    -> Score 0.95
            #      Result: Discard.
            if contrast_score > target_score:
                continue

        # E. Aggregation Rule
        # Keep the highest score seen for this specific ontology label
        # (e.g. "NaCl" (0.85) and "Salt" (0.88) both map to "Salinity Stress" -> Keep 0.88)
        if target_label not in verified_candidates:
            verified_candidates[target_label] = target_score
        else:
            verified_candidates[target_label] = max(verified_candidates[target_label], target_score)

    # 3. Sort and Slice
    # Convert to list of tuples: [('Salinity Stress', 0.88), ('Heat Stress', 0.85)]
    sorted_candidates = sorted(verified_candidates.items(), key=lambda x: x[1], reverse=True)
    
    return sorted_candidates[:top_k]