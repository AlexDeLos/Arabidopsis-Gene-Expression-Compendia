import re
from typing import Dict, List, Any, Optional, Set, Tuple
from sentence_transformers import util # Required for cosine similarity
from src.constants import TissueEnum, TreatmentEnum, MediumEnum
from src.meta_data_processing.utils.classes import GroundingOptimizer

# ... (extract_valid_candidates function remains exactly the same) ...
def extract_valid_candidates(candidate_set: Set[str], optimizer: GroundingOptimizer, target_category: str, top_k: int = 3) -> List[Tuple[str, float]]:
    # ... (Keep existing code from previous steps) ...
    # 1. Define Distractors
    if target_category == 'treatment': contrast_category = 'tissue'
    elif target_category == 'tissue': contrast_category = 'treatment' 
    elif target_category == 'medium': contrast_category = 'tissue' 
    else: contrast_category = None

    NOISE_STOPLIST = {
        "wild type", "wild-type", "columbia", "col-0", "mutant", "genotype", 
        "transgenic", "plants", "analysis", "data", "replicate", "study", 
        "experiment", "grown", "growth", "independent", "agb1", "gpa1",
        "control", "mock", "buffer", "treated", "samples", "using", "total", "rna"
    }

    valid_raw_terms = {} 
    for term in candidate_set:
        term_lower = term.lower()
        if any(bad in term_lower for bad in NOISE_STOPLIST): continue
        if len(term) < 3: continue

        target_label, target_score = optimizer.get_best_match_with_score(term, category=target_category)
        if target_score < 0.82: continue

        if contrast_category:
            _, contrast_score = optimizer.get_best_match_with_score(term, category=contrast_category)
            if contrast_score > target_score: continue

        valid_raw_terms[term] = target_score

    sorted_candidates = sorted(valid_raw_terms.items(), key=lambda x: x[1], reverse=True)
    return sorted_candidates[:top_k]


class UniversalExtractor:
    def __init__(self):
        self.optimizer = GroundingOptimizer()
        self.known_tissues = {t.value.lower() for t in TissueEnum}
        self.known_treatments = {t.value.lower() for t in TreatmentEnum}
        self.known_mediums = {m.value.lower() for m in MediumEnum}

        # --- NEW: Pre-compute Explicit Column Vectors ---
        # These are the "Anchor Terms" defining what a column header looks like
        anchors = {
            'tissue': ['tissue', 'organ', 'source', 'cell type', 'tissue type', 'organism part'],
            'treatment': ['treatment', 'stress', 'condition', 'treatment protocol', 'treated with', 'exposure'],
            'medium': ['medium', 'growth medium', 'substrate', 'growth condition', 'culture medium']
        }
        
        # We store the encoded matrices. 
        # Structure: {'tissue': Tensor(Shape: [N_Anchors, 384]), ...}
        self.explicit_vectors = {
            cat: self.optimizer.model.encode(terms) 
            for cat, terms in anchors.items()
        }

    def extract(self, sample_metadata: Dict, study_metadata: Dict) -> Dict:
        """
        1. Explicit Vector Check (Key Similarity).
        2. Priority Field Scan (Title, Source).
        3. Broad Field Scan (Description).
        """
        extracted_data = {
            "tissue": [],
            "treatment": [],
            "medium": []
        }
        
        # --- PROCESS EACH CATEGORY ---
        for category in ['tissue', 'treatment', 'medium']:
            
            # STEP 1: Vector-Based Column Check
            # Checks if any key (e.g. "growth_cond") semantically matches the category anchor
            explicit_hits = self._check_for_explicit_columns(sample_metadata, category)
            
            if explicit_hits:
                extracted_data[category] = explicit_hits
                continue 

            # STEP 2: Priority Text Scan
            # (We pass the category name to help the text gatherer, though exact keywords aren't used for keys anymore)
            priority_text = self._get_priority_text(sample_metadata, category)
            scanned_hits = self._scan_semantic_match(priority_text, category=category)
            
            if scanned_hits:
                extracted_data[category] = list(set(scanned_hits))
                continue

            # STEP 3: Broad Text Scan
            broad_text = self._get_broad_text(sample_metadata, study_metadata)
            scanned_broad = self._scan_semantic_match(broad_text, category=category)
            extracted_data[category] = list(set(scanned_broad))

        return extracted_data

    def _is_semantic_key_match(self, key_text: str, category: str, threshold: float = 0.85) -> bool:
        """
        Encodes the metadata key and compares it to the pre-computed anchor vectors.
        """
        # Clean the key (remove underscores/camelCase to help the semantic model)
        # e.g. "growth_protocol" -> "growth protocol"
        clean_key = key_text.replace('_', ' ').replace('-', ' ')
        clean_key = re.sub(r'(?<!^)(?=[A-Z])', ' ', clean_key).lower().strip()

        # Encode the specific key from the file
        key_vector = self.optimizer.model.encode(clean_key)

        # Calculate Cosine Similarity against the Anchor Matrix for this category
        # result is a list of scores [score_anchor_1, score_anchor_2, ...]
        scores = util.cos_sim(key_vector, self.explicit_vectors[category])[0]

        # If ANY of the anchors match this key with high confidence, return True
        return bool(scores.max() > threshold)

    def _check_for_explicit_columns(self, metadata: Dict, category: str) -> List[str]:
        """
        Looks for values where the KEY is semantically similar to the category.
        """
        found_values = []
        
        # A. Check keys of the metadata dictionary itself
        # e.g. data['growth_medium'] = 'MS Plates'
        for key, val in metadata.items():
            if self._is_semantic_key_match(key, category,threshold=0.85):
                if isinstance(val, list): found_values.extend([str(v) for v in val])
                else: found_values.append(str(val))
        
        # B. Check inside 'characteristics_ch1' list for "Key: Value" patterns
        target_list_keys = ['characteristics_ch1', 'characteristics', 'source_name_ch1']
        
        for list_key in target_list_keys:
            if list_key in metadata:
                val_list = metadata[list_key]
                if not isinstance(val_list, list): val_list = [str(val_list)]
                
                for item in val_list:
                    item_str = str(item)
                    if ':' in item_str:
                        parts = item_str.split(':', 1)
                        key_part = parts[0].strip() # e.g. "Tissue Type"
                        val_part = parts[1].strip() # e.g. "Leaf"
                        
                        # Use vector check on the key part
                        if self._is_semantic_key_match(key_part, category):
                            val_part = val_part.strip('"').strip("'")
                            if len(val_part) > 1:
                                found_values.append(val_part)

        return list(set(found_values))

    def _get_priority_text(self, metadata: Dict, category: str) -> str:
        """Gathers text from Title, Source, and keys matching the category."""
        text_bits = []
        text_bits.append(str(metadata.get('title', '')))
        text_bits.append(str(metadata.get('source_name_ch1', '')))
        
        # Also grab contents of keys that semantically match
        for key, val in metadata.items():
            if self._is_semantic_key_match(key, category):
                if isinstance(val, list): text_bits.extend([str(v) for v in val])
                else: text_bits.append(str(val))
        return " ".join(text_bits)

    def _get_broad_text(self, sample_meta: Dict, study_meta: Dict) -> str:
        # (Same as before)
        text_bits = []
        text_bits.append(str(sample_meta.get('description', '')))
        text_bits.append(str(sample_meta.get('growth_protocol_ch1', '')))
        text_bits.append(str(sample_meta.get('treatment_protocol_ch1', '')))
        text_bits.append(str(study_meta.get('summary', '')))
        text_bits.append(str(study_meta.get('overall_design', '')))
        return " ".join(text_bits)

    def _scan_semantic_match(self, text: str, category: str) -> List[str]:
        # (Same as before)
        if not self.optimizer or not text.strip(): return []
        clean_text = re.sub(r'[^a-zA-Z0-9\s-]', ' ', text.lower())
        tokens = clean_text.split()
        tokens = tokens[:600] 
        STOPWORDS = {"the", "and", "of", "in", "to", "with", "for", "on", "at", "by", "from", "was", "were", "are", "is", "an", "a", "or", "that", "this", "using", "analysis", "data", "samples", "expression", "profiling", "study", "results", "between", "under", "during", "after", "before", "total", "rna"}
        candidates = set()
        max_n = 3 
        for n in range(1, max_n + 1):
            for i in range(len(tokens) - n + 1):
                chunk_tokens = tokens[i : i + n]
                if chunk_tokens[0] in STOPWORDS or chunk_tokens[-1] in STOPWORDS: continue
                if any(len(t) < 2 for t in chunk_tokens): continue
                candidates.add(" ".join(chunk_tokens))
        if not candidates: return []
        top_candidates = extract_valid_candidates(candidates, self.optimizer, category, top_k=5)
        return [item[0] for item in top_candidates]