from typing import List, Optional
import json
import copy

def load_json(path:str):
    with open(path, 'r') as file:
        object = json.load(file)
    return object


class LabelMap:
    def __init__(self, path:Optional[str]=None):
        self.path = path
        if path is None:
            self.map_treatment = {}
            self.map_tissue = {}
            self.map = {}
        else:
            try:
                self.map_treatment = load_json(path+'/map_treatment.json')
                self.map_tissue = load_json(path+'/map_tissue.json')
                self.map = load_json(path+'/map.json')
            except:
                Warning('Path not found')
                self.map_treatment = {}
                self.map_tissue = {}
                self.map = {}

    def add(self,label:str,id)->None:
        self.map[label] = id

    def add_treatment(self,label:str,id)->None:
        self.map_treatment[label] = id

    def add_tissue(self,label:str,id)->None:
        self.map_tissue[label] = id
    
    def add_mapping(self,og,grounded)->None:
        for el in og:
            if  isinstance(og[el],List):
                for i,term in enumerate(og[el]):
                    if len(og[el]) != len(grounded[el]):
                        self.add_treatment(str(og[el]),grounded[el])
                        break
                    else:
                        self.add_treatment(term,grounded[el][i])
                        
            elif (el =='tissue'):
                self.add_tissue(og[el],grounded[el])
            else:
                self.add(og[el],grounded[el])

    def in_maps_evaluated(self,el)->bool:
        return el in self.map_treatment or el in self.map_tissue

    def in_maps(self,el,key)->bool:
        if key == 'tissue':
            return el in self.map_tissue
        elif key == 'treatment':
            return el in self.map_treatment
        elif key == 'medium':
            return el in self.map
        else:
            raise ValueError(f"unkown key {key}")
    

    def need_new_mapping(self,sample:dict)->bool:
        # will return True is it is in the maps
        run = False
        for key in sample:
            if  isinstance(sample[key],List):
                for term in sample[key]:
                    if not self.in_maps(term,key):
                        if self.in_maps(str(sample[key]),key):
                            pass
                        else:
                            run = True
            elif key !='id':
                if not self.in_maps(sample[key],key):
                    run = True
        return run
    
    def save_map(self) -> None:
        with open(f'{self.path}/map_treatment.json', 'w') as handle:
            json.dump(self.map_treatment, handle)
        with open(f'{self.path}/map_tissue.json', 'w') as handle:
            json.dump(self.map_tissue, handle)
        with open(f'{self.path}/map.json', 'w') as handle:
            json.dump(self.map, handle)
    
    
    def apply_mappings(self,og:dict)-> dict:
        ret = copy.deepcopy(og)
        for el in og:
            if  isinstance(og[el],List):
                for i,_ in enumerate(og[el]):
                    try:
                        ret[el][i] = self.map_treatment[og[el][i]]
                    except:
                        ret[el] =self.map_treatment[str(og[el])]
                        break
   
            elif el =='tissue':
                    ret[el] = self.map_tissue[og[el]]
            else:
                    ret[el] = self.map[og[el]]

        return ret

from sentence_transformers import SentenceTransformer
from sklearn.metrics.pairwise import cosine_similarity
import json
from typing import List, Dict

class GroundingOptimizer:
    def __init__(self, ontology_path='data/ontology.json'):
        # Load a tiny, fast model (runs locally)
        self.model = SentenceTransformer('all-MiniLM-L6-v2')
        
        # Define your standard labels (The Target Ontology)
        self.valid_treatments = [
            "Drought Stress",
            "Dehidration Stress",
            "Salinity Stress",
            "Heat Stress",
            "Cold Stress",
            "Chemical Stress",
            "Nutrient Deficiency",
            "Biotic Stress",
            "Low Light Stress",
            "High Light Stress",
            "Red Light Stress",
            "Other Light Stress",
            "Other stress",
            "No stress"]
        self.valid_tissues = ["root", "leaf", "flower", "shoot", "rosette", "bud", "whole_plant", "silique", "callus", "seed", "seedling", "unknown"]
        # Pre-compute vectors for your ontology
        self.treatment_ontology_vectors = self.model.encode(self.valid_treatments)
        self.tissue_ontology_vectors = self.model.encode(self.valid_tissues)
        
    def find_semantic_match(self, raw_term: str, threshold=0.85) -> Optional[str]:
        """
        Checks if the raw term is semantically identical to a known label
        without using an LLM.
        """
        if not raw_term:
            return None
            
        # Encode the raw term
        term_vector = self.model.encode([raw_term])
        
        # Calculate similarity against all valid ontology terms
        scores = cosine_similarity(term_vector, self.treatment_ontology_vectors)[0]
        
        # Find best match
        best_idx = scores.argmax()
        if scores[best_idx] > threshold:
            return self.valid_treatments[best_idx]
        
        return None

    def batch_process_study(self, data: Dict, extracted_samples: List[Dict], llm_func, label_map_obj:LabelMap):
        """
        TODO: this needs to use the mappings from 'label_map_obj', if not use the find semantic match function, and if not then do the LLM call
        Args:
            label_map_obj: The full 'seen' LabelMap object (contains .map_treatment, .map_tissue, .map)
        """
        
        # We will build a single local cache for this study
        # Structure: { "50mM NaCl": "Salinity Stress", "root": "root", "soil": "soil" }
        mapping_cache_tissue = {}
        mapping_cache_treatment = {}
        mapping_cache_medium = {}
        
        # Set to collect unknowns that need the LLM (Currently only Treatments supported by your LLM prompt)
        unknown_treatments_for_llm = set()
        unknown_tissues_for_llm = set()

        # --- HELPER: Logic to check Map -> then Vector ---
        def resolve_term(term, specific_map):
            # 1. Check Memory (The specific map for this field)
            if specific_map and term in specific_map:
                return specific_map[term]
            
            # 2. Check Vector Match
            match = self.find_semantic_match(term)
            if match:
                return match
            
            return None

        # --- STEP 1: Process Fields ---
        
        for sample in extracted_samples:
            
            # A. Process TREATMENT (List of Strings) -> map_treatment
            treatments = sample.get('treatment', [])
            for t in treatments:
                if t in mapping_cache_treatment: continue # Already processed this string
                
                resolved = resolve_term(t, label_map_obj.map_treatment)
                if resolved:
                    mapping_cache_treatment[t] = resolved
                else:
                    unknown_treatments_for_llm.add(t)

            # B. Process TISSUE (String) -> map_tissue
            tissue = sample.get('tissue')
            if tissue in mapping_cache_tissue: continue # Already processed this string
            
            resolved_tis = resolve_term(tissue, label_map_obj.map_tissue)
            if resolved_tis:
                mapping_cache_tissue[tissue] = resolved_tis
            else:
                unknown_tissues_for_llm.add(tissue)
            # C. Process MEDIUM (String) -> map (generic)
            medium = sample.get('medium')
            if medium and medium not in mapping_cache_medium:
                resolved_m = resolve_term(medium, label_map_obj.map)
                if resolved_m:
                    mapping_cache_medium[medium] = resolved_m
                else:
                    mapping_cache_medium[medium] = medium
        # --- STEP 2: Batch LLM Call (Treatments Only) ---
        data_study = data['study_metadata']
        if unknown_treatments_for_llm:
            print(f"Querying LLM for {len(unknown_treatments_for_llm)} unique treatments...")
            try:
                llm_result = llm_func(list(unknown_treatments_for_llm),data_study)
                mapping_cache_treatment.update(llm_result)
            except Exception as e:
                print(f"LLM Batch Error: {e}")
        
        if unknown_tissues_for_llm:
            print(f"Querying LLM for {len(unknown_tissues_for_llm)} unique tissues...")
            try:
                llm_result = llm_func(list(unknown_tissues_for_llm),data_study)
                mapping_cache_tissue.update(llm_result)
            except Exception as e:
                print(f"LLM Batch Error: {e}")

        # --- STEP 3: Apply mappings back to samples ---
        final_samples = []
        for sample in extracted_samples:
            new_sample = sample.copy()
            
            # 1. Update Treatments
            new_treatments = []
            for t in sample.get('treatment', []):
                # Fallback to "Other stress" only if LLM failed and it's a treatment
                val = mapping_cache_treatment.get(t)
                if not val and t in unknown_treatments_for_llm:
                    val = "Other stress"
                new_treatments.append(val if val else t)
            new_sample['treatment'] = list(set(new_treatments))
            
            # 2. Update Tissue
            if sample.get('tissue') in mapping_cache_tissue:
                new_sample['tissue'] = mapping_cache_tissue[sample['tissue']]
                
            # 3. Update Medium
            if sample.get('medium') in mapping_cache_medium:
                new_sample['medium'] = mapping_cache_medium[sample['medium']]

            final_samples.append(new_sample)
            
        return final_samples