import re

def GSE72949_extractor(sample_metadata: dict) -> dict:
    extracted_data = {
        "tissue": "unknown",
        "treatment": [],
        "medium": "unspecified"
    }

    # Helper to safely get string value from metadata (can be list or single string)
    def _get_value(key, default=""):
        value = sample_metadata.get(key)
        if isinstance(value, list):
            return " ".join(value)
        elif isinstance(value, str):
            return value
        return default

    # Combine relevant text fields for keyword searching
    characteristics_ch1_str = _get_value('characteristics_ch1')
    source_name_ch1_str = _get_value('source_name_ch1')
    title_str = _get_value('title')
    description_str = _get_value('description')
    treatment_protocol_ch1_str = _get_value('treatment_protocol_ch1')
    growth_protocol_ch1_str = _get_value('growth_protocol_ch1')

    # --- Extract Tissue ---
    tissue_found = False
    for char_item in sample_metadata.get('characteristics_ch1', []):
        if 'tissue:' in char_item:
            tissue_val = char_item.split('tissue:', 1)[1].strip().lower()
            if 'seedling' in tissue_val:
                extracted_data['tissue'] = 'seedling'
                tissue_found = True
                break
            elif 'root' in tissue_val:
                extracted_data['tissue'] = 'root'
                tissue_found = True
                break
            elif 'leaf' in tissue_val:
                extracted_data['tissue'] = 'leaf'
                tissue_found = True
                break
            elif 'flower' in tissue_val:
                extracted_data['tissue'] = 'flower'
                tissue_found = True
                break
            elif 'shoot' in tissue_val:
                extracted_data['tissue'] = 'shoot'
                tissue_found = True
                break
            elif 'rosette' in tissue_val:
                extracted_data['tissue'] = 'rosette'
                tissue_found = True
                break
            elif 'bud' in tissue_val:
                extracted_data['tissue'] = 'bud'
                tissue_found = True
                break
            elif 'silique' in tissue_val:
                extracted_data['tissue'] = 'silique'
                tissue_found = True
                break
            elif 'callus' in tissue_val:
                extracted_data['tissue'] = 'callus'
                tissue_found = True
                break
            elif 'seed' in tissue_val:
                extracted_data['tissue'] = 'seed'
                tissue_found = True
                break
            elif 'whole plant' in tissue_val or 'whole-plant' in tissue_val:
                extracted_data['tissue'] = 'whole_plant'
                tissue_found = True
                break

    if not tissue_found:
        combined_text = (source_name_ch1_str + " " + title_str + " " + description_str).lower()
        if 'seedling' in combined_text:
            extracted_data['tissue'] = 'seedling'
        elif 'whole plant' in combined_text or 'whole-plant' in combined_text:
            extracted_data['tissue'] = 'whole_plant'
        elif 'plant' in combined_text and 'seedling' not in combined_text:
            extracted_data['tissue'] = 'whole_plant'

    # --- Extract Treatment ---
    treatments = set()

    # Check characteristics_ch1 for 'stress:'
    for char_item in sample_metadata.get('characteristics_ch1', []):
        if 'stress:' in char_item:
            stress_val = char_item.split('stress:', 1)[1].strip().lower()
            if 'no stress' in stress_val or 'control' in stress_val:
                treatments.add('No stress')
            if 'drought' in stress_val: treatments.add('Drought Stress')
            if 'salinity' in stress_val: treatments.add('Salinity Stress')
            if 'heat' in stress_val: treatments.add('Heat Stress')
            if 'cold' in stress_val: treatments.add('Cold Stress')
            if 'chemical' in stress_val: treatments.add('Chemical Stress')
            if 'nutrient deficiency' in stress_val: treatments.add('Nutrient Deficiency')
            if 'pathogen' in stress_val: treatments.add('Pathogen Attack')
            if 'low light' in stress_val: treatments.add('Low Light Stress')
            if 'high light' in stress_val: treatments.add('High Light Stress')
            if 'red light' in stress_val: treatments.add('Red Light Stress')
            if 'stress' in stress_val and not any(t in stress_val for t in ['drought', 'salinity', 'heat', 'cold', 'chemical', 'nutrient deficiency', 'pathogen', 'light', 'no stress', 'control']):
                treatments.add('Other stress')

    # Check treatment_protocol_ch1
    protocol_text = treatment_protocol_ch1_str.lower()
    if 'heat regime' in protocol_text or 'hs treatment' in protocol_text or 'heat stress' in protocol_text:
        treatments.add('Heat Stress')
    if 'drought' in protocol_text: treatments.add('Drought Stress')
    if 'salinity' in protocol_text: treatments.add('Salinity Stress')
    if 'cold' in protocol_text: treatments.add('Cold Stress')
    if 'chemical' in protocol_text: treatments.add('Chemical Stress')
    if 'nutrient deficiency' in protocol_text: treatments.add('Nutrient Deficiency')
    if 'pathogen' in protocol_text: treatments.add('Pathogen Attack')
    if 'low light' in protocol_text: treatments.add('Low Light Stress')
    if 'high light' in protocol_text: treatments.add('High Light Stress')
    if 'red light' in protocol_text: treatments.add('Red Light Stress')
    if 'stress' in protocol_text and not any(t in protocol_text for t in ['drought', 'salinity', 'heat', 'cold', 'chemical', 'nutrient deficiency', 'pathogen', 'light']):
        treatments.add('Other stress')

    # Check title and description for stress keywords
    combined_title_desc = (title_str + " " + description_str).lower()
    if 'drought' in combined_title_desc: treatments.add('Drought Stress')
    if 'salinity' in combined_title_desc: treatments.add('Salinity Stress')
    if 'heat' in combined_title_desc: treatments.add('Heat Stress')
    if 'cold' in combined_title_desc: treatments.add('Cold Stress')
    if 'chemical' in combined_title_desc: treatments.add('Chemical Stress')
    if 'nutrient deficiency' in combined_title_desc: treatments.add('Nutrient Deficiency')
    if 'pathogen' in combined_title_desc: treatments.add('Pathogen Attack')
    if 'low light' in combined_title_desc: treatments.add('Low Light Stress')
    if 'high light' in combined_title_desc: treatments.add('High Light Stress')
    if 'red light' in combined_title_desc: treatments.add('Red Light Stress')
    if 'stress' in combined_title_desc and not any(t in combined_title_desc for t in ['drought', 'salinity', 'heat', 'cold', 'chemical', 'nutrient deficiency', 'pathogen', 'light']):
        treatments.add('Other stress')

    # If 'No stress' is present, and other stresses are also present, remove 'No stress'
    if 'No stress' in treatments and len(treatments) > 1:
        treatments.remove('No stress')
    elif not treatments: # If no specific stress found, default to 'No stress'
        treatments.add('No stress')

    extracted_data['treatment'] = sorted(list(treatments))

    # --- Extract Medium ---
    import re
    medium_text = growth_protocol_ch1_str.lower()
    
    if 'agar medium' in medium_text or 'agar plates' in medium_text or 'agar-solidified' in medium_text:
        content = ''
        # Try to extract details about the agar medium from parentheses first
        match = re.search(r'(?:agar medium|agar plates|agar-solidified medium)\s*\(?([^\)]*?)\)?', medium_text)
        if match and match.group(1):
            content = match.group(1)
        
        # If no content from parentheses, look for content immediately following "agar medium" or similar
        if not content:
            match = re.search(r'(?:agar medium|agar plates|agar-solidified medium)\s*(?:with|containing|supplemented with)?\s*([^.,;]+)', medium_text)
            if match and match.group(1):
                content = match.group(1)
        
        if content:
            # Clean up content
            content = content.replace('supplemented with', 'with').strip()
            content = content.replace('[w/v]', '').strip()
            content = content.replace('an equal volume of', '').strip()
            content = content.replace('containing', '').strip()
            content = content.replace('petri dishes', '').strip()
            content = content.replace('in', '').strip()
            content = content.replace('equal volume of', '').strip()
            content = content.replace('  ', ' ').strip()
            
            # Specific handling for MS medium to keep capitalization
            if 'ms medium' in content:
                content = content.replace('ms medium', 'MS medium')
            elif 'ms' in content and 'medium' in content:
                content = content.replace('ms', 'MS')

            if content:
                extracted_data['medium'] = content[0].upper() + content[1:] if content[0].islower() else content
            else:
                extracted_data['medium'] = 'Agar medium'
        else:
            extracted_data['medium'] = 'Agar medium'

    elif 'soil' in medium_text:
        extracted_data['medium'] = 'Soil'
    elif 'hydroponic' in medium_text:
        extracted_data['medium'] = 'Hydroponic solution'
    elif 'water' in medium_text:
        extracted_data['medium'] = 'Water'
    
    # Fallback/inference for medium if still unspecified
    if extracted_data['medium'] == 'unspecified':
        if 'ms medium' in medium_text:
            extracted_data['medium'] = 'MS medium'
        elif 'sucrose' in medium_text:
            extracted_data['medium'] = 'Medium with sucrose'
        elif extracted_data['tissue'] == 'whole_plant':
            extracted_data['medium'] = 'Soil'
        elif extracted_data['tissue'] == 'seedling':
            extracted_data['medium'] = 'Agar medium'

    return extracted_data

def GSE71237_extractor(sample_metadata: dict) -> dict:
    extracted_data = {}

    # --- Extract Tissue ---
    tissue_found = "unknown"
    if 'characteristics_ch1' in sample_metadata:
        for char_str in sample_metadata['characteristics_ch1']:
            if char_str.lower().startswith('tissue:'):
                tissue_val = char_str.split(':', 1)[1].strip().lower()
                if 'seedlings' in tissue_val:
                    tissue_found = "seedling"
                elif 'whole plant' in tissue_val:
                    tissue_found = "whole_plant"
                elif 'callus' in tissue_val:
                    tissue_found = "callus"
                elif 'roots' in tissue_val or 'root' in tissue_val:
                    tissue_found = "root"
                elif 'leaves' in tissue_val or 'leaf' in tissue_val:
                    tissue_found = "leaf"
                elif 'flowers' in tissue_val or 'flower' in tissue_val:
                    tissue_found = "flower"
                elif 'shoots' in tissue_val or 'shoot' in tissue_val:
                    tissue_found = "shoot"
                elif 'rosettes' in tissue_val or 'rosette' in tissue_val:
                    tissue_found = "rosette"
                elif 'buds' in tissue_val or 'bud' in tissue_val:
                    tissue_found = "bud"
                elif 'siliques' in tissue_val or 'silique' in tissue_val:
                    tissue_found = "silique"
                elif 'seeds' in tissue_val or 'seed' in tissue_val:
                    tissue_found = "seed"
                break # Found tissue, no need to check other characteristics

    extracted_data['tissue'] = tissue_found

    # --- Extract Treatment ---
    treatments = set()
    text_sources = []
    if 'title' in sample_metadata:
        text_sources.extend(sample_metadata['title'])
    if 'source_name_ch1' in sample_metadata:
        text_sources.extend(sample_metadata['source_name_ch1'])
    if 'description' in sample_metadata:
        text_sources.extend(sample_metadata['description'])
    if 'treatment_protocol_ch1' in sample_metadata:
        text_sources.extend(sample_metadata['treatment_protocol_ch1'])

    combined_text = " ".join(text_sources).lower()

    # Define keywords for treatments
    treatment_keywords = {
        "Drought Stress": ["drought", "water potential", "peg", "low water"],
        "Salinity Stress": ["salinity", "nacl", "salt"],
        "Heat Stress": ["heat", "high temperature", "hot"],
        "Cold Stress": ["cold", "low temperature"],
        "Chemical Stress": ["chemical", "herbicide", "pesticide", "metal", "toxic", "cadmium", "aluminum", "heavy metal"],
        "Nutrient Deficiency": ["nutrient deficiency", "low nitrogen", "low phosphate", "starvation", "nutrient deprivation"],
        "Pathogen Attack": ["pathogen", "infection", "bacteria", "virus", "fungus", "elicitor"],
        "Low Light Stress": ["low light", "darkness", "shade"],
        "High Light Stress": ["high light", "uv-b"],
        "Red Light Stress": ["red light"],
        "Other Light Stress": ["light stress"], # General light stress if not more specific
        "Other stress": ["stress"], # General stress if no specific type found
        "No stress": ["control", "no stress", "untreated", "mock"]
    }

    found_specific_stress = False
    for treatment_type, keywords in treatment_keywords.items():
        if treatment_type == "No stress" or treatment_type == "Other stress" or treatment_type == "Other Light Stress":
            continue # Handle these later to prioritize specific stresses

        for keyword in keywords:
            if keyword in combined_text:
                treatments.add(treatment_type)
                found_specific_stress = True
                break # Move to next treatment type once one keyword is found

    # Handle general light stress if no specific light stress was found
    if not any(t in treatments for t in ["Low Light Stress", "High Light Stress", "Red Light Stress"]):
        for keyword in treatment_keywords["Other Light Stress"]:
            if keyword in combined_text:
                treatments.add("Other Light Stress")
                found_specific_stress = True
                break

    # Handle general "stress" if no specific one was found
    if not found_specific_stress:
        for keyword in treatment_keywords["Other stress"]:
            if keyword in combined_text and not any(k in combined_text for k in treatment_keywords["No stress"]):
                treatments.add("Other stress")
                found_specific_stress = True
                break

    # If no specific stress or general stress was found, check for "No stress" keywords
    if not found_specific_stress:
        for keyword in treatment_keywords["No stress"]:
            if keyword in combined_text:
                treatments.add("No stress")
                break
    
    # If nothing at all was found, default to "No stress"
    if not treatments:
        treatments.add("No stress")

    # Ensure "No stress" is only present if no other stress is found
    if len(treatments) > 1 and "No stress" in treatments:
        treatments.remove("No stress")

    extracted_data['treatment'] = sorted(list(treatments)) # Sort for consistent output

    # --- Extract Medium ---
    medium_found = "unspecified"
    medium_sources = []
    if 'growth_protocol_ch1' in sample_metadata:
        medium_sources.extend(sample_metadata['growth_protocol_ch1'])
    if 'treatment_protocol_ch1' in sample_metadata:
        medium_sources.extend(sample_metadata['treatment_protocol_ch1'])

    combined_medium_text = " ".join(medium_sources).lower()

    if "agar" in combined_medium_text:
        if "ms salts" in combined_medium_text or "ms media" in combined_medium_text:
            medium_found = "agar with MS salts"
        elif "control media" in combined_medium_text:
            medium_found = "agar with control media"
        else:
            medium_found = "agar"
    elif "soil" in combined_medium_text:
        medium_found = "soil"
    elif "hydroponic" in combined_medium_text:
        medium_found = "hydroponic solution"
    elif "liquid culture" in combined_medium_text:
        medium_found = "liquid culture"
    elif "vermiculite" in combined_medium_text:
        medium_found = "vermiculite"
    elif "perlite" in combined_medium_text:
        medium_found = "perlite"
    elif "sand" in combined_medium_text:
        medium_found = "sand"
    elif "water" in combined_medium_text:
        medium_found = "water"
    elif "ms salts" in combined_medium_text or "ms media" in combined_medium_text:
        medium_found = "MS salts medium"
    elif "control media" in combined_medium_text:
        medium_found = "control media"
    
    extracted_data['medium'] = medium_found

    return extracted_data

def GSE40061_extractor(sample_metadata: dict) -> dict:
    """
    Extracts tissue, treatment, and medium information from a sample_metadata dictionary
    following a predefined schema.
    """

    def _get_text(metadata_dict, keys):
        """Helper to concatenate and lowercase text from specified metadata keys."""
        text_parts = []
        for key in keys:
            if key in metadata_dict and isinstance(metadata_dict[key], list):
                text_parts.extend(metadata_dict[key])
        return " ".join(text_parts).lower()

    result = {}

    # --- Extract Tissue ---
    tissue_keywords = {
        "root": ["root"],
        "leaf": ["leaf"],
        "flower": ["flower"],
        "shoot": ["shoot"],
        "rosette": ["rosette"],
        "bud": ["bud"],
        "whole_plant": ["whole plant", "whole-plant", "plant"],
        "silique": ["silique"],
        "callus": ["callus"],
        "seed": ["seed"],
        "seedling": ["seedling"]
    }
    extracted_tissue = "unknown"
    
    # Prioritize 'tissue:' in characteristics_ch1
    for char in sample_metadata.get('characteristics_ch1', []):
        if isinstance(char, str) and char.lower().startswith('tissue:'):
            specific_tissue_val = char.split(':', 1)[1].strip().lower()
            for tissue_enum, keywords in tissue_keywords.items():
                # Check for exact match or if any keyword is a substring of the specific tissue value
                if specific_tissue_val in keywords or any(kw in specific_tissue_val for kw in keywords):
                    extracted_tissue = tissue_enum
                    break
            if extracted_tissue != "unknown":
                break

    # If not found in specific 'tissue:' field, search general text
    if extracted_tissue == "unknown":
        search_text_tissue = _get_text(sample_metadata, ['characteristics_ch1', 'source_name_ch1', 'title', 'description'])
        for tissue_enum, keywords in tissue_keywords.items():
            if any(kw in search_text_tissue for kw in keywords):
                extracted_tissue = tissue_enum
                break
    
    result['tissue'] = extracted_tissue

    # --- Extract Treatment ---
    treatment_keywords = {
        "Drought Stress": ["drought", "water deficit", "water stress", "dehydration"],
        "Salinity Stress": ["salinity", "salt stress", "nacl"],
        "Heat Stress": ["heat stress", "high temperature", "heat shock"],
        "Cold Stress": ["cold stress", "low temperature", "chilling"],
        "Chemical Stress": ["chemical stress", "herbicide", "pesticide", "heavy metal", "cadmium", "aluminum", "toxic", "hormone treatment", "aba treatment", "methyl jasmonate", "meja", "auxin", "cytokinin"],
        "Nutrient Deficiency": ["nutrient deficiency", "starvation", "nitrogen deprivation", "phosphate deprivation", "low nitrogen", "low phosphate"],
        "Pathogen Attack": ["pathogen", "infection", "fungal", "bacterial", "viral", "disease", "biotic stress", "pseudomonas", "botrytis", "powdery mildew"],
        "Low Light Stress": ["low light", "shade", "darkness"],
        "High Light Stress": ["high light", "excess light", "uv stress"],
        "Red Light Stress": ["red light"],
        "Other Light Stress": ["light stress"],
        "Other stress": ["stress", "abiotic stress", "mechanical stress", "wounding", "oxidative stress"],
        "No stress": ["no stress", "control", "untreated", "mock"]
    }
    found_treatments = set()
    search_text_treatment = _get_text(sample_metadata, ['treatment_protocol_ch1', 'characteristics_ch1', 'source_name_ch1', 'title', 'description'])

    # First, find all specific stresses (excluding "No stress", "Other stress", "Other Light Stress")
    for treatment_enum, keywords in treatment_keywords.items():
        if treatment_enum in ["No stress", "Other stress", "Other Light Stress"]:
            continue
        if any(kw in search_text_treatment for kw in keywords):
            found_treatments.add(treatment_enum)

    # Handle "Other Light Stress" if no more specific light stress was found
    if not any(t in found_treatments for t in ["Low Light Stress", "High Light Stress", "Red Light Stress"]):
        if any(kw in search_text_treatment for kw in treatment_keywords["Other Light Stress"]):
            found_treatments.add("Other Light Stress")

    # Handle "Other stress" if no specific stress (including specific light) was found
    if not found_treatments:
        if any(kw in search_text_treatment for kw in treatment_keywords["Other stress"]):
            found_treatments.add("Other stress")

    # Final determination for "No stress"
    if not found_treatments: # If no specific, specific light, or general "Other stress" was found
        if any(kw in search_text_treatment for kw in treatment_keywords["No stress"]):
            found_treatments.add("No stress")
        else:
            found_treatments.add("No stress") # Default if absolutely nothing is found

    result['treatment'] = sorted(list(found_treatments))

    # --- Extract Medium ---
    medium_keywords = {
        "soil": ["soil", "potting mix", "vermiculite", "sand"],
        "agar": ["agar", "agar plate", "gel media"],
        "MS medium": ["ms medium", "murashige skoog", "ms media"],
        "liquid culture": ["liquid culture", "hydroponics", "liquid medium", "hoagland solution"],
    }
    extracted_medium = "unspecified"
    search_text_medium = _get_text(sample_metadata, ['growth_protocol_ch1', 'characteristics_ch1', 'source_name_ch1', 'title', 'description'])

    for medium_enum, keywords in medium_keywords.items():
        if any(kw in search_text_medium for kw in keywords):
            extracted_medium = medium_enum
            break
    
    result['medium'] = extracted_medium

    return result

def GSE6583_extractor(sample_metadata: dict) -> dict:
    """
    Extracts tissue, treatment, and medium information from the sample_metadata
    dictionary according to a predefined schema.
    """

    def _get_text_from_metadata(metadata: dict, keys: list) -> str:
        """
        Extracts and concatenates text from specified keys in the metadata,
        handling list values by joining them.
        """
        text_parts = []
        for key in keys:
            if key in metadata:
                value = metadata[key]
                if isinstance(value, list):
                    text_parts.extend([str(item) for item in value if item is not None])
                else:
                    text_parts.append(str(value))
        return " ".join(text_parts).lower()

    extracted_data = {}

    # Combine relevant fields into a single searchable string
    search_fields = [
        'title', 'source_name_ch1', 'characteristics_ch1',
        'description', 'extract_protocol_ch1', 'hyb_protocol', 'scan_protocol'
    ]
    combined_text = _get_text_from_metadata(sample_metadata, search_fields)

    # --- Extract Tissue ---
    tissue_found = "unknown"
    tissue_keywords = {
        "root": "root", "roots": "root",
        "leaf": "leaf", "leaves": "leaf",
        "flower": "flower", "flowers": "flower",
        "shoot": "shoot", "shoots": "shoot",
        "rosette": "rosette",
        "bud": "bud", "buds": "bud",
        "whole plant": "whole_plant", "whole-plant": "whole_plant",
        "silique": "silique", "siliques": "silique",
        "callus": "callus",
        "seed": "seed", "seeds": "seed",
        "seedling": "seedling", "seedlings": "seedling"
    }

    for keyword, schema_value in tissue_keywords.items():
        if keyword in combined_text:
            tissue_found = schema_value
            break
    
    # Special inference for seedling/whole_plant if age is mentioned and no specific tissue
    if tissue_found == "unknown":
        age_mentioned = any(term in combined_text for term in ["age:", "week", "day", "old"])
        
        if age_mentioned:
            if "seedling" in combined_text or "seed" in combined_text:
                tissue_found = "seedling"
            elif "plant" in combined_text:
                tissue_found = "whole_plant"
            else:
                # "3 weeks" is typically a seedling stage for common model plants like Arabidopsis
                tissue_found = "seedling"
        
    extracted_data["tissue"] = tissue_found

    # --- Extract Treatment ---
    treatments = set()
    treatment_keywords = {
        "drought": "Drought Stress",
        "salinity": "Salinity Stress", "salt": "Salinity Stress",
        "heat": "Heat Stress",
        "cold": "Cold Stress",
        "chemical": "Chemical Stress",
        "nutrient deficiency": "Nutrient Deficiency", "low nitrogen": "Nutrient Deficiency", "low phosphate": "Nutrient Deficiency",
        "pathogen": "Pathogen Attack", "infection": "Pathogen Attack", "bacteria": "Pathogen Attack", "fungus": "Pathogen Attack",
        "low light": "Low Light Stress",
        "high light": "High Light Stress",
        "red light": "Red Light Stress",
    }

    found_specific_treatment = False
    for keyword, schema_value in treatment_keywords.items():
        if keyword in combined_text:
            treatments.add(schema_value)
            found_specific_treatment = True
    
    # Check for "Other Light Stress"
    if "light" in combined_text and not any(s in combined_text for s in ["low light", "high light", "red light"]):
        treatments.add("Other Light Stress")
        found_specific_treatment = True

    # Check for "Other stress" if "stress" is mentioned but no specific type caught
    if "stress" in combined_text and not found_specific_treatment:
        treatments.add("Other stress")
    
    if not treatments:
        treatments.add("No stress")
    
    extracted_data["treatment"] = sorted(list(treatments))

    # --- Extract Medium ---
    medium_found = "unspecified"
    medium_keywords = {
        "soil": "soil",
        "agar": "agar",
        "liquid": "liquid", "liquid medium": "liquid",
        "hydroponic": "hydroponic",
        "ms medium": "MS medium", "murashige skoog": "MS medium",
        "vermiculite": "vermiculite",
        "petri dish": "agar"
    }

    for keyword, schema_value in medium_keywords.items():
        if keyword in combined_text:
            medium_found = schema_value
            break
    
    # Infer "soil" if tissue is whole_plant/seedling and no medium found
    if medium_found == "unspecified" and extracted_data["tissue"] in ["whole_plant", "seedling"]:
        medium_found = "soil"
    
    extracted_data["medium"] = medium_found

    return extracted_data

def GSE19700_extractor(sample_metadata: dict) -> dict:
    extracted_data = {
        "tissue": "unknown",
        "treatment": [],
        "medium": "unspecified"
    }

    # Combine relevant text fields for easier searching, handling potential missing keys or non-string values in lists
    all_text_parts = []
    for key in ['title', 'source_name_ch1', 'characteristics_ch1', 'treatment_protocol_ch1', 'growth_protocol_ch1', 'description']:
        if key in sample_metadata:
            value = sample_metadata[key]
            if isinstance(value, list):
                all_text_parts.extend([item.lower() for item in value if isinstance(item, str)])
            elif isinstance(value, str):
                all_text_parts.append(value.lower())
    
    combined_text = " ".join(all_text_parts)

    # --- Extract Tissue ---
    tissue_mapping = {
        "whole rosette": "rosette",
        "rosette": "rosette",
        "leaf": "leaf",
        "root": "root",
        "flower": "flower",
        "shoot": "shoot",
        "bud": "bud",
        "whole plant": "whole_plant",
        "silique": "silique",
        "callus": "callus",
        "seedling": "seedling",
        "seed": "seed"
    }
    
    # Prioritize 'characteristics_ch1' for explicit "tissue:" tag
    if 'characteristics_ch1' in sample_metadata and isinstance(sample_metadata['characteristics_ch1'], list):
        for char in sample_metadata['characteristics_ch1']:
            char_lower = char.lower()
            if 'tissue:' in char_lower:
                extracted_tissue_str = char_lower.split('tissue:')[1].strip()
                for keyword, mapped_tissue in tissue_mapping.items():
                    if keyword in extracted_tissue_str: # Use 'in' for partial matches like "whole rosette"
                        extracted_data["tissue"] = mapped_tissue
                        break
                if extracted_data["tissue"] != "unknown":
                    break # Found a specific tissue from characteristics, stop searching
    
    # If not found or still 'unknown', search in combined text
    if extracted_data["tissue"] == "unknown":
        for keyword, mapped_tissue in tissue_mapping.items():
            if keyword in combined_text:
                extracted_data["tissue"] = mapped_tissue
                break

    # --- Extract Treatment ---
    found_treatments = set()
    treatment_keywords = {
        "drought stress": ["water-limited", "drought", "no water", "water deficit", "water deprivation", "dry", "desiccation"],
        "salinity stress": ["salinity", "salt stress", "nacl"],
        "heat stress": ["heat stress", "high temperature", "heat shock"],
        "cold stress": ["cold stress", "low temperature", "chilling"],
        "chemical stress": ["chemical stress", "herbicide", "pesticide", "heavy metal", "cadmium", "aluminum", "toxic"],
        "nutrient deficiency": ["nutrient deficiency", "low nitrogen", "phosphate starvation", "nitrogen starvation", "nutrient starvation"],
        "pathogen attack": ["pathogen", "infection", "fungus", "bacteria", "virus", "disease"],
        "low light stress": ["low light", "shade", "darkness"],
        "high light stress": ["high light"],
        "red light stress": ["red light"],
        "no stress": ["control", "no stress", "well-watered", "normal conditions", "untreated", "mock"]
    }

    # Check for explicit "treatment:" in characteristics_ch1 first
    if 'characteristics_ch1' in sample_metadata and isinstance(sample_metadata['characteristics_ch1'], list):
        for char in sample_metadata['characteristics_ch1']:
            char_lower = char.lower()
            if 'treatment:' in char_lower:
                explicit_treatment = char_lower.split('treatment:')[1].strip()
                found_specific_explicit_treatment = False
                for mapped_treat, keywords in treatment_keywords.items():
                    if mapped_treat == "no stress": continue # Handle "no stress" later
                    if any(k in explicit_treatment for k in keywords):
                        found_treatments.add(mapped_treat)
                        found_specific_explicit_treatment = True
                        break
                if not found_specific_explicit_treatment and explicit_treatment not in ["none", "no", "control", "mock"]:
                    found_treatments.add("other stress")
                # If an explicit treatment was found (either specific or "other stress"), we can stop looking in characteristics
                if found_treatments:
                    break

    # Search in combined text for other treatments
    for mapped_treat, keywords in treatment_keywords.items():
        if mapped_treat == "no stress": # Handle "no stress" last, as it's a default
            continue
        for keyword in keywords:
            if keyword in combined_text:
                found_treatments.add(mapped_treat)
                break
    
    # Special handling for "Other Light Stress"
    # Only add if "light stress" is mentioned but no specific light stress (low/high/red) is found
    if "light stress" in combined_text and not any(t in found_treatments for t in ["low light stress", "high light stress", "red light stress"]):
        found_treatments.add("other light stress")

    # If "no stress" keywords are present and no other stress is found
    if not found_treatments and any(k in combined_text for k in treatment_keywords["no stress"]):
        found_treatments.add("no stress")
    
    # If no specific stress or "no stress" is found, default to "No stress"
    if not found_treatments:
        extracted_data["treatment"] = ["No stress"]
    else:
        # Remove "no stress" if other specific stresses are present
        if len(found_treatments) > 1 and "no stress" in found_treatments:
            found_treatments.remove("no stress")
        
        extracted_data["treatment"] = sorted(list(found_treatments)) # Sort for consistent output

    # --- Extract Medium ---
    medium_keywords = {
        "soil": ["soil", "sunshine mix", "potting mix", "compost", "vermiculite", "sand"],
        "agar": ["agar", "ms medium", "murashige and skoog", "gelrite"],
        "hydroponic": ["hydroponic", "nutrient solution"],
        "rockwool": ["rockwool"]
    }

    found_medium = "unspecified"
    if 'growth_protocol_ch1' in sample_metadata and isinstance(sample_metadata['growth_protocol_ch1'], list):
        growth_protocol_text = " ".join([item.lower() for item in sample_metadata['growth_protocol_ch1'] if isinstance(item, str)])
        for medium_type, keywords in medium_keywords.items():
            if any(k in growth_protocol_text for k in keywords):
                found_medium = medium_type
                break
    
    # If not found in growth_protocol_ch1, check combined text
    if found_medium == "unspecified":
        for medium_type, keywords in medium_keywords.items():
            if any(k in combined_text for k in keywords):
                found_medium = medium_type
                break
    
    extracted_data["medium"] = found_medium

    return extracted_data

def GSE19603_extractor(sample_metadata: dict) -> dict:
    extracted_data = {
        "tissue": "unknown",
        "treatment": [],
        "medium": "unspecified"
    }

    # Helper to safely get text from metadata fields
    def get_text(key):
        value = sample_metadata.get(key)
        if isinstance(value, list):
            return " ".join(value).lower()
        elif isinstance(value, str):
            return value.lower()
        return ""

    # Combine relevant text fields for easier searching
    all_text = (
        get_text('title') + " " +
        get_text('source_name_ch1') + " " +
        get_text('characteristics_ch1') + " " +
        get_text('treatment_protocol_ch1') + " " +
        get_text('growth_protocol_ch1') + " " +
        get_text('description')
    )

    # --- Extract Tissue ---
    tissue_keywords = {
        "root": ["root"],
        "leaf": ["leaf", "leaves"],
        "flower": ["flower", "flowers"],
        "shoot": ["shoot", "above-ground parts", "aerial parts"],
        "rosette": ["rosette"],
        "bud": ["bud"],
        "whole_plant": ["whole plant", "plants"],
        "silique": ["silique"],
        "callus": ["callus"],
        "seed": ["seed", "seeds"],
        "seedling": ["seedling", "seedlings"]
    }

    # Prioritize 'tissue:' in characteristics_ch1
    characteristics_ch1_text = get_text('characteristics_ch1')
    import re
    match = re.search(r'tissue:\s*([\w\s-]+)', characteristics_ch1_text)
    if match:
        found_tissue_str = match.group(1).strip()
        for tissue_enum, keywords in tissue_keywords.items():
            if any(keyword in found_tissue_str for keyword in keywords):
                extracted_data["tissue"] = tissue_enum
                break
    
    # If not found in characteristics_ch1, search in all_text
    if extracted_data["tissue"] == "unknown":
        for tissue_enum, keywords in tissue_keywords.items():
            if any(keyword in all_text for keyword in keywords):
                extracted_data["tissue"] = tissue_enum
                break
    
    # Special case for "plants" without specific tissue, often implies "whole_plant"
    if extracted_data["tissue"] == "unknown" and "plants" in all_text and "tissue:" not in characteristics_ch1_text:
        extracted_data["tissue"] = "whole_plant"


    # --- Extract Treatment ---
    found_treatments = set()
    treatment_mapping = {
        "Drought Stress": ["drought", "water deficit", "dehydration"],
        "Salinity Stress": ["salinity", "salt stress", "nacl"],
        "Heat Stress": ["heat stress", "heat treatment", "37oc", "high temperature"],
        "Cold Stress": ["cold stress", "cold treatment", "low temperature", "4oc"],
        "Chemical Stress": ["chemical stress", "herbicide", "pesticide", "heavy metal", "cadmium", "aluminum", "hormone treatment", "auxin", "cytokinin", "gibberellin", "abscisic acid", "ethylene", "jasmonate", "salicylic acid"],
        "Nutrient Deficiency": ["nutrient deficiency", "nitrogen starvation", "phosphate starvation", "low nutrient", "nutrient deprivation"],
        "Pathogen Attack": ["pathogen", "infection", "bacteria", "virus", "fungus", "disease", "elicitor"],
        "Low Light Stress": ["low light", "shade"],
        "High Light Stress": ["high light", "excess light"],
        "Red Light Stress": ["red light"],
        "Other Light Stress": ["light stress", "uv light", "blue light", "far-red light"],
        "Other stress": ["stress", "stress condition", "stress treatment"], # General stress, if specific not found
        "No stress": ["normal conditions", "untreated", "control", "no stress", "ambient conditions", "standard conditions", "wild type untreated"]
    }

    # Search for specific treatments in all_text
    specific_stress_found = False
    for treatment_enum, keywords in treatment_mapping.items():
        if treatment_enum == "No stress" or treatment_enum == "Other stress":
            continue # Handle these later
        for keyword in keywords:
            if keyword in all_text:
                found_treatments.add(treatment_enum)
                specific_stress_found = True
                break
    
    # If no specific stress found, check for "No stress" indicators
    if not specific_stress_found:
        no_stress_indicator_found = False
        for keyword in treatment_mapping["No stress"]:
            if keyword in all_text:
                found_treatments.add("No stress")
                no_stress_indicator_found = True
                break
        
        # If neither specific stress nor "No stress" indicator, default to "Other stress"
        if not no_stress_indicator_found and not found_treatments:
            found_treatments.add("Other stress")
    
    extracted_data["treatment"] = sorted(list(found_treatments)) # Sort for consistent output

    # --- Extract Medium ---
    medium_keywords = {
        "soil": ["soil", "metro mix", "potting mix", "peat", "compost"],
        "agar": ["agar", "gelrite"],
        "ms medium": ["ms medium", "murashige and skoog"],
        "liquid medium": ["liquid medium", "hydroponic", "hoagland solution", "nutrient solution"],
        "vermiculite": ["vermiculite"],
        "perlite": ["perlite"],
        "rockwool": ["rockwool"],
        "sand": ["sand"]
    }

    found_medium = "unspecified"
    for medium_enum, keywords in medium_keywords.items():
        if any(keyword in all_text for keyword in keywords):
            found_medium = medium_enum
            break
    
    # If still unspecified, and tissue is whole_plant or plants, infer soil if growth protocol suggests it
    if found_medium == "unspecified" and extracted_data["tissue"] in ["whole_plant", "seedling"]:
        if "sown directly on soil" in all_text or "kept in growth chambers" in all_text:
            found_medium = "soil"

    extracted_data["medium"] = found_medium

    return extracted_data

import re

def GSE62163_extractor(sample_metadata: dict) -> dict:
    result = {}

    # Helper to safely get and join text from metadata fields
    def _get_text(data, key):
        value = data.get(key)
        if isinstance(value, list):
            return " ".join(value).lower()
        if isinstance(value, str):
            return value.lower()
        return ""

    # Combine relevant text fields for easier searching
    all_text = (
        _get_text(sample_metadata, 'title') + " " +
        _get_text(sample_metadata, 'source_name_ch1') + " " +
        _get_text(sample_metadata, 'characteristics_ch1') + " " +
        _get_text(sample_metadata, 'treatment_protocol_ch1') + " " +
        _get_text(sample_metadata, 'growth_protocol_ch1') + " " +
        _get_text(sample_metadata, 'description')
    )

    # --- Extract Tissue ---
    extracted_tissue = "unknown"
    tissue_keywords = {
        "shoot": ["shoot", "shoots"],
        "root": ["root", "roots"],
        "leaf": ["leaf", "leaves"],
        "flower": ["flower", "flowers"],
        "rosette": ["rosette", "rosettes"],
        "bud": ["bud", "buds"],
        "whole_plant": ["whole plant", "whole_plants", "plant", "plants"],
        "silique": ["silique", "siliques"],
        "callus": ["callus"],
        "seed": ["seed", "seeds"],
        "seedling": ["seedling", "seedlings"]
    }

    # Priority 1: Check characteristics_ch1 for explicit "tissue:"
    char_ch1_text = _get_text(sample_metadata, 'characteristics_ch1')
    match = re.search(r'tissue:\s*([\w\s]+)', char_ch1_text)
    if match:
        found_tissue_phrase = match.group(1).strip()
        for tissue_enum, keywords in tissue_keywords.items():
            if any(kw in found_tissue_phrase for kw in keywords):
                extracted_tissue = tissue_enum
                break
    
    # Priority 2: If not found or still unknown, search in all_text
    if extracted_tissue == "unknown":
        for tissue_enum, keywords in tissue_keywords.items():
            # Use word boundaries for more precise matching
            if any(re.search(r'\b' + re.escape(kw) + r'\b', all_text) for kw in keywords):
                extracted_tissue = tissue_enum
                break
    
    result['tissue'] = extracted_tissue

    # --- Extract Treatment ---
    extracted_treatments = set()
    treatment_keywords = {
        "Drought Stress": ["drought", "water deficit", "water stress"],
        "Salinity Stress": ["salinity", "salt stress", "nacl"],
        "Heat Stress": ["heat stress", "high temperature", "43°c", "hs", "heat shock"],
        "Cold Stress": ["cold stress", "low temperature", "chilling", "4°c"],
        "Chemical Stress": ["chemical stress", "herbicide", "pesticide", "metal", "cadmium", "aluminum", "ebr", "aba", "methyl jasmonate", "meja", "salycilic acid", "sa", "ethylene", "auxin", "cytokinin", "gibberellin", "brassinosteroid", "hormone"],
        "Nutrient Deficiency": ["nutrient deficiency", "nitrogen starvation", "phosphate starvation", "sulfur starvation", "nutrient deprivation"],
        "Pathogen Attack": ["pathogen", "infection", "fungus", "bacteria", "virus", "biotic stress", "botrytis", "pseudomonas"],
        "Low Light Stress": ["low light", "darkness", "shade"],
        "High Light Stress": ["high light", "uv light"],
        "Red Light Stress": ["red light"],
        "Other Light Stress": ["light stress", "blue light", "far-red light"],
        "Other stress": ["stress", "stressed", "abiotic stress"],
        "No stress": ["no stress", "control", "unstressed", "mock"]
    }

    found_specific_stress = False
    for treatment_enum, keywords in treatment_keywords.items():
        if treatment_enum == "No stress": # Handle "No stress" separately
            continue
        if any(re.search(r'\b' + re.escape(kw) + r'\b', all_text) for kw in keywords):
            extracted_treatments.add(treatment_enum)
            found_specific_stress = True
    
    # If no specific stress was found, check for "No stress" or default to it
    if not found_specific_stress:
        if any(re.search(r'\b' + re.escape(kw) + r'\b', all_text) for kw in treatment_keywords["No stress"]):
            extracted_treatments.add("No stress")
        else:
            # Default to "No stress" if no treatment information is found
            extracted_treatments.add("No stress")
    
    # Refine treatments:
    # If "Other stress" was added but a more specific stress was also found, remove "Other stress"
    if len(extracted_treatments) > 1 and "Other stress" in extracted_treatments:
        extracted_treatments.remove("Other stress")
    
    # If "No stress" is present with other specific stresses, remove "No stress"
    # as the sample is likely treated, and "unstressed" refers to a control group in the study.
    if len(extracted_treatments) > 1 and "No stress" in extracted_treatments:
        extracted_treatments.remove("No stress")

    result['treatment'] = sorted(list(extracted_treatments))

    # --- Extract Medium ---
    extracted_medium = "unspecified"
    medium_keywords = {
        "MS medium": ["ms medium", "murashige and skoog", "ms salt", "ms basal"],
        "soil": ["soil", "potting mix", "vermiculite", "sand"],
        "agar": ["agar"],
        "liquid medium": ["liquid medium", "hydroponic", "hoagland solution", "knop solution"]
    }

    growth_protocol_text = _get_text(sample_metadata, 'growth_protocol_ch1')
    
    # Priority 1: Search in growth_protocol_ch1
    for medium_enum, keywords in medium_keywords.items():
        if any(re.search(r'\b' + re.escape(kw) + r'\b', growth_protocol_text) for kw in keywords):
            extracted_medium = medium_enum
            break
    
    # Priority 2: If not found, search in all_text (less specific)
    if extracted_medium == "unspecified":
        for medium_enum, keywords in medium_keywords.items():
            if any(re.search(r'\b' + re.escape(kw) + r'\b', all_text) for kw in keywords):
                extracted_medium = medium_enum
                break
    
    result['medium'] = extracted_medium

    return result


def GSE112161_extractor(sample_metadata: dict) -> dict:
    extracted_data = {}

    # Helper to safely get the first item from a list or an empty string
    def _get_first_or_empty(data, key):
        if key in data and isinstance(data[key], list) and data[key]:
            return data[key][0]
        return ""

    # Helper to safely get all items from a list or an empty list
    def _get_all_or_empty(data, key):
        if key in data and isinstance(data[key], list):
            return data[key]
        return []

    # Combine all relevant text fields for searching, converting to lower case
    all_text = " ".join([
        _get_first_or_empty(sample_metadata, 'title'),
        _get_first_or_empty(sample_metadata, 'description'),
        _get_first_or_empty(sample_metadata, 'source_name_ch1'),
        " ".join(_get_all_or_empty(sample_metadata, 'characteristics_ch1')),
        _get_first_or_empty(sample_metadata, 'treatment_protocol_ch1'),
        _get_first_or_empty(sample_metadata, 'growth_protocol_ch1')
    ]).lower()

    # --- Extract Tissue ---
    tissue_mapping = {
        "seedling": ["seedling", "seedlings", "4-d-old"],
        "root": ["root"],
        "leaf": ["leaf"],
        "flower": ["flower"],
        "shoot": ["shoot"],
        "rosette": ["rosette"],
        "bud": ["bud"],
        "silique": ["silique"],
        "callus": ["callus"],
        "seed": ["seed"],
        "whole_plant": ["whole plant", "plant", "plants"] # General, should be checked after specific parts
    }
    extracted_tissue = "unknown"
    
    # Iterate through tissue types in a specific order (more specific first)
    ordered_tissue_types = [
        "seedling", "root", "leaf", "flower", "shoot", "rosette", "bud", 
        "silique", "callus", "seed", "whole_plant"
    ]

    for tissue_type in ordered_tissue_types:
        keywords = tissue_mapping[tissue_type]
        for keyword in keywords:
            if keyword in all_text:
                extracted_tissue = tissue_type
                break
        if extracted_tissue != "unknown":
            break
    extracted_data["tissue"] = extracted_tissue

    # --- Extract Treatment ---
    extracted_treatments = set()
    treatment_text = " ".join([
        _get_first_or_empty(sample_metadata, 'title'),
        _get_first_or_empty(sample_metadata, 'description'),
        " ".join(_get_all_or_empty(sample_metadata, 'characteristics_ch1')),
        _get_first_or_empty(sample_metadata, 'treatment_protocol_ch1')
    ]).lower()

    # Define treatment keywords and their corresponding enum values
    treatment_keywords_map = {
        "Heat Stress": ["heat stress", "hs treatment", "37°c", "44°c", "high temperature", "heat regime", "acclimation"],
        "Cold Stress": ["cold stress", "low temperature", "4°c", "chilling"],
        "Drought Stress": ["drought", "water deficit", "osmotic stress", "peg treatment"],
        "Salinity Stress": ["salinity", "salt stress", "nacl"],
        "Chemical Stress": ["chemical stress", "herbicide", "pesticide", "heavy metal", "cadmium", "aluminum", "abscisic acid", "aba treatment", "methyl jasmonate", "meja"],
        "Nutrient Deficiency": ["nutrient deficiency", "starvation", "low nitrogen", "low phosphate", "no nitrogen", "no phosphate", "nutrient deprivation"],
        "Pathogen Attack": ["pathogen", "infection", "fungus", "bacteria", "virus", "elicitor", "flg22", "peptidoglycan"],
        "Low Light Stress": ["low light", "shade"],
        "High Light Stress": ["high light", "excess light"],
        "Red Light Stress": ["red light"],
        "Other Light Stress": ["light stress"], # General light stress if not specific
    }

    found_specific_stress = False
    for treatment_type, keywords in treatment_keywords_map.items():
        for keyword in keywords:
            if keyword in treatment_text:
                extracted_treatments.add(treatment_type)
                found_specific_stress = True
                break # Move to next treatment_type after finding one keyword

    # Handle "No stress" and "Other stress"
    if not found_specific_stress:
        if "control" in treatment_text or "non-treated" in treatment_text or "no stress" in treatment_text:
            extracted_treatments.add("No stress")
        elif "treatment" in treatment_text or "stress" in treatment_text: # If general terms are present but no specific stress
            extracted_treatments.add("Other stress")
        else: # If absolutely nothing indicating stress or control
            extracted_treatments.add("No stress")
    
    # Ensure "No stress" is only present if no other stress is found
    if len(extracted_treatments) > 1 and "No stress" in extracted_treatments:
        extracted_treatments.remove("No stress")
    
    # If "Other stress" is present, but a more specific stress is also present, remove "Other stress"
    if "Other stress" in extracted_treatments and len(extracted_treatments) > 1:
        # Check if any other specific stress (not "No stress") is present
        has_more_specific_stress = False
        for t in extracted_treatments:
            if t != "Other stress" and t != "No stress":
                has_more_specific_stress = True
                break
        if has_more_specific_stress:
            extracted_treatments.remove("Other stress")

    extracted_data["treatment"] = sorted(list(extracted_treatments))

    # --- Extract Medium ---
    extracted_medium = "unspecified"
    growth_protocol = _get_first_or_empty(sample_metadata, 'growth_protocol_ch1').lower()

    # Regex to capture common medium names, optionally followed by parenthesized details
    # The order in the regex matters for more specific terms like "MS medium" vs "medium"
    match = re.search(r'(ms medium|murashige skoog|gm medium|agar|soil|hydroponic|liquid medium|water|medium)(?:\s*\(.*?\))?', growth_protocol)
    if match:
        captured_medium = match.group(1).strip()
        if "murashige skoog" in captured_medium or "ms medium" in captured_medium:
            extracted_medium = "MS medium"
        elif "gm medium" in captured_medium:
            extracted_medium = "GM medium"
        elif "liquid medium" in captured_medium:
            extracted_medium = "liquid medium"
        elif "water" in captured_medium:
            extracted_medium = "water"
        elif "agar" in captured_medium:
            extracted_medium = "agar"
        elif "soil" in captured_medium:
            extracted_medium = "soil"
        elif "hydroponic" in captured_medium:
            extracted_medium = "hydroponic"
        elif "medium" in captured_medium: # Generic "medium" if nothing more specific
            extracted_medium = "unspecified" # Default to unspecified if only "medium" is found without context
    
    # If no specific medium found by regex, check for "soil" as a common default for whole plants
    if extracted_medium == "unspecified" and extracted_data["tissue"] == "whole_plant":
        if "soil" in growth_protocol:
            extracted_medium = "soil"

    extracted_data["medium"] = extracted_medium

    return extracted_data

import re

def GSE60960_extractor(sample_metadata: dict) -> dict:
    result = {
        "tissue": "unknown",
        "treatment": [],
        "medium": "unspecified"
    }

    # Helper to get text from metadata fields, handling lists and joining
    def get_text(field_name):
        value = sample_metadata.get(field_name)
        if value is None:
            return ""
        if isinstance(value, list):
            return " ".join(value).lower()
        return str(value).lower()

    # Combine relevant text fields for comprehensive searching
    all_text = (
        get_text('title') + " " +
        get_text('source_name_ch1') + " " +
        get_text('characteristics_ch1') + " " +
        get_text('treatment_protocol_ch1') + " " +
        get_text('growth_protocol_ch1') + " " +
        get_text('description')
    )

    # --- Extract Tissue ---
    tissue_map = {
        "leaf": ["leaf", "leaves"],
        "root": ["root", "roots"],
        "flower": ["flower", "flowers"],
        "shoot": ["shoot", "shoots"],
        "rosette": ["rosette", "rosettes"],
        "bud": ["bud", "buds"],
        "whole_plant": ["whole plant", "whole-plant"],
        "silique": ["silique", "siliques"],
        "callus": ["callus"],
        "seed": ["seed", "seeds"],
        "seedling": ["seedling", "seedlings"]
    }

    # Prioritize 'tissue: X' in characteristics_ch1
    characteristics_ch1_text = get_text('characteristics_ch1')
    match = re.search(r'tissue:\s*([\w\s-]+)', characteristics_ch1_text) # Allow spaces and hyphens
    if match:
        extracted_tissue = match.group(1).strip().lower()
        for enum_val, keywords in tissue_map.items():
            if extracted_tissue in keywords:
                result["tissue"] = enum_val
                break
        if result["tissue"] == "unknown": # If direct match not in map, try general search for 'plant'
            if "plant" in extracted_tissue:
                result["tissue"] = "whole_plant"

    if result["tissue"] == "unknown":
        # Search in all_text for other tissue keywords
        for enum_val, keywords in tissue_map.items():
            for keyword in keywords:
                if re.search(r'\b' + re.escape(keyword) + r'\b', all_text):
                    result["tissue"] = enum_val
                    break
            if result["tissue"] != "unknown":
                break
        # Special case for 'plant' if no specific tissue found
        if result["tissue"] == "unknown" and re.search(r'\bplant\b', all_text) and not any(
            re.search(r'\b' + re.escape(k) + r'\b', all_text) for k in ["leaf", "root", "flower", "seed", "seedling", "shoot", "rosette", "bud", "silique", "callus"]
        ):
             result["tissue"] = "whole_plant"


    # --- Extract Treatment ---
    treatments_found = set()
    treatment_keywords = {
        "Drought Stress": ["drought", "water deficit", "no water", "dry", "dehydration", "water stress"],
        "Salinity Stress": ["salinity", "salt", "nacl"],
        "Heat Stress": ["heat stress", "high temperature", "heat shock"],
        "Cold Stress": ["cold stress", "low temperature", "chilling"],
        "Chemical Stress": ["chemical stress", "herbicide", "pesticide", "heavy metal", "cadmium", "arsenic", "ozone", "oxidative stress", "paraquat", "abscisic acid", "aba", "methyl jasmonate", "meja", "ethylene", "auxin", "cytokinin", "gibberellin", "brassinosteroid", "strigolactone", "salicylic acid", "sa", "jasmonic acid", "ja", "hydrogen peroxide", "h2o2"],
        "Nutrient Deficiency": ["nutrient deficiency", "nitrogen starvation", "phosphate starvation", "iron deficiency", "nutrient limitation", "low nitrogen", "low phosphate"],
        "Pathogen Attack": ["pathogen", "infection", "bacteria", "virus", "fungus", "insect", "herbivory", "elicitor", "flg22", "chitin", "pest"],
        "Low Light Stress": ["low light", "shade"],
        "High Light Stress": ["high light", "excess light"],
        "Red Light Stress": ["red light"],
        "Other Light Stress": ["light stress", "moderate light", "uv-b", "uvb", "darkness", "dark treatment"],
    }

    # Search in specific fields first for stronger signals
    search_fields_for_treatment = [
        get_text('title'),
        get_text('source_name_ch1'),
        get_text('treatment_protocol_ch1'),
        get_text('description'),
        get_text('characteristics_ch1')
    ]

    for field_text in search_fields_for_treatment:
        for enum_val, keywords in treatment_keywords.items():
            for keyword in keywords:
                if re.search(r'\b' + re.escape(keyword) + r'\b', field_text):
                    treatments_found.add(enum_val)

    # Handle "No stress" logic:
    # If no specific stress keywords were found, then it's either a control or standard conditions.
    if not treatments_found:
        # Check for explicit mentions of control/untreated/normal conditions
        if re.search(r'\b(control|untreated|normal conditions|standard conditions)\b', all_text):
            treatments_found.add("No stress")
        else:
            # If no stress keywords and no explicit "control", default to "No stress"
            treatments_found.add("No stress")
    # If specific stresses are found, "No stress" is not added, as the sample itself is stressed.

    result["treatment"] = sorted(list(treatments_found)) # Sort for consistent output

    # --- Extract Medium ---
    medium_keywords = {
        "soil": ["soil", "potting mix", "substrate", "compost"],
        "agar": ["agar", "ms medium", "murashige skoog", "gel", "petri dish"],
        "hydroponic": ["hydroponic", "liquid medium", "nutrient solution", "water culture"],
        "vermiculite": ["vermiculite"],
        "perlite": ["perlite"]
    }

    # Search in growth_protocol_ch1, treatment_protocol_ch1, source_name_ch1, description, characteristics_ch1
    search_fields_for_medium = [
        get_text('growth_protocol_ch1'),
        get_text('treatment_protocol_ch1'),
        get_text('source_name_ch1'),
        get_text('description'),
        get_text('characteristics_ch1')
    ]

    for field_text in search_fields_for_medium:
        for enum_val, keywords in medium_keywords.items():
            for keyword in keywords:
                if re.search(r'\b' + re.escape(keyword) + r'\b', field_text):
                    result["medium"] = enum_val
                    break
            if result["medium"] != "unspecified":
                break
        if result["medium"] != "unspecified":
            break

    return result

import re

def GSE83136_extractor(sample_metadata: dict) -> dict:
    extracted_data = {}

    def _get_first_or_empty(data, key):
        """Safely retrieves the first element of a list from the dictionary, or an empty string."""
        if key in data and data[key]:
            return str(data[key][0])
        return ""

    def _get_joined_list_or_empty(data, key):
        """Safely retrieves and joins elements of a list from the dictionary, or an empty string."""
        if key in data and data[key]:
            return " ".join(map(str, data[key]))
        return ""

    # Combine relevant text fields for easier searching, converting to lowercase
    growth_protocol = _get_first_or_empty(sample_metadata, 'growth_protocol_ch1').lower()
    treatment_protocol = _get_first_or_empty(sample_metadata, 'treatment_protocol_ch1').lower()
    characteristics = _get_joined_list_or_empty(sample_metadata, 'characteristics_ch1').lower()
    title = _get_first_or_empty(sample_metadata, 'title').lower()
    description = _get_first_or_empty(sample_metadata, 'description').lower()
    source_name = _get_first_or_empty(sample_metadata, 'source_name_ch1').lower()

    all_text = f"{growth_protocol} {treatment_protocol} {characteristics} {title} {description} {source_name}"

    # 1. Extract Tissue
    tissue = "unknown"
    tissue_keywords = {
        "seedling": ["seedling", "seedlings"],
        "root": ["root"],
        "leaf": ["leaf", "leaves"],
        "flower": ["flower", "flowers"],
        "shoot": ["shoot"],
        "rosette": ["rosette"],
        "bud": ["bud"],
        "silique": ["silique"],
        "callus": ["callus"],
        "seed": ["seed", "seeds"]
    }

    for t_type, keywords in tissue_keywords.items():
        for keyword in keywords:
            if keyword in all_text:
                tissue = t_type
                break
        if tissue != "unknown":
            break
    
    # Special case for "whole_plant" if "plants" is mentioned without other specific tissues
    if tissue == "unknown" and "plants" in all_text and not any(k in all_text for k_list in tissue_keywords.values() for k in k_list if k != "plants"):
        tissue = "whole_plant"

    extracted_data["tissue"] = tissue

    # 2. Extract Treatment
    treatments = set()
    treatment_mappings = {
        "Drought Stress": ["drought", "water deficit", "water stress"],
        "Salinity Stress": ["salinity", "salt stress", "nacl"],
        "Heat Stress": ["heat stress", "hs treatment", "heat regime", "high temperature", "heat shock"],
        "Cold Stress": ["cold stress", "low temperature", "chilling", "freezing"],
        "Chemical Stress": ["chemical stress", "herbicide", "pesticide", "cadmium", "heavy metal", "paraquat", "hormone treatment", "abscisic acid", "auxin", "cytokinin", "gibberellin", "ethylene"],
        "Nutrient Deficiency": ["nutrient deficiency", "low nitrogen", "low phosphate", "nitrogen starvation", "phosphate starvation", "nutrient starvation"],
        "Pathogen Attack": ["pathogen", "infection", "virus", "bacteria", "fungus", "elicitor", "disease"],
        "Low Light Stress": ["low light"],
        "High Light Stress": ["high light"],
        "Red Light Stress": ["red light"],
        "Other Light Stress": ["light stress"], # Catch-all for light if not specific
    }

    for t_type, keywords in treatment_mappings.items():
        for keyword in keywords:
            if keyword in all_text:
                treatments.add(t_type)
                break
    
    # Handle "Other stress" if a treatment is mentioned but not specific
    if "treatment:" in characteristics and not treatments:
        treatments.add("Other stress")

    if not treatments:
        treatments.add("No stress")
    
    # Remove "Other Light Stress" if a more specific light stress is present
    if "Other Light Stress" in treatments and any(s in treatments for s in ["Low Light Stress", "High Light Stress", "Red Light Stress"]):
        treatments.remove("Other Light Stress")

    extracted_data["treatment"] = sorted(list(treatments)) # Sort for consistent output

    # 3. Extract Medium
    medium = "unspecified"
    
    # Prioritize specific mentions
    if "gm medium" in growth_protocol:
        medium = "GM medium"
    elif "ms medium" in growth_protocol:
        medium = "MS medium"
    elif "agar" in growth_protocol:
        medium = "agar"
    elif "soil" in growth_protocol:
        medium = "soil"
    elif "liquid medium" in growth_protocol:
        medium = "liquid medium"
    elif "hoagland" in growth_protocol:
        medium = "Hoagland solution"
    elif "hydroponic" in growth_protocol:
        medium = "hydroponic solution"
    
    # Fallback to extract any "X medium" if not specifically identified
    if medium == "unspecified":
        match = re.search(r'(\w+\s*medium)', growth_protocol)
        if match:
            medium = match.group(1).strip()
            # Standardize common medium names
            if medium.lower() == "ms medium":
                medium = "MS medium"
            elif medium.lower() == "gm medium":
                medium = "GM medium"
            # Add more specific capitalization rules if needed.

    # If still unspecified and tissue is whole_plant, infer soil as a common default
    if medium == "unspecified" and tissue == "whole_plant":
        medium = "soil"

    extracted_data["medium"] = medium

    return extracted_data

def GSE49418_extractor(sample_metadata: dict) -> dict:
    """
    Extracts tissue, treatment, and medium information from a sample_metadata dictionary
    following a predefined JSON schema.
    """

    extracted_data = {
        "tissue": "unknown",
        "treatment": [],
        "medium": "unspecified"
    }

    # Helper to safely get text from potential list fields and convert to lowercase
    def get_text(key):
        value = sample_metadata.get(key)
        if isinstance(value, list):
            return " ".join(value).lower()
        if isinstance(value, str):
            return value.lower()
        return ""

    # --- Extract Tissue ---
    # Prioritize more specific tissue types or common ones for this dataset
    tissue_found = False
    for field in ['source_name_ch1', 'characteristics_ch1', 'title', 'description', 'treatment_protocol_ch1', 'extract_protocol_ch1']:
        field_text = get_text(field)
        if "seedlings" in field_text or "seedling" in field_text:
            extracted_data["tissue"] = "seedling"
            tissue_found = True
            break
        elif "whole plant" in field_text:
            extracted_data["tissue"] = "whole_plant"
            tissue_found = True
            break
        elif "root" in field_text:
            extracted_data["tissue"] = "root"
            tissue_found = True
            break
        elif "leaf" in field_text:
            extracted_data["tissue"] = "leaf"
            tissue_found = True
            break
        elif "flower" in field_text:
            extracted_data["tissue"] = "flower"
            tissue_found = True
            break
        elif "shoot" in field_text:
            extracted_data["tissue"] = "shoot"
            tissue_found = True
            break
        elif "rosette" in field_text:
            extracted_data["tissue"] = "rosette"
            tissue_found = True
            break
        elif "bud" in field_text:
            extracted_data["tissue"] = "bud"
            tissue_found = True
            break
        elif "silique" in field_text:
            extracted_data["tissue"] = "silique"
            tissue_found = True
            break
        elif "callus" in field_text:
            extracted_data["tissue"] = "callus"
            tissue_found = True
            break
        elif "seed" in field_text:
            extracted_data["tissue"] = "seed"
            tissue_found = True
            break

    # --- Extract Treatment ---
    found_treatments = set()
    
    # Determine if the sample itself is explicitly a control/no stress sample
    is_control_sample = False
    if "ck" in get_text('title') or "control" in get_text('title'):
        is_control_sample = True
    if "without treatment" in get_text('source_name_ch1') or "none" in get_text('characteristics_ch1'):
        is_control_sample = True
    if "no treatment" in get_text('description'):
        is_control_sample = True

    if is_control_sample:
        found_treatments.add("No stress")

    treatment_mapping = {
        "drought": "Drought Stress",
        "dehydration": "Drought Stress",
        "salinity": "Salinity Stress",
        "salt": "Salinity Stress",
        "heat": "Heat Stress",
        "cold": "Cold Stress",
        "chemical": "Chemical Stress",
        "nutrient deficiency": "Nutrient Deficiency",
        "pathogen": "Pathogen Attack",
        "low light": "Low Light Stress",
        "high light": "High Light Stress",
        "red light": "Red Light Stress",
        "light stress": "Other Light Stress",
        "stress": "Other stress" # General stress, should be considered last
    }
    
    # Fields to search for treatments
    treatment_fields = [
        'characteristics_ch1',
        'treatment_protocol_ch1',
        'title',
        'description',
        'source_name_ch1'
    ]

    for field in treatment_fields:
        field_text = get_text(field)
        for keyword, treatment_type in treatment_mapping.items():
            # Ensure "stress" as a general keyword doesn't override more specific ones
            if keyword == "stress" and any(specific_k in field_text for specific_k in treatment_mapping if specific_k != "stress" and specific_k in field_text):
                continue # Skip general "stress" if a more specific stress is also present
            if keyword in field_text:
                found_treatments.add(treatment_type)
    
    # Refine treatment logic based on control status
    if is_control_sample and len(found_treatments) > 1:
        # If the sample is explicitly a control, and other stresses are mentioned (likely in protocols for the experiment),
        # then the treatment for *this specific sample* is "No stress".
        extracted_data["treatment"] = ["No stress"]
    elif found_treatments:
        # If "No stress" is present along with other specific stresses, remove "No stress"
        # unless it's the only treatment found. A sample is either stressed or not.
        if "No stress" in found_treatments and len(found_treatments) > 1:
            found_treatments.remove("No stress")
        extracted_data["treatment"] = sorted(list(found_treatments))
    else:
        # If no treatments are found at all, assume "No stress" as a default for a biological sample.
        extracted_data["treatment"] = ["No stress"]

    # --- Extract Medium ---
    medium_keywords = {
        "soil": "soil",
        "ms medium": "MS medium",
        "murashige skoog": "MS medium",
        "agar": "agar",
        "liquid medium": "liquid medium",
        "hydroponic": "hydroponic solution",
        "vermiculite": "vermiculite",
        "peat": "peat",
        "sand": "sand",
        "rockwool": "rockwool"
    }
    
    medium_found = False
    # Search in growth_protocol_ch1 first, then treatment_protocol_ch1, then characteristics_ch1
    for field in ['growth_protocol_ch1', 'treatment_protocol_ch1', 'characteristics_ch1']:
        field_text = get_text(field)
        for keyword, medium_type in medium_keywords.items():
            if keyword in field_text:
                extracted_data["medium"] = medium_type
                medium_found = True
                break
        if medium_found:
            break
    
    # "Whatman paper" is typically a treatment condition, not a growth medium.
    # If no specific growth medium is found, "unspecified" is the default.

    return extracted_data

def GSE10670_extractor(sample_metadata: dict) -> dict:
    def _get_text_from_metadata(metadata_dict, keys):
        text_parts = []
        for key in keys:
            if key in metadata_dict:
                value = metadata_dict[key]
                if isinstance(value, list):
                    text_parts.extend(value)
                elif isinstance(value, str):
                    text_parts.append(value)
        return " ".join(text_parts).lower()

    result = {}

    # --- Extract Tissue ---
    tissue_enums = ["root", "leaf", "flower", "shoot", "rosette", "bud", "silique", "callus", "seed", "seedling"]
    found_tissue = "unknown"
    search_text = _get_text_from_metadata(sample_metadata, ['source_name_ch1', 'characteristics_ch1', 'title', 'description'])

    # Prioritize specific tissues
    for tissue_type in tissue_enums:
        # Handle "whole_plant" separately as it's a more general term
        if tissue_type == "whole_plant":
            continue
        # Check for exact match or common variations (e.g., "leaf" or "leaves")
        if tissue_type in search_text or tissue_type.replace("_", " ") in search_text or tissue_type + "s" in search_text:
            found_tissue = tissue_type
            break
    
    # If no specific tissue found, check for whole_plant
    if found_tissue == "unknown":
        if "whole plant" in search_text or "whole-plant" in search_text:
            found_tissue = "whole_plant"
        elif "plant" in search_text and not any(t in search_text for t in tissue_enums if t != "whole_plant"):
            # If "plant" is mentioned, and no other specific tissue, assume whole_plant
            found_tissue = "whole_plant"

    result["tissue"] = found_tissue

    # --- Extract Treatment ---
    treatment_keywords_map = {
        "Drought Stress": ["drought", "water deficit", "without water", "no water", "water deprivation"],
        "Salinity Stress": ["salinity", "salt stress", "nacl"],
        "Heat Stress": ["heat", "high temperature", "heat shock"],
        "Cold Stress": ["cold", "low temperature", "chilling"],
        "Chemical Stress": ["chemical", "herbicide", "pesticide", "heavy metal", "cadmium", "aluminum", "hormone treatment", "abscisic acid", "aba", "methyl jasmonate", "meja", "auxin", "cytokinin", "gibberellin", "ethylene", "osmotic stress"],
        "Nutrient Deficiency": ["nutrient deficiency", "nitrogen starvation", "phosphate starvation", "iron deficiency", "low nitrogen", "low phosphate", "nutrient deprivation", "starvation"],
        "Pathogen Attack": ["pathogen", "infection", "fungus", "bacteria", "virus", "insect", "elicitor", "disease", "inoculation"],
        "Low Light Stress": ["low light", "shade"],
        "High Light Stress": ["high light", "excess light"],
        "Red Light Stress": ["red light"],
    }
    
    found_treatments = set()
    search_text_treatment = _get_text_from_metadata(sample_metadata, ['source_name_ch1', 'characteristics_ch1', 'title', 'description', 'treatment_protocol_ch1'])

    for treatment_type, keywords in treatment_keywords_map.items():
        for keyword in keywords:
            if keyword in search_text_treatment:
                found_treatments.add(treatment_type)
                break # Found one keyword for this treatment type, move to next treatment type

    # Handle generic light stress
    if "light stress" in search_text_treatment:
        if not any(t in found_treatments for t in ["Low Light Stress", "High Light Stress", "Red Light Stress"]):
            found_treatments.add("Other Light Stress")
    
    # Handle generic stress
    if "stress" in search_text_treatment:
        # Check if any specific stress (excluding light stresses) is already found
        specific_stresses_found = any(t in found_treatments for t in treatment_keywords_map if t not in ["Low Light Stress", "High Light Stress", "Red Light Stress"])
        if not specific_stresses_found and "Other Light Stress" not in found_treatments:
            found_treatments.add("Other stress")

    if not found_treatments:
        # Check for explicit "no stress", "control", "untreated"
        if "no stress" in search_text_treatment or "control" in search_text_treatment or "untreated" in search_text_treatment:
            found_treatments.add("No stress")
        else:
            found_treatments.add("No stress") # Default if nothing found

    # Remove "No stress" if other stresses are present
    if len(found_treatments) > 1 and "No stress" in found_treatments:
        found_treatments.remove("No stress")

    result["treatment"] = sorted(list(found_treatments))

    # --- Extract Medium ---
    medium_keywords = {
        "soil": ["soil", "potting mix", "compost", "peat moss", "vermiculite-perlite mix"],
        "agar": ["agar", "ms medium", "murashige skoog", "gelrite", "phytoagar"],
        "liquid medium": ["liquid medium", "hydroponic", "hoagland", "knop solution", "nutrient solution", "liquid culture"],
        "vermiculite": ["vermiculite"],
        "peat": ["peat"],
        "sand": ["sand"],
        "water": ["water"] # Less specific, check after others
    }
    found_medium = "unspecified"
    search_text_medium = _get_text_from_metadata(sample_metadata, ['growth_protocol_ch1', 'characteristics_ch1', 'source_name_ch1', 'description'])

    for medium_type, keywords in medium_keywords.items():
        for keyword in keywords:
            if keyword in search_text_medium:
                found_medium = medium_type
                break
        if found_medium != "unspecified":
            break
    
    # Inference for whole_plant if not specified
    if found_medium == "unspecified" and result.get("tissue") == "whole_plant":
        found_medium = "soil" # Common default for whole plants

    result["medium"] = found_medium

    return result

def GSE16765_extractor(sample_metadata: dict) -> dict:
    extracted_data = {
        "tissue": "unknown",
        "treatment": ["No stress"],
        "medium": "unspecified"
    }

    # Helper to safely get text from metadata fields, handling lists and missing keys
    def get_text(key):
        value = sample_metadata.get(key)
        if isinstance(value, list):
            return " ".join(value).lower()
        return str(value).lower() if value is not None else ""

    # Combine all relevant text for keyword searching, converting to lowercase
    all_text = (
        get_text('title') + " " +
        get_text('source_name_ch1') + " " +
        get_text('characteristics_ch1') + " " +
        get_text('treatment_protocol_ch1') + " " +
        get_text('growth_protocol_ch1') + " " +
        get_text('description')
    )

    # --- Extract Tissue ---
    tissue_found_direct = False
    for char_entry in sample_metadata.get('characteristics_ch1', []):
        char_entry_lower = char_entry.lower()
        if char_entry_lower.startswith('tissue:'):
            tissue_value = char_entry_lower.split(':', 1)[1].strip()
            # Map to schema enum values
            if tissue_value in ["root", "leaf", "flower", "shoot", "rosette", "bud", "silique", "callus", "seed", "seedling"]:
                extracted_data["tissue"] = tissue_value
                tissue_found_direct = True
                break
            elif "whole plant" in tissue_value or "plant" in tissue_value:
                extracted_data["tissue"] = "whole_plant"
                tissue_found_direct = True
                break
    
    # If tissue not found directly, try to infer from other text (prioritize specific over general)
    if not tissue_found_direct:
        if "leaf" in all_text:
            extracted_data["tissue"] = "leaf"
        elif "root" in all_text:
            extracted_data["tissue"] = "root"
        elif "flower" in all_text:
            extracted_data["tissue"] = "flower"
        elif "shoot" in all_text:
            extracted_data["tissue"] = "shoot"
        elif "rosette" in all_text:
            extracted_data["tissue"] = "rosette"
        elif "bud" in all_text:
            extracted_data["tissue"] = "bud"
        elif "silique" in all_text:
            extracted_data["tissue"] = "silique"
        elif "callus" in all_text:
            extracted_data["tissue"] = "callus"
        elif "seedling" in all_text:
            extracted_data["tissue"] = "seedling"
        elif "seed" in all_text:
            extracted_data["tissue"] = "seed"
        elif "plant" in all_text: # General "plant" implies whole_plant if nothing more specific
            extracted_data["tissue"] = "whole_plant"


    # --- Extract Treatment ---
    found_treatments = set()
    
    # Define treatment keywords (more specific keywords should be checked first implicitly by the order of the map,
    # or by ensuring specific keywords don't overlap with general ones if possible)
    TREATMENT_KEYWORDS_MAP = {
        "Salinity Stress": ["nacl", "salt stress", "salinity", "salt treatment"],
        "Drought Stress": ["drought", "water deficit", "water stress", "dehydration"],
        "Heat Stress": ["heat stress", "high temperature", "hot stress", "heat shock"],
        "Cold Stress": ["cold stress", "low temperature", "chilling", "cold treatment"],
        "Nutrient Deficiency": ["nutrient deficiency", "starvation", "low nitrogen", "nitrogen deficiency", "low phosphate", "phosphate deficiency", "low potassium", "potassium deficiency", "boron deficiency", "iron deficiency"],
        "Pathogen Attack": ["pathogen", "infection", "virus", "bacteria", "fungus", "disease", "biotic stress"],
        "Low Light Stress": ["low light", "shade", "darkness"],
        "High Light Stress": ["high light", "excess light", "uv-b"],
        "Red Light Stress": ["red light"],
        "Other Light Stress": ["light stress", "blue light", "far-red light"],
        "Chemical Stress": ["chemical stress", "herbicide", "pesticide", "heavy metal", "cadmium", "arsenic", "aluminum", "lead", "mercury", "zinc toxicity", "copper toxicity", "chemical treatment"],
        "Other stress": ["stress", "abiotic stress"] # General stress, last resort
    }

    for stress_type, keywords in TREATMENT_KEYWORDS_MAP.items():
        for keyword in keywords:
            if keyword in all_text:
                found_treatments.add(stress_type)
                # Do not break, as multiple treatments can apply

    # Post-processing for "Other stress" and "No stress"
    # If more specific stresses are found, remove "Other stress"
    if len(found_treatments) > 1 and "Other stress" in found_treatments:
        found_treatments.remove("Other stress")
    
    if not found_treatments:
        extracted_data["treatment"] = ["No stress"]
    else:
        extracted_data["treatment"] = sorted(list(found_treatments)) # Sort for consistent output

    # --- Extract Medium ---
    medium_text = get_text('growth_protocol_ch1')
    
    # Prioritize common growth media
    if "soil" in medium_text or "planting medium" in medium_text or "potting mix" in medium_text or "vermiculite" in medium_text or "peat" in medium_text or "sand" in medium_text:
        extracted_data["medium"] = "soil"
    elif "agar" in medium_text or "ms medium" in medium_text or "murashige skoog" in medium_text or "gelrite" in medium_text or "in vitro" in medium_text:
        extracted_data["medium"] = "agar"
    elif "hydroponics" in medium_text or "hoagland solution" in medium_text or "nutrient solution" in medium_text:
        extracted_data["medium"] = "hydroponics"
    # Default is "unspecified" which is already set

    return extracted_data

def GSE72050_extractor(sample_metadata: dict) -> dict:
    result = {
        "tissue": "unknown",
        "treatment": set(),  # Use a set to store unique treatments
        "medium": "unspecified"
    }

    # Helper to get joined string from a list value in metadata
    def _get_joined_text(key):
        return " ".join(sample_metadata.get(key, [])).lower()

    # Combine relevant text fields for easier searching
    all_text = (
        _get_joined_text('title') + " " +
        _get_joined_text('source_name_ch1') + " " +
        _get_joined_text('characteristics_ch1') + " " +
        _get_joined_text('description') + " " +
        _get_joined_text('treatment_protocol_ch1') + " " +
        _get_joined_text('growth_protocol_ch1')
    )

    # --- Extract Tissue ---
    # Prioritize 'characteristics_ch1'
    for char in sample_metadata.get('characteristics_ch1', []):
        if char.lower().startswith('tissue:'):
            tissue_val = char.split(':', 1)[1].strip().lower()
            if 'leaf' in tissue_val: result['tissue'] = 'leaf'
            elif 'root' in tissue_val: result['tissue'] = 'root'
            elif 'flower' in tissue_val: result['tissue'] = 'flower'
            elif 'shoot' in tissue_val: result['tissue'] = 'shoot'
            elif 'rosette' in tissue_val: result['tissue'] = 'rosette'
            elif 'bud' in tissue_val: result['tissue'] = 'bud'
            elif 'whole plant' in tissue_val: result['tissue'] = 'whole_plant'
            elif 'silique' in tissue_val: result['tissue'] = 'silique'
            elif 'callus' in tissue_val: result['tissue'] = 'callus'
            elif 'seedling' in tissue_val: result['tissue'] = 'seedling'
            elif 'seed' in tissue_val: result['tissue'] = 'seed'
            break # Found tissue, no need to check further in characteristics

    # If not found in characteristics, search in combined text
    if result['tissue'] == 'unknown':
        if 'leaf' in all_text: result['tissue'] = 'leaf'
        elif 'root' in all_text: result['tissue'] = 'root'
        elif 'flower' in all_text: result['tissue'] = 'flower'
        elif 'shoot' in all_text: result['tissue'] = 'shoot'
        elif 'rosette' in all_text: result['tissue'] = 'rosette'
        elif 'bud' in all_text: result['tissue'] = 'bud'
        elif 'whole plant' in all_text: result['tissue'] = 'whole_plant'
        elif 'silique' in all_text: result['tissue'] = 'silique'
        elif 'callus' in all_text: result['tissue'] = 'callus'
        elif 'seedling' in all_text: result['tissue'] = 'seedling'
        elif 'seed' in all_text: result['tissue'] = 'seed'


    # --- Extract Treatment ---
    temp_treatments = set()
    char_stress_info_found = False # Flag to indicate if 'stress:' characteristic was present

    # Prioritize 'characteristics_ch1' for stress
    for char in sample_metadata.get('characteristics_ch1', []):
        if char.lower().startswith('stress:'):
            char_stress_info_found = True
            stress_val = char.split(':', 1)[1].strip().lower()
            if 'drought' in stress_val: temp_treatments.add('Drought Stress')
            elif 'salinity' in stress_val or 'salt' in stress_val: temp_treatments.add('Salinity Stress')
            elif 'heat' in stress_val or 'high temperature' in stress_val: temp_treatments.add('Heat Stress')
            elif 'cold' in stress_val or 'low temperature' in stress_val: temp_treatments.add('Cold Stress')
            elif 'chemical' in stress_val: temp_treatments.add('Chemical Stress')
            elif 'nutrient deficiency' in stress_val: temp_treatments.add('Nutrient Deficiency')
            elif 'pathogen' in stress_val: temp_treatments.add('Pathogen Attack')
            elif 'low light' in stress_val: temp_treatments.add('Low Light Stress')
            elif 'high light' in stress_val: temp_treatments.add('High Light Stress')
            elif 'red light' in stress_val: temp_treatments.add('Red Light Stress')
            elif 'light stress' in stress_val: temp_treatments.add('Other Light Stress')
            elif 'normal condition' in stress_val or 'no stress' in stress_val or 'control' in stress_val:
                temp_treatments.add('No stress')
            else: # If it's a stress but not one of the specific enums
                temp_treatments.add('Other stress')
            break # Found stress info in characteristics, prioritize it.

    # If no specific stress was found in characteristics, or if 'No stress' was the only thing found there,
    # then search in the combined text to add more specific stresses.
    if not char_stress_info_found or ('No stress' in temp_treatments and len(temp_treatments) == 1):
        if 'drought' in all_text or 'without watering' in all_text: temp_treatments.add('Drought Stress')
        if 'salinity' in all_text or 'salt' in all_text: temp_treatments.add('Salinity Stress')
        if 'heat' in all_text or 'high temperature' in all_text: temp_treatments.add('Heat Stress')
        if 'cold' in all_text or 'low temperature' in all_text: temp_treatments.add('Cold Stress')
        if 'chemical' in all_text or 'herbicide' in all_text or 'pesticide' in all_text or 'heavy metal' in all_text: temp_treatments.add('Chemical Stress')
        if 'nutrient deficiency' in all_text or 'low nitrogen' in all_text or 'low phosphate' in all_text: temp_treatments.add('Nutrient Deficiency')
        if 'pathogen' in all_text or 'fungus' in all_text or 'bacteria' in all_text or 'virus' in all_text or 'infection' in all_text: temp_treatments.add('Pathogen Attack')
        if 'low light' in all_text or 'darkness' in all_text: temp_treatments.add('Low Light Stress')
        if 'high light' in all_text: temp_treatments.add('High Light Stress')
        if 'red light' in all_text: temp_treatments.add('Red Light Stress')
        if 'light stress' in all_text and not any(s in all_text for s in ['low light', 'high light', 'red light']): temp_treatments.add('Other Light Stress')
        # General 'stress' keyword, only add if no specific stress was found yet and no other treatments are present
        if 'stress' in all_text and not temp_treatments: temp_treatments.add('Other stress')
        # If 'normal condition' or 'control' is mentioned and no other stress was found
        if not temp_treatments and ('normal condition' in all_text or 'control' in all_text or 'no stress' in all_text):
            temp_treatments.add('No stress')

    # Final cleanup for 'No stress'
    # If other specific stresses are present, 'No stress' is redundant
    if 'No stress' in temp_treatments and len(temp_treatments) > 1:
        temp_treatments.remove('No stress')

    # If no stress was found at all, default to 'No stress'
    if not temp_treatments:
        temp_treatments.add('No stress')

    result['treatment'] = sorted(list(temp_treatments))


    # --- Extract Medium ---
    growth_protocol_text = _get_joined_text('growth_protocol_ch1')
    treatment_protocol_text = _get_joined_text('treatment_protocol_ch1')

    if 'soil' in growth_protocol_text or 'soil' in treatment_protocol_text:
        result['medium'] = 'soil'
    elif 'agar' in growth_protocol_text or 'agar' in treatment_protocol_text or \
         'ms medium' in growth_protocol_text or 'ms medium' in treatment_protocol_text or \
         'murashige and skoog' in growth_protocol_text or 'murashige and skoog' in treatment_protocol_text:
        result['medium'] = 'agar'
    elif 'liquid medium' in growth_protocol_text or 'liquid medium' in treatment_protocol_text or \
         'hydroponics' in growth_protocol_text or 'hydroponics' in treatment_protocol_text or \
         'liquid culture' in growth_protocol_text or 'liquid culture' in treatment_protocol_text:
        result['medium'] = 'liquid medium'
    # Infer 'soil' if tissue is 'whole_plant' and medium is still unspecified
    elif result['tissue'] == 'whole_plant':
        result['medium'] = 'soil'

    return result


def GSE19265_extractor(sample_metadata: dict) -> dict:
    extracted_data = {
        "tissue": "unknown",
        "treatment": [],
        "medium": "unspecified"
    }

    def get_joined_text(key):
        return " ".join(sample_metadata.get(key, [])).lower()

    # Combine relevant text fields for comprehensive searching
    all_text = (
        get_joined_text('title') + " " +
        get_joined_text('source_name_ch1') + " " +
        get_joined_text('characteristics_ch1') + " " +
        get_joined_text('growth_protocol_ch1') + " " +
        get_joined_text('description')
    )

    # --- Extract Tissue ---
    characteristics_ch1_list = sample_metadata.get('characteristics_ch1', [])
    for char_item in characteristics_ch1_list:
        char_item_lower = char_item.lower()
        if char_item_lower.startswith('tissue:'):
            tissue_str = char_item_lower.split(':', 1)[1].strip()
            if "whole seedlings" in tissue_str or "seedlings" in tissue_str:
                extracted_data["tissue"] = "seedling"
                break
            elif "whole plant" in tissue_str:
                extracted_data["tissue"] = "whole_plant"
                break
            elif "root" in tissue_str:
                extracted_data["tissue"] = "root"
                break
            elif "leaf" in tissue_str:
                extracted_data["tissue"] = "leaf"
                break
            elif "flower" in tissue_str:
                extracted_data["tissue"] = "flower"
                break
            elif "shoot" in tissue_str:
                extracted_data["tissue"] = "shoot"
                break
            elif "rosette" in tissue_str:
                extracted_data["tissue"] = "rosette"
                break
            elif "bud" in tissue_str:
                extracted_data["tissue"] = "bud"
                break
            elif "silique" in tissue_str:
                extracted_data["tissue"] = "silique"
                break
            elif "callus" in tissue_str:
                extracted_data["tissue"] = "callus"
                break
            elif "seed" in tissue_str:
                extracted_data["tissue"] = "seed"
                break

    # --- Extract Treatment ---
    found_treatments = set()
    
    # Prioritize explicit "treatment:" in characteristics_ch1
    for char_item in characteristics_ch1_list:
        char_item_lower = char_item.lower()
        if char_item_lower.startswith('treatment:'):
            treatment_str = char_item_lower.split(':', 1)[1].strip()
            if "dexamethasone" in treatment_str or "chemical" in treatment_str:
                found_treatments.add("Chemical Stress")
            elif "drought" in treatment_str:
                found_treatments.add("Drought Stress")
            elif "salinity" in treatment_str or "salt" in treatment_str:
                found_treatments.add("Salinity Stress")
            elif "heat" in treatment_str:
                found_treatments.add("Heat Stress")
            elif "cold" in treatment_str:
                found_treatments.add("Cold Stress")
            elif "nutrient deficiency" in treatment_str:
                found_treatments.add("Nutrient Deficiency")
            elif "biotic" in treatment_str:
                found_treatments.add("Biotic Stress")
            elif "low light" in treatment_str:
                found_treatments.add("Low Light Stress")
            elif "high light" in treatment_str:
                found_treatments.add("High Light Stress")
            elif "red light" in treatment_str:
                found_treatments.add("Red Light Stress")
            elif "light stress" in treatment_str:
                found_treatments.add("Other Light Stress")
            elif "stress" in treatment_str:
                found_treatments.add("Other stress")

    # Check other fields for keywords if not already found or to add more
    if "dexamethasone" in all_text and "Chemical Stress" not in found_treatments:
        found_treatments.add("Chemical Stress")
    if "drought" in all_text and "Drought Stress" not in found_treatments:
        found_treatments.add("Drought Stress")
    if "dehydration" in all_text and "Dehidration Stress" not in found_treatments:
        found_treatments.add("Dehidration Stress")
    if ("salinity" in all_text or "salt" in all_text) and "Salinity Stress" not in found_treatments:
        found_treatments.add("Salinity Stress")
    if "heat" in all_text and "Heat Stress" not in found_treatments:
        found_treatments.add("Heat Stress")
    if "cold" in all_text and "Cold Stress" not in found_treatments:
        found_treatments.add("Cold Stress")
    if "chemical" in all_text and "Chemical Stress" not in found_treatments:
        found_treatments.add("Chemical Stress")
    if "nutrient deficiency" in all_text and "Nutrient Deficiency" not in found_treatments:
        found_treatments.add("Nutrient Deficiency")
    if "biotic stress" in all_text and "Biotic Stress" not in found_treatments:
        found_treatments.add("Biotic Stress")
    if "low light" in all_text and "Low Light Stress" not in found_treatments:
        found_treatments.add("Low Light Stress")
    if "high light" in all_text and "High Light Stress" not in found_treatments:
        found_treatments.add("High Light Stress")
    if "red light" in all_text and "Red Light Stress" not in found_treatments:
        found_treatments.add("Red Light Stress")
    if "light stress" in all_text and not any(s in found_treatments for s in ["Low Light Stress", "High Light Stress", "Red Light Stress", "Other Light Stress"]):
        found_treatments.add("Other Light Stress")
    if "stress" in all_text and not found_treatments:
        found_treatments.add("Other stress")

    if not found_treatments:
        extracted_data["treatment"] = ["No stress"]
    else:
        extracted_data["treatment"] = sorted(list(found_treatments))

    # --- Extract Medium ---
    if "ms plates" in all_text or "liquid ms media" in all_text or "ms medium" in all_text:
        extracted_data["medium"] = "MS medium"
    elif "agar" in all_text:
        extracted_data["medium"] = "agar"
    elif "soil" in all_text:
        extracted_data["medium"] = "soil"
    elif "hydroponic" in all_text:
        extracted_data["medium"] = "hydroponic"

    return extracted_data

import re

def GSE14961_extractor(sample_metadata: dict) -> dict:
    result = {
        "tissue": "unknown",
        "treatment": [],
        "medium": "unspecified"
    }

    # Helper to safely get text from metadata fields
    def get_text(key):
        value = sample_metadata.get(key)
        if isinstance(value, list):
            return " ".join(value).lower()
        elif isinstance(value, str):
            return value.lower()
        return ""

    # Combine relevant text fields for easier searching for tissue
    all_text_for_tissue = (
        get_text('source_name_ch1') + " " +
        get_text('characteristics_ch1') + " " +
        get_text('title') + " " +
        get_text('description')
    )

    # --- Extract Tissue ---
    # Prioritize more specific terms
    if re.search(r'\bseedling(s)?\b', all_text_for_tissue):
        result["tissue"] = "seedling"
    elif re.search(r'\bwhole plant(s)?\b', all_text_for_tissue):
        result["tissue"] = "whole_plant"
    elif re.search(r'\b(root|roots)\b', all_text_for_tissue):
        result["tissue"] = "root"
    elif re.search(r'\b(leaf|leaves)\b', all_text_for_tissue):
        result["tissue"] = "leaf"
    elif re.search(r'\bflower(s)?\b', all_text_for_tissue):
        result["tissue"] = "flower"
    elif re.search(r'\bshoot(s)?\b', all_text_for_tissue):
        result["tissue"] = "shoot"
    elif re.search(r'\brosette(s)?\b', all_text_for_tissue):
        result["tissue"] = "rosette"
    elif re.search(r'\bbud(s)?\b', all_text_for_tissue):
        result["tissue"] = "bud"
    elif re.search(r'\bsilique(s)?\b', all_text_for_tissue):
        result["tissue"] = "silique"
    elif re.search(r'\bcallus\b', all_text_for_tissue):
        result["tissue"] = "callus"
    elif re.search(r'\bseed(s)?\b', all_text_for_tissue):
        result["tissue"] = "seed"
    
    # --- Extract Treatment ---
    treatments = set()
    treatment_text = (
        get_text('title') + " " +
        get_text('treatment_protocol_ch1') + " " +
        get_text('description')
    )

    # Chemical Stress
    if re.search(r'\b(sa|salicylic acid|aba|abscisic acid|chemical|herbicide|pesticide|hormone|auxin|cytokinin|gibberellin|ethylene|jasmonate)\b', treatment_text):
        treatments.add("Chemical Stress")
    
    # Drought/Dehydration Stress
    if re.search(r'\b(drought|water deficit|water deprivation)\b', treatment_text):
        treatments.add("Drought Stress")
    if re.search(r'\b(dehydration)\b', treatment_text): # Specific check for dehydration
        treatments.add("Dehidration Stress")
    
    # Salinity Stress
    if re.search(r'\b(salinity|nacl|salt)\b', treatment_text):
        treatments.add("Salinity Stress")
    
    # Heat Stress
    if re.search(r'\b(heat|high temperature)\b', treatment_text):
        treatments.add("Heat Stress")
    
    # Cold Stress
    if re.search(r'\b(cold|low temperature|chilling)\b', treatment_text):
        treatments.add("Cold Stress")
    
    # Nutrient Deficiency
    if re.search(r'\b(nutrient deficiency|nitrogen starvation|phosphate starvation|potassium starvation|sulfur starvation)\b', treatment_text):
        treatments.add("Nutrient Deficiency")
    
    # Biotic Stress
    if re.search(r'\b(biotic|pathogen|insect|fungus|bacteria|virus|infection)\b', treatment_text):
        treatments.add("Biotic Stress")
    
    # Light Stress (order matters for specific vs. general)
    if re.search(r'\b(red light)\b', treatment_text):
        treatments.add("Red Light Stress")
    elif re.search(r'\b(low light|shade)\b', treatment_text):
        treatments.add("Low Light Stress")
    elif re.search(r'\b(high light|excess light)\b', treatment_text):
        treatments.add("High Light Stress")
    elif re.search(r'\b(light stress|uv light)\b', treatment_text):
        treatments.add("Other Light Stress")

    # If no specific stress is found, and "mock/control" is present, add "No stress".
    # If still no treatments, default to "No stress".
    if not treatments and re.search(r'\b(mock|control|untreated)\b', treatment_text):
        treatments.add("No stress")
    elif not treatments: # If no specific stress and no mock/control, assume no stress
        treatments.add("No stress")

    result["treatment"] = sorted(list(treatments))

    # --- Extract Medium ---
    medium_text = (
        get_text('growth_protocol_ch1') + " " +
        get_text('treatment_protocol_ch1')
    )
    
    if re.search(r'\b(ms medium|murashige and skoog|agar|gelrite)\b', medium_text):
        result["medium"] = "agar"
    elif re.search(r'\b(water|liquid medium|hydroponic)\b', medium_text):
        result["medium"] = "water"
    elif re.search(r'\b(soil|potting mix)\b', medium_text):
        result["medium"] = "soil"
    # Default is "unspecified" if none of the above are found.

    return result

def GSE44405_extractor(sample_metadata: dict) -> dict:
    result = {}

    # Helper to get text from metadata fields, handling lists and missing keys
    def get_text(keys):
        text_parts = []
        for key in keys:
            value = sample_metadata.get(key)
            if value:
                if isinstance(value, list):
                    text_parts.extend(value)
                else:
                    text_parts.append(str(value))
        return " ".join(text_parts).lower()

    # --- Extract Tissue ---
    tissue_text = get_text(['characteristics_ch1', 'source_name_ch1', 'title', 'description'])
    extracted_tissue = "unknown"

    tissue_mapping = {
        "pistil": "flower",
        "flower": "flower",
        "inflorescence": "flower",
        "petal": "flower",
        "sepal": "flower",
        "stamen": "flower",
        "carpel": "flower",
        "ovary": "flower",
        "anther": "flower",
        "root": "root",
        "leaf": "leaf",
        "shoot": "shoot",
        "rosette": "rosette",
        "bud": "bud",
        "whole plant": "whole_plant",
        "whole_plants": "whole_plant",
        "silique": "silique",
        "callus": "callus",
        "seed": "seed",
        "embryo": "seed",
        "seedling": "seedling",
        "cotyledon": "seedling",
        "hypocotyl": "seedling",
        "stem": "shoot", 
        "apex": "shoot", 
        "meristem": "shoot", 
    }

    # Prioritize "tissue: X" pattern
    import re
    match = re.search(r'tissue:\s*([\w\s-]+)', tissue_text)
    if match:
        found_tissue_val = match.group(1).strip()
        for k, v in tissue_mapping.items():
            if k in found_tissue_val: 
                extracted_tissue = v
                break
    
    # If not found via "tissue: X" or if the found value didn't map, search keywords directly in all text
    if extracted_tissue == "unknown":
        for k, v in tissue_mapping.items():
            if k in tissue_text:
                extracted_tissue = v
                break
    
    # Final check against enum values
    valid_tissues = ["root", "leaf", "flower", "shoot", "rosette", "bud", "whole_plant", "silique", "callus", "seed", "seedling", "unknown"]
    if extracted_tissue not in valid_tissues:
        extracted_tissue = "unknown"

    result['tissue'] = extracted_tissue

    # --- Extract Treatment ---
    treatment_text = get_text(['characteristics_ch1', 'treatment_protocol_ch1', 'growth_protocol_ch1', 'title', 'description'])
    extracted_treatments = set()

    treatment_keywords = {
        "drought stress": "Drought Stress",
        "drought": "Drought Stress",
        "water deficit": "Drought Stress",
        "desiccation": "Dehidration Stress",
        "dehydration stress": "Dehidration Stress",
        "salinity stress": "Salinity Stress",
        "salt stress": "Salinity Stress",
        "nacl": "Salinity Stress",
        "heat stress": "Heat Stress",
        "high temperature": "Heat Stress",
        "cold stress": "Cold Stress",
        "low temperature": "Cold Stress",
        "chemical stress": "Chemical Stress",
        "herbicide": "Chemical Stress",
        "pesticide": "Chemical Stress",
        "heavy metal": "Chemical Stress",
        "hormone treatment": "Chemical Stress", 
        "nutrient deficiency": "Nutrient Deficiency",
        "nitrogen starvation": "Nutrient Deficiency",
        "phosphate starvation": "Nutrient Deficiency",
        "low nitrogen": "Nutrient Deficiency",
        "low phosphate": "Nutrient Deficiency",
        "biotic stress": "Biotic Stress",
        "pathogen": "Biotic Stress",
        "insect": "Biotic Stress",
        "fungus": "Biotic Stress",
        "bacteria": "Biotic Stress",
        "virus": "Biotic Stress",
        "infection": "Biotic Stress",
        "low light stress": "Low Light Stress",
        "shade": "Low Light Stress",
        "high light stress": "High Light Stress",
        "excess light": "High Light Stress",
        "red light stress": "Red Light Stress",
        "red light": "Red Light Stress",
        "blue light": "Other Light Stress",
        "uv light": "Other Light Stress",
        "light stress": "Other Light Stress", 
        "pollinated": "Other stress", 
        "wounding": "Other stress",
        "mechanical stress": "Other stress",
        "gravity": "Other stress",
        "osmotic stress": "Other stress", 
    }

    for keyword, treatment_type in treatment_keywords.items():
        if keyword in treatment_text:
            extracted_treatments.add(treatment_type)

    # Check for "no stress" or control conditions
    if not extracted_treatments:
        if "control" in treatment_text or "untreated" in treatment_text or "normal growth" in treatment_text:
            extracted_treatments.add("No stress")
        elif "pollinated" in treatment_text: 
            extracted_treatments.add("Other stress")
        else:
            extracted_treatments.add("No stress") 

    result['treatment'] = sorted(list(extracted_treatments)) 

    # --- Extract Medium ---
    medium_text = get_text(['growth_protocol_ch1', 'characteristics_ch1', 'source_name_ch1', 'description'])
    extracted_medium = "unspecified"

    medium_keywords = {
        "soil": "soil",
        "potting mix": "soil",
        "vermiculite": "soil",
        "rockwool": "soil",
        "agar": "agar",
        "ms medium": "agar",
        "murashige and skoog": "agar",
        "gel": "agar",
        "hydroponic": "hydroponic",
        "liquid culture": "hydroponic",
        "nutrient solution": "hydroponic",
    }

    for keyword, medium_type in medium_keywords.items():
        if keyword in medium_text:
            extracted_medium = medium_type
            break
    
    # Infer soil if whole plant is mentioned and no other medium is specified
    if extracted_medium == "unspecified" and result['tissue'] == "whole_plant":
        extracted_medium = "soil"

    result['medium'] = extracted_medium

    return result

def GSE22836_extractor(sample_metadata: dict) -> dict:
    extracted_data = {}

    # Helper to safely get and join list values, converting to lower case
    def get_combined_text(keys):
        text_parts = []
        for key in keys:
            if key in sample_metadata and isinstance(sample_metadata[key], list):
                text_parts.extend(sample_metadata[key])
        return " ".join(text_parts).lower()

    # --- Extract Tissue ---
    tissue_raw = None
    if 'characteristics_ch1' in sample_metadata and isinstance(sample_metadata['characteristics_ch1'], list):
        for char in sample_metadata['characteristics_ch1']:
            if char.lower().startswith('tissue:'):
                tissue_raw = char.split(':', 1)[1].strip().lower()
                break

    tissue_mapping = {
        'rosette leaves': 'rosette',
        'whole rosettes': 'rosette',
        'rosette': 'rosette',
        'whole plant': 'whole_plant',
        'seedling': 'seedling',
        'root': 'root',
        'leaf': 'leaf',
        'leaves': 'leaf',
        'flower': 'flower',
        'shoot': 'shoot',
        'bud': 'bud',
        'silique': 'silique',
        'callus': 'callus',
        'seed': 'seed',
    }

    extracted_tissue = 'unknown'
    if tissue_raw:
        for key, value in tissue_mapping.items():
            if key in tissue_raw: # Use 'in' for partial matches like "rosette leaves"
                extracted_tissue = value
                break
    
    # Fallback to general text search if not found in characteristics
    if extracted_tissue == 'unknown':
        combined_general_text = get_combined_text(['title', 'source_name_ch1', 'description'])
        for key, value in tissue_mapping.items():
            if key in combined_general_text:
                extracted_tissue = value
                break
        
    extracted_data['tissue'] = extracted_tissue

    # --- Extract Treatment ---
    treatments = set()
    combined_treatment_text = get_combined_text(['title', 'source_name_ch1', 'treatment_protocol_ch1', 'description'])

    treatment_keywords = {
        "Drought Stress": ["drought", "dehydration", "water deficit", "osmotic stress", "peg treatment"],
        "Salinity Stress": ["salinity", "salt stress", "nacl"],
        "Heat Stress": ["heat stress", "high temperature", "heat shock"],
        "Cold Stress": ["cold stress", "low temperature", "chilling"],
        "Chemical Stress": ["dex-induced", "chemical", "herbicide", "pesticide", "heavy metal", "cadmium", "aba", "methyl jasmonate", "meja", "sa", "salicylic acid", "ethephon", "acc", "ga", "gibberellin", "cytokinin", "auxin", "brassinosteroid", "strigolactone", "ethylene", "jasmonate", "oxidative stress", "h2o2", "paraquat", "ozone", "hormone", "drug", "chemical treatment"],
        "Nutrient Deficiency": ["nutrient deficiency", "nitrogen starvation", "phosphate starvation", "iron deficiency", "sulfur deficiency", "boron deficiency", "nutrient deprivation"],
        "Biotic Stress": ["pathogen", "fungus", "bacteria", "virus", "insect", "herbivory", "infection", "elicitor", "pest"],
        "Low Light Stress": ["low light", "shade"],
        "High Light Stress": ["high light", "excess light"],
        "Red Light Stress": ["red light"],
        "Other Light Stress": ["light stress", "uv-b", "uv-a", "blue light", "far-red light", "darkness", "dark treatment"],
        "Other stress": ["stress", "wound", "mechanical stress", "gravity", "anoxia", "hypoxia", "osmotic stress (non-drought)"],
    }

    for treatment_type, keywords in treatment_keywords.items():
        for keyword in keywords:
            if keyword in combined_treatment_text:
                treatments.add(treatment_type)
                break

    is_control = False
    if "control" in combined_treatment_text or "mock" in combined_treatment_text or "untreated" in combined_treatment_text:
        is_control = True
    
    if not treatments:
        if is_control:
            treatments.add("No stress")
        else:
            # If no specific stress keywords and not explicitly a control, assume "No stress"
            treatments.add("No stress")
    
    # If "No stress" is present along with other stresses, remove it.
    # "No stress" should only be present if it's the *only* identified treatment.
    if len(treatments) > 1 and "No stress" in treatments:
        treatments.remove("No stress")

    extracted_data['treatment'] = sorted(list(treatments))

    # --- Extract Medium ---
    medium_text = get_combined_text(['growth_protocol_ch1'])

    extracted_medium = 'unspecified'
    if medium_text:
        if 'soil' in medium_text:
            extracted_medium = 'soil'
        elif 'ms medium' in medium_text or 'murashige skoog' in medium_text:
            extracted_medium = 'MS medium'
        elif 'agar' in medium_text or 'gel' in medium_text:
            extracted_medium = 'agar'
        elif 'liquid medium' in medium_text or 'liquid culture' in medium_text or 'hydroponic' in medium_text:
            extracted_medium = 'liquid medium'
        elif 'vermiculite' in medium_text:
            extracted_medium = 'vermiculite'
        elif 'rockwool' in medium_text:
            extracted_medium = 'rockwool'
        elif 'sand' in medium_text:
            extracted_medium = 'sand'
        elif 'peat' in medium_text:
            extracted_medium = 'peat'
        elif 'perlite' in medium_text:
            extracted_medium = 'perlite'
        # If whole plant is mentioned and no medium, assume soil
        elif extracted_data['tissue'] == 'whole_plant' and extracted_medium == 'unspecified':
            extracted_medium = 'soil'
    
    extracted_data['medium'] = extracted_medium

    return extracted_data

def GSE9996_extractor(sample_metadata: dict) -> dict:
    """
    Extracts tissue, treatment, and medium information from a sample_metadata dictionary
    following a predefined JSON schema.

    Args:
        sample_metadata (dict): A dictionary containing sample metadata, typically
                                from a GEO series entry.

    Returns:
        dict: A dictionary conforming to the specified schema for tissue, treatment, and medium.
    """

    def get_combined_text(metadata: dict, keys: list) -> str:
        """Combines text from specified metadata keys into a single lowercase string."""
        text_parts = []
        for key in keys:
            if key in metadata and isinstance(metadata[key], list):
                text_parts.extend(metadata[key])
        return " ".join(text_parts).lower()

    def extract_tissue(combined_text: str) -> str:
        """Extracts tissue type based on keywords in the combined text."""
        tissue_keywords = {
            "root": ["root"],
            "leaf": ["leaf"],
            "flower": ["flower"],
            "shoot": ["shoot"],
            "rosette": ["rosette"],
            "bud": ["bud"],
            "whole_plant": ["whole plant", "entire plant"],
            "silique": ["silique"],
            "callus": ["callus"],
            "seed": ["seed"],
            "seedling": ["seedling"]
        }

        # Prioritize more specific tissues
        priority_order = [
            "root", "leaf", "flower", "bud", "silique", "callus", "seed", "seedling",
            "rosette", "shoot", "whole_plant"
        ]

        for tissue_type in priority_order:
            for keyword in tissue_keywords[tissue_type]:
                if keyword in combined_text:
                    return tissue_type

        return "unknown"

    def extract_treatment(combined_text: str) -> list:
        """Extracts a list of treatments based on keywords in the combined text."""
        found_treatments = set()

        treatment_keywords = {
            "Drought Stress": ["drought", "water deficit", "water deprivation"],
            "Dehidration Stress": ["dehydration", "desiccation"],
            "Salinity Stress": ["salinity", "salt stress", "nacl"],
            "Heat Stress": ["heat stress", "high temperature"],
            "Cold Stress": ["cold stress", "low temperature", "chilling"],
            "Chemical Stress": ["chemical", "herbicide", "pesticide", "metal", "heavy metal", "toxic", "drug", "hormone"],
            "Nutrient Deficiency": ["nutrient deficiency", "starvation", "low nitrogen", "low phosphate", "low potassium"],
            "Biotic Stress": ["biotic", "pathogen", "fungus", "bacteria", "virus", "insect", "herbivory"],
            "Low Light Stress": ["low light", "shade"],
            "High Light Stress": ["high light", "excess light"],
            "Red Light Stress": ["red light"],
            "Other Light Stress": ["uv light", "blue light", "far-red light", "light stress"],
            "Other stress": ["stress", "wounding", "mechanical stress", "physical stress", "surgical excision"]
        }

        for treatment_type, keywords in treatment_keywords.items():
            for keyword in keywords:
                if keyword in combined_text:
                    found_treatments.add(treatment_type)
                    break # Move to next treatment_type once one keyword is found

        if not found_treatments:
            return ["No stress"]

        # Remove "No stress" if other specific stresses are found
        if len(found_treatments) > 1 and "No stress" in found_treatments:
            found_treatments.remove("No stress")

        return sorted(list(found_treatments))

    def extract_medium(combined_text: str, tissue: str) -> str:
        """Extracts growth medium information based on keywords and inferred tissue."""
        medium_parts = set()

        if "agar" in combined_text or "gel" in combined_text:
            medium_parts.add("agar")
        if "ms salts" in combined_text or "murashige skoog" in combined_text:
            medium_parts.add("MS salts")
        if "sucrose" in combined_text:
            medium_parts.add("sucrose")
        if "soil" in combined_text or "potting mix" in combined_text:
            medium_parts.add("soil")
        if "vermiculite" in combined_text:
            medium_parts.add("vermiculite")
        if "hydroponic" in combined_text:
            medium_parts.add("hydroponic solution")
        if "water" in combined_text and not medium_parts:
            medium_parts.add("water")

        # General medium types, added only if more specific components aren't already found
        if "liquid medium" in combined_text and not any(m in medium_parts for m in ["hydroponic solution", "water"]):
            medium_parts.add("liquid medium")
        if "solid medium" in combined_text and not any(m in medium_parts for m in ["agar", "soil", "vermiculite"]):
            medium_parts.add("solid medium")

        if medium_parts:
            return ", ".join(sorted(list(medium_parts)))

        # Infer based on tissue if no medium found
        if tissue == "whole_plant":
            return "soil"

        return "unspecified"

    # Combine relevant text fields for comprehensive searching
    combined_text = get_combined_text(
        sample_metadata,
        ['title', 'source_name_ch1', 'characteristics_ch1', 'treatment_protocol_ch1', 'growth_protocol_ch1', 'description']
    )

    # Extract information
    extracted_tissue = extract_tissue(combined_text)
    extracted_treatment = extract_treatment(combined_text)
    extracted_medium = extract_medium(combined_text, extracted_tissue)

    return {
        "tissue": extracted_tissue,
        "treatment": extracted_treatment,
        "medium": extracted_medium
    }

def GSE28109_extractor(sample_metadata: dict) -> dict:
    """
    Extracts tissue, treatment, and medium information from a sample_metadata dictionary
    following a predefined schema.

    Args:
        sample_metadata (dict): A dictionary containing sample metadata, typically
                                from a GEO series entry.

    Returns:
        dict: A dictionary with 'tissue', 'treatment', and 'medium' fields
              conforming to the specified JSON schema.
    """

    def _get_all_text(data, keys):
        """Concatenate all text from specified keys, handling lists and converting to lowercase."""
        all_text = []
        for key in keys:
            value = data.get(key)
            if isinstance(value, list):
                all_text.extend(value)
            elif isinstance(value, str):
                all_text.append(value)
        return " ".join(all_text).lower()

    result = {}

    # --- Tissue Extraction ---
    tissue_found = "unknown"
    search_text = _get_all_text(sample_metadata, [
        'source_name_ch1', 'characteristics_ch1', 'title', 'description', 'growth_protocol_ch1'
    ])

    # Prioritize more specific terms or common ones based on typical biological context
    if "whole plant" in search_text or "plants were grown" in search_text:
        tissue_found = "whole_plant"
    if "shoot apices" in search_text or "shoot" in search_text:
        tissue_found = "shoot"
    if "leaf" in search_text:
        tissue_found = "leaf"
    if "root" in search_text:
        tissue_found = "root"
    if "flower" in search_text:
        tissue_found = "flower"
    if "rosette" in search_text:
        tissue_found = "rosette"
    if "bud" in search_text:
        tissue_found = "bud"
    if "silique" in search_text:
        tissue_found = "silique"
    if "callus" in search_text:
        tissue_found = "callus"
    if "seedling" in search_text:
        tissue_found = "seedling"
    if "seed" in search_text:
        tissue_found = "seed"
    
    # Refine based on specific phrases that might override general terms
    # Example: "cell type: Subepidermal" combined with "shoot apices" implies shoot
    if "cell type" in search_text and "subepidermal" in search_text and "shoot apices" in search_text:
        tissue_found = "shoot"

    result["tissue"] = tissue_found

    # --- Treatment Extraction ---
    treatments = set()
    search_text_treatment = _get_all_text(sample_metadata, [
        'treatment_protocol_ch1', 'characteristics_ch1', 'description', 'growth_protocol_ch1'
    ])

    # Map common phrases to schema enum values
    if "drought" in search_text_treatment:
        treatments.add("Drought Stress")
    if "dehydration" in search_text_treatment:
        treatments.add("Dehidration Stress")
    if "salinity" in search_text_treatment or "salt stress" in search_text_treatment:
        treatments.add("Salinity Stress")
    if "heat stress" in search_text_treatment or "high temperature" in search_text_treatment:
        treatments.add("Heat Stress")
    if "cold stress" in search_text_treatment or "low temperature" in search_text_treatment:
        treatments.add("Cold Stress")
    if "chemical" in search_text_treatment or "herbicide" in search_text_treatment or "pesticide" in search_text_treatment:
        treatments.add("Chemical Stress")
    if "nutrient deficiency" in search_text_treatment or "starvation" in search_text_treatment:
        treatments.add("Nutrient Deficiency")
    if "biotic stress" in search_text_treatment or "pathogen" in search_text_treatment or "insect" in search_text_treatment:
        treatments.add("Biotic Stress")
    if "low light" in search_text_treatment:
        treatments.add("Low Light Stress")
    if "high light" in search_text_treatment:
        treatments.add("High Light Stress")
    if "red light" in search_text_treatment:
        treatments.add("Red Light Stress")
    
    # General light stress if not specific
    if "light stress" in search_text_treatment and not any(t in treatments for t in ["Low Light Stress", "High Light Stress", "Red Light Stress"]):
        treatments.add("Other Light Stress")
    
    # General stress if not specific and no other specific stress found
    if "stress" in search_text_treatment and not treatments:
        treatments.add("Other stress")

    if not treatments:
        treatments.add("No stress")
    
    result["treatment"] = sorted(list(treatments))

    # --- Medium Extraction ---
    medium_found = "unspecified"
    search_text_medium = _get_all_text(sample_metadata, [
        'growth_protocol_ch1', 'characteristics_ch1', 'description'
    ])

    if "soil" in search_text_medium:
        medium_found = "soil"
    elif "agar" in search_text_medium:
        medium_found = "agar"
    elif "ms medium" in search_text_medium or "murashige skoog" in search_text_medium:
        medium_found = "MS medium"
    elif "liquid medium" in search_text_medium or "liquid culture" in search_text_medium:
        medium_found = "liquid medium"
    elif "hydroponics" in search_text_medium:
        medium_found = "hydroponics"
    elif "vermiculite" in search_text_medium:
        medium_found = "vermiculite"
    elif "peat" in search_text_medium:
        medium_found = "peat"
    elif "rockwool" in search_text_medium:
        medium_found = "rockwool"
    elif "sand" in search_text_medium:
        medium_found = "sand"
    elif "perlite" in search_text_medium:
        medium_found = "perlite"
    elif "gel" in search_text_medium:
        medium_found = "gel"
    elif "water" in search_text_medium:
        medium_found = "water"
    
    # Inference for medium if not explicitly found
    if medium_found == "unspecified":
        # If the tissue implies a whole plant or a part of it grown in natural conditions, infer soil.
        if result["tissue"] in ["whole_plant", "shoot", "leaf", "flower", "rosette", "bud", "silique"]:
            if "plants were grown" in search_text_medium: # Strong indicator for soil if no other medium
                medium_found = "soil"
        # For other tissues like seedling or callus, "unspecified" is a safer default if no explicit medium.
        # If more specific inference is needed, additional rules would be added here.
    
    result["medium"] = medium_found

    return result

def GSE9728_extractor(sample_metadata: dict) -> dict:
    """
    Extracts tissue, treatment, and medium information from a sample_metadata dictionary
    following a predefined schema.

    Args:
        sample_metadata (dict): A dictionary containing sample metadata, typically
                                from a GEO series entry.

    Returns:
        dict: A dictionary containing the extracted 'tissue', 'treatment', and 'medium'
              information, conforming to the specified JSON schema.
    """

    def _get_text_from_metadata(metadata: dict, keys: list) -> str:
        """Helper to extract and concatenate text from specified metadata keys."""
        texts = []
        for key in keys:
            if key in metadata:
                value = metadata[key]
                if isinstance(value, list):
                    texts.extend(value)
                elif isinstance(value, str):
                    texts.append(value)
        return " ".join(texts).lower()

    extracted_data = {}

    # --- Extract Tissue ---
    tissue = "unknown"
    search_text = _get_text_from_metadata(sample_metadata, ['source_name_ch1', 'characteristics_ch1', 'title', 'description'])

    # Prioritize more specific terms or common terms for this dataset
    if "seedling" in search_text:
        tissue = "seedling"
    elif "whole plant" in search_text:
        tissue = "whole_plant"
    elif "rosette" in search_text:
        tissue = "rosette"
    elif "leaf" in search_text:
        tissue = "leaf"
    elif "root" in search_text:
        tissue = "root"
    elif "flower" in search_text:
        tissue = "flower"
    elif "shoot" in search_text:
        tissue = "shoot"
    elif "bud" in search_text:
        tissue = "bud"
    elif "silique" in search_text:
        tissue = "silique"
    elif "callus" in search_text:
        tissue = "callus"
    elif "seed" in search_text:
        tissue = "seed"
    
    extracted_data["tissue"] = tissue

    # --- Extract Treatment ---
    found_treatments = set()
    search_text_treatment = _get_text_from_metadata(sample_metadata, ['characteristics_ch1', 'title', 'description'])

    if "drought" in search_text_treatment:
        found_treatments.add("Drought Stress")
    if "dehydration" in search_text_treatment:
        found_treatments.add("Dehidration Stress")
    if "salinity" in search_text_treatment or "salt stress" in search_text_treatment:
        found_treatments.add("Salinity Stress")
    if "heat stress" in search_text_treatment:
        found_treatments.add("Heat Stress")
    if "cold stress" in search_text_treatment:
        found_treatments.add("Cold Stress")
    if "chemical" in search_text_treatment:
        found_treatments.add("Chemical Stress")
    if "nutrient deficiency" in search_text_treatment:
        found_treatments.add("Nutient Deficiency") # Typo in schema: "Nutrient Deficiency"
    if "biotic" in search_text_treatment:
        found_treatments.add("Biotic Stress")
    if "low light" in search_text_treatment:
        found_treatments.add("Low Light Stress")
    if "high light" in search_text_treatment:
        found_treatments.add("High Light Stress")
    if "red light" in search_text_treatment:
        found_treatments.add("Red Light Stress")
    if "light stress" in search_text_treatment and not any(s in search_text_treatment for s in ["low light", "high light", "red light"]):
        found_treatments.add("Other Light Stress")
    
    # Handle "No stress" or control conditions
    if not found_treatments:
        if any(s in search_text_treatment for s in ["no stress", "control", "untreated", "normal growth", "light grown"]):
            found_treatments.add("No stress")
        else:
            # Default to "No stress" if no specific stress or control indicator is found
            found_treatments.add("No stress")
    
    # If "No stress" is found along with other specific stresses, remove "No stress"
    if len(found_treatments) > 1 and "No stress" in found_treatments:
        found_treatments.remove("No stress")

    extracted_data["treatment"] = sorted(list(found_treatments))

    # --- Extract Medium ---
    medium = "unspecified"
    search_text_medium = _get_text_from_metadata(sample_metadata, ['characteristics_ch1', 'description', 'extract_protocol_ch1'])

    if "soil" in search_text_medium:
        medium = "soil"
    elif "agar" in search_text_medium:
        medium = "agar"
    elif "ms medium" in search_text_medium:
        medium = "MS medium"
    elif "liquid medium" in search_text_medium:
        medium = "liquid medium"
    elif "hydroponic" in search_text_medium:
        medium = "hydroponic"
    
    extracted_data["medium"] = medium

    return extracted_data

def GSE5615_extractor(sample_metadata: dict) -> dict:
    extracted_data = {}

    # Helper to safely get the first item from a list or None
    def _get_first_item(key, default=None):
        val = sample_metadata.get(key)
        if isinstance(val, list) and val:
            return val[0]
        return default

    # Helper to search for a prefix in a list of strings and return the suffix
    def _search_list_for_prefix(data_list, prefix):
        if isinstance(data_list, list):
            for item in data_list:
                if item.lower().startswith(prefix.lower()):
                    return item[len(prefix):].strip()
        return None

    # --- Extract Tissue ---
    tissue_raw = None
    characteristics_ch1 = sample_metadata.get('characteristics_ch1', [])
    tissue_raw = _search_list_for_prefix(characteristics_ch1, 'Tissue:')

    if tissue_raw:
        tissue_raw_lower = tissue_raw.lower()
        if 'leaves' in tissue_raw_lower or 'leaf' in tissue_raw_lower:
            extracted_data['tissue'] = 'leaf'
        elif 'root' in tissue_raw_lower:
            extracted_data['tissue'] = 'root'
        elif 'flower' in tissue_raw_lower:
            extracted_data['tissue'] = 'flower'
        elif 'shoot' in tissue_raw_lower:
            extracted_data['tissue'] = 'shoot'
        elif 'rosette' in tissue_raw_lower:
            extracted_data['tissue'] = 'rosette'
        elif 'bud' in tissue_raw_lower:
            extracted_data['tissue'] = 'bud'
        elif 'whole plant' in tissue_raw_lower or 'whole_plant' in tissue_raw_lower:
            extracted_data['tissue'] = 'whole_plant'
        elif 'silique' in tissue_raw_lower:
            extracted_data['tissue'] = 'silique'
        elif 'callus' in tissue_raw_lower:
            extracted_data['tissue'] = 'callus'
        elif 'seedling' in tissue_raw_lower:
            extracted_data['tissue'] = 'seedling'
        elif 'seed' in tissue_raw_lower:
            extracted_data['tissue'] = 'seed'
        else:
            extracted_data['tissue'] = 'unknown'
    else:
        extracted_data['tissue'] = 'unknown' # Default if not found

    # --- Extract Treatment ---
    detected_treatments = set()
    treatment_texts = []
    treatment_texts.append(_get_first_item('treatment_protocol_ch1', ''))
    description_list = sample_metadata.get('description', [])
    treatment_texts.extend(description_list)
    
    # Combine all relevant text for treatment detection
    full_text = " ".join(treatment_texts).lower()

    # Keywords for each treatment type
    treatment_keywords = {
        "Drought Stress": ["drought", "water deprivation", "water deficit", "dry", "osmotic stress"],
        "Dehidration Stress": ["dehydration", "desiccation"],
        "Salinity Stress": ["salinity", "salt", "nacl"],
        "Heat Stress": ["heat", "high temperature", "hot", "thermal stress"],
        "Cold Stress": ["cold", "low temperature", "chilling", "freezing"],
        "Chemical Stress": ["chemical", "herbicide", "pesticide", "fungicide", "metal", "heavy metal", "toxic", "gst", "hormone", "aba", "auxin", "cytokinin", "gibberellin", "ethylene", "jasmonate", "salicylic acid", "drug", "inhibitor", "stressor", "cadmium", "aluminum", "arsenic", "lead", "mercury", "pce", "tce", "atrazine", "glyphosate", "paraquat", "methyl viologen", "oxidative stress"],
        "Nutrient Deficiency": ["nutrient deficiency", "starvation", "low nitrogen", "low phosphate", "low potassium", "minus n", "minus p", "minus k", "nitrogen deprivation", "phosphate deprivation", "potassium deprivation"],
        "Biotic Stress": ["pathogen", "fungus", "bacteria", "virus", "insect", "herbivore", "infection", "elicitor", "pest", "microbe", "disease", "wounding", "mechanical damage"],
        "Low Light Stress": ["low light", "darkness", "shade", "etiolation"],
        "High Light Stress": ["high light", "excess light", "uv", "uv-b", "uv-c"],
        "Red Light Stress": ["red light"],
        "Other Light Stress": ["blue light", "far-red light", "light quality", "light spectrum"],
        "Other stress": ["stress", "wound", "mechanical stress", "abiotic stress", "biotic stress"] # Catch-all for general stress if specific type not found
    }

    for treatment_type, keywords in treatment_keywords.items():
        for keyword in keywords:
            if keyword in full_text:
                # Special handling for "Other stress" to avoid over-triggering
                # If a more specific stress is found, don't add "Other stress"
                if treatment_type == "Other stress" and any(t in detected_treatments for t in treatment_keywords.keys() if t != "Other stress"):
                    continue 
                detected_treatments.add(treatment_type)
                break # Move to the next treatment type once a keyword is found

    if not detected_treatments:
        extracted_data['treatment'] = ["No stress"]
    else:
        # Remove "Other stress" if more specific stresses are present
        if len(detected_treatments) > 1 and "Other stress" in detected_treatments:
            detected_treatments.remove("Other stress")
        extracted_data['treatment'] = sorted(list(detected_treatments)) # Sort for consistent output

    # --- Extract Medium ---
    medium_raw = _search_list_for_prefix(description_list, 'Medium:')
    if medium_raw:
        extracted_data['medium'] = medium_raw.strip()
    else:
        # Infer "soil" if whole_plant or root, otherwise "unspecified"
        if extracted_data['tissue'] in ['whole_plant', 'root']:
            extracted_data['medium'] = 'soil'
        else:
            extracted_data['medium'] = 'unspecified'

    return extracted_data

def GSE37130_extractor(sample_metadata: dict) -> dict:
    extracted_data = {}

    # --- Extract Tissue ---
    extracted_tissue = "unknown"
    characteristics_ch1_list = sample_metadata.get('characteristics_ch1', [])
    source_name_ch1 = sample_metadata.get('source_name_ch1', [''])[0].lower()
    title = sample_metadata.get('title', [''])[0].lower()

    tissue_mapping = {
        "root": "root", "leaf": "leaf", "flower": "flower", "shoot": "shoot",
        "rosette": "rosette", "bud": "bud", "silique": "silique", "callus": "callus",
        "seed": "seed", "seedling": "seedling", "whole plant": "whole_plant", "whole_plant": "whole_plant"
    }

    # 1. Prioritize explicit 'tissue:' in characteristics_ch1
    for char in characteristics_ch1_list:
        char_lower = char.lower()
        if 'tissue:' in char_lower:
            tissue_val = char_lower.split('tissue:')[1].strip()
            extracted_tissue = tissue_mapping.get(tissue_val, "unknown")
            if extracted_tissue != "unknown":
                break
    
    # 2. Fallback to other keywords if no explicit 'tissue:' found or mapped
    if extracted_tissue == "unknown":
        # Check for 'seedling' or 'whole plant' in characteristics_ch1
        for char in characteristics_ch1_list:
            char_lower = char.lower()
            if 'seedling' in char_lower:
                extracted_tissue = 'seedling'
                break
            elif 'whole plant' in char_lower or 'whole_plant' in char_lower:
                extracted_tissue = 'whole_plant'
                break
    
    # 3. Further fallback using source_name_ch1 and title
    if extracted_tissue == "unknown":
        if 'seedling' in source_name_ch1 or 'seedling' in title:
            extracted_tissue = 'seedling'
        elif 'plant' in source_name_ch1 or 'plant' in title:
            # If it's a general 'plant' and no specific tissue/age, default to whole_plant
            extracted_tissue = 'whole_plant'

    extracted_data['tissue'] = extracted_tissue

    # --- Extract Treatment ---
    treatments = set()
    text_fields = [
        sample_metadata.get('title', [''])[0],
        sample_metadata.get('source_name_ch1', [''])[0],
        sample_metadata.get('growth_protocol_ch1', [''])[0]
    ]
    text_fields.extend(characteristics_ch1_list)

    combined_text = " ".join(text_fields).lower()

    treatment_keywords_map = {
        "drought": "Drought Stress",
        "dehydration": "Dehidration Stress",
        "salinity": "Salinity Stress", "salt": "Salinity Stress",
        "heat": "Heat Stress",
        "cold": "Cold Stress",
        "chemical": "Chemical Stress",
        "nutrient deficiency": "Nutrient Deficiency", "low nutrient": "Nutrient Deficiency",
        "biotic": "Biotic Stress", "pathogen": "Biotic Stress", "infection": "Biotic Stress",
        "low light": "Low Light Stress",
        "high light": "High Light Stress",
        "red light": "Red Light Stress",
        "stress": "Other stress" # General stress, if not more specific
    }

    for keyword, schema_treatment in treatment_keywords_map.items():
        if keyword in combined_text:
            treatments.add(schema_treatment)

    # Check for "No stress" indicators
    no_stress_indicators = ["room temperature", "control", "no stress", "unstressed", "normal condition"]
    is_no_stress_explicit = False
    for indicator in no_stress_indicators:
        if indicator in combined_text:
            is_no_stress_explicit = True
            break

    if not treatments and is_no_stress_explicit:
        extracted_data['treatment'] = ["No stress"]
    elif treatments:
        extracted_data['treatment'] = sorted(list(treatments))
    else:
        # Default to "No stress" if no specific stress or no-stress indicator is found
        extracted_data['treatment'] = ["No stress"]

    # --- Extract Medium ---
    extracted_medium = "unspecified"
    growth_protocol_ch1 = sample_metadata.get('growth_protocol_ch1', [''])[0].lower()

    medium_keywords_map = {
        "ms plate": "MS plate", "ms medium": "MS plate", "agar plate": "agar plate",
        "soil": "soil", "potting mix": "soil",
        "hydroponic": "hydroponic solution", "hydroponics": "hydroponic solution",
        "liquid medium": "liquid medium", "liquid culture": "liquid medium"
    }

    # Check growth_protocol_ch1 first as it's most likely to contain medium info
    for keyword, schema_medium in medium_keywords_map.items():
        if keyword in growth_protocol_ch1:
            extracted_medium = schema_medium
            break
    
    # If not found, check characteristics_ch1
    if extracted_medium == "unspecified":
        for char in characteristics_ch1_list:
            char_lower = char.lower()
            for keyword, schema_medium in medium_keywords_map.items():
                if keyword in char_lower:
                    extracted_medium = schema_medium
                    break
            if extracted_medium != "unspecified":
                break

    # Special inference for "soil" if tissue is whole_plant and medium is still unspecified
    if extracted_medium == "unspecified" and extracted_data['tissue'] == "whole_plant":
        extracted_medium = "soil"
    
    extracted_data['medium'] = extracted_medium

    return extracted_data

def GSE11505_extractor(sample_metadata: dict) -> dict:
    extracted_data = {}

    # Helper to safely get and join text from a list of strings
    def get_text(field_name):
        if field_name in sample_metadata and isinstance(sample_metadata[field_name], list):
            return " ".join(sample_metadata[field_name]).lower()
        return ""

    # Combine relevant text from various fields for comprehensive searching
    all_text = ""
    for field in ['characteristics_ch1', 'source_name_ch1', 'title', 'description', 'growth_protocol_ch1', 'extract_protocol_ch1']:
        all_text += get_text(field) + " "
    all_text = all_text.strip()

    # --- Extract Tissue ---
    tissue_keywords = {
        "leaf": ["leaf", "leaves"],
        "root": ["root", "roots"],
        "flower": ["flower", "flowers"],
        "shoot": ["shoot", "shoots"],
        "rosette": ["rosette", "rosettes"],
        "bud": ["bud", "buds"],
        "silique": ["silique", "siliques"],
        "callus": ["callus"],
        "seed": ["seed", "seeds"],
        "seedling": ["seedling", "seedlings"],
        "whole_plant": ["whole plant", "plants"]
    }

    extracted_tissue = "unknown"
    found_specific_tissue = False

    # Prioritize specific tissues (e.g., leaf, root) over general terms (seedling, whole_plant)
    for tissue_type, keywords in tissue_keywords.items():
        if tissue_type in ["seedling", "whole_plant"]:
            continue # Handle these with lower priority later
        for keyword in keywords:
            if keyword in all_text:
                extracted_tissue = tissue_type
                found_specific_tissue = True
                break
        if found_specific_tissue:
            break

    # If no specific tissue found, check for seedling or whole_plant
    if not found_specific_tissue:
        if "seedling" in all_text:
            extracted_tissue = "seedling"
        elif "plant" in all_text: # Catch "plant" if no other specific tissue
            extracted_tissue = "whole_plant"
    
    # Special case: if "leaf tissues from seedlings" is mentioned, "leaf" is more specific than "seedling"
    if "leaf tissues from" in all_text and "seedlings" in all_text:
        extracted_tissue = "leaf"

    extracted_data["tissue"] = extracted_tissue

    # --- Extract Treatment ---
    treatment_keywords = {
        "Drought Stress": ["drought", "water deficit", "dehydration"],
        "Dehidration Stress": ["dehydration", "water deficit"],
        "Salinity Stress": ["salinity", "salt stress", "nacl"],
        "Heat Stress": ["heat stress", "high temperature"],
        "Cold Stress": ["cold stress", "low temperature", "chilling", "freezing"],
        "Chemical Stress": ["chemical stress", "herbicide", "pesticide", "heavy metal", "cadmium", "arsenic", "antibiotic", "selection marker"],
        "Nutrient Deficiency": ["nutrient deficiency", "nitrogen starvation", "phosphate starvation", "iron deficiency"],
        "Biotic Stress": ["biotic stress", "pathogen", "fungal", "bacterial", "viral", "insect", "herbivory"],
        "Low Light Stress": ["low light", "shade"],
        "High Light Stress": ["high light", "excess light"],
        "Red Light Stress": ["red light"],
        "Other Light Stress": ["uv light", "blue light", "far-red light"],
        "Other stress": ["stress", "abiotic stress", "oxidative stress", "osmotic stress"]
    }

    found_treatments = set()

    for treatment_type, keywords in treatment_keywords.items():
        for keyword in keywords:
            if keyword in all_text:
                found_treatments.add(treatment_type)
                break # Move to next treatment_type once one keyword is found

    if not found_treatments:
        found_treatments.add("No stress")

    extracted_data["treatment"] = sorted(list(found_treatments))

    # --- Extract Medium ---
    extracted_medium = "unspecified"

    # Prioritize specific and detailed mentions
    if "grown on ms agar plates" in all_text:
        extracted_medium = "MS agar plates"
    elif "ms media" in all_text or "ms medium" in all_text:
        extracted_medium = "MS media"
    elif "soil" in all_text:
        extracted_medium = "soil"
    elif "agar" in all_text:
        extracted_medium = "agar"
    elif "plates" in all_text: # General mention of plates, often implies agar
        extracted_medium = "agar"
    
    # Infer "soil" if tissue is "whole_plant" and medium is still unspecified
    if extracted_medium == "unspecified" and extracted_data["tissue"] == "whole_plant":
        extracted_medium = "soil"

    extracted_data["medium"] = extracted_medium

    return extracted_data

def GSE18112_extractor(sample_metadata: dict) -> dict:
    """
    Extracts tissue, treatment, and medium information from a sample_metadata dictionary
    following a predefined JSON schema.

    Args:
        sample_metadata (dict): A dictionary containing metadata for a sample.

    Returns:
        dict: A dictionary containing the extracted 'tissue', 'treatment', and 'medium'
              information, conforming to the specified schema.
    """
    result = {}

    # Helper to safely get and join text from metadata fields
    def _get_text(data_dict, key):
        value = data_dict.get(key)
        if isinstance(value, list):
            return " ".join(value).lower()
        return str(value).lower() if value is not None else ""

    # Combine relevant text fields for easier searching for tissue and treatment
    combined_text_for_tissue_and_treatment = " ".join([
        _get_text(sample_metadata, 'source_name_ch1'),
        _get_text(sample_metadata, 'characteristics_ch1'),
        _get_text(sample_metadata, 'title'),
        _get_text(sample_metadata, 'description'),
        _get_text(sample_metadata, 'treatment_protocol_ch1')
    ])
    
    # --- Extract Tissue ---
    extracted_tissue = "unknown"
    tissue_keywords = {
        "seed": ["seed", "seeds"],
        "root": ["root"],
        "leaf": ["leaf", "leaves"],
        "flower": ["flower", "flowers"],
        "shoot": ["shoot"],
        "rosette": ["rosette"],
        "bud": ["bud"],
        "whole_plant": ["whole plant", "plant"],
        "silique": ["silique"],
        "callus": ["callus"],
        "seedling": ["seedling"],
    }

    for tissue_type, keywords in tissue_keywords.items():
        for keyword in keywords:
            if keyword in combined_text_for_tissue_and_treatment:
                extracted_tissue = tissue_type
                break
        if extracted_tissue != "unknown":
            break
    result["tissue"] = extracted_tissue

    # --- Extract Treatment ---
    treatments = set()
    
    # Check for "No stress" explicitly in treatment protocol
    treatment_protocol_text = _get_text(sample_metadata, 'treatment_protocol_ch1')
    if "none" in treatment_protocol_text or "no stress" in treatment_protocol_text:
        treatments.add("No stress")
    else:
        # Search for specific stresses in the combined text
        if any(k in combined_text_for_tissue_and_treatment for k in ["drought", "dehydration", "water deficit"]):
            treatments.add("Drought Stress")
        if any(k in combined_text_for_tissue_and_treatment for k in ["salinity", "salt", "nacl"]):
            treatments.add("Salinity Stress")
        if any(k in combined_text_for_tissue_and_treatment for k in ["heat", "high temperature", "hot"]):
            treatments.add("Heat Stress")
        if any(k in combined_text_for_tissue_and_treatment for k in ["cold", "low temperature", "chilling"]):
            treatments.add("Cold Stress")
        # Chemical Stress - broad category, look for common chemical terms or general "chemical"
        if any(k in combined_text_for_tissue_and_treatment for k in ["chemical", "hormone", "pesticide", "herbicide", "aba", "auxin", "cytokinin", "gibberellin", "ethylene", "jasmonate", "salicylic acid"]):
            treatments.add("Chemical Stress")
        if any(k in combined_text_for_tissue_and_treatment for k in ["nutrient deficiency", "low nitrogen", "phosphate starvation", "nitrogen starvation", "potassium deficiency"]):
            treatments.add("Nutrient Deficiency")
        if any(k in combined_text_for_tissue_and_treatment for k in ["pathogen", "insect", "fungal", "bacterial", "virus", "herbivory", "biotic"]):
            treatments.add("Biotic Stress")
        if "low light" in combined_text_for_tissue_and_treatment:
            treatments.add("Low Light Stress")
        if "high light" in combined_text_for_tissue_and_treatment:
            treatments.add("High Light Stress")
        if "red light" in combined_text_for_tissue_and_treatment:
            treatments.add("Red Light Stress")
        # Catch-all for other light stresses, only if more specific light stresses aren't already added
        if any(k in combined_text_for_tissue_and_treatment for k in ["light stress", "darkness", "uv light", "blue light", "far-red light"]):
            if not any(s in treatments for s in ["Low Light Stress", "High Light Stress", "Red Light Stress"]):
                treatments.add("Other Light Stress")
        
        # If no specific stress found and "No stress" wasn't explicitly mentioned, default to "No stress"
        if not treatments:
            treatments.add("No stress")

    result["treatment"] = sorted(list(treatments)) # Sort for consistent output

    # --- Extract Medium ---
    extracted_medium = "unspecified"
    medium_search_text = " ".join([
        _get_text(sample_metadata, 'growth_protocol_ch1'),
        _get_text(sample_metadata, 'description'),
        _get_text(sample_metadata, 'characteristics_ch1')
    ])

    if any(k in medium_search_text for k in ["soil", "greenhouse", "field"]):
        extracted_medium = "soil"
    elif any(k in medium_search_text for k in ["agar", "ms medium", "murashige skoog", "plate", "in vitro", "petri dish"]):
        extracted_medium = "agar"
    elif any(k in medium_search_text for k in ["liquid medium", "hydroponic", "solution culture", "liquid culture"]):
        extracted_medium = "liquid medium"
    # Infer "soil" if the tissue is a whole plant and medium is still unspecified
    elif extracted_tissue == "whole_plant" and extracted_medium == "unspecified":
        extracted_medium = "soil"
        
    result["medium"] = extracted_medium

    return result

def GSE31639_extractor(sample_metadata: dict) -> dict:
    """
    Extracts tissue, treatment, and medium information from a sample_metadata dictionary
    according to a predefined JSON schema.

    Args:
        sample_metadata (dict): A dictionary containing sample metadata, typically
                                from a GEO Series (GSE) entry.

    Returns:
        dict: A dictionary formatted as a JSON instance conforming to the schema:
              {"properties": {
                  "tissue": {"enum": ["root", "leaf", "flower", "shoot", "rosette", "bud", "whole_plant", "silique", "callus", "seed", "seedling", "unknown"]},
                  "treatment": {"items": {"enum": ["Drought Stress", "Dehidration Stress", "Salinity Stress", "Heat Stress", "Cold Stress", "Chemical Stress", "Nutrient Deficiency", "Biotic Stress", "Low Light Stress", "High Light Stress", "Red Light Stress", "Other Light Stress", "Other stress", "No stress"]}},
                  "medium": {"type": "string"}
              }, "required": ["tissue", "treatment", "medium"]}
    """

    extracted_data = {}

    # Helper function to safely get text from metadata fields
    def _get_text(metadata_dict, key):
        value = metadata_dict.get(key)
        if isinstance(value, list):
            return " ".join(value).lower()
        if isinstance(value, str):
            return value.lower()
        return ""

    # --- Extract Tissue ---
    tissue_enum_values = ["root", "leaf", "flower", "shoot", "rosette", "bud", "whole_plant", "silique", "callus", "seed", "seedling", "unknown"]
    tissue_keywords_map = {
        "whole seedling": "seedling",
        "seedlings": "seedling",
        "seedling": "seedling",
        "whole plant": "whole_plant",
        "root": "root",
        "leaf": "leaf",
        "flower": "flower",
        "shoot": "shoot",
        "rosette": "rosette",
        "bud": "bud",
        "silique": "silique",
        "callus": "callus",
        "seed": "seed",
    }
    extracted_tissue = "unknown"

    char_ch1 = _get_text(sample_metadata, 'characteristics_ch1')
    if "tissue:" in char_ch1:
        parts = char_ch1.split("tissue:")
        if len(parts) > 1:
            tissue_str = parts[1].split(',')[0].strip()
            # Try direct match first
            if tissue_str in tissue_enum_values:
                extracted_tissue = tissue_str
            else: # Then try keyword mapping
                for k, v in tissue_keywords_map.items():
                    if k in tissue_str:
                        extracted_tissue = v
                        break

    if extracted_tissue == "unknown":
        source_name = _get_text(sample_metadata, 'source_name_ch1')
        title = _get_text(sample_metadata, 'title')
        combined_text = source_name + " " + title
        for k, v in tissue_keywords_map.items():
            if k in combined_text:
                extracted_tissue = v
                break
    
    extracted_data["tissue"] = extracted_tissue

    # --- Extract Treatment ---
    treatment_enum_values = ["Drought Stress", "Dehidration Stress", "Salinity Stress", "Heat Stress", "Cold Stress", "Chemical Stress", "Nutrient Deficiency", "Biotic Stress", "Low Light Stress", "High Light Stress", "Red Light Stress", "Other Light Stress", "Other stress", "No stress"]
    extracted_treatments = set()

    treatment_text = _get_text(sample_metadata, 'treatment_protocol_ch1') + " " + \
                     _get_text(sample_metadata, 'characteristics_ch1') + " " + \
                     _get_text(sample_metadata, 'title') + " " + \
                     _get_text(sample_metadata, 'description')
    
    # Check for explicit "No stress" or control conditions
    if "none other than" in treatment_text or "no stress" in treatment_text or "control" in treatment_text or "untreated" in treatment_text:
        extracted_treatments.add("No stress")
    else:
        # Look for specific stresses
        if "drought" in treatment_text or "water deprivation" in treatment_text:
            extracted_treatments.add("Drought Stress")
        if "dehydration" in treatment_text:
            extracted_treatments.add("Dehidration Stress")
        if "salinity" in treatment_text or "salt stress" in treatment_text:
            extracted_treatments.add("Salinity Stress")
        if "heat stress" in treatment_text or "high temperature" in treatment_text:
            extracted_treatments.add("Heat Stress")
        if "cold stress" in treatment_text or "low temperature" in treatment_text:
            extracted_treatments.add("Cold Stress")
        if "chemical" in treatment_text or "herbicide" in treatment_text or "pesticide" in treatment_text or "drug" in treatment_text:
            extracted_treatments.add("Chemical Stress")
        if "nutrient deficiency" in treatment_text or "nitrogen starvation" in treatment_text or "phosphate starvation" in treatment_text:
            extracted_treatments.add("Nutrient Deficiency")
        if "biotic stress" in treatment_text or "pathogen" in treatment_text or "insect" in treatment_text or "fungus" in treatment_text:
            extracted_treatments.add("Biotic Stress")
        if "low light" in treatment_text:
            extracted_treatments.add("Low Light Stress")
        if "high light" in treatment_text:
            extracted_treatments.add("High Light Stress")
        if "red light" in treatment_text:
            extracted_treatments.add("Red Light Stress")
        # General "stress" if no specific one was found
        if not extracted_treatments and "stress" in treatment_text:
            extracted_treatments.add("Other stress")

    # If no treatments were found, default to "No stress"
    if not extracted_treatments:
        extracted_treatments.add("No stress")

    # If "No stress" is present with other stresses, remove it (unless it's the only one)
    if len(extracted_treatments) > 1 and "No stress" in extracted_treatments:
        extracted_treatments.remove("No stress")
    
    # Ensure all extracted treatments are valid enum values
    final_treatments = []
    for t in extracted_treatments:
        if t in treatment_enum_values:
            final_treatments.append(t)
        else:
            # Fallback for unrecognized but implied stress, or if a typo occurred
            final_treatments.append("Other stress") 

    extracted_data["treatment"] = sorted(list(set(final_treatments))) # Ensure unique and sorted

    # --- Extract Medium ---
    medium_keywords_map = {
        "ms plates": "MS medium",
        "agar": "MS medium",
        "1/2 ms": "MS medium",
        "murashige skoog": "MS medium",
        "soil": "soil",
        "potting mix": "soil",
        "vermiculite": "soil",
        "peat": "soil",
        "sand": "soil",
        "hydroponic": "hydroponic solution",
        "liquid medium": "hydroponic solution",
        "filter paper": "filter paper",
        "rockwool": "rockwool",
    }
    extracted_medium = "unspecified"

    growth_text = _get_text(sample_metadata, 'growth_protocol_ch1') + " " + \
                  _get_text(sample_metadata, 'characteristics_ch1') + " " + \
                  _get_text(sample_metadata, 'description')

    for k, v in medium_keywords_map.items():
        if k in growth_text:
            extracted_medium = v
            break

    # Inference based on tissue if still unspecified
    if extracted_medium == "unspecified":
        if extracted_tissue == "whole_plant":
            extracted_medium = "soil"
        elif extracted_tissue == "seedling" and ("plates" in growth_text or "agar" in growth_text or "ms" in growth_text):
            extracted_medium = "MS medium"
        elif extracted_tissue == "seed":
            # Seeds can be stored dry, or germinated on various media, so "unspecified" is appropriate
            extracted_medium = "unspecified" 

    extracted_data["medium"] = extracted_medium

    return extracted_data

def GSE35057_extractor(sample_metadata: dict) -> dict:
    """
    Extracts tissue, treatment, and medium information from a sample_metadata dictionary
    following a predefined schema.

    Args:
        sample_metadata (dict): A dictionary containing metadata for a sample.

    Returns:
        dict: A dictionary with 'tissue', 'treatment', and 'medium' keys,
              conforming to the specified JSON schema.
    """
    result = {
        "tissue": "unknown",
        "treatment": [],
        "medium": "unspecified"
    }

    # Helper to safely get and join list values into a single lowercase string
    def get_text(key):
        value = sample_metadata.get(key, [])
        if isinstance(value, list):
            return " ".join(value).lower()
        return str(value).lower()

    # --- Extract Tissue ---
    tissue_found = False
    characteristics_ch1_text = get_text('characteristics_ch1')
    source_name_ch1_text = get_text('source_name_ch1')
    title_text = get_text('title')
    description_text = get_text('description')

    # Priority 1: characteristics_ch1
    match = re.search(r'tissue:\s*([\w\s-]+)', characteristics_ch1_text)
    if match:
        extracted_tissue = match.group(1).strip()
        if "seedling" in extracted_tissue or "seedlings" in extracted_tissue:
            result["tissue"] = "seedling"
            tissue_found = True
        elif "root" in extracted_tissue:
            result["tissue"] = "root"
            tissue_found = True
        elif "leaf" in extracted_tissue:
            result["tissue"] = "leaf"
            tissue_found = True
        elif "flower" in extracted_tissue:
            result["tissue"] = "flower"
            tissue_found = True
        elif "shoot" in extracted_tissue:
            result["tissue"] = "shoot"
            tissue_found = True
        elif "rosette" in extracted_tissue:
            result["tissue"] = "rosette"
            tissue_found = True
        elif "bud" in extracted_tissue:
            result["tissue"] = "bud"
            tissue_found = True
        elif "whole plant" in extracted_tissue or "whole_plant" in extracted_tissue:
            result["tissue"] = "whole_plant"
            tissue_found = True
        elif "silique" in extracted_tissue:
            result["tissue"] = "silique"
            tissue_found = True
        elif "callus" in extracted_tissue:
            result["tissue"] = "callus"
            tissue_found = True
        elif "seed" in extracted_tissue:
            result["tissue"] = "seed"
            tissue_found = True

    # Priority 2: source_name_ch1 if not found yet
    if not tissue_found:
        if "seedling" in source_name_ch1_text or "seedlings" in source_name_ch1_text:
            result["tissue"] = "seedling"
            tissue_found = True
        elif "plant" in source_name_ch1_text and "seedling" not in source_name_ch1_text:
            result["tissue"] = "whole_plant"
            tissue_found = True

    # Priority 3: title/description if not found yet
    if not tissue_found:
        if "seedling" in title_text or "seedlings" in title_text or "seedling" in description_text or "seedlings" in description_text:
            result["tissue"] = "seedling"
            tissue_found = True
        elif "plant" in title_text or "plant" in description_text:
            result["tissue"] = "whole_plant"
            tissue_found = True

    # --- Extract Treatment ---
    extracted_treatments = set()
    treatment_search_text = get_text('characteristics_ch1') + " " + get_text('treatment_protocol_ch1') + " " + get_text('title') + " " + get_text('description')

    # Specific stress keywords
    if re.search(r'dehydration', treatment_search_text):
        extracted_treatments.add("Dehidration Stress")
    elif re.search(r'drought|water deficit', treatment_search_text):
        extracted_treatments.add("Drought Stress")
    
    if re.search(r'salinity|nacl', treatment_search_text):
        extracted_treatments.add("Salinity Stress")

    if re.search(r'heat|high temperature', treatment_search_text):
        extracted_treatments.add("Heat Stress")

    if re.search(r'cold|low temperature|chilling', treatment_search_text):
        extracted_treatments.add("Cold Stress")

    if re.search(r'chemical stress|herbicide|pesticide|heavy metal|toxic|pollutant|cadmium|aluminum|lead|mercury', treatment_search_text):
        extracted_treatments.add("Chemical Stress")

    if re.search(r'nutrient deficiency|low nitrogen|phosphate starvation|sulfur starvation|iron deficiency|boron deficiency|potassium deficiency', treatment_search_text):
        extracted_treatments.add("Nutrient Deficiency")

    if re.search(r'biotic stress|pathogen|insect|fungus|bacteria|virus|herbivory', treatment_search_text):
        extracted_treatments.add("Biotic Stress")

    if re.search(r'low r/fr|supplemental fr|far-red light', treatment_search_text):
        extracted_treatments.add("Red Light Stress")
    elif re.search(r'high light|excess light', treatment_search_text):
        extracted_treatments.add("High Light Stress")
    elif re.search(r'low light|shade', treatment_search_text):
        extracted_treatments.add("Low Light Stress")
    elif re.search(r'uv light|blue light|light stress', treatment_search_text):
        # Only add "Other Light Stress" if no more specific light stress was found
        if not any(s in extracted_treatments for s in ["Red Light Stress", "High Light Stress", "Low Light Stress"]):
            extracted_treatments.add("Other Light Stress")

    # Fallback logic for treatments if no specific stress was identified
    if not extracted_treatments:
        # Check for explicit "no stress" indicators
        if re.search(r'no stress|control|untreated|normal conditions|standard conditions|ambient conditions', treatment_search_text):
            extracted_treatments.add("No stress")
        # Check if any treatment was mentioned at all (even if not mapped to a specific enum)
        elif re.search(r'treatment:', characteristics_ch1_text) or get_text('treatment_protocol_ch1'):
            # If a treatment was mentioned but not mapped and not "no stress" related, it's "Other stress"
            extracted_treatments.add("Other stress")
        else:
            # No treatment mentioned anywhere, assume "No stress"
            extracted_treatments.add("No stress")

    result["treatment"] = sorted(list(extracted_treatments))

    # --- Extract Medium ---
    growth_protocol_text = get_text('growth_protocol_ch1')
    
    # Combine all relevant text for medium inference if not found in growth protocol
    all_relevant_text_for_medium = (
        growth_protocol_text + " " +
        characteristics_ch1_text + " " +
        source_name_ch1_text + " " +
        title_text + " " +
        description_text
    )

    if "half ms" in growth_protocol_text or "ms medium" in growth_protocol_text or "murashige and skoog" in growth_protocol_text:
        result["medium"] = "MS medium"
    elif "agar" in growth_protocol_text or "gelrite" in growth_protocol_text:
        result["medium"] = "agar medium"
    elif "soil" in growth_protocol_text or "potting mix" in growth_protocol_text:
        result["medium"] = "soil"
    elif "hydroponic" in growth_protocol_text or "liquid culture" in growth_protocol_text:
        result["medium"] = "hydroponic solution"
    elif "vermiculite" in growth_protocol_text or "perlite" in growth_protocol_text:
        result["medium"] = "vermiculite/perlite"
    elif "water" in growth_protocol_text:
        result["medium"] = "water"
    else:
        # Infer medium if tissue is whole_plant or seedling and no medium specified
        if result["tissue"] in ["whole_plant", "seedling"]:
            # Check if "soil" is mentioned anywhere in relevant text
            if "soil" in all_relevant_text_for_medium:
                result["medium"] = "soil"
            else:
                result["medium"] = "unspecified" # Default if no strong inference
        else:
            result["medium"] = "unspecified"

    return result

def GSE10247_extractor(sample_metadata: dict) -> dict:
    result = {}

    # Helper to safely get text from metadata fields
    def _get_text(key: str, default: str = "") -> str:
        value = sample_metadata.get(key)
        if isinstance(value, list) and value:
            return value[0].lower()
        elif isinstance(value, str):
            return value.lower()
        return default.lower()

    # Combine relevant text fields for easier searching for tissue and treatment
    # growth_protocol_ch1 is included here as it can provide context for tissue/treatment
    combined_text = (
        _get_text('title') + " " +
        _get_text('source_name_ch1') + " " +
        _get_text('characteristics_ch1') + " " +
        _get_text('treatment_protocol_ch1') + " " +
        _get_text('growth_protocol_ch1') + " " +
        _get_text('description')
    ).strip()

    # --- Extract Tissue ---
    extracted_tissue = "unknown"
    # TISSUE_ENUM is defined in the schema, but not explicitly needed for validation here
    # as we default to "unknown" if no match.

    # Prioritize specific terms. Order matters for overlapping keywords (e.g., "leaves" vs "phloem exudate").
    if "phloem exudate" in combined_text or "adult leaves" in combined_text or "leaves" in combined_text:
        extracted_tissue = "leaf"
    elif "root" in combined_text:
        extracted_tissue = "root"
    elif "seedling" in combined_text:
        extracted_tissue = "seedling"
    elif "seed" in combined_text:
        extracted_tissue = "seed"
    elif "flower" in combined_text:
        extracted_tissue = "flower"
    elif "silique" in combined_text:
        extracted_tissue = "silique"
    elif "callus" in combined_text:
        extracted_tissue = "callus"
    elif "bud" in combined_text:
        extracted_tissue = "bud"
    elif "shoot" in combined_text:
        extracted_tissue = "shoot"
    elif "rosette" in combined_text:
        extracted_tissue = "rosette"
    elif "whole plant" in combined_text:
        extracted_tissue = "whole_plant"
    
    result["tissue"] = extracted_tissue

    # --- Extract Treatment ---
    extracted_treatments = set() # Use a set to avoid duplicates
    # TREATMENT_ENUM is defined in the schema.

    # Check for specific stresses
    if "drought" in combined_text:
        extracted_treatments.add("Drought Stress")
    if "dehydration" in combined_text:
        extracted_treatments.add("Dehidration Stress")
    if "salinity" in combined_text or "salt stress" in combined_text:
        extracted_treatments.add("Salinity Stress")
    if "heat stress" in combined_text or "high temperature" in combined_text:
        extracted_treatments.add("Heat Stress")
    if "cold stress" in combined_text or "low temperature" in combined_text:
        extracted_treatments.add("Cold Stress")
    if "chemical" in combined_text:
        extracted_treatments.add("Chemical Stress")
    if "nutrient deficiency" in combined_text:
        extracted_treatments.add("Nutrient Deficiency")
    if "biotic stress" in combined_text or "pathogen" in combined_text or "insect" in combined_text:
        extracted_treatments.add("Biotic Stress")
    if "low light" in combined_text:
        extracted_treatments.add("Low Light Stress")
    if "high light" in combined_text:
        extracted_treatments.add("High Light Stress")
    if "red light" in combined_text:
        extracted_treatments.add("Red Light Stress")
    
    # General light stress, if not already covered by specific light stresses
    specific_light_stresses_added = any(t in extracted_treatments for t in ["Low Light Stress", "High Light Stress", "Red Light Stress"])
    if "light stress" in combined_text and not specific_light_stresses_added:
        extracted_treatments.add("Other Light Stress")
    
    # General stress, if not already covered by specific stresses
    any_stress_added = bool(extracted_treatments) # True if set is not empty
    if "stress" in combined_text and not any_stress_added:
        extracted_treatments.add("Other stress")

    # If no stress was found at all, assume "No stress"
    if not extracted_treatments:
        extracted_treatments.add("No stress")
    
    result["treatment"] = sorted(list(extracted_treatments)) # Sort for consistent output

    # --- Extract Medium ---
    extracted_medium = "unspecified"
    
    # Use specific fields for medium as they are most likely to contain this info
    growth_protocol_text = _get_text('growth_protocol_ch1')
    characteristics_text = _get_text('characteristics_ch1')
    source_name_text = _get_text('source_name_ch1')

    if "soil" in growth_protocol_text or "soil" in characteristics_text or "soil" in source_name_text:
        extracted_medium = "soil"
    elif "agar" in growth_protocol_text or "agar" in characteristics_text:
        extracted_medium = "agar"
    elif "ms medium" in growth_protocol_text or "murashige skoog" in growth_protocol_text:
        extracted_medium = "MS medium"
    elif "liquid medium" in growth_protocol_text or "hydroponic" in growth_protocol_text:
        extracted_medium = "liquid medium"
    
    result["medium"] = extracted_medium

    return result

import re

def GSE81477_extractor(sample_metadata: dict) -> dict:
    """
    Extracts tissue, treatment, and medium information from sample_metadata
    following a predefined schema for GSE81477 dataset.

    Args:
        sample_metadata (dict): A dictionary containing metadata for a sample.

    Returns:
        dict: A dictionary containing the extracted 'tissue', 'treatment', and 'medium'.
    """
    extracted_data = {}

    # Helper to safely get text from metadata fields, join lists, and convert to lowercase
    def get_text(key):
        value = sample_metadata.get(key)
        if isinstance(value, list):
            return " ".join(value).lower()
        return str(value).lower() if value else ""

    # Combine relevant text fields for easier searching across protocols and characteristics
    all_relevant_text = (
        get_text('title') + " " +
        get_text('source_name_ch1') + " " +
        get_text('characteristics_ch1') + " " +
        get_text('treatment_protocol_ch1') + " " +
        get_text('growth_protocol_ch1')
    )

    # 1. Extract Tissue
    tissue = "unknown"
    source_name_ch1_text = get_text('source_name_ch1')
    characteristics_ch1_text = get_text('characteristics_ch1')

    if "t87 cultured cell" in source_name_ch1_text or \
       "cell line: t87" in characteristics_ch1_text or \
       "suspension cells" in characteristics_ch1_text:
        tissue = "callus"
    # Add more specific tissue extraction logic here if other tissue types are expected
    # in this dataset and can be identified from keywords.

    extracted_data['tissue'] = tissue

    # 2. Extract Treatment
    treatments = set()
    treatment_protocol_text = get_text('treatment_protocol_ch1')
    growth_protocol_text = get_text('growth_protocol_ch1')

    # Check for specific stress keywords in relevant protocol fields
    # Cold Stress
    if re.search(r'cold treated|4°c|4 c', treatment_protocol_text + growth_protocol_text):
        treatments.add("Cold Stress")
    
    # Chemical Stress (cordycepin is a chemical treatment)
    if "cordycepin treatment" in treatment_protocol_text:
        treatments.add("Chemical Stress")

    # General stress keywords (using all_relevant_text for broader search)
    # Drought Stress / Dehydration Stress
    if re.search(r'drought|water deficit|dehydration', all_relevant_text):
        treatments.add("Drought Stress")
    # Salinity Stress
    if re.search(r'salinity|nacl', all_relevant_text):
        treatments.add("Salinity Stress")
    # Heat Stress
    if re.search(r'heat stress|high temperature|37°c|37 c', all_relevant_text):
        treatments.add("Heat Stress")
    # Nutrient Deficiency
    if re.search(r'nutrient deficiency|low nitrogen|low phosphate', all_relevant_text):
        treatments.add("Nutrient Deficiency")
    # Biotic Stress
    if re.search(r'biotic stress|pathogen|fungus|bacteria|virus|insect', all_relevant_text):
        treatments.add("Biotic Stress")
    # Light Stress
    if re.search(r'low light|dim light', all_relevant_text):
        treatments.add("Low Light Stress")
    if re.search(r'high light', all_relevant_text):
        treatments.add("High Light Stress")
    if re.search(r'red light', all_relevant_text):
        treatments.add("Red Light Stress")
    
    # If no specific treatments are identified, default to "No stress"
    if not treatments:
        treatments.add("No stress")

    extracted_data['treatment'] = sorted(list(treatments)) # Sort for consistent output

    # 3. Extract Medium
    medium = "unspecified"
    growth_protocol_text = get_text('growth_protocol_ch1')

    # Specific check for "JPL media" as it's a known medium in the example
    if "jpl media" in growth_protocol_text:
        medium = "JPL media"
    else:
        # Look for patterns like "grown in [NAME] media/medium" or "cultured in [NAME]"
        # This regex attempts to capture the full name of the medium including "media" or "medium" if present.
        match = re.search(r'grown in\s+((?:[\w\s-]+?)\s+(?:media|medium))|cultured in\s+([\w\s-]+)', growth_protocol_text)
        if match:
            # Find the first non-None capturing group
            for i in range(1, match.lastindex + 1):
                if match.group(i):
                    medium = match.group(i).strip()
                    break
    
    # Apply title casing for general medium names, but keep specific names like "JPL media" as is.
    if medium != "unspecified" and medium != "JPL media":
        extracted_data['medium'] = medium.title()
    else:
        extracted_data['medium'] = medium

    return extracted_data

def GSE19258_extractor(sample_metadata: dict) -> dict:
    extracted_data = {}

    # Helper to get text from metadata fields, handling lists and converting to lowercase
    def _get_text_from_metadata(data: dict, keys: list) -> str:
        text_parts = []
        for key in keys:
            if key in data:
                value = data[key]
                if isinstance(value, list):
                    text_parts.extend(value)
                elif isinstance(value, str):
                    text_parts.append(value)
        return " ".join(text_parts).lower()

    # Combine relevant text fields for comprehensive searching
    full_text = _get_text_from_metadata(sample_metadata, [
        'characteristics_ch1', 'source_name_ch1', 'title', 'description', 'growth_protocol_ch1'
    ])

    # 1. Extract Tissue
    tissue = "unknown"
    if "tissue: whole plant excluding roots" in full_text or "whole plant excluding roots" in full_text:
        tissue = "whole_plant"
    elif "whole plant" in full_text:
        tissue = "whole_plant"
    elif "seedling" in full_text:
        tissue = "seedling"
    elif "seed" in full_text:
        tissue = "seed"
    elif "root" in full_text:
        tissue = "root"
    elif "leaf" in full_text:
        tissue = "leaf"
    elif "flower" in full_text:
        tissue = "flower"
    elif "shoot" in full_text:
        tissue = "shoot"
    elif "rosette" in full_text:
        tissue = "rosette"
    elif "bud" in full_text:
        tissue = "bud"
    elif "silique" in full_text:
        tissue = "silique"
    elif "callus" in full_text:
        tissue = "callus"
    extracted_data["tissue"] = tissue

    # 2. Extract Treatment
    treatments = set()
    found_specific_stress = False
    found_light_stress = False

    if "drought" in full_text:
        treatments.add("Drought Stress")
        found_specific_stress = True
    if "dehydration" in full_text:
        treatments.add("Dehidration Stress")
        found_specific_stress = True
    if "salinity" in full_text or "salt stress" in full_text:
        treatments.add("Salinity Stress")
        found_specific_stress = True
    if "heat stress" in full_text or "high temperature" in full_text:
        treatments.add("Heat Stress")
        found_specific_stress = True
    if "cold stress" in full_text or "low temperature" in full_text:
        treatments.add("Cold Stress")
        found_specific_stress = True
    if "chemical stress" in full_text or "herbicide" in full_text or "pesticide" in full_text or "heavy metal" in full_text:
        treatments.add("Chemical Stress")
        found_specific_stress = True
    if "nutrient deficiency" in full_text or "nitrogen starvation" in full_text or "phosphate starvation" in full_text:
        treatments.add("Nutrient Deficiency")
        found_specific_stress = True
    if "biotic stress" in full_text or "pathogen" in full_text or "insect" in full_text:
        treatments.add("Biotic Stress")
        found_specific_stress = True

    if "low light" in full_text:
        treatments.add("Low Light Stress")
        found_light_stress = True
    if "high light" in full_text:
        treatments.add("High Light Stress")
        found_light_stress = True
    if "red light" in full_text:
        treatments.add("Red Light Stress")
        found_light_stress = True
    if "light stress" in full_text and not found_light_stress:
        treatments.add("Other Light Stress")
        found_light_stress = True

    # Add "Other stress" only if no more specific stress (including light stress) was found
    if "stress" in full_text and not found_specific_stress and not found_light_stress:
        treatments.add("Other stress")

    if not treatments:
        treatments.add("No stress")

    extracted_data["treatment"] = sorted(list(treatments))

    # 3. Extract Medium
    medium = "unspecified"
    if "soil" in full_text:
        medium = "soil"
    elif "agar" in full_text:
        medium = "agar"
    elif "ms medium" in full_text or "murashige skoog" in full_text:
        medium = "MS medium"
    elif "liquid medium" in full_text or "liquid culture" in full_text:
        medium = "liquid medium"
    elif "hydroponic" in full_text:
        medium = "hydroponic"
    elif "vermiculite" in full_text: # Example from guidance, not in schema enum but allowed as string
        medium = "vermiculite"
    # Infer soil if tissue is whole_plant and no other medium was explicitly found
    elif extracted_data["tissue"] == "whole_plant":
        medium = "soil"
    extracted_data["medium"] = medium

    return extracted_data

def GSE8248_extractor(sample_metadata: dict) -> dict:
    extracted_data = {}

    # Helper to safely get and join text from potential keys
    def _get_text(data: dict, keys: list) -> str:
        text_parts = []
        for key in keys:
            if key in data:
                value = data[key]
                if isinstance(value, list):
                    text_parts.extend(value)
                elif isinstance(value, str):
                    text_parts.append(value)
        return " ".join(text_parts).lower()

    # --- Extract Tissue ---
    tissue_text = _get_text(sample_metadata, ['source_name_ch1', 'characteristics_ch1', 'title'])
    
    # Define tissue mapping based on schema enum and common synonyms
    tissue_mapping = {
        "root": "root",
        "leaf": "leaf",
        "flower": "flower",
        "shoot": "shoot",
        "rosette": "rosette",
        "bud": "bud",
        "whole plant": "whole_plant",
        "silique": "silique",
        "callus": "callus",
        "seed": "seed",
        "seedling": "seedling",
        "protoplast": "leaf", # Protoplasts are isolated from leaves
        "mesophyll cells": "leaf" # Mesophyll cells are from leaves
    }
    
    found_tissue = "unknown"
    for keyword, schema_value in tissue_mapping.items():
        if keyword in tissue_text:
            found_tissue = schema_value
            break
    
    extracted_data["tissue"] = found_tissue

    # --- Extract Treatment ---
    treatment_text = _get_text(sample_metadata, ['treatment_protocol_ch1', 'title', 'description', 'characteristics_ch1'])
    
    # Define treatment mapping based on schema enum and common synonyms
    treatment_keywords = {
        "drought": "Drought Stress",
        "dehydration": "Dehidration Stress", # Schema uses "Dehidration"
        "salinity": "Salinity Stress",
        "salt": "Salinity Stress",
        "heat": "Heat Stress",
        "cold": "Cold Stress",
        "chemical": "Chemical Stress",
        "nutrient deficiency": "Nutrient Deficiency",
        "biotic": "Biotic Stress",
        "low light": "Low Light Stress",
        "high light": "High Light Stress",
        "red light": "Red Light Stress",
        # "other light stress" and "other stress" are harder to infer generically
    }
    
    found_treatments = set()
    for keyword, schema_value in treatment_keywords.items():
        if keyword in treatment_text:
            found_treatments.add(schema_value)
            
    # Check for "no stress" or control conditions
    if not found_treatments:
        # Look for explicit "control" or "no stress" indicators
        if "control" in treatment_text or "no stress" in treatment_text or "untreated" in treatment_text:
            found_treatments.add("No stress")
        else:
            # If no specific stress or control is mentioned, default to "No stress"
            found_treatments.add("No stress")

    extracted_data["treatment"] = sorted(list(found_treatments)) # Sort for consistent output

    # --- Extract Medium ---
    medium_text = _get_text(sample_metadata, ['growth_protocol_ch1', 'treatment_protocol_ch1', 'characteristics_ch1', 'source_name_ch1'])
    
    # Define medium mapping based on common terms, prioritizing more specific media
    medium_mapping = {
        "mannitol buffer": "mannitol buffer",
        "ms medium": "MS medium",
        "murashige and skoog": "MS medium",
        "agar": "agar",
        "hydroponic": "hydroponic solution",
        "liquid medium": "liquid medium",
        "water": "water",
        "soil": "soil",
        "vermiculite": "vermiculite",
        "perlite": "perlite",
        "rockwool": "rockwool",
        "sand": "sand",
        "peat": "peat",
        "gelrite": "gelrite"
    }
    
    found_medium = "unspecified"
    for keyword, schema_value in medium_mapping.items():
        if keyword in medium_text:
            found_medium = schema_value
            break
            
    # If tissue is whole_plant and no specific medium was found, default to soil
    if found_tissue == "whole_plant" and found_medium == "unspecified":
        found_medium = "soil"
    
    extracted_data["medium"] = found_medium

    return extracted_data

def GSE39236_extractor(sample_metadata: dict) -> dict:
    extracted_data = {}

    # Helper to safely get a string from a list or return an empty string
    def get_string_value(key):
        value = sample_metadata.get(key)
        if isinstance(value, list) and value:
            return " ".join(value)
        elif isinstance(value, str):
            return value
        return ""

    # --- Extract Tissue ---
    tissue_raw = ""
    characteristics_ch1 = get_string_value('characteristics_ch1')
    if characteristics_ch1:
        # Look for "tissue: [value]"
        match = re.search(r'tissue:\s*([^;,\n]+)', characteristics_ch1, re.IGNORECASE)
        if match:
            tissue_raw = match.group(1).strip().lower()

    tissue_mapping = {
        "seedlings": "seedling",
        "seedling": "seedling",
        "root": "root",
        "leaf": "leaf",
        "flower": "flower",
        "shoot": "shoot",
        "rosette": "rosette",
        "bud": "bud",
        "whole plant": "whole_plant",
        "silique": "silique",
        "callus": "callus",
        "seed": "seed",
    }
    extracted_data['tissue'] = tissue_mapping.get(tissue_raw, "unknown")

    # --- Extract Treatment ---
    found_treatments = set()
    all_text_for_treatment = [
        get_string_value('characteristics_ch1'),
        get_string_value('treatment_protocol_ch1'),
        get_string_value('source_name_ch1'),
        get_string_value('title'),
        get_string_value('description')
    ]

    treatment_keywords = {
        "Drought Stress": ["drought", "water deficit"],
        "Dehidration Stress": ["dehydration"], # Note: schema uses "Dehidration"
        "Salinity Stress": ["salt", "salinity", "nacl"],
        "Heat Stress": ["heat", "high temperature", "thermal stress"],
        "Cold Stress": ["cold", "low temperature", "chilling", "freezing"],
        "Chemical Stress": ["chemical", "herbicide", "pesticide", "drug", "cadmium", "hormone", "metal stress", "oxidative stress", "paraquat", "methyl viologen"],
        "Nutrient Deficiency": ["nutrient deficiency", "low nitrogen", "phosphate starvation", "nitrogen starvation", "phosphorus starvation", "sulfur starvation", "iron deficiency"],
        "Biotic Stress": ["biotic", "pathogen", "insect", "fungus", "bacteria", "virus", "infection", "herbivory"],
        "Low Light Stress": ["low light", "shade", "darkness"],
        "High Light Stress": ["high light", "uv-b"],
        "Red Light Stress": ["red light"],
        "Other Light Stress": ["light stress"], # General light stress if not specific
    }

    for text in all_text_for_treatment:
        text_lower = text.lower()
        for schema_treatment, keywords in treatment_keywords.items():
            for keyword in keywords:
                if keyword in text_lower:
                    found_treatments.add(schema_treatment)
                    break # Move to next schema_treatment if one keyword matches

    # Special handling for "No stress" or default
    if not found_treatments:
        no_stress_explicit = False
        for text in all_text_for_treatment:
            text_lower = text.lower()
            if "no stress" in text_lower or "control" in text_lower or "untreated" in text_lower or "wild type" in text_lower and not any(k in text_lower for k_list in treatment_keywords.values() for k in k_list):
                no_stress_explicit = True
                break
        if no_stress_explicit:
            found_treatments.add("No stress")
        else:
            # If no specific stress or "no stress" is found, assume "No stress" as a default
            found_treatments.add("No stress")
            
    extracted_data['treatment'] = sorted(list(found_treatments))

    # --- Extract Medium ---
    medium_raw = "unspecified"
    growth_protocol_ch1 = get_string_value('growth_protocol_ch1')
    
    if growth_protocol_ch1:
        # Prioritize specific medium names like MS medium
        ms_match = re.search(r'(Murashige and Skoog \(MS\) medium(?: agar plates)?(?: \([^)]+\))?)', growth_protocol_ch1, re.IGNORECASE)
        if ms_match:
            medium_raw = ms_match.group(1).strip()
        else:
            # General search for "medium", "agar", "soil", "plates", "solution"
            # This regex tries to capture the phrase describing the medium
            medium_match = re.search(r'(?:grown on|in|cultured on|in)\s+([^,.;]+(?:medium|agar|soil|plates|solution|hydroponics|vermiculite))', growth_protocol_ch1, re.IGNORECASE)
            if medium_match:
                medium_raw = medium_match.group(1).strip()
            else:
                # Fallback for "soil" if whole plant and "soil" is mentioned
                if extracted_data['tissue'] == "whole_plant" and "soil" in growth_protocol_ch1.lower():
                    medium_raw = "soil"
    
    # If still unspecified and tissue is whole_plant, check other fields for "soil"
    if medium_raw == "unspecified" and extracted_data['tissue'] == "whole_plant":
        if "soil" in get_string_value('description').lower() or "soil" in get_string_value('source_name_ch1').lower():
            medium_raw = "soil"

    extracted_data['medium'] = medium_raw

    return extracted_data

import re

def GSE50679_extractor(sample_metadata: dict) -> dict:
    """
    Extracts tissue, treatment, and medium information from a sample_metadata dictionary
    following a specific JSON schema.

    Args:
        sample_metadata (dict): A dictionary containing metadata for a sample.

    Returns:
        dict: A dictionary with 'tissue', 'treatment', and 'medium' fields
              conforming to the specified schema.
    """
    extracted_data = {
        "tissue": "unknown",
        "treatment": [],
        "medium": "unspecified"
    }

    # Helper to get a single string from a list, or an empty string
    def get_first_or_empty(key):
        val = sample_metadata.get(key)
        if isinstance(val, list) and val:
            return val[0]
        return ""

    # Combine relevant text fields for easier searching
    characteristics_ch1_str = " ".join(sample_metadata.get('characteristics_ch1', []))
    source_name_ch1_str = get_first_or_empty('source_name_ch1')
    treatment_protocol_ch1_str = get_first_or_empty('treatment_protocol_ch1')
    growth_protocol_ch1_str = get_first_or_empty('growth_protocol_ch1')
    title_str = get_first_or_empty('title')

    # --- Extract Tissue ---
    tissue_mapping = {
        "seedling": "seedling", "seedlings": "seedling",
        "leaf": "leaf", "leaves": "leaf",
        "root": "root", "roots": "root",
        "flower": "flower", "flowers": "flower",
        "shoot": "shoot", "shoots": "shoot",
        "rosette": "rosette", "rosettes": "rosette",
        "bud": "bud", "buds": "bud",
        "whole plant": "whole_plant", "whole-plant": "whole_plant",
        "silique": "silique", "siliques": "silique",
        "callus": "callus",
        "seed": "seed", "seeds": "seed"
    }

    # 1. From characteristics_ch1 (most specific)
    match = re.search(r'tissue:\s*([^;]+)', characteristics_ch1_str, re.IGNORECASE)
    if match:
        raw_tissue = match.group(1).strip().lower()
        for key, mapped_val in tissue_mapping.items():
            if key in raw_tissue:
                extracted_data["tissue"] = mapped_val
                break
    
    # 2. From source_name_ch1 or title if not found yet
    if extracted_data["tissue"] == "unknown":
        search_text = (source_name_ch1_str + " " + title_str).lower()
        for key, mapped_val in tissue_mapping.items():
            if key in search_text:
                extracted_data["tissue"] = mapped_val
                break

    # --- Extract Treatment ---
    found_treatments = set()
    treatment_enum = [
        "Drought Stress", "Dehidration Stress", "Salinity Stress", "Heat Stress",
        "Cold Stress", "Chemical Stress", "Nutrient Deficiency", "Biotic Stress",
        "Low Light Stress", "High Light Stress", "Red Light Stress",
        "Other Light Stress", "Other stress", "No stress"
    ]

    # 1. From characteristics_ch1 (primary source for specific sample treatment)
    match = re.search(r'treatment:\s*([^;]+)', characteristics_ch1_str, re.IGNORECASE)
    if match:
        raw_treatment = match.group(1).strip().lower()
        if "non-stress" in raw_treatment or "no stress" in raw_treatment:
            found_treatments.add("No stress")
        elif "hypoxic" in raw_treatment or "deprived of air" in raw_treatment:
            found_treatments.add("Other stress")
        elif "drought" in raw_treatment:
            found_treatments.add("Drought Stress")
        elif "salinity" in raw_treatment:
            found_treatments.add("Salinity Stress")
        elif "heat" in raw_treatment:
            found_treatments.add("Heat Stress")
        elif "cold" in raw_treatment:
            found_treatments.add("Cold Stress")

    # 2. From treatment_protocol_ch1 (for more detailed stress types)
    protocol_text = treatment_protocol_ch1_str.lower()
    
    if "drought" in protocol_text or "water deficit" in protocol_text:
        found_treatments.add("Drought Stress")
    if "dehydration" in protocol_text:
        found_treatments.add("Dehidration Stress")
    if "salinity" in protocol_text or "nacl" in protocol_text:
        found_treatments.add("Salinity Stress")
    if "heat" in protocol_text or "high temperature" in protocol_text:
        found_treatments.add("Heat Stress")
    if "cold" in protocol_text or "low temperature" in protocol_text:
        found_treatments.add("Cold Stress")
    if "chemical" in protocol_text or "herbicide" in protocol_text or "pesticide" in protocol_text:
        found_treatments.add("Chemical Stress")
    if "nutrient deficiency" in protocol_text:
        found_treatments.add("Nutrient Deficiency")
    if "biotic stress" in protocol_text or "pathogen" in protocol_text or "insect" in protocol_text:
        found_treatments.add("Biotic Stress")
    if "low light" in protocol_text:
        found_treatments.add("Low Light Stress")
    if "high light" in protocol_text:
        found_treatments.add("High Light Stress")
    if "red light" in protocol_text:
        found_treatments.add("Red Light Stress")
    if "light stress" in protocol_text and not any(s in protocol_text for s in ["low light", "high light", "red light"]):
        found_treatments.add("Other Light Stress")
    if "hypoxic stress" in protocol_text or "deprived of air" in protocol_text:
        found_treatments.add("Other stress")
    
    # 3. Check title for quick hints if no treatment found yet or for confirmation
    if not found_treatments: # Only if no treatments were found from characteristics or protocol
        title_lower = title_str.lower()
        if "non-stress" in title_lower or "no stress" in title_lower or "_ns" in title_lower:
            found_treatments.add("No stress")
        elif "stress" in title_lower:
            found_treatments.add("Other stress")

    # Final treatment logic:
    # If "No stress" was explicitly identified for *this specific sample*
    # (e.g., in characteristics_ch1 or title), it should be the sole treatment.
    if "No stress" in found_treatments and (
        "non-stress" in characteristics_ch1_str.lower() or
        "no stress" in characteristics_ch1_str.lower() or
        "_ns" in title_str.lower() # Specific to this dataset's naming convention
    ):
        extracted_data["treatment"] = ["No stress"]
    elif found_treatments:
        # Filter to ensure only valid enum values are returned and sort for consistency
        extracted_data["treatment"] = sorted([t for t in list(found_treatments) if t in treatment_enum])
    else:
        # Default if no treatment is found at all
        extracted_data["treatment"] = ["No stress"]


    # --- Extract Medium ---
    medium_text = growth_protocol_ch1_str.lower()
    if "ms media" in medium_text or "murashige skoog" in medium_text or "ms salts" in medium_text:
        extracted_data["medium"] = "MS media"
    elif "soil" in medium_text:
        extracted_data["medium"] = "soil"
    elif "liquid medium" in medium_text or "hydroponic" in medium_text:
        extracted_data["medium"] = "liquid medium"
    # If only agar/phytagel is mentioned without MS, it's safer to leave as "unspecified"
    # as "agar medium" is not in the schema and "MS media" is more specific.

    return extracted_data

def GSE72954_extractor(sample_metadata: dict) -> dict:
    """
    Extracts biological experimental information (tissue, treatment, medium)
    from the input sample_metadata dictionary for GSE72954,
    conforming to a predefined JSON schema.

    Args:
        sample_metadata (dict): A dictionary containing metadata for a sample.

    Returns:
        dict: A dictionary containing the extracted 'tissue', 'treatment', and 'medium'
              information, formatted according to the specified schema.
    """

    def get_text(data, key):
        """Helper to safely get and join text from a list value, converting to lowercase."""
        value = data.get(key)
        if value and isinstance(value, list):
            return " ".join(value).lower()
        return ""

    # --- Extract Tissue ---
    # Schema enum for tissue: ["root", "leaf", "flower", "shoot", "rosette", "bud",
    #                          "whole_plant", "silique", "callus", "seed", "seedling", "unknown"]
    tissue_keywords = [
        "root", "leaf", "flower", "shoot", "rosette", "bud",
        "whole_plant", "silique", "callus", "seed", "seedling"
    ]
    extracted_tissue = "unknown"

    # Prioritize 'source_name_ch1' as it often directly states the tissue
    source_name_text = get_text(sample_metadata, 'source_name_ch1')
    if source_name_text:
        for keyword in tissue_keywords:
            if keyword in source_name_text:
                extracted_tissue = keyword
                break
    
    # If not found in source_name_ch1, check other fields (less common for direct tissue info)
    if extracted_tissue == "unknown":
        characteristics_text = get_text(sample_metadata, 'characteristics_ch1')
        title_text = get_text(sample_metadata, 'title')
        description_text = get_text(sample_metadata, 'description')
        combined_other_text = " ".join([characteristics_text, title_text, description_text])
        for keyword in tissue_keywords:
            if keyword in combined_other_text:
                extracted_tissue = keyword
                break

    # --- Extract Treatment ---
    # Schema enum for treatment items: ["Drought Stress", "Dehidration Stress", "Salinity Stress",
    #                                   "Heat Stress", "Cold Stress", "Chemical Stress",
    #                                   "Nutrient Deficiency", "Biotic Stress", "Low Light Stress",
    #                                   "High Light Stress", "Red Light Stress", "Other Light Stress",
    #                                   "Other stress", "No stress"]
    treatment_keywords = [
        "Drought Stress", "Dehidration Stress", "Salinity Stress", "Heat Stress", "Cold Stress",
        "Chemical Stress", "Nutrient Deficiency", "Biotic Stress", "Low Light Stress",
        "High Light Stress", "Red Light Stress", "Other Light Stress", "Other stress"
    ]
    found_treatments = set()

    # Combine text from relevant fields for treatment search
    treatment_text_sources = [
        get_text(sample_metadata, 'treatment_protocol_ch1'),
        get_text(sample_metadata, 'characteristics_ch1'),
        get_text(sample_metadata, 'title'),
        get_text(sample_metadata, 'description')
    ]
    combined_treatment_text = " ".join(treatment_text_sources)

    for treatment_type in treatment_keywords:
        # Perform case-insensitive search
        if treatment_type.lower() in combined_treatment_text:
            found_treatments.add(treatment_type)

    if not found_treatments:
        extracted_treatment = ["No stress"]
    else:
        # Sort for consistent output if multiple treatments are found
        extracted_treatment = sorted(list(found_treatments))

    # --- Extract Medium ---
    # Schema for medium: string
    extracted_medium = "unspecified"
    growth_text = get_text(sample_metadata, 'growth_protocol_ch1')

    if growth_text:
        # Prioritize more specific or common media for plant samples
        if "gm plates" in growth_text or "g.m. plates" in growth_text or "agar plates" in growth_text or "ms plates" in growth_text:
            extracted_medium = "GM plates"
        elif "soil" in growth_text:
            extracted_medium = "soil"
        elif "agar" in growth_text:
            extracted_medium = "agar"
        elif "hydroponic" in growth_text or "hoagland solution" in growth_text or "liquid medium" in growth_text or "ms medium liquid" in growth_text:
            extracted_medium = "liquid medium"
        # Add more specific checks here if other common media types are expected in this dataset.

    return {
        "tissue": extracted_tissue,
        "treatment": extracted_treatment,
        "medium": extracted_medium
    }

import re

def GSE31158_extractor(sample_metadata: dict) -> dict:
    # Constants for enum values
    TISSUE_ENUM = ["root", "leaf", "flower", "shoot", "rosette", "bud", "whole_plant", "silique", "callus", "seed", "seedling", "unknown"]
    TREATMENT_ENUM = ["Drought Stress", "Dehidration Stress", "Salinity Stress", "Heat Stress", "Cold Stress", "Chemical Stress", "Nutrient Deficiency", "Biotic Stress", "Low Light Stress", "High Light Stress", "Red Light Stress", "Other Light Stress", "Other stress", "No stress"]

    extracted_data = {
        "tissue": "unknown",
        "treatment": [],
        "medium": "unspecified"
    }

    # Helper to get the first element of a list or an empty string
    def _get_first_or_empty(data, key):
        val = data.get(key)
        if isinstance(val, list) and val:
            return val[0]
        return ""

    # Helper to get all elements of a list or an empty list
    def _get_list_or_empty(data, key):
        val = data.get(key)
        if isinstance(val, list):
            return val
        return []

    # Combine relevant fields for easier searching (all lowercased)
    all_characteristics_str = " ".join(_get_list_or_empty(sample_metadata, 'characteristics_ch1')).lower()
    source_name_str = _get_first_or_empty(sample_metadata, 'source_name_ch1').lower()
    title_str = _get_first_or_empty(sample_metadata, 'title').lower()
    treatment_protocol_str = _get_first_or_empty(sample_metadata, 'treatment_protocol_ch1').lower()
    growth_protocol_str = _get_first_or_empty(sample_metadata, 'growth_protocol_ch1').lower()

    # Comprehensive search space for tissue and treatment
    search_space = f"{all_characteristics_str} {source_name_str} {title_str} {treatment_protocol_str}"

    # 1. Extract Tissue
    # Prioritize 'characteristics_ch1' as it often explicitly states "tissue: X"
    for char_item in _get_list_or_empty(sample_metadata, 'characteristics_ch1'):
        char_item_lower = char_item.lower()
        if "tissue:" in char_item_lower:
            found_tissue = char_item_lower.split("tissue:", 1)[1].strip()
            # Map to enum values
            for enum_tissue in TISSUE_ENUM:
                if enum_tissue == "unknown": continue
                # Check for exact match or match with spaces instead of underscores
                if found_tissue == enum_tissue or found_tissue == enum_tissue.replace("_", " "):
                    extracted_data["tissue"] = enum_tissue
                    break
            if extracted_data["tissue"] != "unknown":
                break # Found tissue, stop searching
    
    # If not found in characteristics_ch1, search general search_space
    if extracted_data["tissue"] == "unknown":
        for tissue_type in TISSUE_ENUM:
            if tissue_type == "unknown": continue # 'unknown' is a default, not a search term
            # Handle "whole_plant" vs "plant" and "seedling"
            search_term = tissue_type.replace("_", " ")
            if search_term in search_space:
                extracted_data["tissue"] = tissue_type
                break
        # Specific check for "seedling" if not caught by general search
        if extracted_data["tissue"] == "unknown" and "seedling" in search_space:
            extracted_data["tissue"] = "seedling"


    # 2. Extract Treatment
    found_treatments = set()

    # Keywords for treatment mapping (ordered from most specific to least specific for overlapping terms)
    treatment_keywords_map = {
        "red light": "Red Light Stress",
        "low light": "Low Light Stress",
        "high light": "High Light Stress",
        "light stress": "Other Light Stress", # General light stress
        "drought": "Drought Stress",
        "dehydration": "Dehidration Stress",
        "salinity": "Salinity Stress",
        "salt": "Salinity Stress",
        "heat": "Heat Stress",
        "cold": "Cold Stress",
        "chemical": "Chemical Stress",
        "nutrient deficiency": "Nutrient Deficiency",
        "biotic": "Biotic Stress",
        "hypoxia": "Other stress", # Oxygen deprivation
        "normoxia (mock)": "No stress", # Explicit control condition
        "mock": "No stress", # General mock
        "control": "No stress", # General control
        "no stress": "No stress", # Explicit no stress
    }

    # Search for treatments in the combined search space
    for keyword, mapped_treatment in treatment_keywords_map.items():
        if keyword in search_space:
            found_treatments.add(mapped_treatment)
            
    # If no specific treatments (including explicit "No stress" conditions like "mock") were found,
    # then default to "No stress".
    if not found_treatments:
        found_treatments.add("No stress")
    
    # Ensure all found treatments are in the enum and convert to sorted list
    extracted_data["treatment"] = sorted([t for t in list(found_treatments) if t in TREATMENT_ENUM])


    # 3. Extract Medium
    if growth_protocol_str:
        # Look for common medium indicators
        # Example: "grown on 0.5x MS, 1% sucrose, 0.4% phytagel plates grown vertically for 10 days under long day conditions"
        
        # Pattern 1: "grown on X", "on X", "grown in X", "in X"
        # Capture content after the prefix until a common separator or end of relevant phrase
        medium_match = re.search(
            r'(?:grown on|on|grown in|in)\s+([^,;.]+?)(?: plates| grown vertically| under long day| for \d+ days| in \d+ days| under| at| \d+ hr| \d+ days|$)',
            growth_protocol_str
        )
        if medium_match:
            medium_str = medium_match.group(1).strip()
            # Further clean up if necessary, e.g., remove trailing "plates" if not caught by regex
            medium_str = re.sub(r'\s*plates$', '', medium_str).strip()
            if medium_str:
                extracted_data["medium"] = medium_str
        
        # If not found by pattern, try to infer from keywords
        if extracted_data["medium"] == "unspecified":
            if "ms" in growth_protocol_str:
                # Try to extract "X MS" or "MS, Y"
                ms_match = re.search(r'([\d.]+x\s*ms(?:,\s*[^,;.]+)*)', growth_protocol_str)
                if ms_match:
                    extracted_data["medium"] = ms_match.group(1).strip()
                else:
                    extracted_data["medium"] = "MS medium" # Generic if specific not found
            elif "agar" in growth_protocol_str:
                extracted_data["medium"] = "agar"
            elif "soil" in growth_protocol_str:
                extracted_data["medium"] = "soil"
            elif "hydroponic" in growth_protocol_str:
                extracted_data["medium"] = "hydroponic solution"
            elif "water" in growth_protocol_str:
                extracted_data["medium"] = "water"
            elif "vermiculite" in growth_protocol_str:
                extracted_data["medium"] = "vermiculite"
            elif "sand" in growth_protocol_str:
                extracted_data["medium"] = "sand"
            elif "peat" in growth_protocol_str:
                extracted_data["medium"] = "peat"
            elif "rockwool" in growth_protocol_str:
                extracted_data["medium"] = "rockwool"
            elif "gelrite" in growth_protocol_str or "phytagel" in growth_protocol_str:
                extracted_data["medium"] = "gel medium" # Generic for gelling agents

    # Final inference for medium if still unspecified and tissue is whole_plant
    if extracted_data["medium"] == "unspecified" and extracted_data["tissue"] == "whole_plant":
        extracted_data["medium"] = "soil"

    return extracted_data

def GSE34897_extractor(sample_metadata: dict) -> dict:
    extracted_data = {
        "tissue": "unknown",
        "treatment": [],
        "medium": "unspecified"
    }

    def get_text(field_name):
        return " ".join(sample_metadata.get(field_name, [])).lower()

    all_text = (
        get_text('title') + " " +
        get_text('source_name_ch1') + " " +
        get_text('characteristics_ch1') + " " +
        get_text('description')
    )

    # --- Extract Tissue ---
    tissue_keywords = {
        "seedlings": "seedling",
        "seed": "seed",
        "root": "root",
        "leaf": "leaf",
        "flower": "flower",
        "shoot": "shoot",
        "rosette": "rosette",
        "bud": "bud",
        "whole plant": "whole_plant",
        "silique": "silique",
        "callus": "callus"
    }
    found_tissue = None
    for keyword, schema_value in sorted(tissue_keywords.items(), key=lambda item: len(item[0]), reverse=True):
        if keyword in all_text:
            found_tissue = schema_value
            break
    if found_tissue:
        extracted_data["tissue"] = found_tissue

    # --- Extract Treatment ---
    found_specific_treatments = set()
    specific_treatment_keywords = {
        "drought": "Drought Stress",
        "dehydration": "Dehidration Stress",
        "salinity": "Salinity Stress",
        "heat": "Heat Stress",
        "cold": "Cold Stress",
        "chemical": "Chemical Stress",
        "ethanol": "Chemical Stress",
        "mg132": "Chemical Stress",
        "nutrient deficiency": "Nutrient Deficiency",
        "biotic": "Biotic Stress",
        "low light": "Low Light Stress",
        "high light": "High Light Stress",
        "constant light": "High Light Stress",
        "red light": "Red Light Stress",
        "other light": "Other Light Stress",
    }
    
    for keyword, schema_value in sorted(specific_treatment_keywords.items(), key=lambda item: len(item[0]), reverse=True):
        if keyword in all_text:
            found_specific_treatments.add(schema_value)

    if found_specific_treatments:
        extracted_data["treatment"] = sorted(list(found_specific_treatments))
    else:
        # If no specific treatments, check for general "stress"
        # Ensure "no stress" isn't misinterpreted as "Other stress"
        if "stress" in all_text and "no stress" not in all_text:
            extracted_data["treatment"] = ["Other stress"]
        else:
            extracted_data["treatment"] = ["No stress"]

    # --- Extract Medium ---
    medium_keywords = {
        "agar": "agar",
        "ms medium": "MS medium",
        "murashige skoog": "MS medium",
        "soil": "soil",
        "liquid medium": "liquid medium",
        "hydroponic": "hydroponic",
        "vermiculite": "vermiculite",
        "peat": "peat",
        "plate": "agar" # Often implies agar plate
    }
    found_medium = None
    for keyword, schema_value in sorted(medium_keywords.items(), key=lambda item: len(item[0]), reverse=True):
        if keyword in all_text:
            found_medium = schema_value
            break
    if found_medium:
        extracted_data["medium"] = found_medium

    return extracted_data

def GSE5619_extractor(sample_metadata: dict) -> dict:
    extracted_data = {
        "tissue": "unknown",
        "treatment": ["No stress"],
        "medium": "unspecified"
    }

    # Helper to safely get a string from a list or default
    def _get_string_value(key, default=""):
        val = sample_metadata.get(key)
        if isinstance(val, list) and val:
            return " ".join(val) # Join all elements for broader search
        return default

    # Helper to extract value after a keyword in a list of strings
    def _extract_info_from_list(key_list, keyword):
        if not isinstance(key_list, list):
            return None
        for item in key_list:
            if isinstance(item, str) and keyword.lower() in item.lower():
                # Find the exact keyword with correct casing to split
                idx = item.lower().find(keyword.lower())
                if idx != -1:
                    return item[idx + len(keyword):].strip()
        return None

    # --- Extract Tissue ---
    characteristics_ch1 = sample_metadata.get('characteristics_ch1', [])
    tissue_raw = _extract_info_from_list(characteristics_ch1, "Tissue:")

    if tissue_raw:
        tissue_raw_lower = tissue_raw.lower()
        if "mature leaves" in tissue_raw_lower or "leaf" in tissue_raw_lower:
            extracted_data["tissue"] = "leaf"
        elif "whole plant" in tissue_raw_lower:
            extracted_data["tissue"] = "whole_plant"
        elif "seedling" in tissue_raw_lower:
            extracted_data["tissue"] = "seedling"
        elif "seed" in tissue_raw_lower:
            extracted_data["tissue"] = "seed"
        elif "root" in tissue_raw_lower:
            extracted_data["tissue"] = "root"
        elif "flower" in tissue_raw_lower:
            extracted_data["tissue"] = "flower"
        elif "shoot" in tissue_raw_lower:
            extracted_data["tissue"] = "shoot"
        elif "rosette" in tissue_raw_lower:
            extracted_data["tissue"] = "rosette"
        elif "bud" in tissue_raw_lower:
            extracted_data["tissue"] = "bud"
        elif "silique" in tissue_raw_lower:
            extracted_data["tissue"] = "silique"
        elif "callus" in tissue_raw_lower:
            extracted_data["tissue"] = "callus"
        # Default is "unknown" if not matched

    # --- Extract Treatment ---
    treatments_found = set()
    treatment_protocol_ch1_str = _get_string_value('treatment_protocol_ch1')
    description_list = sample_metadata.get('description', [])
    treatment_desc_raw = _extract_info_from_list(description_list, "Treatment:")

    all_treatment_info = []
    if treatment_protocol_ch1_str:
        all_treatment_info.append(treatment_protocol_ch1_str)
    if treatment_desc_raw:
        all_treatment_info.append(treatment_desc_raw)

    for info_str in all_treatment_info:
        info_str_lower = info_str.lower()
        if "drought" in info_str_lower:
            treatments_found.add("Drought Stress")
        if "dehydration" in info_str_lower:
            treatments_found.add("Dehidration Stress")
        if "salinity" in info_str_lower or "salt" in info_str_lower:
            treatments_found.add("Salinity Stress")
        if "heat stress" in info_str_lower:
            treatments_found.add("Heat Stress")
        if "cold stress" in info_str_lower:
            treatments_found.add("Cold Stress")
        if "chemical" in info_str_lower or "ppt" in info_str_lower or "herbicide" in info_str_lower or "pesticide" in info_str_lower:
            treatments_found.add("Chemical Stress")
        if "nutrient deficiency" in info_str_lower or "low nutrient" in info_str_lower:
            treatments_found.add("Nutrient Deficiency")
        if "biotic" in info_str_lower or "pathogen" in info_str_lower or "insect" in info_str_lower:
            treatments_found.add("Biotic Stress")
        if "low light" in info_str_lower:
            treatments_found.add("Low Light Stress")
        if "high light" in info_str_lower:
            treatments_found.add("High Light Stress")
        if "red light" in info_str_lower:
            treatments_found.add("Red Light Stress")
        elif "light stress" in info_str_lower: # General light stress if not specific
            treatments_found.add("Other Light Stress")
        # "Other stress" is not added unless explicitly found, as "No stress" is the default.

    if not treatments_found:
        extracted_data["treatment"] = ["No stress"]
    else:
        extracted_data["treatment"] = sorted(list(treatments_found)) # Sort for consistent output

    # --- Extract Medium ---
    medium_raw = _extract_info_from_list(description_list, "Medium:")
    if medium_raw:
        extracted_data["medium"] = medium_raw.strip()
    else:
        # Infer medium if tissue is whole_plant and not specified
        if extracted_data["tissue"] == "whole_plant":
            extracted_data["medium"] = "soil"
        else:
            extracted_data["medium"] = "unspecified"

    return extracted_data

def GSE115555_extractor(sample_metadata: dict) -> dict:
    def get_text(data, key):
        value = data.get(key)
        if value is None:
            return ""
        if isinstance(value, list):
            return " ".join(map(str, value)).lower()
        return str(value).lower()

    # --- Tissue Extraction ---
    extracted_tissue = "unknown"
    tissue_found_in_char = False
    for char in sample_metadata.get('characteristics_ch1', []):
        if 'tissue:' in char.lower():
            tissue_raw = char.split(':', 1)[1].strip().lower()
            tissue_found_in_char = True
            break

    tissue_mapping = {
        "root": "root",
        "leaf": "leaf",
        "flower": "flower",
        "shoot": "shoot",
        "rosette": "rosette",
        "bud": "bud",
        "whole plant": "whole_plant",
        "silique": "silique",
        "callus": "callus",
        "seed": "seed",
        "seedling": "seedling",
        "hypocotyl": "seedling",
        "embryo": "seed",
        "cotyledon": "seedling",
        "stem": "shoot",
        "apex": "shoot",
        "inflorescence": "flower",
        "cell culture": "callus",
    }

    if tissue_found_in_char:
        extracted_tissue = tissue_mapping.get(tissue_raw, "unknown")
    else:
        # If no explicit 'tissue:' found, look for other clues
        all_text_for_tissue = get_text(sample_metadata, 'characteristics_ch1') + " " + \
                              get_text(sample_metadata, 'growth_protocol_ch1') + " " + \
                              get_text(sample_metadata, 'title')

        if "seedling" in all_text_for_tissue:
            extracted_tissue = "seedling"
        elif "whole plant" in all_text_for_tissue:
            extracted_tissue = "whole_plant"
        elif "cell culture" in all_text_for_tissue:
            extracted_tissue = "callus"

    # --- Medium Extraction ---
    medium_text = get_text(sample_metadata, 'source_name_ch1') + " " + \
                  get_text(sample_metadata, 'growth_protocol_ch1') + " " + \
                  get_text(sample_metadata, 'characteristics_ch1')

    extracted_medium = "unspecified"
    if "solid media" in medium_text or "agar" in medium_text or "phytagel" in medium_text:
        extracted_medium = "solid media"
    elif "ms medium" in medium_text or "murashige skoog" in medium_text:
        extracted_medium = "MS medium"
    elif "liquid media" in medium_text or "hydroponic" in medium_text:
        extracted_medium = "liquid media"
    elif "soil" in medium_text:
        extracted_medium = "soil"
    elif "cell culture" in medium_text:
        extracted_medium = "cell culture medium"

    # --- Treatment Extraction ---
    treatment_keywords = {
        "Drought Stress": ["drought", "water deprivation", "dry condition"],
        "Dehidration Stress": ["dehydration"],
        "Salinity Stress": ["salinity", "salt stress", "nacl"],
        "Heat Stress": ["heat stress", "high temperature", "heat shock"],
        "Cold Stress": ["cold stress", "low temperature", "chilling", "freezing"],
        "Chemical Stress": ["chemical stress", "herbicide", "pesticide", "heavy metal", "cadmium", "arsenic", "aluminum", "paraquat", "osmotic stress", "mannitol", "sorbitol", "oxidative stress", "h2o2"],
        "Nutrient Deficiency": ["nutrient deficiency", "starvation", "low nitrogen", "low phosphate", "low potassium", "nitrogen deprivation", "phosphate deprivation"],
        "Biotic Stress": ["biotic stress", "pathogen", "insect", "fungal", "bacterial", "virus", "herbivory"],
        "Low Light Stress": ["low light stress", "darkness treatment", "shade treatment"],
        "High Light Stress": ["high light stress", "excess light", "uv-b stress"],
        "Red Light Stress": ["red light stress", "red light treatment"],
        "Other Light Stress": ["blue light stress", "blue light treatment", "far-red light stress", "uv-a stress"],
        "Other stress": ["stress", "mechanical stress", "wounding", "gravity", "microgravity", "spaceflight", "radiation", "bric hardware"],
    }

    all_text = get_text(sample_metadata, 'title') + " " + \
               get_text(sample_metadata, 'source_name_ch1') + " " + \
               get_text(sample_metadata, 'characteristics_ch1') + " " + \
               get_text(sample_metadata, 'treatment_protocol_ch1') + " " + \
               get_text(sample_metadata, 'growth_protocol_ch1') + " " + \
               get_text(sample_metadata, 'description') # Include description if it exists

    extracted_treatments = set()
    for treatment_type, keywords in treatment_keywords.items():
        for keyword in keywords:
            if keyword in all_text:
                extracted_treatments.add(treatment_type)
                break # Move to next treatment_type once a keyword is found

    # Handle "No stress"
    if not extracted_treatments:
        extracted_treatments.add("No stress")
    # Remove "No stress" if other stresses are found
    elif "No stress" in extracted_treatments and len(extracted_treatments) > 1:
        extracted_treatments.remove("No stress")

    return {
        "tissue": extracted_tissue,
        "treatment": sorted(list(extracted_treatments)),
        "medium": extracted_medium
    }

def GSE118364_extractor(sample_metadata: dict) -> dict:
    result = {}

    # Helper to safely get and join string values from metadata
    def _get_text(key):
        value = sample_metadata.get(key)
        if value is None:
            return ""
        if isinstance(value, list):
            return " ".join(str(item) for item in value).lower()
        return str(value).lower()

    # Combine relevant fields for easier searching
    all_text = (
        _get_text('title') + " " +
        _get_text('source_name_ch1') + " " +
        _get_text('characteristics_ch1') + " " +
        _get_text('treatment_protocol_ch1') + " " +
        _get_text('growth_protocol_ch1') + " " +
        _get_text('description')
    )

    # --- Extract Tissue ---
    tissue_found = "unknown"
    tissue_enums = ["root", "leaf", "flower", "shoot", "rosette", "bud", "whole_plant", "silique", "callus", "seed", "seedling"]

    # Prioritize characteristics_ch1 for developmental stage
    char_ch1 = _get_text('characteristics_ch1')
    if "developmental stage:" in char_ch1:
        # Extract the stage info, handling potential multiple characteristics
        stage_parts = char_ch1.split("developmental stage:")
        if len(stage_parts) > 1:
            stage_info = stage_parts[1].split(",")[0].strip()
            for tissue_enum in tissue_enums:
                if tissue_enum in stage_info:
                    tissue_found = tissue_enum
                    break
    
    # If not found, search in source_name_ch1
    if tissue_found == "unknown":
        source_name_ch1 = _get_text('source_name_ch1')
        for tissue_enum in tissue_enums:
            if tissue_enum in source_name_ch1:
                tissue_found = tissue_enum
                break
    
    # Fallback to general search in all text
    if tissue_found == "unknown":
        for tissue_enum in tissue_enums:
            if tissue_enum in all_text:
                tissue_found = tissue_enum
                break

    result['tissue'] = tissue_found

    # --- Extract Treatment ---
    treatments = set()
    treatment_protocol_text = _get_text('treatment_protocol_ch1')
    
    # Check for "No stress" explicitly
    if "no treatment" in treatment_protocol_text or "untreated" in treatment_protocol_text or "control" in treatment_protocol_text:
        treatments.add("No stress")
    else:
        # Look for specific stresses
        stress_mapping = {
            "drought": "Drought Stress",
            "dehydration": "Dehidration Stress",
            "salinity": "Salinity Stress",
            "salt stress": "Salinity Stress",
            "heat": "Heat Stress",
            "cold": "Cold Stress",
            "chemical": "Chemical Stress",
            "nutrient deficiency": "Nutrient Deficiency",
            "biotic": "Biotic Stress",
            "low light": "Low Light Stress",
            "high light": "High Light Stress",
            "red light": "Red Light Stress",
        }
        
        found_specific_stress = False
        for keyword, stress_type in stress_mapping.items():
            if keyword in all_text: # Search in all relevant text
                treatments.add(stress_type)
                found_specific_stress = True
        
        # If no specific stress and no "No stress" found yet
        if not found_specific_stress and not treatments:
            if treatment_protocol_text: # If treatment protocol exists but didn't match known stresses
                treatments.add("Other stress")
            else: # If no treatment info at all
                treatments.add("No stress")

    result['treatment'] = sorted(list(treatments)) # Sort for consistent output

    # --- Extract Medium ---
    medium_found = "unspecified"
    growth_protocol_text = _get_text('growth_protocol_ch1')

    if growth_protocol_text:
        # Try to extract the specific medium description
        extracted_medium_phrase = ""
        start_phrases = ["sowed in", "grown in", "cultured in", "on "]
        
        for sp in start_phrases:
            if sp in growth_protocol_text:
                start_idx = growth_protocol_text.find(sp) + len(sp)
                temp_phrase = growth_protocol_text[start_idx:].strip()
                
                # Try to cut off at common conjunctions or end of sentence
                end_phrases = [" and germinated", " and grown", " for ", ".", ","]
                for ep in end_phrases:
                    if ep in temp_phrase:
                        temp_phrase = temp_phrase.split(ep)[0].strip()
                        break
                
                # If the extracted phrase is somewhat descriptive and contains a common medium keyword
                medium_keywords_in_phrase = ["medium", "agar", "soil", "hydroponic", "vermiculite", "peat", "sand", "water"]
                if len(temp_phrase) > 5 and any(m_kw in temp_phrase for m_kw in medium_keywords_in_phrase):
                    extracted_medium_phrase = temp_phrase
                    break # Found a good phrase, stop searching
        
        if extracted_medium_phrase:
            medium_found = extracted_medium_phrase.replace("wth", "with") # Fix typo
        else:
            # Fallback if no specific phrase could be extracted, look for keywords
            if "ms medium" in growth_protocol_text or "murashige skoog" in growth_protocol_text:
                medium_found = "MS medium"
                if "agarized" in growth_protocol_text or "agar" in growth_protocol_text:
                    medium_found += " (agarized)"
                if "liquid" in growth_protocol_text:
                    medium_found += " (liquid)"
            elif "soil" in growth_protocol_text:
                medium_found = "soil"
            elif "hydroponic" in growth_protocol_text:
                medium_found = "hydroponic solution"
            elif "vermiculite" in growth_protocol_text:
                medium_found = "vermiculite"
            elif "peat" in growth_protocol_text:
                medium_found = "peat"
            elif "sand" in growth_protocol_text:
                medium_found = "sand"
            elif "agar" in growth_protocol_text:
                medium_found = "agar medium"
            elif "liquid" in growth_protocol_text and "medium" in growth_protocol_text:
                medium_found = "liquid medium"
            elif "water" in growth_protocol_text and "medium" not in growth_protocol_text:
                medium_found = "water"
            # If still nothing specific, it remains "unspecified"

    # If tissue is whole_plant and medium is unspecified, infer soil
    if result['tissue'] == "whole_plant" and medium_found == "unspecified":
        medium_found = "soil"

    result['medium'] = medium_found.replace("wth", "with") # Final cleanup for typo

    return result

def GSE61897_extractor(sample_metadata: dict) -> dict:
    extracted_data = {
        "tissue": "unknown",
        "treatment": ["No stress"],
        "medium": "unspecified"
    }

    # Helper to safely get the first element of a list or an empty string
    def get_first_or_empty(key):
        val = sample_metadata.get(key)
        if val is None:
            return ""
        if isinstance(val, list) and val:
            return val[0]
        if isinstance(val, str):
            return val
        return ""

    # --- Extract Tissue ---
    characteristics_ch1 = sample_metadata.get('characteristics_ch1', [])
    for char in characteristics_ch1:
        char_lower = char.lower()
        if 'tissue:' in char_lower:
            tissue_val = char_lower.split('tissue:')[1].strip()
            if 'aerial part' in tissue_val:
                extracted_data["tissue"] = "shoot"
            elif 'whole plant' in tissue_val:
                extracted_data["tissue"] = "whole_plant"
            elif 'seedling' in tissue_val:
                extracted_data["tissue"] = "seedling"
            elif 'seed' in tissue_val:
                extracted_data["tissue"] = "seed"
            elif 'root' in tissue_val:
                extracted_data["tissue"] = "root"
            elif 'leaf' in tissue_val:
                extracted_data["tissue"] = "leaf"
            elif 'flower' in tissue_val:
                extracted_data["tissue"] = "flower"
            elif 'silique' in tissue_val:
                extracted_data["tissue"] = "silique"
            elif 'callus' in tissue_val:
                extracted_data["tissue"] = "callus"
            elif 'rosette' in tissue_val:
                extracted_data["tissue"] = "rosette"
            elif 'bud' in tissue_val:
                extracted_data["tissue"] = "bud"
            # Once a tissue is identified, we can stop searching
            break
    
    # --- Extract Treatment ---
    treatments = set()
    
    # Aggregate relevant text fields for treatment search
    treatment_search_text_fields = []
    treatment_search_text_fields.append(get_first_or_empty('growth_protocol_ch1'))
    treatment_search_text_fields.extend(characteristics_ch1)
    treatment_search_text_fields.append(get_first_or_empty('title'))
    treatment_search_text_fields.append(get_first_or_empty('description'))

    for text_field in treatment_search_text_fields:
        text_field_lower = text_field.lower()
        if "drought" in text_field_lower:
            treatments.add("Drought Stress")
        if "dehydration" in text_field_lower:
            treatments.add("Dehidration Stress")
        if "salinity" in text_field_lower or "salt stress" in text_field_lower:
            treatments.add("Salinity Stress")
        if "heat stress" in text_field_lower or "high temperature" in text_field_lower:
            treatments.add("Heat Stress")
        if "cold stress" in text_field_lower or "low temperature" in text_field_lower:
            treatments.add("Cold Stress")
        if "chemical" in text_field_lower and "stress" in text_field_lower:
            treatments.add("Chemical Stress")
        if "nutrient deficiency" in text_field_lower or "low nutrient" in text_field_lower:
            treatments.add("Nutrient Deficiency")
        if "biotic stress" in text_field_lower or "pathogen" in text_field_lower or "insect" in text_field_lower:
            treatments.add("Biotic Stress")
        if "low light" in text_field_lower:
            treatments.add("Low Light Stress")
        if "high light" in text_field_lower:
            treatments.add("High Light Stress")
        if "red light" in text_field_lower:
            treatments.add("Red Light Stress")

    if not treatments:
        extracted_data["treatment"] = ["No stress"]
    else:
        extracted_data["treatment"] = sorted(list(treatments))

    # --- Extract Medium ---
    growth_protocol = get_first_or_empty('growth_protocol_ch1').lower()
    if "murashige and skoog" in growth_protocol or "ms salts" in growth_protocol:
        if "agar" in growth_protocol:
            extracted_data["medium"] = "MS agar"
        elif "liquid" in growth_protocol:
            extracted_data["medium"] = "MS liquid"
        else:
            extracted_data["medium"] = "MS medium"
    elif "agar plates" in growth_protocol or "agar medium" in growth_protocol:
        extracted_data["medium"] = "agar"
    elif "soil" in growth_protocol:
        extracted_data["medium"] = "soil"
    elif "hydroponic" in growth_protocol:
        extracted_data["medium"] = "hydroponic"
    # If no specific medium is found, it remains "unspecified"

    return extracted_data

import re

def GSE5521_extractor(sample_metadata: dict) -> dict:
    def _get_text(metadata_dict, key):
        value = metadata_dict.get(key)
        if isinstance(value, list):
            return " ".join(value).lower()
        return (value or "").lower()

    result = {
        "tissue": "unknown",
        "treatment": ["No stress"],
        "medium": "unspecified"
    }

    all_text = (
        _get_text(sample_metadata, "title") + " " +
        _get_text(sample_metadata, "source_name_ch1") + " " +
        _get_text(sample_metadata, "characteristics_ch1") + " " +
        _get_text(sample_metadata, "treatment_protocol_ch1") + " " +
        _get_text(sample_metadata, "description")
    )

    # --- Extract Tissue ---
    tissue_found = False
    tissue_keywords = {
        "rosette leaves": "rosette",
        "rosette": "rosette",
        "seedling": "seedling",
        "seedlings": "seedling",
        "leaf": "leaf",
        "leaves": "leaf",
        "root": "root",
        "roots": "root",
        "flower": "flower",
        "flowers": "flower",
        "shoot": "shoot",
        "shoots": "shoot",
        "bud": "bud",
        "buds": "bud",
        "whole plant": "whole_plant",
        "silique": "silique",
        "siliques": "silique",
        "callus": "callus",
        "seed": "seed",
        "seeds": "seed"
    }

    char_ch1_text = _get_text(sample_metadata, "characteristics_ch1")
    match = re.search(r"tissue:\s*([a-z\s]+)", char_ch1_text)
    if match:
        extracted_tissue = match.group(1).strip()
        for keyword, mapped_tissue in tissue_keywords.items():
            if keyword in extracted_tissue:
                result["tissue"] = mapped_tissue
                tissue_found = True
                break
    
    if not tissue_found:
        for keyword, mapped_tissue in tissue_keywords.items():
            if keyword in all_text:
                result["tissue"] = mapped_tissue
                tissue_found = True
                break

    # --- Extract Treatment ---
    treatments = set()
    
    if "no stress" in all_text or "control" in all_text or "untreated" in all_text:
        treatments.add("No stress")

    if "aba-treated" in all_text or "chemical" in all_text or "hormone" in all_text or "pesticide" in all_text or "herbicide" in all_text or "fungicide" in all_text:
        treatments.add("Chemical Stress")
    if "dehydration" in all_text:
        treatments.add("Dehidration Stress")
    if "drought" in all_text or "water deficit" in all_text:
        treatments.add("Drought Stress")
    if "salinity" in all_text or "salt" in all_text:
        treatments.add("Salinity Stress")
    if "heat" in all_text or "high temperature" in all_text:
        treatments.add("Heat Stress")
    if "cold" in all_text or "low temperature" in all_text:
        treatments.add("Cold Stress")
    if "nutrient deficiency" in all_text or "nitrogen starvation" in all_text or "phosphate starvation" in all_text:
        treatments.add("Nutrient Deficiency")
    if "biotic" in all_text or "pathogen" in all_text or "insect" in all_text or "fungal infection" in all_text or "bacterial infection" in all_text:
        treatments.add("Biotic Stress")
    if "low light" in all_text or "dark" in all_text or "shade" in all_text:
        treatments.add("Low Light Stress")
    if "high light" in all_text or "uv light" in all_text:
        treatments.add("High Light Stress")
    if "red light" in all_text:
        treatments.add("Red Light Stress")
    elif "light stress" in all_text:
        treatments.add("Other Light Stress")
    
    if len(treatments) > 1 and "No stress" in treatments:
        treatments.remove("No stress")
    
    if treatments:
        result["treatment"] = sorted(list(treatments))
    else:
        result["treatment"] = ["No stress"]

    # --- Extract Medium ---
    medium_found = False
    
    desc_text = _get_text(sample_metadata, "description")
    match = re.search(r"growth substrate:\s*([a-z0-9\s\-\(\):,./]+)", desc_text)
    if match:
        extracted_medium = match.group(1).strip()
        if "peat-vermiculite" in extracted_medium or "peat and vermiculite" in extracted_medium:
            result["medium"] = "peat-vermiculite mixture"
            medium_found = True
        elif "soil" in extracted_medium:
            result["medium"] = "soil"
            medium_found = True
        elif "ms media" in extracted_medium or "murashige skoog" in extracted_medium:
            result["medium"] = "MS media"
            medium_found = True
        elif "agar" in extracted_medium:
            result["medium"] = "agar"
            medium_found = True
    
    if not medium_found:
        if "soil-grown plant" in char_ch1_text:
            result["medium"] = "soil"
            medium_found = True

    if not medium_found:
        if "peat-vermiculite" in all_text or "peat and vermiculite" in all_text:
            result["medium"] = "peat-vermiculite mixture"
            medium_found = True
        elif "ms media" in all_text or "murashige skoog" in all_text:
            result["medium"] = "MS media"
            medium_found = True
        elif "agar" in all_text:
            result["medium"] = "agar"
            medium_found = True
        elif "soil" in all_text:
            result["medium"] = "soil"
            medium_found = True
        elif "hydroponic" in all_text:
            result["medium"] = "hydroponic solution"
            medium_found = True
        elif "liquid culture" in all_text:
            result["medium"] = "liquid culture"
            medium_found = True
        elif "gelrite" in all_text:
            result["medium"] = "gelrite"
            medium_found = True
        elif "water" in all_text:
            result["medium"] = "water"
            medium_found = True
    
    return result

import re

def GSE52013_extractor(sample_metadata: dict) -> dict:
    extracted_data = {
        "tissue": "unknown",
        "treatment": [],
        "medium": "unspecified"
    }

    def get_first_or_none(data, key):
        if key in data and data[key]:
            return data[key][0]
        return None

    def get_joined_or_empty(data, key, separator=' '):
        if key in data and data[key]:
            return separator.join(data[key])
        return ""

    # --- Extract Tissue ---
    characteristics_ch1 = sample_metadata.get('characteristics_ch1', [])
    for char in characteristics_ch1:
        if 'tissue:' in char:
            extracted_tissue = char.split('tissue:')[1].strip().lower()
            if 'seedling' in extracted_tissue:
                extracted_data["tissue"] = "seedling"
            elif 'whole plant' in extracted_tissue:
                extracted_data["tissue"] = "whole_plant"
            elif 'root' in extracted_tissue:
                extracted_data["tissue"] = "root"
            elif 'leaf' in extracted_tissue:
                extracted_data["tissue"] = "leaf"
            elif 'flower' in extracted_tissue:
                extracted_data["tissue"] = "flower"
            elif 'shoot' in extracted_tissue:
                extracted_data["tissue"] = "shoot"
            elif 'rosette' in extracted_tissue:
                extracted_data["tissue"] = "rosette"
            elif 'bud' in extracted_tissue:
                extracted_data["tissue"] = "bud"
            elif 'silique' in extracted_tissue:
                extracted_data["tissue"] = "silique"
            elif 'callus' in extracted_tissue:
                extracted_data["tissue"] = "callus"
            elif 'seed' in extracted_tissue:
                extracted_data["tissue"] = "seed"
            break

    if extracted_data["tissue"] == "unknown":
        source_name = get_first_or_none(sample_metadata, 'source_name_ch1')
        description = get_first_or_none(sample_metadata, 'description')
        if source_name and 'seedling' in source_name.lower():
            extracted_data["tissue"] = "seedling"
        elif description and 'seedling' in description.lower():
            extracted_data["tissue"] = "seedling"
        elif source_name and 'plant' in source_name.lower():
            extracted_data["tissue"] = "whole_plant"
        elif description and 'plant' in description.lower():
            extracted_data["tissue"] = "whole_plant"

    # --- Extract Treatment ---
    treatments_set = set()
    
    combined_text = (
        get_joined_or_empty(sample_metadata, 'title').lower() + " " +
        get_joined_or_empty(sample_metadata, 'source_name_ch1').lower() + " " +
        get_joined_or_empty(sample_metadata, 'characteristics_ch1').lower() + " " +
        get_joined_or_empty(sample_metadata, 'treatment_protocol_ch1').lower() + " " +
        get_joined_or_empty(sample_metadata, 'description').lower()
    )

    if 'drought' in combined_text or 'dehydration' in combined_text:
        treatments_set.add("Drought Stress")
    if 'salinity' in combined_text or 'salt stress' in combined_text:
        treatments_set.add("Salinity Stress")
    if 'heat stress' in combined_text or 'high temperature' in combined_text:
        treatments_set.add("Heat Stress")
    if 'cold stress' in combined_text or 'low temperature' in combined_text:
        treatments_set.add("Cold Stress")
    
    if 'chemical' in combined_text or 'μm' in combined_text or 'compound' in combined_text or 'drug' in combined_text or 'treated with' in combined_text:
        if 'chemically treated' in combined_text or re.search(r'\d+\s*μm', combined_text) or ('dmso' in combined_text and 'control' not in combined_text):
            treatments_set.add("Chemical Stress")
    
    if 'nutrient deficiency' in combined_text or 'low nutrient' in combined_text:
        treatments_set.add("Nutrient Deficiency")
    if 'biotic stress' in combined_text or 'pathogen' in combined_text or 'insect' in combined_text:
        treatments_set.add("Biotic Stress")
    if 'low light' in combined_text:
        treatments_set.add("Low Light Stress")
    if 'high light' in combined_text:
        treatments_set.add("High Light Stress")
    if 'red light' in combined_text:
        treatments_set.add("Red Light Stress")
    if 'light stress' in combined_text and not any(s in combined_text for s in ['low light', 'high light', 'red light']):
        treatments_set.add("Other Light Stress")
    
    if not treatments_set:
        if 'control' in combined_text or 'untreated' in combined_text or 'no stress' in combined_text:
            treatments_set.add("No stress")
        elif 'treated' in combined_text or 'treatment' in combined_text:
            treatments_set.add("Other stress")
        else:
            treatments_set.add("No stress")
    
    if len(treatments_set) > 1 and "No stress" in treatments_set:
        treatments_set.remove("No stress")

    extracted_data["treatment"] = sorted(list(treatments_set))

    # --- Extract Medium ---
    growth_protocol_ch1 = get_joined_or_empty(sample_metadata, 'growth_protocol_ch1')
    if growth_protocol_ch1:
        growth_protocol_lower = growth_protocol_ch1.lower()
        if 'murashige and skoog medium' in growth_protocol_lower or 'ms medium' in growth_protocol_lower:
            extracted_data["medium"] = "Murashige and Skoog medium"
        elif 'agar' in growth_protocol_lower:
            extracted_data["medium"] = "Agar"
        elif 'soil' in growth_protocol_lower:
            extracted_data["medium"] = "soil"

    return extracted_data

def GSE105058_extractor(sample_metadata: dict) -> dict:
    extracted_data = {}

    # --- Extract Tissue ---
    tissue_found = "unknown"
    # Prioritize more specific tissues. Order matters here.
    tissue_priority = [
        ("shoot", ["shoots"]),
        ("root", ["root"]),
        ("leaf", ["leaf"]),
        ("flower", ["flower"]),
        ("rosette", ["rosette"]),
        ("bud", ["bud"]),
        ("silique", ["silique"]),
        ("callus", ["callus"]),
        ("seed", ["seed"]),
        ("whole_plant", ["whole plant", "whole-plant"]),
        ("seedling", ["seedling", "seedlings"]),
    ]

    if 'characteristics_ch1' in sample_metadata and isinstance(sample_metadata['characteristics_ch1'], list):
        for char_str in sample_metadata['characteristics_ch1']:
            if "tissue:" in char_str.lower():
                tissue_desc = char_str.split("tissue:", 1)[1].strip().lower()
                for mapped_value, keywords in tissue_priority:
                    if any(kw in tissue_desc for kw in keywords):
                        tissue_found = mapped_value
                        break # Found the most specific tissue, break from inner loop
            if tissue_found != "unknown":
                break # Found a tissue, break from outer loop

    extracted_data['tissue'] = tissue_found

    # --- Extract Treatment ---
    treatments = set()
    
    # Specific treatments
    specific_treatment_keywords = {
        "drought": "Drought Stress",
        "dehydration": "Dehidration Stress",
        "salinity": "Salinity Stress",
        "heat": "Heat Stress",
        "cold": "Cold Stress",
        "chemical": "Chemical Stress",
        "nutrient deficiency": "Nutrient Deficiency",
        "biotic": "Biotic Stress",
        "low light": "Low Light Stress",
        "high light": "High Light Stress",
        "red light": "Red Light Stress",
        "light stress": "Other Light Stress",
    }
    
    # Explicit "other stress" keywords (e.g., physical stressors)
    other_stress_explicit_keywords = {
        "microgravity": "Other stress",
        "spaceflight": "Other stress",
    }

    search_fields = []
    if 'characteristics_ch1' in sample_metadata and isinstance(sample_metadata['characteristics_ch1'], list):
        search_fields.extend(sample_metadata['characteristics_ch1'])
    if 'title' in sample_metadata and isinstance(sample_metadata['title'], list):
        search_fields.extend(sample_metadata['title'])
    if 'description' in sample_metadata and isinstance(sample_metadata['description'], list):
        search_fields.extend(sample_metadata['description'])
    if 'growth_protocol_ch1' in sample_metadata and isinstance(sample_metadata['growth_protocol_ch1'], list):
        search_fields.extend(sample_metadata['growth_protocol_ch1'])

    search_text = " ".join(search_fields).lower()

    found_any_specific_or_explicit_stress = False

    # Check for specific treatments
    for keyword, mapped_value in specific_treatment_keywords.items():
        # Use regex for whole word matching to avoid false positives (e.g., "chemical" in "biochemical")
        if re.search(r'\b' + re.escape(keyword) + r'\b', search_text):
            treatments.add(mapped_value)
            found_any_specific_or_explicit_stress = True

    # Check for explicit "other stress" keywords
    for keyword, mapped_value in other_stress_explicit_keywords.items():
        if re.search(r'\b' + re.escape(keyword) + r'\b', search_text):
            treatments.add(mapped_value)
            found_any_specific_or_explicit_stress = True

    # If no specific or explicit "other" stress was found, check for general "stress"
    if not found_any_specific_or_explicit_stress:
        if re.search(r'\bstress\b', search_text):
            treatments.add("Other stress")
    
    # If no treatments were identified at all, default to "No stress"
    if not treatments:
        treatments.add("No stress")

    extracted_data['treatment'] = sorted(list(treatments))

    # --- Extract Medium ---
    medium_found = "unspecified"
    medium_keywords = {
        "agar": "agar",
        "ms medium": "MS medium",
        "murashige skoog": "MS medium",
        "soil": "soil",
        "hydroponic": "hydroponic",
        "liquid medium": "liquid medium",
        "vermiculite": "vermiculite",
        "perlite": "perlite",
        "rockwool": "rockwool",
        "gelrite": "gelrite",
    }

    medium_search_text = ""
    if 'growth_protocol_ch1' in sample_metadata and isinstance(sample_metadata['growth_protocol_ch1'], list):
        medium_search_text += " ".join(sample_metadata['growth_protocol_ch1']).lower()
    if 'characteristics_ch1' in sample_metadata and isinstance(sample_metadata['characteristics_ch1'], list):
        medium_search_text += " ".join(sample_metadata['characteristics_ch1']).lower()

    for keyword, mapped_value in medium_keywords.items():
        if re.search(r'\b' + re.escape(keyword) + r'\b', medium_search_text):
            medium_found = mapped_value
            break
    
    # Special case: if "whole plant" tissue and no medium specified, infer "soil"
    if medium_found == "unspecified" and extracted_data['tissue'] == "whole_plant":
        medium_found = "soil"

    extracted_data['medium'] = medium_found

    return extracted_data

def GSE26973_extractor(sample_metadata: dict) -> dict:
    extracted_data = {
        "tissue": "unknown",
        "treatment": [],
        "medium": "unspecified"
    }

    # Access the inner sample_metadata and study_metadata
    inner_sample_metadata = sample_metadata.get('sample_metadata', {})
    study_metadata = sample_metadata.get('study_metadata', {})

    # --- Extract Tissue ---
    characteristics_ch1 = inner_sample_metadata.get('characteristics_ch1', [])
    for char_str in characteristics_ch1:
        if char_str.startswith('tissue:'):
            tissue_raw = char_str.split(':', 1)[1].strip().lower()
            # Map raw tissue string to schema enum values
            if 'leaf' in tissue_raw:
                extracted_data['tissue'] = 'leaf'
            elif 'root' in tissue_raw:
                extracted_data['tissue'] = 'root'
            elif 'flower' in tissue_raw:
                extracted_data['tissue'] = 'flower'
            elif 'shoot' in tissue_raw:
                extracted_data['tissue'] = 'shoot'
            elif 'rosette' in tissue_raw:
                extracted_data['tissue'] = 'rosette'
            elif 'bud' in tissue_raw:
                extracted_data['tissue'] = 'bud'
            elif 'whole plant' in tissue_raw or 'whole_plant' in tissue_raw:
                extracted_data['tissue'] = 'whole_plant'
            elif 'silique' in tissue_raw:
                extracted_data['tissue'] = 'silique'
            elif 'callus' in tissue_raw:
                extracted_data['tissue'] = 'callus'
            elif 'seed' in tissue_raw:
                extracted_data['tissue'] = 'seed'
            elif 'seedling' in tissue_raw:
                extracted_data['tissue'] = 'seedling'
            else:
                extracted_data['tissue'] = 'unknown' # Default if no specific match
            break # Stop after finding the first tissue entry

    # --- Extract Treatment ---
    treatments_set = set() # Use a set to ensure unique treatments

    # Check characteristics_ch1 for 'agent'
    for char_str in characteristics_ch1:
        if char_str.startswith('agent:'):
            agent = char_str.split(':', 1)[1].strip().lower()
            if 'g3p' in agent or 'exudate' in agent: # G3P and Exudate are chemical agents
                treatments_set.add('Chemical Stress')
            # Add more agent-based treatments here if applicable
            # e.g., if 'drought' in agent: treatments_set.add('Drought Stress')

    # Check treatment_protocol_ch1
    treatment_protocol = inner_sample_metadata.get('treatment_protocol_ch1', '').lower()
    if 'g3p' in treatment_protocol or 'exudate' in treatment_protocol:
        treatments_set.add('Chemical Stress')
    if 'drought' in treatment_protocol or 'dehydration' in treatment_protocol:
        treatments_set.add('Drought Stress')
    if 'salinity' in treatment_protocol or 'salt' in treatment_protocol:
        treatments_set.add('Salinity Stress')
    if 'heat' in treatment_protocol:
        treatments_set.add('Heat Stress')
    if 'cold' in treatment_protocol:
        treatments_set.add('Cold Stress')
    if 'nutrient deficiency' in treatment_protocol:
        treatments_set.add('Nutrient Deficiency')
    if 'biotic stress' in treatment_protocol or 'pathogen' in treatment_protocol or 'insect' in treatment_protocol:
        treatments_set.add('Biotic Stress')
    # Add more checks for other stress types from the enum as needed

    # Check description
    description = inner_sample_metadata.get('description', '').lower()
    if 'g3p' in description or 'exudate' in description:
        treatments_set.add('Chemical Stress')
    # Add similar checks for other stress types in description as needed
    
    extracted_data['treatment'] = sorted(list(treatments_set)) # Sort for consistent output

    # --- Extract Medium ---
    # Check growth_protocol_ch1 for explicit medium
    growth_protocol = inner_sample_metadata.get('growth_protocol_ch1', '').lower()
    if 'soil' in growth_protocol:
        extracted_data['medium'] = 'soil'
    elif 'agar' in growth_protocol or 'plate' in growth_protocol:
        extracted_data['medium'] = 'agar'
    elif 'hydroponic' in growth_protocol:
        extracted_data['medium'] = 'hydroponic'
    elif 'liquid' in growth_protocol or ('medium' in growth_protocol and 'ms' in growth_protocol): # e.g., MS liquid medium
        extracted_data['medium'] = 'liquid'
    # Add more specific medium checks if needed

    # Infer 'soil' if context implies whole plants and medium is still unspecified
    if extracted_data['medium'] == 'unspecified':
        study_summary = study_metadata.get('summary', [])
        study_summary_text = ' '.join(study_summary).lower()
        # Look for keywords implying whole plants, typically grown in soil
        if ('arabidopsis' in study_summary_text or 'plant' in study_summary_text) and \
           ('grown' in study_summary_text or 'cultivated' in study_summary_text):
            extracted_data['medium'] = 'soil'
    
    return extracted_data

def GSE44781_extractor(sample_metadata: dict) -> dict:
    """
    Extracts tissue, treatment, and medium information from a sample metadata dictionary
    conforming to the GSE44781 structure.

    Args:
        sample_metadata (dict): A dictionary containing sample and study metadata.

    Returns:
        dict: A dictionary formatted according to the specified schema, containing
              'tissue', 'treatment', and 'medium' information.
    """
    extracted_data = {
        "tissue": "unknown",
        "treatment": ["No stress"],
        "medium": "unspecified"
    }

    # Helper to safely get a value from a list or string, handling potential missing keys
    def _get_value(data, key, default=None):
        if key in data:
            value = data[key]
            if isinstance(value, list):
                return value[0] if value else default
            return value
        return default

    # Access relevant parts of the input dictionary
    sample_info = sample_metadata.get("sample_metadata", {})
    study_info = sample_metadata.get("study_metadata", {})

    # --- Extract Tissue ---
    characteristics_ch1 = sample_info.get("characteristics_ch1", [])
    for char_str in characteristics_ch1:
        if char_str.lower().startswith("tissue:"):
            tissue_raw = char_str.split(":", 1)[1].strip().lower()
            # Map raw tissue to schema enum values
            if "root" in tissue_raw:
                extracted_data["tissue"] = "root"
            elif "leaf" in tissue_raw:
                extracted_data["tissue"] = "leaf"
            elif "flower" in tissue_raw:
                extracted_data["tissue"] = "flower"
            elif "shoot" in tissue_raw:
                extracted_data["tissue"] = "shoot"
            elif "rosette" in tissue_raw:
                extracted_data["tissue"] = "rosette"
            elif "bud" in tissue_raw:
                extracted_data["tissue"] = "bud"
            elif "whole plant" in tissue_raw or "whole_plant" in tissue_raw:
                extracted_data["tissue"] = "whole_plant"
            elif "silique" in tissue_raw:
                extracted_data["tissue"] = "silique"
            elif "callus" in tissue_raw:
                extracted_data["tissue"] = "callus"
            elif "seedling" in tissue_raw:
                extracted_data["tissue"] = "seedling"
            elif "seed" in tissue_raw:
                extracted_data["tissue"] = "seed"
            # "Secondary meristem" from the example does not directly map to the enum,
            # so it will remain "unknown" as initialized.
            break # Assuming only one primary tissue entry

    # --- Extract Treatment ---
    treatments = set()

    # Check characteristics_ch1 for stress information
    for char_str in characteristics_ch1:
        if char_str.lower().startswith("stress:"):
            stress_raw = char_str.split(":", 1)[1].strip().lower()
            if "clipped" in stress_raw or "herbivore" in stress_raw:
                treatments.add("Biotic Stress")
            # Add other specific stress mappings here if needed for other examples

    # Check treatment_protocol_ch1
    treatment_protocol = _get_value(sample_info, "treatment_protocol_ch1", "").lower()
    if "clipped" in treatment_protocol or "herbivore damage" in treatment_protocol:
        treatments.add("Biotic Stress")

    # Check overall_design for additional context on treatments
    overall_design = _get_value(study_info, "overall_design", "").lower()
    if "clipped" in overall_design or "herbivore damage" in overall_design:
        treatments.add("Biotic Stress")
    # Note: "cold treated seeds" in growth_protocol_ch1 or overall_design is considered
    # a pre-treatment for seeds, not a treatment applied to the plant sample itself,
    # so it's not added as "Cold Stress" for a plant sample.

    if treatments:
        extracted_data["treatment"] = sorted(list(treatments))
    else:
        extracted_data["treatment"] = ["No stress"]

    # --- Extract Medium ---
    growth_protocol = _get_value(sample_info, "growth_protocol_ch1", "").lower()
    if "soil" in growth_protocol or "plant mix" in growth_protocol or "pot" in growth_protocol:
        extracted_data["medium"] = "soil"
    elif "hydroponic" in growth_protocol:
        extracted_data["medium"] = "hydroponic"
    elif "agar" in growth_protocol or "plate" in growth_protocol:
        extracted_data["medium"] = "agar"
    # Default is "unspecified" as per initialization if no match is found

    return extracted_data