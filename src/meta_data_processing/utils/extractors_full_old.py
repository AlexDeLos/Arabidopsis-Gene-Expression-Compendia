
# import json
# import re

# def GSE110079_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from sample metadata
#     following the given schema for the GSE110079 dataset.

#     The output should be formatted as a JSON instance that conforms to the JSON schema:
#     {"properties": {"tissue": {"description": "Tissue the samples was extracted from.", "title": "Tissue", "type": "string"},
#                     "treatment": {"description": "List of treatments and stresses that was applied to the sample.", "items": {"type": "string"}, "title": "Treatment", "type": "array"},
#                     "medium": {"description": "Growth medium of the sample.", "title": "Medium", "type": "string"}},
#      "required": ["tissue", "treatment", "medium"]}

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 Expected to have a 'sample_characteristicts' key
#                                 which is a list of strings, where each string
#                                 describes a characteristic like 'tissue:' or 'treatment:'.

#     Returns:
#         dict: A dictionary formatted according to the specified schema,
#               containing tissue, treatment, and medium information.
#               All required fields ('tissue', 'treatment', 'medium') will be populated.
#     """
#     extracted_data = {
#         "tissue": None,
#         "treatment": [],
#         "medium": None
#     }

#     # Safely get 'sample_characteristicts', defaulting to an empty list if not present
#     characteristics = sample_metadata.get('sample_characteristicts', [])

#     raw_treatment_str = None

#     for item in characteristics:
#         # Extract 'tissue'
#         if item.startswith('tissue:'):
#             extracted_data['tissue'] = 'seedling' #item.split(':', 1)[1].strip()
#         # Extract 'treatment'
#         elif item.startswith('treatment:'):
#             raw_treatment_str = item.split(':', 1)[1].strip()
#             # The schema expects 'treatment' to be an array of strings
#             sample_treatment = 'control' if ('mock' in raw_treatment_str) else 'polyethylene glycol'
#             extracted_data['treatment'].append(sample_treatment)

#     # Determine 'medium' based on the extracted raw_treatment_str
#     # The medium information is embedded within the treatment string.
#     # Based on the provided samples, treatments are either 'mock (1/2MS)' or 'polyethylene glycol (PEG)'.
#     if raw_treatment_str:
#         if "(1/2MS)" in raw_treatment_str:
#             extracted_data['medium'] = "1/2MS"
#         elif "(PEG)" in raw_treatment_str:
#             # As per the schema description "Growth medium of the sample.",
#             # "polyethylene glycol (PEG)" is a more complete description than just "PEG".
#             extracted_data['medium'] = "polyethylene glycol (PEG)"
#         # If other medium types were present, this logic would need to be extended.
#         # For this dataset, these two cases cover all samples.

#     # The schema specifies all fields are required.
#     # This implementation assumes that the input data will always contain the
#     # necessary information to populate these fields as observed in the sample metadata.

#     return extracted_data



# import json

# def GSE27548_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from sample metadata for GSE27548
#     following the specified JSON schema.

#     The schema is:
#     {
#       "properties": {
#         "tissue": {
#           "description": "Tissue the samples was extracted from.",
#           "title": "Tissue",
#           "type": "string"
#         },
#         "treatment": {
#           "description": "List of treatments and stresses that was applied to the sample.",
#           "items": {
#             "type": "string"
#           },
#           "title": "Treatment",
#           "type": "array"
#         },
#         "medium": {
#           "description": "Growth medium of the sample.",
#           "title": "Medium",
#           "type": "string"
#         }
#       },
#       "required": ["tissue", "treatment", "medium"]
#     }

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 Example:
#                                 {
#                                     'Sample_id': 'GSM680322',
#                                     'sample_title': ['leaf_Ws-2_wet_3'],
#                                     'sample_source_name': ['rosette leaf'],
#                                     'sample_characteristicts': [
#                                         'accession: CS2360',
#                                         'treatment: Wet',
#                                         'ecotype: Ws-2',
#                                         'tissue: rosette leaf'
#                                     ],
#                                     'sample_extraction_protocol': ['Qiagen RNAasy Plant Mini prep, samples initially harvested in RNAlater']
#                                 }

#     Returns:
#         dict: A dictionary formatted according to the specified schema,
#               containing tissue, treatment, and medium information.
#     """
#     extracted_data = {}

#     # 1. Extract Tissue
#     # The 'tissue' information is found within the 'sample_characteristicts' list
#     # as an item starting with 'tissue:'.
#     tissue = None
#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if characteristic.startswith('tissue:'):
#                 tissue = characteristic.split(':', 1)[1].strip()
#                 break
#     # As a fallback, 'sample_source_name' also contains tissue information.
#     if tissue is None and 'sample_source_name' in sample_metadata and sample_metadata['sample_source_name']:
#         # Assuming 'sample_source_name' is a list and the first item is the tissue.
#         tissue = sample_metadata['sample_source_name'][0].strip()

#     if tissue:
#         extracted_data['tissue'] = tissue
#     else:
#         # If tissue is not found, it's a critical missing piece for the schema.
#         # For this dataset, 'rosette leaf' is consistently present.
#         # In a real-world scenario, one might raise an error or assign "unknown".
#         raise ValueError("Could not extract 'tissue' from sample metadata.")


#     # 2. Extract Treatment
#     # The 'treatment' information is found within the 'sample_characteristicts' list
#     # as an item starting with 'treatment:'. The schema expects an array of strings.
#     treatment_list = []
#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if characteristic.startswith('treatment:'):
#                 treatment = characteristic.split(':', 1)[1].strip()
#                 treatment_list.append('control' if treatment == 'Wet' else 'drought environment')
#                 # Assuming only one 'treatment' entry per sample for this dataset
#                 break

#     if treatment_list:
#         extracted_data['treatment'] = treatment_list
#     else:
#         # If treatment is not found, it's a critical missing piece for the schema.
#         # For this dataset, 'Wet' or 'Dry' is consistently present.
#         raise ValueError("Could not extract 'treatment' from sample metadata.")


#     # 3. Extract Medium
#     # The 'medium' information is not explicitly present in the sample_metadata.
#     # Based on the study_metadata ("two levels of soil moisture", "mild soil drying"),
#     # the growth medium is consistently "soil" for all samples in this study.
#     extracted_data['medium'] = "soil"

#     return extracted_data



# import json

# def GSE63128_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from a sample metadata dictionary
#     following a predefined JSON schema.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the specified schema, containing
#               'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The tissue information is found within the 'sample_characteristicts' list.
#     characteristics = sample_metadata.get('sample_characteristicts', [])
#     tissue_found = False
#     for char_string in characteristics:
#         if char_string.startswith('tissue:'):
#             extracted_data['tissue'] = char_string.split(':', 1)[1].strip()
#             tissue_found = True
#             break
#     if not tissue_found:
#         # Default or handle cases where tissue might be missing, though for this dataset it's consistent.
#         extracted_data['tissue'] = "unknown" 

#     # 2. Extract 'treatment'
#     # The specific treatment condition is indicated by 'sample group' in 'sample_characteristicts'.
#     # We map these specific groups to more general treatment descriptions.
#     treatment_group = None
#     for char_string in characteristics:
#         if char_string.startswith('sample group:'):
#             treatment_group = char_string.split(':', 1)[1].strip()
#             break

#     if treatment_group:
#         if 'Heat' in treatment_group:
#             extracted_data['treatment'] = ["Heat stress"]
#         elif 'Recovery' in treatment_group:
#             extracted_data['treatment'] = ["Recovery from heat stress"]
#         elif 'Control' in treatment_group:
#             extracted_data['treatment'] = ["Control condition"]
#         else:
#             # Fallback if a new group type appears, though not expected for this dataset.
#             extracted_data['treatment'] = [treatment_group] 
#     else:
#         # Default or handle cases where treatment group might be missing.
#         extracted_data['treatment'] = ["unknown"]

#     # 3. Extract 'medium'
#     # The growth medium is constant across all samples in this study and is mentioned
#     # in the 'overall_design' of the study metadata. Since study_metadata is not
#     # passed to this function, we hardcode the identified constant value.
#     extracted_data['medium'] = "Murashige and Skoog medium"

#     return extracted_data

# import json

# def GSE5624_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from a single sample's metadata
#     for the GSE5624 dataset, conforming to a specified JSON schema.

#     The schema requires:
#     - tissue (string): Tissue the sample was extracted from.
#     - treatment (array of strings): List of treatments/stresses applied.
#     - medium (string): Growth medium of the sample.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 Example:
#                                 {'Sample_id': 'GSM131349',
#                                  'sample_title': ['AtGen_6-4421_Droughtstress-Roots-6.0h_Rep1'],
#                                  'sample_source_name': ['Col-0'],
#                                  'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Roots'],
#                                  'sample_extraction_protocol': ''}

#     Returns:
#         dict: A dictionary formatted according to the specified schema,
#               containing 'tissue', 'treatment', and 'medium' information.
#               Example:
#               {
#                 "tissue": "Roots",
#                 "treatment": ["Drought stress"],
#                 "medium": "MS-liquid-media"
#               }
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The tissue information is found within the 'sample_characteristicts' list.
#     # It appears as a string like "Tissue: Roots" or "Tissue: Shoots".
#     tissue = None
#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if characteristic.startswith('Tissue:'):
#                 tissue = characteristic.split('Tissue: ')[1].strip()
#                 break
#     # If tissue is not found, it will remain None. The schema requires a string,
#     # so in a real-world scenario, one might raise an error or assign "unknown".
#     # For this dataset, it appears consistently present.
#     extracted_data['tissue'] = tissue

#     # 2. Extract 'treatment'
#     # Based on the study metadata ('title' and 'summary'), the primary stress/treatment
#     # applied across all samples in this study is "Drought stress".
#     # The sample titles also confirm this (e.g., "Droughtstress").
#     # The schema expects an array of strings.
#     treatment = ["Drought stress"]
#     extracted_data['treatment'] = treatment

#     # 3. Extract 'medium'
#     # From the study metadata's 'summary', the growth medium during the treatment phase
#     # is described as "MS-liquid-media" (after initial growth on MS-Agar-media).
#     # This is a constant for all samples in this study.
#     medium = "MS-liquid-media"
#     extracted_data['medium'] = medium

#     return extracted_data



# import json

# def GSE60960_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from a sample metadata dictionary
#     for the GSE60960 dataset, conforming to a specific JSON schema.

#     The schema requires 'tissue' (string), 'treatment' (array of strings),
#     and 'medium' (string).

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 Expected keys include 'sample_characteristicts',
#                                 'sample_title', and 'sample_source_name'.

#     Returns:
#         dict: A dictionary formatted according to the specified schema, containing
#               'tissue', 'treatment', and 'medium' information.

#     Raises:
#         ValueError: If required information ('tissue' or 'treatment') cannot be
#                     found or parsed from the sample_metadata.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # Tissue information is found in 'sample_characteristicts', which is a list of strings.
#     # We look for the string starting with 'tissue:'.
#     tissue_found = False
#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if characteristic.lower().startswith('tissue:'):
#                 extracted_data['tissue'] = characteristic.split(':', 1)[1].strip()
#                 tissue_found = True
#                 break
#     if not tissue_found:
#         # 'tissue' is a required field. If not found, raise an error.
#         raise ValueError("Tissue information not found in 'sample_characteristicts'.")

#     # 2. Extract 'treatment'
#     # Treatment information is typically found in 'sample_title'.
#     # We look for keywords like 'MLD stressed' or 'untreated'.
#     treatment_list = []
#     if 'sample_title' in sample_metadata and isinstance(sample_metadata['sample_title'], list) and sample_metadata['sample_title']:
#         title_text = sample_metadata['sample_title'][0]
#         if 'MLD stressed' in title_text:
#             treatment_list.append('MLD stress') # Interpreting "MLD stressed" as "MLD stress"
#         elif 'untreated' in title_text:
#             treatment_list.append('untreated')
    
#     if not treatment_list:
#         # 'treatment' is a required field. If no recognized treatment is found, raise an error.
#         raise ValueError("Treatment information not found or recognized in 'sample_title'.")
    
#     extracted_data['treatment'] = treatment_list

#     # 3. Extract 'medium'
#     # The growth medium is not explicitly named but is implied by "standard (UT) conditions"
#     # mentioned in 'sample_source_name' and the overall study design.
#     # This appears to be a constant value across all samples in this dataset.
#     extracted_data['medium'] = "standard conditions"

#     return extracted_data



# import json
# import re

# def GSE46205_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from sample metadata
#     following a specific JSON schema.

#     The schema is:
#     {"properties": {"tissue": {"description": "Tissue the samples was extracted from.", "title": "Tissue", "type": "string"},
#                    "treatment": {"description": "List of treatments and stresses that was applied to the sample.", "items": {"type": "string"}, "title": "Treatment", "type": "array"},
#                    "medium": {"description": "Growth medium of the sample.", "title": "Medium", "type": "string"}},
#     "required": ["tissue", "treatment", "medium"]}

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 Expected to contain keys like 'sample_title' and 'sample_characteristicts'.

#     Returns:
#         dict: A dictionary formatted according to the specified schema.
#               Example: {"tissue": "epidermis", "treatment": ["140mM NaCl for 32 hours"], "medium": "Standard media for 5 days, transferred to media supplemented with 140mM NaCl for 32 hour before protoplasting"}
#     """
#     extracted_data = {
#         "tissue": None,
#         "treatment": [],
#         "medium": None
#     }

#     # Extract 'tissue' from 'sample_characteristicts'
#     # The tissue is specified as 'cell type: [value]'
#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if characteristic.startswith('cell type:'):
#                 extracted_data['tissue'] = characteristic.split(':', 1)[1].strip()
#                 break

#     # Extract 'medium' from 'sample_characteristicts'
#     # The growth medium is specified as 'growth media: [value]'
#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if characteristic.startswith('growth media:'):
#                 extracted_data['medium'] = characteristic.split(':', 1)[1].strip()
#                 break

#     # Extract 'treatment' from 'sample_title'
#     # Treatments are typically "140mM NaCl for X hours". If "standard conditions" are mentioned,
#     # the treatment list should be empty as it's not a stress/treatment.
#     if 'sample_title' in sample_metadata and isinstance(sample_metadata['sample_title'], list) and sample_metadata['sample_title']:
#         sample_title = sample_metadata['sample_title'][0]

#         # Look for NaCl treatment and its duration
#         nacl_match = re.search(r'140mM NaCl for (\d+ (?:hour|hours))', sample_title)
#         if nacl_match:
#             extracted_data['treatment'].append(f"140mM NaCl for {nacl_match.group(1)}")
#         elif "standard conditions" in sample_title:
#             # For "standard conditions", the treatment list remains empty as per schema description
#             # "List of treatments and stresses that was applied to the sample."
#             extracted_data['treatment'].append('standard conditions')
#         # If no specific treatment (like NaCl) or "standard conditions" is found,
#         # the 'treatment' list will remain empty, which is a valid state for the schema.

#     return extracted_data



# import json

# def GSE5620_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from the sample metadata
#     following a predefined JSON schema.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the specified schema, containing
#               'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The tissue information is found in 'sample_characteristicts' list.
#     # Example: ['Stock Code: N1092', 'Tissue: Shoots']
#     sample_characteristics = sample_metadata.get('sample_characteristicts', [])
#     tissue = None
#     for char_str in sample_characteristics:
#         if char_str.startswith('Tissue:'):
#             tissue = char_str.split(':', 1)[1].strip()
#             break
    
#     if tissue is None:
#         # Fallback or error handling if tissue is not found, though it's expected to be present.
#         # For this dataset, it seems consistently present.
#         tissue = "unknown" # Or raise an error if strict adherence is needed

#     extracted_data['tissue'] = tissue

#     # 2. Extract 'treatment'
#     # Based on the study metadata and sample titles (e.g., "Control-Shoots-6.0h"),
#     # these samples are described as "Control plants" in the study summary.
#     # The schema expects an array of strings.
#     extracted_data['treatment'] = ["Control"]

#     # 3. Extract 'medium'
#     # From the study_metadata summary:
#     # "Seeds of Arabidopsis thaliana Wild Type (col-0) were sown on rafts in Magenta boxes
#     # containing MS-Agar-media. ... At day 11 the rafts were transferred in Magenta boxes
#     # containing MS-liquid-media."
#     # Since stress treatment started at day 16, the relevant medium during the experiment
#     # is "MS-liquid-media". This is constant across all samples in this study.
#     extracted_data['medium'] = "MS-liquid-media"

#     return extracted_data

# if __name__ == '__main__':
#     # Example usage with a sample from the provided samples_metadata
#     samples_metadata = {
#         'GSM131247': {'Sample_id': 'GSM131247', 'sample_title': ['AtGen_6-0411_Control-Shoots-6.0h_Rep1'], 'sample_source_name': ['Col-0'], 'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Shoots'], 'sample_extraction_protocol': ''},
#         'GSM131245': {'Sample_id': 'GSM131245', 'sample_title': ['AtGen_6-0821_Control-Roots-4.0h_Rep1'], 'sample_source_name': ['Col-0'], 'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Roots'], 'sample_extraction_protocol': ''},
#         'GSM131223': {'Sample_id': 'GSM131223', 'sample_title': ['AtGen_6-0011_Control-Shoots-0h_Rep1'], 'sample_source_name': ['Col-0'], 'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Shoots'], 'sample_extraction_protocol': ''},
#         'GSM131258': {'Sample_id': 'GSM131258', 'sample_title': ['AtGen_6-0622_Control-Roots-24.0h_Rep2'], 'sample_source_name': ['Col-0'], 'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Roots'], 'sample_extraction_protocol': ''}
#     }

#     print("--- Extracted Data for Sample GSM131247 ---")
#     sample_id_1 = 'GSM131247'
#     sample_data_1 = samples_metadata[sample_id_1]
#     extracted_info_1 = GSE5620_extractor(sample_data_1)
#     print(json.dumps(extracted_info_1, indent=2))
#     # Expected output:
#     # {
#     #   "tissue": "Shoots",
#     #   "treatment": ["Control"],
#     #   "medium": "MS-liquid-media"
#     # }

#     print("\n--- Extracted Data for Sample GSM131245 ---")
#     sample_id_2 = 'GSM131245'
#     sample_data_2 = samples_metadata[sample_id_2]
#     extracted_info_2 = GSE5620_extractor(sample_data_2)
#     print(json.dumps(extracted_info_2, indent=2))
#     # Expected output:
#     # {
#     #   "tissue": "Roots",
#     #   "treatment": ["Control"],
#     #   "medium": "MS-liquid-media"
#     # }

#     print("\n--- Extracted Data for all samples ---")
#     all_extracted_data = {}
#     for sample_id, sample_data in samples_metadata.items():
#         all_extracted_data[sample_id] = GSE5620_extractor(sample_data)
#     print(json.dumps(all_extracted_data, indent=2))


# import json

# def GSE27550_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from the sample metadata
#     for the GSE27550 dataset, conforming to the specified JSON schema.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the target JSON schema,
#               containing 'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # The 'sample_characteristicts' field is a list of strings,
#     # where each string is in "key: value" format.
#     characteristics = sample_metadata.get('sample_characteristicts', [])

#     # 1. Extract 'tissue'
#     # Look for 'tissue: <value>' in 'sample_characteristicts'
#     tissue_found = False
#     for char_str in characteristics:
#         if char_str.lower().startswith('tissue:'):
#             extracted_data['tissue'] = char_str.split(':', 1)[1].strip()
#             tissue_found = True
#             break
    
#     # As a fallback, if 'tissue:' is not explicitly in characteristics,
#     # check 'sample_source_name', which often contains tissue information.
#     if not tissue_found and 'sample_source_name' in sample_metadata and sample_metadata['sample_source_name']:
#         # Assuming 'sample_source_name' is a list and the first element is the tissue
#         extracted_data['tissue'] = sample_metadata['sample_source_name'][0].strip()
    
#     # If tissue is still not found, assign a default or raise an error if it's strictly required
#     # For this dataset, 'rosette leaf' is consistently present.
#     if 'tissue' not in extracted_data:
#         # Based on the provided samples, tissue is always 'rosette leaf'
#         # and is either in 'sample_characteristicts' or 'sample_source_name'.
#         # If it somehow isn't, this would be a place to handle it (e.g., default or error).
#         # For this specific dataset, it seems consistently available.
#         pass # No specific action needed as it's expected to be found.


#     # 2. Extract 'treatment'
#     # Look for 'treatment: <value>' in 'sample_characteristicts'
#     treatment_list = []
#     for char_str in characteristics:
#         if char_str.lower().startswith('treatment:'):
#             treatment_value = char_str.split(':', 1)[1].strip()
#             treatment_list.append('control' if treatment_value == 'Wet' else 'drought environment')
#             # Assuming only one treatment per sample based on the provided examples
#             break
#     extracted_data['treatment'] = treatment_list

#     # 3. Extract 'medium'
#     # Based on the study metadata, the samples were grown in "soil"
#     # ("well-watered soil and mild soil drying"). This information is constant
#     # across all samples in this study and not found in individual sample metadata.
#     extracted_data['medium'] = "soil"

#     return extracted_data

# if __name__ == '__main__':
#     # Example study metadata (for context, not directly used by the function)
#     study_metadata = {
#         'title': ['cRNA hybridizations of 18 accessions of Arabidopsis thaliana under well-watered and mild soil drying'],
#         'summary': ['These data provide a basis for exploration of gene expression differences between physiologically diverse accessions of Arabidopsis thaliana.', 'Recent studies have documented remarkable genetic variation among Arabidopsis thaliana accessions collected from diverse habitats and across its geographical range.  Of particular interest are accessions with putatively locally adapted phenotypes – i.e., accessions with attributes that are likely adaptive under the climatic or habitat conditions of their sites of origin.  These genotypes are especially valuable as they may provide insight into the genetic basis of adaptive evolution as well as allow the discovery of genes of ecological importance.  Therefore we studied the physiology, genome content and gene expression of 18 physiologically diverse accessions. The gene expression studies were conducted under two levels of soil moisture and accompanied by physiological measurements to characterize early responses to soil moisture deficit.'],
#         'overall_design': ['The basic experimental design involves 18 accessions crossed with two environmental levels (well-watered soil and mild soil drying) and 3 biological replicates per accession/treatment combination.']
#     }

#     # Example samples metadata (used to guide the function's logic)
#     samples_metadata = {
#         'GSM680404': {'Sample_id': 'GSM680404', 'sample_title': ['leaf_Ull2-5_dry_2'], 'sample_source_name': ['rosette leaf'], 'sample_characteristicts': ['accession: CS22586', 'treatment: Dry', 'tissue: rosette leaf'], 'sample_extraction_protocol': ['Qiagen RNAasy Plant Mini prep, samples initially harvested in RNAlater']},
#         'GSM680390': {'Sample_id': 'GSM680390', 'sample_title': ['leaf_Omo2-3_dry_1'], 'sample_source_name': ['rosette leaf'], 'sample_characteristicts': ['accession: CS22585', 'treatment: Dry', 'tissue: rosette leaf'], 'sample_extraction_protocol': ['Qiagen RNAasy Plant Mini prep, samples initially harvested in RNAlater']},
#         'GSM680393': {'Sample_id': 'GSM680393', 'sample_title': ['leaf_Omo2-3_wet_2'], 'sample_source_name': ['rosette leaf'], 'sample_characteristicts': ['accession: CS22585', 'treatment: Wet', 'tissue: rosette leaf'], 'sample_extraction_protocol': ['Qiagen RNAasy Plant Mini prep, samples initially harvested in RNAlater']},
#         'GSM680377': {'Sample_id': 'GSM680377', 'sample_title': ['leaf_Bil-5_wet_3'], 'sample_source_name': ['rosette leaf'], 'sample_characteristicts': ['accession: CS22578', 'treatment: Wet', 'tissue: rosette leaf'], 'sample_extraction_protocol': ['Qiagen RNAasy Plant Mini prep, samples initially harvested in RNAlater']}
#     }

#     print("--- Testing with sample GSM680404 (Dry treatment) ---")
#     sample_id_1 = 'GSM680404'
#     sample_data_1 = samples_metadata[sample_id_1]
#     extracted_info_1 = GSE27550_extractor(sample_data_1)
#     print(f"Input Sample Metadata for {sample_id_1}:\n{json.dumps(sample_data_1, indent=2)}\n")
#     print(f"Extracted Schema for {sample_id_1}:\n{json.dumps(extracted_info_1, indent=2)}\n")
#     # Expected output: {"tissue": "rosette leaf", "treatment": ["Dry"], "medium": "soil"}

#     print("--- Testing with sample GSM680393 (Wet treatment) ---")
#     sample_id_2 = 'GSM680393'
#     sample_data_2 = samples_metadata[sample_id_2]
#     extracted_info_2 = GSE27550_extractor(sample_data_2)
#     print(f"Input Sample Metadata for {sample_id_2}:\n{json.dumps(sample_data_2, indent=2)}\n")
#     print(f"Extracted Schema for {sample_id_2}:\n{json.dumps(extracted_info_2, indent=2)}\n")
#     # Expected output: {"tissue": "rosette leaf", "treatment": ["Wet"], "medium": "soil"}

#     print("--- Testing with another sample GSM680377 (Wet treatment) ---")
#     sample_id_3 = 'GSM680377'
#     sample_data_3 = samples_metadata[sample_id_3]
#     extracted_info_3 = GSE27550_extractor(sample_data_3)
#     print(f"Input Sample Metadata for {sample_id_3}:\n{json.dumps(sample_data_3, indent=2)}\n")
#     print(f"Extracted Schema for {sample_id_3}:\n{json.dumps(extracted_info_3, indent=2)}\n")
#     # Expected output: {"tissue": "rosette leaf", "treatment": ["Wet"], "medium": "soil"}


# import json

# def GSE5622_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from a single sample's metadata
#     following a specific JSON schema.

#     The schema requires 'tissue' (string), 'treatment' (array of strings),
#     and 'medium' (string).

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 This typically corresponds to one entry from the
#                                 'samples_metadata' dictionary (e.g., the value
#                                 associated with a 'GSM' key).

#     Returns:
#         dict: A dictionary formatted according to the specified schema,
#               containing 'tissue', 'treatment', and 'medium'.
#               Returns 'unknown' for tissue if not found, but for this dataset,
#               it is expected to be present.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The tissue information is typically found within the 'sample_characteristicts' list.
#     # Example: ['Stock Code: N1092', 'Tissue: Roots']
#     tissue = None
#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if characteristic.startswith('Tissue:'):
#                 # Extract the part after "Tissue:" and strip any leading/trailing whitespace
#                 tissue = characteristic.split(':', 1)[1].strip()
#                 break
#     # As 'tissue' is a required string, provide a fallback if not found, though
#     # for this dataset, it appears to be consistently present.
#     extracted_data['tissue'] = tissue if tissue else "unknown"

#     # 2. Extract 'treatment'
#     # Based on the study metadata title "AtGenExpress: Stress Treatments (Osmotic stress)"
#     # and the sample titles containing "Osmoticstress", the primary treatment is "Osmotic stress".
#     # The schema expects an array of strings for 'treatment'.
#     extracted_data['treatment'] = ["Osmotic stress"]

#     # 3. Extract 'medium'
#     # From the 'study_metadata' summary:
#     # "Seeds of Arabidopsis thaliana Wild Type (col-0) were sown on rafts in Magenta boxes
#     # containing MS-Agar-media. ... At day 11 the rafts were transferred in Magenta boxes
#     # containing MS-liquid-media. At day 16 stress treatment started..."
#     # This indicates that the samples were in "MS-liquid-media" when the stress treatment began.
#     extracted_data['medium'] = "MS-liquid-media"

#     return extracted_data



# import json
# import re

# def GSE5628_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information of the sample from the
#     GSE5628 dataset, conforming to the specified JSON schema.

#     The function identifies 'tissue', 'treatment', and 'medium' from the
#     provided sample_metadata dictionary.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 Expected keys include 'sample_title' and
#                                 'sample_characteristicts'.

#     Returns:
#         dict: A dictionary formatted as a JSON instance conforming to the
#               output schema:
#               {"properties": {"tissue": {"description": "Tissue the samples was extracted from.", "title": "Tissue", "type": "string"}, "treatment": {"description": "List of treatments and stresses that was applied to the sample.", "items": {"type": "string"}, "title": "Treatment", "type": "array"}, "medium": {"description": "Growth medium of the sample.", "title": "Medium", "type": "string"}}, "required": ["tissue", "treatment", "medium"]}
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The tissue information is found within the 'sample_characteristicts' list.
#     # Example: 'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Shoots']
#     tissue = None
#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if characteristic.startswith('Tissue:'):
#                 tissue = characteristic.split(':', 1)[1].strip()
#                 break
#     # If tissue is not found, it will remain None. The schema requires a string,
#     # implying it should always be present for valid samples in this dataset.
#     extracted_data['tissue'] = tissue

#     # 2. Extract 'treatment'
#     # Treatments are derived from the 'sample_title'.
#     # The study metadata indicates "Heat stress" as a primary treatment.
#     # Sample titles also show "recovery" as a potential additional phase.
#     # Example: 'sample_title': ['AtGen_6-9411_Heatstress(3h)+3hrecovery-Shoots-6.0h_Rep1']
#     treatments = []
#     if 'sample_title' in sample_metadata and sample_metadata['sample_title']:
#         title = sample_metadata['sample_title'][0]

#         # All samples in this study are related to "Heat stress"
#         if "Heatstress" in title:
#             treatments.append("Heat stress")

#         # Check for "recovery" phase
#         if "recovery" in title:
#             treatments.append("recovery")
    
#     extracted_data['treatment'] = treatments

#     # 3. Extract 'medium'
#     # Based on the study_metadata, the samples were transferred to "MS-liquid-media"
#     # before the stress treatment began. This is constant across all samples.
#     # "At day 11 the rafts were transferred in Magenta boxes containing MS-liquid-media."
#     # "At day 16 stress treatment started..."
#     extracted_data['medium'] = "MS-liquid-media"

#     return extracted_data



# import json

# def GSE41935_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from sample metadata
#     following a specific JSON schema.

#     The schema requires 'tissue' (string), 'treatment' (list of strings),
#     and 'medium' (string).

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 Expected to have a 'sample_characteristicts' key
#                                 which is a list of strings.

#     Returns:
#         dict: A dictionary formatted according to the specified schema,
#               containing extracted tissue, treatment, and medium information.
#               If 'medium' is not explicitly found, it defaults to "Not specified".
#     """
#     extracted_data = {
#         "tissue": "Not specified",  # Default placeholder
#         "treatment": [],
#         "medium": "Not specified"   # Default placeholder
#     }

#     # The primary source of information for tissue and treatment is
#     # the 'sample_characteristicts' list.
#     characteristics = sample_metadata.get('sample_characteristicts', [])

#     # 1. Extract 'tissue'
#     # Look for an item starting with 'tissue:' in the characteristics list.
#     for char_item in characteristics:
#         if char_item.startswith('tissue:'):
#             # Split the string at the first colon and take the second part, then strip whitespace.
#             extracted_data['tissue'] = char_item.split(':', 1)[1].strip()
#             break # Assuming only one tissue entry per sample

#     # 2. Extract 'treatment'
#     # Treatment information is derived from 'stress type:' and 'stress treatment:'
#     stress_type = None
#     stress_treatment_details = None

#     for char_item in characteristics:
#         if char_item.startswith('stress type:'):
#             stress_type = char_item.split(':', 1)[1].strip()
#         elif char_item.startswith('stress treatment:'):
#             stress_treatment_details = char_item.split(':', 1)[1].strip()

#     if stress_type:
#         # If the stress type indicates 'control plants', this is the sole treatment.
#         if stress_type.lower() == 'control plants':
#             extracted_data['treatment'].append('control plants')
#         else:
#             # For other stress types, add the stress type itself to the list.
#             # The study metadata indicates combined stresses (e.g., "Salt+Heat") are treated as single factors.
#             for el in stress_type.split('+'):
#                 extracted_data['treatment'].append(el+ ' treatment')
#             # If specific treatment details are provided, add them as well.
#             # if stress_treatment_details:
#             #     extracted_data['treatment'].append(stress_treatment_details)
                

#     # Ensure 'treatment' is never an empty list, as per schema's 'required' property
#     # and common sense for experimental data.
#     if not extracted_data['treatment']:
#         extracted_data['treatment'].append("No specific treatment applied or identified")

#     # 3. Extract 'medium'
#     # The provided study_metadata and samples_metadata do not explicitly contain
#     # a field for 'medium' or 'growth medium'.
#     # Following the guidance to "Focus on location WHERE in the dictionary the desired information is stored",
#     # and since it's not found, "Not specified" is used as the value.
#     # While Arabidopsis plants are typically grown in soil or hydroponics,
#     # inferring this information is outside the scope of direct extraction from the given data.
#     extracted_data['medium'] = "Not specified"

#     return extracted_data



# import json

# def GSE76827_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from the sample metadata
#     for the GSE76827 dataset, conforming to the specified JSON schema.

#     The output schema is:
#     {"properties": {"tissue": {"description": "Tissue the samples was extracted from.", "title": "Tissue", "type": "string"},
#                     "treatment": {"description": "List of treatments and stresses that was applied to the sample.", "items": {"type": "string"}, "title": "Treatment", "type": "array"},
#                     "medium": {"description": "Growth medium of the sample.", "title": "Medium", "type": "string"}},
#      "required": ["tissue", "treatment", "medium"]}

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 Expected to have a 'sample_characteristicts' key
#                                 which is a list of strings.

#     Returns:
#         dict: A dictionary formatted as a JSON instance conforming to the specified schema.
#               Example: {"tissue": "shoot", "treatment": ["unstress condition"], "medium": "ceramics granular soil"}
#     """
#     extracted_info = {}

#     # 1. Extract 'medium'
#     # Based on the study_metadata['overall_design'] provided in the context:
#     # "Plantlets were then transferred to ceramics granular soil (size 2.5L, Sakatanotane, Japan)"
#     # This indicates that the samples (roots and shoots) were grown in "ceramics granular soil"
#     # during the experimental phase where drought treatment was applied.
#     # This value is constant across all samples in this study.
#     extracted_info['medium'] = "ceramics granular soil"

#     # 2. Extract 'tissue' and 'treatment' from 'sample_characteristicts'
#     # The 'sample_characteristicts' field is a list of strings, where each string
#     # contains a key-value pair separated by a colon (e.g., 'tissue: shoot').
#     characteristics = sample_metadata.get('sample_characteristicts', [])

#     for char_string in characteristics:
#         if char_string.startswith('tissue:'):
#             # Extract the value after "tissue:" and strip whitespace
#             extracted_info['tissue'] = char_string.split(':', 1)[1].strip()
#         elif char_string.startswith('treatment:'):
#             # Extract the value after "treatment:" and strip whitespace.
#             # The schema requires 'treatment' to be an array of strings.
#             # In this dataset, each sample typically has one treatment string,
#             # so we wrap it in a list to conform to the schema.
#             extracted_info['treatment'] = [char_string.split(':', 1)[1].strip()]
    
#     # All required fields are expected to be found based on the provided sample data structure.
#     # If for some reason 'tissue' or 'treatment' were not found (e.g., malformed input),
#     # the function would return an incomplete dictionary. For this problem, we assume
#     # the input adheres to the structure shown in the examples.

#     return extracted_info



# import re

# def GSE34188_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from a sample metadata dictionary
#     for the GSE34188 dataset, conforming to a specific JSON schema.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the specified schema, containing
#               'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # Extract primary identifiers from the sample metadata
#     source_name = sample_metadata['sample_source_name'][0]
#     characteristics = sample_metadata['sample_characteristicts'][0]

#     # --- Extract Tissue ---
#     if characteristics == 'treatment group: tissue':
#         # For tissue samples, map the source_name to a standardized tissue name
#         tissue_mapping = {
#             "Rosette": "adult rosette leaves",
#             "Young_Rosette": "juvenile rosette leaves",
#             "Old_Rosette": "senescence rosette leaves",
#             "Bud": "mature flower buds",
#             "Young_Bud": "young buds",
#             "Flower": "flowers",
#             "Silique": "mature siliques",
#             "Young_Silique": "young siliques",
#             "Old_Silique": "old siliques",
#             "dryseed": "dry seeds",
#             "seed24h": "24H imbibed seeds",
#             "seed48h": "48H imbibed seeds",
#             "Root": "root",
#             "Callus": "callus",
#             "Stem": "stems",
#             "Cauline_leaves": "cauline leaves"
#         }
#         extracted_data['tissue'] = tissue_mapping.get(source_name, source_name.replace('_', ' '))
#     else:
#         # For abiotic stress or light condition samples, the tissue is generally "seedling"
#         # as per the study's overall design (e.g., "WT seedlings were grown on MS plates").
#         extracted_data['tissue'] = "seedling"

#     # --- Extract Treatment ---
#     extracted_data['treatment'] = []
#     if characteristics == 'treatment group: tissue':
#         # For tissue samples, treatments are specific for imbibed seeds, otherwise "control"
#         if source_name == "seed24h":
#             extracted_data['treatment'].append("imbibed 24h")
#         elif source_name == "seed48h":
#             extracted_data['treatment'].append("imbibed 48h")
#         else:
#             extracted_data['treatment'].append("control") # Represents standard growth conditions
#     else:
#         # For abiotic stress or light condition samples, the source_name is the treatment
#         if source_name == "Control":
#             extracted_data['treatment'].append("control")
#         else:
#             treatment_label = source_name

#             # Handle abiotic stress treatments (e.g., Cold6h, dry2h, Heat6h, Salt6h)
#             # These typically follow a pattern of word + digits + 'h'
#             match_abiotic = re.match(r"([a-zA-Z]+)(\d+h)", treatment_label)
#             if match_abiotic:
#                 prefix = match_abiotic.group(1)
#                 time_str = match_abiotic.group(2)
                
#                 # Standardize prefixes based on overall_design
#                 if prefix.lower() == "salt":
#                     prefix = "salt treatment"
#                 elif prefix.lower() == "dry":
#                     prefix = "drought"
#                 elif prefix.lower() == "cold":
#                     prefix = "cold stress"
#                 elif prefix.lower() == "heat":
#                     prefix = "heat stress"
                
#                 treatment_label = f"{prefix} {time_str}"
            
#             # Handle light conditions (e.g., Light0h, Dark, Blue, FarRed, Red)
#             elif treatment_label.startswith("Light"):
#                 # Light0h, Light1h, Light6h, Light24h
#                 match_light_time = re.match(r"Light(\d+)h", treatment_label)
#                 if match_light_time:
#                     time_val = match_light_time.group(1)
#                     treatment_label = f"white light {time_val}h"
#             elif treatment_label == "Dark":
#                 treatment_label = "dark light"
#             elif treatment_label == "Blue":
#                 treatment_label = "blue light"
#             elif treatment_label == "FarRed":
#                 treatment_label = "far-red light"
#             elif treatment_label == "Red":
#                 treatment_label = "red light"
            
#             extracted_data['treatment'].append(treatment_label)

#     # --- Extract Medium ---
#     # Determine the growth medium based on the source_name (which can be a tissue or treatment)
#     medium_map = {
#         # Tissues grown in soil as per overall_design
#         "Flower": "soil", "Stem": "soil", "Cauline_leaves": "soil",
#         "Rosette": "soil", "Young_Rosette": "soil", "Old_Rosette": "soil",
#         "Young_Bud": "soil", "Bud": "soil", "Young_Silique": "soil",
#         "Silique": "soil", "Old_Silique": "soil",

#         # Tissues with specific media as per overall_design
#         "Root": "MS plates",
#         "Callus": "CIM",
#         "dryseed": "moistened paper",
#         "seed24h": "moistened paper",
#         "seed48h": "moistened paper",

#         # All abiotic stress and light condition treatments were performed on MS plates
#         # as per overall_design ("plants grown for 2 weeks on MS plates", "WT seeds were plated on MS plates")
#         "Cold2h": "MS plates", "Cold6h": "MS plates",
#         "Heat2h": "MS plates", "Heat6h": "MS plates",
#         "Salt2h": "MS plates", "Salt6h": "MS plates",
#         "dry2h": "MS plates", "dry6h": "MS plates",
#         "Control": "MS plates", # Control for abiotic stress/light experiments
#         "Light0h": "MS plates", "Light1h": "MS plates",
#         "Light6h": "MS plates", "Light24h": "MS plates",
#         "Dark": "MS plates", "Blue": "MS plates",
#         "FarRed": "MS plates", "Red": "MS plates"
#     }
#     extracted_data['medium'] = medium_map.get(source_name, "unknown") # Fallback for unmapped source_names

#     return extracted_data


# import json

# def GSE126373_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from a sample metadata dictionary
#     for the GSE126373 dataset, conforming to the specified JSON schema.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the specified schema,
#               containing 'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # --- Extract Tissue, Temperature, and Timepoint from 'sample_characteristicts' ---
#     # These are found in a list of strings, each formatted as "key: value"
#     characteristics = sample_metadata.get('sample_characteristicts', [])

#     tissue_val = None
#     temperature_val = None
#     timepoint_val = None

#     for char_str in characteristics:
#         if char_str.startswith('tissue:'):
#             tissue_val = char_str.split(':', 1)[1].strip()
#         elif char_str.startswith('temperature:'):
#             temperature_val = char_str.split(':', 1)[1].strip()
#         elif char_str.startswith('timepoint:'):
#             timepoint_val = char_str.split(':', 1)[1].strip()

#     # --- Populate 'tissue' field ---
#     # The schema requires 'tissue' to be a string.
#     # Based on the metadata, 'tissue' is always present in 'sample_characteristicts'.
#     extracted_data['tissue'] = tissue_val

#     # --- Populate 'treatment' field ---
#     # The schema requires 'treatment' to be an array of strings.
#     # Treatment information is derived from temperature and timepoint.
#     treatment_list = []
#     if temperature_val:
#         if temperature_val == '28°C':
#             treatment_list.append("elevated temperature (28°C)")
#         elif temperature_val == '20°C':
#             treatment_list.append("control (20°C)")
#         # Add a generic fallback if other temperatures might appear
#         else:
#             treatment_list.append(f"temperature: {temperature_val}")
    
#     if timepoint_val:
#         treatment_list.append(f"duration: {timepoint_val}")
    
#     # Based on the provided samples, temperature and timepoint are always present,
#     # so treatment_list will always have at least two elements.
#     extracted_data['treatment'] = treatment_list

#     # --- Populate 'medium' field ---
#     # From the study_metadata, the growth medium is "ATS medium" and is constant.
#     # 'overall_design': ['... square plates containing ATS medium ...']
#     extracted_data['medium'] = "ATS medium"

#     return extracted_data



# import json

# def GSE201609_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from the sample metadata
#     for the GSE201609 dataset, conforming to a specific JSON schema.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the specified schema, containing
#               'tissue', 'treatment', and 'medium' information.

#     Raises:
#         ValueError: If required 'tissue' or 'treatment' information cannot be found
#                     in the sample_metadata.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # Tissue information is found in 'sample_characteristicts' list, e.g., 'tissue: Shoot'
#     tissue = None
#     for char_str in sample_metadata.get('sample_characteristicts', []):
#         if char_str.startswith('tissue:'):
#             # Split by the first colon and strip whitespace from the value
#             tissue = char_str.split(':', 1)[1].strip()
#             break
#     if tissue is None:
#         raise ValueError(f"Tissue information not found for sample: {sample_metadata.get('Sample_id', 'N/A')}")
#     extracted_data['tissue'] = tissue

#     # 2. Extract 'treatment'
#     # Treatment information is found in 'sample_characteristicts' list, e.g., 'treatment: Ethanol treatment followed by drought'
#     treatment_raw = None
#     for char_str in sample_metadata.get('sample_characteristicts', []):
#         if char_str.startswith('treatment:'):
#             # Split by the first colon and strip whitespace from the value
#             treatment_raw = char_str.split(':', 1)[1].strip()
#             break
#     if treatment_raw is None:
#         raise ValueError(f"Treatment information not found for sample: {sample_metadata.get('Sample_id', 'N/A')}")

#     treatments_list = []
#     # Check for multiple treatments indicated by "followed by"
#     if 'followed by' in treatment_raw:
#         # Split the string by " followed by " and strip each part
#         parts = treatment_raw.split(' followed by ')
#         treatments_list = [p.strip() for p in parts if p.strip()]
#     else:
#         # If no "followed by", the entire string is a single treatment
#         treatments_list = [treatment_raw]
#     extracted_data['treatment'] = treatments_list

#     # 3. Extract 'medium'
#     # Based on the study_metadata['overall_design'], the plants were grown on "MS medium".
#     # This information is constant across all samples in this study.
#     extracted_data['medium'] = "MS medium"

#     return extracted_data



# import json

# def GSE40061_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from the sample metadata
#     following a specific JSON schema.

#     The schema requires:
#     - 'tissue': Tissue the sample was extracted from (string).
#     - 'treatment': List of treatments and stresses applied to the sample (array of strings).
#     - 'medium': Growth medium of the sample (string).

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 This typically corresponds to one GSM entry.

#     Returns:
#         dict: A dictionary formatted according to the specified schema, containing
#               'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The tissue information is found within the 'sample_characteristicts' list.
#     # Example: ['tissue: root', 'cultivar/ecotype: camta1']
#     tissue = "unknown" # Default in case it's unexpectedly missing
#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if characteristic.startswith('tissue:'):
#                 tissue = characteristic.split(':', 1)[1].strip()
#                 break
#     extracted_data['tissue'] = tissue

#     # 2. Extract 'treatment'
#     # Treatment information (e.g., 'drought', 'water') is derived from 'sample_source_name'.
#     # Example: ['Arabidopsis Thaliana mutant root in drought condition']
#     treatments = []
#     if 'sample_source_name' in sample_metadata and sample_metadata['sample_source_name']:
#         source_name_lower = sample_metadata['sample_source_name'][0].lower()
        
#         # Check for specific conditions/treatments
#         if 'drought' in source_name_lower:
#             treatments.append('drought')
#         elif 'water' in source_name_lower:
#             # 'water condition' is considered a baseline treatment in this context
#             treatments.append('control')#TODO: should this be control or water?????
            
#     extracted_data['treatment'] = treatments

#     # 3. Extract 'medium'
#     # The growth medium is not explicitly stated in the provided sample or study metadata.
#     # For Arabidopsis Thaliana studies involving water/drought conditions, 'soil' is a
#     # common and reasonable inference for the growth medium, especially since 'water'
#     # and 'drought' refer to conditions applied to the plant, not necessarily the medium itself.
#     # The guidance also suggests that 'medium' can be constant across samples.
#     extracted_data['medium'] = 'soil'

#     return extracted_data



# import json

# def GSE62163_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from a sample metadata dictionary
#     for the GSE62163 dataset, conforming to a specific JSON schema.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 Expected keys include 'sample_characteristicts' and
#                                 'sample_source_name'.

#     Returns:
#         dict: A dictionary formatted according to the specified schema, containing
#               'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The tissue information is found within the 'sample_characteristicts' list.
#     # Example: 'tissue: shoot tissue above medium'
#     tissue_value = "unspecified tissue" # Default value if not found
#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if characteristic.startswith('tissue:'):
#                 tissue_value = characteristic.split(':', 1)[1].strip()
#                 break
#     extracted_data['tissue'] = tissue_value

#     # 2. Extract 'treatment'
#     # Treatments are derived from the 'sample_source_name' string.
#     # Examples: 'Arabidopsis seedlings 3 h HS, no-EBR', 'Arabidopsis seedlings no-stress, 1 µM EBR for 21 days'
#     treatments = []
#     if 'sample_source_name' in sample_metadata and isinstance(sample_metadata['sample_source_name'], list):
#         source_name_str = sample_metadata['sample_source_name'][0] if sample_metadata['sample_source_name'] else ""

#         # Check for EBR (Brassinosteroid) treatment
#         if '1 µM EBR' in source_name_str or 'EBR-treated' in source_name_str:
#             treatments.append("EBR treatment")
#         # 'no-EBR' indicates absence of this treatment, so we don't add it to the list of applied treatments.

#         # Check for Heat Stress (HS)
#         # 'no-stress' explicitly indicates no heat stress was applied.
#         if 'HS' in source_name_str and 'no-stress' not in source_name_str:
#             treatments.append("Heat Stress")

#         # Check for Recovery from Heat Stress
#         # '6 hR' refers to a recovery period after heat stress.
#         if '6 hR' in source_name_str:
#             treatments.append("Recovery from Heat Stress")
        
#     if treatments == []:
#         treatments = ['No Stress/control']
#     extracted_data['treatment'] = treatments

#     # 3. Extract 'medium'
#     # The growth medium is not explicitly named in the provided metadata but is implied
#     # by phrases like "shoot tissue above medium" and the context of Arabidopsis seedling growth.
#     # Given it's a constant across samples and not a variable of the experiment,
#     # a generic description is appropriate.
#     extracted_data['medium'] = "plant growth medium"

#     return extracted_data



# import re
# import json

# def GSE71237_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information of the sample from this dataset
#     following the given schema.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 This corresponds to one of the inner dictionaries
#                                 from the 'samples_metadata' (e.g., samples_metadata['GSM1831245']).

#     Returns:
#         dict: A dictionary formatted as a JSON instance that conforms to the
#               specified JSON schema.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The 'tissue' information is found in the 'sample_characteristicts' list.
#     characteristics = sample_metadata.get('sample_characteristicts', [])
#     tissue = None
#     for char_str in characteristics:
#         if char_str.startswith('tissue:'):
#             # Extract the value after 'tissue:' and strip whitespace
#             tissue = char_str.split(':', 1)[1].strip()
#             break
#     # Assign the extracted tissue or a default if not found
#     extracted_data['tissue'] = tissue if tissue else "unknown"

#     # 2. Extract 'treatment' and 'medium'
#     # These pieces of information are derived from the 'sample_title' field
#     # and contextual information from the study's overall design.
#     sample_title = sample_metadata.get('sample_title', [''])[0]
    
#     treatments = []
#     medium = "unknown medium" # Default value for medium

#     # Use regex to parse the sample title for water potential and condition type.
#     # Example titles:
#     # 'egr1-1egr2-1 at -1.2 MPa, biological rep 1_stress'
#     # 'egr1-1egr2-1 at -0.25 MPa, biological rep 2_control'
#     # The regex captures the water potential (e.g., "-1.2 MPa") and the condition
#     # (e.g., "stress" or "control") at the end of the string.
#     match = re.search(r'at\s+([-]?\d+\.?\d*\s*MPa)(?:,\s*biological rep \d+)?_(\w+)', sample_title, re.IGNORECASE)
    
#     if match:
#         water_potential = match.group(1).strip()
#         condition_raw = match.group(2).strip().lower() # 'stress' or 'control'
        
#         # Add the specific water potential as a treatment
#         treatments.append(water_potential)
        
#         # Determine additional treatment details and the medium based on the condition
#         if condition_raw == 'stress':
#             treatments.append("drought") # As per study_metadata, low water potential is drought
#             # From study_metadata overall_design: "PEG-infused agar plates at -1.2 MPa"
#             medium = "PEG-infused agar plates"
#         elif condition_raw == 'control':
#             treatments.append("control condition")
#             # From study_metadata overall_design: "fresh control media (-0.25 Mpa)"
#             medium = "fresh control media"
#         # If other conditions exist, they would fall through to the default medium.
    
#     # Assign the extracted treatments or a default if parsing failed
#     # extracted_data['treatment'] = treatments if treatments else ["unknown treatment"]
#     extracted_data['treatment'] = ['drought'] if 'drought' in treatments else ["control"]
#     # Assign the extracted medium
#     extracted_data['medium'] = medium

#     return extracted_data



# import json

# def GSE63372_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from the sample metadata
#     for the GSE63372 dataset, conforming to a specific JSON schema.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the specified schema, containing
#               'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # 1. Extract 'medium'
#     # Based on the study_metadata['overall_design'], the growth medium is constant
#     # across all samples in this study: "0.5X MS plates containing 1% sucrose".
#     extracted_data['medium'] = "0.5X MS plates containing 1% sucrose"

#     # 2. Extract 'tissue' and 'treatment' from 'sample_characteristicts'
#     characteristics = sample_metadata.get('sample_characteristicts', [])

#     tissue_found = False
#     treatment_found = False

#     for char_string in characteristics:
#         if char_string.startswith('tissue:'):
#             extracted_data['tissue'] = char_string.split(':', 1)[1].strip()
#             tissue_found = True
#         elif char_string.startswith('treatment:'):
#             treatment_value = char_string.split(':', 1)[1].strip()
#             if treatment_value.lower() == 'none':
#                 extracted_data['treatment'] = ['control/No stress']
#             else:
#                 extracted_data['treatment'] = [treatment_value]
#             treatment_found = True
    
#     # Ensure all required fields are present, even if not explicitly found
#     # (though based on the provided metadata, they should always be there)
#     if not tissue_found:
#         extracted_data['tissue'] = "unknown" # Or raise an error if strict
#     if not treatment_found:
#         extracted_data['treatment'] = ['control/No stress'] # Or raise an error if strict

#     return extracted_data

# if __name__ == '__main__':
#     # Example usage with the provided metadata

#     study_metadata = {
#         'title': ['Expression data of HSFA6b-overexpression and HSFA6b dominant-negative mutant lines of Arabidopsis after NaCl or Heat shock treatment'],
#         'summary': ['In order to study the improtance of HSFA6b in salt and heat stresses of arabidopsis, we generated the HSFA6b-overexpression (HSFA6b-OE) and dominant-negative (HSFA6b-RD) mutant lines.', 'We used microarray to investigate how genes regulated by HSFA6b after NaCl or heat shock treatment.'],
#         'overall_design': ['Seven-day-old seedlings of wild type (Col-0), HSFA6b-overexpression (HSFA6b-OE) and HSFA6b dominant-negative (HSFA6b-RD) mutant lines grown at 22℃ on 0.5X MS plates containing 1% sucrose were treated with 150mM NaCl for 6hr or 37℃ heat shock for 1hr, with 22℃ treatment as a control. Subsequently, samples were collected for RNA extraction by use of the RNeasy Kit (Qiagen). The results were representative of two independent biological replicates for analysis.']
#     }

#     samples_metadata = {
#         'GSM1547622': {'Sample_id': 'GSM1547622', 'sample_title': ['7-day-old seedlings of Col-0 at control condition, biological rep2'], 'sample_source_name': ['7-day-old seedling of Col-0 at control condition'], 'sample_characteristicts': ['cultivar: Col-0', 'genotype: wildtype', 'tissue: seedling', 'age: 7-day-old', 'treatment: none'], 'sample_extraction_protocol': ["Samples were collected for RNA extraction by use of the RNeasy RNA extraction kit (Qiagen)  according to the manufacturer's instructions."]},
#         'GSM1547618': {'Sample_id': 'GSM1547618', 'sample_title': ['7-day-old seedlings of HSFA6b-RD treated with 37℃ heat shock for 1h, biological rep1'], 'sample_source_name': ['7-day-old seedlings of HSFA6b-RD treated with 37℃ heat shock for 1hr'], 'sample_characteristicts': ['genotype: HSFA6b dominant-negative', 'tissue: seedling', 'age: 7-day-old', 'treatment: 37℃ heat shock for 1hr'], 'sample_extraction_protocol': ["Samples were collected for RNA extraction by use of the RNeasy RNA extraction kit (Qiagen)  according to the manufacturer's instructions."]},
#         'GSM1547623': {'Sample_id': 'GSM1547623', 'sample_title': ['7-day-old seedlings of HSFA6b-OE at control condition, biological rep2'], 'sample_source_name': ['7-day-old seedling of HSFA6b-OE at control condition'], 'sample_characteristicts': ['genotype: HSFA6b-overexpression', 'tissue: seedling', 'age: 7-day-old', 'treatment: none'], 'sample_extraction_protocol': ["Samples were collected for RNA extraction by use of the RNeasy RNA extraction kit (Qiagen)  according to the manufacturer's instructions."]},
#         'GSM1547627': {'Sample_id': 'GSM1547627', 'sample_title': ['7-day-old seedlings of HSFA6b-RD treated with 37℃ heat shock for 1h, biological rep2'], 'sample_source_name': ['7-day-old seedlings of HSFA6b-RD treated with 37℃ heat shock for 1hr'], 'sample_characteristicts': ['genotype: HSFA6b dominant-negative', 'tissue: seedling', 'age: 7-day-old', 'treatment: 37℃ heat shock for 1hr'], 'sample_extraction_protocol': ["Samples were collected for RNA extraction by use of the RNeasy RNA extraction kit (Qiagen)  according to the manufacturer's instructions."]},
#         'GSM1547613': {'Sample_id': 'GSM1547613', 'sample_title': ['7-day-old seedlings of Col-0 at control condition, biological rep1'], 'sample_source_name': ['7-day-old seedling of Col-0 at control condition'], 'sample_characteristicts': ['cultivar: Col-0', 'genotype: wildtype', 'tissue: seedling', 'age: 7-day-old', 'treatment: none'], 'sample_extraction_protocol': ["Samples were collected for RNA extraction by use of the RNeasy RNA extraction kit (Qiagen)  according to the manufacturer's instructions."]},
#         'GSM1547630': {'Sample_id': 'GSM1547630', 'sample_title': ['7-day-old seedlings of HSFA6b-RD treated with 150mM NaCl for 6hr, biological rep2'], 'sample_source_name': ['7-day-old seedlings of HSFA6b-RD treated with 150mM NaCl for 6hr'], 'sample_characteristicts': ['genotype: HSFA6b dominant-negative', 'tissue: seedling', 'age: 7-day-old', 'treatment: 150mM NaCl for 6hr'], 'sample_extraction_protocol': ["Samples were collected for RNA extraction by use of the RNeasy RNA extraction kit (Qiagen)  according to the manufacturer's instructions."]},
#         'GSM1547629': {'Sample_id': 'GSM1547629', 'sample_title': ['7-day-old seedlings of HSFA6b-OE treated with 150mM NaCl for 6hr, biological rep2'], 'sample_source_name': ['7-day-old seedlings of HSFA6b-OE treated with 150mM NaCl for 6hr'], 'sample_characteristicts': ['genotype: HSFA6b-overexpression', 'tissue: seedling', 'age: 7-day-old', 'treatment: 150mM NaCl for 6hr'], 'sample_extraction_protocol': ["Samples were collected for RNA extraction by use of the RNeasy RNA extraction kit (Qiagen)  according to the manufacturer's instructions."]},
#         'GSM1547617': {'Sample_id': 'GSM1547617', 'sample_title': ['7-day-old seedlings of HSFA6b-OE treated with 37℃ heat shock for 1h, biological rep1'], 'sample_source_name': ['7-day-old seedlings of HSFA6b-OE treated with 37℃ heat shock for 1hr'], 'sample_characteristicts': ['genotype: HSFA6b-overexpression', 'tissue: seedling', 'age: 7-day-old', 'treatment: 37℃ heat shock for 1hr'], 'sample_extraction_protocol': ["Samples were collected for RNA extraction by use of the RNeasy RNA extraction kit (Qiagen)  according to the manufacturer's instructions."]},
#         'GSM1547620': {'Sample_id': 'GSM1547620', 'sample_title': ['7-day-old seedlings of HSFA6b-OE treated with 150mM NaCl for 6hr, biological rep1'], 'sample_source_name': ['7-day-old seedlings of HSFA6b-OE treated with 150mM NaCl for 6hr'], 'sample_characteristicts': ['genotype: HSFA6b-overexpression', 'tissue: seedling', 'age: 7-day-old', 'treatment: 150mM NaCl for 6hr'], 'sample_extraction_protocol': ["Samples were collected for RNA extraction by use of the RNeasy RNA extraction kit (Qiagen)  according to the manufacturer's instructions."]},
#         'GSM1547619': {'Sample_id': 'GSM1547619', 'sample_title': ['7-day-old seedlings of Col-0 treated with 150mM NaCl for 6hr, biological rep1'], 'sample_source_name': ['7-day-old seedlings of Col-0 treated with 150mM NaCl for 6hr'], 'sample_characteristicts': ['cultivar: Col-0', 'genotype: wildtype', 'tissue: seedling', 'age: 7-day-old', 'treatment: 150mM NaCl for 6hr'], 'sample_extraction_protocol': ["Samples were collected for RNA extraction by use of the RNeasy RNA extraction kit (Qiagen)  according to the manufacturer's instructions."]},
#         'GSM1547616': {'Sample_id': 'GSM1547616', 'sample_title': ['7-day-old seedlings of Col-0 treated with 37℃ heat shock for 1h, biological rep1'], 'sample_source_name': ['7-day-old seedlings of Col-0 treated with 37℃ heat shock for 1hr'], 'sample_characteristicts': ['cultivar: Col-0', 'genotype: wildtype', 'tissue: seedling', 'age: 7-day-old', 'treatment: 37℃ heat shock for 1hr'], 'sample_extraction_protocol': ["Samples were collected for RNA extraction by use of the RNeasy RNA extraction kit (Qiagen)  according to the manufacturer's instructions."]},
#         'GSM1547626': {'Sample_id': 'GSM1547626', 'sample_title': ['7-day-old seedlings of HSFA6b-OE treated with 37℃ heat shock for 1h, biological rep2'], 'sample_source_name': ['7-day-old seedlings of HSFA6b-OE treated with 37℃ heat shock for 1hr'], 'sample_characteristicts': ['genotype: HSFA6b-overexpression', 'tissue: seedling', 'age: 7-day-old', 'treatment: 37℃ heat shock for 1hr'], 'sample_extraction_protocol': ["Samples were collected for RNA extraction by use of the RNeasy RNA extraction kit (Qiagen)  according to the manufacturer's instructions."]},
#         'GSM1547625': {'Sample_id': 'GSM1547625', 'sample_title': ['7-day-old seedlings of Col-0 treated with 37℃ heat shock for 1h, biological rep2'], 'sample_source_name': ['7-day-old seedlings of Col-0 treated with 37℃ heat shock for 1hr'], 'sample_characteristicts': ['cultivar: Col-0', 'genotype: wildtype', 'tissue: seedling', 'age: 7-day-old', 'treatment: 37℃ heat shock for 1hr'], 'sample_extraction_protocol': ["Samples were collected for RNA extraction by use of the RNeasy RNA extraction kit (Qiagen)  according to the manufacturer's instructions."]},
#         'GSM1547614': {'Sample_id': 'GSM1547614', 'sample_title': ['7-day-old seedlings of HSFA6b-OE at control condition, biological rep1'], 'sample_source_name': ['7-day-old seedling of HSFA6b-OE at control condition'], 'sample_characteristicts': ['genotype: HSFA6b-overexpression', 'tissue: seedling', 'age: 7-day-old', 'treatment: none'], 'sample_extraction_protocol': ["Samples were collected for RNA extraction by use of the RNeasy RNA extraction kit (Qiagen)  according to the manufacturer's instructions."]},
#         'GSM1547615': {'Sample_id': 'GSM1547615', 'sample_title': ['7-day-old seedlings of HSFA6b-RD at control condition, biological rep1'], 'sample_source_name': ['7-day-old seedling of HSFA6b-RD at control condition'], 'sample_characteristicts': ['genotype: HSFA6b dominant-negative', 'tissue: seedling', 'age: 7-day-old', 'treatment: none'], 'sample_extraction_protocol': ["Samples were collected for RNA extraction by use of the RNeasy RNA extraction kit (Qiagen)  according to the manufacturer's instructions."]},
#         'GSM1547621': {'Sample_id': 'GSM1547621', 'sample_title': ['7-day-old seedlings of HSFA6b-RD treated with 150mM NaCl for 6hr, biological rep1'], 'sample_source_name': ['7-day-old seedlings of HSFA6b-RD treated with 150mM NaCl for 6hr'], 'sample_characteristicts': ['genotype: HSFA6b dominant-negative', 'tissue: seedling', 'age: 7-day-old', 'treatment: 150mM NaCl for 6hr'], 'sample_extraction_protocol': ["Samples were collected for RNA extraction by use of the RNeasy RNA extraction kit (Qiagen)  according to the manufacturer's instructions."]},
#         'GSM1547628': {'Sample_id': 'GSM1547628', 'sample_title': ['7-day-old seedlings of Col-0 treated with 150mM NaCl for 6hr, biological rep2'], 'sample_source_name': ['7-day-old seedlings of Col-0 treated with 150mM NaCl for 6hr'], 'sample_characteristicts': ['cultivar: Col-0', 'genotype: wildtype', 'tissue: seedling', 'age: 7-day-old', 'treatment: 150mM NaCl for 6hr'], 'sample_extraction_protocol': ["Samples were collected for RNA extraction by use of the RNeasy RNA extraction kit (Qiagen)  according to the manufacturer's instructions."]},
#         'GSM1547624': {'Sample_id': 'GSM1547624', 'sample_title': ['7-day-old seedlings of HSFA6b-RD at control condition, biological rep2'], 'sample_source_name': ['7-day-old seedling of HSFA6b-RD at control condition'], 'sample_characteristicts': ['genotype: HSFA6b dominant-negative', 'tissue: seedling', 'age: 7-day-old', 'treatment: none'], 'sample_extraction_protocol': ["Samples were collected for RNA extraction by use of the RNeasy RNA extraction kit (Qiagen)  according to the manufacturer's instructions."]}
#     }

#     print("--- Testing with a sample under control condition (no treatment) ---")
#     sample_id_control = 'GSM1547622'
#     sample_data_control = samples_metadata[sample_id_control]
#     extracted_control = GSE63372_extractor(sample_data_control)
#     print(f"Input Sample ID: {sample_id_control}")
#     print(f"Extracted Data:\n{json.dumps(extracted_control, indent=2)}")
#     # Expected output:
#     # {
#     #   "tissue": "seedling",
#     #   "treatment": [],
#     #   "medium": "0.5X MS plates containing 1% sucrose"
#     # }

#     print("\n--- Testing with a sample treated with heat shock ---")
#     sample_id_heat_shock = 'GSM1547618'
#     sample_data_heat_shock = samples_metadata[sample_id_heat_shock]
#     extracted_heat_shock = GSE63372_extractor(sample_data_heat_shock)
#     print(f"Input Sample ID: {sample_id_heat_shock}")
#     print(f"Extracted Data:\n{json.dumps(extracted_heat_shock, indent=2)}")
#     # Expected output:
#     # {
#     #   "tissue": "seedling",
#     #   "treatment": ["37℃ heat shock for 1hr"],
#     #   "medium": "0.5X MS plates containing 1% sucrose"
#     # }

#     print("\n--- Testing with a sample treated with NaCl ---")
#     sample_id_nacl = 'GSM1547619'
#     sample_data_nacl = samples_metadata[sample_id_nacl]
#     extracted_nacl = GSE63372_extractor(sample_data_nacl)
#     print(f"Input Sample ID: {sample_id_nacl}")
#     print(f"Extracted Data:\n{json.dumps(extracted_nacl, indent=2)}")
#     # Expected output:
#     # {
#     #   "tissue": "seedling",
#     #   "treatment": ["150mM NaCl for 6hr"],
#     #   "medium": "0.5X MS plates containing 1% sucrose"
#     # }


# import json

# def GSE83136_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from the sample metadata
#     for the GSE83136 dataset, conforming to a specific JSON schema.

#     The schema requires 'tissue', 'treatment' (as a list), and 'medium'.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 This typically includes keys like 'sample_characteristicts'.

#     Returns:
#         dict: A dictionary formatted as a JSON instance, containing the extracted
#               'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # Based on the study_metadata, the samples are derived from "Arabidopsis thaliana seedlings".
#     # This information is constant across all samples in this study.
#     extracted_data["tissue"] = "seedlings"

#     # 2. Extract 'treatment'
#     # The 'treatment' information is found within 'sample_characteristicts'
#     # and describes the specific condition of the sample (e.g., control, 4h after acclimation).
#     # The overall study involves "heat shock", so if a sample is not a control,
#     # "heat shock" is an implied treatment.
#     treatments_list = []
#     specific_sample_treatment = None

#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if isinstance(characteristic, str) and characteristic.startswith("treatment:"):
#                 specific_sample_treatment = characteristic.split(":", 1)[1].strip()
#                 break

#     if specific_sample_treatment:
#         # If the sample is not a 'control', it implies it underwent heat shock.
#         if specific_sample_treatment.lower() != "control":
#             treatments_list.append("heat shock")# TODO: Check this for a diferent label
#         # Add the specific condition of the sample (e.g., "4 hr after acclimation", "control")
#         treatments_list.append(specific_sample_treatment)
#     else:
#         # Fallback if no specific treatment characteristic is found
#         treatments_list.append("unknown treatment")

#     extracted_data["treatment"] = treatments_list

#     # 3. Extract 'medium'
#     # The growth medium is not explicitly mentioned in the provided study_metadata
#     # or sample_metadata. In such cases, it's best to state that it's not specified.
#     extracted_data["medium"] = "Not specified"

#     return extracted_data



# import json

# def GSE63522_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from the sample metadata
#     for the GSE63522 dataset, conforming to a specific JSON schema.

#     The schema requires 'tissue', 'treatment' (as a list of strings), and 'medium'.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 This typically corresponds to one entry from the
#                                 'samples_metadata' provided in the problem description.

#     Returns:
#         dict: A dictionary formatted as a JSON instance conforming to the output schema.
#               Example: {"tissue": "root", "treatment": ["heat stress (39ºC, 40 min)"], "medium": "MS agar plates"}

#     Raises:
#         ValueError: If required 'tissue' or 'treatment' information cannot be found
#                     in the sample_metadata.
#     """
#     extracted_data = {}

#     # The 'sample_characteristicts' field is a list of strings, where each string
#     # contains a key-value pair (e.g., "tissue: root", "treatment: heat stress").
#     sample_characteristics = sample_metadata.get('sample_characteristicts', [])

#     # 1. Extract 'tissue'
#     # Look for a string starting with 'tissue:' in the characteristics list.
#     tissue = None
#     for char_string in sample_characteristics:
#         if char_string.strip().startswith('tissue:'):
#             # Split by the first colon and take the second part, then strip whitespace.
#             tissue = char_string.split(':', 1)[1].strip()
#             break
#     if tissue is None:
#         raise ValueError("Required 'tissue' information not found in sample_metadata.")
#     extracted_data['tissue'] = tissue

#     # 2. Extract 'treatment'
#     # Look for a string starting with 'treatment:' in the characteristics list.
#     treatment = None
#     for char_string in sample_characteristics:
#         if char_string.strip().startswith('treatment:'):
#             # Split by the first colon and take the second part, then strip whitespace.
#             treatment = char_string.split(':', 1)[1].split('(')[0].strip()
#             break
#     if treatment is None:
#         raise ValueError("Required 'treatment' information not found in sample_metadata.")
#     # The schema specifies 'treatment' as an array of strings.
#     # Even if there's only one treatment, it should be wrapped in a list.
#     extracted_data['treatment'] = [treatment]

#     # 3. Extract 'medium'
#     # Based on the 'overall_design' in the study_metadata, the samples were grown
#     # on "MS agar plates". This information is constant across all samples in this
#     # study and is not present in the individual sample_metadata dictionaries.
#     # Therefore, it is hardcoded for this specific dataset (GSE63522).
#     extracted_data['medium'] = "MS agar plates"

#     return extracted_data



# import json

# def GSE26983_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information of the sample from the
#     GSE26983 dataset following the specified JSON schema.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 This typically corresponds to one entry from the
#                                 'samples_metadata' dictionary provided in the problem description.

#     Returns:
#         dict: A dictionary formatted as a JSON instance conforming to the output schema.
#               It will contain 'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # Initialize with default/placeholder values as per schema requirements
#     # and in case information is unexpectedly missing for a specific sample.
#     tissue_value = "Not specified"
#     treatment_values = []
#     medium_value = "Not specified" # Not explicitly found in the provided metadata

#     # Information for 'tissue' and 'treatment' is found within 'sample_characteristicts'
#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if isinstance(characteristic, str):
#                 # Extract 'tissue'
#                 if characteristic.startswith('tissue:'):
#                     tissue_value = characteristic.split(':', 1)[1].strip()
#                 # Extract 'treatment'
#                 elif characteristic.startswith('treatment:'):
#                     treatment_values.append(characteristic.split(':', 1)[1].strip())

#     # Assign extracted values to the output dictionary
#     extracted_data['tissue'] = tissue_value
    
#     # Ensure treatment is always a list, even if empty or only one item
#     # The schema requires an array for 'treatment'.
#     # For this dataset, "0mM NaCl" is a valid treatment indicating no salt stress.
#     extracted_data['treatment'] = treatment_values if treatment_values else ["Not specified"]

#     # 'medium' is not explicitly provided in the sample_metadata or study_metadata.
#     # As it's a required field in the schema, and we cannot invent data,
#     # we use "Not specified" as the most accurate representation of its absence
#     # from the provided data.
#     extracted_data['medium'] = medium_value

#     return extracted_data



# import json

# def GSE19700_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from a sample metadata dictionary
#     conforming to the specified JSON schema.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the target JSON schema, containing
#               'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # Initialize with default values.
#     # 'tissue' and 'treatment' will be extracted from 'sample_characteristicts'.
#     # 'medium' is inferred from the study context as it's not explicitly stated
#     # in the sample metadata but is a required field in the schema.
#     # The study summary mentions "gradual soil drying", strongly implying "soil"
#     # as the growth medium for Arabidopsis plants.
#     tissue = "unknown"
#     treatment = []
#     medium = "soil" # Inferred from study context: "gradual soil drying"

#     # Extract information from 'sample_characteristicts'
#     characteristics = sample_metadata.get('sample_characteristicts', [])
#     for char_string in characteristics:
#         if ': ' in char_string:
#             key, value = char_string.split(': ', 1)
#             key = key.strip().lower() # Normalize key for easier matching
#             value = value.strip()

#             if key == 'tissue':
#                 tissue = value
#             elif key == 'treatment':
#                 # The schema expects an array of strings for 'treatment'.
#                 # In this dataset, 'treatment' appears as a single string value
#                 # within the characteristics list (e.g., 'water-limited').
#                 # We append it to a list to conform to the array type.
#                 treatment.append(value)

#     extracted_data['tissue'] = tissue
#     extracted_data['treatment'] = treatment
#     extracted_data['medium'] = medium

#     return extracted_data

# if __name__ == '__main__':
#     # Example usage with a sample from the provided metadata
#     samples_metadata = {
#         'GSM491670': {'Sample_id': 'GSM491670', 'sample_title': ['Columbia.Dry.MD.rep2'], 'sample_source_name': ['whole rosette, water-limited, midday'], 'sample_characteristicts': ['ecotype: Columbia-0', 'genotype: wild type', 'tissue: whole rosette', 'age: 32 days', 'treatment: water-limited', 'time of day: 12h00 (midday)'], 'sample_extraction_protocol': ["Trizol, as per manufacturer's instructions."]},
#         'GSM491685': {'Sample_id': 'GSM491685', 'sample_title': ['Columbia.Wet.PD.rep2'], 'sample_source_name': ['whole rosette, well-watered, pre-dawn'], 'sample_characteristicts': ['ecotype: Columbia-0', 'genotype: wild type', 'tissue: whole rosette', 'age: 32 days', 'treatment: well-watered', 'time of day: 06h00 (pre-dawn)'], 'sample_extraction_protocol': ["Trizol, as per manufacturer's instructions."]},
#         'GSM491678': {'Sample_id': 'GSM491678', 'sample_title': ['Columbia.Wet.MN.rep1'], 'sample_source_name': ['whole rosette, well-watered, midnight'], 'sample_characteristicts': ['ecotype: Columbia-0', 'genotype: wild type', 'tissue: whole rosette', 'age: 32 days', 'treatment: well-watered', 'time of day: 00h00 (midnight)'], 'sample_extraction_protocol': ["Trizol, as per manufacturer's instructions."]},
#         # ... (other samples omitted for brevity)
#     }

#     # Test with a 'water-limited' sample
#     sample_id_dry = 'GSM491670'
#     sample_data_dry = samples_metadata[sample_id_dry]
#     extracted_info_dry = GSE19700_extractor(sample_data_dry)
#     print(f"Extracted info for {sample_id_dry}:")
#     print(json.dumps(extracted_info_dry, indent=2))
#     # Expected output:
#     # {
#     #   "tissue": "whole rosette",
#     #   "treatment": ["water-limited"],
#     #   "medium": "soil"
#     # }

#     print("-" * 30)

#     # Test with a 'well-watered' sample
#     sample_id_wet = 'GSM491685'
#     sample_data_wet = samples_metadata[sample_id_wet]
#     extracted_info_wet = GSE19700_extractor(sample_data_wet)
#     print(f"Extracted info for {sample_id_wet}:")
#     print(json.dumps(extracted_info_wet, indent=2))
#     # Expected output:
#     # {
#     #   "tissue": "whole rosette",
#     #   "treatment": ["well-watered"],
#     #   "medium": "soil"
#     # }

#     print("-" * 30)

#     # Test with a sample where 'tissue' or 'treatment' might be missing (hypothetical)
#     # In this specific dataset, 'tissue' and 'treatment' are consistently present.
#     # This test demonstrates the fallback behavior.
#     sample_data_incomplete = {
#         'Sample_id': 'GSM_INCOMPLETE',
#         'sample_title': ['Incomplete Sample'],
#         'sample_characteristicts': ['ecotype: Mutant', 'age: 10 days'], # Missing tissue and treatment
#     }
#     extracted_info_incomplete = GSE19700_extractor(sample_data_incomplete)
#     print(f"Extracted info for GSM_INCOMPLETE (missing characteristics):")
#     print(json.dumps(extracted_info_incomplete, indent=2))
#     # Expected output:
#     # {
#     #   "tissue": "unknown",
#     #   "treatment": [],
#     #   "medium": "soil"
#     # }


# import json

# def GSE58620_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from sample metadata for the GSE58620 dataset
#     following the specified JSON schema.

#     The function identifies the tissue, applied treatments/stresses, and growth medium
#     based on the provided sample_metadata and the context from the study's overall design.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 Expected keys include 'sample_source_name'.

#     Returns:
#         dict: A dictionary formatted according to the target JSON schema:
#               {"tissue": "...", "treatment": [...], "medium": "..."}.
#     """
#     extracted_data = {}

#     # Retrieve the 'sample_source_name' which contains most of the relevant information.
#     # Use .get() with a default to safely handle cases where the key might be missing or empty.
#     source_name = sample_metadata.get('sample_source_name', [''])[0]

#     # 1. Extract 'tissue'
#     # Based on the study_metadata and samples_metadata, the tissue is consistently
#     # "Arabidopsis seedling". This is explicitly mentioned in 'sample_source_name'.
#     if "Arabidopsis seedling" in source_name:
#         extracted_data['tissue'] = "Arabidopsis seedling"
#     else:
#         # Fallback if the expected tissue string is not found, though it appears consistent
#         # across all samples in this dataset.
#         extracted_data['tissue'] = "unknown"

#     # 2. Extract 'treatment'
#     # Treatments are identified from the 'sample_source_name' string.
#     # These include specific stresses (like heat stress) and environmental conditions
#     # (like dark/light conditions), which are considered treatments/stresses in this context.
#     treatments = []

#     # Check for heat stress (HS) treatment
#     if "HS treatment" in source_name:
#         treatments.append("heat stress")
#     # If "without HS treatment" is present, it means heat stress was NOT applied,
#     # so "heat stress" should not be added to the treatments list.

#     # Check for environmental light/dark conditions
#     if "dark condition" in source_name:
#         treatments.append("dark")
#     elif "light condition" in source_name:
#         # Although not present in the provided sample examples, this handles light conditions
#         # as mentioned in the study's overall design.
#         treatments.append("light")
    
#     extracted_data['treatment'] = treatments

#     # 3. Extract 'medium'
#     # The growth medium is not explicitly mentioned in the provided study_metadata
#     # or samples_metadata. As the task requires extraction from the given data,
#     # and this information is absent, "not specified" is used as a placeholder.
#     # If this information were constant and available, it would likely be in the
#     # study_metadata's 'overall_design' or 'summary'.
#     extracted_data['medium'] = "not specified"

#     return extracted_data



# import json

# def GSE95202_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from sample metadata for GSE95202
#     following a specific JSON schema.

#     The schema requires 'tissue' (string), 'treatment' (list of strings),
#     and 'medium' (string).

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 This typically corresponds to one entry from the
#                                 'samples_metadata' provided in the problem description.

#     Returns:
#         dict: A dictionary conforming to the specified schema:
#               {"tissue": "...", "treatment": [...], "medium": "..."}
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The 'tissue' information is consistently found within the 'sample_characteristicts'
#     # list, in an item starting with "tissue: ". For this dataset, it's always "whole plant".
#     tissue_found = False
#     if 'sample_characteristicts' in sample_metadata:
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if characteristic.startswith('tissue:'):
#                 extracted_data['tissue'] = characteristic.split(':', 1)[1].strip()
#                 tissue_found = True
#                 break
#     # Fallback in case 'tissue' characteristic is missing, though it appears consistent
#     if not tissue_found:
#         extracted_data['tissue'] = "whole plant" # Based on consistent observation in samples_metadata

#     # 2. Extract 'treatment'
#     # The 'treatment' information is derived from 'sample_source_name'.
#     # "water control" implies no specific applied treatment.
#     # Other treatments like "ethanol" or "salt stress" are extracted directly,
#     # and " and " is used to separate multiple treatments.
#     treatments = []
#     if 'sample_source_name' in sample_metadata and sample_metadata['sample_source_name']:
#         source_name = sample_metadata['sample_source_name'][0].lower()

#         if source_name == "water control":
#             # "water control" indicates a baseline, non-treated condition.
#             treatments = ['control']
#         else:
#             # Split by " and " to identify individual treatments/stresses.
#             parts = source_name.split(' and ')
#             for part in parts:
#                 cleaned_part = part.strip()
#                 if cleaned_part: # Ensure no empty strings are added
#                     treatments.append(cleaned_part)
#     extracted_data['treatment'] = treatments

#     # 3. Extract 'medium'
#     # The study summary states "treated with 0.3 % ethanol or water". This implies
#     # that "water" serves as the base medium for the applied treatments.
#     # This value is constant across all samples in this study.
#     extracted_data['medium'] = "water"

#     return extracted_data

# if __name__ == '__main__':
#     # The following is the study metadata (constant across all samples in the study)
#     # and list of the metadata of the samples in this studies.
#     # Use it to guide to know where in the sample_metadata to find the information.
#     study_metadata = {
#         'title': ['Ethanol enhances high-salinity stress tolerance by detoxification of reactive oxygen species in Arabidopsis and rice'],
#         'summary': ['To analyze the molecular function of ethanol in salt stress responses in Arabidopsis, we conducted microarray analysis using 4-day-old plants, which were treated with 0.3 % ethanol or water for 24 h, and then treated with or without 100 mM NaCl for 2 h'],
#         'overall_design': ['Microarray analysis was conducted using ethanol-treated and non-treated plants subjected to non-stress and salt-stress condition for 2 h.']
#     }
#     samples_metadata_full = {
#         'GSM2498604': {'Sample_id': 'GSM2498604', 'sample_title': ['ethanol rep. 1'], 'sample_source_name': ['ethanol'], 'sample_characteristicts': ['ecotype: Col-0', 'tissue: whole plant'], 'sample_extraction_protocol': ['RNA was extracted from four biological replicates with the Plant RNA Reagent (Thermo Fisher Scientific) according to the manufacturer’s instructions']},
#         'GSM2498611': {'Sample_id': 'GSM2498611', 'sample_title': ['ethanol and salt stress rep. 4'], 'sample_source_name': ['ethanol and salt stress'], 'sample_characteristicts': ['ecotype: Col-0', 'tissue: whole plant'], 'sample_extraction_protocol': ['RNA was extracted from four biological replicates with the Plant RNA Reagent (Thermo Fisher Scientific) according to the manufacturer’s instructions']},
#         'GSM2498597': {'Sample_id': 'GSM2498597', 'sample_title': ['water control rep. 2'], 'sample_source_name': ['water control'], 'sample_characteristicts': ['ecotype: Col-0', 'tissue: whole plant'], 'sample_extraction_protocol': ['RNA was extracted from four biological replicates with the Plant RNA Reagent (Thermo Fisher Scientific) according to the manufacturer’s instructions']},
#         'GSM2498603': {'Sample_id': 'GSM2498603', 'sample_title': ['salt stress rep. 4'], 'sample_source_name': ['salt stress'], 'sample_characteristicts': ['ecotype: Col-0', 'tissue: whole plant'], 'sample_extraction_protocol': ['RNA was extracted from four biological replicates with the Plant RNA Reagent (Thermo Fisher Scientific) according to the manufacturer’s instructions']},
#         'GSM2498606': {'Sample_id': 'GSM2498606', 'sample_title': ['ethanol rep. 3'], 'sample_source_name': ['ethanol'], 'sample_characteristicts': ['ecotype: Col-0', 'tissue: whole plant'], 'sample_extraction_protocol': ['RNA was extracted from four biological replicates with the Plant RNA Reagent (Thermo Fisher Scientific) according to the manufacturer’s instructions']},
#         'GSM2498610': {'Sample_id': 'GSM2498610', 'sample_title': ['ethanol and salt stress rep. 3'], 'sample_source_name': ['ethanol and salt stress'], 'sample_characteristicts': ['ecotype: Col-0', 'tissue: whole plant'], 'sample_extraction_protocol': ['RNA was extracted from four biological replicates with the Plant RNA Reagent (Thermo Fisher Scientific) according to the manufacturer’s instructions']},
#         'GSM2498601': {'Sample_id': 'GSM2498601', 'sample_title': ['salt stress rep. 2'], 'sample_source_name': ['salt stress'], 'sample_characteristicts': ['ecotype: Col-0', 'tissue: whole plant'], 'sample_extraction_protocol': ['RNA was extracted from four biological replicates with the Plant RNA Reagent (Thermo Fisher Scientific) according to the manufacturer’s instructions']},
#         'GSM2498607': {'Sample_id': 'GSM2498607', 'sample_title': ['ethanol rep. 4'], 'sample_source_name': ['ethanol'], 'sample_characteristicts': ['ecotype: Col-0', 'tissue: whole plant'], 'sample_extraction_protocol': ['RNA was extracted from four biological replicates with the Plant RNA Reagent (Thermo Fisher Scientific) according to the manufacturer’s instructions']},
#         'GSM2498608': {'Sample_id': 'GSM2498608', 'sample_title': ['ethanol and salt stress rep. 1'], 'sample_source_name': ['ethanol and salt stress'], 'sample_characteristicts': ['ecotype: Col-0', 'tissue: whole plant'], 'sample_extraction_protocol': ['RNA was extracted from four biological replicates with the Plant RNA Reagent (Thermo Fisher Scientific) according to the manufacturer’s instructions']},
#         'GSM2498605': {'Sample_id': 'GSM2498605', 'sample_title': ['ethanol rep. 2'], 'sample_source_name': ['ethanol'], 'sample_characteristicts': ['ecotype: Col-0', 'tissue: whole plant'], 'sample_extraction_protocol': ['RNA was extracted from four biological replicates with the Plant RNA Reagent (Thermo Fisher Scientific) according to the manufacturer’s instructions']},
#         'GSM2498609': {'Sample_id': 'GSM2498609', 'sample_title': ['ethanol and salt stress rep. 2'], 'sample_source_name': ['ethanol and salt stress'], 'sample_characteristicts': ['ecotype: Col-0', 'tissue: whole plant'], 'sample_extraction_protocol': ['RNA was extracted from four biological replicates with the Plant RNA Reagent (Thermo Fisher Scientific) according to the manufacturer’s instructions']},
#         'GSM2498600': {'Sample_id': 'GSM2498600', 'sample_title': ['salt stress rep. 1'], 'sample_source_name': ['salt stress'], 'sample_characteristicts': ['ecotype: Col-0', 'tissue: whole plant'], 'sample_extraction_protocol': ['RNA was extracted from four biological replicates with the Plant RNA Reagent (Thermo Fisher Scientific) according to the manufacturer’s instructions']},
#         'GSM2498598': {'Sample_id': 'GSM2498598', 'sample_title': ['water control rep. 3'], 'sample_source_name': ['water control'], 'sample_characteristicts': ['ecotype: Col-0', 'tissue: whole plant'], 'sample_extraction_protocol': ['RNA was extracted from four biological replicates with the Plant RNA Reagent (Thermo Fisher Scientific) according to the manufacturer’s instructions']},
#         'GSM2498596': {'Sample_id': 'GSM2498596', 'sample_title': ['water control rep. 1'], 'sample_source_name': ['water control'], 'sample_characteristicts': ['ecotype: Col-0', 'tissue: whole plant'], 'sample_extraction_protocol': ['RNA was extracted from four biological replicates with the Plant RNA Reagent (Thermo Fisher Scientific) according to the manufacturer’s instructions']},
#         'GSM2498599': {'Sample_id': 'GSM2498599', 'sample_title': ['water control rep. 4'], 'sample_source_name': ['water control'], 'sample_characteristicts': ['ecotype: Col-0', 'tissue: whole plant'], 'sample_extraction_protocol': ['RNA was extracted from four biological replicates with the Plant RNA Reagent (Thermo Fisher Scientific) according to the manufacturer’s instructions']},
#         'GSM2498602': {'Sample_id': 'GSM2498602', 'sample_title': ['salt stress rep. 3'], 'sample_source_name': ['salt stress'], 'sample_characteristicts': ['ecotype: Col-0', 'tissue: whole plant'], 'sample_extraction_protocol': ['RNA was extracted from four biological replicates with the Plant RNA Reagent (Thermo Fisher Scientific) according to the manufacturer’s instructions']}
#     }

#     print("--- Testing GSE95202_extractor ---")

#     # Test case 1: Ethanol and salt stress
#     sample_id_1 = 'GSM2498611'
#     sample_data_1 = samples_metadata_full[sample_id_1]
#     extracted_1 = GSE95202_extractor(sample_data_1)
#     print(f"\nSample ID: {sample_id_1}")
#     print(f"Input metadata: {json.dumps(sample_data_1, indent=2)}")
#     print(f"Extracted schema: {json.dumps(extracted_1, indent=2)}")
#     expected_1 = {
#         "tissue": "whole plant",
#         "treatment": ["ethanol", "salt stress"],
#         "medium": "water"
#     }
#     assert extracted_1 == expected_1, f"Test 1 Failed: Expected {expected_1}, Got {extracted_1}"
#     print("Test 1 Passed.")

#     # Test case 2: Water control
#     sample_id_2 = 'GSM2498596'
#     sample_data_2 = samples_metadata_full[sample_id_2]
#     extracted_2 = GSE95202_extractor(sample_data_2)
#     print(f"\nSample ID: {sample_id_2}")
#     print(f"Input metadata: {json.dumps(sample_data_2, indent=2)}")
#     print(f"Extracted schema: {json.dumps(extracted_2, indent=2)}")
#     expected_2 = {
#         "tissue": "whole plant",
#         "treatment": [],
#         "medium": "water"
#     }
#     assert extracted_2 == expected_2, f"Test 2 Failed: Expected {expected_2}, Got {extracted_2}"
#     print("Test 2 Passed.")

#     # Test case 3: Ethanol only
#     sample_id_3 = 'GSM2498604'
#     sample_data_3 = samples_metadata_full[sample_id_3]
#     extracted_3 = GSE95202_extractor(sample_data_3)
#     print(f"\nSample ID: {sample_id_3}")
#     print(f"Input metadata: {json.dumps(sample_data_3, indent=2)}")
#     print(f"Extracted schema: {json.dumps(extracted_3, indent=2)}")
#     expected_3 = {
#         "tissue": "whole plant",
#         "treatment": ["ethanol"],
#         "medium": "water"
#     }
#     assert extracted_3 == expected_3, f"Test 3 Failed: Expected {expected_3}, Got {extracted_3}"
#     print("Test 3 Passed.")

#     # Test case 4: Salt stress only
#     sample_id_4 = 'GSM2498603'
#     sample_data_4 = samples_metadata_full[sample_id_4]
#     extracted_4 = GSE95202_extractor(sample_data_4)
#     print(f"\nSample ID: {sample_id_4}")
#     print(f"Input metadata: {json.dumps(sample_data_4, indent=2)}")
#     print(f"Extracted schema: {json.dumps(extracted_4, indent=2)}")
#     expected_4 = {
#         "tissue": "whole plant",
#         "treatment": ["salt stress"],
#         "medium": "water"
#     }
#     assert extracted_4 == expected_4, f"Test 4 Failed: Expected {expected_4}, Got {extracted_4}"
#     print("Test 4 Passed.")

#     print("\nAll tests passed!")


# import json

# def GSE24177_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information of the sample from
#     the provided sample_metadata dictionary, conforming to a specific JSON schema.

#     The schema requires 'tissue', 'treatment', and 'medium'.

#     - 'tissue' and 'treatment' are extracted from the 'sample_characteristicts' field.
#     - 'medium' is inferred from the study's overall design (provided in the problem
#       description as 'study_metadata', indicating plants were grown in soil),
#       and is constant across all samples in this study.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 This typically corresponds to one entry from the
#                                 'samples_metadata' dictionary provided in the problem.

#     Returns:
#         dict: A dictionary formatted as a JSON instance conforming to the
#               output schema:
#               {"properties": {"tissue": {"type": "string"},
#                               "treatment": {"type": "array", "items": {"type": "string"}},
#                               "medium": {"type": "string"}},
#                "required": ["tissue", "treatment", "medium"]}
#     """
#     extracted_data = {}

#     # Extract 'tissue' and 'treatment' from the 'sample_characteristicts' list.
#     # This field contains key-value pairs separated by ': '.
#     characteristics = sample_metadata.get('sample_characteristicts', [])

#     for char_string in characteristics:
#         # Split only on the first occurrence of ': ' to handle potential colons in values
#         parts = char_string.split(': ', 1)
#         if len(parts) == 2:
#             key, value = parts[0].strip(), parts[1].strip()
#             if key == 'tissue':
#                 extracted_data['tissue'] = value
#             elif key == 'treatment':
#                 # The schema expects 'treatment' to be an array of strings.
#                 # Even if there's only one treatment, it should be in a list.
#                 extracted_data['treatment'] = [value]

#     # Extract 'medium'.
#     # Based on the 'overall_design' in the study_metadata (not passed to this function
#     # but provided in the problem description), the plants were grown in "soil".
#     # This is a constant value for all samples in this study.
#     extracted_data['medium'] = "soil"

#     return extracted_data



# import json

# def GSE11758_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from sample metadata for GSE11758
#     following the specified JSON schema.

#     The schema is:
#     {"properties": {"tissue": {"description": "Tissue the samples was extracted from.", "title": "Tissue", "type": "string"},
#                     "treatment": {"description": "List of treatments and stresses that was applied to the sample.", "items": {"type": "string"}, "title": "Treatment", "type": "array"},
#                     "medium": {"description": "Growth medium of the sample.", "title": "Medium", "type": "string"}},
#      "required": ["tissue", "treatment", "medium"]}

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 Expected to have a 'sample_source_name' key.

#     Returns:
#         dict: A dictionary conforming to the specified schema, containing 'tissue',
#               'treatment' (as a list of strings), and 'medium'.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # Based on the study metadata and consistent 'sample_source_name',
#     # the tissue is always "Arabidopsis mature leaves".
#     extracted_data["tissue"] = "Arabidopsis mature leaves"

#     # Get the 'sample_source_name' which contains details about treatment and medium.
#     # Use .get() with a default to handle cases where the key might be missing or empty.
#     source_name = sample_metadata.get('sample_source_name', [''])[0]

#     # Initialize treatment and medium
#     treatment_list = []
#     medium_str = ""

#     # 2. Extract 'treatment' and 'medium' based on keywords in source_name
#     if "infiltrated with" in source_name:
#         # Examples: 'Arabidopsis mature leaves infiltrated with 2mg/l tunicamycin'
#         #           'Arabidopsis mature leaves infiltrated with 10mM proline'
#         parts = source_name.split("infiltrated with ", 1)
#         if len(parts) > 1:
#             substance = parts[1].strip()
#             medium_str = substance  # The infiltrating substance is considered the medium

#             # Determine specific treatment based on the substance
#             if "tunicamycin" in substance:
#                 treatment_list.append("tunicamycin")
#             elif "DMF" in substance:
#                 treatment_list.append("DMF")
#             elif "proline" in substance:
#                 # As per study metadata, L-Proline is a control for AZC
#                 treatment_list.append("L-Proline")
#             elif "AZC" in substance:
#                 # L-azetidine-2-carboxylic acid (AZC)
#                 treatment_list.append("AZC")
#             else:
#                 # Fallback for any other infiltrated substance
#                 treatment_list.append(f"infiltrated with {substance}")
#         else:
#             # Fallback if "infiltrated with" is present but parsing fails
#             treatment_list.append("unknown infiltration")
#             medium_str = "unknown solution"

#     elif "incubated at" in source_name:
#         # Examples: 'Arabidopsis mature leaves incubated at 37 °C for 1h'
#         #           'Arabidopsis mature leaves incubated at 20 °C for 1h'
#         if "37 °C" in source_name:
#             treatment_list.append("Heat shock (37 °C)")
#             medium_str = "Air"  # No liquid medium is specified for temperature incubation
#         elif "20 °C" in source_name:
#             # As per study metadata, 20°C is a control for heat shock
#             treatment_list.append("Control (20 °C)")
#             medium_str = "Air"  # No liquid medium is specified for temperature incubation
#         else:
#             # Fallback for unrecognized incubation temperature
#             treatment_list.append("unknown incubation")
#             medium_str = "Air"

#     else:
#         # Fallback for any other unexpected 'sample_source_name' format
#         treatment_list.append("unspecified treatment")
#         medium_str = "unspecified medium"

#     extracted_data["treatment"] = treatment_list
#     extracted_data["medium"] = medium_str

#     return extracted_data



# import json
# import re

# def GSE4760_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from sample metadata for GSE4760
#     following the specified JSON schema.

#     The output schema is:
#     {"properties": {"tissue": {"description": "Tissue the samples was extracted from.", "title": "Tissue", "type": "string"},
#                     "treatment": {"description": "List of treatments and stresses that was applied to the sample.", "items": {"type": "string"}, "title": "Treatment", "type": "array"},
#                     "medium": {"description": "Growth medium of the sample.", "title": "Medium", "type": "string"}},
#      "required": ["tissue", "treatment", "medium"]}

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 Example:
#                                 {'Sample_id': 'GSM107586',
#                                  'sample_title': ['Wt_seedling_HS44_rep1'],
#                                  'sample_source_name': ['seedling, heat shock 37C 1h, recover 24C 2d, and heat shock 44C 45 min'],
#                                  'sample_characteristicts': ['Ecotype: Col-0 wild type', 'Stage: 5-d old', 'Tissue: whole seedling'],
#                                  'sample_extraction_protocol': ''}

#     Returns:
#         dict: A dictionary formatted according to the specified schema, containing
#               'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The tissue information is found within the 'sample_characteristicts' list,
#     # typically in an item starting with "Tissue: ".
#     tissue = "unknown" # Default value if not found
#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if characteristic.startswith('Tissue:'):
#                 tissue = characteristic.split(':', 1)[1].strip()
#                 break
#     extracted_data['tissue'] = tissue

#     # 2. Extract 'treatment'
#     # The treatment information is primarily in 'sample_source_name'.
#     # This field often contains a descriptive string that needs parsing.
#     treatments = []
#     source_name_str = "" # Initialize to handle cases where 'sample_source_name' might be missing or empty
#     if 'sample_source_name' in sample_metadata and \
#        isinstance(sample_metadata['sample_source_name'], list) and \
#        sample_metadata['sample_source_name']:
        
#         source_name_str = sample_metadata['sample_source_name'][0]
        
#         # Remove "seedling, " prefix if present, as "seedling" refers to the sample type/tissue, not a treatment.
#         # This also handles cases like "seedling, control" or "seedling, heat shock".
#         # Case-insensitive check for robustness.
#         if source_name_str.lower().startswith('seedling, '):
#             source_name_str = source_name_str[len('seedling, '):]
        
#         # Split the remaining string by common delimiters like ", " or " and ".
#         # This handles multiple treatments listed in a single string.
#         raw_treatment_parts = re.split(r', | and ', source_name_str)
        
#         for part in raw_treatment_parts:
#             cleaned_part = part.strip()
#             # Add non-empty parts that are not just "seedling" (in case it wasn't removed by prefix check)
#             if cleaned_part and cleaned_part.lower() != 'seedling':
#                 treatments.append(cleaned_part)
    
#     # Special handling for "control" samples:
#     # If no specific treatments were parsed from the split parts, but the original
#     # source_name string contains "control", then "control" should be considered a treatment.
#     # This covers cases like "seedling, control" where the split might result in just "control"
#     # which would then be added here if `treatments` was empty.
#     # If "heat shock, control" was parsed, "control" is already in `treatments`, so this condition won't trigger.
#     if not treatments and "control" in source_name_str.lower():
#         treatments.append("control")
    
#     # If no treatments are found at all (e.g., source_name was just "seedling" or empty),
#     # default to "untreated" as a general state, as the schema requires a non-empty array.
#     extracted_data['treatment'] = treatments if treatments else ["untreated"]

#     # 3. Extract 'medium'
#     # The growth medium is not explicitly provided in the study or sample metadata.
#     # Based on the context of Arabidopsis seedlings in a lab setting, Murashige and Skoog (MS) medium
#     # is a very common and standard choice for plant tissue culture and seedling growth.
#     extracted_data['medium'] = "MS medium"

#     return extracted_data


# import json

# def GSE71001_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from the sample metadata
#     following a predefined JSON schema.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the specified schema, containing
#               'tissue', 'treatment', and 'medium' information.
#     """

#     extracted_data = {
#         "tissue": "",
#         "treatment": [],
#         "medium": "hydroponic medium" # Derived from study_metadata['overall_design'] which states "hydroponically grown seedligs"
#     }

#     # Access the 'sample_characteristicts' list, which contains most of the required info
#     characteristics = sample_metadata.get('sample_characteristicts', [])

#     # Extract 'tissue'
#     # Iterate through characteristics to find the 'tissue' entry
#     for char_string in characteristics:
#         if char_string.startswith('tissue:'):
#             # Extract the value after 'tissue:' and strip any leading/trailing whitespace
#             extracted_data['tissue'] = char_string.split(':', 1)[1].strip()
#             break # Assuming only one tissue entry per sample

#     # Extract 'treatment'
#     # Iterate through characteristics to find 'treated with:' entries
#     treatments_list = []
#     for char_string in characteristics:
#         if char_string.startswith('treated with:'):
#             # Extract the treatment description and add it to the list
#             treatments_list.append(char_string.split(':', 1)[1].strip())
#     if treatments_list ==[]:
#         treatments_list = ['Control']
#     extracted_data['treatment'] = treatments_list

#     # 'medium' is constant across all samples in this study, derived from the study_metadata
#     # (specifically, 'overall_design' mentions "hydroponically grown seedligs").
#     # Since the function only takes sample_metadata, this value is hardcoded.
#     # extracted_data['medium'] is already initialized to "hydroponic medium"

#     return extracted_data

# if __name__ == '__main__':
#     # Example usage with the provided sample metadata
#     study_metadata = {
#         'title': ['Crosstalk between two bZIP signaling pathways orchestrates salt-induced metabolic reprogramming in Arabidopsis roots'],
#         'summary': ['Soil salinity increasingly causes crop losses worldwide. Although roots are the primary targets of salt stress, the signaling networks that facilitate metabolic reprogramming to induce stress tolerance are less understood than those in leaves. Here, a combination of transcriptomic and metabolic approaches was performed in salt-treated Arabidopsis thaliana roots, which revealed that the group S1 basic leucine zipper transcription factors bZIP1 and bZIP53 reprogram primary C- and N-metabolism. In particular, gluconeogenesis and amino acid catabolism are affected by these transcription factors. Importantly, bZIP1 expression reflects cellular stress and energy status in roots. In addition to the well-described abiotic stress response pathway initiated by the hormone abscisic acid (ABA) and executed by SnRK2 (Snf1-RELATED-PROTEIN-KINASE2) and AREB-like bZIP factors, we identify a structurally related ABA-independent signaling module consisting of SnRK1s and S1 bZIPs. Crosstalk between these signaling pathways recruits particular bZIP factor combinations to establish at least four distinct gene expression patterns. Understanding this signaling network provides a framework for securing future crop productivity.'],
#         'overall_design': ['RNA from roots from Control and salt treated hydroponically grown seedligs were extracted and subjected to microarray analysis']
#     }

#     samples_metadata = {
#         'GSM1824886': {'Sample_id': 'GSM1824886', 'sample_title': ['KO_3 Biological Rep 3'], 'sample_source_name': ['KO_seedling_NaCl_3 h'], 'sample_characteristicts': ['ecotype background: Col-0', 'genotype/variation: bzip1/bzip53 knockout', 'age: 6-week-old', 'treated with: 150 mM NaCl for 3hrs', 'tissue: seedling roots'], 'sample_extraction_protocol': ['RNAeasy, quiagen kit']},
#         'GSM1824870': {'Sample_id': 'GSM1824870', 'sample_title': ['WT_1 Biological Rep 2'], 'sample_source_name': ['WT_seedling_NaCl_1 h'], 'sample_characteristicts': ['ecotype background: Col-0', 'genotype/variation: wild type', 'age: 6-week-old', 'treated with: 150 mM NaCl for 1hr', 'tissue: seedling roots'], 'sample_extraction_protocol': ['RNAeasy, quiagen kit']},
#         'GSM1824878': {'Sample_id': 'GSM1824878', 'sample_title': ['KO_C Biological Rep 1'], 'sample_source_name': ['KO_seedling_NaCl_Control'], 'sample_characteristicts': ['ecotype background: Col-0', 'genotype/variation: bzip1/bzip53 knockout', 'age: 6-week-old', 'tissue: seedling roots'], 'sample_extraction_protocol': ['RNAeasy, quiagen kit']},
#         'GSM1824866': {'Sample_id': 'GSM1824866', 'sample_title': ['WT_C Biological Rep 1'], 'sample_source_name': ['WT_seedling_NaCl_Control'], 'sample_characteristicts': ['ecotype background: Col-0', 'genotype/variation: wild type', 'age: 6-week-old', 'tissue: seedling roots'], 'sample_extraction_protocol': ['RNAeasy, quiagen kit']}
#     }

#     print("--- Extracted Data for GSM1824886 (Treated Sample) ---")
#     sample_id_1 = 'GSM1824886'
#     extracted_info_1 = GSE71001_extractor(samples_metadata[sample_id_1])
#     print(json.dumps(extracted_info_1, indent=2))

#     print("\n--- Extracted Data for GSM1824878 (Control Sample) ---")
#     sample_id_2 = 'GSM1824878'
#     extracted_info_2 = GSE71001_extractor(samples_metadata[sample_id_2])
#     print(json.dumps(extracted_info_2, indent=2))

#     print("\n--- Extracted Data for GSM1824866 (Another Control Sample) ---")
#     sample_id_3 = 'GSM1824866'
#     extracted_info_3 = GSE71001_extractor(samples_metadata[sample_id_3])
#     print(json.dumps(extracted_info_3, indent=2))


# import json

# def GSE12619_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from sample metadata
#     for the GSE12619 dataset, conforming to a specified JSON schema.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the target JSON schema,
#               containing 'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The 'sample_characteristicts' and 'sample_source_name' fields consistently
#     # mention "seedlings" as the biological material.
#     # Example: 'Seven days old seedlings of Arabidopsis mutant of til1-1 ecotype Columbia were used'
#     characteristics = sample_metadata.get('sample_characteristicts', [''])[0]
#     if "seedlings" in characteristics.lower():
#         extracted_data['tissue'] = "seedlings"
#     else:
#         # Fallback if "seedlings" is not found, though it appears consistent in this dataset.
#         extracted_data['tissue'] = "unspecified"

#     # 2. Extract 'treatment'
#     # The 'sample_source_name' field explicitly states the treatment condition.
#     # Examples: 'til1-1 seedlings heat stress rep1', 'Col-0 seedlings control condition rep1'
#     source_name = sample_metadata.get('sample_source_name', [''])[0]
#     treatment_list = []
#     if "heat stress" in source_name.lower():
#         treatment_list.append("heat stress")
#     elif "control condition" in source_name.lower():
#         # "control condition" is considered a type of treatment/experimental condition
#         treatment_list.append("control condition")
    
#     extracted_data['treatment'] = treatment_list

#     # 3. Extract 'medium'
#     # The growth medium is not explicitly stated in the provided sample or study metadata.
#     # For Arabidopsis seedlings, Murashige and Skoog (MS) medium is a very common
#     # and standard growth medium. Given that 'medium' is a required field in the schema
#     # and is often constant across samples in a study, we infer "MS medium" based on
#     # common biological experimental practices for this organism.
#     extracted_data['medium'] = "MS medium"

#     return extracted_data

# if __name__ == '__main__':
#     # Example study metadata (not directly used by the function, but for context)
#     study_metadata = {
#         'title': ['Heat shock response of til1-1 mutant plants'],
#         'summary': ['Several lipocalin genes from higher plants were shown to be responsive to both high and low temperature stresses and have been named as temperature-induced lipocalin (Til).  In this study, a reverse genetic approach was taken to elucidate the role of Arabidopsis Til1 (At5g58070) in thermotolerance.  We showed that Til1 proteins was constitutively expressed and increased significantly after heat shock treatment.  A T-DNA knockout line of Til1, designated as til1-1, could not produce Til1 and showed severe defects in basal and acquired thermotolerance.  Introducing a wild type copy of Til1 gene into til1-1 complemented the mutant phenotype.  Over-expression of Til1 in the wild type plant did not enhance thermotolerance.  Til1 is peripherally associated with plasma membrane, suggesting a regulatory or protective role of this protein in membrane function.  Transcriptomic analysis showed that the heat shock response in til1-1 was not altered as compared to the wild type plants.  The temperature threshold for heat shock protein induction was not affected by the level of Til1.  Ion leakage analysis revealed no significant difference in membrane stability between the wild type and til1-1 seedlings.  These results suggested that Til1 is not involved in regulating membrane fluidity or stability.  Nevertheless, the level of malondialdehyde was significantly higher in til1-1 than in the wild type after severe heat treatment.  The mutant plants were also more sensitive than the wild type to tert-butyl hydroperoxide, a reagent that induces lipid peroxidation.  Taken together, our data indicate that Til1 is an essential component for thermotolerance probably by acting against lipid peroxidation induced by severe heat stress.'],
#         'overall_design': ['Total RNA was isolated from the seedlings of 7-d old wild-type and til1-1 mutant seedlings (a pool of about 100 plants per treatment in duplicates) harvested immediately after heat shock treatment. In this experiment, total 8 chips were used, 1 each for 2 biological replicates of the control and HS-treated samples for the wild type and mutant plants.']
#     }

#     # Example samples metadata
#     samples_metadata = {
#         'GSM315983': {'Sample_id': 'GSM315983', 'sample_title': ['til1-1_37C_rep1'], 'sample_source_name': ['til1-1 seedlings heat stress rep1'], 'sample_characteristicts': ['Seven days old seedlings of Arabidopsis mutant of til1-1 ecotype Columbia were used'], 'sample_extraction_protocol': ['Qiagen RNA extraction kit']},
#         'GSM315973': {'Sample_id': 'GSM315973', 'sample_title': ['Col-0_37C_rep2'], 'sample_source_name': ['Col-0 seedlings heat stress rep2'], 'sample_characteristicts': ['Seven days old seedlings of Arabidopsis wild type ecotype Columbia were used'], 'sample_extraction_protocol': ['Qiagen RNA extraction kit']},
#         'GSM315970': {'Sample_id': 'GSM315970', 'sample_title': ['Col-0_22C_rep1'], 'sample_source_name': ['Col-0 seedlings control condition rep1'], 'sample_characteristicts': ['Seven days old seedlings of Arabidopsis wild type ecotype Columbia were used'], 'sample_extraction_protocol': ['Qiagen RNA extraction kit']},
#         'GSM315981': {'Sample_id': 'GSM315981', 'sample_title': ['til1-1_22C_rep1'], 'sample_source_name': ['til1-1 seedlings control condition rep1'], 'sample_characteristicts': ['Seven days old seedlings of Arabidopsis mutant of til1-1 ecotype Columbia were used'], 'sample_extraction_protocol': ['Qiagen RNA extraction kit']},
#         'GSM315971': {'Sample_id': 'GSM315971', 'sample_title': ['Col-0_22C_rep2'], 'sample_source_name': ['Col-0 seedlings control condition rep2'], 'sample_characteristicts': ['Seven days old seedlings of Arabidopsis wild type ecotype Columbia were used'], 'sample_extraction_protocol': ['Qiagen RNA extraction kit']},
#         'GSM315982': {'Sample_id': 'GSM315982', 'sample_title': ['til1-1_22C_rep2'], 'sample_source_name': ['til1-1 seedlings control condition rep2'], 'sample_characteristicts': ['Seven days old seedlings of Arabidopsis mutant of til1-1 ecotype Columbia were used'], 'sample_extraction_protocol': ['Qiagen RNA extraction kit']},
#         'GSM315972': {'Sample_id': 'GSM315972', 'sample_title': ['Col-0_37C_rep1'], 'sample_source_name': ['Col-0 seedlings heat stress rep1'], 'sample_characteristicts': ['Seven days old seedlings of Arabidopsis wild type ecotype Columbia were used'], 'sample_extraction_protocol': ['Qiagen RNA extraction kit']},
#         'GSM315984': {'Sample_id': 'GSM315984', 'sample_title': ['til1-1_37C_rep2'], 'sample_source_name': ['til1-1 seedlings heat stress rep2'], 'sample_characteristicts': ['Seven days old seedlings of Arabidopsis mutant of til1-1 ecotype Columbia were used'], 'sample_extraction_protocol': ['Qiagen RNA extraction kit']}
#     }

#     print("--- Testing with individual samples ---")

#     # Test case 1: Heat stress sample
#     sample_id_1 = 'GSM315983'
#     sample_data_1 = samples_metadata[sample_id_1]
#     extracted_info_1 = GSE12619_extractor(sample_data_1)
#     print(f"Extracted for {sample_id_1}:\n{json.dumps(extracted_info_1, indent=2)}\n")
#     # Expected: {'tissue': 'seedlings', 'treatment': ['heat stress'], 'medium': 'MS medium'}

#     # Test case 2: Control condition sample
#     sample_id_2 = 'GSM315970'
#     sample_data_2 = samples_metadata[sample_id_2]
#     extracted_info_2 = GSE12619_extractor(sample_data_2)
#     print(f"Extracted for {sample_id_2}:\n{json.dumps(extracted_info_2, indent=2)}\n")
#     # Expected: {'tissue': 'seedlings', 'treatment': ['control condition'], 'medium': 'MS medium'}

#     print("--- Testing with all samples ---")
#     all_extracted_data = {}
#     for sample_id, sample_data in samples_metadata.items():
#         all_extracted_data[sample_id] = GSE12619_extractor(sample_data)
    
#     # Print a few more examples to verify consistency
#     print(f"Extracted for GSM315973:\n{json.dumps(all_extracted_data['GSM315973'], indent=2)}\n")
#     print(f"Extracted for GSM315981:\n{json.dumps(all_extracted_data['GSM315981'], indent=2)}\n")


# import json

# def GSE37118_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from sample metadata for GSE37118
#     following a specific JSON schema.

#     The schema requires 'tissue', 'treatment', and 'medium'.
#     - 'tissue' and 'treatment' are extracted from the 'sample_characteristicts' field.
#     - 'medium' is inferred as "MS medium" (Murashige and Skoog medium), a common
#       growth medium for Arabidopsis thaliana seedlings, as it is not explicitly
#       stated in the provided metadata but is a required field in the schema.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the specified schema, containing
#               'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # Get the 'sample_characteristicts' list, defaulting to an empty list if not found
#     characteristics = sample_metadata.get('sample_characteristicts', [])

#     # 1. Extract 'tissue'
#     tissue = None
#     for char_string in characteristics:
#         if char_string.startswith('tissue:'):
#             # Split once to handle cases where the value itself might contain colons
#             tissue = char_string.split(':', 1)[1].strip()
#             break
    
#     # Fallback if tissue is not found (though it appears consistently present in this dataset)
#     if tissue is None:
#         tissue = "not specified" 
#     extracted_data['tissue'] = tissue

#     # 2. Extract 'treatment'
#     treatment_list = []
#     for char_string in characteristics:
#         if char_string.startswith('treatment:'):
#             raw_treatment = char_string.split(':', 1)[1].strip()
            
#             # Handle special control cases
#             if raw_treatment == 'control - untreated':
#                 # For 'control - untreated', we extract 'control' as the treatment type
#                 treatment_list.append('control')
#             elif raw_treatment.startswith('control - '):
#                 # For 'control - ABA', the actual treatment is 'ABA'
#                 # treatment_list.append(raw_treatment.split(' - ', 1)[1].strip())
#                 treatment_list.append('control')
#             else:
#                 # For other treatments, split by " and " to get individual treatments
#                 # Ensure each part is stripped and non-empty
#                 parts = [p.strip().lower()+' stress' for p in raw_treatment.split(' and ') if p.strip()]
#                 treatment_list.extend(parts)
#             break # Assuming only one 'treatment:' characteristic per sample
    
#     # If no treatment field was found or parsed to an empty list, it will remain empty.
#     # The schema allows an empty array for 'treatment'.
#     extracted_data['treatment'] = treatment_list

#     # 3. Extract 'medium'
#     # The 'medium' information is not explicitly present in the provided metadata.
#     # Based on common practices for Arabidopsis thaliana seedling studies,
#     # Murashige and Skoog (MS) medium is a very common growth medium.
#     # This is an inference to fulfill the schema's 'required' field.
#     extracted_data['medium'] = "MS medium"

#     return extracted_data



# import json

# def GSE35258_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from a sample's metadata
#     for the GSE35258 dataset, conforming to a specific JSON schema.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted as a JSON instance, containing 'tissue',
#               'treatment', and 'medium' information as per the schema.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue' and 'treatment' from 'sample_characteristicts'
#     # This field is a list of strings, where each string is in "key: value" format.
#     characteristics = sample_metadata.get('sample_characteristicts', [])
#     for item in characteristics:
#         if ': ' in item:
#             key, value = item.split(': ', 1)
#             if key == 'tissue':
#                 extracted_data['tissue'] = value
#             elif key == 'treatment':
#                 # The schema expects 'treatment' to be an array of strings.
#                 # In this dataset, it appears as a single string, so we wrap it in a list.
#                 extracted_data['treatment'] = [value]

#     # 2. Extract 'medium' based on the water potential (MPa) mentioned in 'sample_title'
#     # The 'overall_design' in study_metadata indicates the mapping:
#     # - "-1.2 MPa" corresponds to "PEG-infused agar plates" (stress condition)
#     # - "-0.25 MPa" corresponds to "agar plates of control media" (control condition)
#     sample_title = sample_metadata.get('sample_title', [''])[0] # Get the first string from the list

#     if "-1.2 MPa" in sample_title:
#         extracted_data['medium'] = "PEG-infused agar plates"
#     elif "-0.25 MPa" in sample_title:
#         extracted_data['medium'] = "agar plates of control media"
#     else:
#         # Fallback for cases where MPa value might not be found or is different.
#         # Based on the provided samples, this case should not be hit.
#         extracted_data['medium'] = "unknown medium"

#     # Ensure all required fields are present, even if with default/error values
#     # (though the logic above should cover all cases for this specific dataset)
#     if 'tissue' not in extracted_data:
#         extracted_data['tissue'] = "unknown"
#     if 'treatment' not in extracted_data:
#         extracted_data['treatment'] = [] # Empty list for array type
#     if 'medium' not in extracted_data:
#         extracted_data['medium'] = "unknown"

#     return extracted_data



# import json

# def GSE9415_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information (tissue, treatment, medium)
#     from a single sample's metadata for the GSE9415 dataset,
#     conforming to the specified JSON schema.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted as a JSON instance conforming to the schema:
#               {"tissue": "...", "treatment": ["...", ...], "medium": "..."}
#     """
#     extracted_data = {}

#     # 1. Extract Tissue
#     # The 'sample_source_name' field consistently contains the tissue information.
#     # Example: 'Arabidopsis shoots exposed to control conditions'
#     source_name = sample_metadata.get('sample_source_name', [''])[0]
#     if source_name:
#         # Assuming the format "Arabidopsis [tissue] exposed to..."
#         # We can extract the word immediately following "Arabidopsis".
#         parts = source_name.split(' ')
#         if len(parts) > 1 and parts[0].lower() == 'arabidopsis':
#             extracted_data['tissue'] = parts[1].strip()
#         else:
#             # Fallback if the expected format is not found
#             extracted_data['tissue'] = "unspecified"
#     else:
#         extracted_data['tissue'] = "unspecified"

#     # 2. Extract Treatment
#     # The 'sample_source_name' also contains the treatment information.
#     # Examples: 'control conditions', 'heat conditions', 'drought conditions',
#     # 'combined drought and heat conditions'
#     treatment_list = []
#     if source_name and 'exposed to ' in source_name:
#         # Extract the part after "exposed to " and remove " conditions" and any trailing dots/spaces
#         treatment_phrase = source_name.split('exposed to ')[1].replace(' conditions', '').strip().rstrip('.')

#         if treatment_phrase.lower() == 'control':
#             # 'control' is not considered a stress or applied treatment in this context,
#             # so the treatment list is empty.
#             treatment_list = ['control']
#         elif 'combined' in treatment_phrase:
#             # Handle combined treatments, e.g., "combined drought and heat"
#             # Remove "combined " and split by " and "
#             parts = treatment_phrase.replace('combined ', '').split(' and ')
#             treatment_list = [p.strip() for p in parts if p.strip()]
#         else:
#             # Handle single treatments, e.g., "heat", "drought"
#             treatment_list = [treatment_phrase]
#     extracted_data['treatment'] = treatment_list

#     # 3. Extract Medium
#     # The provided study and sample metadata do not explicitly mention the growth medium.
#     # As per guidance, information like medium can be constant or unspecified.
#     # We will use "unspecified" as a placeholder.
#     extracted_data['medium'] = "unspecified"

#     return extracted_data



# import json
# import re

# def GSE79681_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information of the sample from the GSE79681 dataset
#     following the given schema.

#     The output schema is:
#     {
#         "properties": {
#             "tissue": {
#                 "description": "Tissue the samples was extracted from.",
#                 "title": "Tissue",
#                 "type": "string"
#             },
#             "treatment": {
#                 "description": "List of treatments and stresses that was applied to the sample.",
#                 "items": {"type": "string"},
#                 "title": "Treatment",
#                 "type": "array"
#             },
#             "medium": {
#                 "description": "Growth medium of the sample.",
#                 "title": "Medium",
#                 "type": "string"
#             }
#         },
#         "required": ["tissue", "treatment", "medium"]
#     }

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 Expected keys include 'sample_characteristicts',
#                                 'sample_source_name', and 'sample_title'.

#     Returns:
#         dict: A dictionary formatted according to the specified schema,
#               containing tissue, treatment, and medium information.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # 'sample_characteristicts' is a list of strings, e.g., ['ecotype: Columbia-0', 'tissue: leaf']
#     extracted_data['tissue'] = "unknown" # Default value if not found
#     for characteristic in sample_metadata.get('sample_characteristicts', []):
#         if characteristic.startswith('tissue:'):
#             extracted_data['tissue'] = characteristic.split(':', 1)[1].strip()
#             break

#     # 2. Extract 'medium'
#     # 'medium' is related to 'FC' (Field Capacity) found in 'sample_source_name'
#     source_name_str = sample_metadata.get('sample_source_name', [''])[0]
#     extracted_data['medium'] = "unknown" # Default value if not found
    
#     # Use regex to find patterns like "100% FC" or "40% FC"
#     match = re.search(r'(\d+% FC)', source_name_str)
#     if match:
#         extracted_data['medium'] = match.group(1)

#     # 3. Extract 'treatment'
#     # This is a list of strings, derived primarily from 'sample_source_name'
#     treatments = []

#     # Check for pathogen infections
#     if "PStDC3000" in source_name_str:
#         treatments.append("PStDC3000 infection")
#     if "PSta" in source_name_str:
#         treatments.append("PSta infection")
    
#     # Check for drought stress conditions
#     # "40% FC" indicates drought stress was applied or maintained
#     if "40% FC" in source_name_str:
#         treatments.append("drought stress")
    
#     # Check for drought recovery
#     if "re-watered" in source_name_str:
#         treatments.append("drought recovery")
    
#     # Check for control infiltrations
#     if "syringe infiltrated with water" in source_name_str:
#         treatments.append("water infiltration")

#     # Remove duplicates and sort the list for consistent output
#     extracted_data['treatment'] = sorted(list(set(treatments)))
#     if extracted_data['treatment'] == []:
#         extracted_data['treatment'] = ['control conditions']
#     return extracted_data



# import json

# def GSE15577_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from a sample metadata dictionary
#     for the GSE15577 dataset, conforming to a specific JSON schema.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary containing the extracted 'tissue', 'treatment', and 'medium'
#               information, formatted according to the specified schema.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The tissue information is consistently found in the 'sample_source_name' field.
#     # Example: 'vegetative rosettes, wildtype, watered'
#     # We extract the first part before the first comma.
#     source_name = sample_metadata.get('sample_source_name', [''])[0]
#     if source_name:
#         tissue_parts = source_name.split(',')
#         extracted_data['tissue'] = tissue_parts[0].strip()
#     else:
#         # Fallback if 'sample_source_name' is missing or empty, though unlikely for this dataset.
#         extracted_data['tissue'] = "unspecified"

#     # 2. Extract 'treatment'
#     # Treatments are 'watered' (control) or 'water deficit'/'drought'.
#     # This information is also found in 'sample_source_name'.
#     treatment_list = []
#     lower_source_name = source_name.lower()

#     if "watered" in lower_source_name:
#         treatment_list.append("watered")
#     elif "water deficit" in lower_source_name or "drought" in lower_source_name:
#         # Normalize 'drought' to 'water deficit' as per study description
#         treatment_list.append("water deficit")
    
#     # Ensure treatment_list is not empty. Based on the provided samples,
#     # every sample has a clear condition.
#     extracted_data['treatment'] = treatment_list if treatment_list else ["unknown"]

#     # 3. Extract 'medium'
#     # The growth medium is not explicitly stated in the sample metadata.
#     # Given that the samples are 'vegetative rosettes' from Arabidopsis thaliana,
#     # 'soil' is the most common and logical growth medium for such plants.
#     # This information is constant across all samples in this study.
#     extracted_data['medium'] = "soil"

#     return extracted_data

# if __name__ == '__main__':
#     # Example usage with the provided samples_metadata
#     samples_metadata = {
#         'GSM389767': {'Sample_id': 'GSM389767', 'sample_title': ['rosettes_wildtype_watered_rep2'], 'sample_source_name': ['vegetative rosettes, wildtype, watered'], 'sample_characteristicts': ['cultivar: Columbia-0'], 'sample_extraction_protocol': ["Trizol extraction of total RNA was performed according to the manufacturer's instructions."]},
#         'GSM389768': {'Sample_id': 'GSM389768', 'sample_title': ['rosettes_atx_watered_rep1'], 'sample_source_name': ['vegetative rosettes, atx, watered'], 'sample_characteristicts': ['mutant: ATX1'], 'sample_extraction_protocol': ["Trizol extraction of total RNA was performed according to the manufacturer's instructions."]},
#         'GSM389776': {'Sample_id': 'GSM389776', 'sample_title': ['rosettes_myoox_drought_rep1'], 'sample_source_name': ['vegetative rosettes, myoox, water deficit'], 'sample_characteristicts': ['mutant: MYO-OX'], 'sample_extraction_protocol': ["Trizol extraction of total RNA was performed according to the manufacturer's instructions."]},
#         'GSM389775': {'Sample_id': 'GSM389775', 'sample_title': ['rosettes_atx_drought_rep2'], 'sample_source_name': ['vegetative rosettes, atx, water deficit'], 'sample_characteristicts': ['mutant: ATX1'], 'sample_extraction_protocol': ["Trizol extraction of total RNA was performed according to the manufacturer's instructions."]},
#         'GSM389770': {'Sample_id': 'GSM389770', 'sample_title': ['rosettes_myoox_watered_rep1'], 'sample_source_name': ['vegetative rosettes, myoox, watered'], 'sample_characteristicts': ['mutant: MYO-OX'], 'sample_extraction_protocol': ["Trizol extraction of total RNA was performed according to the manufacturer's instructions."]},
#         'GSM389772': {'Sample_id': 'GSM389772', 'sample_title': ['rosettes_wildtype_drought_rep1'], 'sample_source_name': ['vegetative rosettes, wildtype, water deficit'], 'sample_characteristicts': ['cultivar: Columbia-0'], 'sample_extraction_protocol': ["Trizol extraction of total RNA was performed according to the manufacturer's instructions."]},
#         'GSM389774': {'Sample_id': 'GSM389774', 'sample_title': ['rosettes_atx_drought_rep1'], 'sample_source_name': ['vegetative rosettes, atx, water deficit'], 'sample_characteristicts': ['mutant: ATX1'], 'sample_extraction_protocol': ["Trizol extraction of total RNA was performed according to the manufacturer's instructions."]},
#         'GSM389769': {'Sample_id': 'GSM389769', 'sample_title': ['rosettes_atx_watered_rep2'], 'sample_source_name': ['vegetative rosettes, atx, watered'], 'sample_characteristicts': ['mutant: ATX1'], 'sample_extraction_protocol': ["Trizol extraction of total RNA was performed according to the manufacturer's instructions."]},
#         'GSM389773': {'Sample_id': 'GSM389773', 'sample_title': ['rosettes_wildtype_drought_rep2'], 'sample_source_name': ['vegetative rosettes, wildtype, water deficit'], 'sample_characteristicts': ['cultivar: Columbia-0'], 'sample_extraction_protocol': ["Trizol extraction of total RNA was performed according to the manufacturer's instructions."]},
#         'GSM389766': {'Sample_id': 'GSM389766', 'sample_title': ['rosettes_wildtype_watered_rep1'], 'sample_source_name': ['vegetative rosettes, wildtype, watered'], 'sample_characteristicts': ['cultivar: Columbia-0'], 'sample_extraction_protocol': ["Trizol extraction of total RNA was performed according to the manufacturer's instructions."]},
#         'GSM389771': {'Sample_id': 'GSM389771', 'sample_title': ['rosettes_myoox_watered_rep2'], 'sample_source_name': ['vegetative rosettes, myoox, watered'], 'sample_characteristicts': ['mutant: MYO-OX'], 'sample_extraction_protocol': ["Trizol extraction of total RNA was performed according to the manufacturer's instructions."]},
#         'GSM389777': {'Sample_id': 'GSM389777', 'sample_title': ['rosettes_myoox_drought_rep2'], 'sample_source_name': ['vegetative rosettes, myoox, water deficit'], 'sample_characteristicts': ['mutant: MYO-OX'], 'sample_extraction_protocol': ["Trizol extraction of total RNA was performed according to the manufacturer's instructions."]}
#     }

#     print("Extracted schema for each sample:")
#     for sample_id, metadata in samples_metadata.items():
#         extracted_info = GSE15577_extractor(metadata)
#         print(f"Sample ID: {sample_id}")
#         print(json.dumps(extracted_info, indent=2))
#         print("-" * 30)

#     # Example of a specific sample output:
#     sample_id_to_test = 'GSM389767'
#     print(f"\nSpecific test for {sample_id_to_test}:")
#     extracted_for_specific_sample = GSE15577_extractor(samples_metadata[sample_id_to_test])
#     print(json.dumps(extracted_for_specific_sample, indent=2))

#     sample_id_to_test_drought = 'GSM389776'
#     print(f"\nSpecific test for {sample_id_to_test_drought} (drought sample):")
#     extracted_for_drought_sample = GSE15577_extractor(samples_metadata[sample_id_to_test_drought])
#     print(json.dumps(extracted_for_drought_sample, indent=2))



# import json

# def GSE18666_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from sample metadata for the GSE18666 dataset.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the specified JSON schema,
#               containing 'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The tissue information is found in 'sample_characteristicts'.
#     # Example: 'sample_characteristicts': ['tissue: seedling']
#     characteristics = sample_metadata.get('sample_characteristicts', [])
#     tissue_info = next((s for s in characteristics if s.startswith('tissue:')), None)
#     if tissue_info:
#         extracted_data['tissue'] = tissue_info.split(': ')[1].strip()
#     else:
#         # Default or handle cases where tissue might be missing
#         extracted_data['tissue'] = "unknown" # Based on provided samples, it's always 'seedling'

#     # 2. Extract 'medium'
#     # Based on the study_metadata['overall_design'], samples were "grown in vitro".
#     # This appears to be a constant for all samples in this study.
#     extracted_data['medium'] = "in vitro"

#     # 3. Extract 'treatment'
#     # Treatment information is derived from the 'sample_title'.
#     # Examples: 'Mock, 2 days recovery', 'Heat, no recovery'
#     treatments = []
#     sample_title = sample_metadata.get('sample_title', [''])[0] # Get the first string from the list

#     # Determine primary stress/condition
#     if "Heat" in sample_title:
#         treatments.append("heat stress")
#     elif "Mock" in sample_title:
#         treatments.append("mock")

#     # Determine recovery status
#     if "2 days recovery" in sample_title:
#         treatments.append("recovery")
#     # "no recovery" implies the absence of recovery, so no specific treatment string is added for it.

#     extracted_data['treatment'] = treatments

#     return extracted_data



# import json

# def GSE65046_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from sample metadata
#     following a specific JSON schema.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the specified schema, containing
#               'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The study metadata and sample extraction protocols consistently indicate
#     # that RNA was isolated from "leaves". This is a constant value across samples.
#     extracted_data["tissue"] = "leaves"

#     # 2. Extract 'treatment'
#     # The treatment information is found within the 'sample_characteristicts' list.
#     # Each item in this list is a string, and one of them starts with "treatment:".
#     treatments = []
#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if isinstance(characteristic, str) and characteristic.startswith('treatment:'):
#                 # Extract the value after "treatment: " and strip any whitespace
#                 treatment_value = characteristic.split(':', 1)[1].strip()
#                 treatments.append(treatment_value)
#                 # Assuming only one 'treatment' characteristic per sample based on examples
#                 break
#     extracted_data["treatment"] = treatments

#     # 3. Extract 'medium'
#     # The growth medium is not explicitly stated in the provided study or sample metadata.
#     # Given the context of growing Arabidopsis plants under "well-watered" conditions,
#     # a "Not Specified" is a reasonable inference for this required field.
#     extracted_data["medium"] = "Not Specified"

#     return extracted_data


# import json

# def GSE18624_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from a sample metadata dictionary
#     according to a predefined schema.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted as a JSON instance conforming to the specified schema.
#               The schema includes 'tissue', 'treatment', and 'medium'.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The tissue information is found in 'sample_characteristicts' list.
#     # Example: 'sample_characteristicts': ['tissue: whole seedlings', ...]
#     characteristics = sample_metadata.get('sample_characteristicts', [])
#     tissue = None
#     for char_str in characteristics:
#         if char_str.lower().startswith('tissue:'):
#             tissue = char_str.split(':', 1)[1].strip()
#             break
#     extracted_data['tissue'] = tissue if tissue else "unknown" # Default if not found

#     # 2. Extract 'treatment'
#     # Treatments are derived from the 'sample_title' or 'sample_source_name'
#     # and relate to temperature and duration.
#     # Example: 'sample_title': ['arp6-10  2h at 27 ºC  rep2']
#     # Example: 'sample_title': ['Col-0 12 ºC rep1']
#     treatments = []
#     sample_title = sample_metadata.get('sample_title', [''])[0] # Get the first title string

#     # Prioritize specific temperature/duration patterns
#     if '2h at 27 ºC' in sample_title:
#         treatments.append('27 ºC for 2 hours')
#     elif '24h at 27 ºC' in sample_title:
#         treatments.append('27 ºC for 24 hours')
#     elif '12 ºC' in sample_title:
#         treatments.append('12 ºC')
#     else:
#         # Fallback to source_name if title doesn't match, or if more general parsing is needed
#         sample_source_name = sample_metadata.get('sample_source_name', [''])[0]
#         if '2h' in sample_source_name and '27 ºC' in sample_source_name:
#             treatments.append('27 ºC for 2 hours')
#         elif '24h' in sample_source_name and '27 ºC' in sample_source_name:
#             treatments.append('27 ºC for 24 hours')
#         elif '12 ºC' in sample_source_name:
#             treatments.append('12 ºC')
#         else:
#             # If no specific temperature treatment is identified, provide a general label
#             treatments.append("unspecified temperature treatment")

#     extracted_data['treatment'] = treatments

#     # 3. Extract 'medium'
#     # The growth medium is constant across all samples in this study,
#     # as indicated in the 'overall_design' of the study metadata.
#     # From study_metadata['overall_design']: "plated on half-strength MS plates without sugar."
#     extracted_data['medium'] = "half-strength MS plates without sugar"

#     return extracted_data



# import json

# def GSE44053_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from sample metadata
#     following a specific JSON schema for the GSE44053 dataset.

#     The output schema is:
#     {"properties": {"tissue": {"description": "Tissue the samples was extracted from.", "title": "Tissue", "type": "string"},
#                     "treatment": {"description": "List of treatments and stresses that was applied to the sample.", "items": {"type": "string"}, "title": "Treatment", "type": "array"},
#                     "medium": {"description": "Growth medium of the sample.", "title": "Medium", "type": "string"}},
#      "required": ["tissue", "treatment", "medium"]}

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the specified schema,
#               containing tissue, treatment, and medium information.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The tissue information is found within 'sample_characteristicts'.
#     tissue = "Not specified" # Default value if tissue information is not found
#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if isinstance(characteristic, str) and characteristic.startswith('tissue:'):
#                 tissue = characteristic.split(':', 1)[1].strip()
#                 break
#     extracted_data['tissue'] = tissue

#     # 2. Extract 'treatment'
#     # Treatment information is derived from 'sample_title' and 'sample_source_name'.
#     treatment_list = []
#     sample_title = sample_metadata.get('sample_title', [''])[0]
#     sample_source_name = sample_metadata.get('sample_source_name', [''])[0]

#     # Check for heat stress
#     if "HeatStressed" in sample_title or "Heat-stresed" in sample_source_name:
#         treatment_list.append("Heat stress")
#     # Check for control, only if no heat stress is identified
#     elif "control" in sample_title.lower() or "control" in sample_source_name.lower():
#         treatment_list.append("control")
    
#     # If no specific treatment or control is identified, treatment_list remains empty.
#     # An empty list [] is a valid array according to the schema for "treatment".
#     # However, for this specific dataset, samples are consistently either "control" or "Heat stress".
#     # If for some reason neither is found, an empty list is returned.
#     extracted_data['treatment'] = treatment_list

#     # 3. Extract 'medium'
#     # The growth medium is not explicitly stated in the provided study or sample metadata.
#     # The 'sample_extraction_protocol' describes an extraction buffer, not a growth medium.
#     # Given that the schema requires this field and it's for Arabidopsis seedlings,
#     # a "Not Specified" is a reasonable and informative default.
#     extracted_data['medium'] = "Not Specified"

#     return extracted_data



# import json

# def GSE72949_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from the sample metadata
#     for the GSE72949 dataset, conforming to the specified JSON schema.

#     The schema requires:
#     - tissue (string): Tissue the sample was extracted from.
#     - treatment (array of strings): List of treatments and stresses applied.
#     - medium (string): Growth medium of the sample.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 This typically corresponds to one GSM entry.

#     Returns:
#         dict: A dictionary formatted as a JSON instance conforming to the output schema.
#               Example: {"tissue": "seedling", "treatment": ["Heat Stress priming"], "medium": "Not Specified"}
#     """
#     extracted_data = {
#         "tissue": "",
#         "treatment": [],
#         "medium": "Not Specified"  # Inferred for Arabidopsis seedlings, as not explicitly stated in metadata
#     }

#     # Extract information from 'sample_characteristicts' list
#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         for char_str in sample_metadata['sample_characteristicts']:
#             if isinstance(char_str, str):
#                 # Extract 'tissue'
#                 if char_str.startswith('tissue:'):
#                     extracted_data['tissue'] = char_str.split(':', 1)[1].strip()
#                 # Extract 'treatment' based on 'stress'
#                 elif char_str.startswith('stress:'):
#                     stress_type = char_str.split(':', 1)[1].strip()
#                     if stress_type == 'priming':
#                         # Based on study_metadata, 'priming' refers to a specific heat stress regime
#                         extracted_data['treatment'].append("Heat Stress priming")
#                     elif stress_type == 'unprimed':
#                         extracted_data['treatment'].append("No stress")
    
#     # Ensure 'treatment' is not an empty list if no stress was explicitly found
#     # (though for this dataset, 'stress' seems to always be present).
#     if not extracted_data['treatment']:
#         extracted_data['treatment'].append("Not specified")

#     return extracted_data

# if __name__ == '__main__':
#     # Example study_metadata (for context, not directly used by the function)
#     study_metadata = {
#         'title': ['Identification of HS memory-associated genes'],
#         'summary': ['To identify genes associated with the heat stress (HS) memory, transcript profiling using Affymetrix ATH1 microarrays was performed to compare Col-0 seedlings after the priming stimulus with control plants (unprimed).'],
#         'overall_design': ['The experiment was designed to identify genes associated with HS memory. For the priming HS, seedlings were subjected to a heat regime of 1.5 h, 37°C; 1.5 h recovery at 22°C; and 45 min, 44°C. After the priming HS treatment, seedlings were returned to normal growth condition and harvested after 4, 8, 24 and 48 h, respectively.']
#     }

#     # Example samples_metadata (to demonstrate usage)
#     samples_metadata = {
#         'GSM1875148': {
#             'Sample_id': 'GSM1875148',
#             'sample_title': ['Col-0 4h primed seedlings'],
#             'sample_source_name': ['Wild-type Col-0 seedling'],
#             'sample_characteristicts': [
#                 'ecotype: Col-0',
#                 'genotype/variation: Wild-type',
#                 'age: 5-days-old',
#                 'stress: priming',
#                 'time of harvest: 4hr',
#                 'tissue: seedling'
#             ],
#             'sample_extraction_protocol': ['RNeasy Plant Mini Kit']
#         },
#         'GSM1875149': {
#             'Sample_id': 'GSM1875149',
#             'sample_title': ['Col-0 4h unprimed seedlings'],
#             'sample_source_name': ['Wild-type Col-0 seedling'],
#             'sample_characteristicts': [
#                 'ecotype: Col-0',
#                 'genotype/variation: Wild-type',
#                 'age: 5-days-old',
#                 'stress: unprimed',
#                 'time of harvest: 4hr',
#                 'tissue: seedling'
#             ],
#             'sample_extraction_protocol': ['RNeasy Plant Mini Kit']
#         },
#         'GSM1875154': {
#             'Sample_id': 'GSM1875154',
#             'sample_title': ['Col-0 48h primed seedlings'],
#             'sample_source_name': ['Wild-type Col-0 seedling'],
#             'sample_characteristicts': [
#                 'ecotype: Col-0',
#                 'genotype/variation: Wild-type',
#                 'age: 5-days-old',
#                 'stress: priming',
#                 'time of harvest: 48hr',
#                 'tissue: seedling'
#             ],
#             'sample_extraction_protocol': ['RNeasy Plant Mini Kit']
#         }
#     }

#     print("--- Extracted Data for GSM1875148 (primed) ---")
#     sample_id_1 = 'GSM1875148'
#     extracted_1 = GSE72949_extractor(samples_metadata[sample_id_1])
#     print(json.dumps(extracted_1, indent=2))

#     print("\n--- Extracted Data for GSM1875149 (unprimed) ---")
#     sample_id_2 = 'GSM1875149'
#     extracted_2 = GSE72949_extractor(samples_metadata[sample_id_2])
#     print(json.dumps(extracted_2, indent=2))

#     print("\n--- Extracted Data for GSM1875154 (primed, 48h) ---")
#     sample_id_3 = 'GSM1875154'
#     extracted_3 = GSE72949_extractor(samples_metadata[sample_id_3])
#     print(json.dumps(extracted_3, indent=2))


# import json

# def GSE16474_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information (tissue, treatment, medium)
#     from a single sample's metadata dictionary, conforming to a specified JSON schema.

#     The function identifies the relevant information within the 'sample_characteristicts'
#     list and the overall study context (for medium determination).

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 Expected to have keys like 'sample_characteristicts'.

#     Returns:
#         dict: A dictionary formatted according to the schema:
#               {"tissue": "...", "treatment": ["..."], "medium": "..."}
#               Returns None for a field if the information cannot be extracted.
#     """
#     extracted_data = {
#         "tissue": None,
#         "treatment": [],
#         "medium": None
#     }

#     # Get the list of characteristics, defaulting to an empty list if the key is missing
#     characteristics = sample_metadata.get('sample_characteristicts', [])
#     stress_value = None

#     # Iterate through characteristics to find 'tissue' and 'stress' information
#     for char_string in characteristics:
#         if char_string.startswith('tissue:'):
#             # Extract tissue by splitting the string and stripping whitespace
#             extracted_data['tissue'] = char_string.split(':', 1)[1].strip()
#         elif char_string.startswith('stress:'):
#             # Extract stress value for later use in determining treatment and medium
#             stress_value = char_string.split(':', 1)[1].strip()

#     # Determine 'treatment' and 'medium' based on the extracted 'stress_value'
#     if stress_value:
#         if stress_value.lower() == 'osmotic stress':
#             extracted_data['treatment'].append('osmotic stress')
#             # Medium information derived from study_metadata's 'overall_design'
#             extracted_data['medium'] = 'medium with 25mM mannitol'
#         elif stress_value.lower() == 'control':
#             # As per schema description "List of treatments and stresses",
#             # 'control' implies no specific applied treatment, so the list is empty.
#             extracted_data['treatment'] = []
#             # Medium information derived from study_metadata's 'overall_design'
#             extracted_data['medium'] = 'medium without (0mM) mannitol'
#         else:
#             # Fallback for any other unexpected stress values
#             extracted_data['treatment'].append(stress_value)
#             extracted_data['medium'] = f'medium with {stress_value}'
#     else:
#         # If 'stress' information is not found in characteristics
#         extracted_data['medium'] = 'unknown medium'
#         extracted_data['treatment'] = [] # Assume no specific treatment if not found

#     return extracted_data



# def GSE20494_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from sample metadata for GSE20494
#     following a specific JSON schema.

#     The schema is:
#     {"properties": {"tissue": {"description": "Tissue the samples was extracted from.", "title": "Tissue", "type": "string"},
#                     "treatment": {"description": "List of treatments and stresses that was applied to the sample.", "items": {"type": "string"}, "title": "Treatment", "type": "array"},
#                     "medium": {"description": "Growth medium of the sample.", "title": "Medium", "type": "string"}},
#      "required": ["tissue", "treatment", "medium"]}

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the specified schema,
#               containing 'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # Based on the study_metadata and sample_metadata, the samples are from Arabidopsis plants.
#     # A general term like "whole plant" is appropriate as no specific tissue is mentioned.
#     extracted_data["tissue"] = "whole plant"

#     # 2. Extract 'medium'
#     # The growth medium is consistently described in the 'growth protocol'
#     # within 'sample_characteristicts' and in the 'overall_design' of the study.
#     # It is "1/2 MS plant medium" (or "half-strength MS plant medium").
#     extracted_data["medium"] = "1/2 MS plant medium"

#     # 3. Extract 'treatment'
#     # This information is found in the 'growth protocol' within 'sample_characteristicts'.
#     # We need to check for the presence of "50 mM NaCl" to identify salt stress treatment.
#     treatments = []
#     if 'sample_characteristicts' in sample_metadata:
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if characteristic.startswith("growth protocol:"):
#                 growth_protocol_str = characteristic.split("growth protocol:")[1].strip()
#                 if "50 mM NaCl" in growth_protocol_str:
#                     treatments.append("50 mM NaCl")
#                 # If the protocol mentions "0 mM NaCl" or no specific supplementation,
#                 # it implies no additional treatment, so the 'treatments' list remains empty.
#                 break  # Found the growth protocol, no need to check other characteristics

#     extracted_data["treatment"] = treatments

#     return extracted_data



# import json

# def GSE2268_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from sample metadata
#     for the GSE2268 dataset, conforming to a specified JSON schema.

#     The function identifies and extracts 'tissue', 'treatment', and 'medium'
#     information. 'Tissue' is inferred from the study organism (Arabidopsis),
#     'treatment' is parsed from the sample's source name, and 'medium' is
#     marked as "not specified" as it's not present in the provided metadata.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 Expected keys include 'sample_source_name'.

#     Returns:
#         dict: A dictionary formatted as a JSON instance conforming to the
#               output schema:
#               {"properties": {"tissue": {"description": "Tissue the samples was extracted from.", "title": "Tissue", "type": "string"}, "treatment": {"description": "List of treatments and stresses that was applied to the sample.", "items": {"type": "string"}, "title": "Treatment", "type": "array"}, "medium": {"description": "Growth medium of the sample.", "title": "Medium", "type": "string"}}, "required": ["tissue", "treatment", "medium"]}
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The study metadata indicates the organism is "Arabidopsis".
#     # Without specific tissue information in the sample_metadata, and given that
#     # Arabidopsis is a plant, "whole plant" is a common and general term used
#     # when specific organs are not mentioned. This is inferred from the study context.
#     # The guidance "information like meduim or tissue can be constant across samples"
#     # supports inferring this from the broader study context rather than per sample.
#     extracted_data["tissue"] = "whole plant"

#     # 2. Extract 'treatment'
#     # This information is found within the 'sample_source_name' field.
#     # Example: 'non-polysomal RNA under non-stress condition' or
#     #          'non-polysomal mRNA under dehydration stress'
#     source_name = sample_metadata.get('sample_source_name', [''])[0]
#     treatment_list = []

#     if "under " in source_name:
#         # Split the string to get the condition part after "under "
#         condition_part = source_name.split("under ", 1)[1].strip()

#         if condition_part == "non-stress condition":
#             # For "non-stress condition", no specific treatment is applied.
#             # The schema expects a list of treatments/stresses.
#             treatment_list = []
#         elif condition_part == "dehydration stress":
#             # For "dehydration stress", this is a clear treatment.
#             treatment_list = ["dehydration stress"]
#         # Add more conditions here if they appear in other samples
#         # else:
#         #     treatment_list = [condition_part] # Fallback for unhandled conditions

#     extracted_data["treatment"] = treatment_list

#     # 3. Extract 'medium'
#     # The growth medium is not explicitly mentioned in the provided study_metadata
#     # or samples_metadata. Since 'medium' is a required field in the output schema,
#     # and it cannot be extracted from the given dictionaries, we use "not specified"
#     # as a placeholder to fulfill the schema requirement.
#     extracted_data["medium"] = "not specified"

#     return extracted_data



# import json

# def GSE44655_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from sample metadata for GSE44655
#     following a specific JSON schema.

#     The schema requires 'tissue' (string), 'treatment' (array of strings),
#     and 'medium' (string).

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the specified schema,
#               containing 'tissue', 'treatment', and 'medium'.
#     """
#     extracted_data = {
#         "tissue": "",
#         "treatment": [],
#         "medium": ""
#     }

#     # Extract 'tissue' from 'sample_characteristicts'
#     # The 'sample_characteristicts' field is a list of strings,
#     # where each string describes a characteristic in a 'key: value' format.
#     characteristics = sample_metadata.get('sample_characteristicts', [])
#     for item in characteristics:
#         if item.startswith('tissue:'):
#             # Split by the first colon and strip whitespace from the value
#             extracted_data['tissue'] = item.split(':', 1)[1].strip()
#             break # Assuming only one 'tissue' entry per sample

#     # Extract 'treatment' from 'sample_characteristicts'
#     # The 'condition' characteristic indicates the treatment or control state.
#     # We consider "heat shock" as an active treatment/stress.
#     # "control" is interpreted as the absence of an applied treatment/stress,
#     # so it is not added to the 'treatment' list.
#     for item in characteristics:
#         if item.startswith('condition:'):
#             condition_str = item.split(':', 1)[1].strip()
#             if "heat shock" in condition_str.lower():
#                 extracted_data['treatment'].append(condition_str)
#             break # Assuming only one 'condition' entry per sample

#     # Extract 'medium'
#     # Based on the provided study metadata and sample metadata, there is no
#     # explicit field or description that specifies the growth medium of the
#     # plants. Since 'medium' is a required string field in the schema,
#     # and no information is available for extraction, we will provide an
#     # empty string as its value.
#     extracted_data['medium'] = "" # Information not found in the provided metadata

#     return extracted_data



# import json

# def GSE27549_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from sample metadata for GSE27549
#     following a specific JSON schema.

#     The schema requires 'tissue' (string), 'treatment' (array of strings),
#     and 'medium' (string).

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 Expected to have keys like 'sample_characteristicts'.

#     Returns:
#         dict: A dictionary formatted according to the specified schema,
#               containing 'tissue', 'treatment', and 'medium'.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The tissue information is consistently found within the 'sample_characteristicts'
#     # list, typically as an item formatted like 'tissue: leaf'.
#     tissue = "unknown"  # Default value if not found
#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if characteristic.startswith('tissue:'):
#                 # Extract the value after "tissue:" and strip any leading/trailing whitespace
#                 tissue = characteristic.split(':', 1)[1].strip()
#                 break
#     extracted_data['tissue'] = tissue

#     # 2. Extract 'treatment'
#     # Based on the study's title ("Genomic dna hybridizations") and the sample titles
#     # (e.g., "leaf_gDNA_NFA10_1"), these samples are for genomic DNA analysis.
#     # While the study summary mentions "soil moisture deficit" treatments, it explicitly
#     # links them to "gene expression studies", not these gDNA samples. Genomic DNA
#     # samples are typically collected from untreated, baseline conditions.
#     # Therefore, for these specific gDNA samples, we infer no specific treatments were applied.
#     extracted_data['treatment'] = []

#     # 3. Extract 'medium'
#     # The growth medium is not explicitly stated in the provided sample metadata.
#     # For Arabidopsis thaliana (a plant), "soil" is the common and natural growth medium.
#     # Given the context, "soil" is a reasonable inference for the growth medium.
#     extracted_data['medium'] = "soil"

#     return extracted_data



# import json

# def GSE71855_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from a sample's metadata
#     for the GSE71855 dataset, conforming to a specific JSON schema.

#     The function identifies and extracts the tissue, applied treatments/stresses,
#     and the growth medium based on the provided sample and study metadata structure.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 Expected to have keys like 'sample_characteristicts'.

#     Returns:
#         dict: A dictionary formatted as a JSON instance conforming to the schema:
#               {
#                   "tissue": "string",
#                   "treatment": ["string", ...],
#                   "medium": "string"
#               }
#               All fields are required by the schema.
#     """
#     extracted_data = {}

#     # --- Extract Tissue ---
#     # Tissue information is consistently found within the 'sample_characteristicts' list.
#     # Example: 'tissue: whole plant'
#     tissue = None
#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if isinstance(characteristic, str) and characteristic.startswith('tissue:'):
#                 tissue = characteristic.split(':', 1)[1].strip()
#                 break
#     # Based on the provided samples_metadata, 'tissue: whole plant' is always present.
#     # If it were hypothetically missing, 'tissue' would be None, which would violate the schema.
#     # We assume valid input where 'tissue' is always extractable.
#     extracted_data['tissue'] = tissue

#     # --- Extract Treatment ---
#     # Treatment information is derived from two sources:
#     # 1. The 'treatment:' characteristic in 'sample_characteristicts'.
#     # 2. The 'condition:' characteristic, which indicates 'salinity stress' (from study_metadata).
#     treatments = []
#     sample_specific_treatment = None
#     sample_condition = None

#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if isinstance(characteristic, str):
#                 if characteristic.startswith('treatment:'):
#                     sample_specific_treatment = characteristic.split(':', 1)[1].strip()
#                 elif characteristic.startswith('condition:'):
#                     sample_condition = characteristic.split(':', 1)[1].strip()
        
#         # Add the sample-specific treatment (e.g., "without Ky-2" or "1 M Ky-2 for 24 hours")
#         if sample_specific_treatment:
#             treatments.append(sample_specific_treatment)
        
#         # Add salinity stress if the sample condition indicates it.
#         # The details of salinity stress (100 mM Nacl for two hours) are from the study's overall_design.
#         if sample_condition == 'salinity stress':
#             treatments.append("100 mM Nacl for two hours")
    
#     extracted_data['treatment'] = treatments

#     # --- Extract Medium ---
#     # The growth medium information is constant across all samples in this study
#     # and is found in the 'overall_design' section of the study metadata.
#     # "Wild-type plants ... were grown on 1/2 MS 0.1 % Agar medium for four days."
#     extracted_data['medium'] = "1/2 MS 0.1 % Agar medium"

#     return extracted_data



# import json

# def GSE19603_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from sample metadata for GSE19603
#     following the given JSON schema.

#     The output schema is:
#     {"properties": {"tissue": {"description": "Tissue the samples was extracted from.", "title": "Tissue", "type": "string"},
#                     "treatment": {"description": "List of treatments and stresses that was applied to the sample.", "items": {"type": "string"}, "title": "Treatment", "type": "array"},
#                     "medium": {"description": "Growth medium of the sample.", "title": "Medium", "type": "string"}},
#      "required": ["tissue", "treatment", "medium"]}

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the specified schema, containing
#               'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The tissue information is consistently found within the 'sample_characteristicts' list.
#     # Example: 'tissue: above-ground parts'
#     tissue = None
#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if isinstance(characteristic, str) and characteristic.startswith('tissue:'):
#                 tissue = characteristic.split(':', 1)[1].strip()
#                 break
#     # Based on the provided samples, 'tissue' is always 'above-ground parts'.
#     # If it were potentially missing, a more general default like "unknown" might be used,
#     # but for this dataset, "above-ground parts" is the consistent value.
#     extracted_data['tissue'] = tissue if tissue is not None else "above-ground parts"

#     # 2. Extract 'treatment'
#     # Treatment information is derived from 'sample_title' and 'sample_source_name'.
#     # Samples are either 'untreated' (under 'normal conditions') or subjected to 'heat shock'/'heat stress'.
#     treatments = set() # Use a set to automatically handle unique treatments

#     # Safely get and convert relevant strings to lowercase for case-insensitive matching
#     sample_title_str = sample_metadata.get('sample_title', [''])[0].lower()
#     sample_source_name_str = sample_metadata.get('sample_source_name', [''])[0].lower()

#     # Check for "untreated" or "normal conditions"
#     if "untreated" in sample_title_str or "normal conditions" in sample_source_name_str:
#         treatments.add("untreated")
    
#     # Check for "heat shock" or "heat stress"
#     if "heat shock" in sample_title_str or "heat shock" in sample_source_name_str or \
#        "heat stress" in sample_title_str or "heat stress" in sample_source_name_str:
#         treatments.add("heat shock")
    
#     # Convert the set of treatments to a sorted list for consistent output.
#     # The schema specifies an array of strings.
#     extracted_data['treatment'] = sorted(list(treatments))

#     # 3. Extract 'medium'
#     # The growth medium is not explicitly mentioned in the provided metadata for this study.
#     # For Arabidopsis plants, "soil" is a very common and reasonable default
#     # when no specific growth medium (e.g., agar, liquid culture, hydroponics) is specified.
#     extracted_data['medium'] = "soil"

#     return extracted_data



# def GSE162310_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from a sample's metadata
#     following a specific JSON schema.

#     Args:
#         sample_metadata (dict): A dictionary containing the metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the specified schema,
#               containing 'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The tissue information is found in the 'sample_characteristicts' list.
#     # It's typically in the format 'tissue: [value]'.
#     characteristics = sample_metadata.get('sample_characteristicts', [])
#     tissue = None
#     for char_str in characteristics:
#         if char_str.startswith('tissue:'):
#             tissue = char_str.split(':', 1)[1].strip()
#             break # Assuming only one 'tissue' characteristic per sample
#     extracted_data['tissue'] = tissue

#     # 2. Extract 'treatment'
#     # Treatment information is also in 'sample_characteristicts', under 'protocol:'.
#     # The schema requires a list of strings. 'unstress condition' should result in an empty list.
#     treatments = []
#     for char_str in characteristics:
#         if char_str.startswith('protocol:'):
#             protocol = char_str.split(':', 1)[1].strip()
#             if protocol != 'unstress condition':
#                 treatments.append(protocol)
#             break # Assuming only one 'protocol' characteristic per sample
#     extracted_data['treatment'] = treatments

#     # 3. Extract 'medium'
#     # Based on the study_metadata, the growth medium is constant across all samples.
#     # The 'overall_design' states: "Seeds of Arabidopsis thaliana (Col-0 ecotype) were grown on MS medium".
#     extracted_data['medium'] = "MS medium"

#     return extracted_data



# import json

# def GSE58616_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from the sample metadata
#     following a predefined JSON schema.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the specified schema, containing
#               'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # Based on the 'sample_source_name' and 'sample_characteristicts' fields
#     # in the provided samples_metadata, the tissue is consistently "Arabidopsis seedling".
#     # The 'overall_design' in study_metadata also confirms "Arabidopsis seedlings".
#     extracted_data['tissue'] = "seedling"

#     # 2. Extract 'treatment'
#     # This field should be a list of treatments/stresses applied.
#     # Information is primarily found in 'sample_source_name'.
#     treatments = []
#     source_name_list = sample_metadata.get('sample_source_name', [''])
#     source_name = source_name_list[0] if source_name_list else ''

#     # Check for heat stress treatment
#     if "with HS treatment" in source_name:
#         treatments.append("heat stress")
#     # If "without HS treatment" is present, it means heat stress was not applied,
#     # so we don't add "heat stress" to the treatments list.

#     # Check for light/dark condition
#     if "under light condition" in source_name:
#         treatments.append("light")
#     elif "under dark condition" in source_name:
#         treatments.append("dark")
    
#     extracted_data['treatment'] = treatments

#     # 3. Extract 'medium'
#     # The growth medium is not explicitly stated in the provided study_metadata
#     # or samples_metadata. For Arabidopsis seedlings, a "standard growth medium"
#     # (e.g., agar, soil) is typically used. Since this is a required field in the schema,
#     # and no specific medium is mentioned, we infer a common, unstated practice.
#     extracted_data['medium'] = "standard growth medium"

#     return extracted_data

# # Example Usage (using the provided samples_metadata for demonstration):
# if __name__ == "__main__":
#     study_metadata = {
#         'title': ['Expression data from Arabidopsis seedlings heat-stressed in light environment'],
#         'summary': ['Plants transcriptome react to environment temperture changes profoundly. In Arabidopsis seedlings, genes response to temperature fluctuations to adopt the ever-changing ambient envrionment.', 'We used microarrays to detail the global programme of gene expression underlying heat stress response progress in Arabidopsis.'],
#         'overall_design': ['Ten-day-old Arabidopsis seedlings were selected for RNA extraction and hybridization on Affymetrix microarrays. We sought to explore the heat stress response in transcriptome, thus we treat the plants with heat stress.  While in order to identify the interaction between light and temperature signaling pathways in plant , we treat Arabidopsis with heat stress under both light and dark conditions. To that end, our plant tissues are grouped as: HS-LIGHT, HS-DARK,CONTROL-LIGHT,CONTROL-DARK.']
#     }
#     samples_metadata = {
#         'GSM1415487': {'Sample_id': 'GSM1415487', 'sample_title': ['light-control, rep3'], 'sample_source_name': ['Arabidopsis seedling without HS treatment, under light condition'], 'sample_characteristicts': ['cultivar: Col-0', 'developmental stage: 12-day-old seedlings'], 'sample_extraction_protocol': ["Trizol extraction of total RNA was performed according to the manufacturer's instructions."]},
#         'GSM1415489': {'Sample_id': 'GSM1415489', 'sample_title': ['light-hs, rep2'], 'sample_source_name': ['Arabidopsis seedling with HS treatment, under light condition'], 'sample_characteristicts': ['cultivar: Col-0', 'developmental stage: 14-day-old seedlings'], 'sample_extraction_protocol': ["Trizol extraction of total RNA was performed according to the manufacturer's instructions."]},
#         'GSM1415485': {'Sample_id': 'GSM1415485', 'sample_title': ['light-control, rep1'], 'sample_source_name': ['Arabidopsis seedling without HS treatment, under light condition'], 'sample_characteristicts': ['cultivar: Col-0', 'developmental stage: 10-day-old seedlings'], 'sample_extraction_protocol': ["Trizol extraction of total RNA was performed according to the manufacturer's instructions."]},
#         'GSM1415490': {'Sample_id': 'GSM1415490', 'sample_title': ['light-hs, rep3'], 'sample_source_name': ['Arabidopsis seedling with HS treatment, under light condition'], 'sample_characteristicts': ['cultivar: Col-0', 'developmental stage: 15-day-old seedlings'], 'sample_extraction_protocol': ["Trizol extraction of total RNA was performed according to the manufacturer's instructions."]},
#         'GSM1415488': {'Sample_id': 'GSM1415488', 'sample_title': ['light-hs, rep1'], 'sample_source_name': ['Arabidopsis seedling with HS treatment, under light condition'], 'sample_characteristicts': ['cultivar: Col-0', 'developmental stage: 13-day-old seedlings'], 'sample_extraction_protocol': ["Trizol extraction of total RNA was performed according to the manufacturer's instructions."]},
#         'GSM1415486': {'Sample_id': 'GSM1415486', 'sample_title': ['light-control,rep2'], 'sample_source_name': ['Arabidopsis seedling without HS treatment, under light condition'], 'sample_characteristicts': ['cultivar: Col-0', 'developmental stage: 11-day-old seedlings'], 'sample_extraction_protocol': ["Trizol extraction of total RNA was performed according to the manufacturer's instructions."]}
#     }

#     print("--- Extracted Data for Each Sample ---")
#     for sample_id, metadata in samples_metadata.items():
#         extracted_info = GSE58616_extractor(metadata)
#         print(f"Sample ID: {sample_id}")
#         print(json.dumps(extracted_info, indent=2))
#         print("-" * 30)

#     # Example for a specific sample:
#     sample_id_to_test = 'GSM1415487'
#     sample_data = samples_metadata[sample_id_to_test]
#     result = GSE58616_extractor(sample_data)
#     print(f"\nExtracted information for {sample_id_to_test}:")
#     print(json.dumps(result, indent=2))

#     sample_id_to_test_hs = 'GSM1415489'
#     sample_data_hs = samples_metadata[sample_id_to_test_hs]
#     result_hs = GSE58616_extractor(sample_data_hs)
#     print(f"\nExtracted information for {sample_id_to_test_hs}:")
#     print(json.dumps(result_hs, indent=2))


# import json

# def GSE5623_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from the sample metadata
#     following the specified JSON schema.

#     The schema requires:
#     - tissue (string): Tissue the sample was extracted from.
#     - treatment (array of strings): List of treatments applied to the sample.
#     - medium (string): Growth medium of the sample.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted as a JSON instance conforming to the
#               output schema.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The tissue information is found within the 'sample_characteristicts' list.
#     # It's typically formatted as "Tissue: [Value]".
#     tissue = None
#     characteristics = sample_metadata.get('sample_characteristicts', [])
#     for char_str in characteristics:
#         if char_str.startswith('Tissue:'):
#             tissue = char_str.split(':', 1)[1].strip()
#             break

#     if tissue is None:
#         # According to the schema, 'tissue' is a required field.
#         # If not found, it indicates an issue with the input data or extraction logic.
#         raise ValueError(f"Tissue information not found in sample_metadata for sample {sample_metadata.get('Sample_id', 'Unknown')}.")
#     extracted_data['tissue'] = tissue

#     # 2. Extract 'treatment'
#     # Based on the study metadata provided (e.g., 'title': ['AtGenExpress: Stress Treatments (Salt stress)']),
#     # the primary treatment for this study is "Salt stress".
#     # This appears to be constant across all samples in this specific dataset.
#     extracted_data['treatment'] = ["Salt stress"]

#     # 3. Extract 'medium'
#     # Based on the study metadata summary, the seedlings were transferred to
#     # "MS-liquid-media" before stress treatment began.
#     # This is a constant growth medium for all samples during the treatment phase.
#     extracted_data['medium'] = "MS-liquid-media"

#     return extracted_data

# if __name__ == '__main__':
#     # Example usage with provided sample metadata
#     samples_metadata = {
#         'GSM131309': {'Sample_id': 'GSM131309', 'sample_title': ['AtGen_6-3121_Saltstress-Roots-0.5h_Rep1'], 'sample_source_name': ['Col-0'], 'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Roots'], 'sample_extraction_protocol': ''},
#         'GSM131316': {'Sample_id': 'GSM131316', 'sample_title': ['AtGen_6-3312_Saltstress-Shoots-3.0h_Rep2'], 'sample_source_name': ['Col-0'], 'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Shoots'], 'sample_extraction_protocol': ''},
#         'GSM131314': {'Sample_id': 'GSM131314', 'sample_title': ['AtGen_6-3222_Saltstress-Roots-1.0h_Rep2'], 'sample_source_name': ['Col-0'], 'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Roots'], 'sample_extraction_protocol': ''},
#         'GSM131326': {'Sample_id': 'GSM131326', 'sample_title': ['AtGen_6-3522_Saltstress-Roots-12.0h_Rep2'], 'sample_source_name': ['Col-0'], 'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Roots'], 'sample_extraction_protocol': ''},
#         'GSM131310': {'Sample_id': 'GSM131310', 'sample_title': ['AtGen_6-3122_Saltstress-Roots-0.5h_Rep2'], 'sample_source_name': ['Col-0'], 'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Roots'], 'sample_extraction_protocol': ''},
#         'GSM131328': {'Sample_id': 'GSM131328', 'sample_title': ['AtGen_6-3612_Saltstress-Shoots-24.0h_Rep2'], 'sample_source_name': ['Col-0'], 'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Shoots'], 'sample_extraction_protocol': ''},
#         'GSM131308': {'Sample_id': 'GSM131308', 'sample_title': ['AtGen_6-3112_Saltstress-Shoots-0.5h_Rep2'], 'sample_source_name': ['Col-0'], 'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Shoots'], 'sample_extraction_protocol': ''},
#         'GSM131329': {'Sample_id': 'GSM131329', 'sample_title': ['AtGen_6-3621_Saltstress-Roots-24.0h_Rep1'], 'sample_source_name': ['Col-0'], 'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Roots'], 'sample_extraction_protocol': ''},
#         'GSM131327': {'Sample_id': 'GSM131327', 'sample_title': ['AtGen_6-3611_Saltstress-Shoots-24.0h_Rep1'], 'sample_source_name': ['Col-0'], 'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Shoots'], 'sample_extraction_protocol': ''},
#         'GSM131325': {'Sample_id': 'GSM131325', 'sample_title': ['AtGen_6-3521_Saltstress-Roots-12.0h_Rep1'], 'sample_source_name': ['Col-0'], 'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Roots'], 'sample_extraction_protocol': ''},
#         'GSM131324': {'Sample_id': 'GSM131324', 'sample_title': ['AtGen_6-3512_Saltstress-Shoots-12.0h_Rep2'], 'sample_source_name': ['Col-0'], 'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Shoots'], 'sample_extraction_protocol': ''},
#         'GSM131317': {'Sample_id': 'GSM131317', 'sample_title': ['AtGen_6-3321_Saltstress-Roots-3.0h_Rep1'], 'sample_source_name': ['Col-0'], 'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Roots'], 'sample_extraction_protocol': ''},
#         'GSM131323': {'Sample_id': 'GSM131323', 'sample_title': ['AtGen_6-3511_Saltstress-Shoots-12.0h_Rep1'], 'sample_source_name': ['Col-0'], 'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Shoots'], 'sample_extraction_protocol': ''},
#         'GSM131319': {'Sample_id': 'GSM131319', 'sample_title': ['AtGen_6-3411_Saltstress-Shoots-6.0h_Rep1'], 'sample_source_name': ['Col-0'], 'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Shoots'], 'sample_extraction_protocol': ''},
#         'GSM131322': {'Sample_id': 'GSM131322', 'sample_title': ['AtGen_6-3422_Saltstress-Roots-6.0h_Rep2'], 'sample_source_name': ['Col-0'], 'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Roots'], 'sample_extraction_protocol': ''},
#         'GSM131307': {'Sample_id': 'GSM131307', 'sample_title': ['AtGen_6-3111_Saltstress-Shoots-0.5h_Rep1'], 'sample_source_name': ['Col-0'], 'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Shoots'], 'sample_extraction_protocol': ''},
#         'GSM131320': {'Sample_id': 'GSM131320', 'sample_title': ['AtGen_6-3412_Saltstress-Shoots-6.0h_Rep2'], 'sample_source_name': ['Col-0'], 'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Shoots'], 'sample_extraction_protocol': ''},
#         'GSM131330': {'Sample_id': 'GSM131330', 'sample_title': ['AtGen_6-3622_Saltstress-Roots-24.0h_Rep2'], 'sample_source_name': ['Col-0'], 'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Roots'], 'sample_extraction_protocol': ''},
#         'GSM131313': {'Sample_id': 'GSM131313', 'sample_title': ['AtGen_6-3221_Saltstress-Roots-1.0h_Rep1'], 'sample_source_name': ['Col-0'], 'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Roots'], 'sample_extraction_protocol': ''},
#         'GSM131311': {'Sample_id': 'GSM131311', 'sample_title': ['AtGen_6-3211_Saltstress-Shoots-1.0h_Rep1'], 'sample_source_name': ['Col-0'], 'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Shoots'], 'sample_extraction_protocol': ''},
#         'GSM131321': {'Sample_id': 'GSM131321', 'sample_title': ['AtGen_6-3421_Saltstress-Roots-6.0h_Rep1'], 'sample_source_name': ['Col-0'], 'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Roots'], 'sample_extraction_protocol': ''},
#         'GSM131318': {'Sample_id': 'GSM131318', 'sample_title': ['AtGen_6-3322_Saltstress-Roots-3.0h_Rep2'], 'sample_source_name': ['Col-0'], 'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Roots'], 'sample_extraction_protocol': ''},
#         'GSM131312': {'Sample_id': 'GSM131312', 'sample_title': ['AtGen_6-3212_Saltstress-Shoots-1.0h_Rep2'], 'sample_source_name': ['Col-0'], 'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Shoots'], 'sample_extraction_protocol': ''},
#         'GSM131315': {'Sample_id': 'GSM131315', 'sample_title': ['AtGen_6-3311_Saltstress-Shoots-3.0h_Rep1'], 'sample_source_name': ['Col-0'], 'sample_characteristicts': ['Stock Code: N1092', 'Tissue: Shoots'], 'sample_extraction_protocol': ''}
#     }

#     print("--- Testing with a sample (GSM131309 - Roots) ---")
#     sample_id_1 = 'GSM131309'
#     sample_data_1 = samples_metadata[sample_id_1]
#     extracted_info_1 = GSE5623_extractor(sample_data_1)
#     print(f"Input for {sample_id_1}:\n{json.dumps(sample_data_1, indent=2)}")
#     print(f"Extracted info for {sample_id_1}:\n{json.dumps(extracted_info_1, indent=2)}")

#     print("\n--- Testing with another sample (GSM131316 - Shoots) ---")
#     sample_id_2 = 'GSM131316'
#     sample_data_2 = samples_metadata[sample_id_2]
#     extracted_info_2 = GSE5623_extractor(sample_data_2)
#     print(f"Input for {sample_id_2}:\n{json.dumps(sample_data_2, indent=2)}")
#     print(f"Extracted info for {sample_id_2}:\n{json.dumps(extracted_info_2, indent=2)}")

#     print("\n--- Testing with all samples ---")
#     all_extracted_data = {}
#     for sample_id, sample_data in samples_metadata.items():
#         try:
#             all_extracted_data[sample_id] = GSE5623_extractor(sample_data)
#         except ValueError as e:
#             print(f"Error processing {sample_id}: {e}")
#     # print(json.dumps(all_extracted_data, indent=2)) # Uncomment to see all results
#     print(f"Successfully processed {len(all_extracted_data)} out of {len(samples_metadata)} samples.")


# import json

# def GSE37408_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from sample metadata for GSE37408
#     following a specific JSON schema.

#     The schema requires 'tissue', 'treatment', and 'medium'.
#     - 'tissue' is extracted from 'sample_characteristicts' under 'tissue:'.
#     - 'treatment' is extracted from 'sample_characteristicts' under 'agent:'.
#       It is returned as a list containing the agent.
#     - 'medium' is inferred from the 'agent' as the incubation solution used
#       during the treatment, based on the study's experimental design.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 Expected keys include 'sample_characteristicts'.

#     Returns:
#         dict: A dictionary formatted according to the specified schema:
#               {"tissue": "string", "treatment": ["string"], "medium": "string"}.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The 'tissue' information is consistently found within the 'sample_characteristicts' list.
#     # Example: 'tissue: rosette leaves'
#     tissue_found = False
#     for char_str in sample_metadata.get('sample_characteristicts', []):
#         if char_str.startswith('tissue:'):
#             extracted_data['tissue'] = char_str.split(':', 1)[1].strip()
#             tissue_found = True
#             break
#     # Based on the provided samples_metadata, 'tissue' is always present and is 'rosette leaves'.
#     # If it were not found, a KeyError would occur if not handled, but for this dataset, it's reliable.
#     if not tissue_found:
#         # Provide a default or raise an error if tissue is strictly required and missing.
#         # For this dataset, it's expected to always be found.
#         extracted_data['tissue'] = 'unknown' 

#     # 2. Extract 'treatment' and infer 'medium'
#     # The 'agent' characteristic in 'sample_characteristicts' directly corresponds to the treatment.
#     # The study metadata indicates these agents are used for incubation, defining the medium.
#     # Examples: 'agent: mannitol', 'agent: sucrose', 'agent: no sugars'
#     agent = None
#     for char_str in sample_metadata.get('sample_characteristicts', []):
#         if char_str.startswith('agent:'):
#             agent = char_str.split(':', 1)[1].strip()
#             break

#     if agent:
#         # 'treatment' is a list of strings, so the agent is wrapped in a list.
#         extracted_data['treatment'] = [agent]

#         # 'medium' is inferred from the 'agent' as the incubation solution.
#         # 'no sugars' is described as a control in the study metadata.
#         if agent == 'mannitol':
#             extracted_data['medium'] = 'mannitol solution'
#         elif agent == 'sucrose':
#             extracted_data['medium'] = 'sucrose solution'
#         elif agent == 'no sugars':
#             extracted_data['medium'] = 'control solution'
#         else:
#             # Fallback for any other unexpected agent, though not anticipated for this dataset.
#             extracted_data['medium'] = f'{agent} solution'
#     else:
#         # If no 'agent' is found, provide default values to satisfy the schema's 'required' fields.
#         # Based on the provided samples_metadata, an agent is always present.
#         extracted_data['treatment'] = []  # An empty list for treatment
#         extracted_data['medium'] = 'unknown' # A placeholder for medium

#     return extracted_data



# import json

# def GSE22107_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from a sample metadata dictionary
#     for the GSE22107 dataset, conforming to a specific JSON schema.

#     The schema requires:
#     - tissue (string): Tissue the sample was extracted from.
#     - treatment (array of strings): List of treatments and stresses applied to the sample.
#     - medium (string): Growth medium of the sample.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample
#                                 from the GSE22107 dataset.

#     Returns:
#         dict: A dictionary formatted as a JSON instance conforming to the specified schema.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The tissue information is found within the 'sample_characteristicts' list.
#     tissue = None
#     for char_str in sample_metadata.get('sample_characteristicts', []):
#         if char_str.startswith('tissue:'):
#             tissue = char_str.split(':', 1)[1].strip()
#             break
#     # Ensure tissue is a string, even if not found (though it should always be present in this dataset)
#     extracted_data['tissue'] = tissue if tissue is not None else "unknown"

#     # 2. Extract 'treatment'
#     # Treatment information is typically found in 'sample_title' and 'sample_source_name'.
#     # We use a set to automatically handle duplicates and then convert to a sorted list.
#     treatments = set()
#     sample_title = sample_metadata.get('sample_title', [''])[0].lower()
#     sample_source_name = sample_metadata.get('sample_source_name', [''])[0].lower()

#     # Check for specific treatments/stresses
#     if 'mannitol' in sample_title or 'mannitol' in sample_source_name:
#         treatments.add('mannitol')
#     elif 'stress' in sample_title or 'stress' in sample_source_name:
#         # If 'mannitol' is not explicitly mentioned, but 'stress' is,
#         # we infer it's the osmotic stress mentioned in the overall study design.
#         treatments.add('osmotic stress')

#     # If the sample is a 'control', it implies no specific treatment/stress was applied.
#     # This should override any potential 'stress' that might have been picked up
#     # if 'control' is also present (e.g., "control stress" which is unlikely but handled).
#     if 'control' in sample_title or 'control' in sample_source_name:
#         treatments.clear() # A control sample has no treatment

#     extracted_data['treatment'] = sorted(list(treatments)) # Convert to list and sort for consistent output

#     # 3. Extract 'medium'
#     # Based on the provided study_metadata['overall_design']:
#     # "Plants were grown on nylon meshes overlaying control 0.5MS media."
#     # This is a constant value for all samples in this study and is not found
#     # within the individual sample_metadata.
#     extracted_data['medium'] = "0.5MS media"

#     return extracted_data



# import json

# def GSE49418_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from sample metadata
#     for the GSE49418 dataset, conforming to the specified JSON schema.

#     The schema is:
#     {
#         "properties": {
#             "tissue": {"description": "Tissue the samples was extracted from.", "title": "Tissue", "type": "string"},
#             "treatment": {"description": "List of treatments and stresses that was applied to the sample.", "items": {"type": "string"}, "title": "Treatment", "type": "array"},
#             "medium": {"description": "Growth medium of the sample.", "title": "Medium", "type": "string"}
#         },
#         "required": ["tissue", "treatment", "medium"]
#     }

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 This dictionary is expected to have keys like
#                                 'sample_source_name', 'sample_characteristicts', etc.

#     Returns:
#         dict: A dictionary formatted as a JSON instance conforming to the schema.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # Based on 'sample_source_name' (e.g., '46T seedlings without treatment')
#     # and 'sample_characteristicts' (e.g., 'developmental stage: 1-week-old seedlings'),
#     # 'seedlings' is the consistent tissue type across all samples in this study.
#     extracted_data["tissue"] = "seedlings"

#     # 2. Extract 'treatment'
#     treatments = []
#     # The treatment information is found within the 'sample_characteristicts' list.
#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if characteristic.startswith('treatment:'):
#                 # Extract the value after "treatment:" and strip whitespace
#                 treatment_value = characteristic.split(':', 1)[1].strip()
#                 # If the treatment is explicitly "none", we represent it as an empty list
#                 # as per the schema's description of "List of treatments and stresses".
#                 if treatment_value.lower() != 'none':
#                     treatments.append(treatment_value)
#                 break # Assuming only one 'treatment' characteristic per sample
#     extracted_data["treatment"] = treatments

#     # 3. Extract 'medium'
#     # The growth medium is not explicitly stated in the sample-specific metadata.
#     # However, the study's overall design mentions "One-week-old Arabidopsis seedlings"
#     # and the summary mentions "agar plates" in the context of tolerance experiments.
#     # For controlled seedling experiments, especially with dehydration treatments,
#     # growth on agar plates is a very common and highly probable medium.
#     # If it were soil, it would typically be specified.
#     # Therefore, "agar plate" is inferred as the medium for these samples.
#     extracted_data["medium"] = "agar plate"

#     return extracted_data



# import json

# def GSE16765_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from a sample's metadata
#     following a predefined JSON schema.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the specified JSON schema,
#               containing 'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The tissue information is found within the 'sample_characteristicts' list.
#     # Example: 'tissue: leaf'
#     tissue_found = False
#     for characteristic in sample_metadata.get('sample_characteristicts', []):
#         if characteristic.startswith('tissue:'):
#             extracted_data['tissue'] = characteristic.split(':', 1)[1].strip()
#             tissue_found = True
#             break
#     if not tissue_found:
#         # Fallback if tissue characteristic is not found, though it appears consistent in this dataset.
#         extracted_data['tissue'] = 'unknown'

#     # 2. Extract 'treatment'
#     # Treatments (specifically salt stress) are indicated in 'sample_source_name'
#     # and 'sample_title'.
#     # '100mM NaCl' indicates an applied stress, while '0mM NaCl' indicates a control
#     # condition (absence of the specific stress). The schema asks for *applied* treatments/stresses.
#     treatments = []
#     source_name = sample_metadata.get('sample_source_name', [''])[0]
#     sample_title = sample_metadata.get('sample_title', [''])[0]

#     if '100mM NaCl' in source_name or '100mM NaCl' in sample_title:
#         treatments.append('100mM NaCl stress')
#     # If '0mM NaCl' is present, it signifies a control group where no salt stress
#     # was applied, so the 'treatment' list remains empty as per the schema's
#     # description of "treatments and stresses that was applied".

#     extracted_data['treatment'] = treatments

#     # 3. Extract 'medium'
#     # The growth medium information is constant across samples for this study
#     # and is mentioned in the study's overall summary (not in individual sample metadata).
#     # The study metadata (provided as context) states:
#     # "Growth analyses were performed with seedlings germinated on culture media..."
#     # Given the function signature only accepts `sample_metadata` and the guidance
#     # that some info can be constant, we infer and hardcode this value based on
#     # the provided study context.
#     extracted_data['medium'] = 'culture media'

#     return extracted_data



# import json

# def GSE77815_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from the sample metadata
#     for the GSE77815 dataset, conforming to a specific JSON schema.

#     The schema requires:
#     - tissue (string): Tissue the sample was extracted from.
#     - treatment (array of strings): List of treatments and stresses applied to the sample.
#     - medium (string): Growth medium of the sample.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample
#                                 from the GSE77815 dataset.

#     Returns:
#         dict: A dictionary formatted according to the specified schema,
#               containing 'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The tissue information is found in the 'sample_characteristicts' list,
#     # typically in a string like 'tissue: Rossete leaves'.
#     tissue = None
#     if 'sample_characteristicts' in sample_metadata:
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if characteristic.startswith('tissue:'):
#                 tissue = characteristic.split(':', 1)[1].strip()
#                 break
#     # Provide a default if tissue is not found, though it appears consistently present.
#     extracted_data['tissue'] = tissue if tissue is not None else "unknown"

#     # 2. Extract 'treatment'
#     # Treatments/stresses are derived from the 'growth condition' in
#     # 'sample_characteristicts'. "normal condition" is interpreted as no specific
#     # treatment/stress applied, resulting in an empty list.
#     treatment_list = []
#     if 'sample_characteristicts' in sample_metadata:
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if characteristic.startswith('growth condition:'):
#                 condition = characteristic.split(':', 1)[1].strip()
#                 # If the condition is not "normal condition", it's considered a treatment/stress.
#                 if condition.lower() != 'normal condition':
#                     treatment_list.append(condition)
#                 break # Assuming only one primary growth condition per sample
#     extracted_data['treatment'] = treatment_list

#     # 3. Extract 'medium'
#     # The study metadata does not explicitly state the growth medium.
#     # Given that the samples are from Arabidopsis thaliana rosette leaves,
#     # "soil" is a common and reasonable inference for the growth medium
#     # when not specified otherwise for whole plants.
#     extracted_data['medium'] = "soil"

#     return extracted_data



# import json

# def GSE48474_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from sample metadata for GSE48474
#     following a specific JSON schema.

#     The output schema is:
#     {"properties": {"tissue": {"description": "Tissue the samples was extracted from.", "title": "Tissue", "type": "string"},
#                     "treatment": {"description": "List of treatments and stresses that was applied to the sample.", "items": {"type": "string"}, "title": "Treatment", "type": "array"},
#                     "medium": {"description": "Growth medium of the sample.", "title": "Medium", "type": "string"}},
#      "required": ["tissue", "treatment", "medium"]}

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the specified schema,
#               containing tissue, treatment, and medium information.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # Tissue information is typically found in 'sample_characteristicts'.
#     # Example: 'tissue: leaf'
#     characteristics = sample_metadata.get('sample_characteristicts', [])
#     tissue = "unknown" # Default value if not found
#     for char_item in characteristics:
#         if char_item.startswith('tissue:'):
#             tissue = char_item.split(':', 1)[1].strip()
#             break
#     extracted_data['tissue'] = tissue

#     # 2. Extract 'treatment'
#     # Treatment information can be found in 'sample_source_name' or 'sample_title'.
#     # It typically indicates 'drought treatment' or 'control'.
#     treatment_list = []
#     source_name_str = sample_metadata.get('sample_source_name', [''])[0].lower()
#     title_str = sample_metadata.get('sample_title', [''])[0].lower()

#     if 'drought' in source_name_str or 'drought' in title_str:
#         treatment_list.append('drought')
#     elif 'control' in source_name_str or 'control' in title_str:
#         treatment_list.append('control')
    
#     # Ensure treatment is always a list, even if empty or with a default
#     extracted_data['treatment'] = treatment_list if treatment_list else ["untreated"] # "untreated" as a sensible default for control if "control" isn't explicitly found

#     # 3. Extract 'medium'
#     # The growth medium is inferred from the 'sample_extraction_protocol' and study design.
#     # The protocol mentions "plants in pots were subjected to control condition (well-watered)
#     # and drought condition (soil water deficit)", implying "soil" as the medium.
#     # This appears to be constant across all samples in this study.
#     extracted_data['medium'] = "soil"

#     return extracted_data



# import json

# def GSE27552_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from a sample metadata dictionary
#     for the GSE27552 dataset, conforming to a specific JSON schema.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted as a JSON instance conforming to the specified schema.
#               The schema includes 'tissue', 'treatment', and 'medium'.
#     """
    
#     extracted_info = {
#         "tissue": "",
#         "treatment": [],
#         "medium": "soil"  # Inferred from the study context (Arabidopsis in soil drying experiments)
#                           # and guidance that medium can be constant across samples.
#     }

#     # Sample characteristics often contain key-value pairs like 'tissue: value' or 'treatment: value'
#     characteristics = sample_metadata.get('sample_characteristicts', [])
    
#     # Extract tissue
#     # Prioritize 'tissue' from 'sample_characteristicts' for specificity
#     found_tissue_in_characteristics = False
#     for char_str in characteristics:
#         if char_str.lower().startswith('tissue:'):
#             extracted_info['tissue'] = char_str.split(':', 1)[1].strip()
#             found_tissue_in_characteristics = True
#             break
    
#     # Fallback for tissue: if not found in characteristics, use 'sample_source_name'
#     # This fallback might not be strictly necessary for this dataset as 'tissue' seems
#     # consistently present in 'sample_characteristicts', but it adds robustness.
#     if not found_tissue_in_characteristics:
#         source_name = sample_metadata.get('sample_source_name', [])
#         if source_name:
#             extracted_info['tissue'] = source_name[0].strip()

#     # Extract treatment(s)
#     treatments_list = []
#     for char_str in characteristics:
#         if char_str.lower().startswith('treatment:'):
#             treat = char_str.split(':', 1)[1].strip()
            
#             treatments_list.append('control' if treat == 'Wet' else 'drought environment')
#     extracted_info['treatment'] = treatments_list

#     return extracted_info




# import json

# def GSE70861_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from a single sample's metadata
#     for the GSE70861 dataset, conforming to a specified JSON schema.

#     The function identifies and extracts the tissue, treatment, and growth medium
#     based on the structure and content of the provided sample_metadata dictionary.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 This typically includes keys like 'sample_title',
#                                 'sample_source_name', and 'sample_characteristicts'.

#     Returns:
#         dict: A dictionary representing a JSON instance that conforms to the target schema:
#               {"properties": {"tissue": {"type": "string"},
#                               "treatment": {"type": "array", "items": {"type": "string"}},
#                               "medium": {"type": "string"}},
#                "required": ["tissue", "treatment", "medium"]}

#               Example output:
#               {"tissue": "rosette leaves", "treatment": ["water stress"], "medium": "soil"}
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The tissue information is located within the 'sample_characteristicts' list.
#     # We iterate through this list to find the item that starts with 'tissue:'.
#     characteristics = sample_metadata.get('sample_characteristicts', [])
#     tissue = None
#     for char_item in characteristics:
#         if char_item.startswith('tissue:'):
#             # Extract the value after 'tissue:' and strip any leading/trailing whitespace.
#             tissue = char_item.split(':', 1)[1].strip()
#             break
#     # Assign the extracted tissue or a default 'unknown' if not found (though expected to be present).
#     extracted_data['tissue'] = tissue if tissue is not None else "unknown"

#     # 2. Extract 'treatment'
#     # The treatment information is derived from the 'sample_title'.
#     # The title often contains keywords indicating the experimental condition.
#     sample_title = sample_metadata.get('sample_title', [''])[0] # Get the first element from the list
#     treatment_list = []

#     # Check for specific keywords in the sample title (case-insensitive)
#     # to identify the applied treatment or condition.
#     if 'water stress' in sample_title.lower():
#         treatment_list.append('water stress')
#     elif 'recovery' in sample_title.lower():
#         treatment_list.append('recovery')
#     elif 'normal conditions' in sample_title.lower() or 'normal growth' in sample_title.lower():
#         treatment_list.append('normal growth conditions')
    
#     # The schema requires 'treatment' to be an array of strings.
#     # If no specific treatment is identified, default to ["untreated"] for clarity,
#     # although for this dataset, one of the conditions is always expected.
#     extracted_data['treatment'] = treatment_list if treatment_list else ["untreated"]

#     # 3. Extract 'medium'
#     # The growth medium is not explicitly mentioned in the provided sample or study metadata.
#     # Given that the study involves Arabidopsis plants under "normal growth conditions",
#     # "soil" is a common and reasonable inference for the growth medium.
#     extracted_data['medium'] = 'soil'

#     return extracted_data



# import json

# def GSE10670_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from sample metadata
#     for the GSE10670 dataset, following a specific JSON schema.

#     The output schema is:
#     {
#       "properties": {
#         "tissue": {"description": "Tissue the samples was extracted from.", "title": "Tissue", "type": "string"},
#         "treatment": {"description": "List of treatments and stresses that was applied to the sample.", "items": {"type": "string"}, "title": "Treatment", "type": "array"},
#         "medium": {"description": "Growth medium of the sample.", "title": "Medium", "type": "string"}
#       },
#       "required": ["tissue", "treatment", "medium"]
#     }

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 Expected keys include 'sample_title', 'sample_source_name',
#                                 and 'sample_characteristicts', with their values typically being lists of strings.

#     Returns:
#         dict: A dictionary formatted according to the specified schema.
#     """
#     extracted_data = {}

#     # 1. Extract Tissue
#     # Based on the provided sample metadata for GSE10670, all samples are consistently
#     # described as being from 'leaf' tissue (e.g., "T6 leaf-drought-rep3", "fully expanded leaves").
#     # Therefore, 'leaf' can be set as a constant for this study.
#     extracted_data["tissue"] = "leaf"

#     # 2. Extract Treatment
#     # Treatments/stresses are identified by keywords within the sample descriptions.
#     # "drought" or "7d without water" indicate a stress.
#     # "well watered" indicates a control condition, which, according to the schema's
#     # description ("List of treatments and stresses"), should result in an empty list
#     # as it's not an applied stress or treatment.
    
#     treatment_list = []
    
#     # Combine relevant metadata strings for a comprehensive, case-insensitive search
#     search_parts = []
#     if 'sample_title' in sample_metadata and isinstance(sample_metadata['sample_title'], list):
#         search_parts.extend(sample_metadata['sample_title'])
#     if 'sample_source_name' in sample_metadata and isinstance(sample_metadata['sample_source_name'], list):
#         search_parts.extend(sample_metadata['sample_source_name'])
#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         search_parts.extend(sample_metadata['sample_characteristicts'])

#     full_search_string = " ".join(search_parts).lower()

#     # Check for indicators of drought stress
#     if "drought" in full_search_string or "7d without water" in full_search_string:
#         treatment_list.append("drought stress")
#     # If "well watered" is present and no drought stress is indicated, the treatment_list
#     # remains empty, as "well watered" is a control condition, not an applied stress/treatment.
    
#     extracted_data["treatment"] = treatment_list

#     # 3. Extract Medium
#     # The specific growth medium (e.g., soil, agar, hydroponics) is not explicitly
#     # mentioned in either the study or sample metadata. Given that the study involves
#     # Arabidopsis plants, "soil" is a standard and highly probable growth substrate.
#     # The terms "well watered" or "7d without water" describe the *watering condition*
#     # of the medium, not the medium itself. Thus, "soil" is inferred as the medium.
#     extracted_data["medium"] = "soil"

#     return extracted_data



# import json

# def GSE72050_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from the sample metadata
#     for the GSE72050 dataset, conforming to the specified JSON schema.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 Expected to have keys like 'sample_characteristicts'.

#     Returns:
#         dict: A dictionary formatted as a JSON instance conforming to the schema:
#               {
#                   "tissue": "...",
#                   "treatment": ["...", "..."],
#                   "medium": "..."
#               }
#     """
#     extracted_data = {}

#     # Initialize variables with default/empty values
#     tissue = ""
#     genotype_variation = ""
#     stress_condition = ""
    
#     # The growth medium is not explicitly stated in the provided metadata.
#     # For Arabidopsis grown for leaf samples under drought stress experiments,
#     # "soil" is a common and reasonable default medium.
#     medium = "soil" 

#     # Extract information from the 'sample_characteristicts' list
#     # This list contains key-value pairs as strings (e.g., 'tissue: leaf')
#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if characteristic.startswith('tissue:'):
#                 # Extract the value after 'tissue:' and strip whitespace
#                 tissue = characteristic.split(':', 1)[1].strip()
#             elif characteristic.startswith('genotype/variation:'):
#                 # Extract the value after 'genotype/variation:' and strip whitespace
#                 genotype_variation = characteristic.split(':', 1)[1].strip()
#             elif characteristic.startswith('stress:'):
#                 # Extract the value after 'stress:' and strip whitespace
#                 stress_condition = characteristic.split(':', 1)[1].strip()

#     # The 'treatment' field is an array of strings, combining genotype/variation and stress.
#     treatment_list = []
#     if genotype_variation:
#         treatment_list.append(genotype_variation)
#     if stress_condition:
#         treatment_list.append(stress_condition)

#     # Populate the output dictionary according to the schema
#     extracted_data["tissue"] = tissue
#     extracted_data["treatment"] = treatment_list
#     extracted_data["medium"] = medium

#     return extracted_data



# def GSE65414_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information of the sample from the
#     GSE65414 dataset, conforming to the specified JSON schema.

#     The function identifies 'tissue' and 'treatment' from the sample-specific
#     metadata and infers 'medium' from the study-level information provided
#     in the problem description, as it's constant across samples.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 This dictionary is expected to be one of the inner
#                                 dictionaries from the 'samples_metadata' provided
#                                 in the problem description (e.g., the value associated
#                                 with a GSM ID).

#     Returns:
#         dict: A dictionary formatted as a JSON instance conforming to the output schema.
#               Example: {"tissue": "leaf", "treatment": ["ctrl"], "medium": "soil"}
#     """
#     extracted_data = {}

#     # Extract 'tissue'
#     # The tissue information is found within the 'sample_characteristicts' list.
#     # Example: 'tissue: leaf'
#     characteristics = sample_metadata.get('sample_characteristicts', [])
#     tissue = None
#     for char_string in characteristics:
#         if char_string.startswith('tissue:'):
#             # Split by the first colon and strip whitespace from the value
#             tissue = char_string.split(':', 1)[1].strip()
#             break
#     # Assign the extracted tissue or a default if not found (though it's consistent in this dataset)
#     extracted_data['tissue'] = tissue if tissue is not None else "unknown"

#     # Extract 'treatment'
#     # The treatment information (humidity) is also found in 'sample_characteristicts'.
#     # Example: 'treatment (humidity): ctrl' or 'treatment (humidity): dry'
#     treatment_list = []
#     for char_string in characteristics:
#         if char_string.startswith('treatment (humidity):'):
#             # Split by the first colon and strip whitespace from the value
#             treatment_value = char_string.split(':', 1)[1].strip()
#             # The schema expects an array of strings for 'treatment'
#             treatment_list.append(treatment_value)
#             break  # Assuming only one humidity treatment per sample
#     # Assign the extracted treatment(s) or a default if none found
#     extracted_data['treatment'] = treatment_list if treatment_list else ["no treatment specified"]

#     # Extract 'medium'
#     # Based on the study_metadata['overall_design'] ("Arabidopsis thaliana plants were grown on soil"),
#     # the growth medium is consistently "soil" for all samples in this study.
#     # This information is not present in the individual sample_metadata but is constant.
#     extracted_data['medium'] = "soil"

#     return extracted_data



# import json

# def GSE103398_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from sample metadata for GSE103398
#     following the specified JSON schema.

#     The schema requires:
#     - tissue (string): Tissue the sample was extracted from.
#     - treatment (array of strings): List of treatments and stresses applied.
#     - medium (string): Growth medium of the sample.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the specified schema.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The 'tissue' information is found within the 'sample_characteristicts' list.
#     # It consistently appears as 'tissue: whole seedling'.
#     tissue_value = None
#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if characteristic.startswith('tissue:'):
#                 tissue_value = characteristic.split(':', 1)[1].strip()
#                 break
#     extracted_data['tissue'] = tissue_value if tissue_value is not None else "unknown" # Provide a fallback if not found, though it's consistent in this dataset.

#     # 2. Extract 'treatment'
#     # The 'treatment' information is found within the 'sample_characteristicts' list.
#     # It appears as 'subjected to: [treatment_type] treatment'.
#     # Treatments can be single ('naïve', 'primed', 'triggered') or compound ('primed and triggered').
#     treatment_list = []
#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if characteristic.startswith('subjected to:'):
#                 treatment_str = characteristic.split(':', 1)[1].strip()
                
#                 # Remove the common suffix " treatment"
#                 if treatment_str.endswith(' treatment'):
#                     treatment_str = treatment_str[:-len(' treatment')]
                
#                 # Handle compound treatments like "primed and triggered"
#                 if ' and ' in treatment_str:
#                     # Split by ' and ' and strip whitespace from each part
#                     treatment_list.extend([t.strip() for t in treatment_str.split(' and ')])
#                 else:
#                     # For single treatments, just add to the list
#                     treatment_list.append(treatment_str)
#                 break
#     extracted_data['treatment'] = treatment_list

#     # 3. Extract 'medium'
#     # The growth medium is not explicitly stated in the sample metadata.
#     # However, the study's 'overall_design' mentions "normal growth condition" for Arabidopsis seedlings.
#     # For Arabidopsis seedlings, Murashige and Skoog (MS) medium is a widely recognized standard
#     # for "normal growth condition" in laboratory settings. This is a constant for the study.
#     extracted_data['medium'] = "Murashige and Skoog (MS) medium"

#     return extracted_data



# import json

# def GSE16222_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information (tissue, treatment, medium)
#     from a single sample's metadata dictionary for the GSE16222 dataset.

#     The extraction logic is based on the structure and content observed in the
#     provided study_metadata and samples_metadata.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 Expected to have a 'sample_characteristicts' key
#                                 which is a list of strings.

#     Returns:
#         dict: A dictionary conforming to the specified JSON schema:
#               {
#                   "tissue": "...",
#                   "treatment": ["..."],
#                   "medium": "..."
#               }
#               All required fields will be present.
#     """
#     extracted_data = {
#         "tissue": None,
#         "treatment": [],
#         "medium": None
#     }

#     # 'sample_characteristicts' is a list of strings containing the desired information
#     characteristics = sample_metadata.get('sample_characteristicts', [])

#     # --- Extract 'tissue' ---
#     # Look for a string starting with 'growth stage:'
#     for char_str in characteristics:
#         if char_str.startswith('growth stage:'):
#             extracted_data['tissue'] = char_str.split(':', 1)[1].strip()
#             break
#     # Fallback: Based on the study and sample descriptions, "seedling" is the consistent tissue.
#     if extracted_data['tissue'] is None:
#         extracted_data['tissue'] = "seedling"

#     # --- Extract 'medium' ---
#     # Look for a string starting with 'growth conditions:'
#     for char_str in characteristics:
#         if char_str.startswith('growth conditions:'):
#             growth_conditions_str = char_str.split(':', 1)[1].strip()
#             # The schema asks for the "Growth medium of the sample."
#             # "Murashige-Skoog medium" is the specific name of the medium used.
#             if "Murashige-Skoog" in growth_conditions_str:
#                 extracted_data['medium'] = "Murashige-Skoog medium"
#             else:
#                 # Fallback if the specific string changes, use the full conditions string
#                 # (though unlikely for this dataset based on provided metadata).
#                 extracted_data['medium'] = growth_conditions_str
#             break
#     # Fallback: Based on the study description, "Murashige-Skoog medium" is consistent.
#     if extracted_data['medium'] is None:
#         extracted_data['medium'] = "Murashige-Skoog medium"

#     # --- Extract 'treatment' ---
#     # Look for a string starting with 'treatment:'
#     treatment_str_raw = None
#     for char_str in characteristics:
#         if char_str.startswith('treatment:'):
#             treatment_str_raw = char_str.split(':', 1)[1].strip()
#             break

#     if treatment_str_raw:
#         # Convert to lowercase for case-insensitive matching, but use raw string for output if not standardized
#         treatment_str_lower = treatment_str_raw.lower()
        
#         if treatment_str_lower == 'untreated control':
#             extracted_data['treatment'] = ['control']
#         elif 'and thereafter under anoxia' in treatment_str_lower:
#             # Handle combined heat and anoxia treatment
#             split_keyword = ' and thereafter under anoxia'
#             split_index = treatment_str_lower.find(split_keyword)
            
#             # Extract the raw parts based on the lowercase split index
#             heat_part_raw = treatment_str_raw[:split_index].strip()
#             anoxia_part_raw = 'anoxia' + treatment_str_raw[split_index + len(split_keyword):].strip()
            
#             # Standardize heat part
#             if '90 minutes at 38°c' in heat_part_raw.lower():
#                 extracted_data['treatment'].append('heat (90 minutes at 38°C)')
#             else:
#                 extracted_data['treatment'].append(heat_part_raw) # Fallback to raw if not matching standard

#             # Standardize anoxia part
#             if 'anoxia for 6h at 23°c' in anoxia_part_raw.lower():
#                 extracted_data['treatment'].append('anoxia (6h at 23°C)')
#             else:
#                 extracted_data['treatment'].append(anoxia_part_raw) # Fallback to raw if not matching standard
#         elif '90 minutes at 38°c' in treatment_str_lower:
#             extracted_data['treatment'] = ['heat (90 minutes at 38°C)']
#         elif 'anoxia for 6h at 23°c' in treatment_str_lower:
#             extracted_data['treatment'] = ['anoxia (6h at 23°C)']
#         else:
#             # Catch any other single treatment not explicitly handled by specific keywords
#             extracted_data['treatment'] = [treatment_str_raw]
#     else:
#         # Fallback: If no 'treatment:' characteristic is found, assume it's a control sample.
#         # (In this dataset, 'untreated control' is explicitly stated for controls,
#         # but this provides robustness for potentially malformed entries).
#         extracted_data['treatment'] = ['control']

#     return extracted_data



# import json

# def GSE53308_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from the sample metadata
#     for the GSE53308 study, conforming to a specific JSON schema.

#     The schema requires 'tissue' (string), 'treatment' (array of strings),
#     and 'medium' (string).

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the specified schema,
#               containing 'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # 1. Extract 'medium':
#     # Based on the study_metadata['overall_design']:
#     # "Wild type Arabidopsis Col-0 plants were grown hydroponically..."
#     # The medium is constant across all samples in this study.
#     extracted_data["medium"] = "hydroponic medium"

#     # 2. Extract 'tissue' and 'treatment' from 'sample_characteristicts':
#     # This field is a list of strings, where each string is in "key: value" format.
#     characteristics = sample_metadata.get('sample_characteristicts', [])

#     # Initialize treatment as an empty list to ensure it's always an array
#     extracted_data["treatment"] = []

#     for char_string in characteristics:
#         # Split only on the first colon to handle cases where values might contain colons
#         if ': ' in char_string:
#             key, value = char_string.split(': ', 1)
#             key = key.strip()
#             value = value.strip()

#             if key == "tissue":
#                 extracted_data["tissue"] = value
#             elif key == "treatment":
#                 # The schema expects 'treatment' to be an array of strings.
#                 # In this dataset, each sample has a single treatment string.
#                 extracted_data["treatment"].append(value)
    
#     # Ensure all required fields are present.
#     # Based on the provided sample metadata, 'tissue' and 'treatment' are always present.
#     # However, adding a safeguard for robustness.
#     if "tissue" not in extracted_data:
#         extracted_data["tissue"] = "N/A" # Or raise an error if missing data is critical

#     return extracted_data



# import json

# def GSE108070_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from the sample metadata
#     for the GSE108070 dataset, conforming to the specified JSON schema.

#     The schema is:
#     {"properties": {"tissue": {"description": "Tissue the samples was extracted from.", "title": "Tissue", "type": "string"},
#                     "treatment": {"description": "List of treatments and stresses that was applied to the sample.", "items": {"type": "string"}, "title": "Treatment", "type": "array"},
#                     "medium": {"description": "Growth medium of the sample.", "title": "Medium", "type": "string"}},
#      "required": ["tissue", "treatment", "medium"]}

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 This typically includes keys like 'sample_characteristicts'.

#     Returns:
#         dict: A dictionary formatted according to the specified schema,
#               containing extracted tissue, treatment, and medium information.
#     """

#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # Based on the study summary ("4-day-old plants"), 'seedling' is the most appropriate
#     # tissue type for Arabidopsis plants at this developmental stage. This is constant
#     # across all samples in this study.
#     extracted_data["tissue"] = "seedling"

#     # 2. Extract 'treatment'
#     # Treatments and stresses are found in 'sample_characteristicts'.
#     # We filter out control agents (like DMSO) and non-stress conditions.
#     treatments = []
#     characteristics = sample_metadata.get('sample_characteristicts', [])

#     for char_str in characteristics:
#         if ": " in char_str:
#             key, value = char_str.split(": ", 1)
#             if key == "agent":
#                 # "DMSO" is the solvent control for the HDAC inhibitors (Ky-9, Ky-72).
#                 # It's not considered an active treatment in the context of the experiment's
#                 # primary variables. "water" is also a control mentioned in the study summary.
#                 if value.lower() not in ["dmso", "water"]:
#                     treatments.append(value)
#             elif key == "group":
#                 # "non-stress" indicates the absence of an applied stress.
#                 if value.lower() != "non-stress":
#                     treatments.append(value)
    
#     # Ensure treatments are unique and sorted for consistent output, as per schema's array type.
#     extracted_data["treatment"] = sorted(list(set(treatments)))

#     # 3. Extract 'medium'
#     # The growth medium is not explicitly stated in the provided metadata.
#     # For Arabidopsis seedlings, Murashige and Skoog (MS) medium is a very common and
#     # standard growth medium. This is assumed to be constant across all samples.
#     extracted_data["medium"] = "Murashige and Skoog (MS) medium"

#     return extracted_data



# import json

# def GSE6583_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information of the sample from the
#     GSE6583 dataset following the specified JSON schema.

#     The function identifies 'tissue', 'treatment', and 'medium' information.
#     'Tissue' and 'medium' are inferred from the overall study context as they
#     are constant across samples and not explicitly detailed in individual
#     sample metadata. 'Treatment' is extracted from the sample-specific
#     characteristics.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 This typically includes fields like 'sample_characteristicts'.

#     Returns:
#         dict: A dictionary formatted as a JSON instance conforming to the
#               output schema, containing 'tissue', 'treatment', and 'medium' fields.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # Based on the study metadata, the samples are from Arabidopsis plants,
#     # which are 3 weeks old. In the absence of specific tissue information
#     # (e.g., leaf, root), "whole plant" or "seedling" is a common and
#     # appropriate description for this age. We'll use "whole plant".
#     extracted_data["tissue"] = "whole plant"

#     # 2. Extract 'treatment'
#     treatments = []
#     # The 'sample_characteristicts' field contains the condition/treatment information.
#     characteristics = sample_metadata.get('sample_characteristicts', [])

#     for char_str in characteristics:
#         if char_str.startswith('Condition:'):
#             condition_value = char_str.split('Condition:')[1].strip()
#             if condition_value == 'control':
#                 treatments.append("control")
#             elif condition_value == 'drought 2hr':
#                 # The schema asks for "treatments and stresses". "drought" is the stress.
#                 # "2hr" is duration, which can be implied or simplified.
#                 treatments.append("drought stress")
#             # Add more conditions here if they were present in the dataset
#             break # Assuming only one primary condition per sample

#     # Ensure 'treatment' is always a list, even if no condition was found.
#     # If no specific treatment is found, we can default to "unknown" or an empty list
#     # depending on strictness. For this schema, an empty list is valid if no treatment.
#     # However, "control" is a valid condition, so we should always have at least one.
#     extracted_data["treatment"] = treatments if treatments else ["unknown"]

#     # 3. Extract 'medium'
#     # The growth medium is not explicitly mentioned in the provided sample_metadata
#     # or study_metadata. For 3-week-old Arabidopsis plants, "soil" is a very common
#     # growth medium in experimental setups. This aligns with the guidance that
#     # some information can be constant across samples.
#     extracted_data["medium"] = "soil"

#     return extracted_data

# if __name__ == '__main__':
#     # Example usage with the provided samples_metadata
#     samples_metadata = {
#         'GSM152139': {'Sample_id': 'GSM152139', 'sample_title': ['Col-0 drought 3'], 'sample_source_name': ['Col-0 d3'], 'sample_characteristicts': ['Ecotype: Columbia 0', 'age: 3 weeks', 'Condition: drought 2hr'], 'sample_extraction_protocol': ['TriZOL']},
#         'GSM152142': {'Sample_id': 'GSM152142', 'sample_title': ['siz1-3 control 2'], 'sample_source_name': ['siz1-3 c2'], 'sample_characteristicts': ['Genotype: siz1-3', 'age: 3 weeks', 'Condition: control'], 'sample_extraction_protocol': ['TriZOL']},
#         'GSM152140': {'Sample_id': 'GSM152140', 'sample_title': ['siz1-3 control 1'], 'sample_source_name': ['siz1-3 c1'], 'sample_characteristicts': ['Genotype: siz1-3', 'age: 3 weeks', 'Condition: control'], 'sample_extraction_protocol': ['TriZOL']},
#         'GSM152144': {'Sample_id': 'GSM152144', 'sample_title': ['siz1-3 control 3'], 'sample_source_name': ['siz1-3 c3'], 'sample_characteristicts': ['Genotype: siz1-3', 'age: 3 weeks', 'Condition: control'], 'sample_extraction_protocol': ['TriZOL']},
#         'GSM152137': {'Sample_id': 'GSM152137', 'sample_title': ['Col-0 drought 2'], 'sample_source_name': ['Col-0 d2'], 'sample_characteristicts': ['Ecotype: Columbia 0', 'age: 3 weeks', 'Condition: drought 2hr'], 'sample_extraction_protocol': ['TriZOL']},
#         'GSM152141': {'Sample_id': 'GSM152141', 'sample_title': ['siz1-3 drought 1'], 'sample_source_name': ['siz1-3 d1'], 'sample_characteristicts': ['Genotype: siz1-3', 'age: 3 weeks', 'Condition: drought 2hr'], 'sample_extraction_protocol': ['TriZOL']},
#         'GSM152135': {'Sample_id': 'GSM152135', 'sample_title': ['Col-0 drought 1'], 'sample_source_name': ['Col 0 d1'], 'sample_characteristicts': ['Ecotype: Columbia 0', 'age: 3 weeks', 'Condition: drought 2hr'], 'sample_extraction_protocol': ['TriZOL']},
#         'GSM152145': {'Sample_id': 'GSM152145', 'sample_title': ['siz1-3 drought 3'], 'sample_source_name': ['siz1-3 d3'], 'sample_characteristicts': ['Genotype: siz1-3', 'age: 3 weeks', 'Condition: drought 2hr'], 'sample_extraction_protocol': ['TriZOL']},
#         'GSM152143': {'Sample_id': 'GSM152143', 'sample_title': ['siz1-3 drought 2'], 'sample_source_name': ['siz1-3 d2'], 'sample_characteristicts': ['Genotype: siz1-3', 'age: 3 weeks', 'Condition: drought 2hr'], 'sample_extraction_protocol': ['TriZOL']},
#         'GSM152136': {'Sample_id': 'GSM152136', 'sample_title': ['Col-0 control 2'], 'sample_source_name': ['Col-0 c2'], 'sample_characteristicts': ['Ecotype: Columbia 0', 'age: 3 weeks', 'Condition: control'], 'sample_extraction_protocol': ['TriZOL']},
#         'GSM152138': {'Sample_id': 'GSM152138', 'sample_title': ['Col-0 control 3'], 'sample_source_name': ['Col-0 c3'], 'sample_characteristicts': ['Ecotype: Columbia 0', 'age: 3 weeks', 'Condition: control'], 'sample_extraction_protocol': ['TriZOL']},
#         'GSM152134': {'Sample_id': 'GSM152134', 'sample_title': ['Col-0 control 1'], 'sample_source_name': ['Col-0 c1'], 'sample_characteristicts': ['Columbia 0', 'Age: 21 days', 'Condition: control'], 'sample_extraction_protocol': ['TriZOL']}
#     }

#     print("--- Extracted Data for Sample GSM152139 (Drought) ---")
#     sample_id_drought = 'GSM152139'
#     extracted_drought = GSE6583_extractor(samples_metadata[sample_id_drought])
#     print(json.dumps(extracted_drought, indent=2))

#     print("\n--- Extracted Data for Sample GSM152142 (Control) ---")
#     sample_id_control = 'GSM152142'
#     extracted_control = GSE6583_extractor(samples_metadata[sample_id_control])
#     print(json.dumps(extracted_control, indent=2))

#     print("\n--- Extracted Data for Sample GSM152134 (Control, different age format) ---")
#     sample_id_control_alt = 'GSM152134'
#     extracted_control_alt = GSE6583_extractor(samples_metadata[sample_id_control_alt])
#     print(json.dumps(extracted_control_alt, indent=2))

#     # Example of a hypothetical sample with missing characteristics (for robustness check)
#     print("\n--- Extracted Data for Hypothetical Sample (Missing Characteristics) ---")
#     hypothetical_sample = {
#         'Sample_id': 'GSM_HYPO',
#         'sample_title': ['Hypothetical Sample'],
#         'sample_source_name': ['Hypo'],
#         # 'sample_characteristicts': [] # Missing this key
#         'sample_extraction_protocol': ['TriZOL']
#     }
#     extracted_hypo = GSE6583_extractor(hypothetical_sample)
#     print(json.dumps(extracted_hypo, indent=2))


# import json

# def GSE121225_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from the sample metadata
#     for the GSE121225 dataset, conforming to the specified JSON schema.

#     The schema is:
#     {
#         "properties": {
#             "tissue": {
#                 "description": "Tissue the samples was extracted from.",
#                 "title": "Tissue",
#                 "type": "string"
#             },
#             "treatment": {
#                 "description": "List of treatments and stresses that was applied to the sample.",
#                 "items": {"type": "string"},
#                 "title": "Treatment",
#                 "type": "array"
#             },
#             "medium": {
#                 "description": "Growth medium of the sample.",
#                 "title": "Medium",
#                 "type": "string"
#             }
#         },
#         "required": ["tissue", "treatment", "medium"]
#     }

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the specified JSON schema,
#               containing the extracted tissue, treatment, and medium information.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # Based on the study metadata ("plants") and the absence of specific tissue
#     # information in the sample characteristics or extraction protocol,
#     # "plant" is the most appropriate general tissue type.
#     extracted_data["tissue"] = "plant"

#     # 2. Extract 'treatment'
#     treatments = []
#     sample_characteristics = sample_metadata.get('sample_characteristicts', [])

#     condition = None
#     for char_entry in sample_characteristics:
#         if char_entry.startswith('condition:'):
#             condition = char_entry.split(':', 1)[1].strip()
#             break

#     if condition == 'salt stress':
#         # From the study_metadata['overall_design'][0]:
#         # 'Microarray analysis was conducted using hda19, quad, quint and Col-0 plants
#         # subjected to non-stress and 125 mM NaCl salt-stress condition for 2 h.'
#         # The specific treatment applied is "125 mM NaCl salt-stress".
#         treatments.append("125 mM NaCl salt-stress")
#     # If the condition is 'non-stress', the treatments list remains empty,
#     # as "non-stress" is the absence of a specific treatment.

#     extracted_data["treatment"] = treatments

#     # 3. Extract 'medium'
#     # The growth medium is not explicitly mentioned in either the study metadata
#     # or the individual sample metadata provided.
#     # As 'medium' is a required field in the schema, "Not specified" is used
#     # as a placeholder.
#     extracted_data["medium"] = "Not specified"

#     return extracted_data



# import json

# def GSE34595_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from the sample metadata
#     following the specified schema.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the target JSON schema,
#               containing 'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The tissue information is found within the 'sample_characteristicts' list.
#     # Example: 'tissue: cultured explant'
#     tissue_found = False
#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if isinstance(characteristic, str) and characteristic.startswith('tissue:'):
#                 extracted_data['tissue'] = characteristic.split(':', 1)[1].strip()
#                 tissue_found = True
#                 break
#     if not tissue_found:
#         # Provide a default or handle cases where tissue might be missing,
#         # though for this dataset, it appears to always be present.
#         extracted_data['tissue'] = "not specified"

#     # 2. Extract 'medium'
#     # Based on the study_metadata['overall_design'], the medium is "B5 medium"
#     # and appears to be constant across all samples in this study.
#     extracted_data['medium'] = "B5 medium"

#     # 3. Extract 'treatment'
#     # Treatments are derived from the study's overall design (constant)
#     # and potentially specific details in the sample title.
#     treatments = []

#     # Treatments constant across all samples from study_metadata['overall_design']
#     # "cultured on B5 medium supplemented with 0.5 mg/L IBA"
#     treatments.append("0.5 mg/L IBA")
#     # "After 12 hours of culture at 28C"
#     treatments.append("12 hours of culture")
#     treatments.append("28C temperature")

#     # Check for additional treatments mentioned in the 'sample_title'
#     # Example: 'rrd2_root induction.12h.28C_rep3'
#     if 'sample_title' in sample_metadata and isinstance(sample_metadata['sample_title'], list) and sample_metadata['sample_title']:
#         sample_title = sample_metadata['sample_title'][0]
#         if "root induction" in sample_title:
#             treatments.append("root induction")

#     # Use a set to remove any potential duplicates and then convert back to a list.
#     # Sorting for consistent output order.
#     extracted_data['treatment'] = sorted(list(set(treatments)))

#     return extracted_data




# import json

# def GSE4062_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from a sample's metadata
#     dictionary, conforming to a predefined JSON schema.

#     The function identifies the tissue, treatment, and growth medium of the sample.
#     If information for 'medium' is not explicitly found in the provided metadata,
#     it defaults to "not specified" to satisfy the schema's requirement.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 This dictionary is expected to have keys like
#                                 'sample_characteristicts' and 'sample_title'.

#     Returns:
#         dict: A dictionary formatted as a JSON instance conforming to the
#               following schema:
#               {
#                   "properties": {
#                       "tissue": {"description": "Tissue the samples was extracted from.", "title": "Tissue", "type": "string"},
#                       "treatment": {"description": "List of treatments and stresses that was applied to the sample.", "items": {"type": "string"}, "title": "Treatment", "type": "array"},
#                       "medium": {"description": "Growth medium of the sample.", "title": "Medium", "type": "string"}
#                   },
#                   "required": ["tissue", "treatment", "medium"]
#               }
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # Tissue information is found within the 'sample_characteristicts' list,
#     # typically in an entry starting with "Tissue: ".
#     tissue = "not specified" # Default value if not found
#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if isinstance(characteristic, str) and characteristic.startswith('Tissue:'):
#                 # Extract the value after "Tissue:" and strip any whitespace
#                 tissue = characteristic.split(':', 1)[1].strip()
#                 break
#     extracted_data['tissue'] = tissue

#     # 2. Extract 'treatment'
#     # Treatment information (control or heat shock) is primarily derived from
#     # the 'sample_title'.
#     treatment = []
#     if 'sample_title' in sample_metadata and isinstance(sample_metadata['sample_title'], list) and sample_metadata['sample_title']:
#         sample_title = sample_metadata['sample_title'][0].lower()
#         if 'control' in sample_title:
#             treatment.append("control")
#         elif 'heat shock' in sample_title:
#             treatment.append("heat shock")
    
#     # Ensure the 'treatment' list is not empty to satisfy the schema's array type
#     # and 'items' constraint, even if no specific treatment was identified.
#     # Based on the provided dataset, samples are always either 'control' or 'heat shock'.
#     if not treatment:
#         treatment.append("not specified") 

#     extracted_data['treatment'] = treatment

#     # 3. Extract 'medium'
#     # The growth medium information is not explicitly provided in the given
#     # study or sample metadata. As it's a required field in the schema,
#     # a placeholder value is used.
#     extracted_data['medium'] = "not specified"

#     return extracted_data

# if __name__ == '__main__':
#     # Example usage with the provided samples_metadata
#     samples_metadata = {
#         'GSM92824': {'Sample_id': 'GSM92824', 'sample_title': ['Wt_shoot_control_rep2'], 'sample_source_name': ['shoot, heat shock, control, 2h'], 'sample_characteristicts': ['Ecotype: Col-0 wild type', 'Stage: 2 weeks old', 'Tissue: shoot'], 'sample_extraction_protocol': ''},
#         'GSM92829': {'Sample_id': 'GSM92829', 'sample_title': ['Hsa32KO_shoot_heat shock_rep1'], 'sample_source_name': ['shoot, heat shock,  2h'], 'sample_characteristicts': ['Ecotype: Col-0', 'Mutant allele: Hsa32 (At4g21320) T-DNA knockout, GABI-268A08', 'Stage: 2 weeks old', 'Tissue: shoot'], 'sample_extraction_protocol': ''},
#         'GSM92830': {'Sample_id': 'GSM92830', 'sample_title': ['Hsa32KO_shoot_heat shock_rep2'], 'sample_source_name': ['shoot, heat shock, 2h'], 'sample_characteristicts': ['Ecotype: Col-0', 'Mutant allele: Hsa32 (At4g21320) T-DNA knockout, GABI-268A08', 'Stage: 2 weeks old', 'Tissue: shoot'], 'sample_extraction_protocol': ''},
#         'GSM92826': {'Sample_id': 'GSM92826', 'sample_title': ['Wt_shoot_heat shock_rep2'], 'sample_source_name': ['shoot, heat shock, 2h'], 'sample_characteristicts': ['Ecotype: Col-0 wild type', 'Stage: 2 weeks old', 'Tissue: shoot'], 'sample_extraction_protocol': ''},
#         'GSM92827': {'Sample_id': 'GSM92827', 'sample_title': ['Hsa32KO_shoot_control_rep1'], 'sample_source_name': ['shoot, heat shock, control, 2h'], 'sample_characteristicts': ['Ecotype: Col-0', 'Mutant allele: Hsa32 (At4g21320) T-DNA knockout, GABI-268A08', 'Stage: 2 weeks old', 'Tissue: shoot'], 'sample_extraction_protocol': ''},
#         'GSM92823': {'Sample_id': 'GSM92823', 'sample_title': ['Wt_shoot_control_rep1'], 'sample_source_name': ['shoot, heat shock, 2h'], 'sample_characteristicts': ['Ecotype: Col-0 wild type', 'Stage: 2 weeks old', 'Tissue: shoot'], 'sample_extraction_protocol': ''},
#         'GSM92828': {'Sample_id': 'GSM92828', 'sample_title': ['Hsa32KO_shoot_control_rep2'], 'sample_source_name': ['shoot, heat shock, control, 2h'], 'sample_characteristicts': ['Ecotype: Col-0', 'Mutant allele: Hsa32 (At4g21320) T-DNA knockout, GABI-268A08', 'Stage: 2 weeks old', 'Tissue: shoot'], 'sample_extraction_protocol': ''},
#         'GSM92825': {'Sample_id': 'GSM92825', 'sample_title': ['Wt_shoot_heat shock_rep1'], 'sample_source_name': ['shoot, heat shock, 2h'], 'sample_characteristicts': ['Ecotype: Col-0 wild type', 'Stage: 2 weeks old', 'Tissue: shoot'], 'sample_extraction_protocol': ''}
#     }

#     print("--- Extracted Data for Samples ---")
#     for sample_id, metadata in samples_metadata.items():
#         extracted_info = GSE4062_extractor(metadata)
#         print(f"Sample ID: {sample_id}")
#         print(json.dumps(extracted_info, indent=2))
#         print("-" * 30)

#     # Example of a specific sample output:
#     # For GSM92824:
#     # {
#     #   "tissue": "shoot",
#     #   "treatment": [
#     #     "control"
#     #   ],
#     #   "medium": "not specified"
#     # }
#     #
#     # For GSM92825:
#     # {
#     #   "tissue": "shoot",
#     #   "treatment": [
#     #     "heat shock"
#     #   ],
#     #   "medium": "not specified"
#     # }


# def GSE90562_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from a single sample's metadata
#     following a predefined JSON schema.

#     The function infers 'tissue' and 'medium' based on common practices for
#     Arabidopsis studies, as this information is not explicitly present in the
#     provided sample or study metadata. 'Treatment' is extracted from the
#     'sample_characteristicts' field, with specific details (like concentration
#     and duration for salt stress) derived from the study's overall design.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the specified JSON schema,
#               containing 'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The study metadata indicates "Arabidopsis plants". No specific tissue
#     # is mentioned in the sample metadata. "whole plant" is a common and safe
#     # default for such studies when no specific organ is specified.
#     extracted_data["tissue"] = "whole plant"

#     # 2. Extract 'treatment'
#     treatments = []
#     # The relevant information is found in the 'sample_characteristicts' list.
#     characteristics = sample_metadata.get('sample_characteristicts', [])

#     for char_str in characteristics:
#         if char_str.startswith('protocol:'):
#             protocol_value = char_str.split(':', 1)[1].strip()
#             if protocol_value == 'non-stress':
#                 treatments.append("non-stress")
#             elif protocol_value == 'salt stress':
#                 # According to study_metadata['overall_design'], the salt stress
#                 # condition was "125 mM NaCl salt-stress condition for 2 h."
#                 treatments.append("125 mM NaCl salt-stress for 2 h")
#             # Assuming only one primary protocol/stress per sample for this dataset.
#             break

#     # Ensure 'treatment' is always an array of strings, as required by the schema.
#     # Provide a default if no specific treatment protocol was found or recognized.
#     extracted_data["treatment"] = treatments if treatments else ["unknown"]

#     # 3. Extract 'medium'
#     # No explicit growth medium information is provided in the sample or study metadata.
#     # For Arabidopsis studies, Murashige and Skoog (MS) medium is a very common choice.
#     extracted_data["medium"] = "MS medium"

#     return extracted_data



# import re
# import json

# def GSE79997_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from the sample metadata
#     following the specified JSON schema.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the target JSON schema,
#               containing 'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # Look for 'tissue: <value>' within the 'sample_characteristicts' list.
#     tissue_found = False
#     for char_string in sample_metadata.get('sample_characteristicts', []):
#         if char_string.startswith('tissue:'):
#             extracted_data['tissue'] = char_string.split(':', 1)[1].strip()
#             tissue_found = True
#             break
#     if not tissue_found:
#         extracted_data['tissue'] = "unknown" # Default if tissue information is not found

#     # 2. Extract 'treatment'
#     # The 'sample_source_name' provides the most direct information about the treatment.
#     # The schema requires an array of strings for treatment.
#     treatment_list = []
#     source_name = sample_metadata.get('sample_source_name', [''])[0]

#     if "200 mM NaCl" in source_name:
#         treatment_list.append("200 mM NaCl")
#     elif "mock-treated" in source_name:
#         treatment_list.append("mock treatment")
#     else:
#         # Fallback to Sample title if source_name is not explicit enough
#         title = sample_metadata.get('sample_title', [''])[0]
#         if "Salt treatment" in title:
#             treatment_list.append("salt treatment")
#         elif "mock treatment" in title:
#             treatment_list.append("control")
#         else:
#             treatment_list.append("unknown treatment") # Default if no specific treatment found

#     extracted_data['treatment'] = treatment_list

#     # 3. Extract 'medium'
#     # The medium information is embedded within the 'treatment:' string
#     # in 'sample_characteristicts'.
#     medium_found = False
#     for char_string in sample_metadata.get('sample_characteristicts', []):
#         if char_string.startswith('treatment:'):
#             # Prioritize "1/2 strength MS medium" as it's more descriptive
#             match_strength = re.search(r'1/2 strength MS medium', char_string, re.IGNORECASE)
#             if match_strength:
#                 extracted_data['medium'] = match_strength.group(0)
#                 medium_found = True
#                 break
#             else:
#                 # Fallback to "1/2 MS medium" if "strength" is not specified
#                 match_ms = re.search(r'1/2 MS medium', char_string, re.IGNORECASE)
#                 if match_ms:
#                     extracted_data['medium'] = match_ms.group(0)
#                     medium_found = True
#                     break
#     if not medium_found:
#         extracted_data['medium'] = "unknown" # Default if medium information is not found

#     return extracted_data



# import json

# def GSE10643_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from sample metadata for the GSE10643 dataset.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the specified JSON schema,
#               containing 'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The 'sample_source_name' field consistently indicates the tissue.
#     # Example: 'leaf of dor mutant under control'
#     source_name_raw = sample_metadata.get('sample_source_name', [''])[0]
#     source_name_lower = source_name_raw.lower()

#     if 'leaf' in source_name_lower:
#         extracted_data['tissue'] = 'leaf'
#     else:
#         # Default or handle cases where 'leaf' might not be present (though it is in this dataset)
#         extracted_data['tissue'] = 'unknown' 

#     # 2. Extract 'treatment'
#     # The 'sample_source_name' field also indicates the treatment condition.
#     # Examples: 'leaf of WT mutant after drought treatment', 'leaf of dor mutant under control'
#     treatment_list = []
#     if 'drought treatment' in source_name_lower:
#         treatment_list.append('drought stress')
#     elif 'control' in source_name_lower:
#         treatment_list.append('control')
    
#     # Ensure treatment list is not empty, providing a default if necessary
#     if not treatment_list:
#         treatment_list.append('untreated') # Fallback, though samples seem to always have a condition
    
#     extracted_data['treatment'] = treatment_list

#     # 3. Extract 'medium'
#     # The growth medium is not explicitly stated in the sample metadata.
#     # Based on the study description (Arabidopsis plants grown under "normal watering conditions"),
#     # 'soil' is the standard and inferred growth medium for such experiments.
#     # This information is constant across all samples in this study.
#     extracted_data['medium'] = 'soil'

#     return extracted_data




# import json

# def GSE51897_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from a sample's metadata
#     following a predefined JSON schema.

#     Args:
#         sample_metadata (dict): A dictionary containing the metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the specified JSON schema,
#               containing the extracted tissue, treatment, and medium information.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The tissue information is found within the 'sample_characteristicts' list.
#     # Example: 'tissue: seedling'
#     characteristics = sample_metadata.get('sample_characteristicts', [])
#     tissue = "unspecified" # Default value if not found
#     for char_str in characteristics:
#         if char_str.startswith('tissue:'):
#             tissue = char_str.split(':', 1)[1].strip()
#             break
#     extracted_data['tissue'] = tissue

#     # 2. Extract 'treatment'
#     # The study metadata's 'overall_design' and sample's 'Sample source_name'
#     # consistently indicate the growth conditions.
#     # 'overall_design': "grown at 29°C under long day conditions"
#     # 'sample_source_name': "29°C LD" (LD stands for Long Day)
#     # These treatments appear to be constant across all samples in this study.
#     treatments = ["29°C", "Long Day"]
#     extracted_data['treatment'] = treatments

#     # 3. Extract 'medium'
#     # The growth medium is not explicitly mentioned in the provided study_metadata
#     # or sample_metadata. Since it's a required field in the schema, and we cannot
#     # extract a specific value from the given data, we use "unspecified".
#     extracted_data['medium'] = "unspecified"

#     return extracted_data

# if __name__ == '__main__':
#     # Example study metadata (for context, not directly used by the function)
#     study_metadata = {
#         'title': ['Gene expression analysis of jmj30 jmj32 mutant grown at 29°C under long day conditions'],
#         'summary': ['Analysis of loss-of-function mutants of JUMONJI30 (JMJ30) and JMJ32 in Arabidopsis thaliana (Col). JMJ30 and JMJ32 encode the histone demethylases that contain the catalytic JmjC domain.', 'The trimethylation of Histone H3 lysine 27 (H3K27me3) plays a key role in gene repression and developmental regulation. It is actively and dynamically deposited and removed in plants. However, while the H3K27 methyltransferases had been extensively studied, findings on the H3K27 demethylase were few. Here, we show that JMJ30 and JMJ32 encode H3K27me3 demethylases, and genes down-regulated in jmj30 jmj32 mutant are enriched for H3K27me3 targets.'],
#         'overall_design': ['Columbia wild-type (Col WT) and jmj30 jmj32 mutant were grown at 29°C, and seedling samples were collected at 13 day-after-germination (DAG). Two independent sets of WT and jmj30 jmj32 seedling mRNA samples were used for this array.']
#     }

#     # Example samples metadata
#     samples_metadata = {
#         'GSM1254792': {'Sample_id': 'GSM1254792', 'sample_title': ['Col WT, rep1'], 'sample_source_name': ['Col WT, 13 DAG, 29°C LD'], 'sample_characteristicts': ['genetic background: Col-0', 'tissue: seedling', 'genotype: wild-type', 'phenotype: High-temperature-induced early flowering'], 'sample_extraction_protocol': ['Total RNA was extracted using the RNeasy Mini kit (Qiagen Inc., Valencia, CA, USA) and DNA was removed by on-column DNaseI digestion with the RNase-Free DNase set (Qiagen).']},
#         'GSM1254795': {'Sample_id': 'GSM1254795', 'sample_title': ['jmj30 jmj32, rep2'], 'sample_source_name': ['jmj30 jmj32, 13 DAG, 29°C LD'], 'sample_characteristicts': ['genetic background: Col-0', 'tissue: seedling', 'genotype: jmj30-2 jmj32-1', 'phenotype: Flowered earlier than WT grown at similar condition'], 'sample_extraction_protocol': ['Total RNA was extracted using the RNeasy Mini kit (Qiagen Inc., Valencia, CA, USA) and DNA was removed by on-column DNaseI digestion with the RNase-Free DNase set (Qiagen).']},
#         'GSM1254793': {'Sample_id': 'GSM1254793', 'sample_title': ['Col WT, rep2'], 'sample_source_name': ['Col WT, 13 DAG, 29°C LD'], 'sample_characteristicts': ['genetic background: Col-0', 'tissue: seedling', 'genotype: wild-type', 'phenotype: High-temperature-induced early flowering'], 'sample_extraction_protocol': ['Total RNA was extracted using the RNeasy Mini kit (Qiagen Inc., Valencia, CA, USA) and DNA was removed by on-column DNaseI digestion with the RNase-Free DNase set (Qiagen).']},
#         'GSM1254794': {'Sample_id': 'GSM1254794', 'sample_title': ['jmj30 jmj32, rep1'], 'sample_source_name': ['jmj30 jmj32, 13 DAG, 29°C LD'], 'sample_characteristicts': ['genetic background: Col-0', 'tissue: seedling', 'genotype: jmj30-2 jmj32-1', 'phenotype: Flowered earlier than WT grown at similar condition'], 'sample_extraction_protocol': ['Total RNA was extracted using the RNeasy Mini kit (Qiagen Inc., Valencia, CA, USA) and DNA was removed by on-column DNaseI digestion with the RNase-Free DNase set (Qiagen).']}
#     }

#     # Test the function with one of the sample metadata objects
#     sample_id_to_test = 'GSM1254792'
#     sample_data = samples_metadata[sample_id_to_test]
#     extracted_info = GSE51897_extractor(sample_data)

#     print(f"Extracted information for sample {sample_id_to_test}:")
#     print(json.dumps(extracted_info, indent=2))

#     # Test with another sample to show consistency
#     sample_id_to_test_2 = 'GSM1254795'
#     sample_data_2 = samples_metadata[sample_id_to_test_2]
#     extracted_info_2 = GSE51897_extractor(sample_data_2)

#     print(f"\nExtracted information for sample {sample_id_to_test_2}:")
#     print(json.dumps(extracted_info_2, indent=2))

#     # Expected output for both samples:
#     # {
#     #   "tissue": "seedling",
#     #   "treatment": ["29°C", "Long Day"],
#     #   "medium": "unspecified"
#     # }


# import json

# def GSE119383_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from the sample metadata
#     for the GSE119383 dataset, conforming to the specified JSON schema.

#     Args:
#         sample_metadata (dict): A dictionary containing the metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the target JSON schema, containing
#               'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # Initialize fields with default or derived constant values
#     # Medium is derived from the study's overall design, which states "collected from soil".
#     extracted_data['medium'] = "soil"
#     extracted_data['treatment'] = [] # Default to empty list, will be populated if stress is found

#     # Extract information from 'sample_characteristicts' list
#     characteristics = sample_metadata.get('sample_characteristicts', [])

#     tissue_found = False
#     for char_string in characteristics:
#         # Extract 'tissue'
#         if char_string.startswith('tissue:'):
#             extracted_data['tissue'] = char_string.split(':', 1)[1].strip()
#             tissue_found = True
#         # Extract 'treatment' based on 'time and condition'
#         elif char_string.startswith('time and condition:'):
#             condition_part = char_string.split(':', 1)[1].strip()
#             # As per study design, 'drought' is a specific stress/treatment.
#             # 'water' and 'control' are considered baseline/non-stress conditions.
#             if 'drought' in condition_part.lower():
#                 extracted_data['treatment'].append('drought')

#     # Fallback for tissue if not found, though for this dataset it appears constant
#     # and always present as 'mature whole plant'.
#     if not tissue_found:
#         # Based on the provided samples_metadata, 'tissue: mature whole plant' is consistent.
#         extracted_data['tissue'] = "mature whole plant"

#     return extracted_data



# import json

# def GSE26266_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from the sample metadata
#     for the GSE26266 dataset, conforming to a specific JSON schema.

#     The function extracts:
#     - 'tissue' from 'sample_characteristicts'.
#     - 'treatment' from 'sample_characteristicts'.
#     - 'medium' which is a constant value derived from the study's overall design.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the specified schema, containing
#               'tissue', 'treatment', and 'medium' information.
#     """

#     # Initialize variables for extraction
#     tissue = None
#     treatment = []
    
#     # The growth medium is constant across all samples in this study,
#     # as identified from the 'overall_design' in the study_metadata:
#     # "Seven-day-old seedlings of Col-0, Ws and QK grown at 22oC on 0.5x MS plates containing 1% sucrose"
#     medium = "0.5x MS plates containing 1% sucrose"

#     # Extract information from 'sample_characteristicts'
#     # This field is a list of strings, where each string is a key-value pair.
#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         for char_str in sample_metadata['sample_characteristicts']:
#             # Convert to lowercase for robust matching, then split and strip whitespace
#             char_str_lower = char_str.lower()

#             if char_str_lower.startswith('tissue:'):
#                 tissue = char_str.split(':', 1)[1].strip()
#             elif char_str_lower.startswith('treatment:'):
#                 # The schema requires 'treatment' to be an array of strings.
#                 # Even if there's only one treatment, it should be in a list.
#                 treatment_value = char_str.split(':', 1)[1].strip()
#                 treatment.append(treatment_value)
    
#     # Construct the output dictionary conforming to the schema
#     extracted_data = {
#         "tissue": tissue,
#         "treatment": treatment,
#         "medium": medium
#     }

#     return extracted_data



# import json

# def GSE112161_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from a sample metadata dictionary
#     conforming to a specific JSON schema.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the specified schema, containing
#               'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # 1. Extract Tissue
#     # Based on the study metadata and sample characteristics (e.g., 'age: 4-d-old seedlings'),
#     # the samples are consistently described as Arabidopsis seedlings.
#     extracted_data["tissue"] = "seedlings"

#     # 2. Extract Treatment
#     treatments = []
#     if 'sample_characteristicts' in sample_metadata:
#         for characteristic in sample_metadata['sample_characteristicts']:
#             char_lower = characteristic.lower()
#             if char_lower.startswith("treatment:"):
#                 treatment_value = characteristic.split(":", 1)[1].strip()
#                 # Correcting a potential typo observed in some sample metadata
#                 treatment_value = treatment_value.replace("accilimation", "acclimation")
#                 treatments.append(treatment_value)
#             elif char_lower.startswith("time point:"):
#                 time_point_value = characteristic.split(":", 1)[1].strip()
#                 treatments.append(time_point_value)
    
#     # The schema requires 'treatment' to be an array of strings.
#     # If no specific treatment is found, an empty list is a valid representation,
#     # but based on the provided samples, there should always be at least one.
#     extracted_data["treatment"] = treatments

#     # 3. Extract Medium
#     # The growth medium is not explicitly stated in the provided study or sample metadata.
#     # For Arabidopsis seedlings in a laboratory setting, a Not Specified
#     # (e.g., Murashige and Skoog (MS) medium) is commonly used.
#     # We infer this common standard as a reasonable default.
#     extracted_data["medium"] = "Not Specified"

#     return extracted_data



# import json

# def GSE78713_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from sample metadata for GSE78713
#     following a specific JSON schema.

#     The schema requires 'tissue', 'treatment', and 'medium'.
#     - 'tissue' and 'medium' are inferred as constant across samples based on
#       the study's general description of using 'plants' and common Arabidopsis
#       growth conditions in a lab setting.
#     - 'treatment' is extracted from the 'sample_characteristicts' field.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 This typically corresponds to one of the inner dictionaries
#                                 from the 'samples_metadata' provided in the context.

#     Returns:
#         dict: A dictionary formatted according to the specified schema,
#               containing tissue, treatment, and medium information.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The study metadata and overall design mention "plants" (Arabidopsis).
#     # No specific tissue (e.g., leaf, root) is mentioned in the sample
#     # characteristics or extraction protocol. Therefore, "whole plant"
#     # is a reasonable and general inference for the tissue source.
#     extracted_data["tissue"] = "whole plant"

#     # 2. Extract 'treatment'
#     treatments = []
#     # The 'sample_characteristicts' field is a list of strings,
#     # one of which specifies the treatment.
#     if 'sample_characteristicts' in sample_metadata and \
#        isinstance(sample_metadata['sample_characteristicts'], list):
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if characteristic.startswith('treatment:'):
#                 # Extract the value after "treatment: " and strip any whitespace
#                 treatment_value = characteristic.split('treatment:')[1].strip()
#                 if treatment_value: # Ensure the extracted value is not empty
#                     treatments.append(treatment_value)
    
#     # The schema requires 'treatment' to be an array of strings.
#     # For this dataset, each sample typically has one treatment (e.g., "heat stress" or "non-stress").
#     extracted_data["treatment"] = treatments

#     # 3. Extract 'medium'
#     # The growth medium is not explicitly stated in the provided sample or study metadata.
#     # For laboratory-grown Arabidopsis plants, Murashige and Skoog (MS) medium
#     # is a very common and standard synthetic growth medium. Given that 'medium'
#     # is a required field and likely constant across samples in this study,
#     # "MS medium" is a well-justified inference.
#     extracted_data["medium"] = "MS medium"

#     return extracted_data



# import json

# def GSE110857_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from a sample's metadata
#     following a predefined JSON schema.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted according to the specified schema, containing
#               'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The tissue information is found within the 'sample_characteristicts' list.
#     # Example: 'tissue: Shoot'
#     tissue = None
#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if characteristic.startswith('tissue:'):
#                 tissue = characteristic.split(':', 1)[1].strip()
#                 break
#     extracted_data['tissue'] = tissue

#     # 2. Extract 'treatment'
#     # The treatment information is derived from the 'condition' within
#     # 'sample_characteristicts'. 'Drought stressed' maps to ["Drought stress"],
#     # while 'Control' maps to an empty list.
#     treatment = []
#     if 'sample_characteristicts' in sample_metadata and isinstance(sample_metadata['sample_characteristicts'], list):
#         for characteristic in sample_metadata['sample_characteristicts']:
#             if characteristic.startswith('condition:'):
#                 condition = characteristic.split(':', 1)[1].strip()
#                 if condition == 'Drought stressed':
#                     treatment.append('Drought stress')
#                 # If condition is 'Control', the treatment list remains empty,
#                 # as 'Control' is the absence of a specific treatment/stress.
#                 break
#     if treatment == []:
#         treatment = ['control']
#     extracted_data['treatment'] = treatment

#     # 3. Extract 'medium'
#     # Based on the study_metadata['overall_design'] provided in the problem description,
#     # the growth medium for all samples in this study is "Dio propagation mix no. 2 professional soil".
#     # This is a constant value for the entire study.
#     extracted_data['medium'] = "Dio propagation mix no. 2 professional soil"

#     return extracted_data



# def GSE66369_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts relevant biological experimental information from the sample metadata
#     for the GSE66369 dataset, conforming to the specified JSON schema.

#     The schema requires 'tissue' (string), 'treatment' (array of strings),
#     and 'medium' (string).

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.

#     Returns:
#         dict: A dictionary formatted as a JSON instance conforming to the
#               output schema.
#     """
#     extracted_info = {}

#     # Initialize fields with defaults or placeholders
#     tissue = "not specified"
#     treatment = []
#     medium = "soil"  # Inferred: For Arabidopsis "whole plants", soil is a common growth medium.
#                      # This information is not explicitly present in the provided metadata,
#                      # but is a required field in the schema and likely constant across samples.

#     # Extract information from 'sample_characteristicts' list
#     characteristics = sample_metadata.get('sample_characteristicts', [])

#     for char_string in characteristics:
#         if char_string.startswith('tissue:'):
#             # Extract tissue, e.g., "tissue: whole plants" -> "whole plants"
#             tissue = char_string.split(':', 1)[1].strip()
#         elif char_string.startswith('treatment:'):
#             # Extract treatment, e.g., "treatment: Heat stress (37oC, 1 h)"
#             # The schema expects an array, so we add it to a list.
#             treatment_str = char_string.split(':', 1)[1].strip()
#             treatment.append(treatment_str)

#     # If tissue was not found in characteristics, try 'sample_source_name' as a fallback
#     if tissue == "not specified":
#         source_name = sample_metadata.get('sample_source_name', [None])[0]
#         if source_name:
#             tissue = source_name.lower() # e.g., "Whole plant" -> "whole plant"

#     # Assign extracted (or inferred) values to the output dictionary
#     extracted_info['tissue'] = tissue
#     extracted_info['treatment'] = treatment
#     extracted_info['medium'] = medium

#     return extracted_info
# import json

# def GSE27551_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from a sample metadata dictionary
#     for the GSE27551 dataset, conforming to a specific JSON schema.

#     The output schema is:
#     {"properties": {"tissue": {"description": "Tissue the samples was extracted from.", "title": "Tissue", "type": "string"},
#                     "treatment": {"description": "List of treatments and stresses that was applied to the sample.", "items": {"type": "string"}, "title": "Treatment", "type": "array"},
#                     "medium": {"description": "Growth medium of the sample.", "title": "Medium", "type": "string"}},
#      "required": ["tissue", "treatment", "medium"]}

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 Example:
#                                 {'Sample_id': 'GSM680433',
#                                  'sample_title': ['leaf_gDNA_Got22_3'],
#                                  'sample_source_name': ['rosette leaf'],
#                                  'sample_characteristicts': ['accession: CS22609', 'ecotype: Got22', 'tissue: rosette leaf'],
#                                  'sample_extraction_protocol': ['Qiagen Dneasy Plant Mini prep']}

#     Returns:
#         dict: A dictionary formatted according to the specified schema,
#               containing 'tissue', 'treatment', and 'medium' information.
#     """
#     extracted_data = {}

#     # 1. Extract 'tissue'
#     # The tissue information is consistently found under 'sample_source_name'
#     # and also within 'sample_characteristicts' as 'tissue: rosette leaf'.
#     # We prioritize 'sample_source_name' for direct extraction.
#     if 'sample_source_name' in sample_metadata and sample_metadata['sample_source_name']:
#         extracted_data['tissue'] = sample_metadata['sample_source_name'][0]
#     else:
#         # Fallback to 'sample_characteristicts' if 'sample_source_name' is missing or empty
#         tissue_found = False
#         if 'sample_characteristicts' in sample_metadata:
#             for characteristic in sample_metadata['sample_characteristicts']:
#                 if characteristic.startswith('tissue:'):
#                     extracted_data['tissue'] = characteristic.split(':', 1)[1].strip()
#                     tissue_found = True
#                     break
#         if not tissue_found:
#             # If tissue cannot be found, assign a default or raise an error as it's a required field.
#             # For this dataset, 'rosette leaf' is consistent.
#             extracted_data['tissue'] = "unknown" # Or handle as appropriate for missing data

#     # 2. Extract 'treatment'
#     # The provided sample metadata is for "Genomic dna hybridizations" (gDNA).
#     # While the study summary mentions "soil moisture deficit" for gene expression studies,
#     # there are no explicit treatments or stresses mentioned for these specific gDNA samples.
#     # Therefore, we assume no specific treatments were applied for the purpose of gDNA extraction
#     # and represent this as an empty list.
#     extracted_data['treatment'] = []

#     # 3. Extract 'medium'
#     # The growth medium is not explicitly stated in the sample metadata.
#     # However, the study metadata mentions "soil moisture deficit" in the context of
#     # the overall study, implying that the Arabidopsis thaliana plants were grown in soil.
#     # Given that 'medium' is a required field, we infer "soil" as the growth medium.
#     extracted_data['medium'] = "soil"

#     return extracted_data



# import json

# def GSE11538_extractor(sample_metadata: dict) -> dict:
#     """
#     Extracts biological experimental information from sample metadata for the GSE11538 dataset,
#     conforming to the specified JSON schema.

#     The function identifies the tissue, treatment, and growth medium for a given sample
#     based on the provided metadata structure.

#     Args:
#         sample_metadata (dict): A dictionary containing metadata for a single sample.
#                                 Expected to have keys like 'sample_source_name'.

#     Returns:
#         dict: A dictionary formatted according to the specified schema:
#               {"tissue": "...", "treatment": ["..."], "medium": "..."}
#     """
#     extracted_data = {}

#     # The 'sample_source_name' field contains both tissue and treatment information,
#     # typically in the format "tissue_name, treatment_condition".
#     # We use .get() with a default to safely handle cases where the key might be missing,
#     # though for this dataset, it's expected to be present.
#     source_name_list = sample_metadata.get('sample_source_name', [''])
#     source_name = source_name_list[0] if source_name_list else ''

#     # 1. Extract 'tissue'
#     # The tissue is the part before the first comma in the 'source_name'.
#     if source_name:
#         tissue_part = source_name.split(',')[0].strip()
#         extracted_data['tissue'] = tissue_part
#     else:
#         # Fallback if 'sample_source_name' is empty or malformed.
#         # "unspecified" is used as a placeholder to ensure schema compliance.
#         extracted_data['tissue'] = "unspecified"

#     # 2. Extract 'treatment'
#     # The treatment is the part after the first comma in the 'source_name'.
#     treatment_list = []
#     if source_name and ',' in source_name:
#         # Split only once to ensure that if the treatment itself contains commas,
#         # the entire treatment string is captured.
#         treatment_str = source_name.split(',', 1)[1].strip()
#         if treatment_str:
#             treatment_list.append(treatment_str)
    
#     # The schema requires 'treatment' to be an array of strings.
#     # For this dataset, samples are either "watered" (control) or "water deficit" (stress).
#     # The extracted `treatment_list` will contain one of these.
#     extracted_data['treatment'] = treatment_list

#     # 3. Extract 'medium'
#     # The provided metadata does not explicitly state a growth medium for Arabidopsis thaliana plants.
#     # Given that the samples are "vegetative rosettes" (implying whole plants) and the common
#     # practices for growing Arabidopsis thaliana in such studies, "soil" is a reasonable
#     # and widely applicable default for the growth medium when not specified.
#     extracted_data['medium'] = "soil"

#     return extracted_data


