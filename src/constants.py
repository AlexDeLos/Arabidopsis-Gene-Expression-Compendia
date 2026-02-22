import os
import pandas as pd
from typing import Dict
#NEED TO CHANGE
GLOBAL_DIR_PATH = F'{os.getcwd()}/'
CLUSTER_RUN = GLOBAL_DIR_PATH !='/home/alex/Documents/GitHub/Dataset_fusion_Microarray/'
EXPERIMENT_NAME = '4.6'

if CLUSTER_RUN:
    STORAGE_DIR ='/tudelft.net/staff-umbrella/GeneExpressionStorage/'
else:
    STORAGE_DIR = f'./new_storage/'


DATA_IMPORT_DIR = './data/downloads_new/'
GEO_DOWNLOAD_DIR = f'{DATA_IMPORT_DIR}geo_downloads_new/'
METADATA_OUTPUT_DIR = f'{DATA_IMPORT_DIR}metadata_new/'
PROCESSED_DATA_FOLDER = f'new_storage/final_data/'

if CLUSTER_RUN:
    PROCESSED_DATA_FOLDER = '/tudelft.net/staff-umbrella/GeneExpressionStorage/final_data/'

COMBINED_DATA_OUTPUT_FILE = f'{PROCESSED_DATA_FOLDER}RMA_Microarray_Combined.csv'
CORE_DATA_DIR = './data/core_data/'
SOFT_PATH = f'{CORE_DATA_DIR}old_geo_downloads/'
FIGURES_DIR = f'./outputs/{EXPERIMENT_NAME}/'
FILTERING_FIGURES = f'{FIGURES_DIR}filter_figures/{EXPERIMENT_NAME}/'
CLUSTER_EXPLORATION_FIGURES_DIR=f'{FIGURES_DIR}exploration_figures/'
MODEL = 'extractor_grouding'
LABELS_PATH = f'{STORAGE_DIR}labels/{MODEL}/{EXPERIMENT_NAME}'
Studies = ['GSE44053', 'GSE77815', 'GSE16474', 'GSE4062', 'GSE9415', 'GSE18624', 'GSE40061', 'GSE112161', 'GSE20494', 'GSE79997', 'GSE27552', 'GSE16222', 'GSE110857', 'GSE4760', 'GSE46205', 'GSE71001', 'GSE58616', 'GSE162310', 'GSE22107', 'GSE62163', 'GSE51897', 'GSE72949', 'GSE90562', 'GSE5628', 'GSE26266', 'GSE34188', 'GSE34595', 'GSE76827', 'GSE119383', 'GSE65046', 'GSE11758', 'GSE65414', 'GSE37408', 'GSE5624', 'GSE10643', 'GSE15577', 'GSE11538', 'GSE70861', 'GSE6583', 'GSE27551', 'GSE12619', 'GSE121225', 'GSE108070', 'GSE78713', 'GSE110079', 'GSE63128', 'GSE60960', 'GSE37118', 'GSE79681', 'GSE63372', 'GSE5622', 'GSE26983', 'GSE27550', 'GSE19603', 'GSE95202', 'GSE53308', 'GSE16765', 'GSE71855', 'GSE58620', 'GSE24177', 'GSE35258', 'GSE10670', 'GSE49418', 'GSE18666', 'GSE83136', 'GSE44655', 'GSE27549', 'GSE19700', 'GSE103398', 'GSE63522', 'GSE201609', 'GSE5620', 'GSE66369', 'GSE2268', 'GSE71237', 'GSE48474', 'GSE41935', 'GSE27548', 'GSE5623', 'GSE72050', 'GSE126373']

# Trackers

STATUS_NOT_TRIED = 0
STATUS_DOWNLOADED = 1
STATUS_PROCESSED = 2
STATUS_IGNORE = 3
STATUS_ERROR = 4

try:
    SAMPLE_STUDY_MAP = pd.read_csv(STORAGE_DIR+'/final_data/RMA_Microarray_Combined_sample_map.csv', index_col=0) # todo: change back
except FileNotFoundError:
    SAMPLE_STUDY_MAP = None

# LABELS
from enum import Enum
from typing import Dict, List, Any

# 1. Define the Labels (The keys for everything)
LABELS = ['tissue', 'treatment', 'medium']

# 2. Define the Enums (The Canonical "Truth")
class TissueEnum(str, Enum):
    ROOT = "root"
    LEAF = "leaf"
    FLOWER = "flower"
    SHOOT = "shoot"
    ROSETTE = "rosette"
    BUD = "bud"
    WHOLE_PLANT = "whole plant"
    SILIQUE = "silique"
    CALLUS = "callus"
    SEEDLING = "seedling"
    SEED = "seed"
    STEM = "stem"
    POLLEN = "pollen"
    CELL_CULTURE = "cell culture"
    UNKNOWN = "unknown"
    UNSPECIFIED = "unspecified"

class TreatmentEnum(str, Enum):
    DROUGHT = "Drought Stress"
    DEHYDRATION = "Dehydration Stress"
    SALINITY = "Salinity Stress"
    HEAT = "Heat Stress"
    COLD = "Cold Stress"
    CHEMICAL = "Chemical Stress"
    NUTRIENT = "Nutrient Deficiency"
    BIOTIC = "Biotic Stress"
    ABIOTIC = "Abiotic Stress"
    LOW_LIGHT = "Low Light Stress"
    HIGH_LIGHT = "High Light Stress"
    OTHER_LIGHT = "Other Light Stress"
    CUT = "Cut"
    NUTRIENT_DEFICIENCY = "Nutrient Deficiency Stress"
    CHEMICAL_STRESS = "Chemical Stress"
    BIOTIC_STRESS = "Biotic Stress"
    OTHER = "Other stress"
    CONTROL = "Control"
    UNKNOWN = "unknown"
    UNSPECIFIED = "unspecified"

class MediumEnum(str, Enum):
    MS = "MS"
    B5 = "Gamborg B5 medium"
    SOIL = "Soil"
    VERMICULITE = "Vermiculite"
    PERLITE = "Perlite"
    SAND = "Sand"
    HYDROPONIC = "Hydroponic"
    LIQUID = "Liquid"
    AGAR = "Agar"
    UNKNOWN = "unknown"
    UNSPECIFIED = "unspecified"
    GAMBORG_B5_MEDIUM = "Gamborg B5"

# 3. Define Synonyms (Map Canonical Enum -> List of Synonyms)
# This replaces TreatmentEnum_alt. You can add as many variations as you want here.
TREATMENT_SYNONYMS = {
    # --- Hydration and Osmotic Stresses ---
    TreatmentEnum.DROUGHT: [
        "Drought", "Water Deficit", "Water withholding", "Water deprivation", 
        "Drought stress", "Dry soil"
    ],
    TreatmentEnum.DEHYDRATION: [
        "Dehydration", "Desiccation", "Dry air", "Air dry"
    ],
    TreatmentEnum.SALINITY: [
        "Salinity", "Salt", "NaCl", "Sodium chloride", "CaCl2", "KCl", 
        "Salt stress", "Osmotic stress"
    ],

    # --- Temperature Stresses ---
    TreatmentEnum.HEAT: [
        "Heat", "High Temperature", "Heat shock", "Elevated temperature", "Warm"
    ],
    TreatmentEnum.COLD: [
        "Cold", "Low Temperature", "Freezing", "Chilling", "Frost", 
        "Cold stress", "Ice", "Frozen"
    ],

    # --- Light and Radiation ---
    TreatmentEnum.LOW_LIGHT: [
        "Dark", "Shade", "Low intensity light", "Darkness", "Etiolated", 
        "Far-red light", "Green light"
    ],
    TreatmentEnum.HIGH_LIGHT: [
        "High Light", "High intensity light", "UV-B", "UV-A", "Ultraviolet", 
        "Photoinhibition", "Long day"
    ],
    TreatmentEnum.OTHER_LIGHT: [
        "Light", "Light quality", "Continuous light", "White light", 
        "Red light", "Blue light", "Photoperiod"
    ],

    # --- Physical / Mechanical ---
    TreatmentEnum.CUT: [
        "Wounded", "Wounding", "Cut", "Excised", "Detached", 
        "Mechanical damage", "Punctured", "Laser capture microdissection"
    ],

    # --- Nutrient Stresses ---
    TreatmentEnum.NUTRIENT_DEFICIENCY: [
        "Nutrient deficiency", "Iron deficiency", "Minus iron", "Nitrogen starvation", 
        "Phosphate deficiency", "Starvation", "Deprivation", "Absence of"
    ],

    # --- Chemical / Atmospheric ---
    TreatmentEnum.CHEMICAL_STRESS: [
        "Hypoxia", "Anoxia", "Low oxygen", "Dexamethasone", "Ozone", 
        "Estradiol", "Ethanol", "Methanol", "Fumigation", "Chemical treatment"
    ],

    # --- Biotic ---
    TreatmentEnum.BIOTIC_STRESS: [
        "Inoculated", "Infected", "Agrobacterium", "Pseudomonas syringae", 
        "Pathogen", "Fungus", "Elicitor", "Chitin", "Flg22"
    ],

    # --- Baseline ---
    TreatmentEnum.CONTROL: [
        "No stress", "Control", "Mock", "Untreated", "Normal conditions", 
        "Standard media", "Ambient", "Vehicle control"
    ]
}

TISSUE_SYNONYMS = {
    TissueEnum.ROOT: [
        "Roots", "Root system", "Radicle", "Root tip", "Lateral root", 
        "Primary root", "Root hair", "Root meristem"
    ],
    TissueEnum.LEAF: [
        "Leaves", "Foliage", "Cotyledon", "Leaf blade", "Leaflet", 
        "Rosette leaf", "Cauline leaf", "True leaves", "Leaf primordia"
    ],
    TissueEnum.SEEDLING: [
        "Seedlings", "Plantlet", "Sprout", "Young plant", "Etiolated seedling"
    ],
    TissueEnum.POLLEN: [
        "Pollen grains", "Pollen tube", "Microspore", "Pollen grain"
    ],
    TissueEnum.FLOWER: [
        "Flowers", "Inflorescence", "Floral", "Petal", "Sepal", 
        "Stamen", "Carpel", "Pistil", "Anther", "Stigma", "Ovary", "Ovule"
    ],
    TissueEnum.CALLUS: [
        "Calli", "Callus culture", "Epidermis", "Epidermal cells"
    ],
    TissueEnum.SEED: [
        "Seeds", "Seed coat", "Endosperm", "Embryo", "Germinating seed", 
        "Dry seed", "Imbibed seed"
    ],
    TissueEnum.SHOOT: [
        "Shoots", "Shoot apex", "Shoot apical meristem", "Aerial parts", 
        "Aerial tissue", "Aboveground parts"
    ],
    TissueEnum.ROSETTE: [
        "Rosettes", "Vegetative rosette", "Whole rosette", "Rosette leaves"
    ],
    TissueEnum.CELL_CULTURE: [
        "Cell culture", "Suspension culture", "Protoplast", "Tissue culture", 
        "Cultured cells", "Liquid culture cells"
    ],
    TissueEnum.BUD: [
        "Buds", "Floral bud", "Flower bud", "Apical bud"
    ],
    TissueEnum.STEM: [
        "Stems", "Hypocotyl", "Stalk", "Inflorescence stem", "Meristem", 
        "Epicotyl", "Internode", "Shoot axis"
    ],
    TissueEnum.SILIQUE: [
        "Siliques", "Pod", "Fruit", "Seedpod", "Valve", "Dehiscence zone"
    ]
}

MEDIUM_SYNONYMS = {
    # --- Synthetic Nutrient Media ---
    MediumEnum.MS: [
        "Murashige and Skoog", "MS", "MS salts", "MS medium", 
        "1/2 MS", "Half-strength MS", "0.5X MS", "MS plates", "MS agar"
    ],
    MediumEnum.GAMBORG_B5_MEDIUM: [
        "Gamborg's", "Gamborg B5", "B5 medium", "Gamborg's B5", "B5 salts"
    ],
    
    # --- Physical Substrates (Solid) ---
    MediumEnum.SOIL: [
        "Soil", "Potting mix", "Earth", "Compost", "Potting soil", 
        "Peat", "Peat moss", "Levington compost"
    ],
    MediumEnum.VERMICULITE: [
        "Vermiculite", "Mica"
    ],
    MediumEnum.PERLITE: [
        "Perlite", "Volcanic glass substrate"
    ],

    # --- Physical States / Gelling Agents ---
    MediumEnum.AGAR: [
        "Agar", "Agar plates", "Solid medium", "Solid media", 
        "Phytagel", "Gelrite", "Agarose", "Gelled medium"
    ],
    MediumEnum.LIQUID: [
        "Liquid", "Liquid culture", "Liquid medium", "Liquid broth", 
        "Suspension medium", "Liquid MS", "Liquid MS medium"
    ],
    
    # --- Specialized Growth Systems ---
    MediumEnum.HYDROPONIC: [
        "Hydroponic", "Hydroponics", "Liquid nutrient solution", 
        "Hydroponic culture", "Aerated liquid culture"
    ]
}

# 4. Master Configuration (The "Registry")
LABEL_CONFIG = {
    'treatment': {
        'enum': TreatmentEnum,
        'synonyms': TREATMENT_SYNONYMS,
        'search_triggers': ['treatment', 'treated', 'stress', 'condition', 'exposure'],
        'priority_cols': ['titel','characteristics_ch1', 'treatment_protocol_ch1','treatment']
    },
    'tissue': {
        'enum': TissueEnum,
        'synonyms': TISSUE_SYNONYMS,
        'search_triggers': ['tissue', 'organ', 'source', 'derived from'],
        'priority_cols': ['titel','source_name_ch1', 'characteristics_ch1','tissue']
    },
    'medium': {
        'enum': MediumEnum,
        'synonyms': MEDIUM_SYNONYMS,
        'search_triggers': ['medium', 'growth medium', 'substrate'],
        'priority_cols': ['titel','growth_protocol_ch1', 'characteristics_ch1','medium']
    }
}

# 5. Auto-Generate Dictionaries
# BUCKET_KEYWORDS: Only the Canonical values (for Vector embedding base)
# EXPLICIT_KEYWORDS: Canonical + All Synonyms (for Search/Extraction)
# CANONICAL_MAP: Synonym String -> Canonical String (for Grounding/Collapsing)

BUCKET_KEYWORDS: Dict[str, List[str]] = {}
EXPLICIT_KEYWORDS: Dict[str, List[str]] = {}
CANONICAL_MAP: Dict[str, Dict[str, str]] = {}
AREA_KEYWORDS: Dict[str, List[str]] = {} 

for label in LABELS:
    config = LABEL_CONFIG.get(label)
    if not config:
        continue

    enum_cls = config['enum']
    synonym_dict = config['synonyms']

    # Initialize lists
    BUCKET_KEYWORDS[label] = []
    EXPLICIT_KEYWORDS[label] = []
    CANONICAL_MAP[label] = {}

    for item in enum_cls:
        canonical_val = item.value
        
        # 1. Add to Bucket (Target)
        BUCKET_KEYWORDS[label].append(canonical_val)
        
        # 2. Add Canonical to Explicit & Map
        EXPLICIT_KEYWORDS[label].append(canonical_val)
        CANONICAL_MAP[label][canonical_val] = canonical_val # Identity map
        
        # 3. Add Synonyms to Explicit & Map
        if item in synonym_dict:
            for syn in synonym_dict[item]:
                if syn not in EXPLICIT_KEYWORDS[label]: # Avoid dupes
                    EXPLICIT_KEYWORDS[label].append(syn)
                CANONICAL_MAP[label][syn] = canonical_val


# Auto-generate trigger mapping for the extractor
ALL_TRIGGERS = {label: config['search_triggers'] + EXPLICIT_KEYWORDS[label] 
                for label, config in LABEL_CONFIG.items()}

# Define which categories should strictly have only ONE label
UNIQUE_LABELS = ['tissue', 'medium'] 

# Map each category to its Control value. 
# Assuming your enums are imported here, use `.value` to get the string.
CONTROL_MAP = {
    'treatment': TreatmentEnum.CONTROL.value, # e.g., 'Control'
    # Add others if needed: 'tissue': TissueEnum.CONTROL.value
}