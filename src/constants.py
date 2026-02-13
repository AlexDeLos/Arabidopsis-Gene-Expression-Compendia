import os
from typing import Dict
#NEED TO CHANGE
GLOBAL_DIR_PATH = F'{os.getcwd()}/'
CLUSTER_RUN = GLOBAL_DIR_PATH !='/home/alex/Documents/GitHub/Dataset_fusion_Microarray/'
EXPERIMENT_NAME = '4.0'

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
CLUSTER_EXPLORATION_FIGURES_DIR=f'{FIGURES_DIR}exploration_figures/{EXPERIMENT_NAME}/'
MODEL = 'extractors_and_gemini'
LABELS_PATH = f'{STORAGE_DIR}labels/{MODEL}/{EXPERIMENT_NAME}'
Studies = ['GSE44053', 'GSE77815', 'GSE16474', 'GSE4062', 'GSE9415', 'GSE18624', 'GSE40061', 'GSE112161', 'GSE20494', 'GSE79997', 'GSE27552', 'GSE16222', 'GSE110857', 'GSE4760', 'GSE46205', 'GSE71001', 'GSE58616', 'GSE162310', 'GSE22107', 'GSE62163', 'GSE51897', 'GSE72949', 'GSE90562', 'GSE5628', 'GSE26266', 'GSE34188', 'GSE34595', 'GSE76827', 'GSE119383', 'GSE65046', 'GSE11758', 'GSE65414', 'GSE37408', 'GSE5624', 'GSE10643', 'GSE15577', 'GSE11538', 'GSE70861', 'GSE6583', 'GSE27551', 'GSE12619', 'GSE121225', 'GSE108070', 'GSE78713', 'GSE110079', 'GSE63128', 'GSE60960', 'GSE37118', 'GSE79681', 'GSE63372', 'GSE5622', 'GSE26983', 'GSE27550', 'GSE19603', 'GSE95202', 'GSE53308', 'GSE16765', 'GSE71855', 'GSE58620', 'GSE24177', 'GSE35258', 'GSE10670', 'GSE49418', 'GSE18666', 'GSE83136', 'GSE44655', 'GSE27549', 'GSE19700', 'GSE103398', 'GSE63522', 'GSE201609', 'GSE5620', 'GSE66369', 'GSE2268', 'GSE71237', 'GSE48474', 'GSE41935', 'GSE27548', 'GSE5623', 'GSE72050', 'GSE126373']

# Trackers

STATUS_NOT_TRIED = 0
STATUS_DOWNLOADED = 1
STATUS_PROCESSED = 2
STATUS_IGNORE = 3
STATUS_ERROR = 4


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
    WHOLE_PLANT = "whole_plant"
    SILIQUE = "silique"
    CALLUS = "callus"
    SEEDLING = "seedling"
    SEED = "seed"
    STEM = "stem"
    POLLEN = "pollen"
    CELL_CULTURE = "cell_culture"
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
    OTHER = "Other stress"
    NONE = "Control"
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

# 3. Define Synonyms (Map Canonical Enum -> List of Synonyms)
# This replaces TreatmentEnum_alt. You can add as many variations as you want here.
TREATMENT_SYNONYMS = {
    TreatmentEnum.DROUGHT: ["Drought", "Water Deficit", "Water withholding"],
    TreatmentEnum.DEHYDRATION: ["Dehydration"],
    TreatmentEnum.SALINITY: ["Salinity", "Salt", "NaCl"],
    TreatmentEnum.HEAT: ["Heat", "High Temperature"],
    TreatmentEnum.COLD: ["Cold", "Low Temperature", "Freezing", "Chilling"],
    TreatmentEnum.LOW_LIGHT: ["Dark", "Shade", "Low intensity light"],
    TreatmentEnum.HIGH_LIGHT: ["High Light", "High intensity light"],
    TreatmentEnum.OTHER_LIGHT: ["Light", "Light quality"],
    TreatmentEnum.NONE: ["No stress", "Control", "Mock"],
    # Add others as needed...
}

TISSUE_SYNONYMS = {
    TissueEnum.ROOT: ["Roots", "Root system", "Radicle"],
    TissueEnum.LEAF: ["Leaves", "Foliage", "Cotyledon"], # Cotyledon often grouped with leaf
    # Add others as needed...
}

MEDIUM_SYNONYMS = {
    MediumEnum.MS: ["Murashige and Skoog", "MS salts", "1/2 MS"],
    MediumEnum.SOIL: ["Potting mix", "Earth", "Compost"],
    # Add others as needed...
}

# 4. Master Configuration (The "Registry")
# To add a new label type (e.g. 'genotype'), add it here with its Enum and Synonyms.
LABEL_CONFIG = {
    'treatment': {
        'enum': TreatmentEnum,
        'synonyms': TREATMENT_SYNONYMS
    },
    'tissue': {
        'enum': TissueEnum,
        'synonyms': TISSUE_SYNONYMS
    },
    'medium': {
        'enum': MediumEnum,
        'synonyms': MEDIUM_SYNONYMS
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

# 6. Legacy / Area Keywords (Kept from your original)
AREA_KEYWORDS = {
    'treatment': ['treatment', 'treated', 'stress', 'condition', 'exposed to', 'exposure', 'incubated', 'temperature', 'growth condition','temperature', 'oC'],
    'tissue': ['tissue', 'organ', 'source', 'derived from', 'cells', 'cell type', 'organism part'],
    'medium': ['medium', 'growth medium', 'grown on', 'cultured in', 'substrate']
}