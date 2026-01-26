import os
#NEED TO CHANGE
GLOBAL_DIR_PATH = F'{os.getcwd()}/'
CLUSTER_RUN = True
EXPERIMENT_NAME = '3.2'

if CLUSTER_RUN:
    STORAGE_DIR ='/tudelft.net/staff-umbrella/AT GE Datasets/'
else:
    STORAGE_DIR = f'./new_storage/'


DATA_IMPORT_DIR = './data/downloads_new/'
GEO_DOWNLOAD_DIR = f'{DATA_IMPORT_DIR}geo_downloads_new/'
METADATA_OUTPUT_DIR = f'{DATA_IMPORT_DIR}metadata_new/'
PROCESSED_DATA_FOLDER = f'new_storage/final_data/'

if CLUSTER_RUN:
    PROCESSED_DATA_FOLDER = '/tudelft.net/staff-umbrella/AT GE Datasets/final_data/'

COMBINED_DATA_OUTPUT_FILE = f'{PROCESSED_DATA_FOLDER}RMA_Microarray_Combined.csv'
CORE_DATA_DIR = './data/core_data/'
SOFT_PATH = f'{CORE_DATA_DIR}old_geo_downloads/'
FIGURES_DIR = f'./outputs/{EXPERIMENT_NAME}/'
FILTERING_FIGURES = f'{FIGURES_DIR}filter_figures/{EXPERIMENT_NAME}/'
CLUSTER_EXPLORATION_FIGURES_DIR=f'{FIGURES_DIR}exploration_figures/{EXPERIMENT_NAME}/'
MODEL = 'extractors_and_gemini'
LABELS_PATH = f'{STORAGE_DIR}labels/{MODEL}/{EXPERIMENT_NAME}'
Studies = ['GSE44053', 'GSE77815', 'GSE16474', 'GSE4062', 'GSE9415', 'GSE18624', 'GSE40061', 'GSE112161', 'GSE20494', 'GSE79997', 'GSE27552', 'GSE16222', 'GSE110857', 'GSE4760', 'GSE46205', 'GSE71001', 'GSE58616', 'GSE162310', 'GSE22107', 'GSE62163', 'GSE51897', 'GSE72949', 'GSE90562', 'GSE5628', 'GSE26266', 'GSE34188', 'GSE34595', 'GSE76827', 'GSE119383', 'GSE65046', 'GSE11758', 'GSE65414', 'GSE37408', 'GSE5624', 'GSE10643', 'GSE15577', 'GSE11538', 'GSE70861', 'GSE6583', 'GSE27551', 'GSE12619', 'GSE121225', 'GSE108070', 'GSE78713', 'GSE110079', 'GSE63128', 'GSE60960', 'GSE37118', 'GSE79681', 'GSE63372', 'GSE5622', 'GSE26983', 'GSE27550', 'GSE19603', 'GSE95202', 'GSE53308', 'GSE16765', 'GSE71855', 'GSE58620', 'GSE24177', 'GSE35258', 'GSE10670', 'GSE49418', 'GSE18666', 'GSE83136', 'GSE44655', 'GSE27549', 'GSE19700', 'GSE103398', 'GSE63522', 'GSE201609', 'GSE5620', 'GSE66369', 'GSE2268', 'GSE71237', 'GSE48474', 'GSE41935', 'GSE27548', 'GSE5623', 'GSE72050', 'GSE126373']

# LABELS
from enum import Enum

# 1. Define the Enums (The Source of Truth)
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
    SEED = "seed"
    SEEDLING = "seedling"
    # PROTOPLASTS = "protoplasts"
    # GUARD_CELLS = "guard_cells"
    POLLEN = "pollen"
    CELL_CULTURE = "cell_culture"
    UNKNOWN = "unknown"

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
    OTHER = "Other stress"
    NONE = "No stress"
class TreatmentEnum_alt(str, Enum):
    DROUGHT = "Drought"
    DEHYDRATION = "Dehydration"
    SALINITY = "Salinity"
    HEAT = "Heat"
    COLD = "Cold"
    CHEMICAL = "Chemical"
    NUTRIENT = "Nutrient Deficiency"
    BIOTIC = "Biotic"
    ABIOTIC = "Abiotic"
    LOW_LIGHT = "Low Light"
    HIGH_LIGHT = "High Light"
    OTHER_LIGHT = "Other Light"
    OTHER = "Other"
    NONE = "Control"
class MediumEnum(str, Enum):
    MS = "MS medium"
    B5 = "Gamborg B5 medium"
    SOIL = "Soil"
    VERMICULITE = "Vermiculite"
    PERLITE = "Perlite"
    SAND = "Sand"
    HYDROPONIC = "Hydroponic"
    LIQUID = "Liquid culture"
    AGAR = "Agar plate"
    UNSPECIFIED = "Unspecified"

#TODO: add devStage, mutant, and cell type

# 2. Maintain your lists for Grounding/Vectors (Backward Compatibility)
# This extracts the values automatically, so you don't need to type them twice.
VALID_TISSUES = [t.value for t in TissueEnum]
VALID_TREATMENTS = [t.value for t in TreatmentEnum]
VALID_TREATMENTS_ALT= [t.value for t in TreatmentEnum_alt]
VALID_MEDIUMS = [m.value for m in MediumEnum]