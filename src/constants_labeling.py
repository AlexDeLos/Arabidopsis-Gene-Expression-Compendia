# LABELS
from enum import Enum
from typing import Dict, List, Any

# 1. Define the Labels (The keys for everything)
LABELS = ['tissue', 'treatment', 'medium', 'genotype']

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

class GenotypeEnum(str, Enum):
    # Wild-type accessions
    COL_0       = "Col-0"
    WS          = "Ws"
    WS_2        = "Ws-2"
    WS_4        = "Ws-4"
    LER         = "Ler"
    C24         = "C24"
    CVI         = "Cvi"
    # Mutant / transgenic classes
    KNOCKOUT    = "Knockout mutant"
    KNOCKDOWN   = "Knockdown mutant"
    OVEREXPRESSOR = "Overexpressor"
    REPORTER    = "Reporter line"
    RNAi        = "RNAi line"
    CRISPR      = "CRISPR mutant"
    # Catch-alls
    WILD_TYPE   = "Wild-type"
    MUTANT      = "Mutant"
    TRANSGENIC  = "Transgenic"
    UNKNOWN     = "unknown"
    UNSPECIFIED = "unspecified"

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

GENOTYPE_SYNONYMS = {
    # Wild-type accessions
    GenotypeEnum.COL_0: [
        "Columbia", "Col", "Col-0", "Columbia-0", "ecotype Col",
        "background Col-0", "WT Col-0",
    ],
    GenotypeEnum.WS: [
        "Wassilewskija", "Ws", "WS",
    ],
    GenotypeEnum.WS_2: ["Ws-2", "WS-2"],
    GenotypeEnum.WS_4: ["Ws-4", "WS-4"],
    GenotypeEnum.LER: [
        "Landsberg erecta", "Ler", "Ler-0", "L. er",
    ],
    GenotypeEnum.C24: ["C24", "ecotype C24"],
    GenotypeEnum.CVI: ["Cape Verde Islands", "Cvi", "Cvi-0"],

    # Functional mutant / transgenic classes
    GenotypeEnum.KNOCKOUT: [
        "T-DNA insertion", "loss-of-function", "null mutant",
        "insertion line", "SALK line", "SAIL line", "GABI-Kat",
    ],
    GenotypeEnum.KNOCKDOWN: [
        "hypomorphic", "partial loss-of-function", "weak allele",
    ],
    GenotypeEnum.OVEREXPRESSOR: [
        "overexpression", "35S::", "35S promoter", "OE", "OX",
        "gain-of-function", "constitutive expression",
    ],
    GenotypeEnum.REPORTER: [
        "GFP fusion", "GUS fusion", "luciferase", "reporter construct",
        "promoter:GFP", "promoter:GUS",
    ],
    GenotypeEnum.RNAi: [
        "RNAi", "RNA interference", "hairpin", "amiRNA",
        "artificial microRNA", "gene silencing",
    ],
    GenotypeEnum.CRISPR: [
        "CRISPR", "CRISPR-Cas9", "Cas9", "genome editing", "edited line",
    ],

    # Catch-alls (used when no specific accession or class is determinable)
    GenotypeEnum.WILD_TYPE: [
        "wild type", "wild-type", "WT", "wt", "unmodified", "non-transgenic",
    ],
    GenotypeEnum.MUTANT: [
        "mutant", "loss of function", "allele", "atxxx mutant",
    ],
    GenotypeEnum.TRANSGENIC: [
        "transgenic", "transformed", "stably transformed",
    ],
}
class IntensityEnum(int, Enum):
    CONTROL = 0
    MILD = 1
    MODERATE = 2
    SEVERE = 3

INTENSITY_DESCRIPTIONS = {
    IntensityEnum.CONTROL: "Control / No stress / Mock treatment",
    IntensityEnum.MILD: "Mild stress (e.g., slight temperature change, low concentration)",
    IntensityEnum.MODERATE: "Moderate stress (e.g., standard stress assay conditions)",
    IntensityEnum.SEVERE: "Severe/Extreme stress (e.g., lethal temperatures, high concentration, prolonged exposure)"
}

# 4. Master Configuration (The "Registry")
LABEL_CONFIG = {
    'treatment': {
        'enum': TreatmentEnum,
        'synonyms': TREATMENT_SYNONYMS,
        'search_triggers': ['treatment', 'treated', 'stress', 'condition', 'exposure'],
        'priority_cols': ['title','characteristics_ch1', 'treatment_protocol_ch1','treatment'],
        'sub_attributes': {
            'intensity': {
                'instruction': "For every treatment extracted, you must assign an intensity score based on the text.",
                'descriptions': INTENSITY_DESCRIPTIONS
            }
        }
    },
    'tissue': {
        'enum': TissueEnum,
        'synonyms': TISSUE_SYNONYMS,
        'search_triggers': ['tissue', 'organ', 'source', 'derived from','cell'],
        'priority_cols': ['title','source_name_ch1', 'characteristics_ch1','tissue']
    },
    'medium': {
        'enum': MediumEnum,
        'synonyms': MEDIUM_SYNONYMS,
        'search_triggers': ['medium', 'growth medium', 'substrate'],
        'priority_cols': ['titel','growth_protocol_ch1', 'characteristics_ch1','medium']
    },
    'genotype': {
        'enum': GenotypeEnum,
        'synonyms': GENOTYPE_SYNONYMS,
        # "ecotype", "background", "genotype" are the most reliable column triggers.
        # "line" is intentionally excluded — too generic and appears in unrelated contexts.
        'search_triggers': ['genotype', 'ecotype', 'background', 'accession', 'strain', 'mutation'],
        'priority_cols': ['characteristics_ch1', 'title', 'source_name_ch1', 'genotype', 'ecotype']
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
UNIQUE_LABELS = ['tissue', 'medium', 'genotype']

# Map each category to its Control value. 
# Assuming your enums are imported here, use `.value` to get the string.
CONTROL_MAP = {
    'treatment': TreatmentEnum.CONTROL.value, # e.g., 'Control'
    # Add others if needed: 'tissue': TissueEnum.CONTROL.value
}