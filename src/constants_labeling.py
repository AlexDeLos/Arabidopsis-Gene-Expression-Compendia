# LABELS
from enum import Enum

# 1. Define the Labels (The keys for everything)
# LABELS = ['tissue', 'treatment', 'medium', 'genotype']
LABELS = ["tissue", "treatment", "medium", "ecotype", "modification", "developmental_stage"]


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
	# SEEDLING = "seedling"
	SEED = "seed"
	STEM = "stem"
	POLLEN = "pollen"
	CELL_CULTURE = "cell culture"
	UNKNOWN = "unknown"
	UNSPECIFIED = "unspecified"


class DevelopmentalStageEnum(str, Enum):
	GERMINATION = "Germination"
	SEEDLING = "Seedling"
	VEGETATIVE = "Vegetative"
	BOLTING = "Bolting"
	FLOWERING = "Flowering"
	FRUITING = "Fruiting"
	SENESCENCE = "Senescence"
	UNKNOWN = "unknown"
	UNSPECIFIED = "unspecified"


class TreatmentEnum(str, Enum):
	DROUGHT = "Drought"
	FLOOD = "Flood"
	DEHYDRATION = "Dehydration"
	SALINITY = "Salinity"
	HEAT = "Heat"
	COLD = "Cold"
	CHEMICAL = "Chemical"
	BIOTIC = "Biotic"
	ABIOTIC = "Abiotic"
	LOW_LIGHT = "Low Light"
	HIGH_LIGHT = "High Light"
	OTHER_LIGHT = "Other Light"
	CUT = "Cut"
	NUTRIENT_DEFICIENCY = "Nutrient Deficiency"
	OTHER = "Other"
	CONTROL = "Control"
	UNKNOWN = "unknown"
	UNSPECIFIED = "unspecified"


class MediumEnum(str, Enum):
	MS = "MS"
	GAMBORG_B5 = "Gamborg B5 medium"
	SOIL = "Soil"
	VERMICULITE = "Vermiculite"
	PERLITE = "Perlite"
	SAND = "Sand"
	HYDROPONIC = "Hydroponic"
	LIQUID = "Liquid"
	AGAR = "Agar"
	UNKNOWN = "unknown"
	UNSPECIFIED = "unspecified"


class EcotypeEnum(str, Enum):
	COL_0 = "Col-0"
	WS = "Ws"
	WS_2 = "Ws-2"
	WS_4 = "Ws-4"
	LER = "Ler"
	C24 = "C24"
	CVI = "Cvi"
	GENERIC_WILD_TYPE = "Wild-type (Unspecified ecotype)"


class ModificationEnum(str, Enum):
	NONE = "None (Wild-type)"
	KNOCKOUT = "Knockout mutant"
	KNOCKDOWN = "Knockdown mutant"
	OVEREXPRESSOR = "Overexpressor"
	REPORTER = "Reporter line"
	RNAI = "RNAi line"
	CRISPR = "CRISPR mutant"
	GENERIC_MUTANT = "Mutant (Unspecified type)"
	GENERIC_TRANSGENIC = "Transgenic (Unspecified type)"


# 3. Define Synonyms (Map Canonical Enum -> List of Synonyms)
# This replaces TreatmentEnum_alt. You can add as many variations as you want here.
TREATMENT_SYNONYMS = {
	# --- Hydration and Osmotic Stresses ---
	TreatmentEnum.DROUGHT: ["Drought", "Water Deficit", "Water withholding", "Water deprivation", "Drought stress", "Dry soil"],
	TreatmentEnum.DEHYDRATION: ["Dehydration", "Desiccation", "Dry air", "Air dry"],
	TreatmentEnum.SALINITY: ["Salinity", "Salt", "NaCl", "Sodium chloride", "CaCl2", "KCl", "Salt stress", "Osmotic stress"],
	# --- Temperature Stresses ---
	TreatmentEnum.HEAT: ["Heat", "High Temperature", "Heat shock", "Elevated temperature", "Warm", "high ambient temperature"],
	TreatmentEnum.COLD: ["Cold", "Low Temperature", "Freezing", "Chilling", "Frost", "Cold stress", "Ice", "Frozen"],
	# --- Light and Radiation ---
	TreatmentEnum.LOW_LIGHT: ["Dark", "Shade", "Low intensity light", "Darkness", "Etiolated", "Far-red light", "Green light", "Short day"],
	TreatmentEnum.HIGH_LIGHT: ["High Light", "High intensity light", "UV-B", "UV-A", "Ultraviolet", "Photoinhibition", "Long day"],
	TreatmentEnum.OTHER_LIGHT: ["Light", "Light quality", "Continuous light", "White light", "Red light", "Blue light", "Photoperiod"],
	# --- Physical / Mechanical ---
	TreatmentEnum.CUT: ["Wounded", "Wounding", "Cut", "Excised", "Detached", "Mechanical damage", "Punctured", "Laser capture microdissection"],
	# --- Nutrient Stresses ---
	TreatmentEnum.NUTRIENT_DEFICIENCY: ["Nutrient deficiency", "Iron deficiency", "Minus iron", "Nitrogen starvation", "Phosphate deficiency", "Starvation", "Deprivation", "Absence of"],
	# --- Baseline ---
	TreatmentEnum.CONTROL: ["No stress", "Control", "Mock", "Untreated", "Normal conditions", "Standard media", "Ambient", "Vehicle control", "Room temperature"],
}

TISSUE_SYNONYMS = {
	TissueEnum.ROOT: ["Roots", "Root system", "Radicle", "Root tip", "Lateral root", "Primary root", "Root hair", "Root meristem"],
	TissueEnum.LEAF: ["Leaves", "Foliage", "Cotyledon", "Leaf blade", "Leaflet", "Rosette leaf", "Cauline leaf", "True leaves", "Leaf primordia"],
	# TissueEnum.SEEDLING: ["Seedlings", "Plantlet", "Sprout", "Young plant", "Etiolated seedling"],
	TissueEnum.POLLEN: ["Pollen grains", "Pollen tube", "Microspore", "Pollen grain"],
	TissueEnum.FLOWER: ["Flowers", "Inflorescence", "Floral", "Petal", "Sepal", "Stamen", "Carpel", "Pistil", "Anther", "Stigma", "Ovary", "Ovule"],
	TissueEnum.CALLUS: ["Calli", "Callus culture", "Epidermis", "Epidermal cells"],
	TissueEnum.SEED: ["Seeds", "Seed coat", "Endosperm", "Embryo", "Germinating seed", "Dry seed", "Imbibed seed"],
	TissueEnum.SHOOT: ["Shoots", "Shoot apex", "Shoot apical meristem", "Aerial parts", "Aerial tissue", "Aboveground parts"],
	TissueEnum.ROSETTE: ["Rosettes", "Vegetative rosette", "Whole rosette", "Rosette leaves"],
	TissueEnum.CELL_CULTURE: ["Cell culture", "Suspension culture", "Protoplast", "Tissue culture", "Cultured cells", "Liquid culture cells"],
	TissueEnum.BUD: ["Buds", "Floral bud", "Flower bud", "Apical bud"],
	TissueEnum.STEM: ["Stems", "Hypocotyl", "Stalk", "Inflorescence stem", "Meristem", "Epicotyl", "Internode", "Shoot axis"],
	TissueEnum.SILIQUE: ["Siliques", "Pod", "Fruit", "Seedpod", "Valve", "Dehiscence zone"],
}

MEDIUM_SYNONYMS = {
	# --- Synthetic Nutrient Media ---
	MediumEnum.MS: ["Murashige and Skoog", "MS", "MS salts", "MS medium", "1/2 MS", "Half-strength MS", "0.5X MS", "MS plates", "MS agar"],
	# --- Physical Substrates (Solid) ---
	MediumEnum.SOIL: ["Soil", "Potting mix", "Earth", "Compost", "Potting soil", "Peat", "Peat moss", "Levington compost"],
	MediumEnum.VERMICULITE: ["Vermiculite", "Mica"],
	MediumEnum.PERLITE: ["Perlite", "Volcanic glass substrate"],
	# --- Physical States / Gelling Agents ---
	MediumEnum.AGAR: ["Agar", "Agar plates", "Solid medium", "Solid media", "Phytagel", "Gelrite", "Agarose", "Gelled medium"],
	MediumEnum.LIQUID: ["Liquid", "Liquid culture", "Liquid medium", "Liquid broth", "Suspension medium", "Liquid MS", "Liquid MS medium"],
	# --- Specialized Growth Systems ---
	MediumEnum.HYDROPONIC: ["Hydroponic", "Hydroponics", "Liquid nutrient solution", "Hydroponic culture", "Aerated liquid culture"],
}

DEVELOPMENTAL_SYNONYMS = {
	DevelopmentalStageEnum.GERMINATION: ["germinating", "imbibed", "stratified", "stage 0"],
	DevelopmentalStageEnum.SEEDLING: ["plantlet", "young plant", "days post germination", "dpg", "stage 1"],
	DevelopmentalStageEnum.VEGETATIVE: ["rosette stage", "leaf production", "pre-flowering", "stage 3"],
	DevelopmentalStageEnum.BOLTING: ["stem emergence", "inflorescence emergence", "stage 5"],
	DevelopmentalStageEnum.FLOWERING: ["anthesis", "blooming", "floral", "stage 6"],
	DevelopmentalStageEnum.FRUITING: ["silique development", "seed filling", "pod development", "stage 8"],
	DevelopmentalStageEnum.SENESCENCE: ["aging", "drying", "yellowing", "terminal", "stage 9"],
	DevelopmentalStageEnum.UNKNOWN: [],
	DevelopmentalStageEnum.UNSPECIFIED: [],
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
	IntensityEnum.SEVERE: "Severe/Extreme stress (e.g., lethal temperatures, high concentration, prolonged exposure)",
}
DEVELOPMENTAL_STAGE_DESCRIPTIONS = {
	DevelopmentalStageEnum.GERMINATION: "The process beginning with dry seed imbibition and ending with radicle emergence.",
	DevelopmentalStageEnum.SEEDLING: "The early growth phase showing expanded cotyledons, prior to true rosette leaf expansion.",
	DevelopmentalStageEnum.VEGETATIVE: "The main growth phase characterized by true rosette leaf production before the transition to flowering.",
	DevelopmentalStageEnum.BOLTING: "The phase of rapid elongation of the primary inflorescence stem.",
	DevelopmentalStageEnum.FLOWERING: "The reproductive period when flowers are open and anthesis occurs.",
	DevelopmentalStageEnum.FRUITING: "The post-fertilization stage focusing on silique expansion and seed development.",
	DevelopmentalStageEnum.SENESCENCE: "The terminal aging phase marked by chlorophyll breakdown, yellowing, and tissue death.",
}

TISSUE_DESCRIPTIONS = {
	TissueEnum.ROOT: "The below-ground portion of the plant responsible for water and nutrient uptake.",
	TissueEnum.LEAF: "The primary photosynthetic organ, including cotyledons, rosette leaves, and cauline leaves.",
	TissueEnum.FLOWER: "The reproductive structure, including petals, sepals, stamens, and carpels.",
	TissueEnum.SHOOT: "The entire above-ground portion of the plant, including stems, leaves, and reproductive organs.",
	TissueEnum.ROSETTE: "The circular, basal arrangement of leaves that forms before the plant bolts.",
	TissueEnum.BUD: "An undeveloped or embryonic shoot, often referring to floral or apical meristems.",
	TissueEnum.WHOLE_PLANT: "The entire organism, including both aerial and subterranean parts.",
	TissueEnum.SILIQUE: "The seed-bearing fruit capsule typical of Arabidopsis and other Brassicaceae.",
	TissueEnum.CALLUS: "An unorganized mass of undifferentiated, actively dividing cells grown in vitro.",
	# TissueEnum.SEEDLING: "A very young plant newly emerged from a seed, typically encompassing the cotyledon stage.",
	TissueEnum.SEED: "The mature fertilized ovule containing the embryo, endosperm, and seed coat.",
	TissueEnum.STEM: "The main structural axis of the plant, including the hypocotyl and the inflorescence stalk.",
	TissueEnum.POLLEN: "The male microgametophytes produced in the anther.",
	TissueEnum.CELL_CULTURE: "Cells grown in liquid suspension or on solid media, often as protoplasts or undifferentiated cell lines.",
}
TREATMENT_DESCRIPTIONS = {
	TreatmentEnum.DROUGHT: "Insufficient water availability in the growth medium (e.g., withholding water from soil).",
	TreatmentEnum.FLOOD: "Excessive water application leading to complete or partial submergence and root hypoxia.",
	TreatmentEnum.DEHYDRATION: "Removal of water from the plant tissue itself, often by placing excised plants in dry air.",
	TreatmentEnum.SALINITY: "Exposure to high concentrations of salts, typically NaCl, causing osmotic and ionic stress.",
	TreatmentEnum.HEAT: "Exposure to temperatures significantly above the optimal growth range (e.g., >30°C).",
	TreatmentEnum.COLD: "Exposure to low, non-freezing chilling temperatures or freezing temperatures.",
	TreatmentEnum.CHEMICAL: "Exposure to exogenous chemicals, hormones, toxins, gases (e.g., ozone, hypoxia), or chemical elicitors.",
	TreatmentEnum.BIOTIC: "Interaction with living organisms, including pathogens, herbivores, or application of pathogen-associated molecular patterns (e.g., flg22).",
	TreatmentEnum.ABIOTIC: "A generic designation for a non-living environmental stressor when the specific type (e.g., heat, cold) is unspecified.",
	TreatmentEnum.LOW_LIGHT: "Exposure to darkness, shading, or suboptimal light intensity.",
	TreatmentEnum.HIGH_LIGHT: "Exposure to excessively high light intensity or damaging UV radiation.",
	TreatmentEnum.OTHER_LIGHT: "Specific light quality treatments (e.g., red, blue, far-red) or photoperiod changes.",
	TreatmentEnum.CUT: "Mechanical damage, wounding, or physical excision of plant tissues.",
	TreatmentEnum.NUTRIENT_DEFICIENCY: "Growth in a medium lacking one or more essential macro- or micro-nutrients (e.g., -N, -P, -Fe).",
	TreatmentEnum.OTHER: "A specified treatment or stress condition that does not fit into any of the standard canonical categories.",
	TreatmentEnum.CONTROL: "Standard, optimal, unperturbed growth conditions. Used as the baseline for comparison.",
}

MEDIUM_DESCRIPTIONS = {
	MediumEnum.MS: "Murashige and Skoog medium, a standard synthetic nutrient mixture for plant tissue culture.",
	MediumEnum.GAMBORG_B5: "Gamborg's B5 medium, a specific synthetic basal salt mixture for in vitro culture.",
	MediumEnum.SOIL: "Natural or commercial potting mixes consisting of organic and inorganic matter.",
	MediumEnum.VERMICULITE: "A hydrous phyllosilicate mineral substrate used for aeration and moisture retention.",
	MediumEnum.PERLITE: "An amorphous volcanic glass substrate used to improve aeration and drainage.",
	MediumEnum.SAND: "Granular material composed of finely divided rock and mineral particles.",
	MediumEnum.HYDROPONIC: "Liquid nutrient solution setups where roots are completely submerged or bathed in liquid without soil.",
	MediumEnum.LIQUID: "General liquid broth or suspension medium without a gelling agent.",
	MediumEnum.AGAR: "Any solid or semi-solid medium gelled with agar, agarose, or phytagel.",
}

ECOTYPE_DESCRIPTIONS = {
	EcotypeEnum.COL_0: "Columbia-0 (Col-0), the standard reference accession.",
	EcotypeEnum.WS: "Wassilewskija (Ws) general ecotype.",
	EcotypeEnum.WS_2: "Wassilewskija-2 (Ws-2) sub-line.",
	EcotypeEnum.WS_4: "Wassilewskija-4 (Ws-4) sub-line.",
	EcotypeEnum.LER: "Landsberg erecta (Ler) ecotype.",
	EcotypeEnum.C24: "C24 ecotype.",
	EcotypeEnum.CVI: "Cape Verde Islands (Cvi) ecotype.",
	EcotypeEnum.GENERIC_WILD_TYPE: "Used when the text mentions 'wild-type' but does not specify which ecotype (e.g., Col-0, Ler) was used.",
}

MODIFICATION_DESCRIPTIONS = {
	ModificationEnum.NONE: "The plant has no genetic modifications; it is a true wild-type baseline.",
	ModificationEnum.KNOCKOUT: "A complete loss-of-function mutation, usually via T-DNA insertion.",
	ModificationEnum.KNOCKDOWN: "Reduced, but not eliminated, gene expression (e.g., weak hypomorphic allele).",
	ModificationEnum.OVEREXPRESSOR: "Engineered to express a gene at high levels (e.g., using a 35S promoter).",
	ModificationEnum.REPORTER: "Expressing a visual/quantifiable marker (e.g., GFP, GUS, Luciferase).",
	ModificationEnum.RNAI: "Utilizing RNA interference or artificial microRNAs to silence transcripts.",
	ModificationEnum.CRISPR: "A targeted mutation or edit generated using CRISPR/Cas technology.",
	ModificationEnum.GENERIC_MUTANT: "A plant with a genetic alteration, but the exact mechanism (e.g., knockout, point mutation) is not stated.",
	ModificationEnum.GENERIC_TRANSGENIC: "A plant carrying introduced foreign DNA, but the specific function (e.g., reporter, overexpressor) is not stated.",
}
# 4. Master Configuration (The "Registry")
LABEL_CONFIG = {
	"treatment": {
		"enum": TreatmentEnum,
		"synonyms": TREATMENT_SYNONYMS,
		"descriptions": TREATMENT_DESCRIPTIONS,
		"search_triggers": ["treatment", "treated", "stress", "condition", "exposure"],
		"priority_cols": ["title", "characteristics_ch1", "treatment_protocol_ch1", "treatment"],
		"sub_attributes": {
			"intensity": {"enum": IntensityEnum, "instruction": "For every treatment extracted, you must assign an intensity score based on the text.", "descriptions": INTENSITY_DESCRIPTIONS}
		},
	},
	"tissue": {
		"enum": TissueEnum,
		"synonyms": TISSUE_SYNONYMS,
		"descriptions": TISSUE_DESCRIPTIONS,
		"search_triggers": ["tissue", "organ", "source", "derived from", "cell"],
		"priority_cols": ["title", "source_name_ch1", "characteristics_ch1", "tissue"],
	},
	"medium": {
		"enum": MediumEnum,
		"synonyms": MEDIUM_SYNONYMS,
		"descriptions": MEDIUM_DESCRIPTIONS,
		"search_triggers": ["medium", "growth medium", "substrate"],
		"priority_cols": ["titel", "growth_protocol_ch1", "characteristics_ch1", "medium"],
	},
	"ecotype": {
		"enum": EcotypeEnum,
		"synonyms": {},
		"descriptions": ECOTYPE_DESCRIPTIONS,
		"search_triggers": ["genotype", "ecotype", "background", "accession", "strain", "mutation"],
		"priority_cols": ["characteristics_ch1", "title", "source_name_ch1", "genotype", "ecotype"],
	},
	"modification": {
		"enum": EcotypeEnum,
		"synonyms": {},
		"descriptions": MODIFICATION_DESCRIPTIONS,
		"search_triggers": ["genotype", "ecotype", "background", "accession", "strain", "mutation"],
		"priority_cols": ["characteristics_ch1", "title", "source_name_ch1", "genotype", "ecotype"],
	},
	"developmental_stage": {
		"enum": DevelopmentalStageEnum,
		"search_triggers": ["stage", "development", "age", "days old", "weeks old", "dpg", "das", "boyes"],
		"synonyms": DEVELOPMENTAL_SYNONYMS,
		"descriptions": DEVELOPMENTAL_STAGE_DESCRIPTIONS,
	},
}

# 5. Auto-Generate Dictionaries
# BUCKET_KEYWORDS: Only the Canonical values (for Vector embedding base)
# EXPLICIT_KEYWORDS: Canonical + All Synonyms (for Search/Extraction)
# CANONICAL_MAP: Synonym String -> Canonical String (for Grounding/Collapsing)

BUCKET_KEYWORDS: dict[str, list[str]] = {}
EXPLICIT_KEYWORDS: dict[str, list[str]] = {}
CANONICAL_MAP: dict[str, dict[str, str]] = {}
AREA_KEYWORDS: dict[str, list[str]] = {}

for label in LABELS:
	config = LABEL_CONFIG.get(label)
	if not config:
		continue

	enum_cls = config["enum"]
	synonym_dict = config["synonyms"]
	# sub_attributes = config['sub_attributes']

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
		CANONICAL_MAP[label][canonical_val] = canonical_val  # Identity map

		# 3. Add Synonyms to Explicit & Map
		if item in synonym_dict:
			for syn in synonym_dict[item]:
				if syn not in EXPLICIT_KEYWORDS[label]:  # Avoid dupes
					EXPLICIT_KEYWORDS[label].append(syn)
				CANONICAL_MAP[label][syn] = canonical_val


# Auto-generate trigger mapping for the extractor
ALL_TRIGGERS = {label: config["search_triggers"] + EXPLICIT_KEYWORDS[label] for label, config in LABEL_CONFIG.items()}

# Define which categories should strictly have only ONE label
# UNIQUE_LABELS = ['tissue', 'medium', 'genotype']
UNIQUE_LABELS = ["tissue", "medium", "genotype", "developmental_stage"]

# Map each category to its Control value.
# Assuming your enums are imported here, use `.value` to get the string.
CONTROL_MAP = {
	"treatment": TreatmentEnum.CONTROL.value,  # e.g., 'Control'
	# Add others if needed: 'tissue': TissueEnum.CONTROL.value
}


STRESS_GO_ROOTS: dict[str, tuple[str, str]] = {
	# treatment_value: (go_id, go_name)
	"Drought":			  ("GO:0009414", "response to water deprivation"),
	"Flood":				("GO:0071456", "cellular response to hypoxia"),
	"Dehydration":		   ("GO:0009414", "response to water deprivation"),
	"Salinity":			 ("GO:0009651", "response to salt stress"),
	"Heat":				 ("GO:0009408", "response to heat"),
	"Cold":				 ("GO:0009409", "response to cold"),
	"Chemical":			 ("GO:0042221", "response to chemical"),
	"Biotic":			   ("GO:0009607", "response to biotic stimulus"),
	"Abiotic":			  ("GO:0009628", "response to abiotic stimulus"),
	"Low Light":			("GO:0009642", "response to light intensity"),
	"High Light":		   ("GO:0009644", "response to high light intensity"),
	"Other Light":		  ("GO:0009416", "response to light stimulus"),
	"Cut":				  ("GO:0009611", "response to wounding"),
	"Nutrient Deficiency":  ("GO:0031667", "response to nutrient levels"),
}