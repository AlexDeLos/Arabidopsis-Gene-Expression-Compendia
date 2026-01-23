from langchain_core.prompts import ChatPromptTemplate
from langchain_core.messages import SystemMessage
from langchain_core.output_parsers import PydanticOutputParser, StrOutputParser
from langchain.chat_models import init_chat_model
from typing import List, Literal
from pydantic import ValidationError, BaseModel, Field
import sys
module_dir = './'
sys.path.append(module_dir)
from src.constants import *

class MetaData(BaseModel):
    tissue: TissueEnum = Field(..., description="Tissue the samples was extracted from.") # "stem", 
    # tissue: str = Field(..., description="Tissue the samples was extracted from.") # "stem", 

    # Field with description explaining its purpose
    treatment: List[TreatmentEnum] = Field(..., description="List of treatments and stresses that was applied to the sample, each unique stress or treatment should have one, and only one, entry in this list")
    
    # Field with description explaining its purpose
    medium: MediumEnum = Field(..., description="Growth medium of the sample.")




def get_metadata(study_info:dict,sample_info:dict,model:str='gemini-2.5-flash',temp:float=0):
    llm = init_chat_model(model=model,
                        model_provider="google_genai",
                        temperature=temp)
    parser = PydanticOutputParser(pydantic_object=MetaData)
    system_message = SystemMessage(content="You are a biology reseracher that has been tasked with creating a script for analysing the metadata of given samples that follows a particular sctructure and returning information on the biological experimental conditions of the sample that are requested. You want to stick to labels that are both relevant but wide enouogh that other studies can fall under the same label. If more than one stress is applied return both.")


    format_instructions = parser.get_format_instructions()

    prompt = ChatPromptTemplate.from_messages([
        ("system", "Extract per schema:\n{format_instructions}"),
        ("human", "{text}"),
    ]).partial(format_instructions=format_instructions)

    parsing_llm = prompt | llm | parser


    parsing_prompt = '''
    <task>
    Extract the relevant biological experimental information of the sample following the given schema:
    </task>
    <metadata>
        <study_metadata>
            {study_info}
        </study_metadata>
        <sample_metadata>
            {sample_info}
        </sample_metadata>
    </metadata>
    '''    
    try:
        result = parsing_llm.invoke({"text": parsing_prompt.format(study_info=study_info,sample_info=sample_info)}) # type: ignore #! may also raise an uncaght JSONDecodeError ->GSM1126269
    except Exception as e:
        print(f'trying a better model-> {e}')
        if model=='gemini-2.5-pro':
            raise e#-> GSM1126262 ->GSM843683
        return get_metadata(study_info,sample_info, model='gemini-2.5-pro', temp=0.1)
    return result

def get_condensed_labels(study_info:dict,sample_info:dict,model:str='gemini-2.5-flash',temp:float=0):
    llm = init_chat_model(model=model,
                        model_provider="google_genai",
                        temperature=temp)
    parser = PydanticOutputParser(pydantic_object=MetaData)


    format_instructions = parser.get_format_instructions()

    prompt = ChatPromptTemplate.from_messages([
        ("system", "You are a biology reseracher that has been tasked with creating a script for analysing the metadata of given samples that follows a particular sctructure and returning information on the biological experimental conditions of the sample that are requested. You want to stick to labels that are both relevant but wide enouogh that other studies can fall under the same label. If more than one stress is applied label it as 'control'. Keep the labels in the following schema:\n{format_instructions}"),
        ("human", "{text}"),
    ]).partial(format_instructions=format_instructions)

    parsing_llm = prompt | llm | parser


    parsing_prompt = '''
    <task>
    Using the study information and the extracted sample labels ground the labels to comonly used term matching the ontology:
    </task>
    <metadata>
        <study_metadata>
            {study_info}
        </study_metadata>
        <samples_label>
            {sample_info}
        </samples_label>
    </metadata>
    '''    
    try:
        result = parsing_llm.invoke({"text": parsing_prompt.format(study_info=study_info,sample_info=sample_info)}) # type: ignore #! may also raise an uncaght JSONDecodeError ->GSM1126269
    except Exception as e:
        print(f'trying a better model-> {e}')
        if model=='gemini-2.5-pro':
            raise e
        return get_metadata(study_info,sample_info, model='gemini-2.5-pro', temp=0.1)
    return result
#### NEW

import re

def clean_metadata_for_context(metadata: dict) -> dict:
    """
    Intelligently prunes nested metadata (Study & Sample) to strictly biological context.
    1. Digs into 'study_metadata' and 'sample_metadata'.
    2. Prioritizes known 'Gold' fields (characteristics, source).
    3. Filters 'Silver' fields (protocols, description) by biological keywords.
    4. Removes technical 'Noise' fields.
    """
    
    final_cleaned_output = {}

    # --- 1. KEYWORD LISTS ---
    BIO_KEYWORDS = [
        'medium', 'agar', 'soil', 'grown', 'growth', 'culture', 'treated', 'treatment',
        'stress', 'incubated', 'light', 'dark', 'cycles', 'conditions', 'tissue', 
        'genotype', 'mutant', 'wild type', 'col-0', 'sample', 'seedling', 'rosette', 
        'leaf', 'root', 'mm', 'celsius', 'hours', 'days', 'h', 'd'
    ]+ VALID_MEDIUMS + VALID_TISSUES + VALID_TREATMENTS
    
    # Keys we ALWAYS keep full (High Density)
    # Fixed typo: 'overall_desing' -> 'overall_design'
    keep_keys = ['title', 'source_name_ch1', 'characteristics_ch1', 'overall_design', 'summary']
    
    # Keys we ALWAYS discard (Technical Noise)
    drop_keys = ['extract', 'hyb', 'scan', 'label', 'data_processing', 'contact', 'taxid', 'platform_id', 'id', 'platform', 'series_id', 'relation']

    # --- 2. HELPER: SENTENCE FILTERING ---
    def filter_text_block(text_list_or_str):
        """
        Takes a list or string, splits into sentences, and only returns 
        sentences containing BIO_KEYWORDS.
        """
        if isinstance(text_list_or_str, list):
            # Join lists (like ["line 1", "line 2"]) into a single block
            text = " ".join([str(x) for x in text_list_or_str])
        else:
            text = str(text_list_or_str)

        # Split by periods to get rough sentences (handles "E. coli" or "Ph.D." edge cases roughly)
        sentences = re.split(r'(?<!\w\.\w.)(?<![A-Z][a-z]\.)(?<=\.|\?)\s', text)
        
        relevant_sentences = []
        for sent in sentences:
            if any(k in sent.lower() for k in BIO_KEYWORDS):
                relevant_sentences.append(sent)
        
        if not relevant_sentences and len(text) < 100:
             return text
        elif not relevant_sentences:
             return text[:200] + "... (truncated)"
             
        return " ".join(relevant_sentences)

    # --- 3. HELPER: DICTIONARY CLEANER ---
    def _clean_sub_dict(sub_dict: dict) -> dict:
        """Applies the filtering logic to a specific sub-dictionary."""
        cleaned_sub = {}
        for key, value in sub_dict.items():
            key_lower = key.lower()
            
            # A. DROP Technical Fields
            if any(bad in key_lower for bad in drop_keys):
                continue

            # B. KEEP "Gold" fields full
            if any(good in key_lower for good in keep_keys):
                cleaned_sub[key] = value
                continue

            # C. PROCESS "Silver" fields (Descriptions, Protocols)
            if any(x in key_lower for x in ['protocol', 'description', 'design', 'notes']):
                filtered_val = filter_text_block(value)
                if filtered_val and len(filtered_val) > 5:
                    cleaned_sub[key] = filtered_val
        return cleaned_sub

    # --- 4. MAIN EXECUTION: Dig into the specific sections ---
    
    # A. Process Study Metadata (General context: Growth conditions, Summary)
    if 'study_metadata' in metadata and isinstance(metadata['study_metadata'], dict):
        final_cleaned_output['study_metadata'] = _clean_sub_dict(metadata['study_metadata'])

    # B. Process Sample Metadata (Specific context: Tissue, Treatment, Genotype)
    if 'sample_metadata' in metadata and isinstance(metadata['sample_metadata'], dict):
        final_cleaned_output['sample_metadata'] = _clean_sub_dict(metadata['sample_metadata'])

    # C. Fallback: If input was flat (no study/sample split), clean the root
    if not final_cleaned_output and metadata:
         final_cleaned_output = _clean_sub_dict(metadata)

    return final_cleaned_output

def get_metadata_script(sample_metadata: dict, study_id: str, model: str = 'gemini-2.5-flash', temp: float = 0):
    """
    Generates a Python extractor function for a specific study based on a sample's metadata.
    """
    
    # --- OPTIMIZATION START ---
    # We create a "lean" version of the metadata just for the LLM prompt.
    # The generated script will still run on the FULL metadata in production,
    # but the LLM only needs to see the structure of the relevant fields to write the code.
    context_metadata = clean_metadata_for_context(sample_metadata)#TODO: improve this
    # context_metadata = sample_metadata
    # --- OPTIMIZATION END ---

    llm = init_chat_model(model=model, model_provider="google_genai", temperature=temp)
    parser = PydanticOutputParser(pydantic_object=MetaData)
    format_instructions = parser.get_format_instructions()

    prompt = ChatPromptTemplate.from_messages([
        ("system", "Produce a python code function that can extract the following schema from the data dictionary object provided by the user. \nReturn ONLY the python code.\n{format_instructions}"),
        ("human", "{text}"),
    ]).partial(format_instructions=format_instructions)

    parsing_llm = prompt | llm | StrOutputParser()

    parsing_prompt = '''
    <task>
    Provide a python function with the signature:
    def {study_id}_extractor(sample_metadata: dict) -> dict:
    
    This function should extract the relevant biological experimental information (tissue, treatment, medium) from the input dictionary.
    </task>

    <guidance>
    1. The input `sample_metadata` provided below is a *representative subset* of the fields containing biological info.
    2. Your code should check for these specific keys (e.g., 'characteristics_ch1', 'source_name_ch1', 'title').
    3. If values are lists, access the first element or join them as appropriate.
    4. If 'medium' is not explicitly found, default to "unspecified" (or infer "soil" if context implies whole plants).
    5. Return ONLY valid Python code.
    </guidance>

    <metadata_structure>
        <sample_metadata>
            {sample_metadata}
        </sample_metadata>
    </metadata_structure>
    '''

    try:
        # Pass the CLEANED metadata to the prompt
        result = parsing_llm.invoke({
            "text": parsing_prompt.format(
                study_id=study_id, 
                sample_metadata=json.dumps(context_metadata, indent=2) # json.dumps saves tokens vs raw dict str
            )
        })
        
        result = result.replace("```python", "").replace("```", "").strip()
        
    except Exception as e:
        print(f'Error generating script with {model}: {e}')
        if 'flash' in model:
            print('Retrying with a Pro model...')
            return get_metadata_script(sample_metadata, study_id, model='gemini-1.5-pro', temp=temp)
        else:
            raise e

    return result
import json
from typing import List, Dict
from langchain_core.output_parsers import StrOutputParser
from langchain_core.prompts import ChatPromptTemplate
# Ensure init_chat_model is imported (it seems to be in your llm_utils already)

def get_batch_labels_treatment(terms: List[str],context, model: str = 'gemini-2.5-flash', temp: float = 0.0) -> Dict[str, str]:
    """
    Maps a list of raw experimental terms to standardized ontology labels using an LLM.
    
    Args:
        terms: List of unique raw strings (e.g. ['50mM NaCl', 'Mock treated'])
        model: LLM model name
        temp: Temperature (0 is best for classification)
        
    Returns:
        Dict mapping input term -> standardized label
    """
    
    # 1. Handle empty input immediately
    if not terms:
        return {}

    # 2. Define your strict ontology (Matching your MetaData class in llm_utils)
    valid_labels = VALID_TREATMENTS

    # 3. Initialize LLM
    llm = init_chat_model(model=model, model_provider="google_genai", temperature=temp)

    # 4. Construct Prompt
    # We ask for strict JSON output and provide the valid list
    prompt_template = ChatPromptTemplate.from_template("""
    You are an expert biological data curator.
    Your task is to map the provided list of raw experimental conditions/treatments to the most appropriate Standardized Label from the specific list below.

    <standard_labels>
    {ontology}
    </standard_labels>

    <rules>
    1. Output MUST be a valid JSON dictionary.
    2. Keys must be the exact strings from the Input List.
    3. Values must be one of the Standard Labels.
    4. If a term implies a control condition (e.g., "Mock", "Water", "DMSO"), map it to "No stress".
    5. If a term is ambiguous or fits none of the specific categories, use "Other stress".
    6. Do not include markdown formatting (like ```json). Just the raw JSON string.
    </rules>

    <condensed_contex>
    {context}
    </condensed_contex>
    <input_list>
    {terms}
    </input_list>
    
    JSON Output:
    """)

    # 5. Build Chain
    chain = prompt_template | llm | StrOutputParser()

    try:
        # 6. Execute
        # We join the valid labels into a string for the prompt
        ontology_str = ", ".join([f'"{x}"' for x in valid_labels])
        context_metadata = clean_metadata_for_context(context)
        response_str = chain.invoke({
            "ontology": ontology_str,
            "context": json.dumps(context_metadata),
            "terms": json.dumps(terms) # Pass terms as a JSON string representation
        })

        # 7. Clean and Parse Response
        # LLMs often wrap JSON in markdown blocks (```json ... ```), remove them.
        cleaned_response = response_str.replace("```json", "").replace("```", "").strip()
        
        result_dict = json.loads(cleaned_response)
        
        return result_dict

    except json.JSONDecodeError:
        print(f"Error: LLM returned invalid JSON for batch.")
        return {} # Return empty dict on failure so pipeline continues
    except Exception as e:
        print(f"Error in get_batch_labels: {e}")
        return {}
def get_batch_labels_tissues(terms: List[str], context, model: str = 'gemini-2.5-flash', temp: float = 0.0) -> Dict[str, str]:
    """
    Maps a list of raw tissue/organ terms to standardized ontology labels using an LLM.
    
    Args:
        terms: List of unique raw strings (e.g. ['rosette leaves', 'whole seedling'])
        context: The study/sample metadata dictionary
        model: LLM model name
        temp: Temperature (0 is best for classification)
        
    Returns:
        Dict mapping input term -> standardized label
    """
    
    # 1. Handle empty input immediately
    if not terms:
        return {}

    # 2. Define your strict ontology (Matching TissueEnum in src.constants)
    valid_tissues = VALID_TISSUES

    # 3. Initialize LLM
    llm = init_chat_model(model=model, model_provider="google_genai", temperature=temp)

    # 4. Construct Prompt
    # Adapted specifically for Tissue mapping rules
    prompt_template = ChatPromptTemplate.from_template("""
    You are an expert biological data curator.
    Your task is to map the provided list of raw tissue descriptions to the most appropriate Standardized Label from the specific list below.

    <standard_labels>
    {ontology}
    </standard_labels>

    <rules>
    1. Output MUST be a valid JSON dictionary.
    2. Keys must be the exact strings from the Input List.
    3. Values must be one of the Standard Labels.
    4. "Aerial parts" or "ground parts" should generally map to "shoot" (unless it is clearly a rosette, then "rosette").
    5. "Cotyledons" should map to "leaf" (or "seedling" if context implies whole young plant).
    6. If the term implies the whole organism (e.g., "whole organism", "entire plant"), map to "whole_plant" (or "seedling" if age < 14 days).
    7. If a term is ambiguous or fits none of the specific categories, use "unknown".
    8. Do not include markdown formatting (like ```json). Just the raw JSON string.
    </rules>

    <condensed_context>
    {context}
    </condensed_context>
    <input_list>
    {terms}
    </input_list>
    
    JSON Output:
    """)

    # 5. Build Chain
    chain = prompt_template | llm | StrOutputParser()

    try:
        # 6. Execute
        ontology_str = ", ".join([f'"{x}"' for x in valid_tissues])
        
        # Use the existing cleaner to reduce token usage
        context_metadata = clean_metadata_for_context(context)
        
        response_str = chain.invoke({
            "ontology": ontology_str,
            "context": json.dumps(context_metadata),
            "terms": json.dumps(terms)
        })

        # 7. Clean and Parse Response
        cleaned_response = response_str.replace("```json", "").replace("```", "").strip()
        result_dict = json.loads(cleaned_response)
        
        return result_dict

    except json.JSONDecodeError:
        print(f"Error: LLM returned invalid JSON for tissue batch.")
        return {} 
    except Exception as e:
        print(f"Error in get_batch_tissues: {e}")
        return {}