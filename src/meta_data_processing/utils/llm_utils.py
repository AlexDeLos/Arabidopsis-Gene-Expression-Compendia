from langchain_core.prompts import ChatPromptTemplate
from langchain_core.messages import SystemMessage
from langchain_core.output_parsers import PydanticOutputParser, StrOutputParser
from langchain.chat_models import init_chat_model
from typing import List, Literal
from pydantic import ValidationError, BaseModel, Field


class MetaData(BaseModel):
    tissue: Literal["root", "leaf", "flower", "shoot", "rosette", "bud", "whole_plant", "silique", "callus", "seed", "seedling", "unknown"] = Field(..., description="Tissue the samples was extracted from.") # "stem", 
    # tissue: str = Field(..., description="Tissue the samples was extracted from.") # "stem", 

    # Field with description explaining its purpose
    treatment: List[Literal[
    "Drought Stress",
    "Salinity Stress",
    "Heat Stress",
    "Cold Stress",
    "Chemical Stress",
    "Nutrient Deficiency",
    "Pathogen Attack",
    "Low Light Stress",
    "High Light Stress",
    "Red Light Stress",
    "Other Light Stress",
    "Other stress",
    "No stress"]] = Field(..., description="List of treatments and stresses that was applied to the sample, each unique stress or treatment should have one, and only one, entry in this list")
    
    # Field with description explaining its purpose
    medium: str = Field(..., description="Growth medium of the sample.")
    # Field with description explaining its purpose
    # age: float = Field(..., description="The age of the plant since germination in days.") 
    # # Field with description explaining its purpose
    # treatment_time: float = Field(..., description="Time in hours that between the application of the treatment to the sample being harvested.") #! this is not as accurate




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

def get_metadata_script(sample_metadata: dict, study_id: str, model: str = 'gemini-1.5-flash', temp: float = 0):
    """
    Generates a Python extractor function for a specific study based on a sample's metadata.
    """
    
    # Initialize the LLM
    llm = init_chat_model(model=model,
                        model_provider="google_genai",
                        temperature=temp)
    
    # Setup the Parser to get the target Schema instructions
    parser = PydanticOutputParser(pydantic_object=MetaData)
    format_instructions = parser.get_format_instructions()

    # Define the System Prompt
    # We ask for a Python function that outputs the specific JSON schema
    prompt = ChatPromptTemplate.from_messages([
        ("system", "Produce a python code function that can extract the following schema from the data dictionary object provided by the user. \nReturn ONLY the python code.\n{format_instructions}"),
        ("human", "{text}"),
    ]).partial(format_instructions=format_instructions)

    parsing_llm = prompt | llm | StrOutputParser()

    # Construct the Task Prompt
    # We removed 'study_info' as it is not passed in the new signature.
    # We focus purely on the sample_metadata structure.
    parsing_prompt = '''
    <task>
    Provide a python function with the signature:
    def {study_id}_extractor(sample_metadata: dict) -> dict:
    
    This function should extract the relevant biological experimental information (tissue, treatment, medium) from the input dictionary following the schema defined in the system prompt.
    </task>

    <guidance>
    1. The input `sample_metadata` will have the structure shown below.
    2. Note that values in the dictionary are often lists of strings (e.g., "title": ["Name"]). You usually want to access the first element or join them.
    3. Look for keywords in fields like 'characteristics_ch1', 'source_name_ch1', 'title', or 'description'.
    4. If specific information (like medium) is missing, infer "unspecified" or a logical default (e.g. "soil" for whole plants) if strongly implied, otherwise "unspecified".
    5. Return only the valid Python code, no markdown backticks.
    </guidance>

    <metadata_structure>
        The following is an example of the sample_metadata dictionary this function will process:
        <sample_metadata>
            {sample_metadata}
        </sample_metadata>
    </metadata_structure>
    '''

    try:
        # Invoke the chain
        result = parsing_llm.invoke({
            "text": parsing_prompt.format(
                study_id=study_id, 
                sample_metadata=sample_metadata
            )
        })
        
        # Simple cleanup to remove markdown code blocks if the LLM adds them
        result = result.replace("```python", "").replace("```", "").strip()
        
    except Exception as e:
        print(f'Error generating script with {model}: {e}')
        
        # Retry logic with a more powerful model if Flash fails
        if 'flash' in model:
            print('Retrying with a Pro model...')
            # Update to use Pro model for retry (using 1.5-pro as standard)
            return get_metadata_script(sample_metadata, study_id, model='gemini-1.5-pro', temp=temp)
        else:
            raise e

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
            raise e#-> GSM1126262 ->GSM843683
        return get_metadata(study_info,sample_info, model='gemini-2.5-pro', temp=0.1)
    return result