
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# import umap
import os
import math
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.impute import KNNImputer
import matplotlib.cm as cm
import seaborn as sns
import polars as pl
import sys
import logging
import json
module_dir = './'
sys.path.append(module_dir)
from src.constants import *

def apply_KNN_impute(df:pd.DataFrame,n_neighbors: int):
    imputer = KNNImputer(n_neighbors=n_neighbors)

    # Fit and transform the dataset
    df_imputed = pd.DataFrame(imputer.fit_transform(df), columns=df.columns, index=df.index)
    return df_imputed

def setup_directories():
    """Create output directories if they don't exist."""
    logging.info("Setting up output directories...")
    os.makedirs(GEO_DOWNLOAD_DIR, exist_ok=True)
    os.makedirs(METADATA_OUTPUT_DIR, exist_ok=True)
    os.makedirs(PROCESSED_DATA_FOLDER, exist_ok=True)

def create_probe_to_gene_map(gpl):
    """
    Creates a mapping from probe IDs to gene symbols by searching for common
    gene symbol column names in the GPL table.
    """
    possible_gene_columns = ['ORF', 'Gene Symbol', 'GENE_SYMBOL', 'ILMN_Gene', 'Symbol']
    for col_name in possible_gene_columns:
        if col_name in gpl.table.columns:
            logging.info(f"Found gene symbol column: '{col_name}'")
            mapping_df = gpl.table[['ID', col_name]].dropna()
            return dict(zip(mapping_df['ID'], mapping_df[col_name]))
    logging.warning(f"Could not find a recognized gene symbol column in GPL {gpl.name}.")
    return None

def process_metadata(geo_accession, gse, gsm, save_path = METADATA_OUTPUT_DIR):
    """Filters and saves metadata for a given sample to the ./metadata folder."""
    excluded_keys = [
        'contact', 'date', '_id', 'taxid', 'data_', 'status', 'type', 
        'contributor', 'relation', 'geo_accession', 'row_count', 'organism', 
        'label', 'supplementary_file'
    ]
    filtered_study_meta = {k: v for k, v in gse.metadata.items() if not any(ex in k for ex in excluded_keys)}
    filtered_sample_meta = {k: v for k, v in gsm.metadata.items() if not any(ex in k for ex in excluded_keys)}
    final_metadata = {
        'study_id': geo_accession,
        'sample_id': gsm.name,
        'platform': gse.metadata['platform_id'][0],
        'study_metadata': filtered_study_meta,
        'sample_metadata': filtered_sample_meta
    }
    output_path = os.path.join(save_path, f"{geo_accession}_{gsm.name}.json")
    try:
        with open(output_path, "w") as fp:
            json.dump(final_metadata, fp, indent=4)
    except Exception as e:
        logging.error(f"Failed to save metadata for {gsm.name}: {e}")

def process_sample_data(geo_accession, gsm, probe_to_gene_map):
    """
    Processes expression data for a sample and returns it as a DataFrame
    with a single, uniquely named column.
    """
    if gsm.table.empty:
        logging.warning(f"Sample {gsm.name} from study {geo_accession} has no data table. Skipping.")
        return None

    df = gsm.table.copy()
    df['gene_symbol'] = df['ID_REF'].map(probe_to_gene_map)
    df.dropna(subset=['gene_symbol'], inplace=True)
    df['gene_symbol'] = df['gene_symbol'].map(lambda x: x.split("_")[0])
    if df.empty:
        logging.warning(f"No genes could be mapped for sample {gsm.name}. Trying a diferent mapping.")
        #TODO: a lot of the entries here already have proper mappings, so just let them through
        df = gsm.table.copy()
        df['gene_symbol'] = df['ID_REF'].map(lambda x: x.split("_")[0])
        # del df['ID_REF']
        # return None

    # Filter for Arabidopsis thaliana genes (AtXgXXXXX format)
    at_gene_regex = r'^At[1-5MC]g\d{5}$'
    df = df[df['gene_symbol'].str.match(at_gene_regex, case=False)]
    if df.empty:
        logging.error(f"No Arabidopsis thaliana genes found in {gsm.name}. Skipping.")
        return None

    df.set_index('gene_symbol', inplace=True)
    df = df[['VALUE']]
    
    # Handle duplicate genes by averaging their expression
    df = df.groupby(df.index).mean()
    
    # Create a unique name for the sample's data column
    unique_sample_id = f"{geo_accession}_{gsm.name}"
    df = df.rename(columns={'VALUE': unique_sample_id})
    
    return df

def handle_duplicates(df_index,sample_df):

    if sample_df.index.duplicated().any():
        # print('Duplicates in df:', in_df[in_df.index.duplicated(keep=False)])

        # Create a dictionary to keep track of the count of duplicates
        duplicate_count = sample_df.index.value_counts().to_dict()

        # Group by the index column and calculate the mean of the values
        sample_df = sample_df.groupby('gene_symbol')[sample_df.columns[0]].mean().reset_index() # assuming there is only one value
    return sample_df


def get_geo_list(path:str):
    read =  pd.read_csv(path)
    read = read.loc[read['depository_source'] == 'GEO']
    read = read.loc[read['species'] == 'Arabidopsis thaliana']
    return list(read['depository_accession'])
def mapping(x):
    if type(x) is str:
        return x.upper()
    else:
        return x
    
def predicate(gene:str, chromosome:str)-> bool:
    return str('AT'+chromosome+'G') in gene

def get_first_indexs(df_index,chromo:list[str]):
    array = []
    for i in chromo:
        gene:str = next(filter(lambda x : predicate(x,str(i)), df_index))
        array.append(df_index.get_loc(gene))
    return array


def plot_sim_matrix(matrix: np.array, indices: list = None, chromosomes: list = None, 
                   name: str = '', save_loc: str = '', title: str = ''):
    """
    Plot similarity matrix and save to specified location, creating directories if needed.
    
    Args:
        matrix: Input data matrix
        indices: List of indices to split the matrix
        chromosomes: List of chromosome names for labeling
        name: Additional name identifier for output file
        save_loc: Base directory to save outputs
        title: Plot title
    """
    # Determine folder structure
    folder = 'Genes/'
    if indices is None:
        indices = [0]
        folder = 'Samples/'

    if chromosomes is None:
        chromosomes = ['']

    # Create directories if they don't exist
    output_dir = os.path.join(save_loc, 'sim_matrix', folder)
    os.makedirs(output_dir, exist_ok=True)

    for i, _ in enumerate(indices):
        # print('Plotting sim matrix', i)
        min_idx = indices[i]
        try:
            max_idx = indices[i+1]
        except IndexError:
            max_idx = len(matrix)
        
        # print('Computing similarity')
        # Compute pairwise cosine similarity
        similarity_matrix = cosine_similarity(matrix[min_idx:max_idx])
        
        # print('Creating plot')
        plt.imshow(similarity_matrix, cmap='hot', interpolation='nearest')
        plt.colorbar()
        plt.title(title)
        
        # Construct output path
        output_path = os.path.join(output_dir, f'sim_{chromosomes[i]}_matrix{name}.svg')
        plt.savefig(output_path)
        plt.close()
        # print(f'Finished plot saved to {output_path}')

    plt.close()
    # print('Done with all similarity plots')



def box_plot(df: pd.DataFrame, cols_per_plot: int, out_path: str):
    """
    Generates and saves boxplots from a DataFrame, with visual separators between studies.
    
    Args:
        df (pd.DataFrame): The input data.
        cols_per_plot (int): The number of columns (samples) to include in each plot.
        out_path (str): The directory to save the output plots.
    """
    num_cols = len(df.columns)
    num_plots = math.ceil(num_cols / cols_per_plot)
    
    # Ensure the output directory exists
    os.makedirs(out_path, exist_ok=True)

    for plot_num in range(num_plots):
        start_idx = plot_num * cols_per_plot
        end_idx = min((plot_num + 1) * cols_per_plot, num_cols)
        
        current_cols = df.iloc[:, start_idx:end_idx]
        
        # Dynamically adjust figure width based on the number of columns
        # This prevents the plots from looking too compressed.
        fig_width = max(20, len(current_cols.columns) * 0.5) 
        plt.figure(figsize=(fig_width, 10))
        
        # Set a fixed y-axis limit for consistent comparison across plots
        plt.ylim(-18, 18)
        
        # Create the boxplot
        plt.boxplot(current_cols, labels=current_cols.columns)
        
        # --- NEW: Add vertical lines between studies ---
        # Extract study IDs from column names (e.g., 'GSE12345' from 'GSE12345_sample1')
        study_ids = [name.split('_')[0] for name in current_cols.columns]
        
        # Iterate through the columns to find where the study ID changes
        for i in range(1, len(study_ids)):
            if study_ids[i] != study_ids[i-1]:
                # Add a vertical line. Positions are 1-based, so the line goes
                # at i + 0.5 to be between box i and box i+1.
                plt.axvline(x=i + 0.5, color='black', linestyle='--', linewidth=1)
        # --- END NEW ---

        # Rotate x-axis labels for better readability
        plt.xticks(rotation=45, ha='right')
        
        # Add title and adjust layout to prevent labels from being cut off
        plt.title(f'Boxplot Group {plot_num + 1} (Columns {start_idx + 1}-{end_idx})')
        plt.tight_layout()
        
        # Save the figure and close it to free up memory
        plt.savefig(os.path.join(out_path, f'boxplot_group_{plot_num + 1}.png'))
        plt.close()

def find_and_plot_missing_genes(present_genes,out_opath, chr):
    """
    Identifies and plots the location of missing AGI gene locus identifiers 
    (assumed to follow the ATxGyzzzz pattern with increments of 10) 
    within each chromosome.

    Args:
        present_genes (list): A list of existing AGI locus identifiers (strings).
    """

    # --- 1. SIMULATION and DATA CLEANING ---
    # Extract chromosome and numeric location from present genes
    parsed_present = []
    for gene_id in present_genes:
        try:
            # Assuming AGI format: AT[1-5]G[00000-99990]
            # Chromosome is the 3rd character (index 2)
            chromosome = gene_id[2]
            # Numeric ID is the last 5 digits (index -5 onwards)
            numeric_id = int(gene_id[-5:])
            parsed_present.append((chromosome, numeric_id, gene_id))
        except (ValueError, IndexError):
            # Skip any IDs that don't match the expected structure
            continue

    if not parsed_present:
        print("Error: No valid AGI gene IDs found in the input list (e.g., AT1G01010).")
        return

    # Find the maximum numeric ID seen for each chromosome to set the simulation range
    max_ids = {}
    for chr, num_id, _ in parsed_present:
        max_ids[chr] = max(max_ids.get(chr, 0), num_id)

    # Convert the present genes back into a set for fast lookup
    present_gene_set = set(g for _, _, g in parsed_present)
    
    # --- 2. IDENTIFY MISSING GENES ---
    missing_genes = []
    
    # Iterate through all chromosomes (1 to 5) and the simulated range
    for chr_num in range(1, 6):
        chr_label = str(chr_num)
        max_id = max_ids.get(chr_label, 0)
        
        # Simulate all expected gene IDs for this chromosome
        # AGI locus identifiers end in zero, e.g. 10010, 10020...
        # We start at the smallest possible ID (01010) and go up to the max seen
        # plus some buffer to ensure we capture the full range.
        
        # Max expected ID is 99990, but we use max_id as a practical limit
        # plus a 100-gene buffer (1000 numeric steps) to cover any potential high-end gaps.
        sim_max = min(99990, max_id + 1000) 

        # Iterate through expected numeric IDs (increments of 10)
        for num_id in range(1010, sim_max + 10, 10): 
            # Format the 5-digit number, e.g., 1010 -> '01010'
            formatted_num = f'{num_id:05d}'
            # Construct the expected AGI ID
            expected_id = f'AT{chr_label}G{formatted_num}'
            
            if expected_id not in present_gene_set:
                # Store the missing gene and its numeric location/chromosome
                missing_genes.append({
                    'chromosome': chr_label,
                    'locus_id': expected_id,
                    'location_index': num_id
                })

    if not missing_genes:
        print("All expected genes within the range of your input list appear to be present.")
        return

    # Convert missing genes to a DataFrame for easier plotting
    missing_df = pd.DataFrame(missing_genes)
    
    # --- 3. PLOTTING ---
    
    print(f"Found {len(missing_genes)} missing gene IDs across chromosomes 1-5.")
    
    # Set up the plot (5 subplots for 5 chromosomes)
    fig, axes = plt.subplots(
        nrows=5, 
        ncols=1, 
        figsize=(12, 10), 
        sharex=True,
        sharey=True
    )
    plt.suptitle('🔍 Missing Arabidopsis thaliana Genes by Chromosome', fontsize=16, y=1.02)
    
    # AGI genes are numbered top (north) to bottom (south)
    plt.xlabel('Gene Location Index (Southward)', fontsize=14) 
    
    # Create the visualization of a karyotype or a chromosome map.
    # 

    for i, chr_num in enumerate(range(1, 6)):
        chr_label = str(chr_num)
        ax = axes[i]
        
        # Filter data for the current chromosome
        chr_data = missing_df[missing_df['chromosome'] == chr_label]
        
        # Plot each missing gene as a point on the chromosome axis
        ax.scatter(
            chr_data['location_index'], 
            # Use a constant y-value to represent the chromosome line
            [1] * len(chr_data), 
            s=5, # size of the marker
            color='red', 
            alpha=0.6,
            label=f'Missing Genes ({len(chr_data)})'
        )
        
        # Set y-axis properties to make it look like a single line
        ax.set_yticks([]) # Remove y-axis ticks
        ax.set_ylim(0.5, 1.5) # A small range to center the line
        
        # Add a horizontal line to represent the chromosome
        ax.axhline(1, color='lightgray', linestyle='-', linewidth=2)
        
        # Add the chromosome label on the left
        ax.set_ylabel(f'Chr {chr_label}', rotation=0, labelpad=30, fontsize=12, weight='bold')
        
        # Set x-axis limit based on the maximum ID for clarity
        max_idx = max_ids.get(chr_label, 1000)
        ax.set_xlim(1000, max_idx + 100)
        
        # Add a title or legend
        ax.legend(loc='upper right', frameon=False)


    plt.tight_layout()
    plt.savefig(f'{out_opath}missingGenesChr{chr}.svg')
    # plt.show()

def plot_platform_usage(df,save_dir):
    """
    Function 1: Plots a histogram (countplot) showing how many samples 
    use a given platform (by id).
    
    Parameters:
    - df: A Polars DataFrame.
    """
    os.makedirs(save_dir,exist_ok=True)
    plt.figure(figsize=(10, 6))
    
    # Extract IDs as a list or numpy array for Seaborn
    # Seaborn accepts sequences (lists/arrays) directly
    platform_ids = df['id'].to_list()
    
    ax = sns.countplot(x=platform_ids, palette='viridis', hue=platform_ids, legend=False)
    
    plt.title('Number of Samples per Platform ID')
    plt.xlabel('Platform ID')
    plt.ylabel('Count of Samples')
    plt.xticks(rotation=45)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(save_dir+'platform_use.svg')
    plt.close()

def plot_platform_intensities(df:pl.DataFrame, save_dir, key,lengend =False):
    """
    Function 2: For each platform (by id), plots an intensity plot of the data.
    Uses a different color for each sample with transparency.
    
    Parameters:
    - df: A Polars DataFrame containing the data.
    - bucket_size: The bin width for the histogram (intensity bucket).
    """
    os.makedirs(save_dir,exist_ok=True)
    # Get unique platforms
    platforms = df['id'].unique().to_list()

    for platform_id in platforms:
        # platform_id ='GPL198'
        # 1. Filter data for this specific platform using Polars syntax
        platform_subset = df.filter(pl.col('id') == platform_id)
        
        
        # 2. Explode the 'data' column using Polars syntax.
        # This expands the list[float] column into individual rows.
        exploded_data = platform_subset.explode('data')
        
        # 3. Plotting
        # Create a dictionary of numpy arrays for Seaborn
        # This avoids converting the entire object to a Pandas DataFrame
        data_array = np.nan_to_num(exploded_data['data'].to_numpy())
        plot_data = {
            'data': data_array,
            'sample_id': exploded_data['sample_id'].to_numpy(),
            'study_id': exploded_data['study_id'].to_numpy()
        }
        bucket_size = data_array.max()/300

        

        plt.figure(figsize=(12, 7))
        
        ax = sns.histplot(
            data=plot_data,
            x='data',
            hue=key,       # Different color for each sample
            binwidth=bucket_size,  # The requested bucket size
            element="step",        # 'step' style is cleaner for overlapping distributions
            fill=True,             # Fill the area under the curve/step
            alpha=0.1,             # Transparency to see overlaps
            # kde=False              # Disable KDE
            legend=lengend
        )
        df.remove(pl.col('id')==platform_id)
        plt.title(f'Intensity Distribution: {platform_id}')
        plt.xlabel('Intensity Value')
        plt.ylabel('Frequency')
        plt.grid(True, linestyle=':', alpha=0.6)
        plt.tight_layout()
        plt.savefig(save_dir+f'platform_intensity_{platform_id}.svg')
        plt.close()
        # return



class tracker():
    def __init__(self,storage) -> None:
        self.schema_dic:dict = {'id':str,'sample_id':str,'study_id':str,'data':list[float]}
        self.platform_counter:pl.DataFrame = pl.DataFrame(schema=self.schema_dic)#.set_index('id')
        os.makedirs(storage,exist_ok=True)
        self.storage_loc:str = storage+'tracker.json'
    
    def update_counter(self,sample_id,sample_info,study_id,data) -> None:
        # new_data = np.array(data)
        new_record:pl.DataFrame = pl.DataFrame([{'id':sample_info.name,'sample_id':sample_id,'study_id':study_id,'data':data}],schema=self.schema_dic)#.set_index('id')
        new_df= pl.concat([self.platform_counter,new_record])
        self.platform_counter = new_df
        pass

    def store(self):
        self.platform_counter.write_json(self.storage_loc)

    def load(self):
        self.platform_counter = pl.read_json(self.storage_loc).match_to_schema(self.schema_dic)

    def log_transform_high_intensity_studies(self):
        df = self.platform_counter
        max = df.explode('data')['data'].to_numpy().max()
        self.platform_counter = df.with_columns(
            pl.when(
                pl.col("data").list.max().max().over("study_id") > 300
            ).then(
                # Apply log2 to every element in the list using list.eval
                pl.col("data").list.eval((pl.element()+1).log(2))
            ).otherwise(
                pl.col("data")
            )
        )
    def plot(self):
        df = self.platform_counter.clone()

        plot_platform_usage(df,'test_plots/')
        # plot_platform_intensities(df,'test_plots/intensity_plots_by_sample/','sample_id')
        plot_platform_intensities(df,'test_plots/intensity_plots_by_study/','study_id',False)



##### FOR THE new_download.py file
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def plot_tracker_results(json_file="tracker_results.json", output_dir="."):
    """
    Reads tracker JSON and saves SVG plots:
    1. Pie Chart: Total Studies (Used vs Skipped)
    2. Pie Chart: Total Samples (Used vs Skipped)
    3. Bar Charts: Platforms (Filtered by >0 samples used)
    """
    # 1. Load Data
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f"Error: {json_file} not found.")
        return

    totals = data['totals']
    platform_data = data['platform_counts']

    # --- PLOT 1: PIE CHART (STUDIES) ---
    studies_used = totals['total_studies_used']
    studies_seen = totals['total_studies_seen']
    studies_skipped = studies_seen - studies_used

    plt.figure(figsize=(6, 6))
    plt.pie(
        [studies_used, studies_skipped], 
        labels=[f'Used ({studies_used})', f'Skipped ({studies_skipped})'],
        autopct='%1.1f%%',
        colors=sns.color_palette('pastel')[0:2],
        startangle=140
    )
    plt.title(f"Total Studies Seen: {studies_seen}")
    plt.savefig(f"{output_dir}/tracker_pie_studies.svg", format='svg')
    plt.close()
    print(f"Saved {output_dir}/tracker_pie_studies.svg")

    # --- PLOT 2: PIE CHART (SAMPLES) ---
    samples_used = totals['total_samples_used']
    samples_seen = totals['total_sample_seen']
    samples_skipped = samples_seen - samples_used

    plt.figure(figsize=(6, 6))
    plt.pie(
        [samples_used, samples_skipped], 
        labels=[f'Used ({samples_used})', f'Skipped ({samples_skipped})'],
        autopct='%1.1f%%',
        colors=sns.color_palette('pastel')[2:4],
        startangle=140
    )
    plt.title(f"Total Samples Seen: {samples_seen}")
    plt.savefig(f"{output_dir}/tracker_pie_samples.svg", format='svg')
    plt.close()
    print(f"Saved {output_dir}/tracker_pie_samples.svg")

    # --- PLOT 3: BAR CHARTS (PLATFORMS) ---
    # Filter and Reshape Data
    records = []
    
    for platform, stats in platform_data.items():
        # CRITICAL FILTER: Only include platforms where we actually used samples
        if stats['samples_used'] > 0:
            
            # Add Seen Data
            records.append({
                'Platform': platform,
                'Count': stats['studies_seen'],
                'Metric': 'Studies',
                'Status': 'Seen'
            })
            records.append({
                'Platform': platform,
                'Count': stats['samples_seen'],
                'Metric': 'Samples',
                'Status': 'Seen'
            })
            
            # Add Used Data
            records.append({
                'Platform': platform,
                'Count': stats['studies_used'],
                'Metric': 'Studies',
                'Status': 'Used'
            })
            records.append({
                'Platform': platform,
                'Count': stats['samples_used'],
                'Metric': 'Samples',
                'Status': 'Used'
            })

    if not records:
        print("No platforms found with used samples. Skipping bar plots.")
        return

    df = pd.DataFrame(records)

    # Sort by Used Samples to make the plot readable
    # We find the order of platforms based on 'Samples' + 'Used' count
    sorter = df[(df['Metric'] == 'Samples') & (df['Status'] == 'Used')].sort_values('Count', ascending=False)
    order = sorter['Platform'].tolist()

    # Create Subplots (1 Row, 2 Columns)
    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(16, 8), sharey=True)

    # Subplot A: Studies
    studies_df = df[df['Metric'] == 'Studies']
    sns.barplot(
        data=studies_df, 
        x='Count', y='Platform', hue='Status', 
        order=order, ax=axes[0], palette="muted"
    )
    axes[0].set_title("Studies per Platform (Filtered)")
    axes[0].set_xlabel("Number of Studies")

    # Subplot B: Samples
    samples_df = df[df['Metric'] == 'Samples']
    sns.barplot(
        data=samples_df, 
        x='Count', y='Platform', hue='Status', 
        order=order, ax=axes[1], palette="muted"
    )
    axes[1].set_title("Samples per Platform (Filtered)")
    axes[1].set_xlabel("Number of Samples")

    plt.tight_layout()
    plt.savefig(f"{output_dir}/tracker_platforms.svg", format='svg')
    plt.close()
    print(f"Saved {output_dir}/tracker_platforms.svg")

def plot_tracker_results_RNA(json_file="tracker_results.json", output_dir="."):
    """
    Reads tracker JSON and saves SVG plots:
    1. Pie Chart: Total Studies (Used vs Skipped)
    2. Pie Chart: Total Samples (Used vs Skipped)
    3. Bar Charts: Platforms (Filtered by >0 samples used)
    """
    # 1. Load Data
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f"Error: {json_file} not found.")
        return

    totals = data['totals']
    platform_data = data['platform_counts']

    # --- PLOT 1: PIE CHART (STUDIES) ---
    studies_used = totals['total_studies_used']
    studies_seen = totals['total_studies_seen']
    studies_skipped = studies_seen - studies_used

    plt.figure(figsize=(6, 6))
    plt.pie(
        [studies_used, studies_skipped], 
        labels=[f'Used ({studies_used})', f'Skipped ({studies_skipped})'],
        autopct='%1.1f%%',
        colors=sns.color_palette('pastel')[0:2],
        startangle=140
    )
    plt.title(f"Total Studies Seen: {studies_seen}")
    plt.savefig(f"{output_dir}/tracker_pie_studies.svg", format='svg')
    plt.close()
    print(f"Saved {output_dir}/tracker_pie_studies.svg")

    # --- PLOT 2: PIE CHART (SAMPLES) ---
    samples_used = totals['total_samples_used']
    samples_seen = totals['total_sample_seen']
    samples_skipped = samples_seen - samples_used

    plt.figure(figsize=(6, 6))
    plt.pie(
        [samples_used, samples_skipped], 
        labels=[f'Used ({samples_used})', f'Skipped ({samples_skipped})'],
        autopct='%1.1f%%',
        colors=sns.color_palette('pastel')[2:4],
        startangle=140
    )
    plt.title(f"Total Samples Seen: {samples_seen}")
    plt.savefig(f"{output_dir}/tracker_pie_samples.svg", format='svg')
    plt.close()
    print(f"Saved {output_dir}/tracker_pie_samples.svg")

    # --- PLOT 3: BAR CHARTS (PLATFORMS) ---
    # Filter and Reshape Data
    records = []
    
    for platform, stats in platform_data.items():
        # CRITICAL FILTER: Only include platforms where we actually used samples
        if stats['samples_with_raw'] > 0:
            
            # Add Seen Data
            records.append({
                'Platform': platform,
                'Count': stats['studies_seen'],
                'Metric': 'Studies',
                'Status': 'Seen'
            })
            records.append({
                'Platform': platform,
                'Count': stats['samples_seen'],
                'Metric': 'Samples',
                'Status': 'Seen'
            })
            
            # Add Used Data
            records.append({
                'Platform': platform,
                'Count': stats['studies_with_raw'],
                'Metric': 'Studies',
                'Status': 'Used'
            })
            records.append({
                'Platform': platform,
                'Count': stats['samples_with_raw'],
                'Metric': 'Samples',
                'Status': 'Used'
            })

    if not records:
        print("No platforms found with used samples. Skipping bar plots.")
        return

    df = pd.DataFrame(records)

    # Sort by Used Samples to make the plot readable
    # We find the order of platforms based on 'Samples' + 'Used' count
    sorter = df[(df['Metric'] == 'Samples') & (df['Status'] == 'Used')].sort_values('Count', ascending=False)
    order = sorter['Platform'].tolist()

    # Create Subplots (1 Row, 2 Columns)
    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(16, 8), sharey=True)

    # Subplot A: Studies
    studies_df = df[df['Metric'] == 'Studies']
    sns.barplot(
        data=studies_df, 
        x='Count', y='Platform', hue='Status', 
        order=order, ax=axes[0], palette="muted"
    )
    axes[0].set_title("Studies per Platform (Filtered)")
    axes[0].set_xlabel("Number of Studies")

    # Subplot B: Samples
    samples_df = df[df['Metric'] == 'Samples']
    sns.barplot(
        data=samples_df, 
        x='Count', y='Platform', hue='Status', 
        order=order, ax=axes[1], palette="muted"
    )
    axes[1].set_title("Samples per Platform (Filtered)")
    axes[1].set_xlabel("Number of Samples")

    plt.tight_layout()
    plt.savefig(f"{output_dir}/tracker_platforms.svg", format='svg')
    plt.close()
    print(f"Saved {output_dir}/tracker_platforms.svg")


import os
import pandas as pd
import sys
from typing import Union, Callable, Tuple, Dict

def combine_files_microarray(
    folder: str, 
    new_file_name: str, 
    new_file_location: str, 
    combine_genes: bool = False, 
    combination_method: Union[str, Callable] = 'mean'
) -> Tuple[pd.DataFrame, Dict[str, str]]:
    """
    Combines individual GSE RMA CSV files into a single large expression matrix
    and generates a mapping of Samples to Study IDs.
    
    Args:
        folder (str): Root directory containing GSE subfolders.
        new_file_name (str): Output filename (e.g. 'combined_data.csv').
        new_file_location (str): Output directory.
        combine_genes (bool): If True, merges 'Gene.1', 'Gene.2' into 'Gene'.
        combination_method (str or callable): Method to merge overlapping values.

    Returns:
        Tuple[pd.DataFrame, Dict]: (The combined expression dataframe, Dictionary {SampleID: StudyID})
    """
    
    # 1. Setup Output Directory
    if not os.path.exists(new_file_location):
        try:
            os.makedirs(new_file_location)
        except OSError as e:
            print(f"Error creating directory {new_file_location}: {e}")
            raise e

    dataframes = []
    sample_to_study_map = {} # Dictionary to store Sample -> Study relationship
    
    print(f"Scanning '{folder}' for processed files...")

    # 2. Iterate over the folder structure
    for root, dirs, files in os.walk(folder):
        for d in dirs:
            # The folder name 'd' is assumed to be the Study ID (e.g., GSE12345)
            study_id = d 
            
            expected_csv_name = f"{d}_RMA_Genes.csv"
            file_path = os.path.join(root, d, expected_csv_name)
            
            if os.path.exists(file_path):
                print(f"  - Found: {expected_csv_name} (Study: {study_id})")
                try:
                    df = pd.read_csv(file_path, index_col=0)
                    
                    # Requirement: Set all genes to full capitalization
                    df.index = df.index.str.upper()
                    
                    # Local deduplication
                    if df.index.duplicated().any():
                        df = df.groupby(df.index).mean()
                    
                    # --- NEW LOGIC: Map Samples to Study ---
                    # df.columns contains the GSM (Sample) IDs
                    for sample_id in df.columns:
                        sample_to_study_map[sample_id] = study_id
                    # ---------------------------------------

                    dataframes.append(df)
                    
                except Exception as e:
                    print(f"    ! Error reading {expected_csv_name}: {e}")

    if not dataframes:
        print("No valid CSV files found to combine.")
        raise ValueError("No valid CSV files found to combine.")

    # 3. Combine Dataframes
    print(f"\nMerging {len(dataframes)} datasets... (This may take a moment)")
    
    combined_df = pd.concat(dataframes, axis=1, join='outer', sort=True)
    
    # 4. Optional: Combine Gene Variants (Gene.1 -> Gene)
    if combine_genes:
        print(f"  - Consolidating gene variants (e.g., 'Gene.1' -> 'Gene') using method: {combination_method}...")
        cleaned_index = combined_df.index.str.split('.').str[0]
        combined_df.index = cleaned_index
        combined_df = combined_df.groupby(combined_df.index).agg(combination_method)

    elif combined_df.index.duplicated().any():
        print("  - resolving duplicate indices in final merge (using mean)...")
        combined_df = combined_df.groupby(combined_df.index).mean()

    # 5. Save Data Matrix
    output_path = os.path.join(new_file_location, new_file_name)
    print(f"Saving combined matrix to: {output_path}")
    
    try:
        combined_df.to_csv(output_path)
        print("SUCCESS: Data matrix saved.")
        print(f"Final Dimensions: {combined_df.shape[0]} Genes x {combined_df.shape[1]} Samples")
    except Exception as e:
        print(f"Error saving data file: {e}")

    # 6. Save Sample Map
    # We generate a filename based on the input name: "data.csv" -> "data_sample_map.csv"
    base_name = os.path.splitext(new_file_name)[0]
    map_filename = f"{base_name}_sample_map.csv"
    map_output_path = os.path.join(new_file_location, map_filename)
    
    print(f"Saving sample map to: {map_output_path}")
    try:
        # Convert dict to DataFrame for easy saving
        map_df = pd.DataFrame.from_dict(sample_to_study_map, orient='index', columns=['StudyID'])
        map_df.index.name = 'SampleID'
        map_df.to_csv(map_output_path)
        print("SUCCESS: Sample map saved.")
    except Exception as e:
        print(f"Error saving map file: {e}")
        
    return combined_df, sample_to_study_map
import matplotlib
# Force matplotlib to not use any Xwindow/GUI backend.
# This must be called BEFORE importing pyplot.
matplotlib.use('Agg') 

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os

import datashader as ds
import datashader.transfer_functions as tf
from datashader.utils import export_image
import pandas as pd
import numpy as np
import os
import datashader as ds
import datashader.transfer_functions as tf
from datashader.utils import export_image
import pandas as pd
import numpy as np
import os

def plot_study_distributions_datashader(df, sample_map, save_file_path, bins=200, width=1200, height=800):
    """
    Uses Datashader to plot density curves aggregated by STUDY.
    Instead of 1 line per sample, it plots 1 line per Study.
    
    Args:
        df (pd.DataFrame): Genes (rows) x Samples (columns).
        sample_map (dict): Dictionary mapping {SampleID: StudyID}.
        save_file_path (str): Output path (without extension).
    """
    # 1. Ensure output directory exists
    directory = os.path.dirname(save_file_path)
    if directory and not os.path.exists(directory):
        os.makedirs(directory)
        
    print("  - [Datashader] Preparing study-level data...")

    # 2. Clean Data
    df_numeric = df.apply(pd.to_numeric, errors='coerce').dropna(axis=1, how='all')
    
    if df_numeric.empty:
        print("Error: No numeric data to plot.")
        return

    # 3. Group samples by Study
    # Invert map for easier access: {StudyID: [Sample1, Sample2, ...]}
    study_groups = {}
    for sample, study in sample_map.items():
        if sample in df_numeric.columns:
            if study not in study_groups:
                study_groups[study] = []
            study_groups[study].append(sample)

    n_studies = len(study_groups)
    print(f"  - [Datashader] Aggregating {df_numeric.shape[1]} samples into {n_studies} study distributions...")

    # 4. Pre-allocate arrays
    # One line per study -> n_studies lines
    total_rows = n_studies * (bins + 1)
    
    xs = np.zeros(total_rows, dtype=np.float32)
    ys = np.zeros(total_rows, dtype=np.float32)
    
    global_min = df_numeric.min().min()
    global_max = df_numeric.max().max()
    bin_edges = np.linspace(global_min, global_max, bins + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # 5. Compute Density per Study
    for i, (study_id, samples) in enumerate(study_groups.items()):
        # Extract all data for this study as a single flat array
        # This treats the study as one giant distribution of values
        study_data = df_numeric[samples].values.flatten()
        
        # Remove NaNs
        study_data = study_data[~np.isnan(study_data)]
        
        if len(study_data) == 0:
            continue
            
        # Compute Histogram
        hist, _ = np.histogram(study_data, bins=bin_edges, density=True)
        
        # Fill arrays
        start_idx = i * (bins + 1)
        end_idx = start_idx + bins
        
        xs[start_idx:end_idx] = bin_centers
        ys[start_idx:end_idx] = hist
        
        # NaN separator
        xs[end_idx] = np.nan
        ys[end_idx] = np.nan

    # 6. Render
    print("  - [Datashader] Rendering...")
    line_df = pd.DataFrame({'x': xs, 'y': ys})
    
    cvs = ds.Canvas(plot_width=width, plot_height=height)
    agg = cvs.line(line_df, 'x', 'y', agg=ds.count())
    img = tf.shade(agg, cmap=['lightblue', 'darkblue', 'red'], how='log')
    
    if save_file_path.endswith('.png'):
        save_file_path = save_file_path[:-4]
        
    export_image(img, save_file_path, background="white", export_path=".")
    print(f"  - [Datashader] Success! Saved to {save_file_path}.png")


import datashader as ds
import datashader.transfer_functions as tf
from datashader.utils import export_image
import pandas as pd
import numpy as np
import os

def plot_study_distributions_incremental(folder_path, save_file_path, bins=200, width=1200, height=800):
    """
    Scans a folder of processed CSVs (e.g., GSE123_RMA_Genes.csv), calculates the 
    density distribution for each study, and incrementally adds it to a Datashader plot.
    
    Args:
        folder_path (str): Directory containing the per-study CSV files.
        save_file_path (str): Output path (without extension).
    """
    
    # 1. Setup Output
    directory = os.path.dirname(save_file_path)
    if directory and not os.path.exists(directory):
        os.makedirs(directory)
        
    print(f"Scanning '{folder_path}' for study files...")
    
    # 2. Initialize Datashader Canvas
    # We need to know the global min/max to set the canvas range correctly.
    # Since we can't load all data to find min/max, we estimate or do a quick pre-scan.
    # For normalized expression (Log2), 0 to 16 is a safe standard range.
    x_range = (0, 16) 
    y_range = (0, 1.0) # Density usually falls 0-1
    
    cvs = ds.Canvas(plot_width=width, plot_height=height, x_range=x_range, y_range=y_range)
    
    # This will hold the aggregated counts (the 'image' data)
    agg = None 

    studies_processed = 0
    
    # 3. Iterate over files
    for root, dirs, files in os.walk(folder_path):
        for d in dirs:
            # Construct path to the expected CSV
            csv_path = os.path.join(root, d, f"{d}_RMA_Genes.csv")
            
            if os.path.exists(csv_path):
                try:
                    # Load ONE study into memory
                    df = pd.read_csv(csv_path, index_col=0)
                    
                    # Ensure numeric
                    df = df.apply(pd.to_numeric, errors='coerce')
                    
                    # --- METHOD: One line per study (Mean of samples) ---
                    # We average all samples to get a "Study Consensus" profile
                    study_profile = df.mean(axis=1).dropna().values
                    
                    if len(study_profile) == 0:
                        continue
                        
                    # Compute Histogram for this study
                    # We map this profile into x,y coordinates for Datashader
                    hist, bin_edges = np.histogram(study_profile, bins=bins, range=x_range, density=True)
                    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
                    
                    # Create a mini dataframe for just this line
                    line_df = pd.DataFrame({'x': bin_centers, 'y': hist})
                    
                    # Add to the global aggregation
                    # If this is the first study, create the agg. Otherwise, add to it.
                    new_agg = cvs.line(line_df, 'x', 'y', agg=ds.count())
                    
                    if agg is None:
                        agg = new_agg
                    else:
                        agg = agg + new_agg
                    
                    studies_processed += 1
                    if studies_processed % 50 == 0:
                        print(f"  - Processed {studies_processed} studies...")
                        
                except Exception as e:
                    print(f"    ! Error reading {d}: {e}")

    if agg is None:
        print("Error: No valid data found to plot.")
        return

    # 4. Render and Save
    print(f"  - Rendering final image for {studies_processed} studies...")
    
    # 'log' mapping works best to see the "main trend" vs "outlier studies"
    img = tf.shade(agg, cmap=['lightblue', 'darkblue', 'red'], how='log')
    
    # if save_file_path.endswith('.png'):
    #     save_file_path = save_file_path[:-4]
        
    export_image(img, save_file_path, background="white", export_path=".")
    print(f"  - Success! Saved to {save_file_path}.png")


import seaborn as sns
import seaborn.objects as so
import matplotlib
# Use Agg backend to prevent display errors/crashes
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

def plot_study_distributions_seaborn(folder_path, save_file_path, bins=200):
    """
    Scans processed CSVs and plots density distributions using the Seaborn Objects API.
    
    Visualization Layers:
    1. Individual Study Lines (Light Blue)
    2. Standard Deviation Band (Gray Area)
    3. Mean Distribution Line (Red Line)
    """
    
    # 1. Setup Output
    directory = os.path.dirname(save_file_path)
    if directory and not os.path.exists(directory):
        os.makedirs(directory)

    print(f"Scanning '{folder_path}' for study files...")
    
    # Standard bounds for Log2 Expression
    x_range = (0, 16) 
    bin_edges = np.linspace(x_range[0], x_range[1], bins + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_width = (x_range[1] - x_range[0]) / bins

    plot_data_list = []
    studies_processed = 0

    # 2. Iterate and Compute Densities
    for root, dirs, files in os.walk(folder_path):
        for d in dirs:
            csv_path = os.path.join(root, d, f"{d}_RMA_Genes.csv")
            
            if os.path.exists(csv_path):
                try:
                    # Load ONE study
                    df = pd.read_csv(csv_path, index_col=0)
                    df = df.apply(pd.to_numeric, errors='coerce')
                    
                    # Calculate Study Consensus (Mean of samples)
                    study_profile = df.mean(axis=1).dropna().values
                    
                    if len(study_profile) == 0:
                        continue
                        
                    # Compute Histogram (Density)
                    hist, _ = np.histogram(study_profile, bins=bin_edges, density=True)
                    
                    # --- ASSERTION CHECK ---
                    integral = np.sum(hist) * bin_width
                    assert np.isclose(integral, 1.0, atol=1e-4), \
                        f"Integration Check Failed for {d}: Area = {integral}"
                    # -----------------------

                    # Store as a lightweight DataFrame
                    temp_df = pd.DataFrame({
                        'Expression': bin_centers,
                        'Density': hist,
                        'Study': d
                    })
                    
                    plot_data_list.append(temp_df)
                    studies_processed += 1
                    
                    if studies_processed % 50 == 0:
                        print(f"  - Processed {studies_processed} studies...")
                        
                except Exception as e:
                    print(f"    ! Error reading {d}: {e}")

    if not plot_data_list:
        print("Error: No valid data found to plot.")
        return

    print(f"  - Aggregating data for plotting...")
    plotting_df = pd.concat(plot_data_list, ignore_index=True)

    # 3. Plotting with Seaborn Objects
    print("  - Generating Seaborn Objects plot...")
    
    # Ensure correct extension
    if not save_file_path.endswith('.svg'):
        save_file_path += '.svg'

    (
        so.Plot(plotting_df, x='Expression', y='Density')
        # Layer 1: Individual Study Lines (Background)
        # We group by 'Study' so each study gets its own thin line
        .add(so.Line(color='cornflowerblue', alpha=0.1, linewidth=0.5), group='Study')
        
        # Layer 2: Standard Deviation (Band)
        # so.Est(errorbar='sd') calculates the standard deviation at each X point
        .add(so.Band(color='red', alpha=0.2), so.Est(errorbar="sd"))
        
        # Layer 3: Mean (Line)
        # so.Agg() calculates the mean at each X point
        .add(so.Line(color='red', linewidth=1.5), so.Agg())
        
        # Layout & Labels
        .label(
            title=f"Gene Expression Density (N={studies_processed})",
            x="Log2 Expression Intensity",
            y="Density"
        )
        .limit(x=x_range, y=(0, 1.2))
        .layout(size=(10, 6))
        
        # Save
        .save(save_file_path, format='svg', dpi=150)
    )
    
    print(f"  - Success! Saved to {save_file_path}")