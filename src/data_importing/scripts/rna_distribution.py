import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

# 1. Configuration
INPUT_FILE = "new_storage/final_data/Salmon_RNAseq_Combined.csv"
OUTPUT_PLOT = "sample_distributions_histogram.png"
MAX_SAMPLES_TO_PLOT = 30  # Change this if you want to plot more/fewer samples

def plot_distributions():
    if not os.path.exists(INPUT_FILE):
        print(f"Error: Could not find {INPUT_FILE}")
        return

    print(f"Loading data from {INPUT_FILE}...")
    # Assuming the first column contains Gene IDs, so we set it as the index
    df = pd.read_csv(INPUT_FILE, index_col=0)
    
    # Drop any metadata columns if they exist (keep only numeric data)
    df = df.select_dtypes(include=[np.number])
    
    print(f"Loaded matrix with {df.shape[0]} genes and {df.shape[1]} samples.")

    # 2. Log Transformation
    # RNA-seq counts are highly skewed. We use log1p (log(x + 1)) for better visualization
    print("Applying log1p transformation...")
    df_log = np.log1p(df)

    # 3. Select Samples to Plot
    # Plotting 100+ samples as individual histograms is impossible to read, 
    # so we select the first N samples (or you can change this to df.columns to do all)
    samples_to_plot = df_log.columns[:MAX_SAMPLES_TO_PLOT]
    
    print(f"Plotting distributions for {len(samples_to_plot)} samples...")

    # 4. Plotting
    plt.figure(figsize=(12, 6))
    
    # Option A: Overlapping Density Plots (Usually best for comparing many RNA-seq samples)
    for sample in samples_to_plot:
        sns.kdeplot(df_log[sample], label=sample, linewidth=1.5, alpha=0.7)
    
    # If you strictly want standard blocky histograms instead, comment out the for-loop above 
    # and uncomment the one below:
    # for sample in samples_to_plot:
    #     sns.histplot(df_log[sample], element="step", fill=False, label=sample, bins=50)

    plt.title("RNA-seq Value Distribution per Sample (Log Transformed)", fontsize=14)
    plt.xlabel("Expression Level (log1p)", fontsize=12)
    plt.ylabel("Density", fontsize=12)
    plt.legend(title="Samples", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    # 5. Save the plot
    plt.savefig(OUTPUT_PLOT, dpi=300)
    print(f"Plot successfully saved to {OUTPUT_PLOT}")

if __name__ == "__main__":
    plot_distributions()