import os
import sys

# Ensure we can import the local modules
module_dir = "./"
sys.path.append(module_dir)

from src.constants import STORAGE_DIR  # noqa: E402
from src.data_importing.RNA_seq_processing_batch import combine_files_rnaseq  # noqa: E402


def main():
    # 1. Define your base paths
    BIG_STORAGE = STORAGE_DIR

    # 2. Set input and output locations
    processed_rnaseq_folder = os.path.join(BIG_STORAGE, "rnaseq_data", "processed_rnaseq")
    final_data_folder = os.path.join(BIG_STORAGE, "final_data", "rnaseq_processed")
    output_filename = "Salmon_RNAseq_Combined_TPM.csv"

    print("Starting RNA-seq merge process...")
    print(f"Input directory: {processed_rnaseq_folder}")
    print(f"Output directory: {final_data_folder}")

    # 3. Call the combining function
    combined_matrix, sample_map = combine_files_rnaseq(
        folder=processed_rnaseq_folder,
        new_file_name=output_filename,
        new_file_location=final_data_folder,
        combination_method="max",  # Matches your microarray logic
    )

    print("\n--- Merge Complete ---")
    print(f"Final Matrix Shape: {combined_matrix.shape}")
    print(f"Total Samples Mapped: {len(sample_map)}")


if __name__ == "__main__":
    main()
