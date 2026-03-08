#!/bin/bash

#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --time=32:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --account=ewi-insy-prb
#SBATCH --job-name=check_fastqs
#SBATCH --output=/home/nfs/alexdelossanto/Dataset_fusion_Microarray/logs_slurm/stdout-%x-%j.txt

module purge
module use /opt/insy/modulefiles
module load miniconda/3.9
eval "$(conda shell.bash hook)"
conda activate download_env

LOG_DIR="/home/nfs/alexdelossanto/Dataset_fusion_Microarray/logs_slurm"
DELETION_LOG="${LOG_DIR}/corrupted_files_$(date +%Y%m%d_%H%M%S).log"

echo "=================================================="
echo " Starting STRICT Python Integrity Scan"
echo " Start time: $(date)"
echo " Deletion log: ${DELETION_LOG}"
echo "=================================================="

python3 - "$DELETION_LOG" << 'EOF'
import os
import sys
import gzip
from glob import glob
from multiprocessing import Pool

DELETION_LOG = sys.argv[1]

def check_gzip_integrity(file_path):
    """
    Reads the entire gzip file through Python's strict gzip parser.
    Returns (file_path, error) if corrupted, None if healthy.
    This catches truncated multi-part streams that `gzip -t` misses.
    """
    try:
        with gzip.open(file_path, 'rb') as f:
            while f.read(1024 * 1024 * 10):  # 10MB chunks — avoids RAM exhaustion
                pass
        return None  # Healthy
    except Exception as e:
        return (file_path, str(e))

if __name__ == '__main__':
    directory = '/tudelft.net/staff-umbrella/GeneExpressionStorage/rnaseq_data/fastq_storage'
    print(f"Locating files in {directory}...")

    files = glob(f"{directory}/**/*.fastq.gz", recursive=True)
    total = len(files)
    print(f"Found {total} files. Scanning with 8 workers...")
    print(f"Progress will update every 100 files.\n")

    corrupted = []
    checked = 0

    # imap_unordered: processes results as they complete rather than waiting for all.
    # More memory efficient and gives real-time progress on large file sets.
    with Pool(8) as pool:
        for result in pool.imap_unordered(check_gzip_integrity, files, chunksize=4):
            checked += 1
            if result is not None:
                file_path, error_msg = result
                corrupted.append((file_path, error_msg))
                print(f"  CORRUPTED [{checked}/{total}]: {file_path}")
                print(f"             Reason: {error_msg}")

            if checked % 100 == 0:
                print(f"  Progress: {checked}/{total} files checked "
                      f"({100*checked//total}%) — {len(corrupted)} corrupted so far.")

    print(f"\nScan complete. {checked} files checked, {len(corrupted)} corrupted.")

    # Write deletion log BEFORE deleting — permanent record even if deletion is interrupted
    with open(DELETION_LOG, 'w') as log:
        import datetime
        log.write(f"Corrupted FASTQ Scan — {datetime.datetime.now()}\n")
        log.write(f"Directory: {directory}\n")
        log.write(f"Total files scanned: {total}\n")
        log.write(f"Corrupted files found: {len(corrupted)}\n\n")
        if corrupted:
            log.write("FILES TO DELETE:\n")
            for file_path, error_msg in corrupted:
                log.write(f"  {file_path}\n")
                log.write(f"    Reason: {error_msg}\n")
        else:
            log.write("No corrupted files found.\n")

    print(f"Deletion log written to: {DELETION_LOG}")

    # Delete with per-file error handling
    if corrupted:
        print("\n==================================================")
        deleted = 0
        failed_to_delete = []
        for file_path, _ in corrupted:
            try:
                os.remove(file_path)
                print(f"  DELETED: {file_path}")
                deleted += 1
            except OSError as e:
                print(f"  FAILED TO DELETE: {file_path} — {e}")
                failed_to_delete.append(file_path)

        print("==================================================")
        print(f"Deleted {deleted}/{len(corrupted)} corrupted files.")

        if failed_to_delete:
            print(f"WARNING: {len(failed_to_delete)} files could not be deleted:")
            for f in failed_to_delete:
                print(f"  {f}")

        # Append deletion outcome to log
        with open(DELETION_LOG, 'a') as log:
            log.write(f"\nDELETION RESULTS:\n")
            log.write(f"  Successfully deleted: {deleted}\n")
            if failed_to_delete:
                log.write(f"  Failed to delete: {len(failed_to_delete)}\n")
                for f in failed_to_delete:
                    log.write(f"    {f}\n")
    else:
        print("\nAll files passed the strict Python integrity check!")
EOF

echo "=================================================="
echo " Scan Complete!"
echo " End time: $(date)"
echo "=================================================="
