import os
import re
import gzip
import shutil
import requests
import pandas as pd
import GEOparse
from typing import Optional

# ---------------------------------------------------------------------------
# Patterns
# ---------------------------------------------------------------------------

# Explicit count keywords — high confidence
_COUNT_KEYWORDS = re.compile(
    r'(count|counts|rawcount|raw_count|gene_count|read_count|'
    r'htseq|featurecount|rsem|salmon|kallisto|star|expression_matrix|'
    r'read_matrix|gene_expression)',
    re.IGNORECASE
)

# Normalised metrics — reject these even if they pass other filters
_NORMALISED = re.compile(
    r'(fpkm|rpkm|tpm|cpm|normalized|normalised|vst|rlog|log2)',
    re.IGNORECASE
)

# File extensions that could be count matrices
_TABULAR_EXT = re.compile(
    r'\.(txt|tsv|csv|tab|xls|xlsx)(\.gz)?$',
    re.IGNORECASE
)

# Hard-skip: archives, known non-count files
_HARD_SKIP = re.compile(
    r'(RAW\.tar|\.tar(\.gz)?$|\.bam$|\.fastq|\.fq|\.bed$|'
    r'README|\.pdf$|\.xml$|\.json$|peak|junction|splice|snp|vcf|'
    r'differential|DE_|DEG_|DESeq_[0-9]|edgeR|fold.?change|padj|pvalue)',
    re.IGNORECASE
)

# ---------------------------------------------------------------------------
# URL helpers
# ---------------------------------------------------------------------------

def _ftp_to_https(url: str) -> str:
    if url.startswith('ftp://ftp.ncbi.nlm.nih.gov/'):
        return url.replace('ftp://ftp.ncbi.nlm.nih.gov/', 'https://ftp.ncbi.nlm.nih.gov/', 1)
    return url


def _list_geo_suppl_dir(gse_id: str) -> list[str]:
    """
    Fetch the GEO series supplementary directory listing via HTTPS and return
    all file URLs found. GEO uses an Apache-style HTML directory index.
    """
    prefix = gse_id[:-3] + 'nnn'
    index_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse_id}/suppl/"
    try:
        r = requests.get(index_url, timeout=30)
        r.raise_for_status()
    except Exception as e:
        print(f"  Could not list FTP directory for {gse_id}: {e}")
        return []
    return [
        index_url + m.group(1)
        for m in re.finditer(r'href="([^"?/][^"]*)"', r.text)
    ]

# ---------------------------------------------------------------------------
# File scoring
# ---------------------------------------------------------------------------

def _score_candidate(filename: str) -> int:
    """
    Score a filename for how likely it is to be a gene-level raw count matrix.
    Higher is better. Returns -1 to disqualify.
    """
    if _HARD_SKIP.search(filename):
        return -1
    if not _TABULAR_EXT.search(filename):
        return -1
    if _NORMALISED.search(filename):
        return -1  # normalised values not usable for DESeq2

    score = 0
    if _COUNT_KEYWORDS.search(filename):
        score += 10
    # Series-wide files (prefixed with GSE ID) are preferred over per-sample
    if re.match(r'GSE\d+', filename, re.IGNORECASE):
        score += 5
    # Compressed files slightly preferred (suggests larger, merged matrices)
    if filename.endswith('.gz'):
        score += 2

    return score


# ---------------------------------------------------------------------------
# Filename extraction from data_processing text
# ---------------------------------------------------------------------------

def _extract_filenames_from_metadata(gse: GEOparse.GEOTypes.GSE) -> list[str]:
    """
    Parse `Sample_data_processing` fields for filenames explicitly mentioned
    by authors (e.g. 'deseq_genes_raw_counts.txt, DESeq_200628_sum_overlap.txt').
    These are sometimes uploaded to FTP without being registered in SOFT metadata.
    """
    filenames = []
    seen = set()
    for gsm in gse.gsms.values():
        for entry in gsm.metadata.get('data_processing', []):
            if 'Supplementary_files_format_and_content' not in entry:
                continue
            # Extract anything that looks like a filename (has an extension)
            for match in re.finditer(r'[\w\-\.]+\.(txt|tsv|csv|tab|gz|xlsx?)', entry, re.IGNORECASE):
                fname = match.group(0)
                if fname not in seen:
                    seen.add(fname)
                    filenames.append(fname)
    return filenames


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def _parse_count_file(path: str) -> Optional[pd.DataFrame]:
    """
    Parse a tabular file into a DataFrame and validate it looks like raw counts.
    Returns None if the file cannot be parsed or fails the count sanity check.
    """
    try:
        compression = 'gzip' if path.endswith('.gz') else None

        # Sniff separator
        opener = gzip.open if compression else open
        with opener(path, 'rt') as f:
            first_line = f.readline()
        sep = '\t' if first_line.count('\t') >= first_line.count(',') else ','

        df = pd.read_csv(
            path, sep=sep, index_col=0,
            compression=compression,
            comment='#',
            low_memory=False
        )

        # Drop HTSeq summary rows (__no_feature etc.)
        df = df[~df.index.astype(str).str.startswith('__')]

        # Keep only numeric columns
        df = df.select_dtypes(include='number')

        if df.empty or df.shape[1] == 0:
            return None

        # Values < 2 max → ratios/fractions, not counts
        if df.values.max() < 2:
            print(f"    Skipping {os.path.basename(path)}: looks normalised (max={df.values.max():.4f})")
            return None

        # Sanity: at least 100 rows (genes) expected
        if df.shape[0] < 100:
            return None

        return df

    except Exception as e:
        print(f"    Could not parse {os.path.basename(path)}: {e}")
        return None


def _download_file(url: str, save_path: str) -> bool:
    """Download a URL to save_path. Returns True on success."""
    if os.path.exists(save_path):
        return True
    try:
        r = requests.get(_ftp_to_https(url), stream=True, timeout=120)
        r.raise_for_status()
        with open(save_path, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
        return True
    except Exception as e:
        print(f"    Download failed for {os.path.basename(save_path)}: {e}")
        if os.path.exists(save_path):
            os.remove(save_path)
        return False


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------

def download_processed_counts(gse_id: str, output_dir: str) -> Optional[pd.DataFrame]:
    """
    Attempt to retrieve a pre-computed gene-level raw count matrix for an
    RNA-seq GEO study. Returns a DataFrame (genes × samples) or None.

    Search order:
      1. Series FTP directory — score all tabular files, try best candidates
         (covers: explicitly named count files AND files with no count keyword
         that happen to be tabular, since many authors don't follow naming
         conventions)
      2. Filenames mentioned in Sample_data_processing metadata text — authors
         often name files there even when supplementary_file_1 = NONE
      3. Sample-level supplementary files via GEOparse — per-sample HTSeq /
         featureCounts files that get merged into one matrix

    Files from studies with no recoverable count matrix are NOT kept on disk.
    """
    study_dir = os.path.join(output_dir, gse_id)
    os.makedirs(study_dir, exist_ok=True)
    downloaded_files: list[str] = []  # track for cleanup on failure

    print(f"\n[{gse_id}] Searching for count matrix...")

    try:
        # ---------------------------------------------------------------
        # Step 1: Score and try all tabular files in series FTP directory
        # ---------------------------------------------------------------
        series_urls = _list_geo_suppl_dir(gse_id)
        print(f"  FTP directory: {len(series_urls)} file(s) found")

        scored = sorted(
            [(url, _score_candidate(url.split('/')[-1])) for url in series_urls],
            key=lambda x: -x[1]
        )
        # Only try files with score >= 0 (not disqualified)
        candidates = [(url, score) for url, score in scored if score >= 0]

        if candidates:
            print(f"  {len(candidates)} tabular candidate(s): "
                  f"{[u.split('/')[-1] for u, _ in candidates]}")

        for url, score in candidates:
            filename = url.split('/')[-1]
            save_path = os.path.join(study_dir, filename)
            print(f"  Trying [{score:+d}] {filename} ...")

            if _download_file(url, save_path):
                downloaded_files.append(save_path)
                df = _parse_count_file(save_path)
                if df is not None:
                    print(f"  ✓ Count matrix found: {df.shape[0]} genes × {df.shape[1]} samples")
                    return df

        # ---------------------------------------------------------------
        # Step 2: Filenames mentioned in data_processing text
        # ---------------------------------------------------------------
        print(f"  Checking data_processing metadata for named files...")
        try:
            gse = GEOparse.get_GEO(geo=gse_id, destdir=study_dir, silent=True)
        except Exception as e:
            print(f"  GEOparse error: {e}")
            gse = None

        if gse is not None:
            mentioned_filenames = _extract_filenames_from_metadata(gse)
            if mentioned_filenames:
                print(f"  Files mentioned in metadata: {mentioned_filenames}")

            # Build candidate URLs: try GSE-prefixed FTP path for each filename
            prefix = gse_id[:-3] + 'nnn'
            base_ftp = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse_id}/suppl/"

            for fname in mentioned_filenames:
                if _score_candidate(fname) < 0:
                    continue
                url = base_ftp + fname
                save_path = os.path.join(study_dir, fname)
                print(f"  Trying metadata-mentioned: {fname} ...")
                if _download_file(url, save_path):
                    downloaded_files.append(save_path)
                    df = _parse_count_file(save_path)
                    if df is not None:
                        print(f"  ✓ Count matrix found via metadata hint: {df.shape[0]} genes × {df.shape[1]} samples")
                        return df

            # ---------------------------------------------------------------
            # Step 3: Sample-level supplementary files — merge per-sample counts
            # ---------------------------------------------------------------
            print(f"  Checking sample-level supplementary files...")
            sample_dfs: list[pd.DataFrame] = []

            for gsm_id, gsm in gse.gsms.items():
                for key, val_list in gsm.metadata.items():
                    if 'supplementary_file' not in key:
                        continue
                    for url in val_list:
                        if not url or url.strip().upper() == 'NONE':
                            continue
                        filename = url.split('/')[-1]
                        score = _score_candidate(filename)
                        if score < 0:
                            continue
                        save_path = os.path.join(study_dir, filename)
                        if _download_file(url, save_path):
                            downloaded_files.append(save_path)
                            df_s = _parse_count_file(save_path)
                            if df_s is not None and df_s.shape[1] == 1:
                                df_s.columns = [gsm_id]
                                sample_dfs.append(df_s)
                        break  # one count file per sample

            if sample_dfs:
                merged = pd.concat(sample_dfs, axis=1)
                print(f"  ✓ Merged {len(sample_dfs)} sample files → {merged.shape[0]} genes × {merged.shape[1]} samples")
                return merged

        # ---------------------------------------------------------------
        # No count matrix found — clean up all downloaded files
        # ---------------------------------------------------------------
        print(f"  ✗ No count matrix found for {gse_id}. Cleaning up.")
        _cleanup(study_dir, downloaded_files)
        return None

    except Exception as e:
        print(f"  Unexpected error for {gse_id}: {e}")
        _cleanup(study_dir, downloaded_files)
        return None


def _cleanup(study_dir: str, files: list[str]) -> None:
    """Remove downloaded files and the study directory if empty."""
    for f in files:
        try:
            if os.path.exists(f):
                os.remove(f)
        except OSError:
            pass
    try:
        # Remove directory only if it's now empty
        if os.path.isdir(study_dir) and not os.listdir(study_dir):
            shutil.rmtree(study_dir)
    except OSError:
        pass