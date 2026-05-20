import json
import os

# Define the path to your json cache file
cache_file_path = "srr_GSM_cache.json"

# 1. Check if the file exists before attempting to open it
if not os.path.exists(cache_file_path):
    print(f"Error: Could not find '{cache_file_path}' in the current directory.")
    exit(1)

# 2. Load the JSON data
print(f"Loading {cache_file_path}...")
with open(cache_file_path, "r") as f:
    try:
        cache_data = json.load(f)
    except json.JSONDecodeError as e:
        print(f"Error: Failed to parse JSON. Please check file formatting. Details: {e}")
        exit(1)

initial_count = len(cache_data)
cleaned_cache = {}

# 3. Process the dictionary
for gsm_key, value_list in cache_data.items():
    # Keep only values that DO NOT start with "GSM"
    filtered_list = [item for item in value_list if not str(item).startswith("GSM")]
    
    # If the array is NOT empty after filtering, retain the key
    if filtered_list:
        cleaned_cache[gsm_key] = filtered_list

final_count = len(cleaned_cache)
removed_keys = initial_count - final_count

# 4. Save the cleaned dictionary back to the JSON file
print("Writing cleaned data back to file...")
with open(cache_file_path, "w") as f:
    json.dump(cleaned_cache, f, indent=4)

print("=" * 50)
print("SUCCESS: Cache cleanup complete!")
print(f"  Initial unique GSM keys: {initial_count}")
print(f"  Final unique GSM keys:   {final_count}")
print(f"  Empty dictionary entries deleted: {removed_keys}")
print("=" * 50)