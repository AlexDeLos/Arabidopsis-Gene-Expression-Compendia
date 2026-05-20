import json
from pathlib import Path
from collections import Counter
import sys
import os

sys.path.append(os.path.abspath("./"))
from src.constants import LABELS_PATH

def load_all_samples(labels_dir):
    samples = []
    for path in sorted(Path(labels_dir).glob("*.json")):
        if path.name.startswith("map_"):
            continue
        data = json.loads(path.read_text())
        if isinstance(data, dict):
            samples.extend(data.values())
        elif isinstance(data, list):
            samples.extend(data)
    return [s for s in samples if isinstance(s, dict)]

def extract_combination(sample, labels):
    """Turn one sample into a frozenset-based hashable combo of its label values."""
    combo = {}
    for label in labels:
        raw = sample.get(label, [])
        if not isinstance(raw, list):
            raw = [raw]
        # For treatment, strip intensity so we group by name only
        vals = []
        for item in raw:
            if isinstance(item, dict):
                vals.append(item.get("val") or item.get("value") or "Unspecified")
            else:
                vals.append(str(item) if item else "Unspecified")
        combo[label] = tuple(sorted(vals))
    return tuple(sorted(combo.items()))  # hashable

groupings = ['treatment','tissue']
samples = load_all_samples(LABELS_PATH)
combos = Counter(extract_combination(s, groupings) for s in samples)

print(f"Unique combinations: {len(combos)}")
print(f"Total samples:       {len(samples)}\n")

# Show the top 10 most common combos
for combo, count in combos.most_common(10):
    readable = {k: list(v) for k, v in combo}
    print(f"  {count}x — {readable}")