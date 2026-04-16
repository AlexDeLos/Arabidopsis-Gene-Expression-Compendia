import json
import os

from src.constants_labeling import LABELS


def load_json(path: str):
    with open(path) as file:
        return json.load(file)


class LabelMap:
    """
    Manages persistent mapping of raw terms to grounded ontology terms.

    Storage format (v2):
        {
          "raw_term": {
              "label":   "Drought Stress",
              "studies": {
                  "GSE001": ["GSM001", "GSM002"],
                  "GSE002": ["GSM099"]
              }
          }
        }

    Automatically migrates older formats on load:
        v0  "raw_term": "grounded_term"
        v1  "raw_term": ["grounded_term", ["GSE001", "GSE002"]]
    """

    # ------------------------------------------------------------------ #
    #  Init / Loading                                                      #
    # ------------------------------------------------------------------ #

    def __init__(self, path: str | None = None):
        self.path = path

        if path is not None and not os.path.exists(path):
            os.makedirs(path, exist_ok=True)

        for label in LABELS:
            attr_name = f"map_{label}"
            if path is None:
                setattr(self, attr_name, {})
            else:
                file_path = os.path.join(path, f"{attr_name}.json")
                try:
                    raw = load_json(file_path)
                    setattr(self, attr_name, self._migrate(raw))
                except Exception:
                    print(f"Warning: {file_path} not found or invalid. Starting fresh.")
                    setattr(self, attr_name, {})
                    with open(file_path, "w") as f:
                        json.dump({}, f)

    # ------------------------------------------------------------------ #
    #  Migration                                                           #
    # ------------------------------------------------------------------ #

    @staticmethod
    def _migrate(raw: dict) -> dict:
        """Upgrade an entire map dict from v0 / v1 to v2."""
        migrated = {}
        for raw_term, entry in raw.items():
            if isinstance(entry, str):
                # v0: plain string value
                migrated[raw_term] = {"label": entry, "studies": {}}

            elif isinstance(entry, list):
                # v1: ["label", ["GSE001", "GSE002"]]
                label_val = entry[0] if entry else "unknown"
                study_ids = entry[1] if (len(entry) > 1 and isinstance(entry[1], list)) else []
                # No sample IDs available in v1 — carry studies with empty sample lists
                migrated[raw_term] = {
                    "label": label_val,
                    "studies": {sid: [] for sid in study_ids},
                }

            elif isinstance(entry, dict):
                # v2 already — ensure required keys exist
                entry.setdefault("label", "unknown")
                entry.setdefault("studies", {})
                migrated[raw_term] = entry

            else:
                migrated[raw_term] = {"label": str(entry), "studies": {}}

        return migrated

    # ------------------------------------------------------------------ #
    #  Internal helpers                                                    #
    # ------------------------------------------------------------------ #

    def _get_map(self, category: str) -> dict | None:
        return getattr(self, f"map_{category}", None)

    def _update_entry(
        self,
        map_dict: dict,
        raw_label: str,
        grounded_label: str,
        study_id: str,
        sample_id: str | None = None,
    ) -> None:
        """
        Insert or update a mapping entry.
        First-winner policy: the label is never overwritten once set.
        Provenance is tracked as studies[GSE_ID] = [GSM_ID, ...].
        """
        if raw_label not in map_dict:
            map_dict[raw_label] = {"label": grounded_label, "studies": {}}

        entry = map_dict[raw_label]

        # Guard against stale in-memory data from a failed migration
        if not isinstance(entry, dict):
            entry = {"label": grounded_label, "studies": {}}
            map_dict[raw_label] = entry

        entry.setdefault("label", grounded_label)  # first-winner
        entry.setdefault("studies", {})

        if study_id:
            if study_id not in entry["studies"]:
                entry["studies"][study_id] = []
            if sample_id and sample_id not in entry["studies"][study_id]:
                entry["studies"][study_id].append(sample_id)

    # ------------------------------------------------------------------ #
    #  Public read                                                         #
    # ------------------------------------------------------------------ #

    def _get_value(self, category: str, raw_label: str) -> str | None:
        """Return just the canonical label string for a raw term, or None if unseen."""
        map_dict = self._get_map(category)
        if map_dict is None or raw_label not in map_dict:
            return None
        entry = map_dict[raw_label]
        if isinstance(entry, dict):
            return entry.get("label")
        if isinstance(entry, list):  # stale v1 still in memory
            return entry[0]
        return entry  # stale v0 still in memory

    # ------------------------------------------------------------------ #
    #  Public write                                                        #
    # ------------------------------------------------------------------ #

    def add_entry(
        self,
        category: str,
        raw_label: str,
        grounded_label: str,
        study_id: str,
        sample_id: str | None = None,
    ) -> None:
        """Add or update a single entry in the correct label map."""
        map_dict = self._get_map(category)
        if map_dict is not None:
            self._update_entry(map_dict, raw_label, grounded_label, study_id, sample_id)
        else:
            print(f"Warning: Category '{category}' is not a valid label map.")

    def add_mapping_dict(
        self,
        mapping: dict[str, str],
        category: str,
        study_id: str,
        sample_id: str | None = None,
    ) -> None:
        """Batch-update from an external {raw: grounded} dict."""
        map_dict = self._get_map(category)
        if map_dict is not None:
            for raw, grounded in mapping.items():
                self._update_entry(map_dict, raw, grounded, study_id, sample_id)
        else:
            print(f"Warning: Category '{category}' is not a valid label map.")

    def add_mapping(
        self,
        og_sample: dict,
        grounded_sample: dict,
        study_id: str,
        sample_id: str | None = None,
    ) -> None:
        """
        Update maps from a (raw_sample, grounded_sample) pair.
        sample_id is read from grounded_sample['sample_id'] if not supplied explicitly.
        Only 1-to-1 raw->grounded pairs are stored to avoid list-alignment ambiguity.
        """
        if sample_id is None:
            sample_id = grounded_sample.get("sample_id") or og_sample.get("sample_id")

        for category in LABELS:
            raw_val = og_sample.get(category)
            ground_val = grounded_sample.get(category)

            if isinstance(raw_val, set):
                raw_val = list(raw_val)
            if isinstance(ground_val, set):
                ground_val = list(ground_val)
            if isinstance(raw_val, str):
                raw_val = [raw_val]
            if isinstance(ground_val, str):
                ground_val = [ground_val]

            if raw_val and ground_val and len(raw_val) == 1 and len(ground_val) == 1:
                self.add_entry(category, str(raw_val[0]), str(ground_val[0]), study_id, sample_id)

    # ------------------------------------------------------------------ #
    #  Persistence                                                         #
    # ------------------------------------------------------------------ #

    def save_map(self) -> None:
        """Write all label maps to disk."""
        if not self.path:
            return
        for label in LABELS:
            file_path = os.path.join(self.path, f"map_{label}.json")
            map_dict = getattr(self, f"map_{label}", {})
            with open(file_path, "w") as f:
                json.dump(map_dict, f, indent=4)
