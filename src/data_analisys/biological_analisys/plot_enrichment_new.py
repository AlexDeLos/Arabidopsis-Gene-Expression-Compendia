"""
plot_enrichment_new.py
----------------------
Two public plotting functions:

  plot_enrichment_scatter_interactive(...)
      – Self-contained HTML scatter plot of GSEA results, with:
          • Switchable Y-axis (NOM / FDR / FWER p-values)
          • Switchable colour metric
          • Searchable multi-select dropdown to highlight gene sets
          • Navigation buttons that preserve highlighted terms across pages
          • Session-storage restoration on reload

  create_gsea_spider_plot(...)
      – Radar / spider chart comparing GSEA statistics across experiment
        configurations for a single GO term.
"""

import json
import os
import sys
import uuid

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from matplotlib.patches import Circle, RegularPolygon
from matplotlib.path import Path
from matplotlib.projections import register_projection
from matplotlib.projections.polar import PolarAxes
from matplotlib.spines import Spine
from matplotlib.transforms import Affine2D

module_dir = "./"
sys.path.append(module_dir)
from src.constants import FIGURES_DIR, GLOBAL_DIR_PATH  # noqa: E402

matplotlib.rc("font", **{"size": 14})


# =============================================================================
# Helper: Pareto frontier labelling
# =============================================================================

def find_pareto_frontier_indices(df: pd.DataFrame) -> pd.Index:
    """
    Return indices of the top-labelling candidates: top-10 by positive NES,
    top-10 by negative NES, and top-10 by -log10_qval.
    Uses a proxy approach (full Pareto is too slow for large GSEA result sets).
    """
    df_sorted = df.sort_values(by=["-log10_qval", "NES"], ascending=[False, False])
    top_pos = df_sorted[df_sorted["NES"] > 0].head(10).index
    top_neg = df_sorted[df_sorted["NES"] < 0].tail(10).index
    top_qval = df_sorted.head(10).index
    return top_pos.union(top_neg).union(top_qval)


# =============================================================================
# 1. Interactive Scatter Plot
# =============================================================================

def plot_enrichment_scatter_interactive(
    enrichment_df: pd.DataFrame,
    title: str = "Gene Set Enrichment Analysis",
    save_path: str = "interactive_plot.html",
    treatments: list | None = None,
    normalizations: list | None = None,
) -> None:
    """
    Generate a self-contained interactive HTML scatter plot from GSEA results.

    Features
    --------
    - Y-axis can be toggled between NOM p-val, FDR q-val, FWER p-val.
    - Marker colour can be toggled between the same three metrics.
    - Searchable multi-select dropdown highlights chosen gene sets.
    - Navigation buttons link to sibling plots (other tissue / normalization /
      stress / threshold / purity combinations) while preserving the selection.
    - Selections are stored in sessionStorage so they survive page reloads.

    Parameters
    ----------
    enrichment_df : pd.DataFrame
        Must contain: Term, NES, NOM p-val, FDR q-val, FWER p-val.
    title : str
        Plot title shown in the HTML.
    save_path : str
        Output HTML file path.  Parent directories are created as needed.
        The path structure
        ``.../plots_enrichment/{version}/{full|sanity}/{tissue}/{norm}/{thresh}/{purity}/{stress}.html``
        is used to build navigation URLs automatically.
    treatments : list or None
        All treatment names — used to build stress-navigation buttons.
    normalizations : list or None
        All normalisation names — used to build normalisation-navigation buttons.
    """
    if treatments is None:
        treatments = []
    if normalizations is None:
        normalizations = []

    # ------------------------------------------------------------------
    # 1. Data preparation
    # ------------------------------------------------------------------
    df = enrichment_df.copy()

    METRICS_MAP = {
        "NOM p-val":  "-log10(NOM p-val)",
        "FDR q-val":  "-log10(FDR q-val)",
        "FWER p-val": "-log10(FWER p-val)",
    }
    p_eps = 1e-10
    for raw_col, log_col in METRICS_MAP.items():
        if raw_col not in df.columns:
            print(f"  Warning: column '{raw_col}' not found — filling with 0.")
            df[log_col] = 0.0
        else:
            df[log_col] = -np.log10(df[raw_col].astype(float).replace(0, p_eps))

    df["hover_text"] = df.apply(
        lambda row: (
            f"<b>{row['Term']}</b><br><br>"
            f"NES: {row['NES']:.3f}<br>"
            f"NOM p-val: {row['NOM p-val']:.3g}<br>"
            f"FDR q-val: {row['FDR q-val']:.3g}<br>"
            f"FWER p-val: {row['FWER p-val']:.3g}"
        ),
        axis=1,
    )

    # For labelling: rename column so find_pareto_frontier_indices can find it
    pareto_df = df.rename(columns={METRICS_MAP["FDR q-val"]: "-log10_qval"})
    pareto_indices = find_pareto_frontier_indices(pareto_df)
    df_to_label = df.loc[pareto_indices]

    # ------------------------------------------------------------------
    # 2. Plotly figure
    # ------------------------------------------------------------------
    fig = go.Figure()
    plot_div_id = f"plotly-graph-{uuid.uuid4()}"

    initial_y = METRICS_MAP["FDR q-val"]
    initial_color = METRICS_MAP["FWER p-val"]

    fig.add_trace(go.Scatter(        # trace 0 — all points
        x=df["NES"],
        y=df[initial_y],
        mode="markers",
        hoverinfo="text",
        hovertext=df["hover_text"],
        name="Gene Sets",
        marker=dict(
            color=df[initial_color],
            colorscale="Viridis",
            showscale=True,
            colorbar=dict(title=initial_color),
            size=8,
            symbol="circle",
        ),
    ))

    fig.add_trace(go.Scatter(        # trace 1 — labels
        x=df_to_label["NES"],
        y=df_to_label[initial_y],
        mode="text",
        text=df_to_label["Term"],
        textposition="top right",
        textfont=dict(size=10, color="#444"),
        hoverinfo="none",
        name="Labels",
        visible=True,
    ))

    fig.add_trace(go.Scatter(        # trace 2 — highlighted points (JS-driven)
        x=[], y=[],
        mode="markers",
        hoverinfo="none",
        name="Selected",
        showlegend=False,
        marker=dict(color="red", size=16, symbol="star",
                    line=dict(width=1, color="black")),
    ))

    # ------------------------------------------------------------------
    # 3. Layout controls
    # ------------------------------------------------------------------
    y_axis_buttons = [
        dict(
            label=f"Y: {key}",
            method="update",
            args=[
                {"y": [df[val], df_to_label[val]]},
                {"yaxis.title.text": val},
                [0, 1],
            ],
        )
        for key, val in METRICS_MAP.items()
    ]

    color_buttons = [
        dict(
            label=f"Color: {key}",
            method="restyle",
            args=[{"marker.color": [df[val]], "marker.colorbar.title.text": val}, [0]],
        )
        for key, val in METRICS_MAP.items()
    ]

    fig.update_layout(
        title=dict(text=f"<b>{title}</b>", x=0.5),
        xaxis_title="Normalized Enrichment Score (NES)",
        yaxis_title=initial_y,
        template="plotly_white",
        height=800,
        hovermode="closest",
        updatemenus=[
            dict(type="buttons", direction="right", active=1,
                 x=0.01, y=1.12, xanchor="left", buttons=y_axis_buttons),
            dict(type="buttons", direction="right", active=2,
                 x=0.35, y=1.12, xanchor="left", buttons=color_buttons),
            dict(type="buttons", direction="right", active=0,
                 x=0.99, y=1.12, xanchor="right",
                 buttons=[
                     dict(label="Show Labels", method="restyle",
                          args=[{"visible": [True, True, True]}, [0, 1, 2]]),
                     dict(label="Hide Labels", method="restyle",
                          args=[{"visible": [True, False, True]}, [0, 1, 2]]),
                 ]),
        ],
    )
    fig.add_vline(x=0, line_width=1, line_dash="dash", line_color="grey")
    fig.add_hline(y=-np.log10(0.05), line_width=1.5, line_dash="dot",
                  line_color="red", annotation_text="p = 0.05",
                  annotation_position="bottom right")

    # ------------------------------------------------------------------
    # 4. JavaScript data payload
    # ------------------------------------------------------------------
    plot_div = fig.to_html(full_html=False, include_plotlyjs="cdn", div_id=plot_div_id)

    coords_cols = ["Term", "NES"] + list(METRICS_MAP.values())
    coords_json = json.dumps(df[coords_cols].to_dict(orient="records"))
    options_html = "".join(
        f'<option value="{t}">{t}</option>' for t in sorted(df["Term"])
    )

    # ------------------------------------------------------------------
    # 5. Navigation URL construction
    # ------------------------------------------------------------------
    path_parts = save_path.split("/")

    dataset_types = ["full", "sanity"]
    tissue_options = ["All-Tissues", "leaf"]
    threshold_options = ["0", "10", "15"]
    purity_options = ["pure", "mixed"]

    # Defaults
    version = "0.0"
    dataset_type = "full"
    tissue = "All-Tissues"
    normalization = normalizations[0] if normalizations else "combat"
    threshold = "0"
    purity = "mixed"
    stress_name = treatments[0] if treatments else "Heat Stress"

    try:
        base_idx = next(
            (i for i, p in enumerate(path_parts) if "plots_enrichment" in p), 2
        )
        if len(path_parts) > base_idx + 1:
            version = path_parts[base_idx + 1]
        if len(path_parts) > base_idx + 2:
            dataset_type = path_parts[base_idx + 2]
        if len(path_parts) > base_idx + 3:
            tissue = path_parts[base_idx + 3]
        if len(path_parts) > base_idx + 4:
            normalization = path_parts[base_idx + 4]
        if len(path_parts) > base_idx + 5:
            threshold = path_parts[base_idx + 5]
        if len(path_parts) > base_idx + 6:
            purity = path_parts[base_idx + 6]
        if path_parts[-1].endswith(".html"):
            stress_name = path_parts[-1].replace(".html", "")
    except (ValueError, IndexError):
        pass

    def _url(v, dt, t, norm, thresh, p, stress):
        enrich_out = path_parts[2] if len(path_parts) > 2 else "plots_enrichment"
        return (
            f"{GLOBAL_DIR_PATH}{FIGURES_DIR.split('.')[1][1:]}"
            f"{enrich_out}/{v}/{dt}/{t}/{norm}/{thresh}/{p}/{stress}.html"
        )

    nav_urls = {
        "dataset_type":  [(dt, _url(version, dt, tissue, normalization, threshold, purity, stress_name))
                          for dt in dataset_types if dt != dataset_type],
        "tissue":        [(t, _url(version, dataset_type, t, normalization, threshold, purity, stress_name))
                          for t in tissue_options if t != tissue],
        "normalization": [(n, _url(version, dataset_type, tissue, n, threshold, purity, stress_name))
                          for n in normalizations if n != normalization],
        "threshold":     [(th, _url(version, dataset_type, tissue, normalization, th, purity, stress_name))
                          for th in threshold_options if th != threshold],
        "purity":        [(p, _url(version, dataset_type, tissue, normalization, threshold, p, stress_name))
                          for p in purity_options if p != purity],
        "stress":        [(s, _url(version, dataset_type, tissue, normalization, threshold, purity, s))
                          for s in treatments if s != stress_name],
    }

    # ------------------------------------------------------------------
    # 6. HTML template
    # ------------------------------------------------------------------
    def _nav_section(title_text: str, div_id: str, key: str) -> str:
        has_buttons = bool(nav_urls.get(key))
        placeholder = (
            "" if has_buttons
            else '<span class="no-buttons">No other options available</span>'
        )
        return f"""
            <div class="nav-section">
                <h3>{title_text}</h3>
                <div class="nav-buttons" id="{div_id}">{placeholder}</div>
            </div>"""

    nav_sections = (
        _nav_section("Navigate by Dataset Type", "dataset-buttons", "dataset_type")
        + _nav_section("Navigate by Tissue", "tissue-buttons", "tissue")
        + _nav_section("Navigate by Normalization", "normalization-buttons", "normalization")
        + _nav_section("Navigate by Threshold", "threshold-buttons", "threshold")
        + _nav_section("Navigate by Purity", "purity-buttons", "purity")
        + _nav_section("Navigate by Stress Type", "stress-buttons", "stress")
    )

    html_template = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8"/>
    <title>{title}</title>
    <link rel="stylesheet"
          href="https://cdn.jsdelivr.net/npm/choices.js/public/assets/styles/choices.min.css"/>
    <script src="https://cdn.jsdelivr.net/npm/choices.js/public/assets/scripts/choices.min.js">
    </script>
    <style>
        body   {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; }}
        .container   {{ max-width: 1400px; margin: 20px auto; text-align: center; }}
        .choices     {{ margin: 0 auto 20px auto; max-width: 80%; text-align: left; }}
        .choices__inner {{ background-color: #f9f9f9; }}
        .nav-buttons {{ margin: 10px 0; display: flex; justify-content: center; gap: 8px; flex-wrap: wrap; }}
        .nav-button  {{ padding: 6px 12px; background-color: #4CAF50; color: white; border: none;
                        border-radius: 4px; cursor: pointer; font-size: 12px; }}
        .nav-button:hover {{ background-color: #45a049; }}
        .nav-section {{ margin: 15px 0; padding: 10px; border: 1px solid #ddd; border-radius: 5px; }}
        .nav-section h3 {{ margin-bottom: 10px; color: #333; font-size: 14px; }}
        .current-info   {{ background-color: #f0f8ff; padding: 10px; border-radius: 5px;
                           margin: 15px 0; font-size: 14px; }}
        .no-buttons     {{ color: #666; font-style: italic; font-size: 12px; }}
        .dropdown-container {{ margin-top: 20px; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="current-info">
            <strong>Current Plot:</strong><br>
            Version: {version} | Dataset: {dataset_type} | Tissue: {tissue}<br>
            Normalization: {normalization} | Threshold: {threshold}
            | Purity: {purity} | Stress: {stress_name}
        </div>

        {nav_sections}

        <div class="dropdown-container">
            <label for="term-select">Search and select gene sets to highlight:</label>
            <select id="term-select" multiple>
                {options_html}
            </select>
        </div>

        {plot_div}
    </div>

    <script>
        const termData   = {coords_json};
        const plotDivId  = '{plot_div_id}';
        const navUrls    = {json.dumps(nav_urls)};
        const coordMap   = new Map(termData.map(item => [item.Term, item]));

        let currentYMetric = "{initial_y}";

        const choices = new Choices('#term-select', {{
            removeItemButton: true,
            searchResultLimit: 150,
            shouldSort: false,
            placeholder: true,
            placeholderValue: 'Type to search...',
        }});

        function createNavButtons(data, containerId) {{
            const container = document.getElementById(containerId);
            if (container.querySelector('.no-buttons')) container.innerHTML = '';
            data.forEach(([label, url]) => {{
                const btn = document.createElement('button');
                btn.className = 'nav-button';
                btn.textContent = label;
                btn.onclick = () => {{
                    const params = new URLSearchParams();
                    choices.getValue(true).forEach(t => params.append('highlight', t));
                    window.location.href = url + (url.includes('?') ? '&' : '?') + params;
                }};
                container.appendChild(btn);
            }});
        }}

        if (navUrls.dataset_type?.length)  createNavButtons(navUrls.dataset_type,  'dataset-buttons');
        if (navUrls.tissue?.length)         createNavButtons(navUrls.tissue,         'tissue-buttons');
        if (navUrls.normalization?.length)  createNavButtons(navUrls.normalization,  'normalization-buttons');
        if (navUrls.threshold?.length)      createNavButtons(navUrls.threshold,      'threshold-buttons');
        if (navUrls.purity?.length)         createNavButtons(navUrls.purity,         'purity-buttons');
        if (navUrls.stress?.length)         createNavButtons(navUrls.stress,         'stress-buttons');

        function updateHighlights() {{
            const selected = choices.getValue(true);
            const xs = [], ys = [];
            selected.forEach(term => {{
                const d = coordMap.get(term);
                if (d) {{ xs.push(d.NES); ys.push(d[currentYMetric]); }}
            }});
            Plotly.restyle(plotDivId, {{x: [xs], y: [ys]}}, [2]);
            sessionStorage.setItem('selectedTerms', JSON.stringify(selected));
        }}

        document.getElementById('term-select').addEventListener('change', updateHighlights);

        // Keep highlights in sync when the Y-axis toggle button is clicked
        document.getElementById(plotDivId).on('plotly_restyle', (data) => {{
            if (data[0]?.['yaxis.title.text']) {{
                currentYMetric = data[0]['yaxis.title.text'];
                updateHighlights();
            }}
        }});

        function restoreSelectedTerms() {{
            const urlParams = new URLSearchParams(window.location.search);
            const fromURL   = urlParams.getAll('highlight');
            const source    = fromURL.length ? fromURL
                            : JSON.parse(sessionStorage.getItem('selectedTerms') || '[]');
            source.forEach(t => choices.setChoiceByValue(t));
            if (source.length) updateHighlights();
        }}
        document.addEventListener('DOMContentLoaded', restoreSelectedTerms);
    </script>
</body>
</html>"""

    # ------------------------------------------------------------------
    # 7. Write file
    # ------------------------------------------------------------------
    if save_path:
        dirpath = os.path.dirname(save_path)
        if dirpath:
            os.makedirs(dirpath, exist_ok=True)
        # os.makedirs(os.path.dirname(save_path), exist_ok=True)
        with open(save_path, "w", encoding="utf-8") as fh:
            fh.write(html_template)
        print(f"Interactive scatter plot saved → {os.path.abspath(save_path)}")


# =============================================================================
# 2. Radar / Spider Plot
# =============================================================================

def _radar_factory(num_vars: int, frame: str = "polygon"):
    """Register a matplotlib 'radar' projection with ``num_vars`` axes."""
    theta = np.linspace(0, 2 * np.pi, num_vars, endpoint=False)

    class RadarTransform(PolarAxes.PolarTransform):
        def transform_path_non_affine(self, path):
            if path._interpolation_steps > 1:
                path = path.interpolated(num_vars)
            return Path(self.transform(path.vertices), path.codes)

    class RadarAxes(PolarAxes):
        name = "radar"
        PolarTransform = RadarTransform

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.set_theta_zero_location("N")

        def fill(self, *args, closed=True, **kwargs):
            return super().fill(closed=closed, *args, **kwargs)

        def plot(self, *args, **kwargs):
            lines = super().plot(*args, **kwargs)
            for line in lines:
                self._close_line(line)
            return lines

        def _close_line(self, line):
            x, y = line.get_data()
            if x[0] != x[-1]:
                line.set_data(np.append(x, x[0]), np.append(y, y[0]))

        def set_varlabels(self, labels):
            self.set_thetagrids(np.degrees(theta), labels)

        def _gen_axes_patch(self):
            if frame == "circle":
                return Circle((0.5, 0.5), 0.5)
            return RegularPolygon((0.5, 0.5), num_vars, radius=0.5, edgecolor="k")

        def _gen_axes_spines(self):
            if frame == "circle":
                return super()._gen_axes_spines()
            spine = Spine(
                axes=self, spine_type="circle",
                path=Path.unit_regular_polygon(num_vars),
            )
            spine.set_transform(
                Affine2D().scale(0.5).translate(0.5, 0.5) + self.transAxes
            )
            return {"polar": spine}

    register_projection(RadarAxes)
    return theta


def create_gsea_spider_plot(df: pd.DataFrame, save_path: str, term: str) -> None:
    """
    Generate and save a radar / spider plot comparing GSEA statistics across
    multiple experiment configurations for a single GO term.

    Each row of `df` represents one experiment configuration.  The 'Name'
    column is used as the legend label for that configuration.

    Statistics plotted: ES, NES, NOM p-val, FDR q-val, FWER p-val.
    P-values are log-transformed (-ln) before plotting.
    All axes are independently normalised to [0, 1] for comparability.
    A dashed threshold line is drawn at p = 0.01 on the p-value axes.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain columns: ES, NES, NOM p-val, FDR q-val, FWER p-val, Name.
    save_path : str
        Output file path (SVG recommended).  Parent dirs are created if needed.
    term : str
        Human-readable GO term / stress name used in the plot title.
    """
    P_VAL_COLS = ["NOM p-val", "FDR q-val", "FWER p-val"]
    OTHER_COLS  = ["ES", "NES"]
    STATS_COLS  = OTHER_COLS + P_VAL_COLS

    try:
        plot_data = df[STATS_COLS].copy().dropna()

        # Log-transform p-values
        for col in P_VAL_COLS:
            plot_data[col] = -np.log(plot_data[col].astype(float) + 1e-10)

        # Per-axis [min, max] ranges with padding; guard against zero range
        ranges: dict[str, tuple[float, float]] = {}

        all_p = plot_data[P_VAL_COLS].values.flatten()
        p_min, p_max = float(all_p.min()), float(all_p.max())
        p_pad = (p_max - p_min) * 0.05 if p_max != p_min else 0.1
        p_range = (p_min - p_pad, p_max + p_pad)
        for col in P_VAL_COLS:
            ranges[col] = p_range

        for col in OTHER_COLS:
            v_min, v_max = float(plot_data[col].min()), float(plot_data[col].max())
            pad = (v_max - v_min) * 0.05 if v_max != v_min else 0.1
            ranges[col] = (v_min - pad, v_max + pad)

        # Normalise to [0, 1]
        scaled = plot_data.copy()
        for col in STATS_COLS:
            lo, hi = ranges[col]
            denom = (hi - lo) if (hi - lo) != 0 else 1.0
            scaled[col] = (plot_data[col] - lo) / denom

        names = df["Name"]

    except Exception as exc:
        print(f"  Error scaling data for spider plot: {exc}")
        return

    num_vars = len(STATS_COLS)
    theta = _radar_factory(num_vars, frame="polygon")

    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection="radar"))

    axis_labels = ["ES", "NES"] + [f"-ln({c})" for c in P_VAL_COLS]
    ax.set_varlabels(axis_labels)
    ax.tick_params(pad=35)

    colors = plt.cm.tab10.colors
    for i, (_, row) in enumerate(scaled.iterrows()):
        vals = row.values.flatten().tolist()
        color = colors[i % len(colors)]
        ax.plot(theta, vals, label=names.iloc[i], linewidth=2, color=color)
        ax.fill(theta, vals, color=color, alpha=0.1)

    # Threshold line at p = 0.01 drawn across p-value axes only
    thresh_val = -np.log(0.01 + 1e-10)
    p_lo, p_hi = ranges[P_VAL_COLS[0]]
    if p_hi != p_lo:
        thresh_norm = (thresh_val - p_lo) / (p_hi - p_lo)
        if 0.0 <= thresh_norm <= 1.0:
            p_indices = [i for i, c in enumerate(STATS_COLS) if c in P_VAL_COLS]
            p_angles  = [theta[i] for i in p_indices]
            p_radii   = [thresh_norm] * len(p_angles)
            line = ax.plot(p_angles, p_radii,
                           color="red", linestyle=":", linewidth=2,
                           label="p = 0.01", zorder=10)[0]
            # Remove the auto-close segment added by RadarAxes.plot
            lx, ly = line.get_data()
            if len(lx) > len(p_angles):
                line.set_data(lx[: len(p_angles)], ly[: len(p_angles)])
            ax.text(p_angles[-1], thresh_norm, "  p=0.01",
                    color="red", ha="left", va="center",
                    fontsize=10, fontweight="bold")

    # Grid annotations
    ax.set_yticklabels([])
    grid_points = [0.0, 0.5, 1.0]
    ax.set_rgrids(grid_points, labels=[], angle=0, color="grey", alpha=0.3)
    for ang, col_name in zip(theta, STATS_COLS):
        lo, hi = ranges[col_name]
        for gp in grid_points:
            real_val = lo + gp * (hi - lo)
            label_r  = 0.12 if gp == 0.0 else gp
            ax.text(ang, label_r, f"{real_val:.2f}",
                    ha="center", va="center", fontsize=11,
                    bbox=dict(facecolor="none", edgecolor="none", pad=1))

    ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.1), fontsize=12)
    plt.title(f"GSEA Results: {term}", size=18, y=1.1)

    dirpath = os.path.dirname(save_path)
    if dirpath:
        os.makedirs(dirpath, exist_ok=True)
    # os.makedirs(os.path.dirname(save_path), exist_ok=True)
    fig.savefig(save_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Spider plot saved → {save_path}")
