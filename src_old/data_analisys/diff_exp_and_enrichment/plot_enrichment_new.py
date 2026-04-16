import pandas as pd
import numpy as np
import plotly.graph_objects as go
import json
import os
import uuid
import sys
module_dir = './'
sys.path.append(module_dir)
from src.constants import *

# --- Helper function (from your original code, required) ---
def find_pareto_frontier_indices(df):
    """
    Finds the indices of the Pareto frontier in a DataFrame.
    Assumes we want to maximize '-log10_qval' and minimize 'NES' (for neg)
    or maximize 'NES' (for pos).
    """
    # A simple proxy: get top N by NES and -log10_qval
    # This is a placeholder. Use your original function.
    df_sorted = df.sort_values(by=['-log10_qval', 'NES'], ascending=[False, False])
    top_pos_nes = df_sorted[df_sorted['NES'] > 0].head(10).index
    top_neg_nes = df_sorted[df_sorted['NES'] < 0].tail(10).index
    top_qval = df_sorted.head(10).index
    return top_pos_nes.union(top_neg_nes).union(top_qval)

# --- Main Function (Improved) ---

def plot_enrichment_scatter_interactive(
    enrichment_df: pd.DataFrame, 
    title: str = "Gene Set Enrichment Analysis", 
    save_path: str = "interactive_plot.html",
    treatments=None,
    normalizations=None,
    global_dir_path: str = "./", # Added for a complete example
    figures_dir: str = "figures/" # Added for a complete example
):
    """
    Generates a self-contained, interactive HTML scatter plot from GSEA results.

    --- NEW FEATURES ---
    - A searchable dropdown to find and highlight multiple gene sets.
    - A button to toggle labels for the most significant terms.
    - Zoom, pan, and hover over points for detailed information.
    - Navigation buttons that change URLs while preserving highlighted terms.
    - Buttons to change the Y-axis metric (NOM p-val, FDR q-val, FWER p-val).
    - Buttons to change the Color metric (NOM p-val, FDR q-val, FWER p-val).

    Args:
        enrichment_df (pd.DataFrame): DataFrame with GSEA results.
                                      MUST contain columns: 'Term', 'NES', 
                                      'NOM p-val', 'FDR q-val', 'FWER p-val'.
        title (str, optional): The title for the plot.
        save_path (str, optional): The path to save the output HTML file.
        ... (other args)
    """
    
    # --- 1. Data Preparation ---
    df = enrichment_df.copy()

    # --- NEW: Define metrics and calculate all -log10 values ---
    # We assume the input DataFrame has these columns
    METRICS_MAP = {
        'NOM p-val': '-log10(NOM p-val)',
        'FDR q-val': '-log10(FDR q-val)',
        'FWER p-val': '-log10(FWER p-val)'
    }
    
    # Use a small epsilon for 0 p-values to avoid -inf
    p_val_epsilon = 1e-10 
    
    for input_col, output_col in METRICS_MAP.items():
        if input_col not in df.columns:
            print(f"Warning: Column '{input_col}' not found. Skipping.")
            # Add a dummy column to avoid errors
            df[output_col] = 0.0
        else:
            df[output_col] = -np.log10(df[input_col].astype(float).replace(0, p_val_epsilon))
    
    # --- MODIFIED: Update hover text to include all p-values ---
    df['hover_text'] = df.apply(
        lambda row: f"""<b>{row['Term']}</b><br><br>
NES: {row['NES']:.3f}<br>
NOM p-val: {row['NOM p-val']:.3g}<br>
FDR q-val: {row['FDR q-val']:.3g}<br>
FWER p-val: {row['FWER p-val']:.3g}
""",
        axis=1
    )
    
    # Identify significant terms to label
    # This will use the default FDR q-val for its calculation
    pareto_indices = find_pareto_frontier_indices(df.rename(columns={METRICS_MAP['FDR q-val']: '-log10_qval'}))
    df_to_label = df.loc[pareto_indices]

    # --- 2. Create the Plotly Figure ---
    fig = go.Figure()
    plot_div_id = f'plotly-graph-{uuid.uuid4()}' # Unique ID for the graph div

    # --- MODIFIED: Set initial state using FDR for Y and FWER for Color ---
    initial_y_metric = METRICS_MAP['FDR q-val']
    initial_color_metric = METRICS_MAP['FWER p-val']

    # Trace 0: Main scatter plot
    fig.add_trace(go.Scatter(
        x=df['NES'], 
        y=df[initial_y_metric],  # Initial Y-axis
        mode='markers', 
        hoverinfo='text',
        hovertext=df['hover_text'], 
        name='Gene Sets',
        marker=dict(
            color=df[initial_color_metric], # Initial Color
            colorscale='Viridis', 
            showscale=True,
            colorbar=dict(title=initial_color_metric), 
            size=8, 
            symbol='circle'
        )
    ))

    # Trace 1: Labels for significant terms
    fig.add_trace(go.Scatter(
        x=df_to_label['NES'], 
        y=df_to_label[initial_y_metric], # Initial Y-axis
        mode='text',
        text=df_to_label['Term'], 
        textposition="top right",
        textfont=dict(size=10, color='#444'), 
        hoverinfo='none', 
        name='Labels',
        visible=True # Initially visible
    ))

    # Trace 2: Placeholder for highlighted points (will be controlled by JavaScript)
    fig.add_trace(go.Scatter(
        x=[], y=[], mode='markers', hoverinfo='none', name='Selected', showlegend=False,
        marker=dict(color='red', size=16, symbol='star', line=dict(width=1, color='black'))
    ))
    
    # --- 3. Configure Layout and Controls ---
    y_axis_buttons = [
        dict(
            label=f"Y: {key}",
            method="update", # <-- CHANGED from 'restyle'
            args=[
                # 1. The data update (for restyle)
                {
                    "y": [df[val], df_to_label[val]]
                },
                # 2. The layout update (for relayout)
                {
                    "yaxis.title.text": val
                },
                # 3. The trace indices for the data update
                [0, 1]
            ]
        ) for key, val in METRICS_MAP.items()
    ]
    
    color_buttons = [
        dict(
            label=f"Color: {key}",
            method="restyle",
            args=[
                {
                    "marker.color": [df[val]],          # Update color data
                    "marker.colorbar.title.text": val # Update colorbar title
                },
                [0] # Target trace 0
            ]
        ) for key, val in METRICS_MAP.items()
    ]

    fig.update_layout(
        title=dict(text=f'<b>{title}</b>', x=0.5),
        xaxis_title="Normalized Enrichment Score (NES)",
        yaxis_title=initial_y_metric, # Initial Y-axis title
        template='plotly_white', height=800, hovermode='closest',
        updatemenus=[
            # --- NEW: Y-Axis Selector ---
            dict(
                type="buttons", direction="right",
                active=1, # Default is FDR (index 1 in METRICS_MAP)
                x=0.01, y=1.1, xanchor='left',
                buttons=y_axis_buttons
            ),
            # --- NEW: Color Selector ---
            dict(
                type="buttons", direction="right",
                active=2, # Default is FWER (index 2 in METRICS_MAP)
                x=0.35, y=1.5, xanchor='left',
                buttons=color_buttons
            ),
            # --- MODIFIED: Original Label Toggle (position adjusted) ---
            dict(
                type="buttons", direction="right", 
                active=0, x=0.99, y=1.5, xanchor='right',
                buttons=[
                    dict(label="Show Labels", method="restyle", args=[{"visible": [True, True, True]}, [0, 1, 2]]),
                    dict(label="Hide Labels", method="restyle", args=[{"visible": [True, False, True]}, [0, 1, 2]]),
                ]
            )
        ]
    )
    fig.add_vline(x=0, line_width=1, line_dash="dash", line_color="grey")
    # This line is now relative to the FDR q-val, which is fine
    fig.add_hline(y=-np.log10(0.05), line_width=1.5, line_dash="dot", line_color="red",
                  annotation_text="p = 0.05", annotation_position="bottom right")

    # --- 4. Prepare Data and HTML for Injection ---
    
    # Generate the core Plotly graph HTML
    plot_div = fig.to_html(
        full_html=False, 
        include_plotlyjs='cdn', 
        div_id=plot_div_id
    )

    # The JS map will now hold NES and all three -log10 p-values
    coords_data_cols = ['Term', 'NES'] + list(METRICS_MAP.values())
    coords_data = df[coords_data_cols].to_dict(orient='records')
    coords_json = json.dumps(coords_data)
    
    # Create HTML <option> elements for the dropdown
    options_html = "".join([f'<option value="{term}">{term}</option>' for term in sorted(df['Term'])])

    # --- 5. Extract path components for navigation ---
    path_parts = save_path.split('/')
    
    # Define possible values for each parameter (you can customize these)
    dataset_types = ['full', 'sanity']
    tissues = ['All-Tissues', 'leaf']  # Add your actual tissues
    thresholds = ['0', '10', '15']
    purities = ['pure', 'mixed']
    stress_names = treatments  # Add your actual stress names
    
    # Initialize with default values
    version = '0.0'
    dataset_type = 'full'
    tissue = 'All-Tissues'
    normalization = '2_way_norm'
    threshold = '14'
    purity = 'pure'
    stress_name = 'Heat'
    
    # Try to extract values from path
    try:
        plots_idx = 2
        if len(path_parts) > plots_idx + 1:
            version = path_parts[plots_idx + 1]
        if len(path_parts) > plots_idx + 2:
            dataset_type = path_parts[plots_idx + 2]
        if len(path_parts) > plots_idx + 3:
            tissue = path_parts[plots_idx + 3]
        if len(path_parts) > plots_idx + 4:
            normalization = path_parts[plots_idx + 4]
        if len(path_parts) > plots_idx + 5:
            threshold = path_parts[plots_idx + 5]
        if len(path_parts) > plots_idx + 6:
            purity = path_parts[plots_idx + 6]
            
        # Extract stress name from filename
        if path_parts[-1].endswith(' Stress.html'):
            stress_name = path_parts[-1].replace('.html', '')
    except (ValueError, IndexError):
        pass  # Keep default values if extraction fails

    # Generate navigation URLs for each category
    nav_urls = {}
    
    # Helper function to generate URL
    def generate_url(v, dt, t, norm, thresh, p, stress):
        enrich_out = save_path.split('/')[2]
        return f"{GLOBAL_DIR_PATH}{FIGURES_DIR.split('.')[1][1:]}{enrich_out}/{v}/{dt}/{t}/{norm}/{thresh}/{p}/{stress}.html"
    
    # Generate URLs for each parameter type
    nav_urls['dataset_type'] = [
        (dt, generate_url(version, dt, tissue, normalization, threshold, purity, stress_name))
        for dt in dataset_types if dt != dataset_type
    ]
    
    nav_urls['tissue'] = [
        (t, generate_url(version, dataset_type, t, normalization, threshold, purity, stress_name))
        for t in tissues if t != tissue
    ]
    
    nav_urls['normalization'] = [
        (norm, generate_url(version, dataset_type, tissue, norm, threshold, purity, stress_name))
        for norm in normalizations if norm != normalization
    ]
    
    nav_urls['threshold'] = [
        (thresh, generate_url(version, dataset_type, tissue, normalization, thresh, purity, stress_name))
        for thresh in thresholds if thresh != threshold
    ]
    
    nav_urls['purity'] = [
        (p, generate_url(version, dataset_type, tissue, normalization, threshold, p, stress_name))
        for p in purities if p != purity
    ]
    
    nav_urls['stress'] = [
        (stress, generate_url(version, dataset_type, tissue, normalization, threshold, purity, stress))
        for stress in stress_names if stress != stress_name
    ]

    # --- 6. Construct the Final HTML Document ---
    html_template = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8" />
        <title>{title}</title>
        <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/choices.js/public/assets/styles/choices.min.css"/>
        <script src="https://cdn.jsdelivr.net/npm/choices.js/public/assets/scripts/choices.min.js"></script>
        <style>
            body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; }}
            .container {{ max-width: 1400px; margin: 20px auto; text-align: center; }}
            .choices {{ margin: 0 auto 20px auto; max-width: 80%; text-align: left; }}
            .choices__inner {{ background-color: #f9f9f9; }}
            .nav-buttons {{ margin: 10px 0; display: flex; justify-content: center; gap: 8px; flex-wrap: wrap; }}
            .nav-button {{ padding: 6px 12px; background-color: #4CAF50; color: white; border: none; border-radius: 4px; cursor: pointer; font-size: 12px; }}
            .nav-button:hover {{ background-color: #45a049; }}
            .nav-section {{ margin: 15px 0; padding: 10px; border: 1px solid #ddd; border-radius: 5px; }}
            .nav-section h3 {{ margin-bottom: 10px; color: #333; font-size: 14px; }}
            .current-info {{ background-color: #f0f8ff; padding: 10px; border-radius: 5px; margin: 15px 0; font-size: 14px; }}
            .no-buttons {{ color: #666; font-style: italic; font-size: 12px; }}
        </style>
    </head>
    <body>
        <div class="container">
            <div class="current-info">
                <strong>Current Plot:</strong><br>
                Version: {version} | Dataset: {dataset_type} | Tissue: {tissue}<br>
                Normalization: {normalization} | Threshold: {threshold} | Purity: {purity} | Stress: {stress_name}
            </div>

            <div class="nav-section">
                <h3>Navigate by Dataset Type</h3>
                <div class="nav-buttons" id="dataset-buttons">
                    {f'<span class="no-buttons">No other dataset types available</span>' if not nav_urls['dataset_type'] else ''}
                </div>
            </div>
            <div class="nav-section">
                <h3>Navigate by Tissue</h3>
                <div class="nav-buttons" id="tissue-buttons">
                    {f'<span class="no-buttons">No other tissues available</span>' if not nav_urls['tissue'] else ''}
                </div>
            </div>
            <div class="nav-section">
                <h3>Navigate by Normalization</h3>
                <div class="nav-buttons" id="normalization-buttons">
                    {f'<span class="no-buttons">No other normalizations available</span>' if not nav_urls['normalization'] else ''}
                </div>
            </div>
            <div class="nav-section">
                <h3>Navigate by Threshold</h3>
                <div class="nav-buttons" id="threshold-buttons">
                    {f'<span class="no-buttons">No other thresholds available</span>' if not nav_urls['threshold'] else ''}
                </div>
            </div>
            <div class="nav-section">
                <h3>Navigate by Purity</h3>
                <div class="nav-buttons" id="purity-buttons">
                    {f'<span class="no-buttons">No other purity types available</span>' if not nav_urls['purity'] else ''}
                </div>
            </div>
            <div class="nav-section">
                <h3>Navigate by Stress Type</h3>
                <div class="nav-buttons" id="stress-buttons">
                    {f'<span class="no-buttons">No other stress types available</span>' if not nav_urls['stress'] else ''}
                </div>
            </div>

            <div class="dropdown-container" style="margin-top: 20px;">
                <label for="term-select">Search and select gene sets to highlight:</label>
                <select id="term-select" multiple>
                    {options_html}
                </select>
            </div>
            
            {plot_div}
            
        </div>

        <script>
            // --- JavaScript to link the dropdown with the Plotly graph ---
            
            // 1. Data passed from Python
            const termData = {coords_json};
            const plotDivId = '{plot_div_id}';
            const navUrls = {json.dumps(nav_urls)};
            
            // --- MODIFIED: Store all coords in the map ---
            // 2. Create a Map for fast coordinate lookups (Term -> {{NES, y1, y2, y3...}})
            const coordMap = new Map(termData.map(item => [item.Term, item]));

            // --- NEW: Keep track of the currently displayed Y-axis metric ---
            let currentYMetric = "{initial_y_metric}"; // Default Y-metric
            const plotDiv = document.getElementById(plotDivId);

            // 3. Initialize the searchable dropdown using Choices.js
            const choices = new Choices('#term-select', {{
                removeItemButton: true,
                searchResultLimit: 150,
                shouldSort: false, // Options are pre-sorted in Python
                placeholder: true,
                placeholderValue: 'Type to search...',
            }});

            // 4. Function to create navigation buttons
            function createNavButtons(buttonData, containerId) {{
                const container = document.getElementById(containerId);
                if (container.querySelector('.no-buttons')) {{
                    container.innerHTML = '';
                }}
                
                buttonData.forEach(([label, url]) => {{
                    const button = document.createElement('button');
                    button.className = 'nav-button';
                    button.textContent = label;
                    button.onclick = function() {{
                        const selectedTerms = choices.getValue(true);
                        const params = new URLSearchParams();
                        selectedTerms.forEach(term => {{
                            params.append('highlight', term);
                        }});
                        const separator = url.includes('?') ? '&' : '?';
                        window.location.href = url + separator + params.toString();
                    }};
                    container.appendChild(button);
                }});
            }}

            // 5. Create all navigation buttons
            if (navUrls.dataset_type && navUrls.dataset_type.length > 0) createNavButtons(navUrls.dataset_type, 'dataset-buttons');
            if (navUrls.tissue && navUrls.tissue.length > 0) createNavButtons(navUrls.tissue, 'tissue-buttons');
            if (navUrls.normalization && navUrls.normalization.length > 0) createNavButtons(navUrls.normalization, 'normalization-buttons');
            if (navUrls.threshold && navUrls.threshold.length > 0) createNavButtons(navUrls.threshold, 'threshold-buttons');
            if (navUrls.purity && navUrls.purity.length > 0) createNavButtons(navUrls.purity, 'purity-buttons');
            if (navUrls.stress && navUrls.stress.length > 0) createNavButtons(navUrls.stress, 'stress-buttons');

            // 6. Listen for changes in the dropdown selection
            document.getElementById('term-select').addEventListener('change', function(event) {{
                const selectedTerms = choices.getValue(true); // Get selected values as an array
                
                const highlightX = [];
                const highlightY = [];

                // Find coordinates for each selected term
                selectedTerms.forEach(term => {{
                    const data = coordMap.get(term);
                    if (data) {{
                        highlightX.push(data.NES);
                        // --- Use the GLOBAL 'currentYMetric' to get the correct Y value ---
                        highlightY.push(data[currentYMetric]); 
                    }}
                }});

                // Update the 'Selected' trace (index 2) on the Plotly graph
                Plotly.restyle(plotDivId, {{
                    x: [highlightX],
                    y: [highlightY]
                }}, [2]);

                sessionStorage.setItem('selectedTerms', JSON.stringify(selectedTerms));
            }});

            // --- NEW: Listen for Plotly button clicks (restyle events) ---
            // This is CRITICAL for keeping the highlights in sync when the Y-axis changes
            plotDiv.on('plotly_restyle', (data) => {{
                // data[0] is the update object, e.g., {{ "yaxis.title.text": "new_title" }}
                // data[1] is an array of trace indices, e.g., [0, 1]
                
                // Check if the y-axis title was part of the update
                if (data[0] && data[0]["yaxis.title.text"]) {{
                    // Update our global variable
                    currentYMetric = data[0]["yaxis.title.text"];
                    console.log('Y-axis metric changed to:', currentYMetric);
                    
                    // The Y-axis has changed, so the highlighted points are now
                    // in the wrong Y position. We must re-calculate and re-draw them.
                    // The easiest way is to re-trigger the dropdown's 'change' event.
                    document.getElementById('term-select').dispatchEvent(new Event('change'));
                }}
            }});

            // 7. Restore selected terms from URL parameters or sessionStorage on page load
            function restoreSelectedTerms() {{
                const urlParams = new URLSearchParams(window.location.search);
                const highlightTerms = urlParams.getAll('highlight');
                
                if (highlightTerms.length > 0) {{
                    highlightTerms.forEach(term => {{
                        choices.setChoiceByValue(term);
                    }});
                    // Trigger the change event to update the plot
                    document.getElementById('term-select').dispatchEvent(new Event('change'));
                }} else {{
                    const storedTerms = sessionStorage.getItem('selectedTerms');
                    if (storedTerms) {{
                        const terms = JSON.parse(storedTerms);
                        terms.forEach(term => {{
                            choices.setChoiceByValue(term);
                        }});
                        document.getElementById('term-select').dispatchEvent(new Event('change'));
                    }}
                }}
            }}
            document.addEventListener('DOMContentLoaded', restoreSelectedTerms);

            // Move the dropdown above the plot
            const plotContainer = plotDiv.parentElement;
            const dropdown = document.querySelector('.dropdown-container');
            plotContainer.parentNode.insertBefore(dropdown, plotContainer);
            
        </script>
    </body>
    </html>
    """

    # --- 7. Save to HTML file ---
    if save_path:
        output_dir = os.path.dirname(save_path)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
        with open(save_path, 'w', encoding='utf-8') as f:
            f.write(html_template)
        print(f"Interactive plot saved to: {os.path.abspath(save_path)}")



import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# --- Imports from matplotlib example ---
from matplotlib.patches import Circle, RegularPolygon
from matplotlib.path import Path
from matplotlib.projections import register_projection
from matplotlib.projections.polar import PolarAxes
from matplotlib.spines import Spine
from matplotlib.transforms import Affine2D

# --- Matplotlib example code for radar_factory ---
# This code is adapted from the matplotlib documentation
# to create a custom 'radar' projection.

def radar_factory(num_vars, frame='circle'):
    """
    Create a radar chart with `num_vars` Axes.

    This function creates a RadarAxes projection and registers it.

    Parameters
    ----------
    num_vars : int
        Number of variables for radar chart.
    frame : {'circle', 'polygon'}
        Shape of frame surrounding Axes.
    """
    # calculate evenly-spaced axis angles
    theta = np.linspace(0, 2*np.pi, num_vars, endpoint=False)

    class RadarTransform(PolarAxes.PolarTransform):

        def transform_path_non_affine(self, path):
            # Paths with non-unit interpolation steps correspond to gridlines,
            # in which case we force interpolation (to defeat PolarTransform's
            # autoconversion to circular arcs).
            if path._interpolation_steps > 1:
                path = path.interpolated(num_vars)
            return Path(self.transform(path.vertices), path.codes)

    class RadarAxes(PolarAxes):

        name = 'radar'
        PolarTransform = RadarTransform

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            # rotate plot such that the first axis is at the top
            self.set_theta_zero_location('N')

        def fill(self, *args, closed=True, **kwargs):
            """Override fill so that line is closed by default"""
            return super().fill(closed=closed, *args, **kwargs)

        def plot(self, *args, **kwargs):
            """Override plot so that line is closed by default"""
            lines = super().plot(*args, **kwargs)
            for line in lines:
                self._close_line(line)

        def _close_line(self, line):
            x, y = line.get_data()
            # FIXME: markers at x[0], y[0] get doubled-up
            if x[0] != x[-1]:
                x = np.append(x, x[0])
                y = np.append(y, y[0])
                line.set_data(x, y)

        def set_varlabels(self, labels):
            self.set_thetagrids(np.degrees(theta), labels)

        def _gen_axes_patch(self):
            # The Axes patch must be centered at (0.5, 0.5) and of radius 0.5
            # in axes coordinates.
            if frame == 'circle':
                return Circle((0.5, 0.5), 0.5)
            elif frame == 'polygon':
                return RegularPolygon((0.5, 0.5), num_vars,
                                      radius=.5, edgecolor="k")
            else:
                raise ValueError("Unknown value for 'frame': %s" % frame)

        def _gen_axes_spines(self):
            if frame == 'circle':
                return super()._gen_axes_spines()
            elif frame == 'polygon':
                # spine_type must be 'left'/'right'/'top'/'bottom'/'circle'.
                spine = Spine(axes=self,
                              spine_type='circle',
                              path=Path.unit_regular_polygon(num_vars))
                # unit_regular_polygon gives a polygon of radius 1 centered at
                # (0, 0) but we want a polygon of radius 0.5 centered at (0.5,
                # 0.5) in axes coordinates.
                spine.set_transform(Affine2D().scale(.5).translate(.5, .5)
                                     + self.transAxes)
                return {'polar': spine}
            else:
                raise ValueError("Unknown value for 'frame': %s" % frame)

    register_projection(RadarAxes)
    return theta


import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import numpy as np
import os
# Source - https://stackoverflow.com/a
# Posted by Herman Schaaf, modified by community. See post 'Timeline' for change history
# Retrieved 2025-11-17, License - CC BY-SA 4.0

font = {'size'   : 14}

matplotlib.rc('font', **font)

# -----------------------------------------------------------------------------
# 2. Main Plotting Function
# -----------------------------------------------------------------------------

def create_gsea_spider_plot(df, save_path, term):
    # ... (Columns setup and initial data processing remains the same) ...
    p_val_cols = ['NOM p-val', 'FDR q-val', 'FWER p-val']
    other_cols = ['ES', 'NES']
    stats_cols = other_cols + p_val_cols
    
    try:
        plot_data = df[stats_cols].copy()
        plot_data = plot_data.dropna()
        # Log transform p-values
        for col in p_val_cols:
            plot_data[col] = -np.log(plot_data[col] + 1e-10)
            # plot_data[col] = plot_data[col].replace(np.nan, 0.1)

        # --- Determine Min/Max Ranges ---
        ranges = {}
        all_p_vals = plot_data[p_val_cols].values.flatten()
        p_min, p_max = all_p_vals.min(), all_p_vals.max()
        
        # 🐛 FIX: Ensure non-zero range for P-values when all are identical
        if p_max == p_min:
            # Enforce a padding of 0.1 if the range is zero. This centers 
            # the constant value at 0.5 (middle) of the scale.
            p_pad_val = 0.1 
            p_range = (p_min - p_pad_val, p_max + p_pad_val)
        else:
            p_pad = (p_max - p_min) * 0.05
            p_range = (p_min - p_pad, p_max + p_pad)
        # END FIX

        for col in p_val_cols:
            ranges[col] = p_range
        # B. Individual scales for ES and NES (already robust)
        for col in other_cols:
            val_min, val_max = plot_data[col].min(), plot_data[col].max()
            # If val_max == val_min, pad is 0.1, creating a non-zero range.
            pad = (val_max - val_min) * 0.05 if val_max != val_min else 0.1
            ranges[col] = (val_min - pad, val_max + pad)

        # --- Normalize Data to [0, 1] ---
        plot_data_scaled = plot_data.copy()
        for col in stats_cols:
            min_v, max_v = ranges[col]
            # Denom will be non-zero thanks to the fix above
            denom = (max_v - min_v) if (max_v - min_v) != 0 else 1.0 
            plot_data_scaled[col] = (plot_data[col] - min_v) / denom

        names = df['Name']
    
    except Exception as e:
        print(f"Error during data scaling: {e}")
        return None, None

    # ... (Plot setup and plotting loop remains the same) ...
    num_vars = len(stats_cols)
    theta = radar_factory(num_vars, frame='polygon')

    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection='radar'))

    labels_display = ['ES', 'NES'] + [f"-ln({c})" for c in p_val_cols]
    ax.set_varlabels(labels_display)
    ax.tick_params(pad=35) 

    colors = plt.cm.tab10.colors 
    for i, row in plot_data_scaled.iterrows():
        values = row.values.flatten().tolist()
        label = names[i]
        color = colors[i % len(colors)]
        
        ax.plot(theta, values, label=label, linewidth=2, color=color)
        ax.fill(theta, values, color=color, alpha=0.1)

    # ---------------------------------------------------------
    # Threshold Line (Including previous fix for open line)
    # ---------------------------------------------------------
    thresh_val = -np.log(0.01 + 1e-10)
    p_min_r, p_max_r = ranges[p_val_cols[0]]
    
    if (p_max_r - p_min_r) != 0:
        thresh_norm = (thresh_val - p_min_r) / (p_max_r - p_min_r)
    else:
        thresh_norm = -1
        
    if 0 <= thresh_norm <= 1.0:
        p_indices = [i for i, col in enumerate(stats_cols) if col in p_val_cols]
        p_angles = [theta[i] for i in p_indices]
        p_radii = [thresh_norm] * len(p_angles)
        
        ax.plot(p_angles, p_radii, color='red', linestyle=':', linewidth=2, label='p=0.01', zorder=10)
        
        line_obj = ax.lines[-1]
        lx, ly = line_obj.get_data()
        if len(lx) > len(p_angles):
            line_obj.set_data(lx[:-1], ly[:-1])
            
        ax.text(p_angles[-1], thresh_norm, ' p=0.01', 
                color='red', ha='left', va='center', fontsize=10, fontweight='bold')
    # ---------------------------------------------------------

    # 5. Custom Grid Labels
    ax.set_yticklabels([]) 
    
    grid_points = [0.0, 0.5, 1.0]
    ax.set_rgrids(grid_points, labels=[], angle=0, color="grey", alpha=0.3)

    # Annotate the values on the axes
    for ang, col_name in zip(theta, stats_cols):
        min_v, max_v = ranges[col_name]
        
        for gp in grid_points:
            real_val = min_v + gp * (max_v - min_v)
            val_str = f"{real_val:.2f}"
            
            label_r = 0.12 if gp == 0.0 else gp

            ax.text(ang, label_r, val_str, 
                    ha='center', va='center', 
                    fontsize=11,  
                    fontweight='normal',
                    color='black', 
                    bbox=dict(facecolor='none', edgecolor='none', pad=1))

    # 6. Add legend and title
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1), fontsize=12)
    plt.title(f'GSEA Results: {term}', size=18, y=1.1)
    
    fig.savefig(f'{save_path}', bbox_inches='tight')
    plt.close()
    return