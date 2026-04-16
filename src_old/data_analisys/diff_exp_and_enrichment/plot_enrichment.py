import pandas as pd
import uuid
import numpy as np
import os
import sys
module_dir = './'
sys.path.append(module_dir)
from src.constants import *


############
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import os
import json

def find_pareto_frontier_indices(df: pd.DataFrame, margin: float = 0.0) -> pd.Index:
    """
    Identifies indices of points on or near the Pareto frontier.
    We aim to maximize: abs(NES), -log10_qval, -log10_FWER_p-val.
    """
    values = df[['NES', '-log10_qval', '-log10_FWER_p-val']].copy()
    values['NES'] = values['NES'].abs()
    points = values.to_numpy()
    
    num_points = points.shape[0]
    is_on_frontier = np.ones(num_points, dtype=bool)

    for i in range(num_points):
        for j in range(num_points):
            if i == j:
                continue
            
            if np.all(points[j] >= points[i] * (1 + margin)) and np.any(points[j] > points[i]):
                is_on_frontier[i] = False
                break
    
    return df.index[is_on_frontier]


##### NEW plot method:
def plot_enrichment_scatter_interactive(enrichment_df: pd.DataFrame, title: str = "Gene Set Enrichment Analysis", save_path: str = "interactive_plot.html",treatments=None,normalizations=None):
    """
    Generates a self-contained, interactive HTML scatter plot from GSEA results.

    This version embeds custom JavaScript to create a searchable, multi-select dropdown
    menu within a single, portable HTML file, plus navigation buttons that change URLs
    while preserving highlighted terms.

    Features:
    - A searchable dropdown to find and highlight multiple gene sets.
    - A button to toggle labels for the most significant terms.
    - Zoom, pan, and hover over points for detailed information.
    - Navigation buttons that change URLs while preserving highlighted terms.

    Args:
        enrichment_df (pd.DataFrame): DataFrame with GSEA results.
        title (str, optional): The title for the plot.
        save_path (str, optional): The path to save the output HTML file.
    """
    # --- 1. Data Preparation ---
    df = enrichment_df.copy()
    df['-log10_qval'] = -np.log10(df['FDR q-val'].astype(float).replace(0, 1e-10))
    df['-log10_FWER_p-val'] = -np.log10(df['FWER p-val'].astype(float).replace(0, 1e-10))
    
    df['hover_text'] = df.apply(
        lambda row: f"""<b>{row['Term']}</b><br><br>
NES: {row['NES']:.3f}<br>
FDR q-val: {row['FDR q-val']:.3g}<br>
FWER p-val: {row['FWER p-val']:.3g}
""",
        axis=1
    )
    
    # Identify significant terms to label
    pareto_indices = find_pareto_frontier_indices(df)
    df_to_label = df.loc[pareto_indices]

    # --- 2. Create the Plotly Figure ---
    fig = go.Figure()
    plot_div_id = f'plotly-graph-{uuid.uuid4()}' # Unique ID for the graph div

    # Trace 0: Main scatter plot
    fig.add_trace(go.Scatter(
        x=df['NES'], y=df['-log10_qval'], mode='markers', hoverinfo='text',
        hovertext=df['hover_text'], name='Gene Sets',
        marker=dict(
            color=df['-log10_FWER_p-val'], colorscale='Viridis', showscale=True,
            colorbar=dict(title="-log10(FWER p-val)"), size=8, symbol='circle'
        )
    ))

    # Trace 1: Labels for significant terms
    fig.add_trace(go.Scatter(
        x=df_to_label['NES'], y=df_to_label['-log10_qval'], mode='text',
        text=df_to_label['Term'], textposition="top right",
        textfont=dict(size=10, color='#444'), hoverinfo='none', name='Labels',
        visible=True # Initially visible
    ))

    # Trace 2: Placeholder for highlighted points (will be controlled by JavaScript)
    fig.add_trace(go.Scatter(
        x=[], y=[], mode='markers', hoverinfo='none', name='Selected', showlegend=False,
        marker=dict(color='red', size=16, symbol='star', line=dict(width=1, color='black'))
    ))
    
    # --- 3. Configure Layout and Controls ---
    fig.update_layout(
        title=dict(text=f'<b>{title}</b>', x=0.5),
        xaxis_title="Normalized Enrichment Score (NES)",
        yaxis_title="-log10(FDR q-value)",
        template='plotly_white', height=700, hovermode='closest',
        updatemenus=[
            dict(
                type="buttons", direction="right", active=0, x=0.9, y=1.1,
                buttons=[
                    dict(label="Show Labels", method="restyle", args=[{"visible": [True, True, True]}, [0, 1, 2]]),
                    dict(label="Hide Labels", method="restyle", args=[{"visible": [True, False, True]}, [0, 1, 2]]),
                ]
            )
        ]
    )
    fig.add_vline(x=0, line_width=1, line_dash="dash", line_color="grey")
    fig.add_hline(y=-np.log10(0.05), line_width=1.5, line_dash="dot", line_color="red",
                  annotation_text="p = 0.05", annotation_position="bottom right")
    fig.add_hline(y=-np.log10(1/2000), line_width=1.5, line_dash="dot", line_color="red",
                  annotation_text=f"p = {1/2000}", annotation_position="bottom right")

    # --- 4. Prepare Data and HTML for Injection ---
    
    # Generate the core Plotly graph HTML
    plot_div = fig.to_html(
        full_html=False, 
        include_plotlyjs='cdn', 
        div_id=plot_div_id
    )

    # Prepare data for JavaScript: a map of Term -> Coordinates
    coords_data = df[['Term', 'NES', '-log10_qval']].to_dict(orient='records')
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
            .container {{ max-width: 1200px; margin: 20px auto; text-align: center; }}
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
            <!-- Current Plot Info -->
            <div class="current-info">
                <strong>Current Plot:</strong><br>
                Version: {version} | Dataset: {dataset_type} | Tissue: {tissue}<br>
                Normalization: {normalization} | Threshold: {threshold} | Purity: {purity} | Stress: {stress_name}
            </div>

            <!-- Navigation Sections -->
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

            <!-- Main Plot -->
            {plot_div}
            
            <label for="term-select">Search and select gene sets to highlight:</label>
            <select id="term-select" multiple>
                {options_html}
            </select>
        </div>

        <script>
            // --- JavaScript to link the dropdown with the Plotly graph ---
            
            // 1. Data passed from Python
            const termData = {coords_json};
            const plotDivId = '{plot_div_id}';
            const navUrls = {json.dumps(nav_urls)};

            // 2. Create a Map for fast coordinate lookups (Term -> {{NES, -log10_qval}})
            const coordMap = new Map(termData.map(item => [item.Term, item]));

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
                // Clear any "no buttons" message
                if (container.querySelector('.no-buttons')) {{
                    container.innerHTML = '';
                }}
                
                buttonData.forEach(([label, url]) => {{
                    const button = document.createElement('button');
                    button.className = 'nav-button';
                    button.textContent = label;
                    button.onclick = function() {{
                        // Get currently selected terms
                        const selectedTerms = choices.getValue(true);
                        
                        // Navigate to new URL with selected terms as URL parameters
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
            if (navUrls.dataset_type && navUrls.dataset_type.length > 0) {{
                createNavButtons(navUrls.dataset_type, 'dataset-buttons');
            }}

            if (navUrls.tissue && navUrls.tissue.length > 0) {{
                createNavButtons(navUrls.tissue, 'tissue-buttons');
            }}

            if (navUrls.normalization && navUrls.normalization.length > 0) {{
                createNavButtons(navUrls.normalization, 'normalization-buttons');
            }}

            if (navUrls.threshold && navUrls.threshold.length > 0) {{
                createNavButtons(navUrls.threshold, 'threshold-buttons');
            }}

            if (navUrls.purity && navUrls.purity.length > 0) {{
                createNavButtons(navUrls.purity, 'purity-buttons');
            }}

            if (navUrls.stress && navUrls.stress.length > 0) {{
                createNavButtons(navUrls.stress, 'stress-buttons');
            }}

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
                        highlightY.push(data['-log10_qval']);
                    }}
                }});

                // Update the 'Selected' trace (index 2) on the Plotly graph
                Plotly.restyle(plotDivId, {{
                    x: [highlightX],
                    y: [highlightY]
                }}, [2]);

                // Store selected terms in sessionStorage to preserve across navigation
                sessionStorage.setItem('selectedTerms', JSON.stringify(selectedTerms));
            }});

            // 7. Restore selected terms from URL parameters or sessionStorage on page load
            function restoreSelectedTerms() {{
                const urlParams = new URLSearchParams(window.location.search);
                const highlightTerms = urlParams.getAll('highlight');
                
                if (highlightTerms.length > 0) {{
                    // Set selections from URL parameters
                    highlightTerms.forEach(term => {{
                        choices.setChoiceByValue(term);
                    }});
                    // Trigger the change event to update the plot
                    document.getElementById('term-select').dispatchEvent(new Event('change'));
                }} else {{
                    // Try to restore from sessionStorage
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

            // Call restore function when page loads
            document.addEventListener('DOMContentLoaded', restoreSelectedTerms);

            // Move the dropdown above the plot for better layout
            const plotContainer = document.getElementById(plotDivId).parentElement;
            const dropdown = document.querySelector('.choices').parentElement;
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