import pandas as pd
import os
import sys
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.cluster import KMeans, MiniBatchKMeans
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import cross_val_score
from sklearn.linear_model import LinearRegression

module_dir = './'
sys.path.append(module_dir)

from src.constants import *
from src.data_analisys.bulk_transformer import get_bulkformer_embeddings
from src.data_analisys.utils.cluster_exploration_utils_2 import (
    prepare_data_structure, align_labels_to_data, 
    run_pca, run_umap, run_tsne
)
from src.data_analisys.utils.cluster_exploration_utils import *
# Ensure plot_metrics_comparison matches your exact function signature
from src.data_analisys.cluster_explotation import plot_metrics_comparison,calculate_asw_batch_within_biology,variance_explained_by_label


# ==========================================
# --- VISUALIZATION FUNCTIONS ---
# ==========================================
def plot_combined_interactive_projections(embeddings_dict, meta_dicts, title, output_path):
    """
    Generates a single HTML file with side-by-side data matrices. 
    Highlights and legend clicks are linked across all displayed plots.
    Includes keyboard shortcuts (1,2,3,4) to highlight specific categories.
    Responsive full-screen layout.
    """
    stages = list(embeddings_dict.keys())
    num_stages = len(stages)
    
    fig = make_subplots(
        rows=1, cols=num_stages, 
        subplot_titles=stages,
        horizontal_spacing=0.02
    )
    
    first_stage = stages[0]
    categories = list(meta_dicts[first_stage].columns)
    
    trace_visibility_by_cat = {cat: [] for cat in categories}
    colors = px.colors.qualitative.Alphabet + px.colors.qualitative.Plotly
    
    for cat in categories:
        all_classes = set()
        for stage in stages:
            all_classes.update(meta_dicts[stage][cat].astype(str).unique())
        all_classes = sorted(list(all_classes))
        
        color_map = {cls: colors[i % len(colors)] for i, cls in enumerate(all_classes)}
        
        for stage_idx, stage in enumerate(stages):
            emb = embeddings_dict[stage]
            meta = meta_dicts[stage]
            
            for cls in all_classes:
                mask = (meta[cat].astype(str) == cls).values
                
                x_data = emb[mask, 0] if mask.any() else []
                y_data = emb[mask, 1] if mask.any() else []
                text_data = meta.index[mask] if mask.any() and isinstance(meta.index, pd.Index) else []
                
                custom_data = meta[categories].values[mask] if mask.any() else []
                    
                fig.add_trace(
                    go.Scatter(
                        x=x_data,
                        y=y_data,
                        mode='markers',
                        marker=dict(size=4, color=color_map[cls], opacity=0.7),
                        name=str(cls),
                        legendgroup=str(cls),
                        showlegend=(stage_idx == 0),
                        text=text_data,
                        customdata=custom_data, 
                        hoverinfo='text+name'
                    ),
                    row=1, col=stage_idx + 1
                )
                
                for c in categories:
                    trace_visibility_by_cat[c].append(c == cat)

    buttons = []
    for cat in categories:
        buttons.append(dict(
            label=cat.capitalize(),
            method='update',
            args=[{'visible': trace_visibility_by_cat[cat]},
                  {'title': f"{title} - Colored by {cat.capitalize()}"}]
        ))
        
    fig.update_layout(
        updatemenus=[dict(
            active=0,
            buttons=buttons,
            x=0.99, y=0.99, # Moved inside the viewport area
            xanchor='right', yanchor='top',
            bgcolor='rgba(255, 255, 255, 0.9)', # Added semi-transparent background so it is readable over points
            bordercolor='gray',
            borderwidth=1
        )],
        title=dict(text=f"{title} - Colored by {categories[0].capitalize()}", x=0.01, y=0.98),
        hovermode='closest',
        autosize=True, # Allow Plotly to dynamically size itself
        margin=dict(l=10, r=10, t=60, b=10), # Tight margins to maximize plot space
        template="plotly_white"
    )
    
    for i, trace in enumerate(fig.data):
        trace.visible = trace_visibility_by_cat[categories[0]][i] # type: ignore
        
    fig.write_html(output_path, include_plotlyjs='cdn', full_html=True)
    
    # --- CSS & JAVASCRIPT QoL ENHANCEMENTS ---
    enhancements_snippet = """
    <style>
        /* Force the HTML and Body to take exactly the screen height and width without scrolling */
        html, body {
            margin: 0;
            padding: 0;
            width: 100%;
            height: 100%;
            overflow: hidden;
            background-color: #ffffff;
            font-family: Arial, sans-serif;
        }
        /* Make the Plotly container stretch to the very edges */
        .plotly-graph-div {
            width: 100vw !important;
            height: 100vh !important;
        }
    </style>
    <script>
    function bindPlotlyEvents() {
        var graph = document.getElementsByClassName('plotly-graph-div')[0];
        if (!graph || !graph.on) {
            setTimeout(bindPlotlyEvents, 500);
            return;
        }

        // Trigger a resize event to make sure Plotly registers the 100vw/vh CSS
        window.dispatchEvent(new Event('resize'));

        var currentHover = null;
        graph.on('plotly_hover', function(data) {
            currentHover = data.points[0];
        });

        document.addEventListener('keydown', function(event) {
            var key = parseInt(event.key);
            
            if (key >= 1 && key <= 4 && currentHover && currentHover.customdata) {
                var catIndex = key - 1;
                var targetValue = currentHover.customdata[catIndex];
                
                var opacities = [];
                for(var i=0; i<graph.data.length; i++) {
                    var trace = graph.data[i];
                    var traceOpacities = [];
                    
                    if (trace.customdata) {
                        for(var j=0; j<trace.customdata.length; j++) {
                            if (trace.customdata[j][catIndex] === targetValue) {
                                traceOpacities.push(1.0); 
                            } else {
                                traceOpacities.push(0.05); 
                            }
                        }
                    }
                    opacities.push(traceOpacities);
                }
                
                Plotly.restyle(graph, {'marker.opacity': opacities});
            }
            
            if (event.key === 'Escape' || event.key === '0') {
                Plotly.restyle(graph, {'marker.opacity': 0.7});
            }
        });
    }

    window.addEventListener('load', bindPlotlyEvents);
    </script>
    """
    
    with open(output_path, 'a') as f:
        f.write(enhancements_snippet)
        
# ==========================================
# --- CORE PIPELINE FUNCTION ---
# ==========================================

def run_exploration_on_dataframe(
    data_df: pd.DataFrame, 
    labels_dict: dict, 
    experiment_name: str,
    output_folder: str,
    use_bulkformer: bool = False,
    gene_length_dict: dict = None,
    target_vocab: list = None,
    ortholog_map: dict = None
):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    if use_bulkformer:
        print(f"  >>> Extracting BulkFormer Embeddings for {experiment_name}...")
        embeddings_df = get_bulkformer_embeddings(
            count_df=data_df.T, 
            gene_length_dict=gene_length_dict,
            target_vocab=target_vocab,
            ortholog_map=ortholog_map
        )
        df_aligned = embeddings_df.T 
    else:
        print(f"  >>> Using standard PCA preprocessing for {experiment_name}...")
        df_aligned = prepare_data_structure(data_df)

    categories = ['treatment', 'tissue', 'medium','study_id']
    
    X_base, _, _, _ = align_labels_to_data(df_aligned, labels_dict, 'study_id')
    
    meta_df = pd.DataFrame({
        c: align_labels_to_data(df_aligned, labels_dict, c)[1] 
        for c in categories
    })

    results_summary = []

    for cat in categories:
        print(f"\n[Metrics: {cat.upper()}]")
        text_labels_np = np.array(meta_df[cat].tolist())
        valid_mask = ~np.isin(text_labels_np, ['unknown', 'unspecified', 'None', 'nan'])

        X_metric = X_base[valid_mask]
        text_labels_metric = text_labels_np[valid_mask]
        batch_text_labels_metric = meta_df['study_id'].values[valid_mask]

        unique_classes, num_labels_metric = np.unique(text_labels_metric, return_inverse=True)

        if X_metric.shape[0] < 5 or len(unique_classes) < 2:
            print(f"  Not enough valid samples/classes for {cat}.")
            sil_score, ari_score, knn_purity, var_explained, batch_asw = [np.nan] * 5
        else:
            if use_bulkformer:
                X_rep_metric = X_metric
            else:
                X_rep_metric, _ = run_pca(X_metric, n_components=min(50, X_metric.shape[0]-1))

            sil_score = silhouette_score(X_rep_metric, num_labels_metric, sample_size=min(5000, X_rep_metric.shape[0]))
            
            kmeans = MiniBatchKMeans(n_clusters=len(unique_classes), random_state=42).fit(X_rep_metric)
            ari_score = adjusted_rand_score(num_labels_metric, kmeans.labels_)

            knn = KNeighborsClassifier(n_neighbors=min(5, X_rep_metric.shape[0] - 1))
            knn_purity = cross_val_score(knn, X_rep_metric, num_labels_metric, cv=2).mean()

            var_explained = variance_explained_by_label(X_rep_metric, text_labels_metric)
            batch_asw = calculate_asw_batch_within_biology(X_rep_metric, batch_text_labels_metric, text_labels_metric)
            
            print(f"  Silhouette: {sil_score:.3f}, ARI: {ari_score:.3f}, KNN Purity: {knn_purity:.3f}, Var Exp: {var_explained:.3f}, Batch ASW: {batch_asw:.3f}")

        results_summary.append({
            'Category': cat, 'Silhouette': sil_score, 'ARI': ari_score, 
            'KNN_Purity': knn_purity, 'Variance_Explained': var_explained, 
            'Batch_ASW_within_Bio': batch_asw
        })

    print(f"\nGenerating standard UMAP & TSNE for {experiment_name}...")
    
    if use_bulkformer:
        X_rep_full = X_base 
    else:
        X_rep_full, _ = run_pca(X_base, n_components=min(50, X_base.shape[0]-1))
        
    embeddings_out = {}
    
    for method, run_func in [("UMAP", run_umap), ("TSNE", run_tsne)]:
        emb = run_func(X_rep_full)
        embeddings_out[method] = emb
        
    res_df = pd.DataFrame(results_summary)
    res_df.to_csv(f'{output_folder}/{experiment_name}_metrics.csv', index=False)
    
    return res_df, embeddings_out, meta_df


# ==========================================
# --- MAIN EXECUTION BLOCK ---
# ==========================================

if __name__ == "__main__":
    all_metrics = {}
    all_umaps = {}
    all_tsnes = {}
    all_metas = {}
    
    print("Loading Labels Map...")
    labels_map = make_df_from_labels(load_labels_study(LABELS_PATH), LABELS).to_dict() 
    
    stages = ['filter', 'imputed', 'study_corrected', 'rankin']
    
    for file in stages:
        data_path = f'{STORAGE_DIR}/final_data/{file}.csv'
        
        if os.path.exists(data_path):
            print(f"\n{'='*50}\nProcessing {file}\n{'='*50}")
            df = pd.read_csv(data_path, index_col=0)
            
            print("  Cleaning sample IDs...")
            df.columns = [c.split('.')[0].upper() for c in df.columns]

            print("  Backfilling missing study_ids using get_study()...")
            if 'study_id' not in labels_map:
                labels_map['study_id'] = {}
                
            count_filled = 0
            for sample in df.columns:
                if sample not in labels_map['study_id']:
                    study_val = get_study(sample)
                    labels_map['study_id'][sample.upper()] = study_val
                    count_filled += 1
            print(f"  -> Added study_id labels for {count_filled} samples.")

            output_dir = f"{CLUSTER_EXPLORATION_FIGURES_DIR}/interactive_plots_3/{file}"
            
            metrics_df, embeddings, meta_df = run_exploration_on_dataframe(
                data_df=df,
                labels_dict=labels_map,
                experiment_name=file,
                output_folder=output_dir,
                use_bulkformer=False
            )
            
            all_metrics[file] = metrics_df
            all_umaps[file] = embeddings['UMAP']
            all_tsnes[file] = embeddings['TSNE']
            all_metas[file] = meta_df
            
        else:
            print(f"Error: Data file not found at {data_path}")

    # Generate the Comparison Plots 
    if len(all_metrics) > 1:
        comparison_output_dir = f"{CLUSTER_EXPLORATION_FIGURES_DIR}/interactive_plots_2/Comparisons"
        os.makedirs(comparison_output_dir, exist_ok=True)
        
        print("\nGenerating Metric Comparisons...")
        # Ironed out the argument name to properly use 'output_dir' as expected by plot_metrics_comparison
        plot_metrics_comparison(
            metrics_dict=all_metrics, 
            metadata_df=pd.DataFrame(labels_map),
            output_folder=comparison_output_dir
        )
        
        print("Generating linked multi-stage UMAP comparison...")
        plot_combined_interactive_projections(
            embeddings_dict=all_umaps, 
            meta_dicts=all_metas, 
            title="UMAP Cross-Stage Comparison", 
            output_path=f"{comparison_output_dir}/Combined_UMAP.html"
        )
        
        print("Generating linked multi-stage t-SNE comparison...")
        plot_combined_interactive_projections(
            embeddings_dict=all_tsnes, 
            meta_dicts=all_metas, 
            title="t-SNE Cross-Stage Comparison", 
            output_path=f"{comparison_output_dir}/Combined_TSNE.html"
        )