import pandas as pd
import os
import sys
import gc
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px

from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.cluster import KMeans, MiniBatchKMeans
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.model_selection import cross_val_score
from scipy.stats import chi2_contingency
from sklearn.linear_model import LinearRegression


module_dir = './'
sys.path.append(module_dir)

from src.constants import *
#TODO: add this
from src.data_analisys.utils.cluster_exploration_utils_2 import (
    prepare_data_structure, align_labels_to_data, 
    run_pca, run_umap, run_tsne
)
from src.data_analisys.utils.cluster_exploration_utils import *


def calculate_asw_batch_within_biology(X_pca, batch_labels, bio_labels):
    scores = []
    bio_labels = np.array(bio_labels)
    batch_labels = np.array(batch_labels)
    
    for bio_class in np.unique(bio_labels):
        if bio_class in ['unspecified', 'unknown', 'None', 'nan']:
            continue 
            
        mask = (bio_labels == bio_class)
        X_sub = X_pca[mask]
        batch_sub = batch_labels[mask]
        
        if len(X_sub) > 2 and len(np.unique(batch_sub)) > 1:
            score = silhouette_score(X_sub, batch_sub)
            scores.append(abs(score))
            
    if not scores:
        return np.nan
        
    return np.mean(scores)


def variance_explained_by_label(X_pca, labels):
    labels_encoded = pd.get_dummies(labels).values
    model = LinearRegression()
    model.fit(labels_encoded, X_pca)
    return model.score(labels_encoded, X_pca)


def plot_interactive_projection(emb, meta_df, title, output_path):
    """
    Generates a highly interactive Plotly HTML file.
    Includes custom JavaScript for keyboard interactions and an 'Unselect All' button.
    """
    fig = go.Figure()
    categories = meta_df.columns.tolist()
    buttons = []
    total_traces = 0
    cat_trace_indices = {}
    
    hover_text = meta_df.apply(lambda row: '<br>'.join([f"<b>{c.capitalize()}</b>: {row[c]}" for c in categories]), axis=1).tolist()
    color_palette = px.colors.qualitative.Alphabet + px.colors.qualitative.Light24 + px.colors.qualitative.Dark24
    
    for cat in categories:
        unique_vals = meta_df[cat].fillna('unspecified').astype(str).unique()
        cat_trace_indices[cat] = []
        
        for i, val in enumerate(unique_vals):
            mask = (meta_df[cat].fillna('unspecified').astype(str) == val).values
            
            # Pass the raw row data to JS via customdata
            custom_data = meta_df.iloc[mask].astype(str).values
            
            fig.add_trace(go.Scatter(
                x=emb[mask, 0], y=emb[mask, 1],
                mode='markers',
                name=str(val),
                legendgroup=str(val),
                marker=dict(color=color_palette[i % len(color_palette)], size=6, opacity=0.8),
                hovertext=[hover_text[j] for j in range(len(mask)) if mask[j]],
                hoverinfo="text",
                customdata=custom_data,
                visible=(cat == categories[0]) 
            ))
            cat_trace_indices[cat].append(total_traces)
            total_traces += 1
            
    for cat in categories:
        visibility = [False] * total_traces
        for idx in cat_trace_indices[cat]:
            visibility[idx] = True
            
        buttons.append(dict(
            label=f"Color by {cat.capitalize()}",
            method="update",
            args=[{"visible": visibility}, {"title": f"{title} (Colored by {cat.capitalize()})"}]
        ))
        
    fig.update_layout(
        updatemenus=[dict(
            active=0,
            buttons=buttons,
            x=1.02, y=1.1,
            xanchor="left", yanchor="top",
            showactive=True
        )],
        title=f"{title} (Colored by {categories[0].capitalize()})",
        legend_title_text="Legend (Double click to isolate)",
        margin=dict(r=200, t=100),
        width=1200, height=800,
        template="plotly_white"
    )
    
    # Save standard plot
    fig.write_html(output_path, include_plotlyjs='cdn')
    
    # Dump categories list to JS so it dynamically maps the column index
    col_names_json = json.dumps(categories)
    
    # Inject Custom Javascript
    js_code = f"""
    <script>
    document.addEventListener('DOMContentLoaded', function() {{
        var plotDivs = document.getElementsByClassName('plotly-graph-div');
        if (plotDivs.length === 0) return;
        var myPlot = plotDivs[0];
        
        // 1. Unselect All Button
        var btn = document.createElement('button');
        btn.innerHTML = 'Unselect All (Active Category)';
        btn.style.position = 'absolute';
        btn.style.top = '10px';
        btn.style.left = '10px';
        btn.style.zIndex = 1000;
        btn.style.padding = '8px';
        btn.style.backgroundColor = '#f8f9fa';
        btn.style.border = '1px solid #ced4da';
        btn.style.borderRadius = '4px';
        btn.style.cursor = 'pointer';
        document.body.appendChild(btn);

        btn.onclick = function() {{
            var numTraces = myPlot.data.length;
            var updates = [];
            for (var i = 0; i < numTraces; i++) {{
                if (myPlot.data[i].visible === true || myPlot.data[i].visible === 'legendonly') {{
                    updates.push('legendonly'); 
                }} else {{
                    updates.push(false);
                }}
            }}
            Plotly.restyle(myPlot, {{'visible': updates}});
        }};

        // 2. Add instructions box
        var inst = document.createElement('div');
        inst.innerHTML = '<b>Interactions (Hover over a point):</b><br>' + 
                         '- Press <b>"1"</b>: Isolate by Study<br>' +
                         '- Press <b>"2"</b>: Isolate by Tissue<br>' +
                         '- Press <b>"3"</b>: Isolate by Treatment<br>' +
                         '- Press <b>"4"</b>: Isolate by Medium<br>' +
                         '- Press <b>"S"</b>: Reset all opacities';
        inst.style.position = 'absolute';
        inst.style.bottom = '10px';
        inst.style.left = '10px';
        inst.style.zIndex = 1000;
        inst.style.backgroundColor = 'rgba(255,255,255,0.9)';
        inst.style.padding = '10px';
        inst.style.border = '1px solid #ccc';
        inst.style.borderRadius = '4px';
        document.body.appendChild(inst);

        // Map keys to the correct DataFrame column index
        var colNames = {col_names_json};
        var keyMap = {{
            '1': colNames.indexOf('study_id'),
            '2': colNames.indexOf('tissue'),
            '3': colNames.indexOf('treatment'),
            '4': colNames.indexOf('medium')
        }};

        // 3. Hover logic
        var currentHover = null;
        myPlot.on('plotly_hover', function(data){{
            if (data.points.length > 0) {{
                currentHover = data.points[0].customdata;
            }}
        }});
        myPlot.on('plotly_unhover', function(data){{
            currentHover = null;
        }});

        // 4. Keyboard Event Listeners
        document.addEventListener('keydown', function(event) {{
            var key = event.key.toLowerCase();
            
            // 'S' or 'R' Key to Reset all points
            if (key === 's' || key === 'r') {{
                var numTraces = myPlot.data.length;
                var resetOpacities = [];
                for (var i = 0; i < numTraces; i++) {{
                    resetOpacities.push(0.8);
                }}
                Plotly.restyle(myPlot, {{'marker.opacity': resetOpacities}});
                return;
            }}

            // 1, 2, 3, or 4 to isolate by specific metadata
            if (['1', '2', '3', '4'].includes(key)) {{
                if (!currentHover) return; // Must be hovering over a point
                
                var targetColIdx = keyMap[key];
                if (targetColIdx === -1) return; // Column not found in metadata

                var refVal = String(currentHover[targetColIdx]).toLowerCase();
                
                // Do not isolate if the label is 'unknown' or 'unspecified'
                if (refVal === 'unspecified' || refVal === 'unknown' || refVal === 'none' || refVal === 'nan') {{
                    return; 
                }}

                var numTraces = myPlot.data.length;
                var opacities = [];
                for (var i = 0; i < numTraces; i++) {{
                    var traceOpacities = [];
                    var traceData = myPlot.data[i].customdata;
                    if (traceData) {{
                        for (var j = 0; j < traceData.length; j++) {{
                            var val = String(traceData[j][targetColIdx]).toLowerCase();
                            var match = (val === refVal);
                            traceOpacities.push(match ? 1.0 : 0.02);
                        }}
                    }}
                    opacities.push(traceOpacities);
                }}
                Plotly.restyle(myPlot, {{'marker.opacity': opacities}});
            }}
        }});
    }});
    </script>
    """
    
    with open(output_path, 'a') as f:
        f.write(js_code)

def run_exploration_on_dataframe(
    data_df: pd.DataFrame, 
    labels_dict: dict, 
    experiment_name: str,
    output_folder: str
):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    df_aligned = prepare_data_structure(data_df)
    categories = ['treatment', 'tissue', 'medium','study_id']
    
    # Align once to guarantee identical row order across all meta categories
    X_base, _, _, _ = align_labels_to_data(df_aligned, labels_dict, 'study_id')
    
    meta_dict = {}
    for c in categories:
        _, t_labels, _, _ = align_labels_to_data(df_aligned, labels_dict, c)
        meta_dict[c] = t_labels
    meta_df = pd.DataFrame(meta_dict)

    results_summary = []

    for cat in categories:
        print(f"\n[Processing Category Metrics: {cat.upper()}]")
        try:
            # We already know X_base row order perfectly matches meta_df
            text_labels = meta_df[cat].tolist()
            text_labels_np = np.array(text_labels)
            
            n_unspecified = text_labels.count('unspecified')
            
            # TASK 1: METRICS IGNORE 'UNKNOWN' / 'UNSPECIFIED'
            valid_mask = ~np.isin(text_labels_np, ['unknown', 'unspecified', 'None', 'nan'])
            
            X_metric = X_base[valid_mask]
            text_labels_metric = text_labels_np[valid_mask]
            batch_text_labels_metric = meta_df['study_id'].values[valid_mask]
            
            # Encode strings to numeric labels for metrics (like ARI)
            unique_classes, num_labels_metric = np.unique(text_labels_metric, return_inverse=True)
            n_classes_metric = len(unique_classes)
            
            print(f"  -> Total Samples: {X_base.shape[0]} | Dropped Unknowns: {np.sum(~valid_mask)} | Valid for metrics: {X_metric.shape[0]}")

            if X_metric.shape[0] < 5 or n_classes_metric < 2:
                print("  -> Skipping Metrics: Not enough valid samples/classes after removing unknowns.")
                sil_score, ari_score, knn_purity, var_explained, batch_asw = np.nan, np.nan, np.nan, np.nan, np.nan
            else:
                # bulk = get_bulkformer_embeddings(X_metric.T.head())
                X_pca_metric, _ = run_pca(X_metric, n_components=min(50, X_metric.shape[0]-1))
                
                try:
                    sil_score = silhouette_score(X_pca_metric, num_labels_metric, metric='euclidean', sample_size=min(5000, X_metric.shape[0]))
                except ValueError:
                    sil_score = -1 
                
                kmeans = MiniBatchKMeans(n_clusters=n_classes_metric, batch_size=256, random_state=42).fit(X_pca_metric)
                ari_score = adjusted_rand_score(num_labels_metric, kmeans.labels_)
                
                try:
                    n_neighbors = min(5, X_metric.shape[0] - 1)
                    cv_splits = max(2, min(5, np.min(np.bincount(num_labels_metric))))
                    knn = KNeighborsClassifier(n_neighbors=n_neighbors)
                    knn_purity = cross_val_score(knn, X_pca_metric, num_labels_metric, cv=cv_splits).mean()
                except:
                    knn_purity = -1
                
                var_explained = variance_explained_by_label(X_pca_metric, text_labels_metric)
                batch_asw = calculate_asw_batch_within_biology(X_pca_metric, batch_text_labels_metric, text_labels_metric)
                
                print(f"     * Silhouette: {sil_score:.4f} | ARI: {ari_score:.4f} | KNN Purity: {knn_purity:.4f}")
            
            results_summary.append({
                'Category': cat,
                'Silhouette': sil_score,
                'ARI': ari_score,
                'KNN_Purity': knn_purity,
                'Variance_Explained': var_explained,
                'Batch_ASW_within_Bio': batch_asw,
                'Num_Classes': n_classes_metric,
                'Unspecified_Count': n_unspecified
            })

        except Exception as e:
            print(f"  -> Error calculating metrics for {cat}: {e}")

    # TASK 2: Output ONLY 2 interactive HTML files per experiment
    print(f"\n[Generating Global Interactive Plots for {experiment_name}]")
    X_pca_full, _ = run_pca(X_base, n_components=min(50, X_base.shape[0]-1))
    
    print("  -> Running UMAP...")
    umap_emb = run_umap(X_pca_full)
    plot_interactive_projection(
        emb=umap_emb, 
        meta_df=meta_df, 
        title=f'{experiment_name.capitalize()} - UMAP Projection', 
        output_path=f'{output_folder}/{experiment_name}_UMAP.html'
    )
    
    print("  -> Running t-SNE...")
    tsne_emb = run_tsne(X_pca_full)
    plot_interactive_projection(
        emb=tsne_emb, 
        meta_df=meta_df, 
        title=f'{experiment_name.capitalize()} - t-SNE Projection', 
        output_path=f'{output_folder}/{experiment_name}_TSNE.html'
    )

    del X_base, X_pca_full, umap_emb, tsne_emb
    gc.collect()

    res_df = pd.DataFrame(results_summary)
    res_df.to_csv(f'{output_folder}/{experiment_name}_metrics.csv', index=False)
    print("Processing Complete.")
    return res_df


def calculate_cramers_v(series1: pd.Series, series2: pd.Series) -> float:
    contingency_table = pd.crosstab(series1, series2)
    if contingency_table.empty or contingency_table.shape[0] < 2 or contingency_table.shape[1] < 2:
        return 0.0
        
    chi2, _, _, _ = chi2_contingency(contingency_table)
    n = contingency_table.sum().sum()
    min_dim = min(contingency_table.shape) - 1
    
    if min_dim == 0 or n == 0:
        return 0.0
    return np.sqrt(chi2 / (n * min_dim))

def calculate_multilabel_association(study_series: pd.Series, multilabel_series: pd.Series) -> float:
    clean_labels = multilabel_series.apply(
        lambda x: list(x) if isinstance(x, (list, tuple, set)) else ([str(x)] if pd.notna(x) else [])
    ).tolist()
    
    mlb = MultiLabelBinarizer()
    binary_matrix = mlb.fit_transform(clean_labels)
    
    if binary_matrix.shape[1] == 0:
        return 0.0
        
    v_scores = []
    for i, label in enumerate(mlb.classes_):
        binary_col = binary_matrix[:, i]
        if len(np.unique(binary_col)) < 2:
            continue
        v = calculate_cramers_v(study_series, binary_col)
        v_scores.append(v)
        
    if not v_scores:
        return 0.0
    return float(np.mean(v_scores))

def plot_metrics_comparison(metrics_dict: dict, 
                            metadata_df: pd.DataFrame, 
                            output_folder: str, 
                            experiment_name: str = "Normalization_Comparison"):
    
    print(f"\n[Generating Metric Comparison Plots...]")
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    combined_data = []
    for stage_name, df in metrics_dict.items():
        df_copy = df.copy()
        df_copy['Stage'] = stage_name
        combined_data.append(df_copy)
    
    plot_df = pd.concat(combined_data, ignore_index=True)

    bio_targets = LABELS
    confounding_data = []
    
    if 'study_id' in metadata_df.columns:
        for target in bio_targets:
            if target in metadata_df.columns:
                target_data = metadata_df[target]
                has_lists = target_data.apply(lambda x: isinstance(x, (list, tuple, set))).any()
                
                if has_lists:
                    v_score = calculate_multilabel_association(metadata_df['study_id'], target_data)
                elif target_data.nunique() > 1:
                    v_score = calculate_cramers_v(metadata_df['study_id'], target_data)
                else:
                    v_score = 0.0
                    
                confounding_data.append({'Variable': target.capitalize(), 'Cramers_V': v_score})
            else:
                confounding_data.append({'Variable': target.capitalize() + " (Missing)", 'Cramers_V': 0})
                
    confounding_df = pd.DataFrame(confounding_data)

    sns.set_theme(style="whitegrid", context="talk")
    fig, axes = plt.subplots(2, 3, figsize=(24, 12))
    fig.suptitle('Batch Correction Evaluation & Confounding Check\n(Calculated exclusively on valid known labels)', fontsize=20, fontweight='bold', y=0.99)

    sns.barplot(data=plot_df, x='Category', y='Variance_Explained', hue='Stage', ax=axes[0, 0], palette='Set2')
    axes[0, 0].set_title('A. Variance Explained (Higher = More Influence)')
    axes[0, 0].set_ylabel('R² Score')

    sns.barplot(data=plot_df, x='Category', y='KNN_Purity', hue='Stage', ax=axes[0, 1], palette='Set2')
    axes[0, 1].set_title('B. KNN Purity (Higher = Better Local Grouping)')
    axes[0, 1].set_ylabel('Purity Score')
    axes[0, 1].set_ylim(0, 1.1)

    bio_only_df = plot_df[plot_df['Category'] != 'study_id']
    sns.barplot(data=bio_only_df, x='Category', y='Batch_ASW_within_Bio', hue='Stage', ax=axes[0, 2], palette='Set2')
    axes[0, 2].set_title('C. Batch ASW within Bio (Lower = +Study Mixing)')
    axes[0, 2].set_ylabel('Silhouette Score of Batch')

    sns.barplot(data=plot_df, x='Category', y='ARI', hue='Stage', ax=axes[1, 0], palette='Set2')
    axes[1, 0].set_title('D. Adjusted Rand Index (Align. w. clutsers)')
    axes[1, 0].set_ylabel('ARI Score')

    sns.barplot(data=plot_df, x='Category', y='Silhouette', hue='Stage', ax=axes[1, 1], palette='Set2')
    axes[1, 1].set_title('E. Silhouette Score (Higher = Tighter Clusters)')
    axes[1, 1].set_ylabel('Silhouette Score')

    ax_conf = axes[1, 2]
    if not confounding_df.empty:
        sns.barplot(data=confounding_df, x='Variable', y='Cramers_V', ax=ax_conf, color=sns.color_palette("deep")[0])
        ax_conf.set_title('F. Inherent Confounding: Study ID vs. Biology')
        ax_conf.set_ylabel("Association (Cramer's V)\n(1.0 = Perfect Confounding)")
        ax_conf.set_ylim(0, 1.1)
        
        for p in ax_conf.patches:
             ax_conf.annotate(f'{p.get_height():.2f}', 
                              (p.get_x() + p.get_width() / 2., p.get_height()), 
                              ha = 'center', va = 'center', 
                              xytext = (0, 9), textcoords = 'offset points', fontsize=12)
    else:
        ax_conf.text(0.5, 0.5, "Metadata missing or insufficient data\nfor confounding check.", ha='center', va='center')
        ax_conf.set_title('F. Inherent Confounding Check')

    active_axes = axes.flatten()
    for i, ax in enumerate(active_axes):
        ax.set_xlabel('')
        if ax.get_legend() is not None:
            if i == 0: 
                ax.legend(title='Pipeline Stage', loc='upper right', bbox_to_anchor=(1.0, 1.05))
            else:
                ax.get_legend().remove()

    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    
    output_path = os.path.join(output_folder, f"{experiment_name}_Summary_with_Confounding.svg")
    try:
        plt.savefig(output_path, format='svg', bbox_inches='tight')
        plt.savefig(output_path.replace('.svg', '.png'), format='png', bbox_inches='tight', dpi=300)
    except Exception as e:
         print(f"[Warning] Could not save with tight bbox layout: {e}. Saving standard layout.")
         plt.savefig(output_path, format='svg')
         plt.savefig(output_path.replace('.svg', '.png'), format='png', dpi=300)

    print(f"  -> Saved comparison plots to {output_path.replace('.svg', '.png')}")
    plt.close()
        

if __name__ == "__main__":
    #TODO: add a baseline random labels in oreder to compare results to.
    all_metrics = {}
    print(f"Loading labels from: {LABELS_PATH}")
    labels = load_labels_study(LABELS_PATH)

    labels_types = LABELS
    labels_df = make_df_from_labels(labels, labels_types)
    labels_map = labels_df.to_dict() 
    
    del labels,labels_df
    for file in ['filter','imputed','study_corrected','rankin']:
        data_path = f'{STORAGE_DIR}/final_data/{file}.csv' 
        
        if os.path.exists(data_path):
            print(f"Loading expression data from: {data_path}")
            df = pd.read_csv(data_path, index_col=0)
            
            print("  Cleaning sample IDs...")
            df.columns = [c.split('.')[0] for c in df.columns]

            print("  Backfilling missing study_ids using get_study()...")
            if 'study_id' not in labels_map:
                labels_map['study_id'] = {}
                
            count_filled = 0
            for sample in df.columns:
                if sample not in labels_map['study_id']:
                    study_val = get_study(sample)
                    labels_map['study_id'][sample.upper()] = study_val
                    count_filled += 1
            df.columns = [c.upper() for c in df.columns]
            print(f"  -> Added study_id labels for {count_filled} samples.")

            output_dir = f"{CLUSTER_EXPLORATION_FIGURES_DIR}/interactive_plots/{file}"
            
            metrics_df = run_exploration_on_dataframe(
                data_df=df,
                labels_dict=labels_map,
                experiment_name=file,
                output_folder=output_dir
            )
            
            all_metrics[file] = metrics_df
            
        else:
            print(f"Error: Data file not found at {data_path}")
            
    if len(all_metrics) > 1:
        comparison_output_dir = f"{CLUSTER_EXPLORATION_FIGURES_DIR}/interactive_plots/Comparisons"
        plot_metrics_comparison(
            metrics_dict=all_metrics, 
            metadata_df=pd.DataFrame(labels_map),
            output_folder=comparison_output_dir
        )