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
    Generates a polished full-screen dark-theme HTML dashboard with side-by-side
    embedding projections. Features:
      - Dark scientific dashboard aesthetic with sidebar controls
      - Category coloring via sidebar buttons (replaces Plotly dropdown)
      - Linked hover-highlighting across all subplots
      - Keyboard shortcuts (1–N) to highlight hovered point's category group
      - Escape / 0 to reset opacity
      - Sample count and category legend in sidebar
      - Responsive, fully fills the browser window
    """
    import json

    stages = list(embeddings_dict.keys())
    num_stages = len(stages)

    first_stage = stages[0]
    categories = list(meta_dicts[first_stage].columns)

    # --- Build a stable, vivid color palette ---
    PALETTE = [
        "#00d4ff", "#ff6b6b", "#51cf66", "#ffd43b", "#cc5de8",
        "#ff922b", "#20c997", "#f06595", "#74c0fc", "#a9e34b",
        "#e599f7", "#66d9e8", "#ffec99", "#ff8787", "#63e6be",
        "#d0bfff", "#ffa94d", "#38d9a9", "#f783ac", "#4dabf7",
        "#ffe066", "#c0eb75", "#e599f7", "#94d82d", "#3bc9db",
    ]

    # Collect all classes per category, build stable color maps
    cat_class_map = {}       # cat -> sorted list of classes
    cat_color_map = {}       # cat -> {cls: color}
    for cat in categories:
        all_classes = set()
        for stage in stages:
            all_classes.update(meta_dicts[stage][cat].astype(str).unique())
        all_classes = sorted(list(all_classes))
        cat_class_map[cat] = all_classes
        cat_color_map[cat] = {cls: PALETTE[i % len(PALETTE)] for i, cls in enumerate(all_classes)}

    # --- Build Plotly figure (one subplot per stage) ---
    fig = make_subplots(
        rows=1, cols=num_stages,
        subplot_titles=[f"<b>{s}</b>" for s in stages],
        horizontal_spacing=0.03,
    )

    trace_visibility_by_cat = {cat: [] for cat in categories}

    for cat in categories:
        all_classes = cat_class_map[cat]
        color_map = cat_color_map[cat]

        for stage_idx, stage in enumerate(stages):
            emb = embeddings_dict[stage]
            meta = meta_dicts[stage]

            for cls in all_classes:
                mask = (meta[cat].astype(str) == cls).values
                x_data = emb[mask, 0] if mask.any() else []
                y_data = emb[mask, 1] if mask.any() else []
                text_data = meta.index[mask].tolist() if mask.any() and isinstance(meta.index, pd.Index) else []
                custom_data = meta[categories].values[mask].tolist() if mask.any() else []

                fig.add_trace(
                    go.Scatter(
                        x=x_data,
                        y=y_data,
                        mode='markers',
                        marker=dict(
                            size=5,
                            color=color_map[cls],
                            opacity=0.82,
                            line=dict(width=0),
                        ),
                        name=str(cls),
                        legendgroup=str(cls),
                        showlegend=False,   # Legend is handled by sidebar HTML
                        text=text_data,
                        customdata=custom_data,
                        hovertemplate=(
                            "<b>%{text}</b><br>"
                            + "<br>".join([f"{c}: %{{customdata[{i}]}}" for i, c in enumerate(categories)])
                            + "<extra></extra>"
                        ),
                    ),
                    row=1, col=stage_idx + 1,
                )

                for c in categories:
                    trace_visibility_by_cat[c].append(c == cat)

    # Dark theme layout — no Plotly dropdown (replaced by sidebar)
    fig.update_layout(
        paper_bgcolor="#0d1117",
        plot_bgcolor="#161b22",
        font=dict(family="'JetBrains Mono', 'Fira Code', monospace", color="#c9d1d9", size=11),
        hovermode='closest',
        autosize=True,
        margin=dict(l=8, r=8, t=42, b=8),
        showlegend=False,
        title=dict(
            text=f"<b>{title}</b>",
            x=0.5, y=0.99,
            xanchor='center', yanchor='top',
            font=dict(size=14, color="#58a6ff"),
        ),
    )

    fig.update_xaxes(
        showgrid=True, gridcolor="#21262d", gridwidth=1,
        zeroline=False, showticklabels=False,
        showline=False,
    )
    fig.update_yaxes(
        showgrid=True, gridcolor="#21262d", gridwidth=1,
        zeroline=False, showticklabels=False,
        showline=False,
    )

    # Style subplot titles
    for ann in fig.layout.annotations:
        ann.font = dict(color="#8b949e", size=12, family="'JetBrains Mono', monospace")
        ann.y = ann.y + 0.01  # type: ignore

    # Set initial visibility to first category
    for i, trace in enumerate(fig.data):
        trace.visible = trace_visibility_by_cat[categories[0]][i]  # type: ignore

    # Serialize data needed by JS sidebar
    js_data = {
        "categories": categories,
        "cat_class_map": cat_class_map,
        "cat_color_map": cat_color_map,
        "trace_visibility": trace_visibility_by_cat,
        "num_traces": len(fig.data),
        "stages": stages,
    }

    fig.write_html(output_path, include_plotlyjs='cdn', full_html=True)

    # --- Inject polished CSS + sidebar + JS ---
    sidebar_html = _build_sidebar_html(js_data, title, categories, cat_class_map, cat_color_map)

    with open(output_path, 'r') as f:
        html_content = f.read()

    # Inject Google Fonts + sidebar CSS before </head>
    font_link = (
        '<link rel="preconnect" href="https://fonts.googleapis.com">'
        '<link href="https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@300;400;600&family=Space+Grotesk:wght@300;400;600&display=swap" rel="stylesheet">'
    )
    html_content = html_content.replace('</head>', font_link + '\n</head>', 1)

    # Wrap existing plotly div in a layout shell, inject sidebar
    shell_open = '<div id="app-shell" style="display:flex;width:100vw;height:100vh;overflow:hidden;background:#0d1117;">'
    shell_close = '</div>'
    sidebar = sidebar_html

    html_content = html_content.replace(
        '<body>',
        f'<body>\n{shell_open}\n{sidebar}\n<div id="plot-area" style="flex:1;min-width:0;height:100vh;">',
        1
    )
    html_content = html_content.replace('</body>', f'</div>\n{shell_close}\n</body>', 1)

    with open(output_path, 'w') as f:
        f.write(html_content)

    # Append JS + global CSS
    with open(output_path, 'a') as f:
        f.write(_build_enhancement_script(js_data, categories))


def _build_sidebar_html(js_data, title, categories, cat_class_map, cat_color_map):
    """Build the HTML for the left sidebar with category buttons and legend."""
    cat_buttons = ""
    for i, cat in enumerate(categories):
        active_cls = "active" if i == 0 else ""
        shortcut = str(i + 1) if i < 9 else ""
        shortcut_badge = f'<span class="shortcut" title="hover + {shortcut} to highlight by this category">{shortcut}</span>' if shortcut else ""
        cat_buttons += (
            f'<button class="cat-btn {active_cls}" data-cat="{cat}" onclick="selectCategory(\'{cat}\')">'
            f'{shortcut_badge}<span class="cat-label">{cat.replace("_", " ").title()}</span>'
            f'</button>\n'
        )

    # Legend for first category
    first_cat = categories[0]
    legend_items = ""
    for cls, color in cat_color_map[first_cat].items():
        count_info = ""  # Will be filled dynamically by JS
        legend_items += (
            f'<div class="legend-item" data-class="{cls}" onclick="highlightClass(\'{cls}\')" title="{cls}">'
            f'<span class="legend-dot" style="background:{color};box-shadow:0 0 6px {color}66;"></span>'
            f'<span class="legend-label">{cls}</span>'
            f'</div>\n'
        )

    return f"""
<div id="sidebar">
  <div id="sidebar-header">
    <div id="sidebar-logo">&#x2022;&#x2022;&#x2022;</div>
    <div id="sidebar-title">Controls</div>
  </div>

  <div class="sidebar-section-label">COLOR BY</div>
  <div id="cat-buttons">
    {cat_buttons}
  </div>

  <div class="sidebar-section-label" style="margin-top:18px;">LEGEND <span id="legend-cat-name"></span></div>
  <div id="legend-scroll">
    <div id="legend-items">
      {legend_items}
    </div>
  </div>

  <div id="sidebar-footer">
    <div class="kbd-hint"><kbd>hover</kbd> + <kbd>1</kbd>–<kbd>N</kbd> Highlight by category</div>
    <div class="kbd-hint"><kbd>Esc</kbd> Reset highlight</div>
  </div>
</div>

<style>
  @import url('https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@300;400;600&family=Space+Grotesk:wght@300;500;600&display=swap');

  *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}

  html, body {{
    width: 100%; height: 100%; overflow: hidden;
    background: #0d1117;
    font-family: 'JetBrains Mono', monospace;
    color: #c9d1d9;
  }}

  #app-shell {{
    display: flex; width: 100vw; height: 100vh; overflow: hidden;
  }}

  #sidebar {{
    width: 210px;
    min-width: 210px;
    height: 100vh;
    background: #010409;
    border-right: 1px solid #21262d;
    display: flex;
    flex-direction: column;
    padding: 0;
    overflow: hidden;
    z-index: 100;
    flex-shrink: 0;
  }}

  #sidebar-header {{
    padding: 16px 14px 12px;
    border-bottom: 1px solid #21262d;
    display: flex;
    align-items: center;
    gap: 10px;
  }}

  #sidebar-logo {{
    font-size: 20px;
    color: #58a6ff;
    letter-spacing: 2px;
    line-height: 1;
  }}

  #sidebar-title {{
    font-family: 'Space Grotesk', sans-serif;
    font-weight: 600;
    font-size: 13px;
    color: #e6edf3;
    letter-spacing: 0.05em;
    text-transform: uppercase;
  }}

  .sidebar-section-label {{
    padding: 12px 14px 6px;
    font-size: 9px;
    font-weight: 600;
    letter-spacing: 0.12em;
    color: #484f58;
    text-transform: uppercase;
    display: flex;
    align-items: center;
    gap: 6px;
  }}

  #legend-cat-name {{
    color: #58a6ff;
    font-size: 9px;
    letter-spacing: 0.08em;
  }}

  #cat-buttons {{
    padding: 4px 10px;
    display: flex;
    flex-direction: column;
    gap: 3px;
  }}

  .cat-btn {{
    display: flex;
    align-items: center;
    gap: 8px;
    width: 100%;
    padding: 7px 10px;
    background: transparent;
    border: 1px solid transparent;
    border-radius: 6px;
    color: #8b949e;
    font-family: 'JetBrains Mono', monospace;
    font-size: 11px;
    cursor: pointer;
    transition: all 0.15s ease;
    text-align: left;
  }}

  .cat-btn:hover {{
    background: #161b22;
    border-color: #30363d;
    color: #c9d1d9;
  }}

  .cat-btn.active {{
    background: #1f2937;
    border-color: #58a6ff44;
    color: #58a6ff;
    font-weight: 600;
  }}

  .shortcut {{
    display: inline-flex;
    align-items: center;
    justify-content: center;
    width: 16px;
    height: 16px;
    background: #21262d;
    border: 1px solid #30363d;
    border-radius: 3px;
    font-size: 9px;
    color: #6e7681;
    flex-shrink: 0;
  }}

  .cat-btn.active .shortcut {{
    background: #1f4068;
    border-color: #58a6ff66;
    color: #58a6ff;
  }}

  /* Shown when this category is currently used as the highlight dimension
     (via keyboard 1-N while hovering), independently of the colour category */
  .cat-btn.highlight-active {{
    border-color: #ffd43b44;
    background: #1f1a0e;
  }}
  .cat-btn.highlight-active .shortcut {{
    background: #3d2e00;
    border-color: #ffd43b66;
    color: #ffd43b;
  }}

  .cat-label {{ flex: 1; }}

  #legend-scroll {{
    flex: 1;
    overflow-y: auto;
    overflow-x: hidden;
    padding: 2px 10px 6px;
    scrollbar-width: thin;
    scrollbar-color: #21262d transparent;
  }}

  #legend-scroll::-webkit-scrollbar {{ width: 4px; }}
  #legend-scroll::-webkit-scrollbar-track {{ background: transparent; }}
  #legend-scroll::-webkit-scrollbar-thumb {{ background: #21262d; border-radius: 2px; }}

  #legend-items {{
    display: flex;
    flex-direction: column;
    gap: 1px;
  }}

  .legend-item {{
    display: flex;
    align-items: center;
    gap: 8px;
    padding: 4px 6px;
    border-radius: 5px;
    cursor: pointer;
    transition: background 0.12s;
    border: 1px solid transparent;
  }}

  .legend-item:hover {{
    background: #161b22;
    border-color: #30363d;
  }}

  .legend-item.dimmed {{ opacity: 0.2; }}
  .legend-item.highlighted {{ border-color: #58a6ff44; background: #1f2937; }}

  .legend-dot {{
    width: 9px; height: 9px;
    border-radius: 50%;
    flex-shrink: 0;
  }}

  .legend-label {{
    font-size: 10px;
    color: #8b949e;
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
    max-width: 148px;
  }}

  .legend-item:hover .legend-label,
  .legend-item.highlighted .legend-label {{ color: #c9d1d9; }}

  #sidebar-footer {{
    padding: 10px 14px 14px;
    border-top: 1px solid #21262d;
    display: flex;
    flex-direction: column;
    gap: 5px;
  }}

  .kbd-hint {{
    font-size: 9px;
    color: #484f58;
    display: flex;
    align-items: center;
    gap: 4px;
  }}

  kbd {{
    display: inline-flex;
    align-items: center;
    padding: 1px 4px;
    background: #161b22;
    border: 1px solid #30363d;
    border-radius: 3px;
    font-family: 'JetBrains Mono', monospace;
    font-size: 8px;
    color: #6e7681;
  }}

  #plot-area {{
    flex: 1;
    min-width: 0;
    height: 100vh;
    position: relative;
  }}

  .plotly-graph-div {{
    width: 100% !important;
    height: 100% !important;
  }}

  /* Style Plotly hover tooltip */
  .hoverlayer .hovertext {{
    font-family: 'JetBrains Mono', monospace !important;
  }}
</style>
"""


def _build_enhancement_script(js_data, categories):
    """Build the JS block that powers sidebar interaction and hover-highlight."""
    import json
    cat_class_map_json = json.dumps(js_data["cat_class_map"])
    cat_color_map_json = json.dumps(js_data["cat_color_map"])
    trace_visibility_json = json.dumps(js_data["trace_visibility"])
    categories_json = json.dumps(categories)

    return f"""
<script>
(function() {{
  var CAT_CLASS_MAP   = {cat_class_map_json};
  var CAT_COLOR_MAP   = {cat_color_map_json};
  var TRACE_VISIBILITY = {trace_visibility_json};
  var CATEGORIES      = {categories_json};

  var currentCat   = CATEGORIES[0];
  var currentHover = null;
  var highlightedClass = null;
  var highlightCatIdx  = null;   // which category dimension the current highlight is in
  var graph        = null;

  function getGraph() {{
    return document.getElementsByClassName('plotly-graph-div')[0];
  }}

  // ---- Category switching ----
  window.selectCategory = function(cat) {{
    currentCat = cat;
    highlightedClass = null;
    highlightCatIdx  = null;

    // Update sidebar active state
    document.querySelectorAll('.cat-btn').forEach(function(btn) {{
      btn.classList.toggle('active', btn.dataset.cat === cat);
      btn.classList.remove('highlight-active');
    }});

    // Rebuild legend
    buildLegend(cat);

    // Update Plotly trace visibility
    var vis = TRACE_VISIBILITY[cat];
    Plotly.restyle(getGraph(), {{ visible: vis }});

    // Reset opacities
    resetOpacity();
  }};

  // ---- Legend builder ----
  function buildLegend(cat) {{
    var container = document.getElementById('legend-items');
    var catNameEl = document.getElementById('legend-cat-name');
    catNameEl.textContent = '· ' + cat.replace(/_/g, ' ').toUpperCase();
    container.innerHTML = '';

    var classes = CAT_CLASS_MAP[cat] || [];
    var colorMap = CAT_COLOR_MAP[cat] || {{}};

    classes.forEach(function(cls) {{
      var color = colorMap[cls] || '#888';
      var item = document.createElement('div');
      item.className = 'legend-item';
      item.dataset.cls = cls;
      item.title = cls;
      item.innerHTML =
        '<span class="legend-dot" style="background:' + color + ';box-shadow:0 0 6px ' + color + '55;"></span>' +
        '<span class="legend-label">' + cls + '</span>';
      item.addEventListener('click', function() {{ window.highlightClass(cls); }});
      container.appendChild(item);
    }});
  }}

  // ---- Highlight a class across all traces (legend click — uses currentCat) ----
  window.highlightClass = function(cls) {{
    var catIdx = CATEGORIES.indexOf(currentCat);

    if (highlightedClass === cls && highlightCatIdx === catIdx) {{
      // Toggle off
      highlightedClass = null;
      highlightCatIdx  = null;
      resetOpacity();
      document.querySelectorAll('.legend-item').forEach(function(el) {{
        el.classList.remove('highlighted', 'dimmed');
      }});
      document.querySelectorAll('.cat-btn').forEach(function(btn) {{
        btn.classList.remove('highlight-active');
      }});
      return;
    }}

    highlightedClass = cls;
    highlightCatIdx  = catIdx;
    applyHighlight(cls, catIdx);

    document.querySelectorAll('.legend-item').forEach(function(el) {{
      var isCls = el.dataset.cls === cls;
      el.classList.toggle('highlighted', isCls);
      el.classList.toggle('dimmed', !isCls);
    }});
    document.querySelectorAll('.cat-btn').forEach(function(btn) {{
      btn.classList.remove('highlight-active');
    }});
  }};

  function applyHighlight(cls, catIdx) {{
    // catIdx: which category column to check in customdata.
    // This is independent of currentCat (the colouring category).
    var g = getGraph();
    var opacities = [];
    for (var i = 0; i < g.data.length; i++) {{
      var trace = g.data[i];
      var traceOp = [];
      if (trace.customdata && trace.customdata.length > 0) {{
        for (var j = 0; j < trace.customdata.length; j++) {{
          var val = trace.customdata[j][catIdx];
          traceOp.push(String(val) === String(cls) ? 1.0 : 0.04);
        }}
      }}
      opacities.push(traceOp);
    }}
    Plotly.restyle(g, {{ 'marker.opacity': opacities }});
  }}

  function resetOpacity() {{
    Plotly.restyle(getGraph(), {{ 'marker.opacity': 0.82 }});
  }}

  // ---- Keyboard shortcuts ----
  // Keys 1–N: while hovering, highlight points sharing the hovered point's
  //           value in category N — WITHOUT changing the active colour category.
  // Esc / 0:  reset all highlights.
  // (H key removed; 1–N now covers that use-case across all categories.)
  document.addEventListener('keydown', function(event) {{
    var key = parseInt(event.key);
    if (!isNaN(key) && key >= 1 && key <= CATEGORIES.length) {{
      // Highlight by category[key-1] using the hovered point's value in that category
      if (currentHover && currentHover.customdata) {{
        var catIdx = key - 1;
        var val = currentHover.customdata[catIdx];
        if (val !== undefined) {{
          var cls = String(val);
          var cat = CATEGORIES[catIdx];

          if (highlightedClass === cls && highlightCatIdx === catIdx) {{
            // Toggle off if same class+category pressed again
            highlightedClass = null;
            highlightCatIdx = null;
            resetOpacity();
            document.querySelectorAll('.legend-item').forEach(function(el) {{
              el.classList.remove('highlighted', 'dimmed');
            }});
          }} else {{
            highlightedClass = cls;
            highlightCatIdx = catIdx;
            applyHighlight(cls, catIdx);

            // Update legend only if the highlight category matches the colour category
            // (legend always reflects the colour category, so we dim/highlight
            //  its items only when they share the same dimension)
            document.querySelectorAll('.legend-item').forEach(function(el) {{
              if (cat === currentCat) {{
                var isCls = el.dataset.cls === cls;
                el.classList.toggle('highlighted', isCls);
                el.classList.toggle('dimmed', !isCls);
              }} else {{
                // Different dimension — show a neutral "highlight active" state
                // so the user sees feedback without misleading legend dimming
                el.classList.remove('highlighted', 'dimmed');
              }}
            }});

            // Flash the category button to give feedback on which dim is active
            document.querySelectorAll('.cat-btn').forEach(function(btn) {{
              btn.classList.toggle('highlight-active', btn.dataset.cat === cat);
            }});
          }}
        }}
      }}
      return;
    }}


    if (event.key === 'Escape' || event.key === '0') {{
      highlightedClass = null;
      highlightCatIdx  = null;
      resetOpacity();
      document.querySelectorAll('.legend-item').forEach(function(el) {{
        el.classList.remove('highlighted', 'dimmed');
      }});
      document.querySelectorAll('.cat-btn').forEach(function(btn) {{
        btn.classList.remove('highlight-active');
      }});
    }}
  }});

  // ---- Init ----
  function init() {{
    graph = getGraph();
    if (!graph || !graph.on) {{
      setTimeout(init, 400);
      return;
    }}

    // Force resize to fill plot-area
    window.dispatchEvent(new Event('resize'));

    graph.on('plotly_hover', function(data) {{
      currentHover = data.points[0];
    }});
    graph.on('plotly_unhover', function() {{
      // Don't clear — keep for keyboard H shortcut
    }});

    // Set legend for first category
    buildLegend(currentCat);
  }}

  window.addEventListener('load', init);
}})();
</script>
"""
        
# ==========================================
# --- CORE PIPELINE FUNCTION ---
# ==========================================

def run_exploration_on_dataframe(
    data_df: pd.DataFrame, 
    labels_dict: dict, 
    experiment_name: str,
    output_folder: str,
    gene_length_dict: dict = None,
    target_vocab: list = None,
    ortholog_map: dict = None
):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

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
                output_folder=output_dir
            )
            
            all_metrics[file] = metrics_df
            all_umaps[file] = embeddings['UMAP']
            all_tsnes[file] = embeddings['TSNE']
            all_metas[file] = meta_df
            
        else:
            print(f"Error: Data file not found at {data_path}")

    # Generate the Comparison Plots 
    if len(all_metrics) > 1:
        comparison_output_dir = f"{CLUSTER_EXPLORATION_FIGURES_DIR}/interactive_plots_3.1/Comparisons"
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