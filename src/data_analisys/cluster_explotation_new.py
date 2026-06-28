import json
import os
import sys
import random
import itertools
import warnings

import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.model_selection import cross_val_score
from sklearn.neighbors import KNeighborsClassifier

from sklearn.metrics import pairwise_distances

module_dir = "./"
sys.path.append(module_dir)

from src.data_analisys.utils.cluster_exploration_utils_final import get_gsm_id  # noqa: E402
from src.constants import CLUSTER_EXPLORATION_FIGURES_DIR, LABELS_PATH, RNA_USED, SAMPLE_STUDY_MAP, STORAGE_DIR	# noqa: E402
from src.constants_labeling import LABELS as LABEL_AXES	# noqa: E402
from src.data_analisys.utils.cluster_exploration_utils_final import (	# noqa: E402
	align_labels_to_data,
	calculate_asw_batch_within_biology,
	drop_singleton_classes,
	load_labels_study,
	make_df_from_labels,
	plot_metrics_comparison,
	run_bulkformer,
	run_pca,
	run_tsne,
	run_umap,
	multinomial_logistic_accuracy_fun,
	find_n_components_for_variance
)
from src.data_analisys.utils.proxy_distance_metric import run_distance_evaluation	# noqa: E402
from src.data_analisys.utils.plot_distance_metric import plot_distance_metrics, plot_similarity_distance_scatter	# noqa: E402
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

	stages = list(embeddings_dict.keys())
	num_stages = len(stages)

	first_stage = stages[0]
	categories = list(meta_dicts[first_stage].columns)

	# Guard: if any column still contains list/dict values (e.g. treatment was not
	# yet flattened upstream), coerce them to their string representation so the
	# visualiser gets clean flat values for colour-by and legend building.
	def _flatten_cell(v) -> str:
		if isinstance(v, list):
			if not v:
				return "unspecified"
			first = v[0]
			if isinstance(first, dict):
				return str(first.get("val", "unspecified"))
			return str(first)
		if isinstance(v, dict):
			return str(v.get("val", "unspecified"))
		return str(v) if v is not None else "unspecified"

	meta_dicts = {stage: df.apply(lambda col: col.map(_flatten_cell)) for stage, df in meta_dicts.items()}

	# --- Build a stable, vivid color palette ---
	PALETTE = [
		"#00d4ff",
		"#ff6b6b",
		"#51cf66",
		"#ffd43b",
		"#cc5de8",
		"#ff922b",
		"#20c997",
		"#f06595",
		"#74c0fc",
		"#a9e34b",
		"#e599f7",
		"#66d9e8",
		"#ffec99",
		"#ff8787",
		"#63e6be",
		"#d0bfff",
		"#ffa94d",
		"#38d9a9",
		"#f783ac",
		"#4dabf7",
		"#ffe066",
		"#c0eb75",
		"#e599f7",
		"#94d82d",
		"#3bc9db",
	]

	# Collect all classes per category, build stable color maps
	cat_class_map = {}	# cat -> sorted list of classes
	cat_color_map = {}	# cat -> {cls: color}
	for cat in categories:
		all_classes = set()
		for stage in stages:
			all_classes.update(meta_dicts[stage][cat].astype(str).unique())
		all_classes = sorted(all_classes)
		cat_class_map[cat] = all_classes
		cat_color_map[cat] = {cls: PALETTE[i % len(PALETTE)] for i, cls in enumerate(all_classes)}

	# --- Build Plotly figure (one subplot per stage) ---
	fig = make_subplots(
		rows=1,
		cols=num_stages,
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
						mode="markers",
						marker={
							"size": 5,
							"color": color_map[cls],
							"opacity": 0.82,
							"line": {"width": 0},
						},
						name=str(cls),
						legendgroup=str(cls),
						showlegend=False,	# Legend is handled by sidebar HTML
						text=text_data,
						customdata=custom_data,
						hovertemplate=("<b>%{text}</b><br>" + "<br>".join([f"{c}: %{{customdata[{i}]}}" for i, c in enumerate(categories)]) + "<extra></extra>"),
					),
					row=1,
					col=stage_idx + 1,
				)

				for c in categories:
					trace_visibility_by_cat[c].append(c == cat)

	# Dark theme layout — no Plotly dropdown (replaced by sidebar)
	fig.update_layout(
		paper_bgcolor="#0d1117",
		plot_bgcolor="#161b22",
		font={"family": "'JetBrains Mono', 'Fira Code', monospace", "color": "#c9d1d9", "size": 11},
		hovermode="closest",
		autosize=True,
		margin={"l": 8, "r": 8, "t": 42, "b": 8},
		showlegend=False,
		title={
			"text": f"<b>{title}</b>",
			"x": 0.5,
			"y": 0.99,
			"xanchor": "center",
			"yanchor": "top",
			"font": {"size": 14, "color": "#58a6ff"},
		},
	)

	fig.update_xaxes(
		showgrid=True,
		gridcolor="#21262d",
		gridwidth=1,
		zeroline=False,
		showticklabels=False,
		showline=False,
	)
	fig.update_yaxes(
		showgrid=True,
		gridcolor="#21262d",
		gridwidth=1,
		zeroline=False,
		showticklabels=False,
		showline=False,
	)

	# Style subplot titles
	for ann in fig.layout.annotations:
		ann.font = {"color": "#8b949e", "size": 12, "family": "'JetBrains Mono', monospace"}
		ann.y = ann.y + 0.01

	# Set initial visibility to first category
	for i, trace in enumerate(fig.data):
		trace.visible = trace_visibility_by_cat[categories[0]][i]

	# Serialize data needed by JS sidebar
	js_data = {
		"categories": categories,
		"cat_class_map": cat_class_map,
		"cat_color_map": cat_color_map,
		"trace_visibility": trace_visibility_by_cat,
		"num_traces": len(fig.data),
		"stages": stages,
	}

	fig.write_html(output_path, include_plotlyjs="cdn", full_html=True)

	# --- Inject polished CSS + sidebar + JS ---
	sidebar_html = _build_sidebar_html(categories, cat_color_map)

	with open(output_path) as f:
		html_content = f.read()

	# Inject Google Fonts + sidebar CSS before </head>
	font_link = (
		'<link rel="preconnect" href="https://fonts.googleapis.com">'
		'<link href="https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@300;400;600&family=Space+Grotesk:wght@300;400;600&display=swap" rel="stylesheet">'
	)
	html_content = html_content.replace("</head>", font_link + "\n</head>", 1)

	# Wrap existing plotly div in a layout shell, inject sidebar
	shell_open = '<div id="app-shell" style="display:flex;width:100vw;height:100vh;overflow:hidden;background:#0d1117;">'
	shell_close = "</div>"
	sidebar = sidebar_html

	html_content = html_content.replace("<body>", f'<body>\n{shell_open}\n{sidebar}\n<div id="plot-area" style="flex:1;min-width:0;height:100vh;">', 1)
	html_content = html_content.replace("</body>", f"</div>\n{shell_close}\n</body>", 1)

	with open(output_path, "w") as f:
		f.write(html_content)

	# Append JS + global CSS
	with open(output_path, "a") as f:
		f.write(_build_enhancement_script(js_data, categories))


def _build_sidebar_html(categories, cat_color_map):
	"""Build the HTML for the left sidebar with category buttons and legend."""
	cat_buttons = ""
	for i, cat in enumerate(categories):
		active_cls = "active" if i == 0 else ""
		shortcut = str(i + 1) if i < 9 else ""
		shortcut_badge = f'<span class="shortcut" title="hover + {shortcut} to highlight by this category">{shortcut}</span>' if shortcut else ""
		cat_buttons += (
			f'<button class="cat-btn {active_cls}" data-cat="{cat}" onclick="selectCategory(\'{cat}\')">{shortcut_badge}<span class="cat-label">{cat.replace("_", " ").title()}</span></button>\n'
		)

	# Legend for first category
	first_cat = categories[0]
	legend_items = ""
	for cls, color in cat_color_map[first_cat].items():
		legend_items += (
			f'<div class="legend-item" data-class="{cls}" onclick="highlightClass(\'{cls}\')" title="{cls}">'
			f'<span class="legend-dot" style="background:{color};box-shadow:0 0 6px {color}66;"></span>'
			f'<span class="legend-label">{cls}</span>'
			f"</div>\n"
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
	cat_class_map_json = json.dumps(js_data["cat_class_map"])
	cat_color_map_json = json.dumps(js_data["cat_color_map"])
	trace_visibility_json = json.dumps(js_data["trace_visibility"])
	categories_json = json.dumps(categories)

	return f"""
<script>
(function() {{
	var CAT_CLASS_MAP	 = {cat_class_map_json};
	var CAT_COLOR_MAP	 = {cat_color_map_json};
	var TRACE_VISIBILITY = {trace_visibility_json};
	var CATEGORIES		= {categories_json};

	var currentCat	 = CATEGORIES[0];
	var currentHover = null;
	var highlightedClass = null;
	var highlightCatIdx	= null;	 // which category dimension the current highlight is in
	var graph		= null;

	function getGraph() {{
	return document.getElementsByClassName('plotly-graph-div')[0];
	}}

	// ---- Category switching ----
	window.selectCategory = function(cat) {{
	currentCat = cat;
	highlightedClass = null;
	highlightCatIdx	= null;

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
		highlightCatIdx	= null;
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
	highlightCatIdx	= catIdx;
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
	//			 value in category N — WITHOUT changing the active colour category.
	// Esc / 0:	reset all highlights.
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
			//	its items only when they share the same dimension)
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
		highlightCatIdx	= null;
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


def compute_mean_pairwise_distance(
	dist_matrix: np.ndarray,
) -> float:
	"""
	Mean distance across all unique sample pairs.
	"""

	n = dist_matrix.shape[0]

	if n < 2:
		warnings.warn(
			"Cannot compute mean pairwise distance "
			"with fewer than 2 samples."
		)
		return float("nan")

	total_distance = dist_matrix.sum() / 2.0

	n_pairs = n * (n - 1) / 2.0

	return float(total_distance / n_pairs)


# ==========================================
# --- CORE PIPELINE FUNCTION ---
# ==========================================


def run_exploration_on_dataframe(data_df: pd.DataFrame, labels_dict: dict, experiment_name: str, output_folder: str, light_weight: bool = False):
	"""
	data_df is (Sample x Genes)
	"""
	print(f"df for run_exploration_on_dataframe is now {data_df.shape}")
	if not os.path.exists(output_folder):
		os.makedirs(output_folder)

	print(f"	>>> Using standard PCA preprocessing for {experiment_name}...")
	# df_aligned = prepare_data_structure(data_df)
	df_aligned = data_df 
	print(f"MATRIX IS {data_df.shape} with head {data_df.head}")
	# Define X_base for metric calculations
	X_base = df_aligned.values

	# Dynamically grab all axes available in the parsed labels_dict
	available_axes = list(labels_dict.keys())
	if "study_id" not in available_axes:
		available_axes.append("study_id")

	# Build meta_df simply and cleanly
	raw_meta = {}
	for axis in available_axes:
		raw_meta[axis] = align_labels_to_data(df_aligned, labels_dict, axis)

	meta_df = pd.DataFrame(raw_meta, index=df_aligned.index)

	results_summary = []

	# Valid-value exclusions — applies to any label axis
	INVALID_VALUES = {"unknown", "unspecified", "none", "nan", "Unknown", "Unspecified", "None", ""}

	# Calculate metrics for every axis
	metric_categories = available_axes

	for cat in metric_categories:
		print(f"\n[Metrics: {cat.upper()}]")
		text_labels_np = np.array(meta_df[cat].tolist(), dtype=str)
		valid_mask = ~np.isin(text_labels_np, list(INVALID_VALUES))

		X_metric = X_base[valid_mask]
		text_labels_metric = text_labels_np[valid_mask]
		batch_text_labels_metric = np.array(meta_df["study_id"].tolist(), dtype=str)[valid_mask]

		unique_classes = np.unique(text_labels_metric)

		# Convert text labels to numeric codes for sklearn metrics
		num_labels_metric = pd.Series(text_labels_metric).astype("category").cat.codes.values if len(unique_classes) > 0 else []

		if X_metric.shape[0] < 5 or len(unique_classes) < 2:
			print(f"	Not enough valid samples/classes for {cat}.")
			sil_score = ari_score = knn_purity = multinomial_logistic_accuracy = batch_asw = np.nan
		else:
			X_rep_metric = X_metric
			sil_score = silhouette_score(X_rep_metric, num_labels_metric, sample_size=min(5000, X_rep_metric.shape[0]))

			kmeans = MiniBatchKMeans(n_clusters=len(unique_classes), random_state=42, n_init="auto").fit(X_rep_metric)
			ari_score = adjusted_rand_score(num_labels_metric, kmeans.labels_)

			# n_neighbors must fit within the *training fold*, not the full
			# dataset: with cv=2 each fold only trains on ~half the samples.
			# Singleton classes are dropped first so StratifiedKFold doesn't
			# warn and arbitrarily assign them to one fold.
			X_knn, y_knn = drop_singleton_classes(X_rep_metric, num_labels_metric)
			knn_cv_splits = 2
			min_train_fold_size = X_knn.shape[0] - (X_knn.shape[0] // knn_cv_splits)
			knn_n_neighbors = max(1, min(5, min_train_fold_size - 1))
			knn = KNeighborsClassifier(n_neighbors=knn_n_neighbors)
			knn_purity = (
				np.nanmean(cross_val_score(knn, X_knn, y_knn, cv=knn_cv_splits))
				if len(np.unique(y_knn)) >= 2
				else np.nan
			)

			multinomial_logistic_accuracy = multinomial_logistic_accuracy_fun(X_rep_metric, text_labels_metric)
			batch_asw = calculate_asw_batch_within_biology(X_rep_metric, batch_text_labels_metric, text_labels_metric)

			print(f"	Silhouette: {sil_score:.3f}, ARI: {ari_score:.3f}, KNN Purity: {knn_purity:.3f}, multinomial_logistic_accuracy: {multinomial_logistic_accuracy:.3f}, Batch ASW: {batch_asw:.3f}")

		# Format strictly for the plot_metrics_comparison
		for metric_name, val in [("Silhouette", sil_score), ("ARI", ari_score), ("KNN_Purity", knn_purity), ("multinomial_logistic_accuracy", multinomial_logistic_accuracy), ("Batch_ASW_within_Bio", batch_asw)]:
			results_summary.append({"Label_Axis": cat, "Metric": metric_name, "Value": val})
	embeddings_out = {}
	if not light_weight:
		print(f"\nGenerating standard UMAP & TSNE for {experiment_name}...")

		# Run PCA first to feed into UMAP/TSNE
		pca_embedding, _ = run_pca(
			df_aligned,
			n_components=min(
				50,
				df_aligned.shape[0] - 1,
				df_aligned.shape[1] - 1
			)
		)

		embeddings_out["UMAP"] = run_umap(pca_embedding)
		embeddings_out["TSNE"] = run_tsne(pca_embedding)

		print("Generating BulkFormer embedding...")
		embeddings_out["bulk"] = run_bulkformer(
			df_aligned,
			experiment_name
		)
	else:
		print(
			"Skipping UMAP, TSNE and BulkFormer embedding generation "
			"(light_weight=True)."
		)
	
	# ---- CONDITIONAL BULK LATENT SPACE METRICS ---
	bulk_results_summary = []
	if not light_weight:
		print("Generating metric values for Bulk latent space data...")
		for cat in metric_categories:
			print(f"\n[Bulk Metrics: {cat.upper()}]")
			text_labels_np = np.array(meta_df[cat].tolist(), dtype=str)
			valid_mask = ~np.isin(text_labels_np, list(INVALID_VALUES))

			X_metric = embeddings_out["bulk"][valid_mask]
			text_labels_metric = text_labels_np[valid_mask]
			batch_text_labels_metric = np.array(meta_df["study_id"].tolist(), dtype=str)[valid_mask]

			unique_classes = np.unique(text_labels_metric)

			num_labels_metric = pd.Series(text_labels_metric).astype("category").cat.codes.values if len(unique_classes) > 0 else []

			if X_metric.shape[0] < 5 or len(unique_classes) < 2:
				print(f"	Not enough valid samples/classes for {cat}.")
				sil_score = ari_score = knn_purity = multinomial_logistic_accuracy = batch_asw = np.nan
			else:
				X_rep_metric = X_metric
				sil_score = silhouette_score(X_rep_metric, num_labels_metric, sample_size=min(5000, X_rep_metric.shape[0]))

				kmeans = MiniBatchKMeans(n_clusters=len(unique_classes), random_state=42, n_init="auto").fit(X_rep_metric)
				ari_score = adjusted_rand_score(num_labels_metric, kmeans.labels_)

				X_knn, y_knn = drop_singleton_classes(X_rep_metric, num_labels_metric)
				knn_cv_splits = 2
				min_train_fold_size = X_knn.shape[0] - (X_knn.shape[0] // knn_cv_splits)
				knn_n_neighbors = max(1, min(5, min_train_fold_size - 1))
				knn = KNeighborsClassifier(n_neighbors=knn_n_neighbors)
				knn_purity = (
					np.nanmean(cross_val_score(knn, X_knn, y_knn, cv=knn_cv_splits))
					if len(np.unique(y_knn)) >= 2
					else np.nan
				)

				multinomial_logistic_accuracy = multinomial_logistic_accuracy_fun(X_rep_metric, text_labels_metric)
				batch_asw = calculate_asw_batch_within_biology(X_rep_metric, batch_text_labels_metric, text_labels_metric)

				print(f"	Silhouette: {sil_score:.3f}, ARI: {ari_score:.3f}, KNN Purity: {knn_purity:.3f}, Var Exp: {multinomial_logistic_accuracy:.3f}, Batch ASW: {batch_asw:.3f}")

			for metric_name, val in [
				("Silhouette", sil_score),
				("ARI", ari_score),
				("KNN_Purity", knn_purity),
				("multinomial_logistic_accuracy", multinomial_logistic_accuracy),
				("Batch_ASW_within_Bio", batch_asw),
			]:
				bulk_results_summary.append({"Label_Axis": cat, "Metric": metric_name, "Value": val})
	else:
		print("Skipping BulkFormer embedding and metrics generation (light_weight=True).")

	# --- Build and persist both DataFrames separately
	res_df = pd.DataFrame(results_summary)
	res_df.to_csv(f"{output_folder}/{experiment_name}_metrics.csv", index=False)

	bulk_res_df = pd.DataFrame(bulk_results_summary)
	bulk_res_df.to_csv(f"{output_folder}/{experiment_name}_bulk_metrics.csv", index=False)

	return res_df, bulk_res_df, embeddings_out, meta_df

# ==========================================
# --- MAIN EXECUTION BLOCK ---
# ==========================================

if __name__ == "__main__":
	# parser = argparse.ArgumentParser()
	# parser.add_argument("--rna", action="store_true", default=False)
	# args = parser.parse_args()
	N_SAMPLES = 1000
	FULL = True
	LIGHT_WEIGHT = True
	
	all_metrics = {}
	all_bulk_metrics = {}
	all_umaps = {}
	all_tsnes = {}
	all_metas = {}
	all_bulk = {}
	all_dist_metrics = {}
	best_weights_per_stage = {}  # Tracks the optimal weights for each stage

# --- [WEIGHT TESTING CONSTRAINTS] ---
	# Set to 1 or 0 to LOCK the axis to that specific value.
	# Set to None to leave it UNLOCKED (the code will test both 0 and 1).
	WEIGHT_CONSTRAINTS = {
		"tissue": 1,				  # LOCKED: Always active
		"developmental_stage": 1,	 # LOCKED: Always active
		"treatment": None,			# UNLOCKED: Test both 0 and 1
		"ecotype": 1,				 # LOCKED: Always disabled
		"modification": 0,		 # UNLOCKED: Test both 0 and 1
		"medium": 0,			   # UNLOCKED: Test both 0 and 1
		"treatment_intensity": None,  # UNLOCKED: Test both 0 and 1
	}
	
	weight_keys = list(WEIGHT_CONSTRAINTS.keys())
	
	# 1. Dynamically build the search pool per category
	# If locked, list holds 1 element (e.g., [1]). If unlocked, list holds 2 elements ([0, 1]).
	choices = [
		[WEIGHT_CONSTRAINTS[k]] if WEIGHT_CONSTRAINTS[k] is not None else [0, 1]
		for k in weight_keys
	]
	
	# 2. Unpack choices (*) to generate ONLY the permitted combinations
	binary_combinations = []
	for comb in itertools.product(*choices):
		# Create a temporary mapping to perform a safety check
		temp_weights = {weight_keys[i]: comb[i] for i in range(len(weight_keys))}
		if any(temp_weights.values()):  # Prevents testing an all-zero vector
			binary_combinations.append(comb)
			
	print(f"🔒 Applied constraints. Generated {len(binary_combinations)} targeted combinations to test.")

	print("Running a light version of the code")
	print(f"Loading Labels Map from {LABELS_PATH}...")
	labels_map = make_df_from_labels(load_labels_study(LABELS_PATH)).to_dict()
	
	stages = ["filter_norm", "combat_norm", "rankin"] if RNA_USED else ["filter_norm", "combat_norm", "rankin"]
	
	for file in stages:
		data_path = f"{STORAGE_DIR}final_data/rnaseq_processed/{file}.csv" if RNA_USED else f"{STORAGE_DIR}final_data/{file}.csv"

		if os.path.exists(data_path):
			print(f"\n{'=' * 50}\nProcessing {file}\n{'=' * 50}")
			
			# --- [DATA LOADING & PCA PROCESSING (Done ONCE per stage)] ---
			if N_SAMPLES is not None:
				print(f"Loading a random selection of {N_SAMPLES} samples (memory-saving mode)...")
				header = pd.read_csv(data_path, nrows=0)
				total_columns = len(header.columns)
				n_to_sample = min(N_SAMPLES, total_columns - 1)
				random_col_indices = random.sample(range(1, total_columns), n_to_sample)
				final_usecols = [0] + random_col_indices
				df = pd.read_csv(data_path, index_col=0, usecols=final_usecols)
			else:
				df = pd.read_csv(data_path, index_col=0)

			print(" Cleaning sample IDs...")
			df.columns = [c.split(".")[0].upper() for c in df.columns]

			print(" Backfilling missing study_ids using get_study()...")
			if "study_id" not in labels_map:
				labels_map["study_id"] = {}

			count_filled = 0
			for sample in df.columns:
				if sample not in labels_map["study_id"]:
					study_val = str(SAMPLE_STUDY_MAP.at[sample, "StudyID"]) if sample in SAMPLE_STUDY_MAP.index else "Unknown_Study"
					labels_map["study_id"][sample.upper()] = study_val
					count_filled += 1
			print(f"	-> Added study_id labels for {count_filled} samples.")

			output_dir = f"{CLUSTER_EXPLORATION_FIGURES_DIR}/interactive_plots/{file}"
			os.makedirs(output_dir, exist_ok=True)

			df = df.T
			print(f"df is now {df.shape}")
			n_components, cumulative_variance, pca, pca_embedding = find_n_components_for_variance(
				df, variance_threshold=0.90, save_path=output_dir
			)
			pca_embedding = pca_embedding[:, :n_components]
			df = pd.DataFrame(pca_embedding, index=df.index, columns=[f"PC{i + 1}" for i in range(pca_embedding.shape[1])])
			print(f"  Sliced PCA embedding to {pca_embedding.shape[1]} components (90% variance subset)")
			
			if FULL:
				metrics_df, bulk_metrics_df, embeddings, meta_df = run_exploration_on_dataframe(
					data_df=df, labels_dict=labels_map, experiment_name=file, output_folder=output_dir, light_weight=LIGHT_WEIGHT
				)
				all_metrics[file] = metrics_df
				if not LIGHT_WEIGHT:
					all_bulk_metrics[file] = bulk_metrics_df
					all_umaps[file] = embeddings["UMAP"]
					all_tsnes[file] = embeddings["TSNE"]
					all_bulk[file] = embeddings["bulk"]
				all_metas[file] = meta_df

			# --- [COMBINATIONS SWEEP FOR WEIGHTS CORRELATION] ---
			print(f"\n--> Starting grid search over {len(binary_combinations)} binary weight vectors...")
			stage_results = []

			print(f"\n[DistMetrics] Running for stage: '{file}'")
			if RNA_USED:
				df.index = [get_gsm_id(col.split('_')[1]) for col in df.index]
			if RNA_USED:
				try:
					SAMPLE_STUDY_MAP.index = [get_gsm_id(ind.split('_')[1]) for ind in SAMPLE_STUDY_MAP.index]
				except IndexError:
					pass
			nan_count = np.isnan(df).sum()
			inf_count = np.isinf(df).sum()
			print(f"DEBUG: Found {nan_count} NaNs and {inf_count} Infs in the matrix!")
			# --------------------------------------------------
			# PCA distance matrix
			# --------------------------------------------------
			# Samples × Principal_components convention:
			# rows = samples
			# columns = Principal_components
			ordered_samples = list(df.index)

			dist_matrix = pairwise_distances(
				df,
				metric="euclidean",
			)
			print("distance matrix built")

			dist_bar = compute_mean_pairwise_distance(
				dist_matrix
			)
			print("mean distance calculated")

			for combo in binary_combinations:
				# Map the current binary tuple back to its descriptive dictionary keys
				current_weights = {weight_keys[i]: combo[i] for i in range(len(weight_keys))}
				
				# Execute the evaluation on the already-computed PCA data frame
				dist_metrics = run_distance_evaluation(
					dist_matrix=dist_matrix,
					dist_bar=dist_bar,
					ordered_samples=ordered_samples,
					labels_dict=labels_map,
					sample_study_map=SAMPLE_STUDY_MAP,
					experiment_name=file,
					axis_weights=current_weights
				)
				
				try:
					current_corr = dist_metrics.get("SimilarityDistanceSpearman").values[0]
				except Exception:
					current_corr = 0 
				
				# Store the correlation, the weights, and the full metric payload together
				stage_results.append({
					"correlation": current_corr,
					"weights": current_weights,
					"dist_metrics": dist_metrics
				})

			# Sort all combinations by correlation in descending order (highest first)
			stage_results.sort(key=lambda x: x["correlation"], reverse=False)

			# Print out the Top 10 directly to your stdout logs
			print(f"\n📊 Top 10 weight combinations for stage '{file}':")
			for rank, res in enumerate(stage_results[:10], 1):
				# Extra formatting trick: show which keys are turned "ON" (equal to 1)
				active_axes = [k for k, v in res["weights"].items() if v == 1]
				print(f"  Rank {rank:2d} | Corr: {res['correlation']:.4f} | Active Axes: {active_axes}")
			
			# Extract the absolute best run to pass seamlessly into your downstream plotting assets
			all_dist_metrics[file] = stage_results[0]["dist_metrics"]
			best_weights_per_stage[file] = stage_results[0]["weights"]
			
			# Save the clean, sorted list of records (minus the heavy metric objects) to export to text
			if 'all_ranked_runs' not in locals() and 'all_ranked_runs' not in globals():
				all_ranked_runs = {}
			all_ranked_runs[file] = [
				{"rank": i+1, "correlation": round(res["correlation"], 4), "weights": res["weights"]} 
				for i, res in enumerate(stage_results)
			]

		else:
			print(f"Error: Data file not found at {data_path}")

	# --- [GENERATING COMPARISON PLOTS & SAVING COMPLETE RANKINGS] ---
	comparison_output_dir = f"{CLUSTER_EXPLORATION_FIGURES_DIR}/Comparisons"
	os.makedirs(comparison_output_dir, exist_ok=True)
	
	try:
		# Instead of just dumping a single dict, write out a beautifully formatted ranking report
		with open(os.path.join(comparison_output_dir, 'geekyfile.txt'), 'wt') as geeky_file:
			import json
			geeky_file.write("=====================================================================\n")
			geeky_file.write("	  COMPLETE SIMILARITY WEIGHTS COMBINATION RANKING REPORT		\n")
			geeky_file.write("=====================================================================\n\n")
			# Uses json.dumps for clean human-readable indentation of the 128 runs per stage
			geeky_file.write(json.dumps(all_ranked_runs, indent=4))
		print(f"💾 Full ranked report successfully written to {comparison_output_dir}/geekyfile.txt")
	except Exception as e:
		print(f"Unable to write to file because {e}")
	for el in all_dist_metrics:
		plot_similarity_distance_scatter(
			all_dist_metrics[el]["PairwiseSimilarityDistanceDF"].iloc[0],
			output_folder=comparison_output_dir,
			experiment_name=f"dist-sim-plot_{el}"
		)
	
	plot_distance_metrics(
		all_dist_metrics=all_dist_metrics,
		output_folder=comparison_output_dir,
		experiment_name="Distance_Metrics_Comparison_Optimal_Grid_Search",
		plot_ratio=True
	)
	
	if FULL:
		print("\nGenerating Metric Comparisons (gene expression space)...")
		combined_meta = pd.concat(all_metas.values())
		try:
			combined_meta = combined_meta.drop_duplicates()
		except TypeError:
			combined_meta = combined_meta.drop_duplicates(subset=[c for c in combined_meta.columns if combined_meta[c].apply(lambda x: isinstance(x, str)).all()])
		plot_metrics_comparison(metrics_dict=all_metrics, metadata_df=combined_meta, bio_targets=LABEL_AXES, output_folder=comparison_output_dir)

		if not LIGHT_WEIGHT:
			print("\nGenerating Metric Comparisons (BulkFormer latent space)...")
			plot_metrics_comparison(metrics_dict=all_bulk_metrics, metadata_df=combined_meta, bio_targets=LABEL_AXES, output_folder=comparison_output_dir, experiment_name="Bulk_Latent_Comparison")
			print("Generating linked multi-stage UMAP comparison...")
			plot_combined_interactive_projections(embeddings_dict=all_umaps, meta_dicts=all_metas, title="UMAP Cross-Stage Comparison", output_path=f"{comparison_output_dir}/Combined_UMAP.html")
			print("Generating linked multi-stage t-SNE comparison...")
			plot_combined_interactive_projections(embeddings_dict=all_tsnes, meta_dicts=all_metas, title="t-SNE Cross-Stage Comparison", output_path=f"{comparison_output_dir}/Combined_TSNE.html")
			print("Generating Bulk comparison...")
			plot_combined_interactive_projections(embeddings_dict=all_bulk, meta_dicts=all_metas, title="Bulk Cross-Stage Comparison", output_path=f"{comparison_output_dir}/Combined_bulk.html")