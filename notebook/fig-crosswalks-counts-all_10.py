from __future__ import annotations

from pathlib import Path
import matplotlib as mpl
import numpy as np
import requests
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from functools import lru_cache

import warnings
warnings.filterwarnings("ignore")

def load_crosswalk(url: str, source: str):
    df = pd.read_csv(url, skiprows=10)
    df["source"] = source
    return df

@lru_cache(maxsize=None)
def get_label(uberon_id):
    iri = f'http://purl.obolibrary.org/obo/{uberon_id.replace(":", "_")}'
    url = "https://www.ebi.ac.uk/ols4/api/ontologies/uberon/terms"
    resp = requests.get(url, params={"iri": iri}, timeout=30)
    if resp.status_code == 200 and resp.json()["_embedded"]["terms"]:
        return resp.json()["_embedded"]["terms"][0]["label"]
    return None

LOCAL_CROSSWALKS = {
    "azimuth": "https://cdn.humanatlas.io/digital-objects/ctann/azimuth/v1.4/assets/azimuth-crosswalk.csv",
    "celltypist": "https://cdn.humanatlas.io/digital-objects/ctann/celltypist/v1.3/assets/celltypist-crosswalk.csv",
    "deepcelltypes": "https://cdn.humanatlas.io/digital-objects/ctann/deepcelltypes/v1.2/assets/deepcelltypes-crosswalk.csv",
    "deepcelltypes_hubmap": "https://cdn.humanatlas.io/digital-objects/ctann/deepcelltypes-hubmap/v1.2/assets/deepcelltypes-hubmap-crosswalk.csv",
    "frmatch": "https://cdn.humanatlas.io/digital-objects/ctann/frmatch/v1.0/assets/frmatch-crosswalk.csv",
    "pan_human_azimuth": "https://cdn.humanatlas.io/digital-objects/ctann/pan-human-azimuth/v1.2/assets/pan-human-azimuth-crosswalk.csv",
    "popv": "https://cdn.humanatlas.io/digital-objects/ctann/popv/v1.4/assets/popv-crosswalk.csv",
    "ribca": "https://cdn.humanatlas.io/digital-objects/ctann/ribca/v1.0/assets/ribca-crosswalk.csv",
    "stellar": "https://cdn.humanatlas.io/digital-objects/ctann/stellar/v1.0/assets/stellar-crosswalk.csv",
    "vccf": "https://cdn.humanatlas.io/digital-objects/ctann/vccf/v1.2/assets/vccf-crosswalk.csv",
}

SOURCE_LABELS = {
    "azimuth": "Azimuth",
    "celltypist": "CellTypist",
    "deepcelltypes": "deepcelltypes",
    "deepcelltypes_hubmap": "deepcelltypes-HuBMAP",
    "frmatch": "FRMatch",
    "pan_human_azimuth": "Pan-human Azimuth",
    "popv": "popV",
    "ribca": "RIBCA",
    "stellar": "STELLAR",
    "vccf": "VCCF",
}

# Figure / layout settings
SORT_MODE = "alphabetical"  # or "coverage"
FIGURE_FONT_SIZE = 7.0
ROW_SPACING_SCALE = 0.50
BOTTOM_SPACE_SCALE = 1.65
FIGURE_WIDTH_MM = 235
ORIGINAL_FIGURE_HEIGHT_MM = 135
ORIGINAL_BOTTOM_FRACTION = 0.39
ORIGINAL_TOP_FRACTION = 0.975
MM_TO_INCH = 1 / 25.4

OUTPUT_DIR = Path("vis")
OUTPUT_STEM = "crosswalk-fig-all-10"


def load_crosswalk(path: str, source: str) -> pd.DataFrame:
    df = pd.read_csv(path, skiprows=10)
    df["source"] = source
    return df


def marker_area(count_like) -> np.ndarray:
    return 6 + 0.42 * np.asarray(count_like, dtype=float)


def format_legend_value(value: float) -> str:
    return f"{value:.0f}" if float(value).is_integer() else f"{value:.1f}"


def main() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": FIGURE_FONT_SIZE,
        "axes.labelsize": FIGURE_FONT_SIZE,
        "xtick.labelsize": FIGURE_FONT_SIZE,
        "ytick.labelsize": FIGURE_FONT_SIZE,
        "legend.fontsize": FIGURE_FONT_SIZE,
        "legend.title_fontsize": FIGURE_FONT_SIZE,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })

    print("Loading local crosswalk data...")
    df_combined = pd.concat(
        [load_crosswalk(path, source) for source, path in LOCAL_CROSSWALKS.items()],
        ignore_index=True,
    )

    # Size metric: unique mapped CL labels per tool-organ combination.
    df_size = (
        df_combined.groupby(["source", "Organ_ID"], as_index=False)["CL_Label"]
        .nunique()
        .rename(columns={"CL_Label": "n_unique_cell_label"})
    )

    # Color metric: exact_match_rate = n_exact_matches / n_mappings.
    mapping_keys = [
        "source",
        "Organ_ID",
        "Annotation_Label_ID",
        "CL_ID",
        "CL_Match",
    ]
    df_mappings = df_combined.drop_duplicates(mapping_keys).copy()

    df_match_quality = (
        df_mappings.groupby(["source", "Organ_ID"], as_index=False)
        .agg(
            n_mappings=("CL_Match", "size"),
            n_exact_matches=("CL_Match", lambda s: s.eq("skos:exactMatch").sum()),
            n_narrow_matches=("CL_Match", lambda s: s.eq("skos:narrowMatch").sum()),
            n_broad_matches=("CL_Match", lambda s: s.eq("skos:broadMatch").sum()),
        )
    )
    df_match_quality["exact_match_rate"] = (
        df_match_quality["n_exact_matches"] / df_match_quality["n_mappings"]
    )

    df_plot = df_size.merge(
        df_match_quality,
        on=["source", "Organ_ID"],
        how="left",
        validate="one_to_one",
    )

    df_plot["source"] = df_plot["source"].map(SOURCE_LABELS)
    df_plot["source"] = pd.Categorical(
        df_plot["source"],
        categories=list(SOURCE_LABELS.values()),
        ordered=True,
    )
    # Resolve each identifier once, then reuse its human-readable label for all tools.
    organ_labels = {
        organ_id: get_label(organ_id)
        for organ_id in df_plot["Organ_ID"].dropna().unique()
    }
    df_plot["organ_label"] = df_plot["Organ_ID"].map(organ_labels)
    df_plot["organ_label"] = (
        df_plot["organ_label"]
        .fillna("Unspecified")
        .astype(str)
        .replace({"None": "Unspecified", "nan": "Unspecified"})
    )

    tool_order = [
        label for label in SOURCE_LABELS.values()
        if label in set(df_plot["source"].astype(str))
    ]

    organ_summary = (
        df_plot.groupby("organ_label", observed=True)
        .agg(
            n_tools=("source", "nunique"),
            total_cell_types=("n_unique_cell_label", "sum"),
        )
    )

    if SORT_MODE == "alphabetical":
        organ_order = sorted(organ_summary.index, key=str.casefold)
    elif SORT_MODE == "coverage":
        organ_order = (
            organ_summary
            .sort_values(["n_tools", "total_cell_types"], ascending=[False, False], kind="stable")
            .index
            .tolist()
        )
    else:
        raise ValueError("SORT_MODE must be 'alphabetical' or 'coverage'.")

    x_position = {organ: i for i, organ in enumerate(organ_order)}
    y_position = {tool: i for i, tool in enumerate(tool_order)}

    bottom_margin_mm = ORIGINAL_FIGURE_HEIGHT_MM * ORIGINAL_BOTTOM_FRACTION * BOTTOM_SPACE_SCALE
    top_margin_mm = ORIGINAL_FIGURE_HEIGHT_MM * (1 - ORIGINAL_TOP_FRACTION)
    original_panel_height_mm = ORIGINAL_FIGURE_HEIGHT_MM * (ORIGINAL_TOP_FRACTION - ORIGINAL_BOTTOM_FRACTION)
    panel_height_mm = original_panel_height_mm * ROW_SPACING_SCALE
    figure_height_mm = bottom_margin_mm + panel_height_mm + top_margin_mm

    figure_size = (
        FIGURE_WIDTH_MM * MM_TO_INCH,
        figure_height_mm * MM_TO_INCH,
    )

    print(f"Preparing figure, height={figure_height_mm:.2f} mm")
    fig, ax = plt.subplots(figsize=figure_size)

    norm = mpl.colors.Normalize(vmin=0.0, vmax=1.0)
    cmap = plt.cm.Reds
    scatter = ax.scatter(
        df_plot["organ_label"].map(x_position),
        df_plot["source"].astype(str).map(y_position),
        s=marker_area(df_plot["n_unique_cell_label"]),
        c=df_plot["exact_match_rate"],
        cmap=cmap,
        norm=norm,
        alpha=0.95,
        linewidths=0.15,
        edgecolors="#8a8a8a",
    )

    ax.set_xticks(np.arange(len(organ_order)))
    ax.set_xticklabels(
        organ_order,
        rotation=45,
        ha="right",
        va="top",
        fontsize=FIGURE_FONT_SIZE,
        fontfamily="sans-serif",
        rotation_mode="anchor",
    )

    ax.set_yticks(np.arange(len(tool_order)))
    ax.set_yticklabels(
        tool_order,
        fontsize=FIGURE_FONT_SIZE,
        fontfamily="sans-serif",
    )

    ax.invert_yaxis()
    ax.set_xlabel("Organ", labelpad=7)
    ax.set_ylabel("Tool", labelpad=7)
    ax.set_xlim(-0.6, len(organ_order) - 0.4)
    ax.set_ylim(len(tool_order) - 0.4, -0.6)
    ax.grid(False)

    # Show all four spines so the axes form a box.
    for side in ["top", "right", "bottom", "left"]:
        ax.spines[side].set_visible(True)
        ax.spines[side].set_linewidth(0.6)

    ax.tick_params(axis="x", pad=1, labelsize=FIGURE_FONT_SIZE)
    ax.tick_params(axis="y", pad=2, labelsize=FIGURE_FONT_SIZE)

    # Size legend.
    count_values = df_plot["n_unique_cell_label"].astype(float)
    min_count = float(count_values.min())
    max_count = float(count_values.max())
    midpoint_target = (min_count + max_count) / 2
    observed_counts = np.sort(count_values.dropna().unique())
    midpoint_count = float(observed_counts[np.argmin(np.abs(observed_counts - midpoint_target))])
    legend_values = [min_count, midpoint_count, max_count]
    legend_handles = [
        ax.scatter([], [], s=marker_area(value), color="#d7301f", alpha=0.70, linewidths=0)
        for value in legend_values
    ]
    legend_labels = [format_legend_value(v) for v in legend_values]

    # Leave a dedicated right-hand column for the two legends.
    fig.subplots_adjust(
        left=0.15,
        right=0.74,
        top=1 - top_margin_mm / figure_height_mm,
        bottom=bottom_margin_mm / figure_height_mm,
    )

    # Position the colorbar and size legend in a stacked column to the right of the main axes.
    ax_pos = ax.get_position()
    right_col_x = ax_pos.x1 + 0.035
    cbar_width = min(0.18, 0.98 - right_col_x)
    cbar_height = 0.035
    cbar_y = ax_pos.y0 + 0.64 * ax_pos.height

    cbar_ax = fig.add_axes([right_col_x, cbar_y, cbar_width, cbar_height])
    cbar = fig.colorbar(scatter, cax=cbar_ax, orientation="horizontal")
    ticks = np.linspace(0, 1, 6)
    cbar.set_ticks(ticks)
    cbar.ax.set_xticklabels([f"{int(v * 100)}%" for v in ticks])
    cbar.ax.set_title("Exact match rate", fontsize=FIGURE_FONT_SIZE, pad=6)

    # Place the size legend below the colorbar but keep it above the x-axis level.
    legend_y = ax_pos.y0 + 0.46 * ax_pos.height
    fig.legend(
        legend_handles,
        legend_labels,
        title="#Cell types",
        loc="upper left",
        bbox_to_anchor=(right_col_x, legend_y),
        bbox_transform=fig.transFigure,
        frameon=False,
        borderaxespad=0,
        labelspacing=1.0,
        handleheight=1.0,
        handletextpad=0.9,
    )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    svg_path = OUTPUT_DIR / f"{OUTPUT_STEM}.svg"
    png_path = OUTPUT_DIR / f"{OUTPUT_STEM}.png"

    fig.savefig(png_path, dpi=600, facecolor="white", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(svg_path, facecolor="white", bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)

if __name__ == "__main__":
    main()
