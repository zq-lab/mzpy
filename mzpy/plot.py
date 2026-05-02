"""
mzpy.plot
=========

Unified plotting utilities for mass-spectrometry and metabolomics workflows.

All *plotnine*-based functions share a single colour scheme and theme.
*Venn* diagrams and a few specialty plots (centroid spectrum, molecule grid,
chemical network) keep their original backends for technical reasons.

Colour scheme
-------------
Three families are predefined and can be changed globally via
:func:`set_color_scheme`:

* ``sequential`` – blue gradient for continuous data.
* ``diverging``  – centred red-white-blue palette (ColorBrewer2).
* ``qualitative`` – category-safe 10-colour set.

The diverging palette defaults to
``['#ca0020','#f4a582','#f7f7f7','#92c5de','#0571b0']``.
"""

from __future__ import annotations

import warnings
from itertools import chain
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

# ── matplotlib / seaborn ────────────────────────────────────────
import matplotlib
matplotlib.use("Agg")  # 非交互式后端，适用于服务器环境
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, to_hex
import seaborn as sns

# ── plotnine ────────────────────────────────────────────────────
from plotnine import (
    aes,
    coord_flip,
    element_blank,
    element_rect,
    element_text,
    facet_wrap,
    geom_boxplot,
    geom_col,
    geom_hline,
    geom_line,
    geom_point,
    geom_ribbon,
    geom_segment,
    geom_text,
    geom_vline,
    ggplot,
    labs,
    scale_color_brewer,
    scale_color_manual,
    scale_fill_brewer,
    scale_fill_manual,
    scale_shape_manual,
    stat_ellipse,
    theme,
    theme_bw,
    theme_classic,
    theme_matplotlib,
    xlim,
    ylim,
)

# ── sklearn / rdkit / igraph (optional) ─────────────────────────
try:
    from sklearn.decomposition import PCA
    from sklearn.linear_model import LogisticRegression
    _HAS_SKLEARN = True
except Exception:  # pragma: no cover
    _HAS_SKLEARN = False

try:
    from rdkit import Chem, DataStructs  # type: ignore[import-untyped]
    from rdkit.Chem import Draw, rdFingerprintGenerator  # type: ignore[import-untyped]
    _HAS_RDKIT = True
except Exception:  # pragma: no cover
    _HAS_RDKIT = False

try:
    from igraph import Graph, plot as igplot  # type: ignore[import-untyped]
    _HAS_IGRAPH = True
except Exception:  # pragma: no cover
    _HAS_IGRAPH = False

# ═══════════════════════════════════════════════════════════════
#  Global defaults
# ═══════════════════════════════════════════════════════════════

_DEFAULT_COLORS: Dict[str, List[str]] = {
    "sequential": [
        "#f7fbff", "#deebf7", "#c6dbef", "#9ecae1",
        "#6baed6", "#4292c6", "#2171b5", "#08519c", "#08306b",
    ],
    "diverging": [
        "#ca0020", "#f4a582", "#f7f7f7", "#92c5de", "#0571b0",
    ],
    "qualitative": [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
        "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    ],
}

_CURRENT_COLORS = {k: v.copy() for k, v in _DEFAULT_COLORS.items()}

_DEFAULT_FIGSIZE: Tuple[float, float] = (5.0, 5.0)
_DEFAULT_DPI: int = 300
_DEFAULT_FONTSIZE: int = 12


# ═══════════════════════════════════════════════════════════════
#  Public helpers – colour & size
# ═══════════════════════════════════════════════════════════════

def set_color_scheme(scheme: Dict[str, Sequence[str]]) -> None:
    """Replace the global colour palette.

    Parameters
    ----------
    scheme : dict
        Must contain keys ``'sequential'``, ``'diverging'``,
        ``'qualitative'``.  Values are lists of hex colour strings.
    """
    global _CURRENT_COLORS
    for key in ("sequential", "diverging", "qualitative"):
        if key in scheme:
            _CURRENT_COLORS[key] = list(scheme[key])


def set_figure_size(
    size: Tuple[float, float] = (5.0, 5.0),
    dpi: int = 300,
    fontsize: int = 12,
) -> None:
    """Set the default figure size, resolution and font size."""
    global _DEFAULT_FIGSIZE, _DEFAULT_DPI, _DEFAULT_FONTSIZE
    _DEFAULT_FIGSIZE = size
    _DEFAULT_DPI = dpi
    _DEFAULT_FONTSIZE = fontsize


def get_color_scheme() -> Dict[str, List[str]]:
    """Return the currently active colour scheme (deep copy)."""
    return {k: v.copy() for k, v in _CURRENT_COLORS.items()}


# ═══════════════════════════════════════════════════════════════
#  Internal helpers
# ═══════════════════════════════════════════════════════════════

def _get_theme(
    figure_size: Optional[Tuple[float, float]] = None,
    dpi: Optional[int] = None,
    fontsize: Optional[int] = None,
):
    """Build a reusable plotnine theme."""
    fs = fontsize or _DEFAULT_FONTSIZE
    fig = figure_size or _DEFAULT_FIGSIZE
    dp = dpi or _DEFAULT_DPI
    return (
        theme_matplotlib()
        + theme(
            axis_title=element_text(size=fs * 1.3),
            axis_text=element_text(size=fs),
            legend_title=element_text(size=fs),
            legend_text=element_text(size=fs * 0.85),
            legend_background=element_rect(fill=None, color=None),
            legend_position=(0.95, 0.95),
            figure_size=fig,
            dpi=dp,
        )
    )


def _palette(color_type: str = "qualitative") -> List[str]:
    """Return the current palette for a given colour family."""
    return _CURRENT_COLORS.get(color_type, _CURRENT_COLORS["qualitative"])


def _save_plotnine(p, save_to: Optional[str], dpi: Optional[int] = None) -> None:
    if save_to is not None:
        p.save(save_to, dpi=dpi or _DEFAULT_DPI)


def _save_mpl(fig, save_to: Optional[str], dpi: Optional[int] = None) -> None:
    if save_to is not None:
        fig.savefig(save_to, dpi=dpi or _DEFAULT_DPI, bbox_inches="tight")


# ═══════════════════════════════════════════════════════════════
#  1. Bar chart
# ═══════════════════════════════════════════════════════════════

def bar(
    df: pd.DataFrame,
    x: str,
    y: Optional[str] = None,
    fill: Optional[str] = None,
    title: str = "",
    na_replace: str = "other",
    color_type: str = "qualitative",
    figure_size: Optional[Tuple[float, float]] = None,
    save_to: Optional[str] = None,
):
    """Frequency or value bar chart (plotnine).

    If *y* is omitted the frequency of *x* is computed automatically.
    """
    if x not in df.columns:
        raise ValueError(f"Column '{x}' not found in DataFrame.")

    data = df.copy()
    data[x] = data[x].fillna(na_replace).astype(str)

    if y is None:
        freq = data[x].value_counts().reset_index()
        freq.columns = [x, "count"]
        freq = freq.sort_values("count", ascending=True).reset_index(drop=True)
        freq[x] = pd.Categorical(freq[x], categories=freq[x].tolist(), ordered=True)
        y_col = "count"
        y_label = "Frequency"
    else:
        if y not in data.columns:
            raise ValueError(f"Column '{y}' not found in DataFrame.")
        freq = data.sort_values(y, ascending=True).reset_index(drop=True)
        freq[x] = pd.Categorical(freq[x], categories=freq[x].tolist(), ordered=True)
        y_col = y
        y_label = y

    colors = _palette(color_type)

    if fill is not None and fill in freq.columns:
        p = (
            ggplot(freq, aes(x=x, y=y_col, fill=fill))
            + geom_col(show_legend=True)
            + scale_fill_manual(values=colors)
            + coord_flip()
        )
    else:
        p = (
            ggplot(freq, aes(x=x, y=y_col))
            + geom_col(fill=colors[0], show_legend=False)
            + coord_flip()
        )

    p = (
        p
        + labs(title=title, x=x, y=y_label)
        + _get_theme(figure_size=figure_size)
        + theme(
            axis_text_x=element_text(size=_DEFAULT_FONTSIZE * 0.7),
            axis_text_y=element_text(size=_DEFAULT_FONTSIZE * 0.7),
        )
    )

    _save_plotnine(p, save_to)
    return p


# ═══════════════════════════════════════════════════════════════
#  2. Box plot
# ═══════════════════════════════════════════════════════════════

def box(
    df: pd.DataFrame,
    x: str,
    y: str,
    fill: Optional[str] = None,
    facet: Optional[str] = None,
    ncol: int = 4,
    group_order: Optional[List[str]] = None,
    color_type: str = "qualitative",
    show_trend: bool = True,
    log_transform: bool = False,
    figure_size: Optional[Tuple[float, float]] = None,
    save_to: Optional[str] = None,
):
    """Faceted boxplot from long-format data (plotnine)."""
    df_plot = df.copy()
    df_plot[x] = df_plot[x].astype(str)

    if log_transform:
        df_plot = df_plot[df_plot[y] > 0].copy()
        df_plot[y] = np.log10(df_plot[y])
        y_label = f"log10({y})"
    else:
        y_label = y

    if group_order is not None:
        df_plot[x] = pd.Categorical(df_plot[x], categories=group_order, ordered=True)

    colors = _palette(color_type)
    fill_col = fill if fill is not None else x

    p = (
        ggplot(df_plot, aes(x=x, y=y))
        + geom_boxplot(aes(fill=fill_col), outlier_shape=None)
        + geom_point(alpha=0.4, size=1.0)
        + scale_fill_manual(values=colors)
        + labs(x="Group", y=y_label, fill=fill_col)
        + theme_bw()
        + theme(axis_text_x=element_text(rotation=45, hjust=1))
        + _get_theme(figure_size=figure_size)
    )

    if facet is not None:
        p = p + facet_wrap(f"~{facet}", scales="free_y", ncol=ncol)
        if show_trend:
            summary_df = (
                df_plot.groupby([facet, x], observed=True)[y]
                .median()
                .reset_index()
                .rename(columns={y: "summary_value"})
            )
            p = p + geom_line(
                data=summary_df,
                mapping=aes(x=x, y="summary_value", group=1),
                color="grey",
                alpha=0.6,
                size=0.7,
            )

    _save_plotnine(p, save_to)
    return p


# ═══════════════════════════════════════════════════════════════
#  3. Donut chart  (matplotlib)
# ═══════════════════════════════════════════════════════════════

def donut(
    df: pd.DataFrame,
    column: str,
    title: str = "",
    na_replace: str = "other",
    color_type: str = "qualitative",
    show_percent: bool = True,
    figure_size: Optional[Tuple[float, float]] = None,
    save_to: Optional[str] = None,
):
    """Donut / ring chart (matplotlib)."""
    if column not in df.columns:
        raise ValueError(f"Column '{column}' not found in DataFrame")

    data = df.copy()
    data[column] = data[column].fillna(na_replace).astype(str)

    freq = data[column].value_counts().reset_index()
    freq.columns = [column, "count"]
    freq["percent"] = freq["count"] / freq["count"].sum() * 100
    freq["label"] = freq.apply(
        lambda r: f"{r[column]} ({r['percent']:.1f}%)" if show_percent else r[column],
        axis=1,
    )

    colors = _palette(color_type)
    colors = colors * (len(freq) // len(colors) + 1)

    fig, ax = plt.subplots(
        figsize=figure_size or _DEFAULT_FIGSIZE,
        dpi=_DEFAULT_DPI,
        constrained_layout=True,
    )
    wedges, texts = ax.pie(
        freq["count"],
        labels=freq["label"],
        startangle=90,
        counterclock=False,
        colors=colors[: len(freq)],
        textprops={"fontsize": _DEFAULT_FONTSIZE * 0.7},
        wedgeprops=dict(width=0.3, edgecolor="white"),
    )
    centre_circle = plt.Circle((0, 0), 0.70, fc="white")
    ax.add_artist(centre_circle)
    ax.set_title(title, fontsize=_DEFAULT_FONTSIZE * 1.2)
    ax.axis("equal")

    _save_mpl(fig, save_to)
    return fig


# ═══════════════════════════════════════════════════════════════
#  4. Heatmap  (seaborn)
# ═══════════════════════════════════════════════════════════════

def heatmap(
    df: pd.DataFrame,
    data_transfer: Optional[str] = "log10",
    color_type: str = "diverging",
    title: str = r"Heatmap of Log$_{10}$ (Peak Area)",
    xlab: str = "Sample",
    ylab: str = "Metabolite",
    figure_size: Optional[Tuple[float, float]] = None,
    save_to: Optional[str] = None,
):
    """Clustered heatmap (seaborn)."""
    if data_transfer == "log10":
        data = df.apply(lambda x: np.log(x + 1))
        cbar_label = r"log$_{10}$ (peak area)"
    elif data_transfer == "relative":
        base = df.max().max()
        data = (df / base) * 100 if base > 0 else df
        cbar_label = "Relative intensity (%)"
    elif data_transfer is None or data_transfer == "":
        data = df
        cbar_label = "Intensity"
    else:
        raise ValueError(
            f"Unknown data_transfer: {data_transfer}. "
            "Use 'log10', 'relative', or None."
        )

    n_rows, n_cols = data.shape
    figsize = figure_size or (max(5, n_cols * 0.5), max(5, n_rows * 0.5))

    colors = _palette(color_type)
    if len(colors) >= 5:
        cmap = sns.blend_palette(colors, as_cmap=True)
    else:
        cmap = sns.color_palette(colors, as_cmap=True)

    sns.set_style("white")
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        data,
        cmap=cmap,
        center=np.median(data.values.flatten()),
        linewidths=0.05,
        linecolor="white",
        square=True,
        cbar_kws={"shrink": 0.8, "label": cbar_label},
        ax=ax,
    )
    ax.set_title(title, fontsize=16)
    ax.set_xlabel(xlab, fontsize=12)
    ax.set_ylabel(ylab, fontsize=12)
    plt.setp(ax.get_xticklabels(), rotation=90, ha="center")
    plt.setp(ax.get_yticklabels(), rotation=0, ha="right")

    _save_mpl(fig, save_to)
    return fig


# ═══════════════════════════════════════════════════════════════
#  5. Line with error band  (plotnine)
# ═══════════════════════════════════════════════════════════════

def line(
    df: pd.DataFrame,
    x: str,
    y: str,
    group: Optional[str] = None,
    color: Optional[str] = None,
    color_type: str = "qualitative",
    show_error: bool = True,
    figure_size: Optional[Tuple[float, float]] = None,
    save_to: Optional[str] = None,
):
    """Line plot with optional error ribbon (plotnine).

    *df* is expected to be a **summary** table with columns *x*, *y*,
    ``Mean``, ``StdDev`` (or ``Error``).  If you have raw long-format data,
    pre-aggregate it with ``df.groupby(...).agg(Mean=(..., 'mean'),
    StdDev=(..., 'std')).reset_index()``.
    """
    colors = _palette(color_type)
    line_color = colors[0]
    ribbon_color = colors[1] if len(colors) > 1 else colors[0]

    mapping = aes(x=x, y=y)
    if group is not None:
        mapping = aes(x=x, y=y, group=group)
    if color is not None:
        mapping = aes(x=x, y=y, group=group or color, color=color)

    p = (
        ggplot(df, mapping)
        + geom_line(color=line_color if color is None else None)
    )

    if color is not None:
        p = p + scale_color_manual(values=colors)

    if show_error and "Error" in df.columns:
        p = p + geom_ribbon(
            aes(ymin=f"{y} - Error", ymax=f"{y} + Error"),
            alpha=0.15,
            fill=ribbon_color,
        )

    p = (
        p
        + labs(title=f"{y} by {x}", x=x, y=y)
        + theme_bw()
        + theme(
            legend_position="none",
            panel_grid_major=element_blank(),
            panel_grid_minor=element_blank(),
        )
        + _get_theme(figure_size=figure_size)
    )

    _save_plotnine(p, save_to)
    return p


# ═══════════════════════════════════════════════════════════════
#  6. Lollipop chart  (plotnine)
# ═══════════════════════════════════════════════════════════════

def lollipop(
    df: pd.DataFrame,
    x: str,
    y: str,
    fill: Optional[str] = None,
    color_type: str = "qualitative",
    figure_size: Optional[Tuple[float, float]] = None,
    save_to: Optional[str] = None,
):
    """Horizontal lollipop chart (plotnine)."""
    data = df.sort_values(by=x, ascending=True).copy()
    data[y] = pd.Categorical(data[y], categories=data[y].unique(), ordered=True)

    colors = _palette(color_type)

    if fill is not None and fill in data.columns:
        p = ggplot(data, aes(x=x, y=y, fill=fill))
        p = p + scale_fill_manual(values=colors)
    else:
        p = ggplot(data, aes(x=x, y=y))

    p = (
        p
        + geom_segment(aes(x=0, xend=x, y=y, yend=y))
        + geom_point(shape="o", size=3, color="black")
        + _get_theme(figure_size=figure_size)
    )

    _save_plotnine(p, save_to)
    return p


# ═══════════════════════════════════════════════════════════════
#  7. PCA score plot  (plotnine)
# ═══════════════════════════════════════════════════════════════

def pca(
    df: pd.DataFrame,
    groups: Optional[Sequence] = None,
    labels: Optional[Sequence[str]] = None,
    color_type: str = "qualitative",
    palette: Optional[str] = None,  # backward-compat alias
    add_ellipse: bool = True,
    data_transfer: Optional[str] = "log10",
    figure_size: Optional[Tuple[float, float]] = None,
    save_to: Optional[str] = None,
):
    if palette is not None:
        color_type = palette
    """PCA score scatter plot (plotnine).

    *df* should be a sample-by-feature matrix (rows = samples).
    """
    if not _HAS_SKLEARN:
        raise ImportError("scikit-learn is required for PCA plots.")

    if data_transfer == "log10":
        df_plot = np.log10(df + 1)
    elif data_transfer == "relative":
        base = df.max().max()
        df_plot = (df / base) * 100 if base > 0 else df.copy()
    elif data_transfer is None or data_transfer == "":
        df_plot = df.copy()
    else:
        raise ValueError(
            f"Unknown data_transfer: {data_transfer}. "
            "Use 'log10', 'relative', or None."
        )

    pca_model = PCA(n_components=2).fit(df_plot)
    scores = pca_model.transform(df_plot)
    pca_df = pd.DataFrame(scores, columns=["PC1", "PC2"])

    if groups is not None:
        pca_df["group"] = pd.Categorical(groups)
    else:
        pca_df["group"] = pd.Categorical(["Group"] * len(pca_df))

    colors = _palette(color_type)
    p = (
        ggplot(pca_df, aes("PC1", "PC2", fill="group"))
        + labs(
            x=f"PC1: {100 * pca_model.explained_variance_ratio_[0]:.1f} %",
            y=f"PC2: {100 * pca_model.explained_variance_ratio_[1]:.1f} %",
        )
        + scale_fill_manual(values=colors)
        + _get_theme(figure_size=figure_size)
    )

    if add_ellipse:
        p = p + stat_ellipse(geom="polygon", level=0.95, alpha=0.2)

    p = p + geom_point(alpha=0.6, size=3, shape="o", stroke=0)

    if labels is not None:
        pca_df["label"] = labels
        p = p + geom_text(
            label=pca_df.label, nudge_x=0.1, nudge_y=0.1,
            size=_DEFAULT_FONTSIZE * 0.6,
        )

    _save_plotnine(p, save_to)
    return p


# ═══════════════════════════════════════════════════════════════
#  8. PLS-DA score plot  (matplotlib)
# ═══════════════════════════════════════════════════════════════

def plsda(
    T_scores: np.ndarray,
    y: Sequence,
    color_type: str = "qualitative",
    palette: Optional[str] = None,  # backward-compat alias
    data_transfer: Optional[str] = "log10",
    figure_size: Tuple[float, float] = (4, 4),
    save_to: Optional[str] = None,
):
    if palette is not None:
        color_type = palette
    """PLS-DA score scatter with decision regions (matplotlib)."""
    if not _HAS_SKLEARN:
        raise ImportError("scikit-learn is required for PLS-DA plots.")

    if data_transfer == "log10":
        # signed log10 to handle both positive and negative scores safely
        T_plot = np.sign(T_scores) * np.log10(np.abs(T_scores) + 1)
    elif data_transfer == "relative":
        base = np.max(np.abs(T_scores))
        T_plot = (T_scores / base) * 100 if base > 0 else T_scores.copy()
    elif data_transfer is None or data_transfer == "":
        T_plot = T_scores.copy()
    else:
        raise ValueError(
            f"Unknown data_transfer: {data_transfer}. "
            "Use 'log10', 'relative', or None."
        )

    class_names = np.unique(y)
    df_scores = pd.DataFrame(T_plot, columns=["LV1", "LV2"])
    df_scores["group"] = y

    clf = LogisticRegression(max_iter=1000)
    clf.fit(df_scores[["LV1", "LV2"]].to_numpy(), y)

    pad = 0.10
    x_min, x_max = df_scores["LV1"].min(), df_scores["LV1"].max()
    y_min, y_max = df_scores["LV2"].min(), df_scores["LV2"].max()
    x_pad = (x_max - x_min) * pad
    y_pad = (y_max - y_min) * pad
    x_min, x_max = x_min - x_pad, x_max + x_pad
    y_min, y_max = y_min - y_pad, y_max + y_pad

    xx, yy = np.meshgrid(
        np.linspace(x_min, x_max, 400),
        np.linspace(y_min, y_max, 400),
    )
    grid = np.c_[xx.ravel(), yy.ravel()]
    Z = clf.predict(grid).reshape(xx.shape)

    fig, ax = plt.subplots(figsize=figure_size, dpi=_DEFAULT_DPI)
    colors = _palette(color_type)
    palette_map = {c: colors[i % len(colors)] for i, c in enumerate(class_names)}
    idx_map = {c: i for i, c in enumerate(class_names)}
    Z_idx = np.vectorize(idx_map.get)(Z)
    cmap_light = ListedColormap(
        [sns.desaturate(palette_map[c], 0.9) for c in class_names]
    )

    ax.contourf(xx, yy, Z_idx, alpha=0.25, cmap=cmap_light, levels=len(class_names))
    ax.contour(
        xx, yy, Z_idx, levels=np.arange(len(class_names)),
        colors="k", alpha=0.2, linewidths=0.8,
    )

    sns.scatterplot(
        data=df_scores, x="LV1", y="LV2", hue="group",
        s=30, edgecolor="none", linewidth=1.0,
        palette=[palette_map[c] for c in class_names],
        ax=ax,
    )

    ax.set_xlabel("LV1 (PLS-DA)")
    ax.set_ylabel("LV2 (PLS-DA)")
    ax.set_title("PLS-DA scores with decision regions")
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_linewidth(1.2)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.legend(title="Group", frameon=False)
    plt.tight_layout()

    _save_mpl(fig, save_to)
    return fig


# ═══════════════════════════════════════════════════════════════
#  9. Volcano plot  (plotnine)
# ═══════════════════════════════════════════════════════════════

def volcano(
    df: pd.DataFrame,
    x: str,
    y: str,
    fill: str,
    xcut: float = 1.0,
    ycut: float = 2.0,
    title: str = "",
    color_type: str = "diverging",
    xlab: Optional[str] = None,
    ylab: Optional[str] = None,
    palette: Optional[str] = None,  # backward-compat alias
    figure_size: Optional[Tuple[float, float]] = None,
    save_to: Optional[str] = None,
):
    if palette is not None:
        color_type = palette
    """Volcano plot (plotnine)."""
    colors = _palette(color_type)
    # Pick two extremes + neutral grey
    if len(colors) >= 2:
        up_color = colors[0]
        dn_color = colors[-1]
    else:
        up_color = dn_color = colors[0]
    no_color = "#D3D3D3"

    # Guard against all-NaN columns (e.g. zero-variance data)
    x_vals = df[x].replace([np.inf, -np.inf], np.nan).dropna()
    x_limit = max(abs(x_vals.min()), abs(x_vals.max())) if not x_vals.empty else 1.0

    y_vals = df[y].replace([np.inf, -np.inf], np.nan).dropna()
    y_max = max(2.5, y_vals.max()) if not y_vals.empty else 2.5

    p = (
        ggplot(df, aes(x=x, y=y, color=fill))
        + geom_point(alpha=0.5, size=2, shape="o", stroke=0)
        + labs(
            title=f"{title}, Volcano Plot" if title else "Volcano Plot",
            x=xlab if xlab is not None else x,
            y=ylab if ylab is not None else y,
        )
        + xlim(-x_limit, x_limit)
        + ylim(0, y_max)
        + scale_color_manual(
            values={"up": up_color, "dn": dn_color, "no": no_color}
        )
        + geom_vline(xintercept=-xcut, linetype="dashed", color="grey")
        + geom_vline(xintercept=xcut, linetype="dashed", color="grey")
        + geom_hline(yintercept=ycut, linetype="dashed", color="grey")
        + _get_theme(figure_size=figure_size)
    )

    _save_plotnine(p, save_to)
    return p


# ═══════════════════════════════════════════════════════════════
#  10. Scatter / TIC-RT-mz style  (plotnine)
# ═══════════════════════════════════════════════════════════════

def scatter(
    df: pd.DataFrame,
    x: str,
    y: str,
    size: Optional[str] = None,
    color: Optional[str] = None,
    shape: Optional[str] = None,
    alpha: float = 0.55,
    color_type: str = "qualitative",
    figure_size: Optional[Tuple[float, float]] = None,
    save_to: Optional[str] = None,
):
    """Generic scatter plot (plotnine)."""
    colors = _palette(color_type)

    aes_map = {"x": x, "y": y}
    if size is not None:
        aes_map["size"] = size
    if color is not None:
        aes_map["color"] = color
    if shape is not None:
        aes_map["shape"] = shape

    p = ggplot(df, aes(**aes_map)) + geom_point(alpha=alpha)

    if color is not None:
        p = p + scale_color_manual(values=colors)
    else:
        p = p + geom_point(color=colors[0], alpha=alpha)

    if shape is not None:
        df[shape] = df[shape].astype(str)
        # Provide a default manual shape mapping if the user didn't supply one
        # (plotnine will auto-pick shapes otherwise)

    p = (
        p
        + labs(x=x, y=y)
        + theme_classic()
        + _get_theme(figure_size=figure_size)
    )

    _save_plotnine(p, save_to)
    return p


# ═══════════════════════════════════════════════════════════════
#  11. Venn diagram  (matplotlib)
# ═══════════════════════════════════════════════════════════════

def venn(
    data: Dict[str, Sequence],
    fill_mode: Optional[List[str]] = None,
    color_type: str = "qualitative",
    palette: Optional[str] = None,  # backward-compat alias
    alpha: float = 0.65,
    fontsize: Optional[int] = None,
    figure_size: Optional[Tuple[float, float]] = None,
    save_to: Optional[str] = None,
):
    if palette is not None:
        color_type = palette
    """Draw a 2-, 3- or 4-set Venn diagram (matplotlib).

    The number of circles / ellipses is inferred automatically from
    ``len(data)`` (must be 2, 3 or 4).

    Parameters
    ----------
    data : dict
        Mapping from set names to element sequences.
    fill_mode : list of str, optional
        Label modes, e.g. ``["number"]``, ``["logic", "percent"]``.
    """
    if fill_mode is None:
        fill_mode = ["number"]
    if not isinstance(data, dict):
        raise ValueError("Input data must be a dictionary.")
    n_sets = len(data)
    if n_sets not in (2, 3, 4):
        raise ValueError("Venn diagram supports 2, 3 or 4 sets only.")

    colors = _palette(color_type)
    cmap = ListedColormap(colors)
    fs = fontsize or _DEFAULT_FONTSIZE

    fig, ax = plt.subplots(figsize=figure_size or (9, 7), dpi=_DEFAULT_DPI)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_ylim(0.0, 1.0)
    ax.set_xlim(0.0, 1.0)

    values = list(data.values())
    keys = list(data.keys())
    labels = _venn_divide(values, fill=fill_mode)

    draw_func = globals()[f"_venn_draw{n_sets}"]
    draw_func(ax, labels, keys, cmap, alpha, fs)

    if save_to:
        fig.savefig(save_to, transparent=True)
    return fig


def _venn_divide(data, fill=None):
    """Generate region labels for a Venn diagram."""
    if fill is None:
        fill = ["number"]
    N = len(data)
    sets_data = [set(data[i]) for i in range(N)]
    s_all = set(chain(*data))
    set_collections = {}

    for n in range(1, 2 ** N):
        key = bin(n).split("0b")[-1].zfill(N)
        value = s_all
        inter = [sets_data[i] for i in range(N) if key[i] == "1"]
        diff = [sets_data[i] for i in range(N) if key[i] == "0"]
        for s in inter:
            value = value & s
        for s in diff:
            value = value - s
        set_collections[key] = value

    labels = {k: "" for k in set_collections}
    if "logic" in fill:
        for k in set_collections:
            labels[k] = k + ": "
    if "number" in fill:
        for k in set_collections:
            labels[k] += str(len(set_collections[k]))
    if "percent" in fill:
        total = len(s_all)
        for k in set_collections:
            labels[k] += "(%.1f%%)" % (100.0 * len(set_collections[k]) / total)
    return labels


def _venn_draw_ellipse(ax, xy, width, height, angle, color, alpha):
    ax.add_patch(
        patches.Ellipse(
            xy=xy, width=width, height=height, angle=angle,
            alpha=alpha, color=color,
        )
    )


def _venn_draw_text(ax, x, y, text, fs, ha="center", va="center"):
    ax.text(x, y, text, fontsize=fs, ha=ha, va=va)


def _venn_draw2(ax, labels, names, cmap, alpha, fs):
    ax.set_ylim(0.0, 0.7)
    _venn_draw_ellipse(ax, (0.375, 0.3), 0.5, 0.5, 0, cmap(0), alpha)
    _venn_draw_ellipse(ax, (0.625, 0.3), 0.5, 0.5, 0, cmap(1), alpha)
    _venn_draw_text(ax, 0.74, 0.30, labels.get("01", ""), fs)
    _venn_draw_text(ax, 0.26, 0.30, labels.get("10", ""), fs)
    _venn_draw_text(ax, 0.50, 0.30, labels.get("11", ""), fs)
    _venn_draw_text(ax, 0.20, 0.56, names[0], fs)
    _venn_draw_text(ax, 0.80, 0.56, names[1], fs)


def _venn_draw3(ax, labels, names, cmap, alpha, fs):
    _venn_draw_ellipse(ax, (0.333, 0.633), 0.55, 0.55, 0, cmap(0), alpha)
    _venn_draw_ellipse(ax, (0.666, 0.633), 0.55, 0.55, 0, cmap(1), alpha)
    _venn_draw_ellipse(ax, (0.500, 0.310), 0.55, 0.55, 0, cmap(2), alpha)
    _venn_draw_text(ax, 0.50, 0.27, labels.get("001", ""), fs)
    _venn_draw_text(ax, 0.73, 0.65, labels.get("010", ""), fs)
    _venn_draw_text(ax, 0.61, 0.46, labels.get("011", ""), fs)
    _venn_draw_text(ax, 0.27, 0.65, labels.get("100", ""), fs)
    _venn_draw_text(ax, 0.39, 0.46, labels.get("101", ""), fs)
    _venn_draw_text(ax, 0.50, 0.65, labels.get("110", ""), fs)
    _venn_draw_text(ax, 0.50, 0.51, labels.get("111", ""), fs)
    _venn_draw_text(ax, 0.15, 0.87, names[0], fs)
    _venn_draw_text(ax, 0.85, 0.87, names[1], fs)
    _venn_draw_text(ax, 0.50, 0.02, names[2], fs)


def _venn_draw4(ax, labels, names, cmap, alpha, fs):
    o, dx, dy = 0.500, 0.18, 0.08
    _venn_draw_ellipse(ax, (o - dx, o - dy), 4 * dx, 2 * dx, 135, cmap(0), alpha)
    _venn_draw_ellipse(ax, (o, o), 4 * dx, 2 * dx, 135, cmap(1), alpha)
    _venn_draw_ellipse(ax, (o, o), 4 * dx, 2 * dx, 45, cmap(2), alpha)
    _venn_draw_ellipse(ax, (o + dx, o - dy), 4 * dx, 2 * dx, 45, cmap(3), alpha)
    positions = [
        (o + dx * 2.00, o + dy * 0.50, labels.get("0001", "")),
        (o + dx * 0.75, o + dy * 2.50, labels.get("0010", "")),
        (o + dx * 1.25, o + dy * 1.25, labels.get("0011", "")),
        (o - dx * 0.75, o + dy * 2.50, labels.get("0100", "")),
        (o + dx, o - dy * 2.00, labels.get("0101", "")),
        (o, o + dy * 1.25, labels.get("0110", "")),
        (o + dx * 0.75, o - dy * 0.25, labels.get("0111", "")),
        (o - dx * 2.00, o + dy * 0.50, labels.get("1000", "")),
        (o, o - dy * 3.75, labels.get("1001", "")),
        (o - dx, o - dy * 2.00, labels.get("1010", "")),
        (o - dx * 0.25, o - dy * 2.75, labels.get("1011", "")),
        (o - dx * 1.25, o + dy * 1.25, labels.get("1100", "")),
        (o + dx * 0.25, o - dy * 2.75, labels.get("1101", "")),
        (o - dx * 0.75, o - dy * 0.25, labels.get("1110", "")),
        (o, o - dy * 1.75, labels.get("1111", "")),
        (o - dx * 2.25, o + dy * 2.75, names[0]),
        (o - dx * 1.00, o + dy * 3.75, names[1]),
        (o + dx * 1.00, o + dy * 3.75, names[2]),
        (o + dx * 2.25, o + dy * 2.75, names[3]),
    ]
    for x, y, text in positions:
        _venn_draw_text(ax, x, y, text, fs)


# ═══════════════════════════════════════════════════════════════
#  12. Colour swatch  (matplotlib)
# ═══════════════════════════════════════════════════════════════

def swatch(
    colors: Optional[Sequence[str]] = None,
    color_type: str = "qualitative",
    figure_size: Optional[Tuple[float, float]] = None,
    save_to: Optional[str] = None,
):
    """Draw colour swatches (matplotlib)."""
    if colors is None:
        colors = _palette(color_type)
    colors = list(colors)
    block = 2.5 / 2.54  # 2.5 cm in inch
    fig, ax = plt.subplots(figsize=figure_size or (len(colors) * block, block))
    for i, c in enumerate(colors):
        ax.add_patch(plt.Rectangle((i, 0), 1, 1, color=c))
    ax.set_xlim(0, len(colors))
    ax.set_ylim(0, 1)
    ax.set_xticks([i + 0.5 for i in range(len(colors))])
    ax.set_xticklabels(colors, rotation=45, ha="right", fontsize=8)
    ax.set_yticks([])
    ax.set_title("Color Swatches")
    plt.tight_layout()
    _save_mpl(fig, save_to)
    return fig


# ═══════════════════════════════════════════════════════════════
#  13. Chemical network  (igraph + rdkit)
# ═══════════════════════════════════════════════════════════════

def chemicals_net(
    smiles_list: List[str],
    threshold: float = 0.6,
    fp_type: str = "morgan",
    radius: int = 2,
    n_bits: int = 2048,
    use_sanitization: bool = True,
    save_to: Optional[str] = None,
):
    """Build a chemical-similarity network from SMILES (igraph).

    Returns the igraph plot object; also writes to *save_to* when given.
    """
    if not _HAS_RDKIT:
        raise ImportError("rdkit is required for chemical network plots.")
    if not _HAS_IGRAPH:
        raise ImportError("igraph is required for chemical network plots.")
    if not smiles_list:
        raise ValueError("smiles_list cannot be empty")

    mols = []
    fps = []
    if fp_type.lower() == "morgan":
        morgan_gen = rdFingerprintGenerator.GetMorganGenerator(
            radius=radius, fpSize=n_bits
        )
    for s in smiles_list:
        mol = Chem.MolFromSmiles(s, sanitize=use_sanitization)
        if mol is None:
            raise ValueError(f"Cannot parse SMILES: {s}")
        if fp_type.lower() == "morgan":
            fp = morgan_gen.GetFingerprint(mol)
        elif fp_type.lower() in {"maccs", "rdkitmaccs"}:
            from rdkit.Chem import MACCSkeys  # type: ignore[import-untyped]
            fp = MACCSkeys.GenMACCSKeys(mol)
        else:
            raise ValueError("fp_type must be 'morgan' or 'maccs'")
        mols.append(mol)
        fps.append(fp)

    n = len(fps)
    sim_matrix = np.zeros((n, n), dtype=float)
    for i in range(n):
        sim_matrix[i, i] = 1.0
        for j in range(i + 1, n):
            sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            sim_matrix[i, j] = sim_matrix[j, i] = sim

    g = Graph()
    g.add_vertices(n)
    g.vs["label"] = [str(i) for i in range(n)]
    g.vs["smiles"] = smiles_list

    edges, weights = [], []
    for i in range(n):
        for j in range(i + 1, n):
            if sim_matrix[i, j] >= threshold:
                edges.append((i, j))
                weights.append(float(sim_matrix[i, j]))
    if edges:
        g.add_edges(edges)
        g.es["weight"] = weights

    layout = g.layout("fr")
    img = igplot(
        g,
        target=save_to,
        layout=layout,
        vertex_size=12,
        vertex_label=g.vs["label"],
        vertex_label_size=8,
        vertex_color="#d7191c",
        edge_color="grey",
        background=None,
    )
    return img


# ═══════════════════════════════════════════════════════════════
#  14. Molecule grid  (rdkit)
# ═══════════════════════════════════════════════════════════════

def molecules_grid(
    smiles_list: List[str],
    cols: int = 4,
    mols_per_row: Optional[int] = None,
    sub_img_size: Tuple[int, int] = (250, 250),
    legends: bool = False,
    save_to: Optional[str] = None,
):
    """Render a grid of 2-D molecular structures (rdkit).

    Returns a PIL Image; if *save_to* ends with ``.svg`` the SVG string is
    written to disk.
    """
    if not _HAS_RDKIT:
        raise ImportError("rdkit is required for molecule grid plots.")
    if not isinstance(smiles_list, list):
        raise TypeError("smiles_list must be a list")
    if mols_per_row is None:
        mols_per_row = cols
    if mols_per_row <= 0:
        raise ValueError("cols / mols_per_row must be positive")

    mols = []
    legends_list = []
    for i, smi in enumerate(smiles_list):
        if smi is None:
            mol = None
        else:
            mol = Chem.MolFromSmiles(str(smi).strip())
        if mol is None:
            mol = Chem.MolFromSmiles("")
        mols.append(mol)
        legends_list.append(str(i))

    img = Draw.MolsToGridImage(
        mols,
        molsPerRow=mols_per_row,
        subImgSize=sub_img_size,
        legends=legends_list if legends else None,
        useSVG=True,
    )
    if save_to:
        with open(save_to, "w", encoding="utf-8") as f:
            f.write(img.data)  # type: ignore[attr-defined]
    return img


# ═══════════════════════════════════════════════════════════════
#  15. Centroid spectrum  (matplotlib)
# ═══════════════════════════════════════════════════════════════

def spectrum(
    peaks: Sequence[Sequence[float]],
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 4),
    color: Optional[str] = None,
    line_width: float = 1.6,
    normalize: bool = True,
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    show: bool = True,
    save_to: Optional[str] = None,
    show_mz_labels: bool = True,
    mz_label_fmt: str = "{:.4f}",
    mz_label_fontsize: int = 10,
    mz_label_offset_ratio: float = 0.001,
    mz_label_horizontal_offset: float = 0.0,
    mz_label_intensity_threshold: float = 0.1,
):
    """Draw a centroid mass spectrum (matplotlib stem plot)."""
    if color is None:
        color = _palette("qualitative")[0]

    if peaks is None or len(peaks) == 0:
        raise ValueError("peaks cannot be empty")
    arr = np.asarray(peaks, dtype=float)
    if arr.ndim != 2 or arr.shape[1] < 2:
        raise ValueError("peaks must be [[mz, intensity], ...]")

    mz = arr[:, 0]
    intensity = arr[:, 1]
    if normalize:
        max_i = float(np.max(intensity)) if len(intensity) else 0.0
        if max_i > 0:
            intensity = intensity / max_i

    if xlim is None:
        xmin, xmax = float(np.min(mz)), float(np.max(mz))
        dx = xmax - xmin
        pad = dx * 0.02 if dx > 0 else 1.0
        xlim_eff = (xmin - pad, xmax + pad)
    else:
        xlim_eff = xlim
    if ylim is None:
        ymin = float(np.min(intensity)) if len(intensity) else 0.0
        ymax = float(np.max(intensity)) if len(intensity) else 1.0
        pad = (ymax - ymin) * 0.05 if ymax > ymin else 0.1
        ylim_eff = (ymin - pad, ymax + pad)
    else:
        ylim_eff = ylim

    xlim_eff = (min(xlim_eff[0], 0.0), max(xlim_eff[1], 0.0))
    ylim_eff = (min(ylim_eff[0], 0.0), max(ylim_eff[1], 0.0))

    fig, ax = plt.subplots(figsize=figsize)
    _, stemlines, _ = ax.stem(
        mz, intensity, linefmt=color, markerfmt=" ", basefmt=" ", bottom=0.0
    )
    try:
        stemlines.set_linewidth(line_width)
    except Exception:
        pass

    if show_mz_labels and len(intensity) > 0:
        threshold = mz_label_intensity_threshold * float(np.max(intensity))
        y_range = float(np.max(intensity) - np.min(intensity))
        y_offset = y_range * mz_label_offset_ratio if y_range > 0 else 0.01
        for xi, yi in zip(mz, intensity):
            if xi < xlim_eff[0] or xi > xlim_eff[1]:
                continue
            if yi < threshold:
                continue
            ax.text(
                xi + mz_label_horizontal_offset,
                yi + y_offset,
                mz_label_fmt.format(xi),
                ha="center", va="bottom", rotation=0,
                fontsize=mz_label_fontsize, color=color,
            )

    ax.set_xlim(xlim_eff)
    ax.set_ylim(ylim_eff)
    ax.spines["left"].set_position(("data", 0.0))
    ax.spines["bottom"].set_position(("data", 0.0))
    ax.set_xlabel("m/z")
    ax.set_ylabel("Intensity" if not normalize else "Normalized Intensity")
    if title:
        ax.set_title(title)
    plt.tight_layout()

    if save_to is not None:
        plt.savefig(save_to, dpi=300, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close()
    return fig


# ═══════════════════════════════════════════════════════════════
#  Quick smoke-test when executed directly
# ═══════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("mzpy.plot smoke test")
    print("-" * 40)
    print("Active colour scheme:")
    for k, v in get_color_scheme().items():
        print(f"  {k:12s}: {v[:3]} ...")
    print("\nTry:  from mzpy.plot import bar, venn, ...")
