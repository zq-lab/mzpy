"""
mzpy.plotly
===========

Plotly-based plotting utilities for web rendering.
All functions return a ``plotly.graph_objects.Figure`` object that can be
embedded directly in Dash / Streamlit / Jupyter or exported to HTML/JSON.

Colour scheme
-------------
The same three families as ``mzpy.plot`` are predefined:
``sequential``, ``diverging``, ``qualitative``.
Change them globally via :func:`set_color_scheme`.
"""

from __future__ import annotations

import base64
import io
import warnings
from itertools import chain
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

# ── plotly ──────────────────────────────────────────────────────
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# ── sklearn / rdkit (optional) ─────────────────────────────────
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
    from igraph import Graph  # type: ignore[import-untyped]
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
_DEFAULT_DPI: int = 96  # screen dpi for px sizing
_DEFAULT_FONTSIZE: int = 12


def _inch_to_px(inches: float) -> int:
    return int(inches * _DEFAULT_DPI)


def _figsize_to_px(
    figure_size: Optional[Tuple[float, float]] = None,
) -> Tuple[int, int]:
    fig = figure_size or _DEFAULT_FIGSIZE
    return (_inch_to_px(fig[0]), _inch_to_px(fig[1]))


# ═══════════════════════════════════════════════════════════════
#  Public helpers – colour & size
# ═══════════════════════════════════════════════════════════════

def set_color_scheme(scheme: Dict[str, Sequence[str]]) -> None:
    """Replace the global colour palette."""
    global _CURRENT_COLORS
    for key in ("sequential", "diverging", "qualitative"):
        if key in scheme:
            _CURRENT_COLORS[key] = list(scheme[key])


def set_figure_size(
    size: Tuple[float, float] = (5.0, 5.0),
    dpi: int = 96,
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

def _palette(color_type: str = "qualitative") -> List[str]:
    """Return the current palette for a given colour family."""
    return _CURRENT_COLORS.get(color_type, _CURRENT_COLORS["qualitative"])


def _base_layout(
    title: str = "",
    x_title: str = "",
    y_title: str = "",
    figure_size: Optional[Tuple[float, float]] = None,
    legend_title: str = "",
    show_legend: bool = True,
) -> go.Layout:
    """Build a reusable Plotly layout."""
    w, h = _figsize_to_px(figure_size)
    return go.Layout(
        title=dict(text=title, font=dict(size=_DEFAULT_FONTSIZE * 1.3)),
        width=w,
        height=h,
        font=dict(size=_DEFAULT_FONTSIZE),
        plot_bgcolor="white",
        paper_bgcolor="white",
        xaxis=dict(
            title=x_title,
            showgrid=True,
            gridcolor="#eeeeee",
            linecolor="black",
            linewidth=1,
            ticks="outside",
        ),
        yaxis=dict(
            title=y_title,
            showgrid=True,
            gridcolor="#eeeeee",
            linecolor="black",
            linewidth=1,
            ticks="outside",
        ),
        legend=dict(
            title=dict(text=legend_title) if legend_title else None,
            visible=show_legend,
            xanchor="left",
            yanchor="top",
            x=0.98,
            y=0.98,
            bgcolor="rgba(0,0,0,0)",
        ),
        margin=dict(l=60, r=40, t=60, b=60),
    )


def _save_plotly(fig: go.Figure, save_to: Optional[str]) -> None:
    if save_to is None:
        return
    path = str(save_to).lower()
    if path.endswith(".html"):
        fig.write_html(save_to)
    elif path.endswith(".json"):
        fig.write_json(save_to)
    else:
        fig.write_image(save_to)


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
    """Horizontal bar chart (Plotly)."""
    if x not in df.columns:
        raise ValueError(f"Column '{x}' not found in DataFrame.")

    data = df.copy()
    data[x] = data[x].fillna(na_replace).astype(str)

    if y is None:
        freq = data[x].value_counts().reset_index()
        freq.columns = [x, "count"]
        freq = freq.sort_values("count", ascending=True)
        y_col = "count"
        y_label = "Frequency"
    else:
        if y not in data.columns:
            raise ValueError(f"Column '{y}' not found in DataFrame.")
        freq = data.sort_values(y, ascending=True)
        y_col = y
        y_label = y

    colors = _palette(color_type)

    if fill is not None and fill in freq.columns:
        fig = px.bar(
            freq,
            x=y_col,
            y=x,
            color=fill,
            color_discrete_sequence=colors,
            orientation="h",
        )
    else:
        fig = px.bar(
            freq,
            x=y_col,
            y=x,
            color_discrete_sequence=[colors[0]],
            orientation="h",
        )

    fig.update_layout(
        _base_layout(
            title=title,
            x_title=y_label,
            y_title=x,
            figure_size=figure_size,
        )
    )
    fig.update_yaxes(categoryorder="total ascending")
    _save_plotly(fig, save_to)
    return fig


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
    """Faceted boxplot from long-format data (Plotly)."""
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

    if facet is not None:
        n_facets = df_plot[facet].nunique()
        n_rows = int(np.ceil(n_facets / ncol))
        w, h = _figsize_to_px(figure_size)
        facet_size = (w * ncol, h * n_rows)
        fig = px.box(
            df_plot,
            x=x,
            y=y,
            color=fill_col,
            facet_col=facet,
            facet_col_wrap=ncol,
            color_discrete_sequence=colors,
            points="all",
        )
        fig.update_layout(
            _base_layout(
                x_title="Group",
                y_title=y_label,
                figure_size=facet_size,
            )
        )
        if show_trend:
            summary_df = (
                df_plot.groupby([facet, x], observed=True)[y]
                .median()
                .reset_index()
                .rename(columns={y: "summary_value"})
            )
            for facet_val in summary_df[facet].unique():
                sub = summary_df[summary_df[facet] == facet_val]
                fig.add_trace(
                    go.Scatter(
                        x=sub[x],
                        y=sub["summary_value"],
                        mode="lines",
                        line=dict(color="grey", width=1),
                        opacity=0.6,
                        showlegend=False,
                        xaxis=fig.layout.xaxis.domain if hasattr(fig.layout, "xaxis") else None,
                        yaxis=fig.layout.yaxis.domain if hasattr(fig.layout, "yaxis") else None,
                    )
                )
    else:
        fig = px.box(
            df_plot,
            x=x,
            y=y,
            color=fill_col,
            color_discrete_sequence=colors,
            points="all",
        )
        fig.update_layout(
            _base_layout(
                x_title="Group",
                y_title=y_label,
                figure_size=figure_size,
            )
        )

    fig.update_xaxes(tickangle=45)
    _save_plotly(fig, save_to)
    return fig


# ═══════════════════════════════════════════════════════════════
#  3. Donut chart
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
    """Donut / ring chart (Plotly)."""
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

    fig = go.Figure(
        data=[
            go.Pie(
                labels=freq["label"],
                values=freq["count"],
                hole=0.4,
                marker=dict(colors=colors[: len(freq)], line=dict(color="white", width=2)),
                textinfo="label+percent" if show_percent else "label",
                textposition="outside",
                rotation=90,
            )
        ]
    )
    fig.update_layout(
        _base_layout(title=title, figure_size=figure_size, show_legend=False)
    )
    _save_plotly(fig, save_to)
    return fig


# ═══════════════════════════════════════════════════════════════
#  4. Heatmap
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
    """Clustered heatmap (Plotly)."""
    if data_transfer == "log10":
        data = df.apply(lambda x: np.log(x + 1))
        cbar_label = "log10 (peak area)"
    elif data_transfer == "relative":
        base = df.max().max()
        data = (df / base) * 100 if base > 0 else df
        cbar_label = "Relative intensity (%)"
    elif data_transfer is None or data_transfer == "":
        data = df
        cbar_label = "Intensity"
    else:
        raise ValueError(
            f"Unknown data_transfer: {data_transfer}. Use 'log10', 'relative', or None."
        )

    n_rows, n_cols = data.shape
    w, h = _figsize_to_px(figure_size) or (_inch_to_px(max(5, n_cols * 0.5)), _inch_to_px(max(5, n_rows * 0.5)))

    colors = _palette(color_type)

    fig = go.Figure(
        data=go.Heatmap(
            z=data.values,
            x=data.columns.tolist(),
            y=data.index.tolist(),
            colorscale=[(i / (len(colors) - 1), c) for i, c in enumerate(colors)] if len(colors) > 1 else [(0, colors[0]), (1, colors[0])],
            colorbar=dict(title=cbar_label, len=0.8),
            hoverongaps=False,
        )
    )
    fig.update_layout(
        _base_layout(title=title, x_title=xlab, y_title=ylab, figure_size=(w, h), show_legend=False)
    )
    fig.update_xaxes(side="bottom")
    _save_plotly(fig, save_to)
    return fig


# ═══════════════════════════════════════════════════════════════
#  5. Line with error band
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
    """Line plot with optional error ribbon (Plotly)."""
    colors = _palette(color_type)
    line_color = colors[0]
    ribbon_color = colors[1] if len(colors) > 1 else colors[0]

    if color is not None:
        fig = px.line(df, x=x, y=y, color=color, color_discrete_sequence=colors)
    elif group is not None:
        fig = px.line(df, x=x, y=y, line_group=group, color_discrete_sequence=[line_color])
    else:
        fig = px.line(df, x=x, y=y, color_discrete_sequence=[line_color])

    if show_error and "Error" in df.columns:
        fig.add_trace(
            go.Scatter(
                x=df[x].tolist() + df[x].tolist()[::-1],
                y=(df[y] + df["Error"]).tolist() + (df[y] - df["Error"]).tolist()[::-1],
                fill="toself",
                fillcolor=ribbon_color,
                line=dict(color="rgba(0,0,0,0)"),
                opacity=0.15,
                showlegend=False,
                hoverinfo="skip",
            )
        )

    fig.update_layout(
        _base_layout(
            title=f"{y} by {x}",
            x_title=x,
            y_title=y,
            figure_size=figure_size,
            show_legend=color is not None,
        )
    )
    fig.update_layout(plot_bgcolor="white", xaxis_showgrid=False, yaxis_showgrid=False)
    _save_plotly(fig, save_to)
    return fig


# ═══════════════════════════════════════════════════════════════
#  6. Lollipop chart
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
    """Horizontal lollipop chart (Plotly)."""
    data = df.sort_values(by=x, ascending=True).copy()

    colors = _palette(color_type)
    color_col = fill if fill is not None and fill in data.columns else None

    fig = go.Figure()
    for i, row in data.iterrows():
        c = colors[int(row.get(color_col, 0)) % len(colors)] if color_col is not None else colors[0]
        fig.add_trace(
            go.Scatter(
                x=[0, row[x]],
                y=[row[y], row[y]],
                mode="lines",
                line=dict(color="grey", width=1),
                showlegend=False,
                hoverinfo="skip",
            )
        )
        fig.add_trace(
            go.Scatter(
                x=[row[x]],
                y=[row[y]],
                mode="markers",
                marker=dict(color=c, size=10, line=dict(color="black", width=1)),
                showlegend=False,
                hovertemplate=f"{y}=%{{y}}<br>{x}=%{{x}}<extra></extra>",
            )
        )

    fig.update_layout(
        _base_layout(
            x_title=x,
            y_title=y,
            figure_size=figure_size,
            show_legend=False,
        )
    )
    fig.update_yaxes(type="category", categoryorder="array", categoryarray=data[y].tolist())
    _save_plotly(fig, save_to)
    return fig


# ═══════════════════════════════════════════════════════════════
#  7. PCA score plot
# ═══════════════════════════════════════════════════════════════

def pca(
    df: pd.DataFrame,
    groups: Optional[Sequence] = None,
    labels: Optional[Sequence[str]] = None,
    color_type: str = "qualitative",
    palette: Optional[str] = None,
    add_ellipse: bool = True,
    data_transfer: Optional[str] = "log10",
    figure_size: Optional[Tuple[float, float]] = None,
    save_to: Optional[str] = None,
):
    """PCA score scatter plot (Plotly)."""
    if palette is not None:
        color_type = palette

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
        raise ValueError(f"Unknown data_transfer: {data_transfer}. Use 'log10', 'relative', or None.")

    pca_model = PCA(n_components=2).fit(df_plot)
    scores = pca_model.transform(df_plot)
    pca_df = pd.DataFrame(scores, columns=["PC1", "PC2"])

    if groups is not None:
        pca_df["group"] = pd.Categorical(groups).astype(str)
    else:
        pca_df["group"] = "Group"

    if labels is not None:
        pca_df["label"] = labels

    colors = _palette(color_type)
    fig = px.scatter(
        pca_df,
        x="PC1",
        y="PC2",
        color="group",
        text="label" if labels is not None else None,
        color_discrete_sequence=colors,
    )

    if add_ellipse and groups is not None:
        for grp in pca_df["group"].unique():
            sub = pca_df[pca_df["group"] == grp][["PC1", "PC2"]]
            if len(sub) < 3:
                continue
            cov = np.cov(sub.T)
            vals, vecs = np.linalg.eigh(cov)
            order = vals.argsort()[::-1]
            vals, vecs = vals[order], vecs[:, order]
            theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
            width, height = 2 * np.sqrt(vals) * 2.448  # 95 % confidence
            fig.add_shape(
                type="ellipse",
                x0=sub["PC1"].mean() - width / 2,
                y0=sub["PC2"].mean() - height / 2,
                x1=sub["PC1"].mean() + width / 2,
                y1=sub["PC2"].mean() + height / 2,
                xref="x",
                yref="y",
                line_color=colors[list(pca_df["group"].unique()).index(grp) % len(colors)],
                fillcolor=colors[list(pca_df["group"].unique()).index(grp) % len(colors)],
                opacity=0.15,
            )

    fig.update_layout(
        _base_layout(
            x_title=f"PC1: {100 * pca_model.explained_variance_ratio_[0]:.1f} %",
            y_title=f"PC2: {100 * pca_model.explained_variance_ratio_[1]:.1f} %",
            figure_size=figure_size,
        )
    )
    fig.update_traces(textposition="top center")
    _save_plotly(fig, save_to)
    return fig


# ═══════════════════════════════════════════════════════════════
#  8. PLS-DA score plot
# ═══════════════════════════════════════════════════════════════

def plsda(
    T_scores: np.ndarray,
    y: Sequence,
    color_type: str = "qualitative",
    palette: Optional[str] = None,
    data_transfer: Optional[str] = "log10",
    figure_size: Tuple[float, float] = (4, 4),
    save_to: Optional[str] = None,
):
    """PLS-DA score scatter with decision regions (Plotly)."""
    if palette is not None:
        color_type = palette

    if not _HAS_SKLEARN:
        raise ImportError("scikit-learn is required for PLS-DA plots.")

    if data_transfer == "log10":
        T_plot = np.sign(T_scores) * np.log10(np.abs(T_scores) + 1)
    elif data_transfer == "relative":
        base = np.max(np.abs(T_scores))
        T_plot = (T_scores / base) * 100 if base > 0 else T_scores.copy()
    elif data_transfer is None or data_transfer == "":
        T_plot = T_scores.copy()
    else:
        raise ValueError(f"Unknown data_transfer: {data_transfer}. Use 'log10', 'relative', or None.")

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
        np.linspace(x_min, x_max, 200),
        np.linspace(y_min, y_max, 200),
    )
    grid = np.c_[xx.ravel(), yy.ravel()]
    Z = clf.predict(grid).reshape(xx.shape)

    colors = _palette(color_type)
    palette_map = {c: colors[i % len(colors)] for i, c in enumerate(class_names)}

    fig = go.Figure()

    # decision-region contours
    fig.add_trace(
        go.Contour(
            x=np.linspace(x_min, x_max, 200),
            y=np.linspace(y_min, y_max, 200),
            z=Z,
            colorscale=[(i / (len(class_names) - 1), palette_map[c]) for i, c in enumerate(class_names)] if len(class_names) > 1 else [(0, palette_map[class_names[0]]), (1, palette_map[class_names[0]])],
            opacity=0.25,
            showscale=False,
            contours=dict(coloring="fill"),
            hoverinfo="skip",
        )
    )

    for grp in class_names:
        sub = df_scores[df_scores["group"] == grp]
        fig.add_trace(
            go.Scatter(
                x=sub["LV1"],
                y=sub["LV2"],
                mode="markers",
                name=str(grp),
                marker=dict(color=palette_map[grp], size=8),
            )
        )

    fig.update_layout(
        _base_layout(
            title="PLS-DA scores with decision regions",
            x_title="LV1 (PLS-DA)",
            y_title="LV2 (PLS-DA)",
            figure_size=figure_size,
        )
    )
    fig.update_xaxes(range=[x_min, x_max])
    fig.update_yaxes(range=[y_min, y_max])
    _save_plotly(fig, save_to)
    return fig


# ═══════════════════════════════════════════════════════════════
#  9. Volcano plot
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
    palette: Optional[str] = None,
    figure_size: Optional[Tuple[float, float]] = None,
    save_to: Optional[str] = None,
):
    """Volcano plot (Plotly)."""
    if palette is not None:
        color_type = palette

    colors = _palette(color_type)
    if len(colors) >= 2:
        up_color = colors[0]
        dn_color = colors[-1]
    else:
        up_color = dn_color = colors[0]
    no_color = "#D3D3D3"

    x_vals = df[x].replace([np.inf, -np.inf], np.nan).dropna()
    x_limit = max(abs(x_vals.min()), abs(x_vals.max())) if not x_vals.empty else 1.0

    y_vals = df[y].replace([np.inf, -np.inf], np.nan).dropna()
    y_max = max(2.5, y_vals.max()) if not y_vals.empty else 2.5

    color_map = {"up": up_color, "dn": dn_color, "no": no_color}
    df["_color"] = df[fill].map(color_map)

    fig = go.Figure()
    for category in df[fill].unique():
        sub = df[df[fill] == category]
        fig.add_trace(
            go.Scatter(
                x=sub[x],
                y=sub[y],
                mode="markers",
                name=str(category),
                marker=dict(color=color_map.get(str(category), no_color), size=6, opacity=0.7),
            )
        )

    fig.add_vline(x=-xcut, line=dict(dash="dash", color="grey", width=1))
    fig.add_vline(x=xcut, line=dict(dash="dash", color="grey", width=1))
    fig.add_hline(y=ycut, line=dict(dash="dash", color="grey", width=1))

    fig.update_layout(
        _base_layout(
            title=f"{title}, Volcano Plot" if title else "Volcano Plot",
            x_title=xlab if xlab is not None else x,
            y_title=ylab if ylab is not None else y,
            figure_size=figure_size,
        )
    )
    fig.update_xaxes(range=[-x_limit, x_limit])
    fig.update_yaxes(range=[0, y_max])
    _save_plotly(fig, save_to)
    return fig


# ═══════════════════════════════════════════════════════════════
#  10. Scatter / TIC-RT-mz style
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
    """Generic scatter plot (Plotly)."""
    colors = _palette(color_type)

    if shape is not None:
        df = df.copy()
        df[shape] = df[shape].astype(str)

    fig = px.scatter(
        df,
        x=x,
        y=y,
        size=size,
        color=color,
        symbol=shape,
        color_discrete_sequence=colors,
        opacity=alpha,
    )

    if color is None:
        fig.update_traces(marker=dict(color=colors[0]))

    fig.update_layout(
        _base_layout(
            x_title=x,
            y_title=y,
            figure_size=figure_size,
            show_legend=color is not None or shape is not None,
        )
    )
    _save_plotly(fig, save_to)
    return fig


# ═══════════════════════════════════════════════════════════════
#  11. Venn diagram
# ═══════════════════════════════════════════════════════════════

def venn(
    data: Dict[str, Sequence],
    fill_mode: Optional[List[str]] = None,
    color_type: str = "qualitative",
    palette: Optional[str] = None,
    alpha: float = 0.65,
    fontsize: Optional[int] = None,
    figure_size: Optional[Tuple[float, float]] = None,
    save_to: Optional[str] = None,
):
    """Draw a 2-, 3- or 4-set Venn diagram (Plotly shapes)."""
    if palette is not None:
        color_type = palette

    if fill_mode is None:
        fill_mode = ["number"]
    if not isinstance(data, dict):
        raise ValueError("Input data must be a dictionary.")
    n_sets = len(data)
    if n_sets not in (2, 3, 4):
        raise ValueError("Venn diagram supports 2, 3 or 4 sets only.")

    colors = _palette(color_type)
    fs = fontsize or _DEFAULT_FONTSIZE

    labels = _venn_divide(list(data.values()), fill=fill_mode)
    keys = list(data.keys())

    fig = go.Figure()
    fig.update_layout(
        _base_layout(figure_size=figure_size, show_legend=False),
        xaxis=dict(range=[0, 1], showgrid=False, zeroline=False, visible=False),
        yaxis=dict(range=[0, 1], showgrid=False, zeroline=False, visible=False),
        margin=dict(l=20, r=20, t=20, b=20),
    )

    if n_sets == 2:
        _venn_draw2_plotly(fig, labels, keys, colors, alpha, fs)
    elif n_sets == 3:
        _venn_draw3_plotly(fig, labels, keys, colors, alpha, fs)
    else:
        _venn_draw4_plotly(fig, labels, keys, colors, alpha, fs)

    _save_plotly(fig, save_to)
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


def _venn_add_ellipse(fig, x, y, width, height, angle, color, alpha):
    fig.add_shape(
        type="ellipse",
        x0=x - width / 2,
        y0=y - height / 2,
        x1=x + width / 2,
        y1=y + height / 2,
        xref="x",
        yref="y",
        fillcolor=color,
        line_color=color,
        opacity=alpha,
    )


def _venn_add_text(fig, x, y, text, fs):
    fig.add_annotation(x=x, y=y, text=text, showarrow=False, font=dict(size=fs))


def _venn_draw2_plotly(fig, labels, names, colors, alpha, fs):
    _venn_add_ellipse(fig, 0.375, 0.30, 0.5, 0.5, 0, colors[0], alpha)
    _venn_add_ellipse(fig, 0.625, 0.30, 0.5, 0.5, 0, colors[1], alpha)
    _venn_add_text(fig, 0.74, 0.30, labels.get("01", ""), fs)
    _venn_add_text(fig, 0.26, 0.30, labels.get("10", ""), fs)
    _venn_add_text(fig, 0.50, 0.30, labels.get("11", ""), fs)
    _venn_add_text(fig, 0.20, 0.56, names[0], fs)
    _venn_add_text(fig, 0.80, 0.56, names[1], fs)


def _venn_draw3_plotly(fig, labels, names, colors, alpha, fs):
    _venn_add_ellipse(fig, 0.333, 0.633, 0.55, 0.55, 0, colors[0], alpha)
    _venn_add_ellipse(fig, 0.666, 0.633, 0.55, 0.55, 0, colors[1], alpha)
    _venn_add_ellipse(fig, 0.500, 0.310, 0.55, 0.55, 0, colors[2], alpha)
    _venn_add_text(fig, 0.50, 0.27, labels.get("001", ""), fs)
    _venn_add_text(fig, 0.73, 0.65, labels.get("010", ""), fs)
    _venn_add_text(fig, 0.61, 0.46, labels.get("011", ""), fs)
    _venn_add_text(fig, 0.27, 0.65, labels.get("100", ""), fs)
    _venn_add_text(fig, 0.39, 0.46, labels.get("101", ""), fs)
    _venn_add_text(fig, 0.50, 0.65, labels.get("110", ""), fs)
    _venn_add_text(fig, 0.50, 0.51, labels.get("111", ""), fs)
    _venn_add_text(fig, 0.15, 0.87, names[0], fs)
    _venn_add_text(fig, 0.85, 0.87, names[1], fs)
    _venn_add_text(fig, 0.50, 0.02, names[2], fs)


def _venn_draw4_plotly(fig, labels, names, colors, alpha, fs):
    o, dx, dy = 0.500, 0.18, 0.08
    _venn_add_ellipse(fig, o - dx, o - dy, 4 * dx, 2 * dx, 135, colors[0], alpha)
    _venn_add_ellipse(fig, o, o, 4 * dx, 2 * dx, 135, colors[1], alpha)
    _venn_add_ellipse(fig, o, o, 4 * dx, 2 * dx, 45, colors[2], alpha)
    _venn_add_ellipse(fig, o + dx, o - dy, 4 * dx, 2 * dx, 45, colors[3], alpha)
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
        _venn_add_text(fig, x, y, text, fs)


# ═══════════════════════════════════════════════════════════════
#  12. Colour swatch
# ═══════════════════════════════════════════════════════════════

def swatch(
    colors: Optional[Sequence[str]] = None,
    color_type: str = "qualitative",
    figure_size: Optional[Tuple[float, float]] = None,
    save_to: Optional[str] = None,
):
    """Draw colour swatches (Plotly)."""
    if colors is None:
        colors = _palette(color_type)
    colors = list(colors)

    fig = go.Figure()
    for i, c in enumerate(colors):
        fig.add_trace(
            go.Bar(
                x=[c],
                y=[1],
                marker_color=c,
                showlegend=False,
                hovertemplate=f"{c}<extra></extra>",
            )
        )

    w = _inch_to_px(len(colors) * 0.5)
    h = _inch_to_px(1.0)
    fig.update_layout(
        _base_layout(title="Color Swatches", figure_size=figure_size or (w, h), show_legend=False),
        barmode="group",
        bargap=0.1,
        xaxis=dict(type="category", categoryorder="array", categoryarray=colors),
        yaxis=dict(visible=False, range=[0, 1]),
    )
    _save_plotly(fig, save_to)
    return fig


# ═══════════════════════════════════════════════════════════════
#  13. Chemical network
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
    """Build a chemical-similarity network from SMILES (Plotly)."""
    if not _HAS_RDKIT:
        raise ImportError("rdkit is required for chemical network plots.")

    if not smiles_list:
        raise ValueError("smiles_list cannot be empty")

    mols = []
    fps = []
    if fp_type.lower() == "morgan":
        morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
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

    # Fruchterman-Reingold layout via igraph
    if _HAS_IGRAPH:
        g = Graph()
        g.add_vertices(n)
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
        pos = {i: (layout[i][0], layout[i][1]) for i in range(n)}
    else:
        # fallback: random circular layout
        np.random.seed(42)
        angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
        pos = {i: (0.5 + 0.4 * np.cos(a), 0.5 + 0.4 * np.sin(a)) for i, a in enumerate(angles)}

    fig = go.Figure()

    # edges
    for i in range(n):
        for j in range(i + 1, n):
            if sim_matrix[i, j] >= threshold:
                x0, y0 = pos[i]
                x1, y1 = pos[j]
                fig.add_trace(
                    go.Scatter(
                        x=[x0, x1, None],
                        y=[y0, y1, None],
                        mode="lines",
                        line=dict(color="grey", width=sim_matrix[i, j] * 3),
                        showlegend=False,
                        hoverinfo="skip",
                    )
                )

    # nodes
    node_x = [pos[i][0] for i in range(n)]
    node_y = [pos[i][1] for i in range(n)]
    fig.add_trace(
        go.Scatter(
            x=node_x,
            y=node_y,
            mode="markers+text",
            marker=dict(size=12, color="#d7191c"),
            text=[str(i) for i in range(n)],
            textposition="top center",
            hovertemplate="index: %{text}<extra></extra>",
            showlegend=False,
        )
    )

    fig.update_layout(
        _base_layout(title="Chemical Similarity Network", show_legend=False),
        xaxis=dict(visible=False, range=[-0.1, 1.1]),
        yaxis=dict(visible=False, range=[-0.1, 1.1], scaleanchor="x"),
        margin=dict(l=20, r=20, t=40, b=20),
    )
    _save_plotly(fig, save_to)
    return fig


# ═══════════════════════════════════════════════════════════════
#  14. Molecule grid
# ═══════════════════════════════════════════════════════════════

def molecules_grid(
    smiles_list: List[str],
    cols: int = 4,
    mols_per_row: Optional[int] = None,
    sub_img_size: Tuple[int, int] = (250, 250),
    legends: bool = False,
    save_to: Optional[str] = None,
):
    """Render a grid of 2-D molecular structures (RDKit + Plotly)."""
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

    # Return a Plotly figure containing the SVG as an image
    svg_str = img.data if hasattr(img, "data") else str(img)
    encoded = base64.b64encode(svg_str.encode("utf-8")).decode("utf-8")

    n_rows = int(np.ceil(len(mols) / mols_per_row))
    fig_w = mols_per_row * sub_img_size[0]
    fig_h = n_rows * sub_img_size[1]

    fig = go.Figure()
    fig.add_layout_image(
        dict(
            source=f"data:image/svg+xml;base64,{encoded}",
            x=0,
            y=1,
            sizex=1,
            sizey=1,
            sizing="contain",
            xref="paper",
            yref="paper",
            layer="below",
        )
    )
    fig.update_layout(
        _base_layout(show_legend=False),
        width=fig_w,
        height=fig_h,
        xaxis=dict(visible=False, range=[0, 1]),
        yaxis=dict(visible=False, range=[0, 1]),
        margin=dict(l=0, r=0, t=0, b=0),
    )

    if save_to:
        with open(save_to, "w", encoding="utf-8") as f:
            f.write(svg_str)
    return fig


# ═══════════════════════════════════════════════════════════════
#  15. Centroid spectrum
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
    """Draw a centroid mass spectrum (Plotly lollipop / stem plot)."""
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

    fig = go.Figure()

    # stem lines
    for xi, yi in zip(mz, intensity):
        fig.add_trace(
            go.Scatter(
                x=[xi, xi],
                y=[0, yi],
                mode="lines",
                line=dict(color=color, width=line_width),
                showlegend=False,
                hoverinfo="skip",
            )
        )

    # markers at top
    fig.add_trace(
        go.Scatter(
            x=mz,
            y=intensity,
            mode="markers",
            marker=dict(color=color, size=4),
            showlegend=False,
            hovertemplate="m/z: %{x:.4f}<br>Intensity: %{y:.4f}<extra></extra>",
        )
    )

    if show_mz_labels and len(intensity) > 0:
        threshold = mz_label_intensity_threshold * float(np.max(intensity))
        y_range = float(np.max(intensity) - np.min(intensity))
        y_offset = y_range * mz_label_offset_ratio if y_range > 0 else 0.01
        annotations = []
        for xi, yi in zip(mz, intensity):
            if xi < xlim_eff[0] or xi > xlim_eff[1]:
                continue
            if yi < threshold:
                continue
            annotations.append(
                dict(
                    x=xi + mz_label_horizontal_offset,
                    y=yi + y_offset,
                    text=mz_label_fmt.format(xi),
                    showarrow=False,
                    font=dict(size=mz_label_fontsize, color=color),
                )
            )
        fig.update_layout(annotations=annotations)

    fig.update_layout(
        _base_layout(
            title=title or "",
            x_title="m/z",
            y_title="Intensity" if not normalize else "Normalized Intensity",
            figure_size=figsize,
            show_legend=False,
        )
    )
    fig.update_xaxes(range=xlim_eff)
    fig.update_yaxes(range=ylim_eff)

    if save_to is not None:
        _save_plotly(fig, save_to)
    if not show:
        fig.show = lambda *a, **k: None  # type: ignore[method-assign]
    return fig


# ═══════════════════════════════════════════════════════════════
#  Quick smoke-test when executed directly
# ═══════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("mzpy.plotly smoke test")
    print("-" * 40)
    print("Active colour scheme:")
    for k, v in get_color_scheme().items():
        print(f"  {k:12s}: {v[:3]} ...")
    print("\nTry:  from mzpy.plotly import bar, venn, ...")
