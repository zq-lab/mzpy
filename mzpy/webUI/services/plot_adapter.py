"""
Adapter layer: wrap mzpy.plot functions so they write PNGs into the
session plot folder and return a URL path.
"""

from __future__ import annotations

from typing import Optional

import pandas as pd

from .session_store import save_plot, save_plotnine


def plot_bar(
    df: pd.DataFrame,
    x: str,
    y: Optional[str] = None,
    fill: Optional[str] = None,
    title: str = "",
    color_type: str = "qualitative",
) -> str:
    from ...plot import plot_bar as _fn
    p = _fn(df, x=x, y=y, fill=fill, title=title, color_type=color_type)
    return save_plotnine("bar", p)


def plot_box(
    df: pd.DataFrame,
    x: str,
    y: str,
    fill: Optional[str] = None,
    facet: Optional[str] = None,
    color_type: str = "qualitative",
    show_trend: bool = True,
) -> str:
    from ...plot import plot_box as _fn
    p = _fn(df, x=x, y=y, fill=fill, facet=facet, color_type=color_type,
            show_trend=show_trend)
    return save_plotnine("box", p)


def plot_volcano(
    df: pd.DataFrame,
    x: str,
    y: str,
    fill: str,
    xcut: float = 1.0,
    ycut: float = 2.0,
    title: str = "",
    color_type: str = "diverging",
) -> str:
    from ...plot import plot_volcano as _fn
    p = _fn(df, x=x, y=y, fill=fill, xcut=xcut, ycut=ycut,
            title=title, color_type=color_type)
    return save_plotnine("volcano", p)


def plot_pca(
    df: pd.DataFrame,
    groups=None,
    labels=None,
    color_type: str = "qualitative",
    add_ellipse: bool = True,
    data_transfer: Optional[str] = "log10",
) -> str:
    from ...plot import plot_pca as _fn
    p = _fn(df, groups=groups, labels=labels,
            color_type=color_type, add_ellipse=add_ellipse,
            data_transfer=data_transfer)
    return save_plotnine("pca", p)


def plot_plsda(
    T_scores,
    y,
    color_type: str = "qualitative",
    data_transfer: Optional[str] = "log10",
) -> str:
    from ...plot import plot_plsda as _fn
    fig = _fn(T_scores, y, color_type=color_type,
              data_transfer=data_transfer)
    return save_plot("plsda", fig)


def plot_venn(
    data: dict,
    color_type: str = "qualitative",
    alpha: float = 0.65,
) -> str:
    from ...plot import plot_venn as _fn
    fig = _fn(data, color_type=color_type, alpha=alpha)
    return save_plot("venn", fig)


def plot_heatmap(
    df: pd.DataFrame,
    color_type: str = "diverging",
    title: str = "Heatmap",
) -> str:
    from ...plot import plot_heatmap as _fn
    fig = _fn(df, color_type=color_type, title=title)
    return save_plot("heatmap", fig)
