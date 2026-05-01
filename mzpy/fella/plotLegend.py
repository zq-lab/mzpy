"""
Python port of FELLA:::plotLegend

Adds a legend to a graph plot.
"""

from typing import Optional

import matplotlib.pyplot as plt


def plotLegend(
    GO_annot: bool = False,
    cex: float = 1.2,
    markersize: float = 12,
    frameon: bool = True,
    fancybox: bool = True,
    ax: Optional[plt.Axes] = None,
):
    """
    Add legend to the plot.

    Parameters
    ----------
    GO_annot : bool
        Whether GO annotations are included.
    cex : float
        Relative font scale.
    markersize : float
        Marker size in the legend.
    frameon : bool
        Whether to draw a frame around the legend.
    fancybox : bool
        Whether to use a round box for the legend.
    ax : matplotlib Axes or None
        Axes to draw on.
    """
    if ax is None:
        ax = plt.gca()

    legend_elements = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#CD0000", markersize=markersize, label="Pathway"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#CD96CD", markersize=markersize, label="Module"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#FFD500", markersize=markersize, label="Enzyme"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#8DB6CD", markersize=markersize, label="Reaction"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#548B54", markersize=markersize, label="Compound"),
        plt.Line2D([0], [0], marker="s", color="w", markerfacecolor="#548B54", markersize=markersize, label="Input compound"),
    ]
    leg1 = ax.legend(
        handles=legend_elements,
        loc="lower left",
        title="Categories for each node",
        fontsize=cex * 10,
        title_fontsize=cex * 10,
        ncol=3,
        frameon=frameon,
        fancybox=fancybox,
    )
    ax.add_artist(leg1)

    if GO_annot:
        go_elements = [
            plt.Line2D([0], [0], marker="^", color="w", markerfacecolor="#FFD500", markersize=markersize, label="Simil < 0.5"),
            plt.Line2D([0], [0], marker="^", color="w", markerfacecolor="#FF5500", markersize=markersize, label="Simil < 0.7"),
            plt.Line2D([0], [0], marker="^", color="w", markerfacecolor="#FF0000", markersize=markersize, label="Simil < 0.9"),
            plt.Line2D([0], [0], marker="^", color="w", markerfacecolor="#B300FF", markersize=markersize, label="Simil <= 1"),
        ]
        ax.legend(
            handles=go_elements,
            loc="lower right",
            title="Enzymes with CC similarity",
            fontsize=cex * 10,
            title_fontsize=cex * 10,
            ncol=2,
            frameon=frameon,
            fancybox=fancybox,
        )
