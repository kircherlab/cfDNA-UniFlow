#!/usr/bin/env python

import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.ndimage import median_filter
from matplotlib.gridspec import GridSpec


def load_table(path: str):
    df = pd.read_csv(path, sep="\t", index_col=[0, 1])
    return df


def GC_plot_with_read_distribution(
    df, sample_name, figsize=(35, 30), filter=None, quantile_threshold=None
):
    fig = plt.figure(figsize=figsize)
    R_GC = df.loc["R_gc"]
    # apply filter for smoothing
    if filter:
        R_GC_index = R_GC.index
        R_GC = pd.DataFrame(median_filter(R_GC, size=filter))
        R_GC.index = R_GC_index

    # calculate min,max values for coloring
    if quantile_threshold:
        R_GC_min, R_GC_max = np.nanmin(R_GC), np.nanquantile(R_GC, quantile_threshold)
    else:
        R_GC_min, R_GC_max = np.nanmin(R_GC), np.nanmax(R_GC)

    F_GC = df.loc["F_gc"]
    # reads = F_GC.sum().sum()
    percent_reads_per_length = (
        F_GC.sum(axis=1) / F_GC.sum(axis=1).sum()
    ) * 100  # .sort_index(ascending=False)
    percent_reads_per_gc = (
        F_GC.sum(axis=0) / F_GC.sum(axis=0).sum()
    ) * 100  # .sort_index(ascending=False)
    fig.suptitle(f"{sample_name}")

    # create axes

    gs = GridSpec(nrows=2, ncols=3, width_ratios=(0.3, 4, 1), height_ratios=(1, 4))
    ax = fig.add_subplot(gs[1, 1])
    ax_cbar = fig.add_subplot(gs[1, 0])
    ax_marginx = fig.add_subplot(gs[0, 1], sharex=ax)
    ax_marginy = fig.add_subplot(gs[1, 2], sharey=ax)

    # create plots
    bias = ax.pcolormesh(
        R_GC.columns,
        R_GC.index,
        R_GC.values,
        cmap="viridis",
        vmin=R_GC_min,
        vmax=R_GC_max,
    )
    ax.set_ylabel("fragment length")
    ax.set_xticks(
        ticks=np.arange(
            min(R_GC.columns.astype(int)), max(R_GC.columns.astype(int)) + 1, 10
        )
    )
    ax.tick_params(axis="x", labelrotation=45)
    ax.set_xlabel("GC content [%]")

    gc_reads = ax_marginx.stackplot(
        percent_reads_per_gc.index, percent_reads_per_gc.values
    )
    ax_marginx.tick_params(axis="x", labelbottom=False)
    ax_marginx.set_ylabel("reads [%]")

    length_reads = ax_marginy.stackplot(
        percent_reads_per_length.values, percent_reads_per_length.index
    )
    # ax_marginy.grid()
    ax_marginy.tick_params(axis="y", labelleft=False)
    ax_marginy.set_ylim([30, 250])
    ax_marginy.set_xlabel("reads [%]")

    # plt.ylabel("fragment length")
    # plt.xticks(ticks=np.arange(min(R_GC.columns.astype(int)), max(R_GC.columns.astype(int)) + 1, 5))
    # plt.xlabel("GC content [%]")
    # plt.colorbar()
    # fig.suptitle(f"GC correction values for multiple samplesizes")

    fig.colorbar(bias, cax=fig.axes[1], ticklocation="left", label="GC bias ratio")
    fig.tight_layout()  # rect=[0, 0.03, 1, 0.95])
    # plt.tight_layout()
    return fig


@click.command()
@click.argument("input_file", type=click.Path(readable=True))
@click.option(
    "--sample_name",
    "-s",
    "sample_name",
    required=True,
    type=click.STRING,
    help="Name of the sample",
)
@click.option(
    "--quantile_threshold",
    "-q",
    "quantile_threshold",
    type=click.FLOAT,
    default=None,
    help="Quantile threshold for heatmap scale (e.g. 0.95)",
)
@click.option(
    "--filter",
    "filter",
    type=click.INT,
    default=None,
    show_default=True,
    help="""Filtersize for applying median filter on bias values.""",
)
@click.option(
    "--output",
    "-o",
    "output_file",
    required=True,
    type=click.Path(writable=True),
    help="Output_file (image format).",
)
@click.option(
    "--figsize",
    "figsize",
    type=(int, int),
    default=(12, 9),
    show_default=True,
    help="""Figsize of the output plot.""",
)
def main(input_file, sample_name, quantile_threshold, filter, output_file, figsize):
    """Loads GC_bias file (TSV with first two rows being a multiindex for category (N_GC,F_GC,R_GC) and fragmentlength),
       creates a heatmap for the calculated GC bias values and read distributions and saves the results as a combined plot.

    Args:
        input_file (tab separated file): TSV containing submatrices with expected fragments, measured fragments and bias value.
        sample_name (string): Name of the sample displayed in the plot.
        quantile_threshold (float): Quantile threshold for heatmap scale (e.g. 0.95), removing visual outliers skewing the scale.
        filter (integer): Filtersize for applying median filter on bias values.
        output_file (string): Output_file (image format understood by matplotlib).
        figsize (two integers): Figsize of the output plot.
    """
    # set figure details
    cm = 1 / 2.54
    sns.set_theme(
#        context="poster",
        style="darkgrid",
        palette="deep",
        font_scale=1.2,
        color_codes=True,
#        rc={"figure.figsize": (38 * cm, 32 * cm)},
    )
    # load data
    df = load_table(path=input_file)
    # make figure
    fig = GC_plot_with_read_distribution(
        df=df,
        sample_name=sample_name,
        figsize=figsize,
        filter=filter,
        quantile_threshold=quantile_threshold,
    )
    # save figure
    fig.get_figure().savefig(output_file, bbox_inches="tight")


if __name__ == "__main__":
    main()
