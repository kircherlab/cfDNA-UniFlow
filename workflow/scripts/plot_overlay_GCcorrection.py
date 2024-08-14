#!/usr/bin/env python

import click
import numpy as np
import pandas as pd
import seaborn as sns
from loguru import logger
from matplotlib import pyplot as plt
from scipy import stats
from scipy.signal import savgol_filter
import re


def file_name_func(path):
    """Takes a path as string and returns the filename without file extension.

    Args:
        path (str): Path to file.

    Returns:
        str: file name as ID used for downstream labeling
    """
    return path.split("/")[-1].split(".")[0]


def make_name_regex_func(name_regex):
    """Takes a regex pattern for extracting a label from the filename. Compiles pattern and returns closure function for searching in strings.

    Args:
        name_regex (str): Regex pattern.

    Returns:
        function: name_regex_function that extracts a substring from a file path.
    """
    pattern = re.compile(name_regex)

    def name_regex_func(text):
        return pattern.search(text).group(1)

    return name_regex_func


def load_table(path: str):
    """Loads table and sets index according to expected columns "ID" and optional column "coordinate"

    Args:
        path (str): path to the input file in .csv/.tsv format. Supports zipped files.

    Returns:
        pandas.DataFrame: DataFrame containing loaded data.
    """
    ## determine dtypes
    if any([x in path for x in [".tsv", ".tsv.gz"]]):
        dtype_df = pd.read_csv(path, sep="\t", header=None, nrows=1000)
    else:
        dtype_df = pd.read_csv(path, header=None, nrows=1000)

    non_numeric = dtype_df.select_dtypes(exclude=["float", "int"]).columns
    numeric = dtype_df.select_dtypes(include=["float", "int"]).columns
    dtype_dict = {key: str for key in non_numeric.values}
    dtype_dict.update({col: np.float32 for col in numeric if col not in non_numeric})
    del dtype_df
    ## load actual data
    if any([x in path for x in [".tsv", ".tsv.gz"]]):
        df = pd.read_csv(path, sep="\t", header=None, dtype=dtype_dict)
    else:
        df = pd.read_csv(path, header=None, dtype=dtype_dict)
    if pd.api.types.is_string_dtype(df[0]):
        if pd.api.types.is_string_dtype(df[1]):
            df = df.set_index([0, 1])
            df.index.names = ["ID", "coordinate"]
        else:
            df = df.set_index(0)
            df.index.name = "ID"
    return df


def calculate_flanking_regions(val: int):
    """Calculates flanking regions for point of interest.
    Args:
        val (int): should be length of value vector
    Raises:
        TypeError: Only integers are allowed
    Returns:
        [iterator]: range of values around center point (e.g. range(-1000,1000))
    """

    if not isinstance(val, int):
        raise TypeError("Only integers are allowed")

    if val % 2 == 0:
        flank = int(val / 2)
        region = range(-flank, flank)
    elif val % 2 == 1:
        flank_l = int(val / 2 - 0.5)
        flank_r = int(val / 2 + 0.5)
        region = range(-flank_l, flank_r)
    return region


def get_window_slice(window, width=1):
    """Generates upper and lower limits of a specified window from the center point for slicing data.

    Args:
        window (int): size of a window (e.g. 5000)
        width (int, optional): width of the window around the window center. Defaults to 1.

    Returns:
        tuple: upper and lower limits of a slice
    """
    halfwidth = width / 2
    center = window / 2
    start = center - halfwidth
    stop = center + halfwidth
    return int(start), int(stop)


def process_sample(
    sample: pd.DataFrame,
    overlay_mode: str = "mean",
    smoothing: bool = False,
    smooth_window: int = 21,
    smooth_polyorder: int = 2,
    scaling: bool = False,
    rolling: bool = False,
    rolling_window: int = 1000,
    detrend: str = "substract",
    edge_norm: bool = False,
    window: int = 1000,
    flank: int = 500,
) -> pd.DataFrame:
    """Process sample for plotting. This includes smoothing, normalization and the removal of edge regions to mitigate edge effects.

    Args:
        sample (pd.DataFrame): Sample signal extracted from cfDNA.
        overlay_mode (str, optional): Mode for aggregating signal over multiple regions. Can be mean, median or confidence interval. Defaults to "mean".
        smoothing (bool, optional): Switch on whether the signal should be smoothed. Defaults to False.
        smooth_window (int, optional): Smoothing window for Savitzky-Golay filter. Defaults to 21.
        smooth_polyorder (int, optional): Order of polynomial used in Savitzky-Golay filter. Defaults to 2.
        scaling (bool, optional): Switch on whether the signal should be standardized. Defaults to False.
        rolling (bool, optional): Switch on whether the signal should be demeaned by substracting the rolling median. Defaults to False.
        rolling_window (int, optional): Window size being used in rolling_window demeaning.
        detrend (str, optional): Mode for detrending. Can be substract or divide. Defaults to "substract".
        edge_norm (bool, optional): Switch on whether the signal should be demeaned by substracting the average signal in the flanking regions. Defaults to False.
        window (int, optional): Size of window centered on region of interest. Defaults to 1000.
        flank (int, optional): Size of flanking regions for edge_norm. Defaults to 500.

    Raises:
        ValueError: Invalid overlay_mode keyword.

    Returns:
        pd.DataFrame: Processed sample.
    """

    if edge_norm:
        flank_start, flank_end = get_window_slice(len(sample.columns), flank * 2)

        helper = abs(sample.iloc[:, flank_start:flank_end].mean(axis=1))
        if (helper == 0).any():
            zero_mask = helper == 0
            logger.info(
                f"Regions with zero mean in Normalization slice [{flank_start}, {flank_end}] encountered. Removing respective regions."
            )
            for zero_region in zero_mask[zero_mask].index:
                logger.info(f"Removing region: {zero_region}")
            helper.drop(helper[zero_mask].index, inplace=True)
            sample.drop(sample[zero_mask].index, inplace=True)
        sample = sample.div(helper, axis=0)

    if smoothing:
        sample = sample.apply(
            lambda x: savgol_filter(
                x, window_length=smooth_window, polyorder=smooth_polyorder
            ),
            axis=1,
            result_type="expand",
        )

    if overlay_mode.lower() == "mean":
        sample = pd.DataFrame(sample.mean(numeric_only=True))
    elif overlay_mode.lower() == "median":
        sample = pd.DataFrame(sample.median(numeric_only=True))
    else:
        raise ValueError(f"{overlay_mode} is not a valid keyword.")


    if rolling:
        trend = sample.rolling(rolling_window, center=True, min_periods=1).median()
        if detrend == "substract":
            sample = sample.sub(trend)
        elif detrend == "divide":
            sample = sample.div(trend)
        else:
            raise ValueError("{detrend} is not a valid keyword.".format(detrend=detrend))

    sample["position"] = calculate_flanking_regions(len(sample))
    sample = sample.set_index("position")

    return sample


def plot_correction_overlay(
    uncorrected: pd.DataFrame,
    corrected: pd.DataFrame,
    figsize: tuple = (12, 4),
    signal: str = "Coverage",
    target: str = "ROI",
    cmap: str = "tab20",
    lower_limit: float = None,
    upper_limit: float = None,
):
    ylab = f"Normalized {signal}"

    title = f"{target} {signal}"

    fig, ax = plt.subplots(1, 2, figsize=figsize, sharey="row")

    for col in uncorrected:
        ax[0].plot(uncorrected[col], alpha=0.7)

    for col in corrected:
        ax[1].plot(corrected[col], label=col, alpha=0.7)

    ax[0].set_title("Uncorrected coverage", fontsize=14)
    ax[1].set_title("GC corrected coverage", fontsize=14)
    ax[0].set_xlabel("Position relative to target site [bp]", fontsize=14)
    ax[1].set_xlabel("Position relative to target site [bp]", fontsize=14)
    ax[0].set_ylabel(ylab, fontsize=14)

    if lower_limit or upper_limit:
        ax0_ylim = ax[0].get_ylim()
        ax1_ylim = ax[1].get_ylim()
        shared_ylim = (min(ax0_ylim[0], ax1_ylim[0]), max(ax0_ylim[1], ax1_ylim[1]))
        if lower_limit:
            shared_ylim = (min(shared_ylim[0], lower_limit), shared_ylim[1])
        if upper_limit:
            shared_ylim = (shared_ylim[0], max(shared_ylim[1], upper_limit))

        ax[0].set_ylim(shared_ylim)
        ax[1].set_ylim(shared_ylim)

    fig.suptitle(title, fontsize=16)

    if len(ax[-1].get_legend_handles_labels()[1]) < 15:
        fig.legend(
            bbox_to_anchor=(1, 0.90),
            loc="upper left",
            title="Sample names [uncorrected, corrected]:",
        )
    else:
        fig.legend(
            bbox_to_anchor=(1, 0.90),
            loc="upper left",
            ncol=2,
            title="Sample names [uncorrected, corrected]:",
        )

    fig.tight_layout()

    return fig


@click.command()
@click.option(
    "--uncorrected_samples",
    "-us",
    "uncorrected_samples",
    required=True,
    type=click.Path(readable=True),
    multiple=True,
    help="""Path to the uncorrected signal files in .tsv format. 
            Supports gzipped formats by adding .gz file extension.
            Multiple files can be supplied by repeating the option (e.g. -us a.tsv -us b.tsv).""",
)
@click.option(
    "--corrected_samples",
    "-cs",
    "corrected_samples",
    required=True,
    type=click.Path(readable=True),
    multiple=True,
    help="""Path to the  GC corrected signal files in .tsv format. 
            Supports gzipped formats by adding .gz file extension.
            Multiple files can be supplied by repeating the option (e.g. -cs a.tsv -cs b.tsv).""",
)
@click.option(
    "--regex",
    "-r",
    "name_regex",
    type=click.STRING,
    help="""Simple regex for extracting label from filename. For example 'path/to/(.*?).tsv.gz' for file 'path/to/sample_corrected.tsv.gz'.""",
)
@click.option(
    "--output",
    "-o",
    "output",
    required=True,
    type=click.Path(writable=True),
    help="Path to the outputfile in .tsv format. Supports gzipped formats by adding .gz file extension.",
)
@click.option(
    "--signal",
    "signal",
    default="Coverage",
    show_default=True,
    type=click.STRING,
    help="Signal type (e.g. Coverage or WPS).",
)
@click.option(
    "--target",
    "-t",
    "target",
    default="Region of Interest",
    show_default=True,
    type=click.STRING,
    help="Name of the target region that will be used in the plot title.",
)
@click.option(
    "--overlay_mode",
    "-m",
    "overlay_mode",
    default="mean",
    show_default=True,
    type=click.Choice(["mean", "median"], case_sensitive=False),
    help="Sets overlay mode, specifying how regions should be aggregated for each sample.",
)
@click.option(
    "--smoothing",
    "-s",
    "smoothing",
    is_flag=True,
    show_default=True,
    default=False,
    help="Activates smoothing with Savitzky-Golay filter.",
)
@click.option(
    "--smooth_window",
    "-sw",
    "smooth_window",
    type=click.INT,
    default=21,
    show_default=True,
    help="Sets windows size used for smoothing with Savitzky-Golay filter.",
)
@click.option(
    "--smooth_polyorder",
    "-sp",
    "smooth_polyorder",
    type=click.INT,
    default=2,
    show_default=True,
    help="Sets order of polynomial used for smoothing with Savitzky-Golay filter.",
)
@click.option(
    "--rolling",
    "-r",
    "rolling",
    is_flag=True,
    show_default=True,
    default=False,
    help="Activates trend removal with a rolling median filter.",
)
@click.option(
    "--rolling_window",
    "-rw",
    "rolling_window",
    type=click.INT,
    default=200,
    show_default=True,
    help="Sets window size used in rolling median filter.",
)
@click.option(
    "--flank_norm",
    "-fn",
    "flank_norm",
    is_flag=True,
    show_default=True,
    default=False,
    help="Activates normalization by dividing the signals by the mean coverage in flanking intervals around the region of interest",
)
@click.option(
    "--flank",
    "-f",
    "flank",
    type=click.INT,
    default=2500,
    show_default=True,
    help="Sets the size of the flanking intervals around the region of interest. Should be <= 0.5 of the extracted signals.",
)
@click.option(
    "--display_window",
    "display_window",
    type=(int, int),
    default=(-1500, 1500),
    show_default=True,
    help="""Window centered around center of region of interest to plot.""",
)
@click.option(
    "--figsize",
    "figsize",
    type=(int, int),
    default=(12, 6),
    show_default=True,
    help="""Figsize of the output plot.""",
)
@click.option(
    "--lower_limit",
    "lower_limit",
    type=click.FLOAT,
    default=None,
    show_default=True,
    help="Sets the lower limit of the Y axis displayed in plotting.",
)
@click.option(
    "--upper_limit",
    "upper_limit",
    type=click.FLOAT,
    default=None,
    show_default=True,
    help="Sets the upper limit of the Y axis displayed in plotting.",
)
def main(
    uncorrected_samples,
    corrected_samples,
    name_regex,
    output,
    signal,
    target,
    overlay_mode,
    smoothing,
    smooth_window,
    smooth_polyorder,
    rolling,
    rolling_window,
    flank_norm,
    flank,
    display_window,
    figsize,
    lower_limit,
    upper_limit,
):
    if not name_regex:
        name_func = file_name_func
    else:
        name_func = make_name_regex_func(name_regex)

    if signal == "Coverage":
        detrend = "divide"
    elif signal == "WPS":
        detrend = "substract"

    logger.info("Loading uncorrected samples.")
    uncorrected_df = pd.DataFrame()
    for sample in uncorrected_samples:
        name = name_func(sample)
        logger.info(f"Loading uncorrected sample: {name}")
        tmpdf = load_table(sample)
        tmpdf = process_sample(
            sample=tmpdf,
            overlay_mode=overlay_mode,
            smoothing=smoothing,
            smooth_window=smooth_window,
            smooth_polyorder=smooth_polyorder,
            rolling=rolling,
            rolling_window=rolling_window,
            detrend=detrend,
            edge_norm=flank_norm,
            flank=flank,
        )
        uncorrected_df[name] = tmpdf
        del tmpdf

    logger.info("Loading corrected samples.")
    corrected_df = pd.DataFrame()
    for sample in corrected_samples:
        name = name_func(sample)
        logger.info(f"loading corrected sample: {name}")
        tmpdf = load_table(sample)
        tmpdf = process_sample(
            sample=tmpdf,
            overlay_mode=overlay_mode,
            smoothing=smoothing,
            smooth_window=smooth_window,
            smooth_polyorder=smooth_polyorder,
            rolling=rolling,
            detrend=detrend,
            edge_norm=flank_norm,
            flank=flank,
        )
        corrected_df[name] = tmpdf
        del tmpdf

    logger.info("Adjusting display window.")
    if display_window:
        uncorr_min = uncorrected_df.index.min()
        uncorr_max = uncorrected_df.index.max()
        corr_min = corrected_df.index.min()
        corr_max = corrected_df.index.max()

        uncorr_lower = max(uncorr_min, display_window[0])
        uncorr_upper = min(uncorr_max, display_window[1])
        corr_lower = max(corr_min, display_window[0])
        corr_upper = min(corr_max, display_window[1])

        uncorrected_df = uncorrected_df.loc[uncorr_lower:uncorr_upper]
        corrected_df = corrected_df.loc[corr_lower:corr_upper]

    sns.set_palette("hls", len(corrected_df.columns))

    logger.info("Plotting figure.")
    Fig = plot_correction_overlay(
        uncorrected=uncorrected_df,
        corrected=corrected_df,
        figsize=figsize,
        signal=signal,
        target=target,
        lower_limit=lower_limit,
        upper_limit=upper_limit,
    )
    logger.info("Saving figure.")
    Fig.savefig(output, bbox_inches="tight")


if __name__ == "__main__":
    main()
