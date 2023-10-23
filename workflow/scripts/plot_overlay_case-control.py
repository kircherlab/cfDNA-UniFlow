#!/usr/bin/env python

import click
import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
from loguru import logger
from matplotlib import pyplot as plt
from scipy import stats
from scipy.signal import savgol_filter


def load_controls(path: str):
    if any([x in path for x in [".tsv", ".tsv.gz"]]):
        df = pd.read_csv(path, sep="\t")
    else:
        df = pd.read_csv(path)

    if "position" in df.columns:
        df = df.set_index("position", drop=True)

    return df


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
    edge_norm: bool = False,
    window: int = 1000,
    flank:int = 500,
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
        edge_norm (bool, optional): Switch on whether the signal should be demeaned by substracting the average signal in the flanking regions. Defaults to False.
        window (int, optional): Size of window centered on region of interest. Defaults to 1000.
        flank (int, optional): Size of flanking regions for edge_norm. Defaults to 500.

    Raises:
        ValueError: Invalid overlay_mode keyword.

    Returns:
        pd.DataFrame: Processed sample.
    """

    #print(locals())
    if window % 2 == 0:
        fstart = int(window / 2 + 1)
        fstop = int(-window / 2)
    elif window % 2 == 1:
        fstart = int(window / 2 - 0.5 + 1)
        fstop = int(-window / 2 + 0.5)



    if edge_norm:
        flank_start, flank_end = get_window_slice(len(sample.columns),flank*2)

        helper = abs(
            sample.iloc[:, flank_start:flank_end].mean(axis=1)
        )
        if (helper==0).any():
            zero_mask = (helper==0)
            print(f"Regions with zero mean in Normalization slice [{flank_start}, {flank_end}] encountered. Removing respective regions.")
            for zero_region in zero_mask[zero_mask].index:
                print(f"Removing region: {zero_region}")
            helper.drop(helper[zero_mask].index, inplace=True)
            sample.drop(sample[zero_mask].index, inplace=True)
        sample = sample.div(helper, axis=0)

    if smoothing:
        sample = sample.apply(
            lambda x: savgol_filter(
                x, window_length=smooth_window, polyorder=smooth_polyorder
            ),
            axis=1,
            result_type="broadcast",
        )
             
    if overlay_mode.lower() == "mean":
        sample = pd.DataFrame(sample.mean(numeric_only=True))
    elif overlay_mode.lower() == "median":
        sample = pd.DataFrame(sample.median(numeric_only=True))
    else:
        raise ValueError(f"{overlay_mode} is not a valid keyword.")
   
        
    if rolling:
        if edge_norm:
            trend = sample.rolling(rolling_window, center=True, min_periods=1).median()
            sample = sample.div(trend)
        else:
            flank_start, flank_end = get_window_slice(len(sample),flank*2)
            norm_sample = sample.div(sample.iloc[flank_start:flank_end].mean())
            trend = norm_sample.rolling(rolling_window, center=True, min_periods=1).median()
            sample = sample.div(trend)
        
            
    sample["position"] = calculate_flanking_regions(len(sample))
    sample = sample.set_index("position")
    
    #sample = sample.iloc[fstop:fstart, :]

    return sample


def plot_case_control_overlay(
    cases: pd.DataFrame,
    controls: pd.DataFrame,
    aggregate_controls: bool = False,
    figsize: tuple = (12, 9),
    signal: str = "Coverage",
    target: str = "ROI",
    corrected: bool = False,
    lower_limit:float = None,
    upper_limit: float = None,
):
    ylab = f"Normalized {signal}"

    if corrected:
        signal_correction = "GC corrected"
    else:
        signal_correction = "GC uncorrected"

    title = f"{signal_correction} {signal} at {target}"

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    if aggregate_controls:
        ax.plot(controls["control_mean"], label="Controls", color="black", alpha=0.7)
        ax.fill_between(
            controls["control_mean"].index.astype(int),
            (
                controls["control_mean"].values.flatten()
                - 2 * controls["control_std"].values
            ),
            (
                controls["control_mean"].values.flatten()
                + 2 * controls["control_std"].values
            ),
            alpha=0.2,
            color="black",
        )

    else:
        for col in controls:
            ax.plot(controls[col], label=col, color="black", alpha=0.7)

    for col in cases:
        ax.plot(cases[col], label=col, alpha=0.7)

    ax.set_xlabel("Position relative to target site [bp]", fontsize=14)
    ax.set_ylabel(ylab, fontsize=14)

    if lower_limit or upper_limit:
        shared_ylim = ax.get_ylim()
        if lower_limit:
            shared_ylim = ( min(shared_ylim[0], lower_limit), shared_ylim[1])
        if upper_limit:
            shared_ylim = ( shared_ylim[0], max(shared_ylim[1], upper_limit) )
        
        ax.set_ylim(shared_ylim)

    fig.suptitle(title, fontsize=16)

    if len(ax.get_legend_handles_labels()[1]) < 15:
        fig.legend(bbox_to_anchor=(1, 0.95), loc="upper left")
    else:
        fig.legend(bbox_to_anchor=(1, 0.95), loc="upper left", ncol=2)

    fig.tight_layout()

    return fig


@click.command()
@click.argument("cases", type=click.Path(readable=True), nargs=-1)
@click.option(
    "--control_samples",
    "-c",
    "control_samples",
    required=True,
    type=click.Path(readable=True),
    help="Path to the file containing processed and aggregated control samples in .tsv format. Supports gzipped formats by adding .gz file extension.",
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
    "--case_name",
    "case_name",
    default="Case",
    show_default=True,
    type=click.STRING,
    help="Name of the case samples that will be used for naming.",
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
    "--aggregate_controls",
    "aggregate_controls",
    is_flag=True,
    show_default=True,
    default=False,
    help="Switch between showing controls as separate lines [switch off] or as average with 95% confidence interval [switch on].",
)
@click.option(
    "--GC_corrected",
    "GC_corrected",
    is_flag=True,
    show_default=True,
    default=False,
    help="Switch for showing the correct plot descriptions on whether the data was GC corrected or not.",
)
@click.option(
    "--figsize",
    "figsize",
    type=(int, int),
    default=(12, 9),
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
    cases,
    control_samples,
    output,
    signal,
    target,
    case_name,
    overlay_mode,
    smoothing,
    smooth_window,
    smooth_polyorder,
    rolling,
    rolling_window,
    flank_norm,
    flank,
    display_window,
    aggregate_controls,
    GC_corrected,
    figsize,
    lower_limit,
    upper_limit,
):
    """Takes aggregated control samples and case samples. Processes the cases and plots them together with the controls and saves the plot as image."""

    controls_df = load_controls(control_samples)
    if aggregate_controls:
        control_conf_df = pd.DataFrame()
        control_mean = controls_df.mean(axis=1, numeric_only=True)
        control_std = controls_df.std(axis=1, numeric_only=True)
        control_conf_df["control_mean"] = control_mean
        control_conf_df["control_std"] = control_std
        controls_df = control_conf_df

    cases_df = pd.DataFrame()
    for i, case in enumerate(cases):
        print(f"{case_name}_{i}", case)
        ID = f"{case_name}_{i}"
        tmpdf = load_table(case)
        tmpdf = process_sample(
            tmpdf,
            overlay_mode=overlay_mode,
            smoothing=smoothing,
            smooth_window=smooth_window,
            smooth_polyorder=smooth_polyorder,
            rolling=rolling,
            rolling_window=rolling_window,
            edge_norm=flank_norm,
            flank=flank,
        )
        cases_df[ID] = tmpdf

    if display_window:
        control_min = controls_df.index.min()
        control_max = controls_df.index.max()
        case_min = cases_df.index.min()
        case_max = cases_df.index.max()

        control_lower = max(control_min, display_window[0])
        control_upper = min(control_max, display_window[1])
        case_lower = max(case_min, display_window[0])
        case_upper = min(case_max, display_window[1])

        controls_df = controls_df.loc[control_lower:control_upper]
        cases_df = cases_df.loc[case_lower:case_upper]

    Fig = plot_case_control_overlay(
        controls=controls_df,
        cases=cases_df,
        aggregate_controls=aggregate_controls,
        target=target,
        figsize=figsize,
        signal=signal,
        corrected=GC_corrected,
        lower_limit=lower_limit,
        upper_limit=upper_limit,
    )
    Fig.savefig(output, bbox_inches="tight")


if __name__ == "__main__":
    main()
