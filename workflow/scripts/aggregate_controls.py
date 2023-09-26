#!/usr/bin/env python

import click
import numpy as np
import pandas as pd
from scipy import stats
from scipy.signal import savgol_filter
from loguru import logger


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
        edge_norm (bool, optional): Switch on whether the signal should be demeaned by substracting the average signal in the flanking regions. Defaults to False.
        window (int, optional): Size of window centered on region of interest. Defaults to 1000.
        flank (int, optional): Size of flanking regions for edge_norm. Defaults to 500.

    Raises:
        ValueError: Invalid overlay_mode keyword.

    Returns:
        pd.DataFrame: Processed sample.
    """

    # print(locals())
    if window % 2 == 0:
        fstart = int(window / 2 + 1)
        fstop = int(-window / 2)
    elif window % 2 == 1:
        fstart = int(window / 2 - 0.5 + 1)
        fstop = int(-window / 2 + 0.5)

    if edge_norm:
        flank_start, flank_end = get_window_slice(len(sample.columns), flank * 2)

        helper = abs(sample.iloc[:, flank_start:flank_end].mean(axis=1))
        if (helper == 0).any():
            zero_mask = helper == 0
            print(
                f"Regions with zero mean in Normalization slice [{flank_start}, {flank_end}] encountered. Removing respective regions."
            )
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
            flank_start, flank_end = get_window_slice(len(sample), flank * 2)
            norm_sample = sample.div(sample.iloc[flank_start:flank_end].mean())
            trend = norm_sample.rolling(
                rolling_window, center=True, min_periods=1
            ).median()
            sample = sample.div(trend)

    sample["position"] = calculate_flanking_regions(len(sample))
    sample = sample.set_index("position")

    # sample = sample.iloc[fstop:fstart, :]

    return sample


@click.command()
@click.argument("samples", type=click.Path(readable=True), nargs=-1)
@click.option(
    "--output",
    "-o",
    "output",
    required=True,
    type=click.Path(writable=True),
    help="Path to the outputfile in .tsv format. Supports gzipped formats by adding .gz file extension.",
)
@click.option(
    "--control_name",
    "-c",
    "control_name",
    default="control",
    show_default=True,
    type=click.STRING,
    help="Name of the control samples (e.g. control or healthy). Will be used for naming aggregated data.",
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
def main(
    samples,
    output,
    control_name,
    overlay_mode,
    smoothing,
    smooth_window,
    smooth_polyorder,
    rolling,
    rolling_window,
    flank_norm,
    flank,
):
    """
    Takes an arbitrary number of control signal files (.csv/.tsv or .csv.gz/.tsv.gz) with rows as individual sites and columns as positions.
    Signals are loaded, processed and aggregated per sample and exported as a combined file (.tsv or .tsv.gz).
    """

    logger.info("Parameters")
    logger.info(locals())

    logger.info("Starting aggregation...")
    control_df = pd.DataFrame()
    for i, control in enumerate(samples):
        logger.info(f"Loading {control_name}_{i}: {control}")
        ID = f"{control_name}_{i}"
        tmpdf = load_table(control)
        logger.info(f"Processing {control_name}_{i}: {control}")
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
        control_df[ID] = tmpdf
        del tmpdf
    logger.info(f"Saving aggregated results to: {output}")
    control_df.to_csv(output, sep="\t")
    logger.info("Done.")


if __name__ == "__main__":
    main()
