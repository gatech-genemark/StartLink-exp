# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 18/05/2019

import logging
import argparse
import os

import pandas as pd
from typing import *

# noinspection PyUnresolvedReferences
import pathmagic  # add path to custom library

# Custom library imports
import sbsp_general
from sbsp_general import Environment
from sbsp_general.general import get_value
from sbsp_viz.general import FigureOptions, plot_scatter_for_dataframe_columns, plot_scatter_matrix, \
    plot_scatter_matrix_for_dataframe_columns

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #

parser = argparse.ArgumentParser("Plot a scatter plots of columns in a delimited file.")

parser.add_argument('--pf-data', required=True, type=str, help="Data file")

parser.add_argument('--delimiter', required=False, default=",", help="Delimiter in data file")

column_input_group = parser.add_mutually_exclusive_group(required=True)
column_input_group.add_argument('--column-names', type=str, nargs="+",
                    help="Column name(s) for which to construct pairwise scatter plots")
column_input_group.add_argument('--pf-column-pairs',
                                help="File containing (comma separated) pairs of columns (one pair per line)")


parser.add_argument('--filter-by-equal', required=False, nargs=2, help="E.g. q-is-true 1")
parser.add_argument('--plot-type', required=False, default="separate", choices=["separate", "matrix"])
parser.add_argument('--scatter-in-separate-files', required=False, default=False, action="store_true")
parser.add_argument('--color-by-value', required=False,  help="E.g. q-is-true")
parser.add_argument('--limit-x-axis-features', required=False, default=None, type=str, nargs="+",
                    help="If set, use those features as the x-axis, and all rest as the y-axis. "
                         "Note, this only works in scatter-in-separate-files mode")

parser.add_argument('--jitter', default=False, action="store_true")
parser.add_argument('--title-name', required=False, help="Title name to be added to title")
parser.add_argument('--filter-in-range', required=False, nargs="+")

parser.add_argument('--pd-work', required=False, default=None, help="Path to working directory")
parser.add_argument('--pd-data', required=False, default=None, help="Path to data directory")
parser.add_argument('--pd-results', required=False, default=None, help="Path to results directory")
parser.add_argument("-l", "--log", dest="loglevel", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                    help="Set the logging level", default='WARNING')

parsed_args = parser.parse_args()

# ------------------------------ #
#           Main Code            #
# ------------------------------ #

# Load environment variables
my_env = Environment(pd_data=parsed_args.pd_data,
                     pd_work=parsed_args.pd_work,
                     pd_results=parsed_args.pd_results)

# Setup logger
logging.basicConfig(level=parsed_args.loglevel)
logger = logging.getLogger("logger")                    # type: logging.Logger


def filter_dataframe_by_equal(df, column_name, value):
    # type: (pd.DataFrame, str, str) -> pd.DataFrame

    if column_name not in df.columns.values:
        raise ValueError("Invalid filter column name ({})".format(column_name))

    column_type = df[column_name].dtype.type

    try:
        value = column_type(value)          # try casting
    except ValueError:
        raise ValueError("Filter value ({}) not of valid type ({})".format(value, column_type))

    return df[df[column_name] == value]


def plot_scatter_for_columns_from_files(env, pf_data, column_names, delimiter=",", **kwargs):
    # type: (Environment, str, list[str], str, **str) -> None

    filter_by_equal = get_value(kwargs, "filter_by_equal", None)
    scatter_separately = get_value(kwargs, "scatter_in_separate_files", False)
    limit_x_axis_features = get_value(kwargs, "limit_x_axis_features", None)
    color_by_value = get_value(kwargs, "color_by_value", None)

    title = get_value(kwargs, "title", None)
    df = pd.read_csv(pf_data, delimiter=delimiter)

    if filter_by_equal is not None:
        filter_column_name, value = filter_by_equal
        df = filter_dataframe_by_equal(df, filter_column_name, value)

    if scatter_separately:

        x_axis_column_names = column_names
        if limit_x_axis_features is not None:
            x_axis_column_names = limit_x_axis_features

        for f1 in x_axis_column_names:
            for f2 in column_names:
                plot_scatter_for_dataframe_columns(df, [f1, f2], color_by_value=color_by_value, figure_options=FigureOptions(
                    title=title,
                    save_fig=os.path.join(env["pd-work-results"], "scatter_{}_{}".format(f1, f2))
                ))
    else:
        if color_by_value is not None:
            plot_scatter_matrix(df, column_names, color_by=color_by_value, figure_options=FigureOptions(
                save_fig=os.path.join(env["pd-work-results"], "scatter.pdf")
            ))
        else:
            plot_scatter_matrix_for_dataframe_columns(df, column_names, figure_options=FigureOptions(
                save_fig=os.path.join(env["pd-work-results"], "scatter.pdf")
            ))


def df_plot_scatter_separate(env, df, column_pairs, **kwargs):
    # type: (Environment, pd.DataFrame, List[List], Dict[str, Any]) -> None

    color_by_value = get_value(kwargs, "color_by_value", None)
    limit_x_axis_features = get_value(kwargs, "limit_x_axis_features", None)
    jitter = get_value(kwargs, "jitter", None)
    title = get_value(kwargs, "title", None)


    if limit_x_axis_features is not None:
        column_pairs = [x for x in column_pairs if x[0] in limit_x_axis_features]

    for f1, f2 in column_pairs:

        plot_scatter_for_dataframe_columns(df, [f1, f2], color_by_value=color_by_value, figure_options=FigureOptions(
            save_fig=os.path.join(env["pd-work-results"], "scatter_{}_{}.pdf".format(f1, f2)),
            xlabel=f1,
            ylabel=f2,
            title=title
        ))


def df_plot_scatter_matrix(env, df, column_names, **kwargs):
    # type: (Environment, pd.DataFrame, Union[List, Set], Dict[str, Any]) -> None

    color_by_value = get_value(kwargs, "color_by_value", None)

    if color_by_value is not None:
        plot_scatter_matrix(df, column_names, color_by=color_by_value, figure_options=FigureOptions(
            save_fig=os.path.join(env["pd-work-results"], "scatter.pdf")
        ), **kwargs)
    else:

        plot_scatter_matrix(df, column_names, color_by=color_by_value, figure_options=FigureOptions(
            save_fig=os.path.join(env["pd-work-results"], "scatter.pdf")
        ), **kwargs)
        #plot_scatter_matrix_for_dataframe_columns(df, column_names, figure_options=FigureOptions(
        #    save_fig=os.path.join(env["pd-work-results"], "scatter")
        #))


def read_pairs_from_file(pf_pairs, delimiter=","):
    # type: (str, str) -> List[List]

    with open(pf_pairs, "r") as f:
        pairs = list()
        for line in f:
            items = line.strip().split(delimiter)

            pairs.append([x for x in items])
        return pairs


def all_combinations(items):
    # type: (List) -> List[List]

    pairs = list()

    for a in items:
        for b in items:
            if a == b:
                continue
            pairs.append([a, b])
    return pairs


def df_plot_scatter(env, df, **kwargs):
    # type: (Environment, pd.DataFrame, Dict[str, Any]) -> None

    # Steps:
    # 1) Get plot information
    #   - Type of plot (matrix versus individual scatters)
    #   - Columns or pairs of columns (based on type of plot)
    # 2) Plot

    plot_type = get_value(kwargs, "plot_type", "separate", valid={"separate", "matrix"})
    filter_by_equal = get_value(kwargs, "filter_by_equal", None)
    pf_column_pairs = get_value(kwargs, "pf_column_pairs", None)
    column_names = get_value(kwargs, "column_names", None)

    if filter_by_equal is not None:
        filter_column_name, value = filter_by_equal
        df = filter_dataframe_by_equal(df, filter_column_name, value)

    if plot_type == "separate":

        column_pairs = None
        if column_names:
            column_pairs = all_combinations(column_names)

        if pf_column_pairs:
            column_pairs = read_pairs_from_file(pf_column_pairs)

        df_plot_scatter_separate(env, df, column_pairs=column_pairs, **kwargs)
    else:

        df_plot_scatter_matrix(env, df, **kwargs)


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    # check custom args
    if args.limit_x_axis_features is not None:
        if not args.scatter_in_separate_files:
            raise argparse.ArgumentError("Scatter-in-separate-files should be enabled for --limit-x-axis-features")

    df = pd.read_csv(args.pf_data, header=0)

    if args.column_names is not None and    "aa-mismatch-fraction" in args.column_names:
        df["aa-mismatch-fraction"] = 1 - df["aa-match-fraction"]

    if args.filter_in_range:
        if len(args.filter_in_range) % 3 != 0:
            raise ValueError("Filter-in-range should be a multiple of 3")

        num_filters = len(args.filter_in_range) / 3
        for i in range(int(num_filters)):
            filter_info = args.filter_in_range[i * 3: i * 3 + 3]
            name = filter_info[0]
            min_val = float(filter_info[1])
            max_val = float(filter_info[2])

            df = df[(df[name] > min_val) & (df[name] < max_val)]

    df_plot_scatter(env,
                    df,
                    column_names=args.column_names,
                    plot_type=args.plot_type,
                    filter_by_equal=args.filter_by_equal,
                    limit_x_axis_features=args.limit_x_axis_features,
                    color_by_value=args.color_by_value,
                    pf_column_pairs=args.pf_column_pairs,
                    jitter=args.jitter,
                    title=args.title_name
                    )



if __name__ == "__main__":
    main(my_env, parsed_args)
