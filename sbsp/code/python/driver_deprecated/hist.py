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
from sbsp_viz.general import FigureOptions
from sbsp_viz.hist import plot_hist_by_group

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #

parser = argparse.ArgumentParser("Plot a histogram of columns in a delimited file.")

parser.add_argument('--pf-data', required=True, type=str, help="Data file")
parser.add_argument('--column-names', required=True, type=str, nargs="+",
                    help="Column name(s) for which to build histograms")
parser.add_argument('--delimiter', required=False, default=",", help="Delimiter in data file")

parser.add_argument('--filter-by-equal', required=False, nargs=2, help="E.g. q-is-true 1")
parser.add_argument('--group-by', required=False, type=str)
parser.add_argument('--title-name', required=False, help="Title name to be added to title")

parser.add_argument('--xlim', nargs="+", type=float)
parser.add_argument('--bins', type=int, required=False, default=20, help="Number of bins")

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


def plot_histograms_for_columns(env, df_data, column_names, **kwargs):
    # type: (Environment, pd.DataFrame, List[str], Dict[str, Any]) -> None

    group_by = get_value(kwargs, "group_by", None)

    title_name = get_value(kwargs, "title_name", "", default_if_none=True)
    xlim = get_value(kwargs, "xlim", None)

    for c in column_names:
        plot_hist_by_group(
            df_data,
            c,
            group_by,
            figure_options=FigureOptions(
                xlabel=c, title="{}".format(title_name), ylabel="Frequency",
                save_fig=os.path.join(
                    env["pd-work"],
                    "hist{}.pdf".format(c.replace(" ", "_").replace("(", "").replace(")", ""))
                ),
                xlim=xlim
            ),
            bins=get_value(kwargs, "bins", 10)
        )


def filter_dataframe_by_equal(df, column_name, value):
    # type: (pd.DataFrame, str, str) -> pd.DataFrame

    if column_name not in df.columns.values:
        raise ValueError("Invalid filter column name ({})".format(column_name))

    column_type = df[column_name].dtype.type

    try:
        value = column_type(value)  # try casting
    except ValueError:
        raise ValueError("Filter value ({}) not of valid type ({})".format(value, column_type))

    return df[df[column_name] == value]


def plot_histograms_for_columns_from_files(env, pf_data, column_names, delimiter=",", **kwargs):
    # type: (Environment, str, list[str], str, **str) -> None

    filter_by_equal = get_value(kwargs, "filter_by_equal", None)

    df = pd.read_csv(pf_data, delimiter=delimiter)

    if filter_by_equal is not None:
        filter_column_name, value = filter_by_equal
        df = filter_dataframe_by_equal(df, filter_column_name, value)

    plot_histograms_for_columns(env, df, column_names, **kwargs)


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    plot_histograms_for_columns_from_files(env, args.pf_data, args.column_names,
                                           delimiter=args.delimiter,
                                           filter_by_equal=args.filter_by_equal,
                                           group_by=args.group_by,
                                           title_name=args.title_name,
                                           bins=args.bins,
                                           xlim=args.xlim)


if __name__ == "__main__":
    main(my_env, parsed_args)

