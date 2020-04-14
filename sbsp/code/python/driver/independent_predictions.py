# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20
import logging
import argparse
import numpy as np
import pandas as pd
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_general import Environment
import sbsp_viz.sns as sns
from sbsp_general.shelf import next_name
from sbsp_viz.colormap import ColorMap as CM

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_viz.general import FigureOptions, save_figure

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--max-candidates', required=False, type=int, default=100)
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
logger = logging.getLogger("logger")  # type: logging.Logger

def sensitivity_random(num_candidates, sen_a, sen_b):
    # type: (int, float, float) -> float
    return 1.0 / num_candidates

def sensitivity_independent(num_candidates, sen_a, sen_b):
    # type: (int, float, float) -> float
    err_a = 1-sen_a
    err_b = 1-sen_b

    if sen_a == 0 or sen_b == 0:
        return 0

    try:
        prob = 1.0 / (1.0 + (num_candidates - 1) * ((err_a*err_b)/(sen_a*sen_b)))
    except ValueError:
        prob = 0

    return prob


def sensitivity_fully_dependent(num_candidates, sen_a, sen_b):
    # type: (int, float, float) -> float
    err_a = 1 - sen_a
    err_b = 1 - sen_b

    if sen_a == 0 or sen_b == 0:
        return 0

    try:
        prob = 1.0 / (1.0 + (num_candidates - 1) * (err_b / sen_b))
    except ValueError:
        prob = 0

    return prob

def agree_given_pred_random(num_candidates, sen_a, sen_b):
    # type: (int, float, float) -> float
    return num_candidates / pow(float(num_candidates),  3.0)

def agree_given_pred_independent(num_candidates, sen_a, sen_b):
    # type: (int, float, float) -> float
    if num_candidates == 1:
        return 1.0

    err_a = 1-sen_a
    err_b = 1-sen_b

    return (sen_a * sen_b) / (num_candidates) + (num_candidates-1) * err_b * err_a / (num_candidates * (num_candidates-1)*(num_candidates-1))


def agree_given_pred_fully_dependent(num_candidates, sen_a, sen_b):
    # type: (int, float, float) -> float
    if num_candidates == 1:
        return 1.0

    err_a = 1 - sen_a
    err_b = 1 - sen_b

    return sen_a / float(num_candidates) + err_a / float(num_candidates)


def stack_columns_as_rows(df, stack, new_col, labels, label_col="class"):
    # type: (pd.DataFrame, List[str], str, List[str], str) -> pd.DataFrame

    if labels is not None:
        if len(labels) != len(stack):
            raise ValueError("Labels and stack must have same number of elements: {} != {}".format(
                len(labels), len(stack))
            )
    else:
        labels = stack


    list_df = list()

    for c, l in zip(stack, labels):
        if c not in df.columns:
            logger.warning("Cannot find column {} in data frame.".format(c))
            continue

        df_curr = df.copy()
        df_curr.rename(columns={c: new_col}, inplace=True)            # new column name
        df_curr[label_col] = l

        list_df.append(df_curr)

    return pd.concat(list_df, axis=0, sort=False)


def plot_sensitivities_vs_num_candidates(sensitiviies_func, max_candidates, sen_a, sen_b):
    # type: (Dict[str, Callable], int, float, float) -> None

    list_entries = list()
    for i in range(1, max_candidates+1):
        curr_sensitivities = {
            name: sensitiviies_func[name](i, sen_a, sen_b) for name in sensitiviies_func.keys()
        }
        list_entries.append({
            "Number of candidates": i,
            **curr_sensitivities
        })

    df = pd.DataFrame(list_entries)
    conditions = sorted(list(sensitiviies_func.keys()))
    df_stacked = stack_columns_as_rows(df, conditions, "Probability", conditions, "Condition")

    sns.lineplot(df_stacked, "Number of candidates", "Probability", hue="Condition")


def compute_data(sensitivity_funcs, agree_given_pred_funcs, max_candidates):
    # type: (Dict[str, Callable[[int, float, float], float]], Dict[str, Callable[[int, float, float], float]], int) -> pd.DataFrame

    min_sen = 0.0
    max_sen = 1.0
    step_sen = 0.05

    list_entries = list()
    for c in range(1, max_candidates+1):
        for sen_a in np.arange(min_sen, max_sen+0.001, step_sen):
            for sen_b in np.arange(min_sen, max_sen+0.001, step_sen):
                for condition in sensitivity_funcs.keys():

                    sen_a = round(sen_a, 3)
                    sen_b = round(sen_b, 3)

                    list_entries.append({
                        "Number of candidates": c,
                        "Condition": condition,
                        "Sensitivity A": sen_a,
                        "Sensitivity B": sen_b,
                        "Probability": sensitivity_funcs[condition](c, sen_a, sen_b),
                        "Agree given prediction": agree_given_pred_funcs[condition](c, sen_a, sen_b)
                    })



    df = pd.DataFrame(list_entries)
    return df










def analyze_independent_predictions(max_candidates, sen_a, sen_b):
    # type: (int, float, float) -> None

    sensitivities = {
        "Random": sensitivity_random,
        "Independent": sensitivity_independent,
        "Fully dependent": sensitivity_fully_dependent
    }

    agree_given_pred = {
        "Random": agree_given_pred_random,
        "Independent": agree_given_pred_independent,
        "Fully dependent": agree_given_pred_fully_dependent
    }

    df = compute_data(sensitivities, agree_given_pred, max_candidates)

    plot_sensitivities_vs_num_candidates(sensitivities, max_candidates, sen_a, sen_b)

    sns.lineplot(
        df[(df["Sensitivity A"] == 0.9) & (df["Sensitivity B"] == 0.9)],
        "Number of candidates", "Probability", hue="Condition",
        sns_kwargs={"palette": CM.get_map("independence-conditions")},
        legend_loc="best",
        figure_options=FigureOptions(
            save_fig=next_name("."),
            ylabel=r"$P(y=s|x_1=y, x_2=y)$"

        )
    )

    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, sharey="all", figsize=(10,4))

    sns.lineplot(
        df[(df["Sensitivity A"] == 0.9) & (df["Sensitivity B"] == 0.9)],
        "Number of candidates", "Probability", hue="Condition",
        sns_kwargs={"palette": CM.get_map("independence-conditions")},
        ax=axes[0],
        legend=False,
        figure_options=FigureOptions(
            title="Sensitivity = 0.9",
        )
    )

    sns.lineplot(
        df[(df["Sensitivity A"] == df["Sensitivity B"]) & (df["Number of candidates"] == 25)],
        "Sensitivity A", "Probability", hue="Condition",
        ax=axes[1],
        sns_kwargs={"palette": CM.get_map("independence-conditions")}, figure_options=FigureOptions(
            ylim=[0,1.05],
            xlim=[0,1],
            xlabel="Sensitivity",
            title="Number of candidates = 25",
        )
    )

    save_figure(FigureOptions(save_fig=next_name(".")), fig)
    plt.show()

    df_tmp =  df[
            (df["Sensitivity A"] == df["Sensitivity B"]) &
            (df["Condition"] == "Independent") &
            (df["Sensitivity A"].isin({0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}))
            ]
    df_tmp.rename(columns={"Sensitivity A": "Sensitivity"}, inplace=True)

    sns.lineplot(
       df_tmp,
        "Number of candidates", "Probability", hue="Sensitivity", figure_options=FigureOptions(
            # ylim=[0, 1.05],
            # xlim=[0, 1],
            title="Independent algorithms",
            save_fig=next_name(".")
        ),
    )

    # for condition in set(df["Condition"]):
    #
    #     sns.kdeplot(
    #         df[(df["Condition"] == condition) & (df["Sensitivity A"] == df["Sensitivity B"])],
    #         "Sensitivity A", "Number of candidates", "Probability",
    #         figure_options=FigureOptions(
    #             title=condition
    #         ))
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, sharey="all", figsize=(10, 4))

    sns.lineplot(
        df[(df["Sensitivity A"] == 0.9) & (df["Sensitivity B"] == 0.9)],
        "Number of candidates", "Agree given prediction", hue="Condition",
        sns_kwargs={"palette": CM.get_map("independence-conditions")},
        ax=axes[0],
        legend=False,
        figure_options=FigureOptions(
            title="Sensitivity = 0.9",
        )
    )

    sns.lineplot(
        df[(df["Sensitivity A"] == df["Sensitivity B"]) & (df["Number of candidates"] == 25)],
        "Sensitivity A", "Agree given prediction", hue="Condition",
        ax=axes[1],
        sns_kwargs={"palette": CM.get_map("independence-conditions")}, figure_options=FigureOptions(
            ylim=[0, 1.05],
            xlim=[0, 1],
            xlabel="Sensitivity",
            title="Number of targets = 25",
        )
    )

    save_figure(FigureOptions(save_fig=next_name(".")), fig)
    plt.show()


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    analyze_independent_predictions(args.max_candidates, 0.9, 0.9)


if __name__ == "__main__":
    main(my_env, parsed_args)
