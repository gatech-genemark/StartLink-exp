from typing import *
import pandas as pd
from sbsp_general.general import get_value
from sbsp_viz.general import FigureOptions
import matplotlib.pyplot as plt
import seaborn as sns


def plot_hist_by_group(df_data, column_x, column_group=None, figure_options=None, **kwargs):
    # type: (pd.DataFrame, str, Union[str, None], FigureOptions, Dict[str, Any]) -> None

    bins = get_value(kwargs, "bins", 10)

    _, ax = plt.subplots()

    cumulative = get_value(kwargs, "cumulative", False)
    shade = False if cumulative else True
    cut = [min(df_data[column_x]), max(df_data[column_x])]
    if column_group is not None:
        for name, df_group in df_data.groupby(column_group):
            sns.distplot(df_group[column_x], hist=False, kde_kws={"shade": shade, "cumulative": cumulative}, label=name)
    else:
        # sns.distplot(df_data[column_x], hist=True, kde_kws={"shade": shade, "cumulative": cumulative, "clip": cut})
        sns.distplot(df_data[column_x], bins=bins, hist=True, kde=False, hist_kws={"edgecolor": "black"})

    FigureOptions.set_properties_for_axis(ax, figure_options)

    # plt.xlim([min(df_data[column_x]), max(df_data[column_x])])
    if figure_options is not None and figure_options.save_fig is not None:
        plt.savefig(figure_options.save_fig, bbox_index="tight")

    plt.show()

def plot_catplot(df, column_x, column_y, figure_options=None):
    _, ax = plt.subplots()
    sns.catplot(x=column_x, y=column_y, kind="bar", data=df)

    FigureOptions.set_properties_for_axis(ax, figure_options)
    if figure_options is not None and figure_options.save_fig is not None:
        plt.savefig(figure_options.save_fig, bbox_index="tight")

    plt.show()
