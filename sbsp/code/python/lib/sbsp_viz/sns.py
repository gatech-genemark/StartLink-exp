import logging
from typing import *
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


from sbsp_general.general import get_value
from sbsp_viz.general import FigureOptions, add_identity, save_figure

logger = logging.getLogger(__name__)


def scatterplot(df, x, y, hue=None, figure_options=None, **kwargs):
    # type: (pd.DataFrame, str, str, Union[str, None], FigureOptions, Dict[str, Any]) -> None

    sns_kwargs = get_value(kwargs, "sns_kwargs", dict())

    identity = get_value(kwargs, "identity", False)

    _, ax = plt.subplots()

    g = sns.scatterplot(x=x, y=y, hue=hue, data=df, linewidth=0, **sns_kwargs)

    if identity:
        add_identity(ax, color="r", ls="--")

    if hue is not None:
        plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))

    FigureOptions.set_properties_for_axis(ax, figure_options)
    save_figure(figure_options)
    plt.show()


def lineplot(df, x, y, hue=None, figure_options=None, **kwargs):
    # type: (pd.DataFrame, str, str, Union[str, None], FigureOptions, Dict[str, Any]) -> None

    sns_kwargs = get_value(kwargs, "sns_kwargs", dict())

    identity = get_value(kwargs, "identity", False)

    _, ax = plt.subplots()

    g = sns.lineplot(x=x, y=y, hue=hue, data=df, ax=ax,  **sns_kwargs)

    if identity:
        add_identity(ax, color="r", ls="--")

    if hue is not None:
        plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))

    FigureOptions.set_properties_for_axis(ax, figure_options)
    save_figure(figure_options)
    plt.show()



def catplot(df, x, y, hue=None, kind=None, figure_options=None, **kwargs):
    # type: (pd.DataFrame, str, str, str, Union[str, None], FigureOptions, Dict[str, Any]) -> None
    sns_kwargs = get_value(kwargs, "sns_kwargs", dict())

    g = sns.catplot(x=x, y=y, data=df, kind=kind, hue=hue, legend=False, aspect=1.5, **sns_kwargs)

    if hue is not None:
        plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))

    if kind == "point":
        plt.setp(g.ax.lines, linewidth=1.5)  # set lw for all lines of g axes
        # plt.setp(g.ax.lines, markersize=0)  # set lw for all lines of g axes

    FigureOptions.set_properties_for_axis(g.axes[0][0], figure_options)
    save_figure(figure_options)
    plt.show()


def lmplot(df, x, y, hue=None, figure_options=None, **kwargs):
    # type: (pd.DataFrame, str, str, Union[str, None], FigureOptions, Dict[str, Any]) -> None

    sns_kwargs = get_value(kwargs, "sns_kwargs", dict())

    g = sns.lmplot(x=x, y=y, hue=hue, data=df, aspect=2, legend=False, **sns_kwargs)

    if hue is not None:
        g.axes[0][0].legend(loc='center left', bbox_to_anchor=(1.05, 0.5))

    FigureOptions.set_properties_for_axis(g.axes[0][0], figure_options)
    save_figure(figure_options, fig=g.fig)
    plt.subplots_adjust(right=1)
    plt.show()


def distplot(df, x, figure_options=None):
    _, ax = plt.subplots()

    g = sns.distplot(df[x], bins=50, kde=False)

    FigureOptions.set_properties_for_axis(g.axes, figure_options)
    save_figure(figure_options)
    plt.show()


def jointplot(df, x, y, hue=None, figure_options=None, **kwargs):
    _, ax = plt.subplots()
    sns_kwargs = get_value(kwargs, "sns_kwargs", dict())
    # g = sns.lmplot(x=x, y=y, hue=hue, data=df, aspect=2, legend=False, ci=None)
    g = sns.jointplot(x, y, data=df, ax=ax, **sns_kwargs)

    if hue is not None:
        plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))

    # FigureOptions.set_properties_for_axis(g.axes[0][0], figure_options)
    save_figure(figure_options)
    plt.show()


def tsplot(df, x, y, hue=None, figure_options=None):
    _, ax = plt.subplots()

    # g = sns.lmplot(x=x, y=y, hue=hue, data=df, aspect=2, legend=False, ci=None)
    sns.tsplot(df[y].values, df[x].values)

    if hue is not None:
        plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))

    # FigureOptions.set_properties_for_axis(g.axes[0][0], figure_options)
    save_figure(figure_options)
    plt.show()