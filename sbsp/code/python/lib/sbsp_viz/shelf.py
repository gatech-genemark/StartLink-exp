import logging
import math

import matplotlib.cm
import matplotlib.colors
import pandas as pd
from typing import *

import seaborn
from statsmodels import api as sm

from sbsp_general import Environment
from sbsp_general.general import get_value
from sbsp_viz.general import FigureOptions

log = logging.getLogger(__name__)


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
            log.warning("Cannot find column {} in data frame.".format(c))
            continue

        df_curr = df.copy()
        df_curr.rename(columns={c: new_col}, inplace=True)            # new column name
        df_curr[label_col] = l

        df_curr.drop(set(labels).difference({l}), inplace=True, axis=1)

        list_df.append(df_curr)

    return pd.concat(list_df, axis=0, sort=True)


def get_loess(local_x, local_y):
    # l = loess(local_x, local_y)
    # l.fit()
    lowess = sm.nonparametric.lowess(local_y, local_x)
    return lowess[:,1]
    pred = l.predict(local_x, stderror=False)
    return pred.values

import matplotlib.pyplot as plt
import numpy as np; np.random.seed(10)
import seaborn.distributions as sd
from seaborn.palettes import color_palette, blend_palette
from six import string_types

def _bivariate_kdeplot(x, y, filled, fill_lowest,
                       kernel, bw, gridsize, cut, clip,
                       axlabel, cbar, cbar_ax, cbar_kws, ax, **kwargs):
    """Plot a joint KDE estimate as a bivariate contour plot."""
    # Determine the clipping
    if clip is None:
        clip = [(-np.inf, np.inf), (-np.inf, np.inf)]
    elif np.ndim(clip) == 1:
        clip = [clip, clip]

    # Calculate the KDE
    if sd._has_statsmodels:
        xx, yy, z = sd._statsmodels_bivariate_kde(x, y, bw, gridsize, cut, clip)
    else:
        xx, yy, z = sd._scipy_bivariate_kde(x, y, bw, gridsize, cut, clip)

    # Plot the contours
    n_levels = kwargs.pop("n_levels", 10)
    cmap = kwargs.get("cmap", "BuGn" if filled else "BuGn_d")
    if isinstance(cmap, string_types):
        if cmap.endswith("_d"):
            pal = ["#333333"]
            pal.extend(color_palette(cmap.replace("_d", "_r"), 2))
            cmap = blend_palette(pal, as_cmap=True)
        else:
            cmap = plt.cm.get_cmap(cmap)

    kwargs["cmap"] = cmap
    contour_func = ax.contourf if filled else ax.contour
    cset = contour_func(xx, yy, z, n_levels, **kwargs)
    if filled and not fill_lowest:
        cset.collections[0].set_alpha(0)
    kwargs["n_levels"] = n_levels

    if cbar:
        cbar_kws = {} if cbar_kws is None else cbar_kws
        ax.figure.colorbar(cset, cbar_ax, ax, **cbar_kws)

    # Label the axes
    if hasattr(x, "name") and axlabel:
        ax.set_xlabel(x.name)
    if hasattr(y, "name") and axlabel:
        ax.set_ylabel(y.name)

    return ax, cset


def loess_with_stde(df, xcol, ycol, ax, label, **kwargs):

    xlim = get_value(kwargs, "xlim", None)
    ylim = get_value(kwargs, "ylim", None)
    x = df[xcol].values

    df.set_index(xcol, inplace=True, drop=False)
    w = 30

    y = df[ycol].values #df[ycol].rolling(window=w, min_periods=1).mean().values
    std = df[ycol].rolling(window=w, min_periods=1).std().values
    std[0] = 0


    y = get_loess(x, y)
    std = get_loess(x, std)
    y_u = y + std
    y_l = y - std
    import numpy as np
    # ax.plot(x, df[ycol], ".", alpha=0.1)

    # seaborn.kdeplot(df[xcol], df[ycol], cmap="Reds", ax=ax,
    #                 **kwargs)
    # heatmap1_data = pd.pivot_table(df, values=xcol,
    #                                index=['continent'],
                                   # columns=ycol)
    # seaborn.heatmap(heatmap1_data, ax=ax)
    heatmap_grid_data_single(None, df, xcol, ycol, ax=ax, figure_options=None,  **kwargs)

    ax.set_xlabel(None)
    ax.set_ylabel(None)

    ax2 = ax.twinx().twiny()
    ax2.plot(x, y, label=label, color="blue")
    if xlim is not None:
        ax2.set_xlim(*xlim)
    if ylim is not None:
        ax2.set_ylim(*ylim)
    ax2.set_xticks([])
    ax2.set_yticks([])
    # ax2.axes.xaxis.set_visible(False)
    # ax2.axes.yaxis.set_visible(False)

    # ax.fill_between(x, y_l, y_u,  alpha=0.33, color="orange", zorder=4)

    return x, y, y_l, y_u


def create_mappable_for_colorbar(values, cmap):
    import matplotlib.colors
    import matplotlib.cm
    vmin=min(values)
    vmax = max(values)

    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    mappable = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array(values)
    return mappable


def heatmap_grid_data_single(env, df_group, x, y, figure_options, ax, **kwargs):
    # type: (Environment, pd.DataFrame, str, str, FigureOptions, plt.Axes, Dict[str, Any]) -> None
    import matplotlib.pyplot as plt

    xlabel = get_value(kwargs, "xlabel", None)
    ylabel = get_value(kwargs, "ylabel", None)
    hue = get_value(kwargs, "hue", None)
    num_steps = get_value(kwargs, "num_steps", 20)
    cbar_max = get_value(kwargs, "cbar_max", None)
    cbar_ax = get_value(kwargs, "cbar_ax", None)
    cbar = get_value(kwargs, "cbar", False)

    xlim = get_value(kwargs, "xlim", None)
    ylim = get_value(kwargs, "ylim", None)
    balance = get_value(kwargs, "balance", False)

    if xlabel is None:
        xlabel = x
    if ylabel is None:
        ylabel = y



    axis_idx = 0
    curr_row = 0
    curr_col = 0
    mappable = None



    axis_idx += 1
    curr_col += 1

    min_x = min(df_group[x])
    max_x = max(df_group[x]) + 0.000000001

    min_y = min(df_group[y])
    max_y = max(df_group[y]) + 0.000000001

    if xlim is not None:
        min_x, max_x = xlim
    if ylim is not None:
        min_y, max_y = ylim

    if balance:
        min_x = min_y = min(min_x, min_y)
        max_x = max_y = max(max_x, max_y)

    ss_x = (max_x - min_x) / float(num_steps)
    ss_y = (max_y - min_y) / float(num_steps)

    num_col = num_steps
    num_row = num_steps
    import numpy as np
    gms2_eq_sbsp_and_ncbi = np.zeros([num_row, num_col], dtype=float)
    gms2_eq_sbsp_eq_ncbi = np.zeros([num_row, num_col], dtype=float)

    total = 0
    df_group.reset_index(inplace=True, drop=True)
    for index in df_group.index:
        x_val = df_group.at[index, x]
        y_val = df_group.at[index, y]

        x_pos = int((x_val - min_x) / ss_x)
        y_pos = int((y_val - min_y) / ss_y)

        y_pos = min(y_pos, gms2_eq_sbsp_eq_ncbi.shape[1]-1)

        gms2_eq_sbsp_and_ncbi[x_pos][y_pos] += 1
        total += 1

    accuracy = 100 * np.divide(gms2_eq_sbsp_and_ncbi, total)

    import seaborn
    import matplotlib.pyplot as plt

    xticks = [int(z) for z in list(np.arange(0, num_steps+0.001, int(num_steps / 5)))]
    yticks = [int(z) for z in list(np.arange(0, num_steps+0.001, int(num_steps / 5)))]

    l_x = np.arange(min_x, max_x+0.001, ss_x)
    l_y = np.arange(min_y, max_y+0.001, ss_y)
    xticklabels = [round(l_x[i], 2) for i in xticks]
    yticklabels = [round(l_y[i], 2) for i in yticks]
    g = seaborn.heatmap(
        accuracy.transpose(), vmin=0, vmax=np.amax(accuracy) if cbar_max is None else cbar_max,
        xticklabels=xticklabels, yticklabels=yticklabels, ax=ax,
        # cmap=seaborn.light_palette("green"),
        cmap="Reds",
        cbar=cbar,
        cbar_ax=cbar_ax
    )
    # cbar_ax=None if axis_idx != 0 else cbar_ax, cbar=axis_idx==0)

    # cbar=g.cbar

    g.invert_yaxis()
    g.set_xticks(xticks)
    g.set_yticks(yticks)
    g.set_xticklabels(xticklabels, rotation=0)



def heatmap_grid_data(env, df, x, y, figure_options, **kwargs):
    # type: (Environment, pd.DataFrame, str, str, FigureOptions, Dict[str, Any]) -> None
    import matplotlib.pyplot as plt

    xlabel = get_value(kwargs, "xlabel", None)
    ylabel = get_value(kwargs, "ylabel", None)
    hue = get_value(kwargs, "hue", None)
    num_steps = get_value(kwargs, "num_steps", 20)

    xlim = get_value(kwargs, "xlim", None)
    ylim = get_value(kwargs, "ylim", None)
    balance = get_value(kwargs, "balance", False)

    if xlabel is None:
        xlabel = x
    if ylabel is None:
        ylabel = y

    hue_values = sorted(list(set(df[hue])))
    fig, axes = plt.subplots(2, math.ceil(len(hue_values) / 2), sharex="all", sharey="all")
    cbar_ax = fig.add_axes([.91, .3, .03, .4])

    # fig = plt.figure()
    num_rows = 2
    num_cols = math.ceil(len(hue_values) / 2)


    axis_idx = 0
    curr_row = 0
    curr_col = 0
    mappable = None
    for ancestor, df_group in df.groupby(hue, as_index=False):
        ax = axes.ravel()[axis_idx]

        axis_idx += 1
        curr_col += 1
        if curr_col == math.ceil(len(hue_values) / 2):
            curr_row += 1
            curr_col = 0

        min_x = min(df_group[x])
        max_x = max(df_group[x]) + 0.000000001

        min_y = min(df_group[y])
        max_y = max(df_group[y]) + 0.000000001

        if balance:
            min_x = min_y = min(min_x, min_y)
            max_x = max_y = max(max_x, max_y)

        ss_x = (max_x - min_x) / float(num_steps)
        ss_y = (max_y - min_y) / float(num_steps)

        num_col = num_steps
        num_row = num_steps
        import numpy as np
        gms2_eq_sbsp_and_ncbi = np.zeros([num_row, num_col], dtype=float)
        gms2_eq_sbsp_eq_ncbi = np.zeros([num_row, num_col], dtype=float)

        total = 0
        for index in df_group.index:

            x_val = df_group.at[index, x]
            y_val = df_group.at[index, y]

            x_pos = int((x_val-min_x) / ss_x)
            y_pos = int((y_val-min_y) / ss_y)

            gms2_eq_sbsp_and_ncbi[x_pos][y_pos] += 1
            total += 1

        accuracy = 100 * np.divide(gms2_eq_sbsp_eq_ncbi, total)

        import seaborn
        import matplotlib.pyplot as plt

        xticks = list(range(0, num_steps, int(num_steps/5)))
        yticks = list(range(0, num_steps, int(num_steps/5)))

        l_x = np.arange(min_x, max_x, ss_x)
        l_y = np.arange(min_y, max_y, ss_y)
        xticklabels = [round(l_x[i], 2) for i in xticks]
        yticklabels = [round(l_y[i], 2) for i in yticks]
        g = seaborn.heatmap(
            accuracy.transpose(), vmin=0, vmax=100,
            xticklabels=xticklabels, yticklabels=yticklabels, ax=ax,
            # cmap=seaborn.light_palette("green"),
            cmap="Blues",
        cbar=False)
        # cbar_ax=None if axis_idx != 0 else cbar_ax, cbar=axis_idx==0)

        # cbar=g.cbar

        g.invert_yaxis()
        g.set_xticks(xticks)
        g.set_yticks(yticks)
        g.set_xticklabels(xticklabels, rotation=0)

        g.set_title(ancestor)
        mappable = ax.collections[0]

    cbar_ax = fig.axes[-1]


    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(top=False, bottom=False, left=False, right=False, which="both",
                    labelbottom=False, labeltop=False, labelleft=False, labelright=False)
    plt.xlabel(xlabel, labelpad=20)
    plt.ylabel(ylabel, labelpad=30)

    # ax3 = plt.subplot2grid((num_rows, num_cols), (0, num_cols - 1), rowspan=num_rows,
    #                        )
    plt.colorbar(mappable, cax=cbar_ax)
    fig.tight_layout(rect=[0, 0, .9, 1])

    plt.show()
