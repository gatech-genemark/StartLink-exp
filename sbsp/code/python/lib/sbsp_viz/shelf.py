import logging
import pandas as pd
from typing import *

import seaborn
from statsmodels import api as sm

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


def loess_with_stde(df, xcol, ycol, ax, label):

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

    # ax.plot(x, df[ycol], ".", alpha=0.1)
    seaborn.kdeplot(df[xcol], df[ycol], cmap="Reds", ax=ax)
    ax.set_xlabel(None)
    ax.set_ylabel(None)

    ax.plot(x, y, label=label, color="blue")
    # ax.fill_between(x, y_l, y_u,  alpha=0.33, color="orange", zorder=4)

    return x, y, y_l, y_u