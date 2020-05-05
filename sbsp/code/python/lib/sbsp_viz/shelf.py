import logging
import pandas as pd
from typing import *

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