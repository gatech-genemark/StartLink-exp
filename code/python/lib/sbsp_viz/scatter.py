import logging
import pandas as pd

from typing import *
import matplotlib.pyplot as plt
import seaborn as sns

from sbsp_general.general import get_value
from sbsp_viz.general import FigureOptions, add_identity

logger = logging.getLogger(__name__)


def scatter(df, column_x, column_y, figure_options=None, **kwargs):
    # type: (pd.DataFrame, str, str, FigureOptions, Dict[str, Any]) -> None

    column_z = get_value(kwargs, "column_z", None)
    identity = get_value(kwargs, "identity", False)
    hue = df[column_z] if column_z is not None else None

    _, ax = plt.subplots()


    sns.jointplot(df[column_x], df[column_y], kind="scatter", alpha=0.3, s=10, linewidth=0)
    #sns.scatterplot(df[column_x], df[column_y], hue=hue, alpha=0.3, s=10, linewidth=0)
    if identity:
        add_identity(ax, color="r", ls="--")

    FigureOptions.set_properties_for_axis(ax, figure_options)
    if figure_options is not None and figure_options.save_fig is not None:
        plt.savefig(figure_options.save_fig, bbox_index="tight")

    plt.show()
