import logging
import os
import pandas as pd
from typing import *

from sbsp_general.labels_comparison_detailed import LabelsComparisonDetailed
from sbsp_io.general import mkdir_p

logger = logging.getLogger(__name__)


class LabelsComparisonDetailedViz:

    def __init__(self, lcd, **kwargs):
        # type: (LabelsComparisonDetailed, Dict[str, Any]) -> None

        self.lcd = lcd

    def run(self, pd_output):
        # type: (str) -> None

        comparison = self.lcd.comparison

        self._run_helper(comparison["all"], os.path.join(pd_output, "all"))

        for attr_name in comparison["attribute"]:
            self._run_helper_for_attribute(comparison["attribute"][attr_name],
                                           os.path.join(pd_output, attr_name))

    def _run_helper(self, comparison, pd_output):
        # type: (Dict[str, Any], str) -> None

        mkdir_p(pd_output)
        df = self._stats_summary_to_df(comparison["stats"])
        df.to_csv(os.path.join(pd_output, "summary_stats.csv"), index=False)

    def _run_helper_for_attribute(self, value_to_comparison, pd_output):
        # type: (Dict[Any, Dict[str, Any]], str) -> None

        mkdir_p(pd_output)

        list_df = list()
        for value, comparison in value_to_comparison.items():
            df = self._stats_summary_to_df(comparison["stats"])
            list_df.append((value, df))

        df_numbers = self._merge_multiple_stats_summary(list_df, ["Common 3'", "Common 5'"])
        df_percentages = self._merge_multiple_stats_summary(list_df, ["% Common 3'", "% Common 5'"])

        df_numbers.to_csv(os.path.join(pd_output, "numbers.csv"))
        df_percentages.to_csv(os.path.join(pd_output, "percentages.csv"))

    def _merge_multiple_stats_summary(self, list_df, column_names):
        # type: (List[Tuple[str, pd.DataFrame]], List[str]) -> pd.DataFrame

        total_a_tag = "Total {}".format(self.lcd.name_a)

        df = pd.DataFrame({
            "Species": list_df[0][1]["Species"],
            total_a_tag: list_df[0][1][total_a_tag]},
            index=[0]
        )

        for cn in column_names:
            for item in list_df:
                value, df = item
                df["{} - {}".format(cn, value)] = df[cn]

        return df


    def _stats_summary_to_df(self, stats):
        # type: (Dict[str, Any]) -> pd.DataFrame

        total_a_tag = "Total {}".format(self.lcd.name_a)
        total_b_tag = "Total {}".format(self.lcd.name_b)

        dict_summary = {
            "Species": self.lcd.tag,
            total_a_tag: stats["total-a"],
            total_b_tag: stats["total-b"],
            "Common 3'": stats["match-3p"],
            "Common 5'": stats["match-3p-5p"],
        }

        df = pd.DataFrame(dict_summary, index=[0])

        df["% Common 3'"] = 0
        if df[total_a_tag] > 0:
            df["% Common 3'"] = df["Common 3'"] / df[total_a_tag]

        df["Common 5'"] = 0
        if df["Common 3'"] > 0:
            df["% Common 5'"] = df["Common 5'"] / df["Common 3'"]

        df["% Common 3'"].map('${:,.2f}'.format)
        df["% Common 5'"].map('${:,.2f}'.format)

        return df

















