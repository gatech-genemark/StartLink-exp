import logging
from typing import *

from sbsp_general.labels import Labels
from sbsp_general.labels_comparison_detailed import LabelsComparisonDetailed

logger = logging.getLogger(__name__)


class MultipleLCD:

    def __init__(self, list_lcd, **kwargs):
        # type: (List[LabelsComparisonDetailed], Dict[str, Any]) -> None

        self.list_lcd = list_lcd

    def get_comparisons_per_tag(self, attribute_info=None):
        # type: (Union[Tuple[str, Any], None]) -> Dict[str, Dict[str, Any]]

        dict_comparisons = dict()

        for tag, lcd in self.tag_to_lcd.items():
            if attribute_info is None:
                dict_comparisons[tag] = lcd.comparison["all"]

            else:
                name, value = attribute_info

                if name in lcd.comparison["attribute"] and value in lcd.comparison["attribute"][name]:
                    dict_comparisons[tag] = lcd.comparison["attribute"][name][value]
                else:
                    dict_comparisons[tag] = LabelsComparisonDetailed(Labels(), Labels())

        return dict_comparisons






    def _parse(self):
        # type: () -> None

        tag_to_lcd = dict()
        tag_unique = 0
        for lcd in self.list_lcd:
            tag = lcd.tag
            if tag is None:
                tag = tag_unique
                tag_unique += 1

            tag_to_lcd[tag] = lcd

        self.tag_to_lcd = tag_to_lcd

