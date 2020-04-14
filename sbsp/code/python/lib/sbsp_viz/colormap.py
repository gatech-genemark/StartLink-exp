import logging
from typing import *
import seaborn

logger = logging.getLogger(__name__)


seaborn.set_palette(seaborn.xkcd_palette(["windows blue", "amber", "faded green", "dusty purple"]))


def _init_mapping_ancestors():
    colors = ["windows blue", "amber", "faded green", "dusty purple"]
    ancestors = ["Archaea", "Actinobacteria", "Enterobacterales", "FCB group"]
    palette = seaborn.xkcd_palette(colors)
    return {x[0]: x[1] for x in zip(ancestors, palette)}


def _init_mapping_independence_conditions():
    colors = ["windows blue", "amber", "faded green"]
    conditions = ["Random", "Independent", "Fully dependent"]
    palette = seaborn.xkcd_palette(colors)
    return {x[0]: x[1] for x in zip(conditions, palette)}


class ColorMap:

    _mappings = {
        "ancestor": _init_mapping_ancestors(),
        "independence-conditions": _init_mapping_independence_conditions()
    }

    @staticmethod
    def get_map(name):
        # type: (str) -> Dict[str, Any]

        if name not in ColorMap._mappings:
            raise ValueError("Unknown color mapping for: {}".format(name))

        return ColorMap._mappings[name]


