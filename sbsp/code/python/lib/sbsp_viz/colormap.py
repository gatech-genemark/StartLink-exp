import logging
from typing import *
import seaborn

logger = logging.getLogger(__name__)


class ColorMap:
    colors = ["windows blue", "amber", "faded green", "dusty purple"]
    ancestors = ["Archaea", "Actinobacteria", "Enterobacterales", "FCB group"]
    palette = seaborn.xkcd_palette(colors)
    color_mapping = {x[0]: x[1] for x in zip(ancestors, palette)}

    _mappings = {
        "ancestor": color_mapping
    }

    seaborn.set_palette(palette)

    @staticmethod
    def get_map(name):
        # type: (str) -> Dict[str, Any]

        if name not in ColorMap._mappings:
            raise ValueError("Unknown color mapping for: {}".format(name))

        return ColorMap._mappings[name]
