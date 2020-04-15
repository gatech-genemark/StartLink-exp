import logging
import timeit
import pandas as pd
from typing import *

from sbsp_general.general import except_if_not_in_set, get_value
from sbsp_io.general import write_string_to_file

logger = logging.getLogger(__name__)

class Timer:

    seconds_to_unit_multiplier = {
        "ms": 60.0,
        "s": 1.0,
        "m": 1/60.0,
        "h": 1/3600.0
    }

    def __init__(self):
        self._timers = dict()

    def start(self, name):
        # type: (str) -> None

        if self._has_finished(name):
            logger.warning("Timer '{}' has already completed. Restarting...".format(name))
        elif self._has_started(name):
            logger.warning("Timer '{}' has already started. Restarting...".format(name))

        self._timers[name] = {
            "begin": timeit.default_timer(),
            "end": None
        }

    def finish(self, name, unit="s"):
        # type: (str, str) -> float
        except_if_not_in_set(unit, Timer.seconds_to_unit_multiplier.keys())

        if not self._has_started(name):
            raise ValueError("Cannot finish timer '{}'. It hasn't started.".format(name))

        self._timers[name]["end"] = timeit.default_timer()
        return self.elapsed_time(name, unit)

    def elapsed_time(self, name, unit):
        # type: (str, str) -> float
        """Checks if a timer has already been completed"""
        if not self._has_started(name):
            raise ValueError("Cannot compute elapsed time. Timer '{}' has not started.".format(name))
        if not self._has_finished(name):
            raise ValueError("Cannot compute elapsed time. Timer '{}' has not finished.".format(name))

        try:
            return Timer._from_seconds_to_unit(self._timers["end"] - self._timers["begin"], unit)
        except ValueError as e:
            raise e

    def to_csv(self, pf_csv, **kwargs):
        # type: (str) -> None

        out = self.to_string(**kwargs)
        write_string_to_file(out, pf_csv)

    def to_string(self, **kwargs):
        # type: (Dict[str, Any]) -> str

        order = get_value(kwargs, "order", "ascending", valid_choices={"ascending", "descending"})
        order_by = get_value(kwargs, "order_by", "start", valid_choices={"start", "end"})

        df = pd.DataFrame(
            [{
                "Name": name, **self._timers[name]
            } for name in self._timers.keys()]
        )

        ascending = True if order == "ascending" else False
        df.sort_values(order_by, ascending=ascending, inplace=True)

        return df.to_string(index=False)










    def _has_started(self, name):
        # type: (str) -> bool
        """Checks if a timer has already been started (whether or not it has completed)"""
        return name in self._timers

    def _has_finished(self, name):
        # type: (str) -> bool
        """Checks if a timer has already been completed"""
        return name in self._timers and self._timers[name]["end"] is not None




    @staticmethod
    def _from_seconds_to_unit(seconds, unit):
        # type: (float, str) -> float
        except_if_not_in_set(unit, Timer.seconds_to_unit_multiplier.keys())
        return seconds * Timer.seconds_to_unit_multiplier[unit]




