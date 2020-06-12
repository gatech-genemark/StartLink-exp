import os
import logging
from typing import *
from sbsp_general import Environment
from sbsp_general.general import os_join, get_value
from sbsp_options.options import Options

log = logging.getLogger(__name__)


class ParallelizationOptions(Options):
    """Options for parallelization of program executions, including PBS systems"""

    def __init__(self, env, pf_options_custom, **kwargs):
        # type: (Environment, str, Dict[str, Any]) -> None
        super(ParallelizationOptions, self).__init__(env, pf_options_custom, **kwargs)

        self._check_pbs_paths_for_directories()

    def path_to_default_options_file(self, env):
        # type: (Environment) -> str
        return os_join(env["pd-config"], "parallelization_defaults.conf")

    def required(self):  # type: () -> Union[Set[str], None]
        return None

    def _check_pbs_paths_for_directories(self):
        # type: () -> None

        # if pbs head node directory not specified, use current working directory
        if self._options["pbs-pd-head"] is None:
            self._options["pbs-pd-head"] = self.env["pd-work"]

        # make sure it's an absolute path
        self._options["pbs-pd-head"] = os.path.abspath(self._options["pbs-pd-head"])

        # if path to compute not specified, use PBS head directory
        if self._options["pbs-pd-root-compute"] is None:
            self._options["pbs-pd-root-compute"] = self._options["pbs-pd-head"]

    @classmethod
    def init_from_dict(cls, env, dict_options):
        # type: (Environment, Dict[str, Any]) -> TypeVar("T", bound=Options)

        # start by reading options from custom file
        pf_custom_options = get_value(dict_options, "pf_parallelization_options", None)

        options = ParallelizationOptions(env, pf_custom_options)

        # check if any valid keys in dictionary to update

        valid_keys = options._options.keys()
        for k in valid_keys:
            key_in_dict = k.replace("-", "_")
            if key_in_dict in dict_options and dict_options[key_in_dict] is not None:
                options[k] = dict_options[key_in_dict]

        return options
