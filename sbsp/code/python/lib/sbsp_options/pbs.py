from __future__ import print_function
import os
import yaml
from typing import *

import sbsp_io.general

import logging
logger = logging.getLogger(__name__)


class PBSOptions:

    def __init__(self, env, pf_options_custom=None):

        pf_default = os.path.join(env["pd-config"], "pbs-defaults.conf")
        self._default_options = PBSOptions.read_from_defaults_file(pf_default)

        self._custom_options = dict()
        if pf_options_custom is not None:
            self._custom_options = PBSOptions.read_from_file(pf_options_custom)

        self._options = PBSOptions.merge_custom_with_default(self._default_options, self._custom_options)

        if self._options["pd-work-pbs"] is None:
            self._options["pd-work-pbs"] = env["pd-work"]

    def __getitem__(self, item):
        # type: (str) -> Any
        return self._options[item]

    def __setitem__(self, item, value):
        # type: (str, Any) -> None
        self._options[item] = value

    def option_names(self):
        # type: () -> KeysView[str]
        return self._options.keys()

    def to_string(self):
        # type: () -> str
        return yaml.dump(self._options)

    def to_file(self, pf_options):
        # type: (str) -> None
        sbsp_io.general.write_string_to_file(self.to_string(), pf_options)

    @staticmethod
    def read_from_file(pf_options):
        # type: (str) -> Dict[str, Any]
        try:
            f = open(pf_options, "r")
            return yaml.load(f)
        except IOError:
            logger.warning("Options File Not Found: {}".format(pf_options))
            return dict()


    @staticmethod
    def read_from_defaults_file(pf_default):
        # type: (str) -> Dict[str, Any]

        try:
            f = open(pf_default, "r")
            return yaml.load(f)
        except IOError:
            logger.warning("Defaults File Not Found: {}".format(pf_default))
            return dict()

    @staticmethod
    def merge_custom_with_default(default, custom):
        # type: (Dict[str, Any], Dict[str, Any]) -> Dict[str, Any]

        if default is None and custom is None:
            return dict()

        if default is None:
            return custom
        if custom is None:
            return default

        import copy
        combined = copy.deepcopy(default)
        combined.update(custom)
        return combined

    @staticmethod
    def init_from_dict(env, dict_options):
        # type: (Dict[str, Any], Dict[str, Any]) -> PBSOptions

        pf_custom_options = None
        if "pf_pbs_options" in dict_options:
            pf_custom_options = dict_options["pf_pbs_options"]
            print (pf_custom_options)

        options = PBSOptions(env, pf_custom_options)

        valid_keys = options._options.keys()

        for k in valid_keys:
            k_in_args = k.replace("-", "_")
            if k_in_args in dict_options and dict_options[k_in_args] is not None:
                options[k] = dict_options[k_in_args]

        return options












