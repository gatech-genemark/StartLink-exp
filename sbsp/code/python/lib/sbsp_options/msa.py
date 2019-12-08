from __future__ import print_function
import os
import copy
import yaml
import pandas as pd
from typing import *

import sbsp_io.general

import logging

from sbsp_general import Environment

logger = logging.getLogger(__name__)


class MSAOptions:

    def __init__(self, env, pf_options_custom = None):

        pf_default = os.path.join(env["pd-config"], "msa-defaults.conf")
        self._default_options = MSAOptions.read_from_defaults_file(pf_default)

        self._custom_options = dict()
        if pf_options_custom is not None:
            self._custom_options = MSAOptions.read_from_file(pf_options_custom)

        self._options = MSAOptions.merge_custom_with_default(self._default_options, self._custom_options)

    def __getitem__(self, item):
        # type: (str) -> Any
        return self._options[item]

    def __setitem__(self, item, value):
        # type: (str, Any) -> None
        self._options[item] = value

    def __contains__(self, item):
        return item in self._options.keys()

    def option_names(self):
        # type: () -> KeysView[str]
        return self._options.keys()

    def to_string(self):
        # type: () -> str
        return yaml.dump(self._options)

    def to_series(self):
        # type: () -> pd.Series
        return pd.Series(self._options)

    def to_dict(self):
        # type: () -> Dict[str, Any]
        return copy.deepcopy(self._options)

    def to_file(self, pf_options):
        # type: (str) -> None
        sbsp_io.general.write_string_to_file(self.to_string(), pf_options)

    def safe_get(self, item):
        # type: (str) -> Any
        if item in self:
            return self[item]
        return None

    def update(self, other):
        # type: (Dict[str, Any]) -> None

        for key, value in other.items():
            self._options[key] = value

    @staticmethod
    def read_from_file(pf_options):
        # type: (str) -> Dict[str, Any]
        try:
            f = open(pf_options, "r")
            # return yaml.load(f, Loader=yaml.FullLoader)
            return yaml.load(f)
        except IOError:
            logger.warning("MSA-Options File Not Found: {}".format(pf_options))
            return dict()


    @staticmethod
    def read_from_defaults_file(pf_default):
        # type: (str) -> Dict[str, Any]

        try:
            f = open(pf_default, "r")
            return yaml.load(f)
        except IOError:
            logger.warning("MSA-Defaults File Not Found: {}".format(pf_default))
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
        # type: (Environment, Dict[str, Any]) -> MSAOptions

        pf_custom_options = None
        if "pf_msa_options" in dict_options:
            pf_custom_options = dict_options["pf_msa_options"]
            print (pf_custom_options)

        msa_options = MSAOptions(env, pf_custom_options)

        valid_keys = msa_options._options.keys()

        for k in valid_keys:
            k_in_args = k.replace("-", "_")
            if k_in_args in dict_options and dict_options[k_in_args] is not None:
                msa_options[k] = dict_options[k_in_args]

        return msa_options












