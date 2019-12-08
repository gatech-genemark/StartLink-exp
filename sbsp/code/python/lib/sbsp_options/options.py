import os
import yaml
import logging
from typing import *
from abc import abstractmethod

import sbsp_io.general
from sbsp_general import Environment

logger = logging.getLogger(__name__)


class Options:
    """Generic class for holding configuration options
    All options have a dash (-) separator instead of an underscore (_)
    """

    def __init__(self, env, pf_options_custom=None, **kwargs):
        # type: (Environment, str, Dict[str, Any]) -> None

        self.env = env.duplicate()

        # get default options
        self._default_options = self.get_default_options(env)

        # get custom options
        self._custom_options = Options.read_from_file(pf_options_custom)

        # update defaults with custom options
        self._options = Options.merge_custom_with_default(self._default_options, self._custom_options)

        # update options with kwargs
        Options.merge_kwargs_with_options_dict(self._options, **kwargs)

        # check that required parameters have been set
        self._check_requirements()

    def get_default_options(self, env):
        # type: (Environment) -> Dict[str, Any]
        pf_default = self.path_to_default_options_file(env)
        return Options.read_from_defaults_file(pf_default)

    def __getitem__(self, item):
        # type: (str) -> Any
        return self._options[item]

    def __setitem__(self, item, value):
        # type: (str, Any) -> None
        self._options[item] = value

    def __contains__(self, item):
        # type: (Any) -> bool
        return item in self._options

    def option_names(self):
        # type: () -> KeysView[str]
        return self._options.keys()

    def to_string(self):
        # type: () -> str
        return yaml.dump(self._options)

    def to_file(self, pf_options):
        # type: (str) -> None
        sbsp_io.general.write_string_to_file(self.to_string(), pf_options)

    @abstractmethod
    def path_to_default_options_file(self, env):
        # type: (Environment) -> str
        # Needs to be implemented by child class, and returns the path to the file containing
        # the class's default options
        pass

    @staticmethod
    def read_from_file(pf_options):
        # type: (Union[str, None]) -> Dict[str, Any]
        if pf_options is None:
            return dict()

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
        # type: (Environment, Dict[str, Any]) -> TypeVar('T', bound=Options)

        pf_custom_options = None
        if "pf_msa_options" in dict_options:
            pf_custom_options = dict_options["pf_msa_options"]

        msa_options = Options(env, pf_custom_options)

        valid_keys = msa_options._options.keys()

        for k in valid_keys:
            k_in_args = k.replace("-", "_")
            if k_in_args in dict_options and dict_options[k_in_args] is not None:
                msa_options[k] = dict_options[k_in_args]

        return msa_options

    @staticmethod
    def merge_kwargs_with_options_dict(options, **kwargs):
        # type: (Dict[str, Any], Dict[str, Any]) -> None
        for k, v in kwargs.items():
            options[k.replace("_", "-")] = v

    def _check_requirements(self):
        # type: () -> None
        requirements = self.required()
        if requirements is None:
            return

        for r in requirements:
            if r not in self._options or self._options[r] is None:
                raise ValueError("Option required: {}".format(r))

    def required(self):
        # type: () -> Union[Set[str], None]
        """Returns a set of required options to be set"""
        return None









