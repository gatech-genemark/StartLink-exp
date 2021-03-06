import os
import copy
from typing import *
from sbsp_general import general

TOOL = "StartLink"
TOOLp = f"{TOOL}+"


class Environment:
    """A class representing the project environment variables, including
    paths to useful directories (data, code, etc...)
    """

    def __init__(self, pd_data=None, pd_work=None, **kwargs):
        # type: (str, str, Dict[str, Any]) -> None

        self._env = Environment.load_environment_variables(pd_data, pd_work, **kwargs)

    def __getitem__(self, item):
        # type: (str) -> Any
        return self._env[item]

    def __setitem__(self, key, value):
        # type: (str, Any) -> None
        self._env[key] = value

    def duplicate(self, new_values=None):
        # type: (Dict[str, Any]) -> Environment
        """Creates a copy of the environment, with update variables
        """
        new_env = copy.deepcopy(self)
        if new_values is not None:
            for item in new_values.keys():
                new_env[item] = new_values[item]

        return new_env

    @staticmethod
    def load_environment_variables(pd_data=None, pd_work=None, **kwargs):
        # type: (str, str, Dict[str, Any]) -> Dict[str, str]

        import os

        if pd_data is not None:
            pd_data = os.path.abspath(pd_data)
        if pd_work is not None:
            pd_work = os.path.abspath(pd_work)

        pd_results = kwargs["pd_results"] if "pd_results" in kwargs else None

        # path to current file
        pd_curr = os.path.abspath(general.get_file_dir(__file__))

        pd_base = os.path.join(pd_curr, "../../../../")
        pd_bin = os.path.join(pd_base, "bin")
        pd_code = os.path.join(pd_base, "code")
        pd_tmp = os.path.join(pd_base, "tmp")
        pd_mat = os.path.join(pd_base, "mat")
        pd_config = os.path.join(pd_base, "config")

        if pd_results is None:
            pd_results = os.path.join(pd_base, "results")

        if pd_data is None:
            pd_data = os.path.join(pd_base, "data/all")
#            pd_data = os.path.join("/storage4/karl/sbsp/similarity-based-start-prediction/data/large")

        if pd_work is None:
            pd_work = os.path.abspath(".")

        pd_work_results = os.path.join(pd_results, os.path.basename(pd_work))
        if pd_work == os.path.abspath("."):
            pd_work_results = pd_work

        if not os.path.exists(pd_work):
            os.makedirs(pd_work)

        if not os.path.exists(pd_work_results):
            os.makedirs(pd_work_results)

        pd_bash = os.path.join(pd_code, "bash")
        pd_bash_lib = os.path.join(pd_bash, "lib")
        pd_bash_driver = os.path.join(pd_bash, "lib")

        env = {
            "base": pd_base,
            "bin": pd_bin,  # deprecated
            "code": pd_code,
            "data": pd_data,
            "tmp": pd_tmp,
            "mat": pd_mat,
            "bash": pd_bash,
            "bash-lib": pd_bash_lib,
            "bash-driver": pd_bash_driver,
            "working-dir": pd_work,  # deprecated
            "pd-bash-lib": pd_bash_lib,
            "pd-base": pd_base,
            "pd-results": pd_results,
            "pd-work": pd_work,
            "pd-code": pd_code,
            "pd-data": pd_data,
            "pd-bin": pd_bin,
            "pd-config": pd_config,
            "pd-work-results": pd_work_results,
            "pd-bin-external": os.path.join(pd_base, "bin_external"),
            "pd-runs": os.path.join(pd_base, "runs") 
        }

        import copy
        global ENV
        ENV = copy.deepcopy(env)

        return env


ENV = Environment()
