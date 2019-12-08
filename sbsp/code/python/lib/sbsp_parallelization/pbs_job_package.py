from typing import *
import cloudpickle as dill


def save_obj(obj, name):
    # type: (object, str) -> None
    with open(name + '.pkl', 'wb') as f:
        dill.dump(obj, f)


def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        loaded_data = dill.load(f)
        return loaded_data


class PBSJobPackage:

    @staticmethod
    def save(pbs_job_package, pf_save_to):
        # type: (Dict[str, Any], str) -> None
        save_obj(pbs_job_package, pf_save_to)

    @staticmethod
    def load(pf_load_from):
        # type: (str) -> PBSJobPackage
        return load_obj(pf_load_from)


class PBSJobInputPackage:
    """A class that represents an input to a PBS Job"""

    def __init__(self, func, func_arguments):
        self._package = {
            "func": func,
            "func_arguments": func_arguments
        }


class PBSJobOutputPackage:
    """A class that represents an output from a PBS Job"""

    def __init__(self):
        pass