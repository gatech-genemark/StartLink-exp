import os

from typing import *
from sbsp_general import Environment
from sbsp_options.options import Options


class PBSOptions(Options):
    """Options for PBS scheduler"""

    def path_to_default_options_file(self, env):
        # type: (Environment) -> str
        return os.path.join(env["pd-config"], "pbs-defaults.conf")

    def required(self):
        # type: () -> Union[Set[str], None]
        return {
            # directories for pbs
            "pd-root-compute", "dn-compute", "pd-head",
            # node computation configuration
            "num-processors", "num-nodes-per-job", "walltime",
            "num-jobs",
            "split-tag",
            "use-pbs"
        }

    @staticmethod
    def init_from_dict(env, dict_options):
        # type: (Environment, Dict[str, Any]) -> TypeVar('T', bound=Options)

        pf_custom_options = None
        if "pf_msa_options" in dict_options:
            pf_custom_options = dict_options["pf_msa_options"]

        msa_options = PBSOptions(env, pf_custom_options)

        valid_keys = msa_options._options.keys()

        for k in valid_keys:
            k_in_args = k.replace("-", "_")
            if k_in_args in dict_options and dict_options[k_in_args] is not None:
                msa_options[k] = dict_options[k_in_args]

        return msa_options