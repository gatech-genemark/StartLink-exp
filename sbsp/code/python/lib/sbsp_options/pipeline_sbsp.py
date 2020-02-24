import os
from typing import *
from sbsp_general import Environment
from sbsp_options.sbsp import SBSPOptions
from sbsp_options.options import Options


class PipelineSBSPOptions(Options):
    """Options for SBSP pipeline"""

    def __init__(self, env, pf_q_list, pf_t_db, pf_output, sbsp_optionos, pf_options_custom=None, **kwargs):
        # type: (Environment, str, str, str, SBSPOptions, str, Dict[str, Any]) -> None
        super(PipelineSBSPOptions, self).__init__(env, pf_options_custom,
                                                  sbsp_options=sbsp_optionos,
                                                  pf_output=pf_output,
                                                  pf_q_list=pf_q_list,
                                                  pf_t_db=pf_t_db,
                                                  **kwargs)

    def perform_step(self, step):
        if "steps" not in self._options or self["steps"] is None:
            return True
        return step in self["steps"]

    def path_to_default_options_file(self, env):
        # type: (Environment) -> str
        return os.path.join(env["pd-config"], "pipeline_sbsp_defaults.conf")

    def required(self):
        # type: () -> Union[Set[str], None]
        return {
            # output files for sbsp steps
            "fn-orthologs", "fn-compute-features", "fn-filter", "fn-msa", "fn-accuracy", "pf-output",
            # input files
            "pf-q-list", "pf-t-db", "fn-q-labels"
        }

    def use_pbs(self):
        # type: () -> bool
        return "pbs-options" in self and self["pbs-options"]["use-pbs"]

    @staticmethod
    def init_from_dict(env, dict_options):
        # type: (Environment, Dict[str, Any]) -> TypeVar('T', bound=Options)

        pf_custom_options = None
        if "pf_pipeline_sbsp_options" in dict_options:
            pf_custom_options = dict_options["pf_pipeline_sbsp_options"]

        options = PipelineSBSPOptions(env, pf_custom_options=pf_custom_options, **dict_options)

        valid_keys = options._options.keys()

        for k in valid_keys:
            k_in_args = k.replace("-", "_")
            if k_in_args in dict_options and dict_options[k_in_args] is not None:
                options[k] = dict_options[k_in_args]

        return options
