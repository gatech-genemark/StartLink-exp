import os
from typing import *
from sbsp_general import Environment
from sbsp_options.msa import MSAOptions
from sbsp_options.options import Options


class PipelineSBSPOptions(Options):
    """Options for SBSP pipeline"""

    def __init__(self, env, pf_q_list, pf_t_list, pf_output, msa_options, pf_options_custom=None, **kwargs):
        # type: (Environment, str, str, str, MSAOptions, str, Dict[str, Any]) -> None
        super(PipelineSBSPOptions, self).__init__(env, pf_options_custom,
                                                  msa_options=msa_options,
                                                  pf_output=pf_output,
                                                  pf_q_list=pf_q_list,
                                                  pf_t_list=pf_t_list,
                                                  **kwargs)

    def perform_step(self, step):
        if "steps" not in self._options or self["steps"] is None:
            return True
        return step in self["steps"]

    def path_to_default_options_file(self, env):
        # type: (Environment) -> str
        return os.path.join(env["pd-config"], "pipeline-sbsp-defaults.conf")

    def required(self):
        # type: () -> Union[Set[str], None]
        return {
            # output files for sbsp steps
            "fn-orthologs", "fn-compute-features", "fn-filter", "fn-msa", "fn-accuracy", "pf-output",
            # input files
            "pf-q-list", "pf-t-list", "fn-q-labels", "fn-t-labels"
        }

    def use_pbs(self):
        # type: () -> bool
        return "pbs-options" in self and self["pbs-options"]["use-pbs"]

    @staticmethod
    def init_from_dict(env, dict_options):
        # type: (Environment, Dict[str, Any]) -> TypeVar('T', bound=Options)

        pf_custom_options = None
        if "pf_msa_options" in dict_options:
            pf_custom_options = dict_options["pf_msa_options"]

        msa_options = PipelineSBSPOptions(env,
                                          pf_custom_options=pf_custom_options, **dict_options)

        # valid_keys = msa_options._options.keys()
        #
        # for k in valid_keys:
        #     k_in_args = k.replace("-", "_")
        #     if k_in_args in dict_options and dict_options[k_in_args] is not None:
        #         msa_options[k] = dict_options[k_in_args]

        return msa_options


#
# class PipelineMSAOptions:
#
#     def __init__(self, env, **kwargs):
#
#         # column names
#         self.suffix_coordinates = get_value(kwargs, "suffix_coordinates", None)
#         self.suffix_gene_sequence = get_value(kwargs, "suffix_gene_sequence", "gene-sequence")
#         self.column_k2p_distance = get_value(kwargs, "k2p_distance", "k2p-distance")
#         self.tag_msa = get_value(kwargs, "tag_msa", "msa")
#
#         # pbs options
#         self.pbs_options = get_value(kwargs, "pbs_options", PBSOptions(env))
#
#         # directories
#         self.dn_msa_output = get_value(kwargs, "dn_msa_output", None)
#
#         # paths to files
#         self.pf_orthologs = get_value(kwargs, "pf_orthologs", os.path.join(env["pd-work"], "orthologs.csv"))
#         self.pf_features = get_value(kwargs, "pf_features", os.path.join(env["pd-work"], "features.csv"))
#         self.pf_filter = get_value(kwargs, "pf_filter", os.path.join(env["pd-work"], "filter.csv"))
#         self.pf_msa = get_value(kwargs, "pf_msa", os.path.join(env["pd-work"], "msa.csv"))
#         self.pf_accuracy = get_value(kwargs, "pf_accuracy", os.path.join(env["pd-work"], "accuracy.csv"))
#
#         self.pf_q_list = get_value(kwargs, "pf_q_list", None)
#         self.pf_t_list = get_value(kwargs, "pf_t_list", None)
#
#         # file names
#         self.fn_q_labels = get_value(kwargs, "fn_q_labels", "verified.gff")
#         self.fn_t_labels = get_value(kwargs, "fn_t_labels", "ncbi.gff")
#         self.fn_q_labels_true = get_value(kwargs, "fn_q_labels_true", None)
#         self.suffix_labels = get_value(kwargs, "suffix_labels", "")
#
#         # options
#         self.upstream_length_nt = get_value(kwargs, "upstream_length_nt", None)
#         self.downstream_length_nt = get_value(kwargs, "downstream_length_nt", None)
#         self.steps = get_value(kwargs, "steps", None)
#         self.multiple_msa_start_selection_threshold = get_value(kwargs, "multiple_msa_start_selection_thresholds", None)
#
#         self.msa_options = get_value(kwargs, "msa_options",
#                                 sbsp_options.msa.MSAOptions(env))  # type: sbsp_alg.msaoptions.MSAOptions
#
#     def perform_step(self, step):
#         # type: (str) -> bool
#
#         if self.steps is None:
#             return True
#
#         return step in self.steps
