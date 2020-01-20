import os
import timeit
from typing import *

from sbsp_alg.sbsp_steps import sbsp_step_get_orthologs, sbsp_step_compute_features, sbsp_step_filter, sbsp_step_msa, \
    sbsp_step_accuracy
from sbsp_general import Environment
from sbsp_io.general import read_rows_to_list
from sbsp_options.msa import MSAOptions
from sbsp_options.pipeline_sbsp import PipelineSBSPOptions


class PipelineMSA:
    class PipelineState:
        def __init__(self, list_pf_data):
            self.use_pbs = True  # FIXME: allow non-pbs option
            self.list_pf_data = list_pf_data

            self._keep_only_exists()

        def _keep_only_exists(self):
            # type: () -> ()
            self.list_pf_data = [curr for curr in self.list_pf_data if os.path.exists(curr)]

        @classmethod
        def from_file(cls, pf_list_pf_data):
            # type: (str) -> PipelineMSA.PipelineState

            list_pf_data = read_rows_to_list(pf_list_pf_data)
            return cls(list_pf_data)

    def __init__(self, env, pipeline_options, **kwargs):
        # type: (Environment, PipelineSBSPOptions, Dict[str, Any]) -> None

        self.env = env
        self.pipeline_options = pipeline_options

    def run(self):
        # type: () -> None

        elapsed_times = dict()

        curr_time = timeit.default_timer()
        state = self._run_get_orthologs()
        elapsed_times["1-orthologs"] = timeit.default_timer() - curr_time

        curr_time = timeit.default_timer()
        state = self._run_compute_features(state)
        elapsed_times["2-features"] = timeit.default_timer() - curr_time

        curr_time = timeit.default_timer()
        state = self._run_filter(state)
        elapsed_times["3-filter"] = timeit.default_timer() - curr_time

        curr_time = timeit.default_timer()
        state = self._run_msa(state)
        elapsed_times["4-msa"] = timeit.default_timer() - curr_time

        curr_time = timeit.default_timer()
        state = self._accuracy(state)
        elapsed_times["5-accuracy"] = timeit.default_timer() - curr_time

        time_string = "\n".join([
                "{},{}".format(key, value) for key, value in elapsed_times.items()
        ])

        pf_time = os.path.join(self.env["pd-work"], "time.csv")
        with open(pf_time, "w") as f:
            f.write(time_string)
            f.close()

    def _run_get_orthologs(self):
        # type: () -> PipelineState
        curr_env = self.env.duplicate({
            "pd-work": os.path.join(self.env["pd-work"], "orthologs")
        })
        result = sbsp_step_get_orthologs(curr_env, self.pipeline_options)

        return PipelineMSA.PipelineState(result)

    def _run_compute_features(self, state):
        # type: (PipelineState) -> PipelineState
        curr_env = self.env.duplicate({
            "pd-work": os.path.join(self.env["pd-work"], "features")
        })
        result = sbsp_step_compute_features(curr_env, self.pipeline_options, state.list_pf_data)

        return PipelineMSA.PipelineState(result)

    def _run_filter(self, state):
        # type: (PipelineState) -> PipelineState
        curr_env = self.env.duplicate({
            "pd-work": os.path.join(self.env["pd-work"], "filter")
        })
        result = sbsp_step_filter(curr_env, self.pipeline_options, state.list_pf_data)

        return PipelineMSA.PipelineState(result)

    def _run_msa(self, state):
        # type: (PipelineState) -> PipelineState
        curr_env = self.env.duplicate({
            "pd-work": os.path.join(self.env["pd-work"], "msa")
        })
        result = sbsp_step_msa(curr_env, self.pipeline_options, state.list_pf_data)

        return PipelineMSA.PipelineState(result)

    def _accuracy(self, state):
        # type: (PipelineState) -> PipelineState
        curr_env = self.env.duplicate({
            "pd-work": os.path.join(self.env["pd-work"], "accuracy")
        })
        result = sbsp_step_accuracy(curr_env, self.pipeline_options, state.list_pf_data)

        return PipelineMSA.PipelineState(result)
