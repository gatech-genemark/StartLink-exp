# Karl Gemayel
# Georgia Institute of Technology
#
# Modified: 06/19/2020

import os
from shutil import copyfile
from typing import *

from sbsp_general import Environment
from sbsp_general.general import os_join
from sbsp_io.general import read_rows_to_list, mkdir_p
from sbsp_options.pipeline_sbsp import PipelineSBSPOptions
from sbsp_alg.sbsp_steps import sbsp_step_compare, sbsp_steps


class PipelineSBSP:
    """Class that can run the SBSP pipeline"""

    class PipelineState:
        def __init__(self, list_pf_data):
            self.use_pbs = True  # FIXME: allow non-pbs option
            self.list_pf_data = list_pf_data

            self._only_keep_files_that_exist()

        def _only_keep_files_that_exist(self):
            # type: () -> ()
            self.list_pf_data = [curr for curr in self.list_pf_data if os.path.exists(curr)]

        @classmethod
        def from_file(cls, pf_list_pf_data):
            # type: (str) -> PipelineSBSP.PipelineState

            list_pf_data = read_rows_to_list(pf_list_pf_data)
            return cls(list_pf_data)

    def __init__(self, env, pipeline_options, **kwargs):
        # type: (Environment, PipelineSBSPOptions, Dict[str, Any]) -> None

        self.env = env
        self.pipeline_options = pipeline_options

    def run(self):
        # type: () -> None
        pd_work = self.env["pd-work"]

        # make sure working directory is up and running
        mkdir_p(pd_work)

        # Copy genome file to local directory, and write sbsp options
        copyfile(self.pipeline_options["pf-q-list"], os_join(pd_work, "run.list"))
        self.pipeline_options["sbsp_options"].to_file(os_join(pd_work, "sbsp_options.conf"))

        state = self._run_helper()          # run compute steps
        self._compare(state)                # run comparison

    def _compare(self, state):
        # type: (PipelineState) -> PipelineState
        curr_env = self.env
        result = sbsp_step_compare(curr_env, self.pipeline_options, state.list_pf_data)

        return PipelineSBSP.PipelineState(result)

    def _run_helper(self, ):
        # type: () -> PipelineState
        """
        Runs the steps for prediction of gene-starts
        :return: State of pipeline after successful run
        """

        curr_env = self.env.duplicate({
            "pd-work": os.path.join(self.env["pd-work"], "steps")
        })
        result = sbsp_steps(curr_env, self.pipeline_options)

        return PipelineSBSP.PipelineState(result)
