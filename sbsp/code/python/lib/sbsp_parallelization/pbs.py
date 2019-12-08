import os
import copy
from typing import *

from sbsp_general import Environment
from sbsp_io.general import mkdir_p
from sbsp_options.pbs import PBSOptions

from sbsp_general.general import get_value, run_shell_cmd
import sbsp_io.general
from sbsp_parallelization.pbs_job_package import PBSJobPackage


class FunctionArguments:

    def __init__(self, **kwargs):

        self._kwargs = kwargs

    def get_arguments(self, data):
        # type: (Dict[str, Any]) -> Dict[str, Any]
        new_kwargs = copy.deepcopy(self._kwargs)
        new_kwargs.update(data)

        return new_kwargs


class PBS:
    """Runs any function on input data using PBS scheduler"""

    def __init__(self, env, pbs_options, splitter, merger, **kwargs):
        # type: (Environment, PBSOptions, Callable, Callable, Dict[str, Any]) -> None
        """Create a PBS instance that can run a function on
        :param env:
        :param pbs_options:
        :param data_kind:
        :param data_format:
        :param kwargs:
        """

        self._dry_run = get_value(kwargs, "dry_run", False)
        self._env = env

        if pbs_options is None:
            raise ValueError("PBSOptions cannot be None")

        self._pbs_options = copy.deepcopy(pbs_options)
        self._pbs_options["pd-head"] = os.path.abspath(self._pbs_options["pd-head"])
        self._pbs_options["pd-root-compute"] = os.path.abspath(self._pbs_options["pd-root-compute"])

        self._splitter = splitter
        self._merger = merger

    def _setup_pbs_run(self):

        mkdir_p(self._pbs_options["pd-head"])

    def run(self, data, func, func_kwargs, **kwargs):
        # type: (Dict[str, Any], Callable, Dict[str, Any], Dict[str, Any]) -> Any
        """
        Run function on data using PBS scheduler
        :param data: Dictionary containing data name and value
        :param func: function to call on data with arguments in func_kwargs
        :param func_kwargs: Additional arguments to be passed
        :param kwargs:
        :return: output of merger function
        """



        job_name = get_value(kwargs, "job_name", "JOBNAME")
        pf_input_package_template_formatted = os.path.join(
            self._pbs_options["pd-head"], "input_package_{}"
        )

        num_jobs = self._pbs_options["num-jobs"]

        # 1) Parse all PBS arguments
        self._setup_pbs_run()

        # 2) Create input packages files, one for every PBS run
        self.create_input_package_files(
            data, func, func_kwargs, num_jobs,
            pf_package_template_formatted=pf_input_package_template_formatted,
            **kwargs
        )

        # 3) Run all
        list_pf_output_job_packages = self.execute_function_on_input_packages(
            pf_input_package_template_formatted,
            job_name=job_name, num_jobs=num_jobs
        )

        # 4) Merge end-results
        data_output = None
        if not self._dry_run:
            data_output = self.merge_output_package_files(list_pf_output_job_packages)

        return data_output

    def create_input_package_files(self, data, func, func_kwargs, num_splits, **kwargs):
        """
        Run a function on the data using PBS
        :param data: the entire data
        :type data: DataHandler.D
        :param data_arg_name: the name of the data argument in func
        :type data_arg_name: str
        :param func: the function to execute on the (split) data
        :type func: Callable
        :param func_kwargs: the remaining arguments (i.e. not data) to be passed to the function
        :type func_kwargs: Dict[str, Any]
        :param num_splits: number of job splits
        :type num_splits: int
        :param kwargs:
        :return: List of paths to input package files
        :rtype: List[str]
        """

        pd_work_pbs = self._pbs_options["pd-head"]

        pf_package_template_formatted = get_value(
            kwargs, "pf_package_template_formatted", os.path.join(pd_work_pbs, "input_package_{}")
        )

        # Split data
        list_split_data = self._splitter(data, num_splits, pd_work_pbs)

        # Write package to disk
        list_pf_data = self._package_and_save_list_data(list_split_data, func, func_kwargs,
                                                        pf_package_template_formatted)

        # return list of filenames
        return list_pf_data

    def execute_function_on_input_packages(self, pf_input_package_template_formatted, job_name, num_jobs):
        """
        Create PBS file for run and execute it, returning the paths to all the job output packages
        :param pf_input_package_template_formatted:
        :param job_name:
        :param num_jobs:
        :returns: list of paths to output file packages
        :rtype: str
        """

        pf_pbs = os.path.join(self._pbs_options["pd-head"], "run.pbs")
        pf_input_package_template = pf_input_package_template_formatted.format("${PBS_ARRAYID}")

        # create pbs file
        pf_output_package_template = "{}_output".format(pf_input_package_template)
        self._create_pbs_file(job_name, num_jobs, pf_pbs, pf_input_package_template, pf_output_package_template)

        # run
        if not self._dry_run:
            run_shell_cmd("qsub {}".format(pf_pbs))

        # write summary file
        list_pf_outputs = [PBS.create_concrete_from_template(pf_output_package_template, x) for x in range(num_jobs)]
        pf_pbs_summary = os.path.join(self._pbs_options["pd-head"], self._pbs_options["fn-pbs-summary"])
        sbsp_io.general.write_string_to_file("\n".join(list_pf_outputs), pf_pbs_summary)

        return list_pf_outputs

    def _read_data_from_output_packages(self, list_pf_output_packages):

        list_data = list()

        for pf_output_package in list_pf_output_packages:
            list_data.append(PBSJobPackage.load(pf_output_package)["data"])

        return list_data

    def merge_output_package_files(self, list_pf_output_packages):

        list_output_data = self._read_data_from_output_packages(list_pf_output_packages)

        # 4-a) Merge data while loading packages one by one
        data_output = self._merger(list_output_data)

        return data_output

    def _package_and_save_data(self, data, func, func_kwargs, pf_package):
        # type: (Dict[str, Any], Callable, Dict[str, Any], str) -> None

        complete_func_kwargs = FunctionArguments(**func_kwargs).get_arguments(data)

        PBSJobPackage.save(
            {
                "func": func,
                "func_kwargs":  complete_func_kwargs
            },
            pf_package
        )

    def _package_and_save_list_data(self, list_data, func, func_kwargs, pf_package_template_formatted):
        # type: (List[Dict[str, Any]], Callable, Dict[str, Any], str) -> List[str]

        list_pf = list()
        file_number = 0

        for data in list_data:
            pf_save = pf_package_template_formatted.format(file_number)

            self._package_and_save_data(data, func, func_kwargs, pf_save)
            list_pf.append(pf_save)

            file_number += 1

        return list_pf

    def _create_pbs_file(self, jobname, num_jobs, pf_pbs, pf_input_package_template, pf_output_package_template):
        """
        Create PBS file for runnning all input jobs
        :param jobname: Name of job
        :param num_jobs:
        :param pf_pbs:
        :param pf_input_package_template:
        :return:
        """

        pbs_text = PBS._generate_pbs_header_array(num_jobs, jobname, self._pbs_options)

        pbs_text += "\n{}\n".format(
            PBS._generate_call_command(self._env,
                                       pf_input_package_template,
                                       pf_output_package_template,
                                       self._pbs_options
            )
        )

        # write to file
        from sbsp_io.general import write_string_to_file
        write_string_to_file(pbs_text, pf_pbs)

    @staticmethod
    def _generate_pbs_header_array(num_jobs, job_name, pbs_options):
        """

        :param num_jobs:
        :param job_name:
        :param pbs_options:
        :type pbs_options: PBSOptions
        :return:
        """

        num_nodes = pbs_options["num-nodes-per-job"]
        ppn = pbs_options["num-processors"]
        walltime = pbs_options["walltime"]

        pd_compute = os.path.abspath(os.path.join(pbs_options["pd-root-compute"], pbs_options["dn-compute"]))

        pd_job_template = os.path.join(pd_compute, "job_${PBS_ARRAYID}")


        pbs_text = " "

        pbs_text += "#PBS -N " + str(job_name) + "\n"
        pbs_text += "#PBS -o " + "{}/{}".format(pd_job_template, "error_${PBS_ARRAYID}") + "\n"
        pbs_text += "#PBS -j oe" + "\n"
        pbs_text += "#PBS -l nodes=" + str(num_nodes) + ":ppn=" + str(ppn) + "\n"
        pbs_text += "#PBS -l walltime=" + str(walltime) + "\n"

        if pbs_options:
            array_param = "1-{}".format(num_jobs)
            if pbs_options["max-concurrent-nodes"]:
                total_concurrent_jobs = pbs_options["max-concurrent-nodes"] * int(8 / ppn)
                array_param = "{}%{}".format(array_param, total_concurrent_jobs)

            pbs_text += "#PBS -t {}".format(array_param) + "\n"

        pbs_text += "#PBS -W umask=002" + "\n"

        pbs_text += "export PATH=\"/home/karl/anaconda/envs/sbsp/bin:$PATH\"\n"

        pbs_text += "mkdir -p {}".format(pd_job_template) + "\n"

        pbs_text += "set PBS_O_WORKDIR = " + pd_job_template + "\n"
        pbs_text += "cd $PBS_O_WORKDIR \n"

        pbs_text += "echo The working directory is `echo $PBS_O_WORKDIR`" + "\n"
        pbs_text += "echo This job runs on the following nodes:" + "\n"
        pbs_text += "echo `cat $PBS_NODEFILE`" + "\n"

        return pbs_text

    @staticmethod
    def _generate_call_command(env, pf_job_input, pf_job_output, pbs_options):

        cmd = "{} --pf-job-input {} --pf-job-output {}".format(
            "{}/{}".format(env["pd-bin"], "run-pbs-job_py.sh"),
            pf_job_input,
            pf_job_output
        )

        return cmd

    @staticmethod
    def create_concrete_from_template(pf_template, file_number):
        """Create a concrete file name based on template and file number
        e.g. Calling the function with filename_${PBS_ARRAYID}.txt, 5 returns
        filename_5.txt

        :param pf_template: template of file
        :type pf_template: str
        :param file_number: the file's number
        :returns: a concrete filename
        """

        return pf_template.replace("${PBS_ARRAYID}", str(file_number))
