# Karl Gemayel
# Georgia Institute of Technology
#
# Modified: 06/19/2020

import random
import logging
import argparse
import pandas as pd
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
import sbsp_argparse.sbsp
import sbsp_argparse.parallelization
from sbsp_general import Environment
from sbsp_io.general import mkdir_p
from sbsp_options.sbsp import SBSPOptions
from sbsp_general.general import os_join, get_value
from sbsp_parallelization.generic_threading import run_one_per_thread
from sbsp_pipeline.pipeline_msa import PipelineSBSP
from sbsp_options.pipeline_sbsp import PipelineSBSPOptions
from sbsp_options.parallelization import ParallelizationOptions
from sbsp_container.genome_list import GenomeInfoList, GenomeInfo

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #


parser = argparse.ArgumentParser("Run SBSP on a list of genomes.")

parser.add_argument('--pf-q-list', required=True, help="File of query genomes")
parser.add_argument('--pf-db-index', required=True,
                    help="Path to file containing database information for each ancestor clade")

parser.add_argument('--simultaneous-genomes', type=int, default=1, help="Number of genomes to run on simultaneously.")
parser.add_argument('--dn-run', default="sbsp", help="Name of directory with SBSP run.")

parser.add_argument('--fn-q-labels', default="ncbi.gff", required=False, type=Union[str],
                    help="Name of query file(s) containing gene labels")
parser.add_argument('--fn-t-labels', default="ncbi.gff", required=False, type=Union[str],
                    help="Name of target file(s) containing gene labels")

parser.add_argument("--fn-q-labels-compare", default="ncbi.gff", required=False, type=Union[str],
                    help="Name of true labels file. If set, accuracy is computed after MSA.")

parser.add_argument('--steps', nargs="+", required=False,
                    choices=["prediction", "comparison"],
                    default=None)

sbsp_argparse.parallelization.add_parallelization_options(parser)
sbsp_argparse.sbsp.add_sbsp_options(parser)

parser.add_argument('--pd-work', required=False, default=None, help="Path to working directory")
parser.add_argument('--pd-data', required=False, default=None, help="Path to data directory")
parser.add_argument('--pd-results', required=False, default=None, help="Path to results directory")
parser.add_argument("-l", "--log", dest="loglevel", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                    help="Set the logging level", default='WARNING')

parsed_args = parser.parse_args()

# ------------------------------ #
#           Main Code            #
# ------------------------------ #

# Load environment variables
my_env = Environment(pd_data=parsed_args.pd_data,
                     pd_work=parsed_args.pd_work,
                     pd_results=parsed_args.pd_results)

# Setup logger
logging.basicConfig(level=parsed_args.loglevel)
logger = logging.getLogger("logger")  # type: logging.Logger


def sbsp_on_gi(gi, pipeline_options):
    # type: (GenomeInfo, PipelineSBSPOptions) -> None
    random.seed(1)
    logger.debug("Running for {}".format(gi.name))
    PipelineSBSP(pipeline_options.env, pipeline_options).run()
    logger.debug("Done for {}".format(gi.name))


def get_clade_to_pf_db(pf_db_index):
    # type: (str) -> Dict[str, str]
    df = pd.read_csv(pf_db_index)
    return {
        r["Clade"]: r["pf-db"] for _, r in df.iterrows()
    }


def setup_gi_and_run(env, gi, sbsp_options, prl_options, clade_to_pf_db, **kwargs):
    # type: (Environment, GenomeInfo, SBSPOptions, ParallelizationOptions, Dict[str, str], Dict[str, Any]) -> None

    dn_run = get_value(kwargs, "dn_run", "sbsp")

    # Check if clade is known
    try:
        pf_t_db = clade_to_pf_db[gi.attributes["ancestor"]]
    except KeyError:
        raise ValueError("Unknown clade {}".format(gi.attributes["ancestor"]))

    logger.info("Scheduling: {}".format(gi.name))

    pd_work = os_join(env["pd-work"], gi.name, dn_run)  # genome working environment
    curr_env = env.duplicate({"pd-work": pd_work})  # create environment for genome
    pf_output = os_join(pd_work, "output.csv")  # output file

    mkdir_p(pd_work)  # create working directory

    # write genome name to file list (for running)
    pf_list = os_join(pd_work, "query.list")
    GenomeInfoList([gi]).to_file(pf_list)

    # create options for pipeline for current genome
    po = PipelineSBSPOptions(
        curr_env, pf_list, pf_t_db=pf_t_db, pf_output=pf_output, sbsp_options=sbsp_options,
        prl_options=prl_options, **kwargs
    )
    sbsp_on_gi(gi, po)


def run_sbsp_on_genome_list(env, gil, sbsp_options, prl_options, clade_to_pf_db, **kwargs):
    # type: (Environment, GenomeInfoList, SBSPOptions, ParallelizationOptions, Dict[str, str], Dict[str, str]) -> None
    """
    Runs SBSP on list of genomes using specified options.
    :param env: General environment
    :param gil: list of genomes
    :param sbsp_options: Options for controlling algorithm behavior
    :param prl_options: Options for controlling parallelization of runs
    :param clade_to_pf_db: map of clade to file containing target database
    :param kwargs: Optional arguments:
        simultaneous_genomes: Number of genomes to run simultaneously
        dn_run: Name of directory in which to put run
    :return: None
    """

    simultaneous_genomes = get_value(kwargs, "simultaneous_genomes", 1, default_if_none=True)
    dn_run = get_value(kwargs, "dn_run", "sbsp")

    run_one_per_thread(
        gil, setup_gi_and_run, data_arg_name="gi",
        func_kwargs={
            "env": env, "sbsp_options": sbsp_options, "prl_options": prl_options, "clade_to_pf_db": clade_to_pf_db,
            **kwargs,
        },
        simultaneous_runs=simultaneous_genomes
    )


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    gil = GenomeInfoList.init_from_file(args.pf_q_list)  # read genome list
    sbsp_options = SBSPOptions.init_from_dict(env, vars(args))  # read tool options
    prl_options = ParallelizationOptions.init_from_dict(env, vars(args))  # read parallelization options

    # read database index file: shows locations of databases used as targets for each clade
    clade_to_pf_db = get_clade_to_pf_db(args.pf_db_index)

    run_sbsp_on_genome_list(
        env, gil, sbsp_options, prl_options, clade_to_pf_db,
        simultaneous_genomes=args.simultaneous_genomes,
        dn_run=args.dn_run,
        steps=args.steps,
        fn_q_labels=args.fn_q_labels,
        fn_t_labels=args.fn_t_labels,
        fn_q_labels_compare=args.fn_q_labels_compare
    )


if __name__ == "__main__":
    main(my_env, parsed_args)
