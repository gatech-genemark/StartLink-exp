# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20
import logging
import argparse
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_container.genome_list import GenomeInfoList, GenomeInfo
from sbsp_general import Environment
from sbsp_general.general import get_value, os_join, run_shell_cmd
from sbsp_io.general import mkdir_p
from sbsp_options.pbs import PBSOptions
from sbsp_parallelization.pbs import PBS
from sbsp_pbs_data.mergers import merge_identity
from sbsp_pbs_data.splitters import split_genome_info_list
import sbsp_argparse.parallelization
# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #

parser = argparse.ArgumentParser("Run external prediction tools on genome list.")

parser.add_argument('--pf-genome-list', required=True, help="List of genomes")
parser.add_argument('--type', required=True, choices=["archaea", "bacteria"], help="Is the list archaea or bacteria")
parser.add_argument('--dn-run', required=False, help="Name of directory that will contain the run")
parser.add_argument('--tool', choices=["gms2", "prodigal"], required=True, help="Tool used for prediction")
sbsp_argparse.parallelization.add_pbs_options(parser)

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


def run_gms2(env, gi, **kwargs):
    # type: (Environment, GenomeInfo, Dict[str, Any]) -> None

    genome_type = get_value(kwargs, "genome_type", "auto")
    pd_data = env["pd-data"]
    pd_work = env["pd-work"]
    pe_tool = os_join(env["pd-bin-external"], "gms2", "gms2.pl")

    pf_sequence = os_join(pd_data, gi.name, "sequence.fasta")

    # FIXME: put in genetic code
    run_shell_cmd(
        "cd {}; {} --gcode 11 --format gff --out gms2.gff --seq {}  --v --genome-type {} "
        "--fgio-dist-thresh 25".format(
            pd_work, pe_tool, pf_sequence, genome_type
        )
    )


def run_prodigal(env, gi, **kwargs):
    # type: (Environment, GenomeInfo, Dict[str, Any]) -> None
    pd_data = env["pd-data"]
    pd_work = env["pd-work"]
    pe_tool = os_join(env["pd-bin-external"], "prodigal", "prodigal")

    pf_sequence = os_join(pd_data, gi.name, "sequence.fasta")

    # FIXME: put in genetic code
    run_shell_cmd(
        "cd {}; {}  -i {}  -g 11  -o prodigal.gff  -f gff  -t prodigal.parameters  -q \n".format(
            pd_work, pe_tool, pf_sequence
        )
    )


def run_tool_on_gil(env, gil, tool, **kwargs):
    # type: (Environment, GenomeInfoList, str, Dict[str, Any]) -> None

    dn_run = get_value(kwargs, "dn_run", tool, default_if_none=True)
    func = {
        "gms2": run_gms2,
        "prodigal": run_prodigal,
    }[tool]

    for gi in gil:
        pd_work = os_join(env["pd-work"], gi.name, dn_run)
        mkdir_p(pd_work)
        curr_env = env.duplicate({"pd-work": pd_work})

        func(curr_env, gi, **kwargs)


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    gil = GenomeInfoList.init_from_file(args.pf_genome_list)
    pbs_options = PBSOptions.init_from_dict(env, vars(args))

    pbs = PBS(env, pbs_options,
              splitter=split_genome_info_list,
              merger=merge_identity
              )

    output = pbs.run(
        data={"gil": gil},
        func=run_tool_on_gil,
        func_kwargs={
            "env": env,
            "tool": args.tool,
            "dn_run": args.dn_run,
        }
    )


if __name__ == "__main__":
    main(my_env, parsed_args)
