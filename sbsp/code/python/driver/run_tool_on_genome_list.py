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

def create_pbs_file(env, cmd_run, pf_pbs, **kwargs):
        
        job_name = get_value(kwargs, "job_name", "JOB")
        num_nodes = get_value(kwargs, "num_nodes", 1)
        ppn = get_value(kwargs, "ppn", 1)
        node_property = get_value(kwargs, "node_property", "")
        walltime = get_value(kwargs, "pbs-walltime", "07:00:00")

        pd_work = env["pd-work"]

        pbs_text = ""

        pbs_text += "#PBS -N " + str(job_name) + "\n"
        pbs_text += "#PBS -o " + "{}/{}".format(pd_work, "error") + "\n"
        pbs_text += "#PBS -j oe" + "\n"
        pbs_text += "#PBS -l nodes=" + str(num_nodes) + ":ppn=" + str(ppn) + "{}\n".format(node_property)
        pbs_text += "#PBS -l walltime=" + str(walltime) + "\n"

        pbs_text += "#PBS -W umask=002" + "\n"

        pbs_text += "export PATH=\"/home/karl/anaconda/envs/sbsp/bin:$PATH\"\n"

        pbs_text += "PBS_O_WORKDIR=" + pd_work + "\n"
        pbs_text += "cd $PBS_O_WORKDIR \n"

        pbs_text += "echo The working directory is `echo $PBS_O_WORKDIR`" + "\n"
        pbs_text += "echo This job runs on the following nodes:" + "\n"
        pbs_text += "echo `cat $PBS_NODEFILE`" + "\n"

        pbs_text += "\n{}\n".format(cmd_run)

        from sbsp_io.general import write_string_to_file
        write_string_to_file(pbs_text, pf_pbs)


def run_gms2(env, gi, **kwargs):
    # type: (Environment, GenomeInfo, Dict[str, Any]) -> None

    genome_type = get_value(kwargs, "genome_type", "auto")
    pd_data = env["pd-data"]
    pd_work = env["pd-work"]
    pe_tool = os_join(env["pd-bin-external"], "gms2", "gms2.pl")

    pf_sequence = os_join(pd_data, gi.name, "sequence.fasta")

    # FIXME: put in genetic code
    cmd_run = "{} --gcode 11 --format gff --out gms2.gff --seq {}  --v --genome-type {} --fgio-dist-thresh 25".format(
            pe_tool, pf_sequence, genome_type
        )

    pf_pbs = os_join(pd_work, "run.pbs")
    create_pbs_file(env, cmd_run, pf_pbs, job_name=gi.name, **kwargs)

    run_shell_cmd("qsub {} &".format(pf_pbs))


def run_prodigal(env, gi, **kwargs):
    # type: (Environment, GenomeInfo, Dict[str, Any]) -> None
    pd_data = env["pd-data"]
    pd_work = env["pd-work"]
    pe_tool = os_join(env["pd-bin-external"], "prodigal", "prodigal")

    pf_sequence = os_join(pd_data, gi.name, "sequence.fasta")

    # FIXME: put in genetic code
    cmd_run ="{}  -i {}  -g 11  -o prodigal.gff  -f gff  -t prodigal.parameters  -q \n".format(
            pe_tool, pf_sequence
        )
    pf_pbs = os_join(pd_work, "run.pbs")
    create_pbs_file(env, cmd_run, pf_pbs, job_name=gi.name, **kwargs)

    run_shell_cmd("qsub {} &".format(pf_pbs))


def run_tool_on_gil(env, gil, tool, **kwargs):
    # type: (Environment, GenomeInfoList, str, Dict[str, Any]) -> None

    logger.info("Running tool {} on {} genomes".format(tool, len(gil)))
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
    
    run_tool_on_gil(env, gil, args.tool, dn_run=args.dn_run, genome_type=args.type)


if __name__ == "__main__":
    main(my_env, parsed_args)
