# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 2/16/20

import logging
import argparse
import os
from os import remove
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_argparse.parallelization
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_alg.ortholog_finder import extract_labeled_sequences_for_genome, extract_labeled_sequences_for_genomes
from sbsp_alg.sbsp_steps import duplicate_parallelization_options_with_updated_paths
from sbsp_container.genome_list import GenomeInfoList, GenomeInfo
from sbsp_general import Environment

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_io.sequences import read_fasta_into_hash
from sbsp_options.parallelization import ParallelizationOptions
from sbsp_parallelization.pbs import PBS
from sbsp_pbs_data.mergers import merge_identity
from sbsp_pbs_data.splitters import split_genome_info_list

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-genome-list', required=True, help="Genome list")
parser.add_argument('--pf-output', required=True, help="Output file")
sbsp_argparse.parallelization.add_parallelization_options(parser)


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


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    logger.debug("Running: sbsp_step_get_orthologs")

    prl_options = ParallelizationOptions.init_from_dict(env, vars(args))
    prl_options = duplicate_parallelization_options_with_updated_paths(env, prl_options)

    if prl_options.safe_get("pd-data-compute"):
        env = env.duplicate({"pd-data": prl_options["pd-data-compute"]})

    pbs = PBS(env, prl_options,
              splitter=split_genome_info_list,
              merger=merge_identity
              )

    gil = GenomeInfoList.init_from_file(args.pf_genome_list)

    output = pbs.run(
        data={"gil": gil,
              "pf_output_template": os.path.join(prl_options["pbs-pd-head"], "sequences_{}.faa")},
        func=extract_labeled_sequences_for_genomes,
        func_kwargs={
            "env": env,
            "fn_labels": "ncbi.gff",
            "reverse_complement": True,
            "ignore_frameshifted": True,
            "ignore_partial": True
        }
    )

    with open(args.pf_output, "w") as f_output:

        for pf_tmp in output:
            if pf_tmp is None or not os.path.isfile(pf_tmp):
                continue

            sequences = read_fasta_into_hash(pf_tmp, stop_at_first_space=False)

            for k, v in sequences.items():
                f_output.write(">{}\n{}\n".format(k, v))

            remove(pf_tmp)

        f_output.close()



if __name__ == "__main__":
    main(my_env, parsed_args)
