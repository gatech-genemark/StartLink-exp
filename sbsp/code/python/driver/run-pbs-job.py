# Karl Gemayel
# Georgia Institute of Technology
#
# Created:

import logging
import argparse
from typing import *

# noinspection PyUnresolvedReferences
import pathmagic                        # add path to custom library

# Custom library imports
from sbsp_general import Environment
import sbsp_general
from sbsp_parallelization.pbs import PBSJobPackage

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-job-input', required=True)
parser.add_argument('--pf-job-output', required=True)


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
logger = logging.getLogger("logger")                    # type: logging.Logger


def main(env, args):
    # type: (Dict[str, Any], argparse.Namespace) -> None

    pbs_package = PBSJobPackage.load(args.pf_job_input)
    func = pbs_package["func"]
    func_args = pbs_package["func_kwargs"]
    logger.critical("{}\n{}".format(func, func_args))
    output = {
        "data": func(**func_args)
    }

    PBSJobPackage.save(output, args.pf_job_output)


if __name__ == "__main__":
    main(my_env, parsed_args)
