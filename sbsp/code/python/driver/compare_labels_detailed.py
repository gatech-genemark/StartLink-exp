# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 12/19/19

import logging
import argparse
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_general import Environment

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_general.labels import Labels
from sbsp_general.labels_comparison_detailed import LabelsComparisonDetailed
from sbsp_io.labels import read_labels_from_file
from sbsp_viz.labels_comparison_detailed import LabelsComparisonDetailedViz

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-a', required=True, help="Path to file containing labels")
parser.add_argument('--pf-b', required=True, help="Path to file containing labels")

parser.add_argument('--name-a', required=False, default="A", help="Name of first set")
parser.add_argument('--name-b', required=False, default="B", help="Name of second set")

parser.add_argument('--split-on-attributes', default=None, nargs="+", help="If set, labels also split by attribute (GFF)")

parser.add_argument('--tag', required=False, help="Name of genome/species for label set")

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

    labels_a = read_labels_from_file(args.pf_a)
    labels_b = read_labels_from_file(args.pf_b)

    lcd = LabelsComparisonDetailed(labels_a, labels_b,
                             name_a=args.name_a,
                             name_b=args.name_b,
                             tag=args.tag,
                                   split_on_attributes=args.split_on_attributes)

    LabelsComparisonDetailedViz(lcd).run(env["pd-work"])


if __name__ == "__main__":
    main(my_env, parsed_args)
