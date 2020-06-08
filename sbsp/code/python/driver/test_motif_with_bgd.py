# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20
import logging
import argparse
from typing import *
import logomaker as lm
import matplotlib.pyplot as plt

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_container.gms2_mod import GMS2Mod
from sbsp_general import Environment

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_general.GMS2Noncoding import GMS2Noncoding
from sbsp_general.MotifModel import MotifModel

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-gms2-mod', required=True)

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

    mod = GMS2Mod.init_from_file(args.pf_gms2_mod)

    motif_mod = MotifModel(mod.items["RBS_MAT"], None)
    nonc_mod = GMS2Noncoding(mod.items["NON_MAT"])

    df_motif = motif_mod.pwm_to_df()
    bgd = nonc_mod.pwm_to_array(0)

    print(df_motif.to_csv())
    print(bgd)


    fig, axes = plt.subplots(1, 2, sharex="all", sharey="all", figsize=(8,4))

    # relative
    rel_mat = lm.transform_matrix(df_motif, from_type="probability", to_type="information", background=bgd)
    lm.Logo(rel_mat, color_scheme="classic", ax=axes[0])
    axes[0].set_ylim(*[0, 2])

    # shannon
    sha_mat = lm.transform_matrix(df_motif, from_type="probability", to_type="information")
    lm.Logo(sha_mat, color_scheme="classic", ax=axes[1])
    axes[1].set_ylim(*[0, 2])
    plt.show()






if __name__ == "__main__":
    main(my_env, parsed_args)
