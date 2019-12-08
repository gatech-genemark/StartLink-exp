import os

from typing import *
from sbsp_general import Environment
from sbsp_options.options import Options


class PBSOptions(Options):
    """Options for PBS scheduler"""

    @staticmethod
    def path_to_default_options_file(env):
        # type: (Environment) -> str
        return os.path.join(env["pd-config"], "pbs-defaults.conf")

    @staticmethod
    def required():
        # type: () -> Union[Set[str], None]
        return {
            # directories for pbs
            "pd-root-compute", "dn-compute", "pd-head",
            # node computation configuration
            "num-processors", "num-nodes", "walltime",
            "num-jobs",
            "max-concurrent-nodes",
            "split-tag",
            "group-by",
            "use-pbs"
        }
