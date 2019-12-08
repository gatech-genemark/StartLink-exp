import os

from sbsp_general import Environment
from sbsp_options.options import Options


class MSAOptions(Options):
    """Options for multiple sequence alignment"""

    @staticmethod
    def path_to_default_options_file(env):
        # type: (Environment) -> str
        return os.path.join(env["pd-config"], "msa-defaults.conf")
