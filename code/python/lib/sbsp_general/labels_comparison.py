from typing import *
import sbsp_general.general
from sbsp_general import Environment


class LabelsComparison:

    def __init__(self, env, pf_labels_a, pf_labels_b):
        # type: (Environment, str, str) -> None

        self.pf_labels_a = pf_labels_a
        self.pf_labels_b = pf_labels_b

        cmd = "compp -a {} -b {} -v".format(
            pf_labels_a,
            pf_labels_b
        )

        compp_output = sbsp_general.general.run_shell_cmd(cmd)

        self._parse_compp_output(compp_output)

    def _parse_compp_output(self, output):
        # type: (str) -> None
        d = self._split_compp_output(output)

        # put into reasonably-named keys
        good_and_bad = [
            ("found", "found_in_A"),
            ("identical", "identical_in_A")
        ]

        for letter in [("a", "A"), ("b", "B")]:
            g = letter[0]       # good
            b = letter[1]       # bad

            good_and_bad += [
                ("in-{}".format(g), "in_{}".format(b)),
                ("long-in-{}".format(g), "long_in_{}".format(b)),
                ("short-in-{}".format(g), "short_in_{}".format(b)),
                ("unique-in-{}".format(g), "unique_in_{}".format(b)),
            ]

        self.stats = {goodname: d[badname] for goodname, badname in good_and_bad}

    @staticmethod
    def _split_compp_output(output):
        # type: (str) -> Dict[str, Any]

        lines = output.strip().split("\n")

        import re

        p = re.compile(r"^\s*([^:]+):\s*(\d+)$")

        d = {}

        for line in lines:
            x = p.match(line)
            if x:
                key = x.group(1)
                val = x.group(2)

                d[key] = val

        return d


    @staticmethod
    def stringify_genome_accuracies(genome_to_comparison, delimiter="\t"):
        # type: (Dict[str, LabelsComparison], str) -> str

        out = "Genome{}Total{}Attempted{}Correct".format(delimiter,delimiter,delimiter)

        for genome in sorted(genome_to_comparison.keys()):

            out += "\n{}{}{}{}{}{}{}".format(
                genome, delimiter,
                genome_to_comparison[genome].stats["in-a"], delimiter,
                genome_to_comparison[genome].stats["found"], delimiter,
                genome_to_comparison[genome].stats["identical"]
            )

        out += "\n"

        return out



