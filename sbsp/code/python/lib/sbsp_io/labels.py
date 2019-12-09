import os
import re
import sys

sys.path.append(os.path.dirname(__file__) + "/..")       # add custom library directory to path

import sbsp_general.labels

def write_string_to_file(astring, filename):

    with open(filename, "w") as f:
        f.write(astring)


def write_labels_to_file(labels, filename):

    out = labels.to_string(shift_coordinates_by=1)
    write_string_to_file(out, filename)


def create_attribute_dict(attribute_string, delimiter=";"):
    # type: (str, str) -> Dict[str, Any]

    attributes = dict()

    for current in attribute_string.strip().split(sep=delimiter):
        if len(current.strip().split("=")) == 2:
            k, v = current.strip().split("=")
            attributes[k] = v

    return attributes

def read_labels_from_file(filename, shift=-1):
    # type: (str) -> sbsp_general.labels.Labels
    # FIXME: only supports gff

    labels = sbsp_general.labels.Labels()

    pattern = re.compile(r"([^\t]+)\t([^\t]+)\t(CDS)\t(\d+)\t(\d+)\t([^\t]+)\t([+-])\t([^\t]+)\t([^\t]+)")

    with open(filename, "r") as f:

        for line in f:

            line = line.strip()

            m = pattern.match(line)
            if m:

                attributes = create_attribute_dict(m.group(9))

                label = sbsp_general.labels.Label.from_fields(
                    {
                        "left" : int(m.group(4)) + shift,
                        "right" : int(m.group(5)) + shift,
                        "strand" : m.group(7),
                        "seqname" : m.group(1),
                    },
                    attributes=attributes
                )

                labels.add(label)
    return labels



