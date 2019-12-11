

def read_assembly_summary(fname):
    # type: (str) -> dict

    # Output structure:
    # {
    #    "column_names": [column1, column2, ..., columnN],
    #    "data": [
    #        {"column1": value1, "column2": value2, ..., "columnN": valueN},
    #        {"column1": value1, "column2": value2, ..., "columnN": valueN},
    #        ...
    #    ]
    # }

    with open(fname, "r") as f:

        columns = list()
        data = list()

        for line in f:

            line = line.rstrip("\n")
            line = line.lstrip("\t")

            if len(line) == 0:
                continue

            if line.strip()[0] == "#":
                # check for header
                if "assembly_accession" in line:
                    columns = line[1:].strip().split("\t")      # skip hash and then split into column names
            else:
                # data line
                data_arr = line.split("\t")

                name_data_pairs = dict()

                for pair in zip(columns, data_arr):
                    name_data_pairs[pair[0]] = pair[1]

                data += [name_data_pairs]

        return {"column_names": columns, "data": data}


def write_assembly_summary(summary, fname):
    # type: (dict, str) -> None

    # Output structure:
    # {
    #    "column_names": [column1, column2, ..., columnN],
    #    "data": [
    #        {"column1": value1, "column2": value2, ..., "columnN": valueN},
    #        {"column1": value1, "column2": value2, ..., "columnN": valueN},
    #        ...
    #    ]
    # }

    with open (fname, "w") as f:

        column_names = summary["column_names"]
        data = summary["data"]

        # create header line
        out = "#\t" + "\t".join(column_names)

        f.write(out + "\n")

        for d in data:

            out = ""
            prefix = ""
            for name in column_names:

                out += prefix + str(d[name])
                prefix = "\t"

            f.write(out + "\n")


def get_rows_by_key(pf_assembly_summary, key):
    # type: (str, str) -> Dict[int, List[Dict[str, Any]]]

    assembly_info = read_assembly_summary(pf_assembly_summary)

    result = dict()

    for d in assembly_info["data"]:
        taxid = int(d["taxid"])

        if d["assembly_level"] not in {"Complete Genome", "Scaffold", "Contig"}:
            continue

        if taxid not in result.keys():
            result[taxid] = list()
        result[taxid].append(d)

    return result