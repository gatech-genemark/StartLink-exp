from __future__ import absolute_import


from . import labels
import pandas as pd

import sbsp_io.sequences
import sbsp_general.sequences


def create_key(genome_name, accession, left, right, strand):
    return "{};{};{};{};{}".format(genome_name, accession, left, right, strand)


def create_key_for_label(genome_name, label):
    return create_key(genome_name, label.seqname(), label.coordinates().left,
                           label.coordinates().right, label.coordinates().strand)


def get_label_for_row(df_row, field_prefix, field_suffix):

    col_left = "{}-{}".format(field_prefix, "left")
    col_right = "{}-{}".format(field_prefix, "right" )
    col_strand = "{}-{}".format(field_prefix, "strand")
    col_accession = "{}-{}".format(field_prefix, "accession")

    if field_suffix:
        col_left = "{}-{}".format(col_left, field_suffix)
        col_right = "{}-{}".format(col_right, field_suffix)
        col_strand = "{}-{}".format(col_strand, field_suffix)

    label = labels.Label.from_fields({
        "left": int(df_row[col_left]) - 1,
        "right": int(df_row[col_right]) - 1,
        "strand": df_row[col_strand],
        "seqname": df_row[col_accession],
    })

    return label

def get_labels_per_genome(df, source, **kwargs):
    # type: (pd.DataFrame, str, **str) -> Dict[str, labels.Labels]

    source_prefix_pair = {"query": "q", "target": "t"}

    if source not in source_prefix_pair:
        raise ValueError("Source should be", ", ".join(source_prefix_pair))

    prefix = source_prefix_pair[source]

    suffix = ""
    if "coordinates_suffix" in kwargs and kwargs["coordinates_suffix"] is not None:
        suffix = ("-" + str(kwargs["coordinates_suffix"]))



    col_left = "{}-{}".format(prefix, "left" + suffix)
    col_right = "{}-{}".format(prefix, "right" + suffix)
    col_strand = "{}-{}".format(prefix, "strand" + suffix)
    col_genome = "{}-{}".format(prefix, "genome")
    col_accession = "{}-{}".format(prefix, "accession")

    labels_per_genome = dict()



    for index, row in df.iterrows():
        label = labels.Label.from_fields({
            "left": int(row[col_left])-1,
            "right": int(row[col_right])-1,
            "strand": row[col_strand],
            "seqname": row[col_accession],
        })

        genome = row[col_genome]

        if genome not in labels_per_genome:
            labels_per_genome[genome] = labels.Labels()

        labels_per_genome[genome].add(label)

    return labels_per_genome


def get_regions_per_key(labels_per_genome, p_dir_data, upstream_len, downstream_len, seq_type="nucl"):
    # type: (Dict[str, labels.Labels], str, int, int, str) -> Dict[str, str]

    if seq_type == "prot":
        if (upstream_len + downstream_len) % 3 != 0:
            raise ValueError("Region length should be a multiple of three when converting to protein sequences")

    key_to_sequence = dict()

    # do it genome after genome, so as to only read each sequence file once
    for genome_name  in labels_per_genome:

        # read in sequence file
        p_seq = p_dir_data + "/" + genome_name + "/" + "sequence.fasta"

        genome_seqs = sbsp_io.sequences.read_fasta_into_hash(p_seq)

        # for each label for that genome
        for label in labels_per_genome[genome_name]:

            left = label.coordinates().left
            right = label.coordinates().right
            strand = label.coordinates().strand
            seqname = label.seqname()

            frag = sbsp_general.sequences.get_region_around_coordinate(genome_seqs[seqname], label.get_5prime(), strand,
                                                                upstream_len, downstream_len, seq_type)

            key = create_key(genome_name, seqname, left, right, strand)

            key_to_sequence[key] = frag

    return key_to_sequence


def parse_data(df, features, ground_truth_col_name):
    def get_values_for_features(df, features, ground_truth_col_name):
        X = df.loc[:, df.columns != ground_truth_col_name]

        if len(features) > 0:
            X = X[features]

        return X

    def get_ground_truth(df, ground_truth_col_name):
        y = None
        if ground_truth_col_name in df:
            y = df[ground_truth_col_name]
            y = y.replace(1, 1)
            y = y.replace(0, -1)

        return y

    return [get_values_for_features(df, features, ground_truth_col_name),
            get_ground_truth(df, ground_truth_col_name)]


def read_data(pf_data, features, ground_truth_col_name):

    df = pd.read_csv(pf_data, header=0, delimiter=",")
    if features is None:
        features = list(df.columns.values)

    [X, y] = parse_data(df, features, ground_truth_col_name)

    return [df, X, y]



def sort_clustal_by_tag(tag_func_):
    pass







