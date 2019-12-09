import os

from Bio import AlignIO

import sbsp_general
import sbsp_general.dataframe
from sbsp_general.general import get_value
from sbsp_general.labels import Label
from sbsp_general.msa_2 import MSAType, MSASinglePointMarker




def create_key_3prime_from_label(label, genome_name=None):
    # type: (Label, str) -> str
    if label.strand() == "+":
        return "{};{};{};{};{}".format(genome_name, label.seqname(), "",
                                       label.coordinates().right, label.strand())
    else:
        return "{};{};{};{};{}".format(genome_name, label.seqname(), label.coordinates().left,
                                       "", label.strand())



def add_gene_labels_from_file(env, df, **kwargs):
    # type: (Dict[str, Any], pd.DataFrame, Dict[str, Any]) -> None

    fn_q_labels = get_value(kwargs, "fn_q_labels", "verified.gff")
    source = get_value(kwargs, "source", "q")
    suffix_coordinates = get_value(kwargs, "suffix_corrdinates", "ref")

    from sbsp_io.labels import read_labels_from_file
    import sbsp_general.dataframe

    all_genomes = set(df["{}-genome".format(source)])
    genome_to_genekey_to_label = dict()

    for genome_name in all_genomes:
        pf_q_labels = os.path.join(env["pd-data"], genome_name, fn_q_labels)

        labels = read_labels_from_file(pf_q_labels)

        key_3prime_to_label = dict()
        for l in labels:
            key_3prime = create_key_3prime_from_label(l)
            key_3prime_to_label[key_3prime] = l

        genome_to_genekey_to_label[genome_name] = key_3prime_to_label


    # now add to data frame
    column_left = "{}-left-{}".format(source, suffix_coordinates)
    column_right = "{}-right-{}".format(source, suffix_coordinates)
    column_strand = "{}-strand-{}".format(source, suffix_coordinates)

    df[column_left] = -1
    df[column_right] = -1
    df[column_strand] = ""

    for index, row in df.iterrows():

        curr_genome = row["{}-genome".format(source)]
        curr_label = sbsp_general.dataframe.df_get_label_from_row(df, index, source)
        curr_key = create_key_3prime_from_label(curr_label)

        if curr_key in genome_to_genekey_to_label[curr_genome].keys():

            sbsp_general.dataframe.df_coordinates_to_row(df, index, curr_label, source, suffix_coordinates=suffix_coordinates)


def get_reference_position_in_msa(msa_t, df, **kwargs):
    # type: (MSAType, pd.DataFrame, Dict[str, Any]) -> Union[int, None]

    suffix_msa_coordinates = get_value(kwargs, "suffix_msa_coordinates", "msa")
    suffix_ref_coordinates = get_value(kwargs, "suffix_ref_coordinates", "ref")
    msa_nt = get_value(kwargs, "msa_nt", False)
    fn_q_labels_true = get_value(kwargs, "fn_q_labels_true", "verified.gff")

    idx_0 = df.index.values[0]
    msa_label = sbsp_general.dataframe.df_get_label_from_row(df, idx_0, "q",
                                                             suffix_coordinates=suffix_msa_coordinates)
    ref_label = sbsp_general.dataframe.df_get_label_from_row(df, idx_0, "q",
                                                             suffix_coordinates=suffix_ref_coordinates)

    msa_5prime = msa_label.get_5prime()
    ref_5prime = ref_label.get_5prime()

    if ref_5prime < 0:
        return None


    # distance between msa and ref; negative means ref is upstream
    if msa_label.strand() == "+":
        distance = ref_5prime - msa_5prime
    else:
        distance = msa_5prime - ref_5prime

    if distance % 3 != 0:
        raise ValueError("Distance between MSA and Ref must be a multiple of 3")

    distance_aa = int(distance / 3)
    if msa_nt:
        distance_aa *= 3

    msa_pos_in_frag = msa_t.get_mark_position("selected")
    ref_pos_in_frag = msa_pos_in_frag

    if ref_pos_in_frag is None:
        return None

    if distance_aa > 0:
        # add that many non-gapped
        for r in range(distance_aa):
            ref_pos_in_frag += 1
            if ref_pos_in_frag > msa_t.alignment_length():
                ref_pos_in_frag = None
                break

            while msa_t[0].seq._data[ref_pos_in_frag] == "-":
                ref_pos_in_frag += 1

                if ref_pos_in_frag > msa_t.alignment_length():
                    ref_pos_in_frag = None
                    break

            if ref_pos_in_frag is None:
                return None

    elif distance_aa < 0:

        for r in range(abs(distance_aa)):
            ref_pos_in_frag -= 1
            if ref_pos_in_frag < 0:
                ref_pos_in_frag = None
                break

            while msa_t[0].seq._data[ref_pos_in_frag] == "-":
                ref_pos_in_frag -= 1

                if ref_pos_in_frag < 0:
                    ref_pos_in_frag = None
                    break

            if ref_pos_in_frag is None:
                return None

    return ref_pos_in_frag


def add_true_starts_to_msa_output(env, df, **kwargs):
    # type: (Dict[str, Any], pd.DataFrame, Dict[str, Any]) -> None

    msa_nt = get_value(kwargs, "msa_nt", False)
    fn_q_labels_true = get_value(kwargs, "fn_q_labels_true", "verified.gff")

    add_gene_labels_from_file(env, df, fn_q_labels=fn_q_labels_true)

    column_pf_msa_output = "pf-msa-output"

    for pf_msa_output, df_group in df.groupby(column_pf_msa_output):
        if msa_nt:
            pf_msa_output += "_nt"

        msa_t = MSAType.init_from_file(pf_msa_output)

        ref_position_in_msa = get_reference_position_in_msa(msa_t, df_group, **kwargs)

        marker = MSASinglePointMarker(ref_position_in_msa, msa_t.alignment_length(), name="ref")

        msa_t.add_marker(marker, unique=True)

        msa_t.to_file(pf_msa_output)


def update_msa_outputs(df, **kwargs):
    # type: (pd.DataFrame, Dict[str, Any]) -> None

    column_q_upstream_distance = get_value(kwargs, "column_q_upstream_distance", None)
    column_t_upstream_distance = get_value(kwargs, "column_t_upstream_distance", None)

    column_pf_msa_output = "pf-msa-output"
    column_t_pos_in_msa_output = "t-pos-in-msa-output"

    groups = df.groupby(column_pf_msa_output)

    for pf_msa_output, df_group in groups:

        msa_t = MSAType.init_from_file(pf_msa_output)

        # for each alignment, update id with upstream distance

        for index, row in df_group.iterrows():

            pos_in_alignments = int(row[column_t_pos_in_msa_output])

            msa_t[pos_in_alignments].id = msa_t[pos_in_alignments].id.split("-")[0]       # get id

            # add distance to upstream
            if "t-upstream-pos-in-frag-msa" in row:
                msa_t[pos_in_alignments].id += ";{}".format(row["t-upstream-pos-in-frag-msa"])

            # add all distances
            for dt in ["kimura", "kimura3", "ds", "dn"]:
                if dt in row:
                    msa_t[pos_in_alignments].id += ";{}".format(round(float(row[dt]), 2))

            if column_t_upstream_distance is not None and column_t_upstream_distance in row:
                dist_upstream = row[column_t_upstream_distance]

                msa_t[pos_in_alignments].id += ";{}".format(int(dist_upstream))

            gene_length = row["t-right-msa"] - row["t-left-msa"] + 1
            msa_t[pos_in_alignments].id += ";{}".format(int(gene_length))

        # add to query
        msa_t[0].id = msa_t[0].id.split("-")[0]
        if "q-upstream-pos-in-frag-msa" in df_group:
            msa_t[0].id += ";{}".format(df_group["q-upstream-pos-in-frag-msa"].iloc[0])

        if column_q_upstream_distance is not None and column_q_upstream_distance in df_group.columns.values:
            msa_t[0].id += ";{}".format(int(df_group[column_q_upstream_distance].iloc[0]))


        msa_t[0].id += ";{}".format(int(df_group["q-right-msa"].iloc[0])-int(df_group["q-left-msa"].iloc[0])+1)

        if "q-upstream-pos-in-frag-msa" in df_group:
            msa_t.add_marker(
                # MSASinglePointMarker(df_group["q-upstream-pos-in-frag-msa"].iloc, msa_t.alignment_length(),
                #                      name="q3prime", mark="*")
                MSASinglePointMarker(
                    df_group["q-upstream-pos-in-frag-msa"].iloc[0],
                    msa_t.alignment_length(),
                    name="q3prime", mark="*")
            )

        msa_t.to_file(pf_msa_output)

