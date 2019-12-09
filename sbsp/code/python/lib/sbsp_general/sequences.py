

from Bio import Seq


def get_region_around_coordinate(sequence, coord, strand, upstream_len, downstream_len, seq_type="nucl"):
    # type: (Seq, int, str, int, int) -> Seq

    if strand not in ["+", "-"]:
        raise ValueError("Strand must be:", " or ".join(["+", "-"]))

    if strand == "+":
        region_left = coord - upstream_len
        region_right = coord - 1 + downstream_len

        region = sequence[region_left:region_right+1]
    else:
        region_left = coord - downstream_len + 1
        region_right = coord + upstream_len

        region = sequence[region_left:region_right+1]
        region = region.reverse_complement()

    if seq_type == "prot":
        region = region.translate()

    return region


def count_number_of_candidate_starts_upstream(sequence, coord, strand):
    # type: (Seq, int, str) -> dict

    if strand not in ["+", "-"]:
        raise ValueError("Strand must be:", " or ".join(["+", "-"]))

    valid_starts_pos = ["ATG", "GTG", "TTG"]
    valid_starts_neq = ["CAT", "CAC", "CAA"]

    valid_stops_pos = ["TAA", "TGA", "TAG"]
    valid_stops_neg = ["TTA", "TCA", "CTA"]

    def is_valid_start(codon, valid_starts):
        return codon in valid_starts

    def is_valid_stop(codon, valid_stops):
        return codon in valid_stops

    start_counts = {
        "total" : 0,
        "atg" : 0,
        "gtg" : 0,
        "ttg" : 0,
    }

    if strand == "+":

        i = coord - 3

        while i >= 0:

            codon = sequence[i:i+3]

            if is_valid_stop(codon, valid_stops_pos):
                break

            if is_valid_start(codon, valid_starts_pos):
                start_counts[codon.lower()] += 1
                start_counts["total"] += 1

            i -= 3
    else:

        i = coord + 3

        while i < len(sequence):

            codon = sequence[i-2:i+1]       # type: Seq

            if is_valid_stop(codon, valid_stops_neg):
                break

            if is_valid_start(codon, valid_starts_neq):
                codon_rc = codon.reverse_complement()

                start_counts[codon_rc.lower()] += 1
                start_counts["total"] += 1

            i += 3

    return start_counts



