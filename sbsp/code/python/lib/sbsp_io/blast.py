from Bio.Blast import NCBIXML


def read_hits(fname):
    # type: (str) -> list[NCBIXML.Record.Blast]

    fresult = open(fname, "r")
    records = NCBIXML.parse(fresult)

    return records


def split_blast_output_by_query(pf_hits, tag, num_splits):
    # type: (str, str, int) -> list[str]

    def count_records(pf_hits):
        # type: (str) -> int
        records = read_hits(pf_hits)
        num_records = 0
        for _ in records:
            num_records += 1

        return num_records

    # count number of queries
    num_queries = count_records(pf_hits)

    # queries per split
    import math
    num_queries_per_split = math.ceil(num_queries / float(num_splits))

    from sbsp_io.xml_breaker import XMLBreaker

    return XMLBreaker.run_breaker(pf_hits, tag, num_queries_per_split)







