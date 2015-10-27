#!/usr/bin/env python
# Filter out UTRs and, if a gene has multiple UTRs, return a single UTR with minimum start and maximum end coordinate
# usage: python filter_utr.py input.gtf output.gtf

import sys
from collections import OrderedDict

def get_gtf_fields():
    return [ "chr", "source", "feature", "start", "end", "score", "strand", "frame", "group" ]

def get_feature_dict(line, gtf_fields, utr_dict, feature="UTR"):
    """
    Return a dictionary with lines of a GTF if the line describes a UTR.
    Key is the first attribute of the group.
    """
    if line.split("\t")[2] == feature:
        fields = line.strip().split("\t")
        fields[3] = int(fields[3])
        fields[4] = int(fields[4])
        gene_id = fields[-1].split("\"")[1]
        if not gene_id in utr_dict:
            utr_dict[gene_id] = [OrderedDict(zip(gtf_fields, fields))]
        else:
            utr_dict[gene_id].append(OrderedDict(zip(gtf_fields, fields)))

def remove_five_prime_utrs(utr_dict, start_codon_dict):
    for gene, features in start_codon_dict.iteritems():
        is_reverse = features[0]["strand"] == "-"
        if is_reverse:
            stop_codon_end = min([feature["start"] for feature in features])
        else:
            stop_codon_end = max([feature["end"] for feature in features])# get last start codon
        to_remove = []
        if gene in utr_dict:
            for utr in utr_dict[gene]:
                start = utr["start"]
                end = utr["end"]
                if is_reverse:
                    if end >= stop_codon_end:
                        to_remove.append(utr)
                else:
                    if start <= stop_codon_end:
                        to_remove.append(utr)
            [utr_dict[gene].remove(utr) for utr in to_remove ]


def get_longest_utr(utr_dict):
    """
    Start of the composite utr is the most 5p start, end is the most 3p end.
    """
    gtf = []
    for gene_id, values  in utr_dict.iteritems():
        if not values:
            continue
        if len(values) == 1:
            values[0]["start"] = str(values[0]["start"])
            values[0]["end"] = str(values[0]["end"])
            gtf.append( "\t".join( values[0].values() ) )
        else:
            start = min( [fields["start"] for fields in values] )
            end = max( [fields["end"] for fields in values] )
            values[0]["start"] = str(start)
            values[0]["end"] = str(end)
            gtf.append( "\t".join( values[0].values() ) )
    return gtf

def main():
    utr_dict = OrderedDict()
    start_codon_dict = {}
    header = []
    gtf_fields = get_gtf_fields()
    with open(sys.argv[1]) as input:
        for line in input:
            if line.startswith("#"):
                header.append( line.strip() )
            else:
                get_feature_dict(line, gtf_fields, utr_dict, feature="UTR")
                get_feature_dict(line, gtf_fields, start_codon_dict, feature="CDS")
    remove_five_prime_utrs(utr_dict, start_codon_dict)
    gtf = header + get_longest_utr(utr_dict)
    if len(sys.argv) == 3:
        with open(sys.argv[2], "w") as output:
            [output.write(line + "\n") for line in gtf]
    else:
        for line in gtf:
            print(line)

if __name__ == "__main__":
    main()







