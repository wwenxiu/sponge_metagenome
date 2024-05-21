#!/usr/bin/env python
# coding: utf-8
import pandas as pd
from collections import defaultdict
import sys, io

fn1=sys.argv[1] 
fn2=sys.argv[2] 

if __name__ == "__main__":
    aligned_contig = pd.read_csv(fn1, sep='\t', header=None, names=['read_id', 'contig_id'])
    contig_counts = aligned_contig['read_id'].value_counts()
    read_ids = contig_counts[contig_counts > 1].index
    grouped_contigs = defaultdict(list)

    for read_id, contig_id in aligned_contig.values:
        if read_id in read_ids:
            grouped_contigs[read_id].append(contig_id)

    with io.open(fn2, 'w', encoding='utf-8') as f:
        for read_id, contig_ids in grouped_contigs.items():
            sorted_contig_ids = sorted(contig_ids)
            if sorted_contig_ids[0] != sorted_contig_ids[1]:
                f.write(f"{read_id}\t{sorted_contig_ids[0]}\t{sorted_contig_ids[1]}\t1\n")