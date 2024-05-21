##!/usr/bin/env python
# coding: utf-8

import pandas as pd
import sys

fn1=sys.argv[1] 
fn2=sys.argv[2] 
fn3=sys.argv[3] 

if __name__ == "__main__":
    alllink_df = pd.read_csv(fn1,sep='\t',header=None)
    alllink_df.rename(columns={0:'read_name',1:'contig_1',2:'contig_2', 3:'number'},inplace=True)
    
    filter_reads = pd.read_csv(fn2, sep='\t', header=None)
    filter_reads.rename(columns={0:'read_name'},inplace=True)
    filter_reads['filter']='filter'

    merge_df = pd.merge(alllink_df,filter_reads,how='left')
    merge_df = merge_df[merge_df['filter'] != 'filter']
    merge_df = merge_df[['read_name','contig_1','contig_2','number']]
    merge_df.to_csv(fn3, sep='\t',index=0,header=None)