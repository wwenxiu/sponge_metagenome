#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import sys

fn1=sys.argv[1] 
fn2=sys.argv[2] 

def get_df_V(df):
    filter_mask = ((df['contig_1'].str.contains('bin', na=False)) ^ (df['contig_2'].str.contains('bin', na=False)))
    df_filtered = df[filter_mask].copy()
    df_filtered.reset_index(drop=True, inplace=True)
    return df_filtered

def change(df):
    mask = df['contig_2'].str.contains('bin', na=False)
    df.loc[mask, ['contig_1', 'contig_2']] = df.loc[mask, ['contig_2', 'contig_1']].values
    return df

if __name__ == "__main__":
    link_table = pd.read_csv(fn1, sep='\t', header=None, names=['read', 'contig_1', 'contig_2', 'number'])
    group_link_table = link_table.groupby(['contig_1', 'contig_2'])['number'].sum().reset_index()
    group_link_table = group_link_table[(group_link_table['number'] >= 0)] 
    group_link_table_V = get_df_V(group_link_table)
    group_link_table_V_reorder = change(group_link_table_V)
    group_link_table_V_reorder['bin_id'] = group_link_table_V_reorder['contig_1'].str.rsplit("_", 1, expand=True)[0]
    bin_V_table = group_link_table_V_reorder.groupby(['bin_id', 'contig_2'])['number'].sum().reset_index()

    writer = pd.ExcelWriter(fn2)
    group_link_table.to_excel(writer, sheet_name='filter', index=False)
    group_link_table_V_reorder.to_excel(writer, sheet_name='filter_V', index=False)
    bin_V_table.to_excel(writer, sheet_name='bin_V', index=False)
    writer.save()