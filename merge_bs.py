import numpy as np
import json
import yaml
from argparse import ArgumentParser
import os
import ROOT
import re
from unfold_utils import *
from arg_parsing import *
from configs import ObsConfig
import pandas as pd
from GOF.binned import *
import uproot
import awkward as ak
import subprocess
import ast
from multiprocessing.pool import Pool
from itertools import combinations_with_replacement

def append_df(df_origin,df_append):
  if df_origin is None:
    df_origin = df_append.copy()
  else:
    df_origin = df_origin.append(df_append,ignore_index = True)
  df_origin.drop_duplicates(ignore_index = True, inplace = True)
  return df_origin

def get_bs_entries(df,label_row):
  for ikey, key in enumerate(list(label_row.keys())):
    if ikey == 0:
      sel = (df[key]==label_row[key])
    else:
      sel = sel & (df[key]==label_row[key])
  return df[sel]

def df_column_to_ak(column):
  array = column.replace('inf','2e308')
  array = array.apply(ast.literal_eval)
  array = array.values
  array = list(array)
  array = ak.Array(array)
  return array

def write_mean_std(entry_merge,entries_origin):
  bin_values = entries_origin['bin_values']
  bin_values = df_column_to_ak(bin_values)
  entry_merge['bin_values'] = ak.to_list(ak.mean(bin_values,axis = 0))
  entry_merge['bin_errors'] = ak.to_list(np.sqrt((ak.mean(bin_values*bin_values,axis = 0) - np.square(ak.mean(bin_values,axis = 0)))*len(bin_values)/(len(bin_values)-1)))

def write_cov_corr(entry_merge_cov,entry_merge_corr,entries_origin):
  bin_values_x = entries_origin['bin_values_x']
  bin_values_y = entries_origin['bin_values_y']
  bin_values_x = df_column_to_ak(bin_values_x)
  bin_values_y = df_column_to_ak(bin_values_y)
  bin_values_x =   bin_values_x[:,:,1:-1] #remove the underflow and overflow part in dim2
  bin_values_y =   bin_values_y[:,:,1:-1]
  bin_values_x = ak.flatten(bin_values_x,axis=-1)
  bin_values_y = ak.flatten(bin_values_y,axis=-1)
  cov = np.cov(bin_values_x,bin_values_y,rowvar=False)[:np.shape(bin_values_x)[-1],np.shape(bin_values_x)[-1]:]
  entry_merge_cov['bin_values'] = cov.to_list()
  entry_merge_cov['bin_errors'] = None
  corr = np.nan_to_num(np.corrcoef(bin_values_x,bin_values_y,rowvar=False))[:np.shape(bin_values_x)[-1],np.shape(bin_values_x)[-1]:]
  entry_merge_corr['bin_values'] = corr.to_list()
  entry_merge_corr['bin_errors'] = None

if __name__=="__main__":

    parser = ArgumentParser()
    parser.add_argument('-i','--input', default="config/merge_unfolddata.txt", help="The text file including the names of the csv files for merging")
    parser.add_argument('-o','--output', default="merge.csv", help="The name of the csv file for output")
    args = parser.parse_args()

    with open(args.input) as file:
      lines = [line.rstrip() for line in file]

    labels = ['obs1','obs2','histtype','iter','dataset','method','from_step1','do_eff_acc','dim1_isgen','dim2_isgen','bin_edges_dim1','bin_edges_dim2']

    df_merge = None
    df_sysbs_merge = None
    for name in lines:
      df = pd.read_csv(name)
      df_non_bs = df[~df['datatype'].str.contains('bs')]
      df_sysbs = df[df['datatype'].str.contains('sysbs')]
      df_merge = append_df(df_merge,df_non_bs)
      df_sysbs_merge = append_df(df_sysbs_merge,df_sysbs)

    for datatype,df_bs_merge in zip(['unfold_sys'],[df_sysbs_merge]):
      df_labels = df_bs_merge[labels]
      df_labels = df_labels.drop_duplicates(ignore_index = True)
      df_labels_cov = df_labels[(df_labels['histtype'] == 'inclusive')&(df_labels['dim1_isgen'])&(df_labels['dim2_isgen'])]
      for _, row in df_labels.iterrows():
        bs_rows = get_bs_entries(df_bs_merge,row)
        bs_row_merge = bs_rows.iloc[0].copy()
        write_mean_std(bs_row_merge,bs_rows)
        bs_row_merge['datatype'] = datatype
        df_merge = df_merge.append(bs_row_merge.copy(),ignore_index = True)
      for rows in combinations_with_replacement(df_labels_cov.iterrows(),2):
        _,row1 = rows[0]
        _,row2 = rows[1]
        if row1['iter'] != row2['iter']: continue
        bs_rows_1 = get_bs_entries(df_bs_merge,row1)
        bs_rows_2 = get_bs_entries(df_bs_merge,row2)
        cov_row  = bs_rows_1.iloc[0].copy()
        for column_obs in ['obs1','obs2','bin_edges_dim1','bin_edges_dim2']:
          cov_row[column_obs] = (row1[column_obs],row2[column_obs])
        cov_row['datatype'] = datatype
        corr_row = cov_row.copy()
        cov_row['histtype'] = 'covariance'
        corr_row['histtype'] = 'correlation'
        bs_rows_combine = bs_rows_1.merge(bs_rows_2,on='datatype')
        write_cov_corr(cov_row,corr_row,bs_rows_combine)
        df_merge = df_merge.append(cov_row.copy(),ignore_index = True)
        df_merge = df_merge.append(corr_row.copy(),ignore_index = True)
    df_merge.to_csv(args.output,mode='w',index=False)

