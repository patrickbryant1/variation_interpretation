#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import numpy as np
import pandas as pd
import requests
import json
import pdb
#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Write the format needed for dbNSFP.''')

parser.add_argument('--pred', nargs=1, type= str,
                  default=sys.stdin, help = 'Predicted pathogenic and neutral variants.')
parser.add_argument('--true', nargs=1, type= str,
                  default=sys.stdin, help = 'True pathogenic and neutral variants.')
parser.add_argument('--outdir', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to outfile')


#####FUNCTIONS#####
def evaluate(true_sig,pred_sig):
    '''Compare the true and pred sig
    '''

#####MAIN#####
args = parser.parse_args()
pred = pd.read_csv(args.pred[0], sep='\t')
true = pd.read_csv(args.true[0])
outdir = args.outdir[0]

#Merge pred and true
merged_df = pd.merge(pred, true,  how='left', left_on=['chr','pos','ref','alt'], right_on = ['chr','pos','DNAref','DNAalt'])

predictors = {'SIFT_pred':['D:Pathogenic','T:Benign'],
            'FATHMM_pred':['D:Pathogenic','T:Benign'],
            'MetaSVM_pred':['D:Pathogenic','T:Benign'],
            'M_CAP_pred':['D:Pathogenic','T:Benign'],
            'DEOGEN2_pred':['D:Pathogenic','T:Benign'],
            'BayesDel_addAF_pred':['D:Pathogenic','T:Benign'],
            'ClinPred_pred':['D:Pathogenic','T:Benign']
            }

sel_cols = ['clinical_significance']
sel_cols.extend([*predictors.keys()])
selection = merged_df[sel_cols]
#Convert the predictor scores to Pathogenic or Benign
for key in predictors:
    conversions = predictors[key]
    for item in conversions:
        old,new = item.split(':')
        inds = selection[selection[key]==old].index
        selection.at[inds,key]=new
pdb.set_trace()
#Add info to df

#Save df
variants.to_csv('mapped_variants.csv',index=None)
pdb.set_trace()
