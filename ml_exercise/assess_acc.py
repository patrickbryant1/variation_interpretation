#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve, roc_auc_score
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
    Sensitivity = TP/P
    Specificity = TN/N
    PPV = TP/(TP+FP)
    NPV = TN/(TN+FN)
    Accuracy = pred correct/all
    MCC = (TP*TN-FP*FN)/np.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    nMCC = (1+MCC)/2
    OPM = ((PPV+NPV)*(Sensitivity+Specificity)*(Accuracy+nMCC))/8
    AUC - from scikit
    '''
    #Get the positive and neg
    P = np.argwhere(true_sig==1)[:,0]
    N = np.argwhere(true_sig==0)[:,0]
    #Get the pred pos and neg
    pred_P = np.argwhere(pred_sig==1)[:,0]
    pred_N = np.argwhere(pred_sig==0)[:,0]
    #TP and TN
    TP = np.intersect1d(P,pred_P).shape[0]
    FP = len(pred_P)-TP
    TN = np.intersect1d(N,pred_N).shape[0]
    FN= len(pred_N)-TN
    #Metrics
    Sensitivity = TP/P.shape[0]
    Specificity = TN/N.shape[0]
    PPV = TP/(TP+FP)
    NPV = TN/(TN+FN)
    Accuracy = np.argwhere(true_sig==pred_sig).shape[0]/len(pred_sig)
    MCC = (TP*TN-FP*FN)/np.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    nMCC = (1+MCC)/2
    OPM = ((PPV+NPV)*(Sensitivity+Specificity)*(Accuracy+nMCC))/8
    AUC  = roc_auc_score(true_sig,pred_sig)

    return Sensitivity, Specificity, PPV, NPV, Accuracy, MCC, OPM, AUC


#####MAIN#####
args = parser.parse_args()
pred = pd.read_csv(args.pred[0], sep='\t')
true = pd.read_csv(args.true[0])
outdir = args.outdir[0]

#Pred contains duplicates
keep_inds = pred[['chr','pos','ref','alt']].drop_duplicates().index
pred = pred.loc[keep_inds]
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

#Replace pathogenic with 1 and benighn with 0
selection = selection.replace('Pathogenic',1)
selection = selection.replace('Benign',0)

Sensitivities = []
Specificities = []
PPVs =[]
NPVs = []
Accuracies = []
MCCs = []
OPMs = []
AUCs = []
N_predicted = []
#Evaluate preds
for predictor in selection.columns[1:]:

    #Get only the indices where something has been predicted
    pred_inds = selection[(selection[predictor]==0) | (selection[predictor]==1)].index
    print(predictor,len(pred_inds),'out of',len(selection))
    Sensitivity, Specificity, PPV, NPV, Accuracy, MCC, OPM, AUC = evaluate(selection.loc[pred_inds,'clinical_significance'].values,selection.loc[pred_inds,predictor].values)
    #Save
    Sensitivities.append(Sensitivity)
    Specificities.append(Specificity)
    PPVs.append(PPV)
    NPVs.append(NPV)
    Accuracies.append(Accuracy)
    MCCs.append(MCC)
    OPMs.append(OPM)
    AUCs.append(AUC)
    N_predicted.append(len(pred_inds))

#Create df
score_df = pd.DataFrame()
score_df['Model']=selection.columns[1:]
score_df['Sensitivity']=Sensitivities
score_df['Specificity']=Specificities
score_df['PPV']=PPVs
score_df['NPV']=NPVs
score_df['Accuracy']=Accuracies
score_df['MCC']=MCCs
score_df['OPM']=OPMs
score_df['AUC']=AUCs
score_df['Number of predictions']=N_predicted

#Save
score_df = score_df.round(2)
score_df.to_csv(outdir+'score_df.csv')
