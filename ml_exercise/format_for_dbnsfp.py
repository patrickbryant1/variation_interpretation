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

parser.add_argument('--tsv', nargs=1, type= str,
                  default=sys.stdin, help = 'Pathogenic and neutral variants.')
parser.add_argument('--outfile', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to outfile')


#####MAIN#####
args = parser.parse_args()
variants = pd.read_csv(args.tsv[0], sep='\t')
outfile = args.outfile[0]

#See if the positions are available
try:
    hgvs_c_to_pos = pd.read_csv('hgvs_c_to_pos.tsv',sep='\t')
    hgvs_c_to_pos = hgvs_c_to_pos[['hgvs_c','Chromosomal Variant']]
except:
    variants.hgvs_c.to_csv('variants_hgvs_c.csv',index=None)
    print('Need to convert hgvs_c to chromosomal positions')
    sys.exit()
#Merge the variants with the chr pposition
variants = pd.merge(variants,hgvs_c_to_pos,on='hgvs_c',how='left')
#Remove missing
keep_indices = variants['Chromosomal Variant'].dropna().index
print(str(len(variants)-len(keep_indices))+' out of '+str(len(variants))+' variants w/o chr mapping')
variants = variants.loc[keep_indices]
variants = variants.reset_index()

#dbNSFP needs chromosome, pos, DNAref, DNAalt
chromosomes = []
positions = []
DNAref = []
DNAalt = []
with open(outfile,'w') as file:
    for index,row in variants.iterrows():
        info = row['Chromosomal Variant'].split('.')
        chr = int(info[0].split('_')[1])
        pos = info[2][:-3]
        ref = info[2][-3]
        alt = info[2][-1]
        file.write(str(chr)+' '+pos+' '+ref+' '+alt+'\n')
        #Save
        chromosomes.append(chr)
        positions.append(pos)
        DNAref.append(ref)
        DNAalt.append(alt)

#Add info to df
variants['chr']=chromosomes
variants['pos']=positions
variants['DNAref']=DNAref
variants['DNAalt']=DNAalt
#Save df
variants.to_csv('mapped_variants.csv',index=None)
pdb.set_trace()
