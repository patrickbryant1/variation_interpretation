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
    hgvc_c_to_pos = pd.read_csv('hgvc_c_to_pos.csv')
except:
    variants.hgvs_c.to_csv('variants_hgvs_c.csv',index=None)
    print('Need to convert hgvs_c to chromosomal positions')
    sys.exit()
# pdb.set_trace()
# #dbNSFP needs chromosome, pos, DNAref, DNAalt
# with open(outfile,'w') as file:
#     for index,row in variants.iterrows():
