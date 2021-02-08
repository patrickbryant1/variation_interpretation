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


#####MAIN#####
args = parser.parse_args()
pred = pd.read_csv(args.pred[0], sep='\t')
true = pd.read_csv(args.true[0])
outdir = args.outdir[0]

pdb.set_trace()
#Add info to df
variants['chr']=chromosomes
variants['pos']=positions
variants['DNAref']=DNAref
variants['DNAalt']=DNAalt
#Save df
variants.to_csv('mapped_variants.csv',index=None)
pdb.set_trace()
