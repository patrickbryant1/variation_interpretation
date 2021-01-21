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
parser = argparse.ArgumentParser(description = '''Convert csv files to fasta format.''')

parser.add_argument('--pathogenic', nargs=1, type= str,
                  default=sys.stdin, help = 'Pathogenic variants.')
parser.add_argument('--neutral', nargs=1, type= str,
                  default=sys.stdin, help = 'Neutral variants.')
parser.add_argument('--id_map', nargs=1, type= str,
                  default=sys.stdin, help = 'Map from uniprot id to protein name.')
parser.add_argument('--outfile', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to outfile')


#####MAIN#####
args = parser.parse_args()
pathogenic = pd.read_csv(args.pathogenic[0], sep='\t')
neutral = pd.read_csv(args.neutral[0], sep='\t')
id_map = pd.read_csv(args.id_map[0], sep='\t')
id_map = id_map.rename(columns={"From": "ID_UNIPROT", "To": "PROTEIN_NAME"})
outfile = args.outfile[0]

#Assign pathogenicity
pathogenic['pathogenicity']=1
neutral['pathogenicity']=0
#Mrge
all = pathogenic.append(neutral)
#Merge to get protein name
all = pd.merge(all,id_map,on='ID_UNIPROT',how='left')
#Format for ensembl variant recoder
all['variant_id']=all['PROTEIN_NAME']+':p.'+all['MUTATION']

server = "http://rest.ensembl.org"
ext = "/variant_recoder/homo_sapiens"
headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
#Write request str
ids = all['variant_id'].values
request_str = '{ "ids" : ['
for id in ids:
    request_str +='"'+id+'",'
request_str = request_str[:-1]+'] }'
#data='{ "ids" : ["A0PJY2:p.H278Y" ] }'
result = requests.post(server+ext, headers=headers, data=request_str)
decoded = result.json()
#Save
file = open("test.json",'w')
json.dump(decoded,file)

#dbNSFP needs chromosome, pos, DNAref, DNAalt
chr = []
pos = []
DNAref = []
DNAalt = []
for d in decoded:
    d_hgvsg = d['hgvsg'][0]
    d_hgvsg = d_hgvsg.split('_')[1]
    chr.append(int(d_hgvsg.split('.')[0]))
    d_hgvsg = d_hgvsg.split('.')[2]
    pos = int(d_hgvsg[:-3])
    DNAref.append(d_hgvsg[-3])
    DNAalt.append(d_hgvsg[-1])

# 'hgvsg': ['NC_000007.14:g.122303281G>A']
#chr 7 band 14 position 122303281 subst G to A
pdb.set_trace()
