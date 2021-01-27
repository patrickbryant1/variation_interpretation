#!/usr/bin/env bash

TSV=VTD_SCR_TESTD_D1_F1_ClinVar_one_Star.tsv
OUT=query.txt
./format_for_dbnsfp.py --tsv $TSV --outfile $OUT
