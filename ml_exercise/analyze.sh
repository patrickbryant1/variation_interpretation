#!/usr/bin/env bash

PRED=./q100_result.tsv
TRUE=./mapped_variants.csv
OUTDIR=./
./assess_acc.py --pred $PRED --true $TRUE --outdir $OUTDIR
