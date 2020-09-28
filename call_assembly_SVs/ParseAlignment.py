#!/usr/bin/env python
import sys
import pandas as pd
import argparse

parser = argparse.ArgumentParser(prog='Parse Alignments')
parser.add_argument('--inputsam',  help='Input sorted sam file', required=True)
parser.add_argument('--parsedsam',  help='Output parsed sam file', required=True)
args = parser.parse_args()

maternal = pd.read_csv(args.inputsam, sep ="\t", header=None)
remove = [0] * maternal.shape[0]
current_end = maternal.iloc[0][8] + maternal.iloc[0][3]
for i in range(1, maternal.shape[0]):
    if maternal.iloc[i-1][2] == maternal.iloc[i][2]:  # check chr is the same
        if maternal.iloc[i][8] + maternal.iloc[i][3] <= current_end: 
            remove[i] = 1
        else:
            current_end = maternal.iloc[i][8] + maternal.iloc[i][3]
    else:
        current_end = maternal.iloc[i][8] + maternal.iloc[i][3]
maternal['remove'] = remove
maternal_filtered = maternal[maternal['remove'] == 0] 
del maternal_filtered['remove']
maternal_filtered.to_csv(args.parsedsam, header=False,  index=False, sep="\t")