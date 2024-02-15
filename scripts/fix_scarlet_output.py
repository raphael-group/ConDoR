#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 2022

@author: Palash Sashittal
"""

import pandas as pd
import sys
import argparse
import itertools
import math
import numpy as np

import networkx as nx

def main(args):

    df_B_ancestor = pd.read_csv(args.i, index_col = 0)
    df_character_matrix = pd.read_csv(args.d, index_col = 0)
    
    df_B_ancestor_fixed = df_B_ancestor.copy()
    for cell in df_character_matrix.index:
        if cell not in df_B_ancestor.index:
            cluster_id = df_character_matrix.loc[cell]['cluster_id']
            ancestor_data = df_B_ancestor.loc[f'ANC:{cluster_id}']
            ancestor_data['CN'] = cluster_id
            df_B_ancestor_fixed.loc[cell] = ancestor_data

    df_B_ancestor_fixed.to_csv(f'{args.o}.B_ancestor_fixed')
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help='possibly incorrect scarlet file', required=True)
    parser.add_argument('-d', type=str, help='input character matrix', required=True)
    parser.add_argument('-o', type=str, help='output prefix', required=True)
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)