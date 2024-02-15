#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 28 2022

@author: Palash Sashittal
"""

import pandas as pd
import sys
import argparse
import math
import numpy as np

def main(args):

    n = args.n
    m = args.m
    prefix = args.o
    cell_names_string = ' '.join([f'sc{idx}' for idx in range(1,n+1)])

    with open(f'{prefix}_cell_names.txt', 'w') as out:
        out.write(cell_names_string)
    
    with open(f'{prefix}_gene_names.txt', 'w') as out:
        for idx in range(m):
            out.write(f'{idx}\n')
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', type=int, help='number of samples', required=True)
    parser.add_argument('-m', type=int, help='number of SNV mutations', required=True)
    parser.add_argument('-o', type=str, help='output prefix', required=True)
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)
