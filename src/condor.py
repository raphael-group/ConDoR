#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 5 2021

@author: Palash Sashittal
"""

import pandas as pd
import sys
import argparse
import itertools
import math
import numpy as np

from solveConstrainedDollo import solveConstrainedDollo

def tree_to_newick(T, root=None):
    if root is None:
        roots = list(filter(lambda p: p[1] == 0, T.in_degree()))
        assert 1 == len(roots)
        root = roots[0][0]
    subgs = []
    while len(T[root]) == 1:
        root = list(T[root])[0]
    for child in T[root]:
        while len(T[child]) == 1:
            child = list(T[child])[0]
        if len(T[child]) > 0:
            child_newick = tree_to_newick(T, root=child)
            if child_newick != '()':
                subgs.append(child_newick)
        else:
            #if child.startswith('s'):
            subgs.append(child)
    # return "(" + ','.join(map(str, subgs)) + ")"
    if len(subgs) == 1:
        return str(subgs[0])
    else:
        return "(" + ','.join(map(str, subgs)) + ")"

def main(args):

    df_character_matrix = pd.read_csv(f'{args.i}', index_col = 0)
    if args.r is not None:
        df_total_readcounts = pd.read_csv(f'{args.r}', index_col = 0)
    if args.v is not None:
        df_variant_readcounts = pd.read_csv(f'{args.v}', index_col = 0)
    
    snp_list = []
    if args.s is not None:
        with open(args.s, 'r') as inp:
            for line in inp:
                snp_list.append(line.rstrip('\n'))
    
    k = args.k
    fp = args.a
    fn = args.b
    ado = args.ado
    
    if args.r is not None:
        solver = solveConstrainedDollo(df_character_matrix, df_total_readcounts=df_total_readcounts,
                                       df_variant_readcounts=df_variant_readcounts, k=k, fp=fp, fn=fn,
                                       ado_precision = ado, snp_list=snp_list)
    else:
        solver = solveConstrainedDollo(df_character_matrix, k=k, fp=fp, fn=fn)
        
    solver.solveSetInclusion()

    prefix = args.o
    solver.writeSolution(f'{prefix}_B.csv')
    solver.writeDOT(f'{prefix}_tree.dot')
    solver.writeDOT(f'{prefix}_tree_without_cells.dot', withcells=False)
    
    if solver.solT_cell is not None:
        with open(f'{prefix}_tree.newick', 'w') as out:
                out.write(tree_to_newick(solver.solT_cell) + ';')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help='csv file with mutation matrix and cluster id')
    parser.add_argument('-r', type=str, help='csv file with total read count matrix')
    parser.add_argument('-v', type=str, help='csv file with variant read count matrix')
    # parser.add_argument('--ado', type=float, help='allelic dropout rate [0.15]', default=0.15)
    parser.add_argument('-s', type=str, help='file containing list of SNPs')
    parser.add_argument('-a', type=float, help='false positive error rate [0.001]', default = 0.001)
    parser.add_argument('-b', type=float, help='false negative error rate [0.001]', default = 0.001)
    parser.add_argument('--ado', type=float, help='precision parameter for ADO', default=15)
    parser.add_argument('-k', type=int, help='maximum number of losses for an SNV', default = 0)
    parser.add_argument('-o', type=str, help='output prefix', required=True)
    parser.add_argument('-t', type=int, help='time limit in seconds [1800]', default = 1800)
    
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    if args.i is None and (args.r is None and args.v is None):
        raise Exception("please provide either the binarized mutation matrix, or the total and variant readcount matrices!")
    
    main(args)
