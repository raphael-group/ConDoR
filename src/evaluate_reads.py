#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 5 2022

@author: Palash Sashittal
"""

import pandas as pd
import sys
import argparse
import itertools
import math
import numpy as np

import networkx as nx
import ete3
from solveConstrainedDollo import solveConstrainedDollo

sys.setrecursionlimit(1500)

def get_scarlet_tree_dataframe(scarlet_fname, scarlet_edgelist_fname):
    
    df_corrected = pd.read_csv(scarlet_fname, index_col = 0)
    
    T_scarlet_nx = nx.DiGraph()
    with open(scarlet_edgelist_fname, 'r') as inp:
        for line in inp:
            nodes = line.rstrip('\n').split(',')
            parent_node = nodes[0].split(' ')[0]
            child_node = nodes[1].split(' ')[0]

            parent_type = parent_node.split(':')[0]
            parent_name = parent_node.split(':')[1]
            child_type = child_node.split(':')[0]
            child_name = child_node.split(':')[1]

            if parent_type == 'ROOT':
                parent_node_name = f'r{parent_name}'
            elif parent_type == 'MUT':
                parent_node_name = f'{parent_name}_1'
            else:
                parent_node_name = parent_name

            if child_type == 'ROOT':
                child_node_name = f'r{child_name}'
            elif child_type == 'MUT':
                child_node_name = f'{child_name}_1'
            else:
                child_node_name = child_name

            T_scarlet_nx.add_edge(parent_node_name, child_node_name)
    
    Tsol = ete3.Tree(tree_to_newick(T_scarlet_nx) + ';')
    
    return Tsol, df_corrected, T_scarlet_nx

def get_condor_tree_dataframe(condor_fname, condor_newick):
    df_multi = pd.read_csv(condor_fname, index_col = 0)

    df_corrected = df_multi.copy()
    df_corrected[df_corrected > 1] = 0

    Tsol = ete3.Tree(condor_newick, format=1)
    
    df_binary = solveConstrainedDollo.expand_multi_state_to_binary(df_multi)
    _, nxtree_condor = solveConstrainedDollo.generate_perfect_phylogeny(df_binary)    
    
    return Tsol, df_corrected, nxtree_condor

def read_sphyr(fname):
    
    with open(fname, 'r') as inp:
        idx = 0
        data = []
        for line in inp:
            if idx == 0:
                n = int(line.split(' ')[0])
            elif idx == 1:
                m = int(line.split(' ')[0])
            else:
                data.append(list(map(int, line.split(' '))))
            idx += 1
    
    return pd.DataFrame(data, columns = [f'c{idx}' for idx in range(m)], index = [f's{idx}' for idx in range(n)])

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
            if child.startswith('s'):
                subgs.append(child)
    # return "(" + ','.join(map(str, subgs)) + ")"
    if len(subgs) == 1:
        return str(subgs[0])
    else:
        return "(" + ','.join(map(str, subgs)) + ")"

def get_sphyr_tree_dataframe(sphyr_fname):

    df_multi = read_sphyr(sphyr_fname)
    df_binary = solveConstrainedDollo.expand_multi_state_to_binary(df_multi)
    _, nxtree_sphyr = solveConstrainedDollo.generate_perfect_phylogeny(df_binary)

    df_corrected = df_multi.copy()
    df_corrected[df_corrected > 1] = 0

    T_sphyr = ete3.Tree(tree_to_newick(nxtree_sphyr)+';', format=1)
    
    return T_sphyr, df_corrected, nxtree_sphyr

def get_sifit_tree_dataframe(sifit_fname, sifit_newick):
    df_corrected = pd.read_csv(sifit_fname, sep='\t', header=None, index_col = 0)
    n = len(df_corrected)
    m = len(df_corrected.columns)
    
    df_corrected.rename(columns={idx: f'c{idx-1}' for idx in range(1,m+1)}, index={f'sc{idx}': f's{idx-1}' for idx in range(1, n+1)}, inplace=True)
    df_corrected = df_corrected.loc[[f's{idx}' for idx in range(n)]]
    df_corrected = df_corrected.sort_index()
    
    sifit_newick_string = ''
    with open(sifit_newick, 'r') as inp:
        for line in inp:
            sifit_newick_string += line.rstrip('\n')    

    mod_sifit_newick_list = []
    for idx, string in enumerate(sifit_newick_string.split('sc')):
        if idx == 0:
            mod_sifit_newick_list.append(string)
        else:
            splited_string = string.split(':')
            splited_string[0] = str(int(splited_string[0]) - 1)
            mod_string = ':'.join(splited_string)
            mod_sifit_newick_list.append('s' + mod_string)    

    T_sifit = ete3.Tree(''.join(mod_sifit_newick_list))
    
    return T_sifit, df_corrected

def get_scite_tree_dataframe(n, m, gv_fname):

    # scite_nx_tree = nx.DiGraph()
    # with open(gv_fname, 'r') as inp:
    #     for line in inp:
    #         if not line.startswith('digraph') and not line.startswith('node') and not line.startswith('}'):
    #             data = line.rstrip(';\n').split(' -> ')
    #             scite_nx_tree.add_edge(data[0], data[1])
    # T_scite = ete3.Tree(tree_to_newick(scite_nx_tree) + ';')
    
    gv_file = open(gv_fname, 'r')
    gv_file.readline()
    gv_file.readline()    

    pi = [-2 for i in range(m+1)]
    for _ in range(m):
        line = gv_file.readline()
        data = line.rstrip(";\n").split()
        source = int(data[0]) - 1
        target = int(data[2]) - 1
        assert 0 <= target < m
        parent = -2
        if source == m:
            parent = -1
        else:
            parent = source
        pi[target] = parent    

    samples = [-2 for i in range(n)]
    gv_file.readline()

    for _ in range(n):
        line = gv_file.readline()
        data = line.rstrip(";\n").split()
        source = int(data[0]) - 1
        target = int(data[2][1:])
        parent = -2
        if source == n:
            parent = -1
        else:
            parent = source
        samples[target] = source        
        
    scite_corrected = []
    for p in range(len(samples)):
        states = [ 0 for c in range(len(pi)) ]
        parent = samples[p]
        while parent != -1:
            #print p, parent
            states[parent] = 1
            parent = pi[parent]
        scite_corrected.append(states[:len(pi) - 1])
        
    df_corrected = pd.DataFrame(scite_corrected, columns = [f'c{idx}' for idx in range(m)], index = [f's{idx}' for idx in range(n)])
    
    df_corrected_modified = df_corrected.copy()
    df_corrected_modified.columns = [f'{x}_1' for x in df_corrected.columns]
    _, nxtree_scite = solveConstrainedDollo.generate_perfect_phylogeny(df_corrected_modified)
    T_scite = ete3.Tree(tree_to_newick(nxtree_scite) + ';')
    
    return T_scite, df_corrected, nxtree_scite

def get_descendant_mutations(T, node):
    if node not in T.nodes:
        return []
    if node.startswith('c'):
        descendants = [node]
    else:
        return []
    for child in T[node]:
        descendants += get_descendant_mutations(T, child)
    return descendants

def get_clustered_mutations(T, node):
    if node not in T.nodes:
        return []
    if node.startswith('c'):
        clustered = [node]
        if len(list(T[node])) == 1:
            child = list(T[node])[0]
            if child.startswith('c'):
                clustered += get_clustered_mutations(T, child)
    else:
        return []
    return clustered

def main(args):

    n = args.n
    m = args.m
    p = args.p
    k = args.k
    d = args.d
    if d == 0:
        d = int(d)
    s = args.s
    
    output_fname = args.o
    
    method = args.method
    simulation_dir = args.dir
    
    if method == 'condor':
        
        condor_fname = f'{simulation_dir}/condor/n{n}_m{m}_p{p}_k{k}_d{d}_s{s}_B.csv'
        condor_newick = f'{simulation_dir}/condor/n{n}_m{m}_p{p}_k{k}_d{d}_s{s}_tree.newick'
        Tsol, df_corrected, nxtree_sol = get_condor_tree_dataframe(condor_fname, condor_newick)
        
    if method == 'condorreads':

        condor_fname = f'{simulation_dir}/condor_reads/n{n}_m{m}_p{p}_k{k}_d{d}_s{s}_B.csv'
        condor_newick = f'{simulation_dir}/condor_reads/n{n}_m{m}_p{p}_k{k}_d{d}_s{s}_tree.newick'
        Tsol, df_corrected, nxtree_sol = get_condor_tree_dataframe(condor_fname, condor_newick)
        
    elif method == 'sphyr':
        
        sphyr_fname = f'{simulation_dir}/sphyr/n{n}_m{m}_p{p}_k{k}_d{d}_s{s}.out'
        Tsol, df_corrected, nxtree_sol = get_sphyr_tree_dataframe(sphyr_fname)
        
    elif method == 'sifit':
        
        sifit_fname = f'{simulation_dir}/sifit/n{n}_m{m}_p{p}_k{k}_d{d}_s{s}_character_matrix.tsv'        
        sifit_newick = f'{simulation_dir}/sifit/n{n}_m{m}_p{p}_k{k}_d{d}_s{s}_mlTree.newick'

        Tsol, df_corrected = get_sifit_tree_dataframe(sifit_fname, sifit_newick)
        
    elif method == 'scite':
        
        # scite_newick = f'{simulation_dir}/scite/n{n}_m{m}_p{p}_k{k}_d{d}_s{s}_ml0.newick'
        scite_gv_fname = f'{simulation_dir}/scite/n{n}_m{m}_p{p}_k{k}_d{d}_s{s}_ml0.gv'
        Tsol, df_corrected, nxtree_sol = get_scite_tree_dataframe(n, m, scite_gv_fname)
    
    elif method == 'scarlet':
        
        scarlet_fname = f'{simulation_dir}/scarlet/n{n}_m{m}_p{p}_k{k}_d{d}_s{s}.B'
        scarlet_edgelist_fname = f'{simulation_dir}/scarlet/n{n}_m{m}_p{p}_k{k}_d{d}_s{s}.edgelist'
        Tsol, df_corrected, nxtree_sol = get_scarlet_tree_dataframe(scarlet_fname, scarlet_edgelist_fname)
    
    # evaluate
    gt_fname = f'{simulation_dir}/ground_truth/n{n}_m{m}_p{p}_k{k}_d{d}_s{s}_character_matrix_without_noise.csv'
    gt_newick = f'{simulation_dir}/ground_truth/n{n}_m{m}_p{p}_k{k}_d{d}_s{s}_tree.newick'
    gt_multi = f'{simulation_dir}/ground_truth/n{n}_m{m}_p{p}_k{k}_d{d}_s{s}_multi_state_character_matrix.csv'
    df_character_matrix = pd.read_csv(gt_fname, index_col = 0)
    df_character_matrix = df_character_matrix[df_character_matrix.columns[:-1]]
    Tground = ete3.Tree(gt_newick, format=1)
    df_Bcell = pd.read_csv(gt_multi, index_col = 0)
    _, nxtree_gt = solveConstrainedDollo.generate_perfect_phylogeny(solveConstrainedDollo.expand_multi_state_to_binary(df_Bcell[df_Bcell.columns[:-1]]))
    
    # RF metric
    (rf, rf_max, names, edges_t1, edges_t2, discarded_edges_t1, discarded_edges_t2) = Tground.robinson_foulds(Tsol, unrooted_trees=True)
    # mutation error
    nerror = np.sum(abs(df_corrected - df_character_matrix).values)
    
    # recall for ancestry and incomparibility
    if method != 'sifit':
        sol_descendant_dictionary = {}
        for mutation_idx in range(m):
            sol_descendant_dictionary[f'c{mutation_idx}_1'] = get_descendant_mutations(nxtree_sol, f'c{mutation_idx}_1')

        sol_clustered_dictionary = {}
        for mutation_idx in range(m):
            sol_clustered_dictionary[f'c{mutation_idx}_1'] = get_clustered_mutations(nxtree_sol, f'c{mutation_idx}_1')            
            
        gt_descendant_dictionary = {}
        for mutation_idx in range(m):
            gt_descendant_dictionary[f'c{mutation_idx}_1'] = get_descendant_mutations(nxtree_gt, f'c{mutation_idx}_1')    

        gt_clustered_dictionary = {}
        for mutation_idx in range(m):
            gt_clustered_dictionary[f'c{mutation_idx}_1'] = get_clustered_mutations(nxtree_gt, f'c{mutation_idx}_1')
            
#         confusion_mat = np.zeros((3,3))
#         for mutation_idx1, mutation_idx2 in itertools.combinations(range(m), 2):
#             mutation1 = f'c{mutation_idx1}_1'
#             mutation2 = f'c{mutation_idx2}_1'

#             if mutation2 in sol_descendant_dictionary[mutation1] and mutation1 in sol_descendant_dictionary[mutation2]:
#                 print('problem sol')
#                 print(mutation1, mutation2)
#                 raise ValueError
#                 break    

#             x_idx = 2
#             if mutation2 in sol_descendant_dictionary[mutation1]:
#                 x_idx = 0
#             if mutation1 in sol_descendant_dictionary[mutation2]:
#                 x_idx = 1

#             if mutation2 in gt_descendant_dictionary[mutation1] and mutation1 in gt_descendant_dictionary[mutation2]:
#                 print('problem gt')
#                 print(mutation1, mutation2)
#                 raise ValueError
#                 break

#             y_idx = 2
#             if mutation2 in gt_descendant_dictionary[mutation1]:
#                 y_idx = 0
#             if mutation1 in gt_descendant_dictionary[mutation2]:
#                 y_idx = 1    

#             confusion_mat[x_idx, y_idx] += 1

        confusion_mat = np.zeros((4,4))
        for mutation_idx1, mutation_idx2 in itertools.combinations(range(m), 2):
            mutation1 = f'c{mutation_idx1}_1'
            mutation2 = f'c{mutation_idx2}_1'

            if mutation2 in sol_descendant_dictionary[mutation1] and mutation1 in sol_descendant_dictionary[mutation2]:
                print('problem sol')
                print(mutation1, mutation2)
                break    

            x_idx = 3
            if mutation1 in sol_clustered_dictionary[mutation2] or mutation2 in sol_clustered_dictionary[mutation1]:
                x_idx = 2
            else:
                if mutation2 in sol_descendant_dictionary[mutation1]:
                    x_idx = 0
                if mutation1 in sol_descendant_dictionary[mutation2]:
                    x_idx = 1

            if mutation2 in gt_descendant_dictionary[mutation1] and mutation1 in gt_descendant_dictionary[mutation2]:
                print('problem gt')
                print(mutation1, mutation2)
                break

            y_idx = 3
            if mutation1 in gt_clustered_dictionary[mutation2] or mutation2 in gt_clustered_dictionary[mutation1]:
                y_idx = 2
            else:
                if mutation2 in gt_descendant_dictionary[mutation1]:
                    y_idx = 0
                if mutation1 in gt_descendant_dictionary[mutation2]:
                    y_idx = 1

            confusion_mat[x_idx, y_idx] += 1

        ancestry_recall = (confusion_mat[0,0] + confusion_mat[1,1]) /  np.sum(confusion_mat[:,:2])
        # incomparability_recall = confusion_mat[2,2] / np.sum(confusion_mat[:,2])
        incomparability_recall = np.nan
        accuracy = np.trace(confusion_mat)/np.sum(confusion_mat)                
    else:
        ancestry_recall = np.nan
        incomparability_recall = np.nan
        accuracy = np.nan
    
    df_result = pd.DataFrame([[n, m, p, k, d, s, method, rf, rf_max, nerror, ancestry_recall, incomparability_recall, accuracy]],
                             columns = ['ncells', 'ncharacters', 'nclusters', 'k', 'dropout', 'seed', 'method',
                                        'RF', 'RF_max', 'nerror', 'ancestry_recall', 'incomparibility_recall', 'relation_accuracy'])
    
    df_result.to_csv(f'{output_fname}')
    
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', type=int, help='number of samples [5]', default = 5)
    parser.add_argument('-m', type=int, help='number of SNV mutations [5]', default = 5)
    parser.add_argument('-p', type=int, help='number of clusters [1]', default = 1)
    parser.add_argument('-k', type=int, help='number of SNV losses per character [0]', default = 0)
    parser.add_argument('-o', type=str, help='output filename', default='sample.csv')
    parser.add_argument('-s', type=int, help='seed [0]', default = 0)
    parser.add_argument('-d', type=float, help='missing data rate [0.0]', default=0)
    parser.add_argument('--dir', type=str, help='simulation directory', required=True)
    parser.add_argument('--method', type=str, help='method for comparison with ground truth', required=True)
    parser.add_argument('-v', action='store_true', default=False)
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)