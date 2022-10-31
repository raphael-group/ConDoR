#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 5 2022

@author: Palash Sashittal
"""

import gurobipy as gp
import numpy as np
import pandas as pd
import networkx as nx
import itertools
from scipy.stats import betabinom

# minimum correction tree parsimonious clone reconciliation problem
class solveConstrainedDollo():

    def __init__(self, df_character_matrix, df_total_readcounts = None, df_variant_readcounts = None, snp_list = [],
                 k = None, fp = None, fn = None, ado_precision = 15, threads = 1, timelimit = None, verbose = True):
        
        # input character matrix and clustering
        self.df_character_matrix = df_character_matrix
        self.clustering = self.df_character_matrix['cluster_id'].values
        self.A = df_character_matrix.values[:, :-1]
        self.snp_list = snp_list
        self.mutation_list = list(df_character_matrix.columns[:-1])
        
        # input data parameters
        self.ncells = len(self.df_character_matrix)
        self.nclusters = len(self.df_character_matrix['cluster_id'].unique())
        self.nmutations = len(self.df_character_matrix.columns) - 1
        self.k = k
        self.fp = fp
        self.fn = fn        
        
        # read count matrices
        self.df_total_readcounts = df_total_readcounts        
        if df_total_readcounts is not None:
            self.cell_list = list(df_total_readcounts.index)
            self.mutation_list = list(df_total_readcounts.columns)

            bb_alpha = fp * ado_precision
            bb_beta = (1 - fp) * ado_precision
            
            coeff_mat = np.zeros((self.ncells, self.nmutations))
            for cell_idx, cell in enumerate(self.cell_list):
                for mut_idx, mutation in enumerate(self.mutation_list):
                    total_reads = df_total_readcounts.loc[cell][mutation]
                    variant_reads = df_variant_readcounts.loc[cell][mutation]
                    if total_reads > 0:
                        coeff = -np.log(total_reads) - betabinom.logpmf(variant_reads, total_reads, bb_alpha, bb_beta)
                        #coeff = betabinom.logpmf(variant_reads, total_reads, 1, 1) - betabinom.logpmf(variant_reads, total_reads, bb_alpha, bb_beta)
                        coeff_mat[cell_idx, mut_idx] = coeff
            self.coeff_mat = coeff_mat                    
                    
#             for cell in range(self.ncells):
#                 for mutation in range(self.nmutations):
#                     total_reads = df_total_readcounts.values[cell, mutation]
#                     variant_reads = df_variant_readcounts.values[cell, mutation]
#                     if total_reads > 0:
#                         coeff = -np.log(total_reads) - betabinom.logpmf(variant_reads, total_reads, bb_alpha, bb_beta)
#                         coeff_mat[cell, mutation] = coeff
#             self.coeff_mat = coeff_mat
        else:
            self.coeff_mat = None
        
        # gurobi parameters
        self.threads = threads
        self.worklimit = timelimit
        self.verbose = verbose


        self.fpweight = np.log(1 - fp) - np.log(fp)
        self.fnweight = np.log(fn) - np.log(1 - fn)
        
        # solution
        # self.B = np.zeros((self.ncells, self.nmutations))
        self.solB = None
        self.solT_cell = None
        self.solT_mut = None

    def solveSetInclusion(self):

        ncells = self.ncells
        nmutations = self.nmutations
        nclusters = self.nclusters
        
        print(f'n = {ncells}, m = {nmutations}, p = {nclusters}')
        
        clustering = self.clustering
        k = self.k

        model = gp.Model('solveConstrainedDollo')

        # character matrix variables
        b = model.addVars(ncells, nmutations, vtype=gp.GRB.BINARY, name='b')
        c = model.addVars(nclusters, nmutations, k, vtype=gp.GRB.BINARY, name='c')

        # 1+ indicator
        x = model.addVars(ncells, nmutations, vtype=gp.GRB.CONTINUOUS, lb = 0, ub = 1, name='x')

        # row clustering consistency variables
        g = model.addVars(nclusters, nmutations, vtype=gp.GRB.CONTINUOUS, lb=0, ub=1, name='g')
        g0 = model.addVars(nclusters, nmutations, vtype=gp.GRB.CONTINUOUS, lb=0, ub=1, name='g0')
        g1 = model.addVars(nclusters, nmutations, vtype=gp.GRB.CONTINUOUS, lb=0, ub=1, name='g1')

        # column compatibility variables
        y0 = model.addVars(nmutations, nmutations, k, k, vtype=gp.GRB.CONTINUOUS, lb=0, ub=1, name='y0')
        y1 = model.addVars(nmutations, nmutations, k, vtype=gp.GRB.CONTINUOUS, lb=0, ub=1, name='y1')
        y2 = model.addVars(nmutations, nmutations, k, vtype=gp.GRB.CONTINUOUS, lb=0, ub=1, name='y2')
        y3 = model.addVars(nmutations, nmutations, vtype=gp.GRB.CONTINUOUS, lb=0, ub=1, name='y3')
        z0 = model.addVars(nmutations, nmutations, k, k, vtype=gp.GRB.CONTINUOUS, lb=0, ub=1, name='z0')
        z1 = model.addVars(nmutations, nmutations, k, vtype=gp.GRB.CONTINUOUS, lb=0, ub=1, name='z1')
        z2 = model.addVars(nmutations, nmutations, vtype=gp.GRB.CONTINUOUS, lb=0, ub=1, name='z2')

        # encode one-hot-like constraint on b and c
        for i in range(ncells):
            for j in range(nmutations):
                csum = gp.LinExpr()
                for s in range(k):
                    cluster = clustering[i]
                    csum += c[cluster,j,s]

                if self.mutation_list[j] not in self.snp_list:
                    model.addConstr(csum + b[i,j] <= 1)
                else:
                    model.addConstr(csum + b[i,j] == 1)
                model.addConstr(x[i,j] == b[i,j] + csum)
        
        # symmetry breaking
        for s in range(1,k):
            csum1 = gp.LinExpr()
            csum2 = gp.LinExpr()
            for l in range(nclusters):
                for j in range(nmutations):
                    csum1 += c[l,j,s]
                    csum2 += c[l,j,s-1]

            model.addConstr(csum2 >= csum1)        
        
        # encode consistency constraints
        ## g0 constraints
        for l in range(nclusters):
            for j in range(nmutations):
                xsum = gp.LinExpr()
                cluster_size = 0
                for i in np.where(clustering == l)[0]:
                    xsum += x[i,j]
                    model.addConstr(g0[l,j] >= 1 - x[i,j])
                    cluster_size += 1
                model.addConstr(g0[l,j] <= cluster_size - xsum)

        ## g1 constraints
        for l in range(nclusters):
            for j in range(nmutations):
                bsum = gp.LinExpr()
                for i in np.where(clustering == l)[0]:
                    bsum += b[i,j]
                    model.addConstr(g1[l,j] >= b[i,j])
                model.addConstr(g1[l,j] <= bsum)

        ## g constraints
        for l in range(nclusters):
            for j in range(nmutations):
                model.addConstr(g[l,j] <= g0[l,j])
                model.addConstr(g[l,j] <= g1[l,j])
                model.addConstr(g[l,j] >= g0[l,j] + g1[l,j] - 1)        
        
        for j in range(nmutations):
            gsum = gp.LinExpr()
            for l in range(nclusters):
                gsum += g[l,j]
            model.addConstr(gsum <= 1)        
        

        # encode set inclusion constraints
        ## y constraints (containment)
        for j1 in range(nmutations):
            for j2 in range(nmutations):
                if j1 == j2:
                    continue                
                for s1 in range(k):
                    for s2 in range(k):
                        for l in range(nclusters):
                            model.addConstr(y0[j1,j2,s1,s2] <= 1 - c[l,j2,s2] + c[l,j1,s1])

        for j1 in range(nmutations):
            for j2 in range(nmutations):
                if j1 == j2:
                    continue                
                for s2 in range(k):
                    for l in range(nclusters):
                        model.addConstr(y1[j1,j2,s2] <= 1 - c[l,j2,s2] + (1 - g0[l,j1]))

        for j1 in range(nmutations):
            for j2 in range(nmutations):
                if j1 == j2:
                    continue                
                for s1 in range(k):
                    for l in range(nclusters):
                        csum = gp.LinExpr()
                        for s in range(k):
                            csum += c[l,j2,s]
                        model.addConstr(y2[j1,j2,s1] <= 1 - (g1[l,j2] + csum) + c[l,j1,s1])
                        # model.addConstr(y2[j1,j2,s1] <= 1 - (1 - g0[l,j2]) + c[l,j1,s1])


        for j1 in range(nmutations):
            for j2 in range(nmutations):
                if j1 == j2:
                    continue
                for i in range(ncells):
                    model.addConstr(y3[j1,j2] <= 1 - x[i,j2] + x[i,j1])
        
        ## z constraints (disjointness)
        for j1 in range(nmutations):
            for j2 in range(nmutations):
                if j1 >= j2:
                    continue                
                for s1 in range(k):
                    for s2 in range(k):
                        for l in range(nclusters):
                            model.addConstr(z0[j1,j2,s1,s2] <= 2 - c[l,j1,s1] - c[l,j2,s2])
                            # model.addConstr(z1[j1,j2,s1,s2] <= 2 - c[l,j1,s1] - c[l,j2,s2] - (1 - g0[l,j1]))

                for i in range(ncells):
                    l = clustering[i]
                    model.addConstr(z2[j1,j2] <= 2 - x[i,j2] - x[i,j1])

        for j1 in range(nmutations):
            for j2 in range(nmutations):
                if j1 == j2:
                    continue                
                for s2 in range(k):
                    for l in range(nclusters):
                        csum = gp.LinExpr()
                        for s in range(k):
                            csum += c[l,j1,s]
                        model.addConstr(z1[j1,j2,s2] <= 2 - c[l,j2,s2] - (g1[l,j1] + csum))
                        # model.addConstr(z1[j1,j2,s2] <= 2 - c[l,j2,s2] - (1 - g0[l,j1]))
        
        ## avoid conflict
        for j1 in range(nmutations):
            for j2 in range(nmutations):
                if j1 >= j2:
                    continue                
                for s1 in range(k):
                    for s2 in range(k):
                        model.addConstr(y0[j1,j2,s1,s2] + y0[j2,j1,s2,s1] + z0[j1,j2,s1,s2] >= 1)
                        # model.addConstr(y1[j1,j2,s1,s2] + y1[j2,j1,s1,s2] + z1[j1,j2,s1,s2] >= 1)
                        # model.addConstr(y2[j1,j2,s1,s2] + y2[j2,j1,s1,s2] + z1[j1,j2,s1,s2] >= 1)
                model.addConstr(y3[j1,j2] + y3[j2,j1] + z2[j1,j2] >= 1)


        for j1 in range(nmutations):
            for j2 in range(nmutations):
                if j1 == j2:
                    continue                
                for s2 in range(k):
                    model.addConstr(y1[j1,j2,s2] + y2[j2,j1,s2] + z1[j1,j2,s2] >= 1)
        

        # set objective function
        obj_sum = gp.LinExpr()
        if self.coeff_mat is not None:
            for i in range(ncells):
                for j in range(nmutations):
                    if self.df_total_readcounts.values[i,j] > 0:
                        obj_sum += self.coeff_mat[i,j] * b[i,j]
        else:
            for i in range(ncells):
                for j in range(nmutations):
                    if self.df_character_matrix.values[i,j] == 0:
                        obj_sum += self.fnweight * b[i,j]
                    elif self.df_character_matrix.values[i,j] == 1:
                        obj_sum += self.fpweight * b[i,j]
            
        model.setObjective(obj_sum, gp.GRB.MAXIMIZE)
        
        model.setParam(gp.GRB.Param.Threads, self.threads)
        model.setParam(gp.GRB.Param.Method, 4)
        
        model.setParam(gp.GRB.Param.FeasibilityTol, 1e-6)
        model.setParam(gp.GRB.Param.IntFeasTol, 1e-6)
        model.setParam(gp.GRB.Param.OptimalityTol, 1e-6)
        
        # model.write('test_model_turing.lp')
        
        model.optimize()
        if model.status == gp.GRB.OPTIMAL:
            nzero_entries = np.sum(self.A == 0)
            none_entries = np.sum(self.A == 1)            
            opt_obj_value = model.getObjective().getValue()            
            print(f'{nzero_entries}, {none_entries}, {opt_obj_value}')
            #print(f'log likelihood: {opt_obj_value + nzero_entries * np.log(1 - self.fn) + none_entries * np.log(self.fp)}')

#             solb = np.zeros((ncells, nmutations))
            
#             for i in range(ncells):
#                 for j in range(nmutations):
#                     solb[i,j] = np.abs(model.getAttr('x', b)[i,j])
            solb = np.rint(np.reshape(model.getAttr('x', b).values(), (ncells, nmutations)))

            solc = model.getAttr('x', c)
            for l in range(nclusters):
                for j in range(nmutations):
                    for s in range(k):
                        if solc[l,j,s] > 0.5:
                            for i in np.where(clustering == l)[0]:
                                solb[i,j] = s + 2

            # df_solb = pd.DataFrame(solb, index = self.df_character_matrix.index,
            #                        columns = self.df_character_matrix.columns[:-1])
            
            # df_solb.to_csv('new_test_condor_output.csv')
            # df_solb.to_csv('newest_test_condor_output.csv')
            
            df_solb = pd.DataFrame(solb, index = self.df_character_matrix.index,
                                   columns = self.df_character_matrix.columns[:-1], dtype=int)
            # print(model.getAttr('x', b))
            # print('-'*50)
            # print(model.getAttr('x', g))
            # print('-'*50)
            # print(solc)
            
            self.solB = df_solb
            df_solb_binary = solveConstrainedDollo.expand_multi_state_to_binary(df_solb)
            self.solT_mut, self.solT_cell = solveConstrainedDollo.generate_perfect_phylogeny(df_solb_binary)
    
    def writeSolution(self, fname):
        if self.solB is not None:
            self.solB.to_csv(fname)
    
    @staticmethod
    def expand_multi_state_to_binary(df_multistate):

        ncells = len(df_multistate)
        binarized_mat = None
        binary_col_dict = {}
        for column in df_multistate.columns:
            max_state = df_multistate[column].max()
            for s in range(1, max_state+1):
                state_col = np.zeros((ncells))
                if s == 1:
                    state_col[df_multistate[column] > 0] = 1
                else:
                    state_col[df_multistate[column] == s] = 1

                binary_col_dict[f'{column}_{s}'] = state_col

        df_binary = pd.DataFrame(binary_col_dict, index = df_multistate.index, dtype=int)
        return df_binary    
    
    @staticmethod
    def generate_perfect_phylogeny(df_binary):

        solT_mut = nx.DiGraph()
        solT_mut.add_node('root')

        solT_cell = nx.DiGraph()
        solT_cell.add_node('root')

        df_binary = df_binary[df_binary.sum().sort_values(ascending=False).index]    

        for cell_id, row in df_binary.iterrows():
            if cell_id == 'root':
                continue

            curr_node = 'root'
            for column in df_binary.columns[row.values == 1]:
                if column in solT_mut[curr_node]:
                    curr_node = column
                else:
                    if column in solT_mut.nodes:
                        raise NameError(f'{column} is being repeated')
                    solT_mut.add_edge(curr_node, column)
                    solT_cell.add_edge(curr_node, column)
                    curr_node = column

            solT_cell.add_edge(curr_node, cell_id)   

        return solT_mut, solT_cell

    def writeDOT(self, dot_file, withcells=True):
        if withcells is True:
            writeTree = self.solT_cell
        else:
            writeTree = self.solT_mut
        
        with open(dot_file, 'w') as output:

            output.write(f'digraph N {{\n')
            output.write(f"\toverlap=\"false\"\n")
            output.write(f"\trankdir=\"TB\"\n")

            idx_dict = {}
            idx = 0
            if writeTree is not None:
                for node in writeTree.nodes:
                    idx_dict[node] = idx
                    output.write(f'\t{idx} [label=\"{node}\", style=\"bold\"];\n')
                    idx += 1

                for edge in writeTree.edges:
                    output.write(f"\t{idx_dict[edge[0]]} -> {idx_dict[edge[1]]} [style=\"bold\"];\n")

            output.write(f'}}')
