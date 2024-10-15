from typing import List, Any
import os
import numpy as np
import pandas as pd
import math
from ete3 import Tree
import random
from statistics import mean
from Bio import Phylo
from Bio.Phylo import draw_ascii
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio.Phylo import draw_ascii
import time
from itertools import groupby, product, filterfalse, starmap, islice, count
from scipy.optimize import fsolve

def onlyfirst(l):
    return l[0]

def onlysecond(l):
    return l[1]

def cogonly(l):
    return list(map(onlyfirst, l))



class TDDTR:
#_______________________________________________________________________________________________________________________
# Simulations parameters
#_______________________________________________________________________________________________________________________
    trials = 10                     # Number of repetitions
    mean_edge_length = 0.05         # Mean edge length. Each directional edge length is from exponential distribution
                                    # with that mean.
    leaves_number = 10              # Number of leaves (taxa) in the random tree
    n0 = 2000                       # Number of genes in at the root genome
    with_unimog = False             # Should be true only if UNIMOG is downloaded
#_______________________________________________________________________________________________________________________
    genome0 = list(range(n0))       # The root genome
    genomes = [genome0]             # List of genome; one for each node at the rooted tree.
    leaf_genomes = []               # The genomes that are represented by leaves. This is the inpyt for all methods.
    random_tree = [[0, 0, 0, -1]]   # The random tree; for each node, the serial number of the node, the serial number
                                    # of the childs and the serial number of the parent.
    leaves = [0]                    # The serial number of the leaves at the random tree
    labels = []                     # Labels for the newick tree
    test_tree = []                  # The ete tree of the random tree
    gene_serial = n0                # Each gene in a repetition have a serial number.
    node_serial_number = 1          # The serial numer of the next created node on the random tree
    rfd_unimog = 0                  # Sum the Robinson - Foulds distance for UNIMOG
    rfd_sgc = 0                     # Sum the Robinson - Foulds distance for SGC
    rfd_tddtr = 0                   # Sum the Robinson - Foulds distance for TDDTR
    rfd_max = 2 * leaves_number - 6 # Maximal possible Robinson - Foulds distance
    rf_factor = trials * rfd_max
    time_unimog =0                  # Mean time for UNIMOG each cycle
    time_sgc = 0                    # Mean the time for SCG each cycle
    time_tddtr = 0                  # Mean the time for TDDR each cycle
    mean_success_rate = 0           # Mean rate of success for TDDR
    dm = []                         # Distance matrix for UNIMOG and SGC

    leavesgnm = []                  # Set of COGs for each taxa. Without repetitions
    pairs = []                      # Set intersections betweens each pair of taxa
    log_leaves_size = []            # list of log of number of distinct COGS for each taxa

    nmax = 0                        # hold the number of leaves of the tree that we look for his cherry. Starts equal
                                    # to leaves_number and then reduce by one for each repetition of the main loop.
    new_node_number = 0             # Each node is given a serial number. This is the serial number of the current
                                    # new node.
    pair1 = 0                       # Leaf number of the cherry first leaf
    pair2 = 0                       # Leaf number of the cherry second leaf

    newick_tree = []                # Newick tree
    ete_tree = []                   # ete3 tree reference tree
    ete_tree_tddtr = []             # ete3 tree tddtr tree

    nodename = []                   # list of node serial number of the current leaves
    intersect_res = []              # List of to matrices: lpairs, ltriplets. These matrices contain the log of the size
                                    # of pairs intersections and triplets intersections
    leaf_out = []                   # For every leaf a list for all triplet with that leaf that containes distance
                                    # from the center and the other two leaves
    leaf_in = []                    # For every leaf a list for all triplet with that leaf that containes distance
                                    # to the center and the other two leaves
    mean_out = []                   #  For each leaf estimator for the distance from the closest node to the leaf
    mean_in = []                    #  For each leaf estimator for the distance to the closest node from the leaf

    pairs_dis = []                  # Structure that containes all directional distances for each pair

    pairs_0_in = []                 # For every pair of leaves m<n, all the triplets m,n,p the distance from m to the
                                    # center
    pairs_0_out = []                # For every pair of leaves m<n, all the triplets m,n,p the distance from center to m
    pairs_1_in = []                 # For every pair of leaves m<n, all the triplets m,n,p the distance from n to the
                                    # center
    pairs_1_out = []                # For every pair of leaves m<n, all the triplets m,n,p the distance from center to n
    pairs_2_in = []                 # For every pair of leaves m<n, all the triplets m,n,p the distance from p to the
                                    # center
    pairs_2_out = []                # For every pair of leaves m<n, all the triplets m,n,p the distance from center to p
    pairs_2_id = []                 # For every pair of leaves m<n, all the triplets m,n,p the node number of p
    constructed_tree_leaves = []    # Holds the true serial numbers of nodes in case the reconstructed tree is correct
                                    # For the check of the correctness of the reconstructed tree


    pairs_0_mean_in = []
    pairs_0_mean_out = []
    pairs_1_mean_in = []
    pairs_1_mean_out = []
    success = 1
    number_of_successes = 0
    errors = []                     # For each taxa the number of distances of triplets with that leaf with negative
                                    # value
    sum_errors = []                 # For each taxa the sum of negative values in triplets that containes this leaf.

    def input(self):

        self.nmax = self.leaves_number
        self.constructed_tree_leaves = list(self.leaves)
        return

    def intersect(self):
        # Computes intersection sets for all pair of taxa (pairs) and compute the log of the intersection size (lpairs).
        # Compute the log of the intersection size of triplets of taxa (ltriplets).
        self.leavesgnm = []
        #print("size:")
        for m in range(self.leaves_number):
            self.leavesgnm.append(set(self.leaf_genomes[m]))
            self.log_leaves_size.append(np.log(len(self.leavesgnm[m])))
            #print(m, len(self.leavesgnm[m]))
        self.pairs = []
        #print("pairs:")
        for m in range(self.leaves_number):
            self.pairs.append([])
            for n in range(self.leaves_number):
                self.pairs[m].append([])
        lpairs = np.zeros([self.leaves_number, self.leaves_number])

        for m in range(self.leaves_number - 1):
            for n in range(m + 1, self.leaves_number):
                self.pairs[m][n] = self.leavesgnm[m] & self.leavesgnm[n]
                lpairs[m][n] = np.log(len(self.leavesgnm[m] & self.leavesgnm[n]))
                #print(m,n,round(np.exp(lpairs[m][n]),0),lpairs[m][n])
        #print("triplets:")
        ltriplets = np.zeros([self.leaves_number, self.leaves_number, self.leaves_number])
        for m in range(self.leaves_number - 2):
            for n in range(m + 1, self.leaves_number - 1):
                for p in range(n + 1, self.leaves_number):
                    ltriplets[m][n][p] = np.log(len(self.pairs[m][n] & self.leavesgnm[p]))
                    #print(m,n,p,round(np.exp(ltriplets[m][n][p]),0),ltriplets[m][n][p])
        self.intersect_res = [lpairs, ltriplets]
        return
    def add_to_pairs_struct(self):
        # Prepairs structure
        self.pairs_0_in.append([])
        self.pairs_0_out.append([])
        self.pairs_1_in.append([])
        self.pairs_1_out.append([])
        self.pairs_2_in.append([])
        self.pairs_2_out.append([])
        self.pairs_2_id.append([])
        return
    def append_to_pairs_mn_data(self, m,n,p,dlin: object, dlout: object, drin: object, drout: object, dmin: object, dmout: object) -> object:
        #Creates lists of data for pairs
        self.pairs_0_in[m][n].append(dlin)
        self.pairs_0_out[m][n].append(dlout)
        self.pairs_1_in[m][n].append(drin)
        self.pairs_1_out[m][n].append(drout)
        self.pairs_2_in[m][n].append(dmin)
        self.pairs_2_out[m][n].append(dmout)
        self.pairs_2_id[m][n].append(self.nodename[p])
        return
    def initial_pairs_struct(self):
        #Prepairs structure
        self.errors = []
        self.sum_errors = []
        self.pairs_dis = []
        self.leaf_in = []
        self.leaf_out = []
        self.pairs_0_in = []
        self.pairs_0_out = []
        self.pairs_1_in = []
        self.pairs_1_out = []
        self.pairs_2_in = []
        self.pairs_2_out = []
        self.pairs_2_id = []
        self.pairs_0_mean_in = []
        self.pairs_0_mean_out = []
        self.pairs_1_mean_in = []
        self.pairs_1_mean_out = []
        for m in range(self.leaves_number):
            self.pairs_dis.append([])
            self.add_to_pairs_struct()
            self.leaf_in.append([])
            self.leaf_out.append([])
            self.errors.append(0)
            self.sum_errors.append(0)
            self.pairs_0_mean_in.append([])
            self.pairs_0_mean_out.append([])
            self.pairs_1_mean_in.append([])
            self.pairs_1_mean_out.append([])

            for n in range(self.leaves_number):
                self.pairs_dis[m].append([])
                self.pairs_0_in[m].append([])
                self.pairs_0_out[m].append([])
                self.pairs_1_in[m].append([])
                self.pairs_1_out[m].append([])
                self.pairs_2_in[m].append([])
                self.pairs_2_out[m].append([])
                self.pairs_2_id[m].append([])

                self.pairs_0_mean_in[m].append([])
                self.pairs_0_mean_out[m].append([])
                self.pairs_1_mean_in[m].append([])
                self.pairs_1_mean_out[m].append([])

        return
    def preparing_loop(self):
        #Prepairing lists for the loop

        self.newick_tree = list(self.labels)
        self.nodename = list(range(self.leaves_number))
        self.new_node_number = self.leaves_number
        return
    def initial_values(self):
    # Filling directional distances triplet of taxa
        for m in range(self.leaves_number - 2):
            for n in range(m + 1, self.leaves_number - 1):
                for p in range(n + 1, self.leaves_number):

                    dmout = self.intersect_res[0][n][p] - self.intersect_res[1][m][n][p]
                    dnout = self.intersect_res[0][m][p] - self.intersect_res[1][m][n][p]
                    dpout = self.intersect_res[0][m][n] - self.intersect_res[1][m][n][p]

                    dmin = (self.intersect_res[1][m][n][p] + self.log_leaves_size[m] - self.intersect_res[0][m][n]
                           - self.intersect_res[0][m][p])
                    dnin = (self.intersect_res[1][m][n][p] + self.log_leaves_size[n] - self.intersect_res[0][m][n]
                           - self.intersect_res[0][n][p])
                    dpin = (self.intersect_res[1][m][n][p] + self.log_leaves_size[p] - self.intersect_res[0][m][p]
                           - self.intersect_res[0][n][p])
                    #print(m,n,p,dmout,dnout,dpout,dmin,dnin,dpin)
                    if dmin < 0:
                        self.errors[m] += 1
                        self.errors[n] += 1
                        self.errors[p] += 1
                        self.sum_errors[m] += dmin
                        self.sum_errors[n] += dmin
                        self.sum_errors[p] += dmin
                        dmin = 0
                    if dnin < 0:
                        self.errors[m] += 1
                        self.errors[n] += 1
                        self.errors[p] += 1
                        self.sum_errors[m] += dnin
                        self.sum_errors[n] += dnin
                        self.sum_errors[p] += dnin
                        dnin = 0
                    if dpin < 0:
                        self.errors[m] += 1
                        self.errors[n] += 1
                        self.errors[p] += 1
                        self.sum_errors[m] += dpin
                        self.sum_errors[n] += dpin
                        self.sum_errors[p] += dpin
                        dpin = 0

                    self.leaf_out[m].append([dmout, n, p])
                    self.leaf_out[n].append([dnout, m, p])
                    self.leaf_out[p].append([dpout, m, n])

                    self.leaf_in[m].append([dmin, n, p])
                    self.leaf_in[n].append([dnin, m, p])
                    self.leaf_in[p].append([dpin, m, n])

                    self.pairs_dis[m][n].append([p, [dmout, dnout, dpout], [dmin, dnin, dpin]])
                    self.append_to_pairs_mn_data(m,n,p,dmin, dmout, dnin, dnout, dpin, dpout)

                    self.pairs_dis[m][p].append([n, [dmout, dpout, dnout], [dmin, dpin, dnin]])
                    self.append_to_pairs_mn_data(m,p,n, dmin, dmout, dpin, dpout, dnin, dnout)

                    self.pairs_dis[n][p].append([m, [dnout, dpout, dmout], [dnin, dpin, dmin]])
                    self.append_to_pairs_mn_data(n,p,m,dnin, dnout, dpin, dpout, dmin, dmout)
        return



    def choosing_cherry(self):
        #Choosing the fittes couple of leaves to be cherry in the reconstruct tree
        for m in range(self.nmax - 1):
            for n in range(self.nmax):
                self.pairs_dis[m][n] = sorted(self.pairs_dis[m][n], key=onlyfirst)
        self.mean_out = []
        self.mean_in = []
        for m in range(self.nmax):
            self.leaf_out[m] = sorted(self.leaf_out[m])
            smeanout = []
            for n in range(self.nmax - 2):
                smeanout.append(self.leaf_out[m][n][0])
            result = mean(smeanout)
            #print(m,"out:",result)
            self.mean_out.append(result)
            self.leaf_in[m] = sorted(self.leaf_in[m])
            smeanin = []
            for n in range(self.nmax - 2):
                smeanin.append(self.leaf_in[m][n][0])
            result = mean(smeanin)
            #print(m, "in:", result)
            self.mean_in.append(result)

        critlist = []
        for m in range(self.nmax - 1):
            for n in range(m + 1, self.nmax):
                #print("m",m,"n",n)
                self.pairs_0_mean_in[m][n] = mean(self.pairs_0_in[m][n])
                self.pairs_0_mean_out[m][n] = mean(self.pairs_0_out[m][n])
                self.pairs_1_mean_in[m][n] = mean(self.pairs_1_in[m][n])
                self.pairs_1_mean_out[m][n] = mean(self.pairs_1_out[m][n])
                #print("0in",self.pairs_0_mean_in[m][n])
                #print("0out", self.pairs_0_mean_out[m][n])
                #print("1in",self.pairs_1_mean_in[m][n])
                #print("1out", self.pairs_1_mean_out[m][n])
                pairs_d_0_in = self.pairs_0_mean_in[m][n] - self.mean_in[m]
                pairs_d_0_out = self.pairs_0_mean_out[m][n] - self.mean_out[m]
                pairs_d_1_in = self.pairs_1_mean_in[m][n] - self.mean_in[n]
                pairs_d_1_out = self.pairs_1_mean_out[m][n] - self.mean_out[n]
                #crit = np.max([pairs_d_0_out, pairs_d_1_out])
                crit = np.max([pairs_d_0_in, pairs_d_0_out, pairs_d_1_in, pairs_d_1_out])
                #print(m,n,"crit:",crit)
                critlist.append([crit, m, n])

        critlist = sorted(critlist)
        self.pair1 = critlist[0][1]
        self.pair2 = critlist[0][2]

        return

    def check(self):
        if self.success == 1:
            parent = self.random_tree[self.constructed_tree_leaves[self.pair1]][3]
            child1 = self.random_tree[self.constructed_tree_leaves[self.pair1]][1]
            child2 = self.random_tree[self.constructed_tree_leaves[self.pair1]][2]
            nghbrs1 = {parent}
            if child1 > 0:
                nghbrs1.add(child1)
            if child2 > 0:
                nghbrs1.add(child2)
            parent = self.random_tree[self.constructed_tree_leaves[self.pair2]][3]
            child1 = self.random_tree[self.constructed_tree_leaves[self.pair2]][1]
            child2 = self.random_tree[self.constructed_tree_leaves[self.pair2]][2]
            nghbrs2 = {parent}
            if child1 > 0:
                nghbrs2.add(child1)
            if child2 > 0:
                nghbrs2.add(child2)
            connect = nghbrs1 & nghbrs2
            connect_list = list(connect)
            if len(connect_list) == 1:
                self.number_of_successes += 1
                self.constructed_tree_leaves.pop(self.pair2)
                self.constructed_tree_leaves.pop(self.pair1)
                self.constructed_tree_leaves.append(connect_list[0])
            else:
                self.success = 0

            return

    def updating_trees_data(self):
        # Updating trees structures
        self.nodename.append(self.new_node_number)

        self.newick_tree.append("(" + self.newick_tree[self.pair1] + ":" + str(round(max([0, self.mean_in[self.pair1]
            + self.mean_out[self.pair1]]), 3)) + "," + self.newick_tree[ self.pair2] + ":"
            + str(round(max([0, self.mean_in[self.pair2] + self.mean_out[self.pair2]]), 3)) + ")")

        self.newick_tree.pop(self.pair2)
        self.newick_tree.pop(self.pair1)
        return

    def extending_struct_for_new_node(self):
        #Extending the structure for the new node data
        self.leaf_out.append([])
        self.leaf_in.append([])
        self.pairs_dis.append([])
        self.add_to_pairs_struct()

        for m in range(self.nmax + 1):
            self.pairs_dis[m].append([])
            self.pairs_0_in[m].append([])
            self.pairs_0_out[m].append([])
            self.pairs_1_in[m].append([])
            self.pairs_1_out[m].append([])
            self.pairs_2_in[m].append([])
            self.pairs_2_out[m].append([])
            self.pairs_2_id[m].append([])

        for m in range(self.nmax):
            self.pairs_dis[self.nmax].append([])
            self.pairs_0_in[self.nmax].append([])
            self.pairs_0_out[self.nmax].append([])
            self.pairs_1_in[self.nmax].append([])
            self.pairs_1_out[self.nmax].append([])
            self.pairs_2_in[self.nmax].append([])
            self.pairs_2_out[self.nmax].append([])
            self.pairs_2_id[self.nmax].append([])
        return
    def computing_data_for_new_node(self):
        #Computes the data for the new node and updating the structure
        for m in range(self.nmax - 1):
            if not (m in [self.pair1, self.pair2]):
                locm = m - (self.pair1 < m) - (self.pair2 < m)
                for n in range(m + 1, self.nmax):
                    if not (n in [self.pair1, self.pair2]):
                        locn = n - (self.pair1 < n) - (self.pair2 < n)
                        locpair1 = self.pair1 - (m < self.pair1) - (n < self.pair1)
                        locpair2 = self.pair2 - (m < self.pair2) - (n < self.pair2)

                        newtrpmout = (self.pairs_dis[m][n][locpair1][1][0]
                                      + self.pairs_dis[m][n][locpair2][1][0]) / 2
                        newtrpnout = (self.pairs_dis[m][n][locpair1][1][1]
                                      + self.pairs_dis[m][n][locpair2][1][1]) / 2

                        newtrpmin = (self.pairs_dis[m][n][locpair1][2][0]
                                     + self.pairs_dis[m][n][locpair2][2][0]) / 2
                        newtrpnin = (self.pairs_dis[m][n][locpair1][2][1]
                                     + self.pairs_dis[m][n][locpair2][2][1]) / 2
                        #print("newtrpmout",newtrpmout, "newtrpnout",newtrpnout, "newtrpmin",newtrpmin, "newtrpnin",newtrpnin)
                        newnodeout = (
                                             (self.pairs_dis[m][n][locpair1][1][2] - self.mean_out[self.pair1]) +
                                             (self.pairs_dis[m][n][locpair2][1][2] - self.mean_out[self.pair2]) +
                                             (self.pairs_dis[self.pair1][self.pair2][locm][2][2] - newtrpmin) +
                                             (self.pairs_dis[self.pair1][self.pair2][locn][2][2] - newtrpnin)) / 4

                        newnodein = (
                                            (self.pairs_dis[m][n][locpair1][2][2] - self.mean_in[self.pair1]) +
                                            (self.pairs_dis[m][n][locpair2][2][2] - self.mean_in[self.pair2]) +
                                            (self.pairs_dis[self.pair1][self.pair2][locm][1][2] - newtrpmout) +
                                            (self.pairs_dis[self.pair1][self.pair2][locn][1][2] - newtrpnout)) / 4
                        #print("newnodeout",newnodeout, "newnodein",newnodein)
                        """
                        newnodeout = (
                                             (self.pairs_dis[m][n][locpair1][1][2] - self.pairs_0_mean_out[self.pair1][self.pair2]) +
                                             (self.pairs_dis[m][n][locpair2][1][2] - self.pairs_1_mean_out[self.pair1][self.pair2]) +
                                             (self.pairs_dis[self.pair1][self.pair2][locm][2][2] - newtrpmin) +
                                             (self.pairs_dis[self.pair1][self.pair2][locn][2][2] - newtrpnin)) / 4

                        newnodein = (
                                            (self.pairs_dis[m][n][locpair1][2][2] - self.pairs_0_mean_in[self.pair1][self.pair2]) +
                                            (self.pairs_dis[m][n][locpair2][2][2] - self.pairs_1_mean_in[self.pair1][self.pair2]) +
                                            (self.pairs_dis[self.pair1][self.pair2][locm][1][2] - newtrpmout) +
                                            (self.pairs_dis[self.pair1][self.pair2][locn][1][2] - newtrpnout)) / 4
                        """
                        self.leaf_out[m].append([newtrpmout, self.nodename[n], self.nodename[self.nmax]])
                        self.leaf_out[n].append([newtrpnout, self.nodename[m], self.nodename[self.nmax]])
                        self.leaf_out[self.nmax].append([newnodeout, self.nodename[m], self.nodename[n]])

                        self.leaf_in[m].append([newtrpmin, self.nodename[n], self.nodename[self.nmax]])
                        self.leaf_in[n].append([newtrpnin, self.nodename[m], self.nodename[self.nmax]])
                        self.leaf_in[self.nmax].append([newnodein, self.nodename[m], self.nodename[n]])

                        self.pairs_dis[m][n].append(
                            [self.nodename[self.nmax], [newtrpmout, newtrpnout, newnodeout],
                            [newtrpmin, newtrpnin, newnodein]])
                        self.append_to_pairs_mn_data(m, n, self.nmax,
                            newtrpmin, newtrpmout, newtrpnin, newtrpnout, newnodein, newnodeout)

                        self.pairs_dis[m][self.nmax].append(
                            [self.nodename[n], [newtrpmout, newnodeout, newtrpnout],
                            [newtrpmin, newnodein, newtrpnin]])
                        self.append_to_pairs_mn_data(m, self.nmax,n ,
                            newtrpmin, newtrpmout,newnodein, newnodeout, newtrpnin, newtrpnout )

                        self.pairs_dis[n][self.nmax].append(
                            [self.nodename[m], [newtrpnout, newnodeout, newtrpmout],
                            [newtrpnin, newnodein, newtrpmin]])
                        self.append_to_pairs_mn_data(n,self.nmax, m,
                            newtrpnin, newtrpnout, newnodein, newnodeout,newtrpmin, newtrpmout)
        return

    def poping_pair1_pair2(self):
        #Deleting the data of the cherry leaves from the structures

        lg = (self.nmax - 1) * (self.nmax - 2) / 2
        pairset = {self.nodename[self.pair1], self.nodename[self.pair2]}
        for m in range(self.nmax):
            if not (m == self.pair1 or m == self.pair2):
                lgtemp = lg
                n = 0
                while n < lgtemp:
                    if len({self.leaf_out[m][n][1], self.leaf_out[m][n][2]} & pairset) > 0:
                        self.leaf_out[m].pop(n)
                        lgtemp += -1
                    else:
                        n += 1
                lgtemp = lg
                n = 0
                while n < lgtemp:
                    if len({self.leaf_in[m][n][1], self.leaf_in[m][n][2]} & pairset) > 0:
                        self.leaf_in[m].pop(n)
                        lgtemp += -1
                    else:
                        n += 1
        for m in range(self.nmax - 1):
            if not (m in [self.pair1, self.pair2]):
                for n in range(m + 1, self.nmax):
                    if not (n in [self.pair1, self.pair2]):
                        locpair1 = self.pair1 - (m < self.pair1) - (n < self.pair1)
                        locpair2 = self.pair2 - (m < self.pair2) - (n < self.pair2)

                        self.pairs_dis[m][n].pop(locpair2)
                        self.pairs_dis[m][n].pop(locpair1)

        for m in range(self.nmax - 1):
            if not (m in [self.pair1, self.pair2]):
                for n in range(m + 1, self.nmax):
                    ltemp = self.nmax - 2
                    if not (n in [self.pair1, self.pair2]):
                        p = 0
                        while p < ltemp:
                            if self.pairs_2_id[m][n][p] in pairset:

                                self.pairs_0_in[m][n].pop(p)
                                self.pairs_0_out[m][n].pop(p)
                                self.pairs_1_in[m][n].pop(p)
                                self.pairs_1_out[m][n].pop(p)
                                self.pairs_2_in[m][n].pop(p)
                                self.pairs_2_out[m][n].pop(p)
                                self.pairs_2_id[m][n].pop(p)
                                ltemp += -1
                            else:
                                p += 1

        for m in range(self.nmax):
            if not (self.pair2 == m):
                self.pairs_dis[m].pop(self.pair2)
                self.pairs_0_in[m].pop(self.pair2)
                self.pairs_0_out[m].pop(self.pair2)
                self.pairs_1_in[m].pop(self.pair2)
                self.pairs_1_out[m].pop(self.pair2)
                self.pairs_2_in[m].pop(self.pair2)
                self.pairs_2_out[m].pop(self.pair2)
                self.pairs_2_id[m].pop(self.pair2)

            if not (self.pair1 == m):
                self.pairs_dis[m].pop(self.pair1)
                self.pairs_0_in[m].pop(self.pair1)
                self.pairs_0_out[m].pop(self.pair1)
                self.pairs_1_in[m].pop(self.pair1)
                self.pairs_1_out[m].pop(self.pair1)
                self.pairs_2_in[m].pop(self.pair1)
                self.pairs_2_out[m].pop(self.pair1)
                self.pairs_2_id[m].pop(self.pair1)

        self.pairs_dis.pop(self.pair2)
        self.pairs_0_in.pop(self.pair2)
        self.pairs_0_out.pop(self.pair2)
        self.pairs_1_in.pop(self.pair2)
        self.pairs_1_out.pop(self.pair2)
        self.pairs_2_in.pop(self.pair2)
        self.pairs_2_out.pop(self.pair2)
        self.pairs_2_id.pop(self.pair2)

        self.pairs_dis.pop(self.pair1)
        self.pairs_0_in.pop(self.pair1)
        self.pairs_0_out.pop(self.pair1)
        self.pairs_1_in.pop(self.pair1)
        self.pairs_1_out.pop(self.pair1)
        self.pairs_2_in.pop(self.pair1)
        self.pairs_2_out.pop(self.pair1)
        self.pairs_2_id.pop(self.pair1)

        self.leaf_out.pop(self.pair2)
        self.leaf_out.pop(self.pair1)
        self.leaf_in.pop(self.pair2)
        self.leaf_in.pop(self.pair1)

        self.nodename.pop(self.pair2)
        self.nodename.pop(self.pair1)
        return

    def updating_data(self):
        #Manages all the updating process
        self.updating_trees_data()
        self.extending_struct_for_new_node()
        self.computing_data_for_new_node()

        if self.nmax > 4:
            self.poping_pair1_pair2()
            self.new_node_number += 1
        return

    def program_end(self):
        self.new_node_number += 1
        self.nodename.append(self.new_node_number)

        self.newick_tree.append("(" + self.newick_tree[0] + ":" + str(round(max([0, self.mean_in[0]
            + self.mean_out[0]]), 3)) + "," + self.newick_tree[1] + ":" + str(round(max([0, self.mean_in[1]
            + self.mean_out[1]]), 3)) + "," + self.newick_tree[2] + ":" + str(round(max([0, self.mean_in[2]
            + self.mean_out[2]]), 3)) + ")" + ";")
        self.ete_tree_tddtr =Tree(self.newick_tree[3])
        return

    def tddtr_prog(self):
        self.input()
        self.preparing_loop()
        self.initial_pairs_struct()
        self.intersect()
        self.initial_values()


        while self.nmax > 3:
            self.choosing_cherry()
            self.check()
            self.updating_data()
            self.nmax += -1
        if self.success == 1:
            self.mean_success_rate += 1
        self.program_end()
        return

    def edge_process(self, directed_distances_pair, parent_genome):
        temp_genome = list(parent_genome)
        n0 = len(parent_genome)
        d1 = directed_distances_pair[0]
        d2 = directed_distances_pair[1]
        #print("d1",d1,"d2",d2)
        sum_et = 0
        nc = n0
        sd = d1 + d2
        while sum_et < sd:
            et = 1 / nc
            sum_et += np.random.exponential(scale=et, size=1)[0]
            lgnm = len(temp_genome)
            rgen = random.randrange(lgnm)
            if np.random.uniform(0, 1) > d1 / sd:
                temp_genome.pop(rgen)
                nc += -1
            else:
                rins = rgen + random.choice([0, 1])
                temp_genome.insert(rins, self.gene_serial)
                self.gene_serial += 1
                nc += 1
        return list(temp_genome)

    def create_child(self,second,random_new_parent):
        parent_genome = list(self.genomes[random_new_parent])
        edge_directed_distances: list[Any] = list(np.random.exponential(scale=self.mean_edge_length, size=2))
        child_genome = self.edge_process(edge_directed_distances, parent_genome)
        self.genomes.append(list(child_genome))
        self.random_tree[random_new_parent][1+second] = self.node_serial_number

        if self.node_serial_number < 3:
            self.random_tree.append([self.node_serial_number, 0, 0, 2 - second])
        else:
            self.random_tree.append([self.node_serial_number , 0, 0, random_new_parent])
        self.leaves.append(self.node_serial_number)
        self.node_serial_number +=1
        return

    def create_random_tree(self):

        current_number_of_leaves = 1
        while current_number_of_leaves < self.leaves_number: #construct rando tree with nmax leaves
        #Eliminating parent from leaves list
            loc = random.randrange(current_number_of_leaves)
            random_new_parent = self.leaves[loc]
            self.leaves.pop(loc)

            # creating two random childs and adding them to leaves list
            second = 0
            self.create_child(second, random_new_parent)
            second = 1
            self.create_child(second, random_new_parent)
            current_number_of_leaves = len(self.leaves)
        self.leaf_genomes = []

        for m in range(self.leaves_number):
            self.leaf_genomes.append(self.genomes[self.leaves[m]])
        #print(self.random_tree)
        return

    #Constructing newick string for the created random tree
    def newick(self ):

        self.labels = []
        for m in range(self.leaves_number):
            self.labels.append("L" + str(m))
        newick_tree = list(self.labels)
        node_serials = list(range(self.leaves_number))
        new_node_serial = self.leaves_number
        constructed_nodes = list(self.leaves)
        components_number = len(node_serials)
        while components_number  > 3:
            pair_found = 0
            for m in range(components_number - 1):
                for n in range(m + 1 , components_number):
                    parent = self.random_tree[constructed_nodes[node_serials[m]]][3]
                    child1 = self.random_tree[constructed_nodes[node_serials[m]]][1]
                    child2 = self.random_tree[constructed_nodes[node_serials[m]]][2]
                    nghbrs1 = {parent}
                    if child1 > 0:
                        nghbrs1.add(child1)
                    if child2 > 0:
                        nghbrs1.add(child2)
                    parent = self.random_tree[constructed_nodes[node_serials[n]]][3]
                    child1 = self.random_tree[constructed_nodes[node_serials[n]]][1]
                    child2 = self.random_tree[constructed_nodes[node_serials[n]]][2]

                    nghbrs2 = {parent}
                    if child1 > 0:
                        nghbrs2.add(child1)
                    if child2 > 0:
                        nghbrs2.add(child2)
                    connect = nghbrs1 & nghbrs2
                    connect_list = list(connect)

                    if len(connect_list) == 1:
                        if not(connect_list[0] in constructed_nodes):
                            constructed_nodes.append(connect_list[0])
                            pair_found = 1
                            newick_tree.append("(" + newick_tree[m] + "," + newick_tree[n]+")")
                            newick_tree.pop(n)
                            newick_tree.pop(m)
                            node_serials.pop(n)
                            node_serials.pop(m)
                            node_serials.append(new_node_serial)
                            new_node_serial += 1
                            break
                if pair_found == 1:
                    break
            if pair_found == 0:
                print("error")
                exit()
            components_number = len(node_serials)
        newick_tree_string = "(" + newick_tree[0] + "," + newick_tree[1] + "," + newick_tree[2] + ");"
        self.test_tree =Tree(newick_tree_string)
        return

    def dm_to_tree(self):
        m = DistanceMatrix(list(self.labels), list(self.dm))
        constructor = DistanceTreeConstructor()
        njtree = constructor.nj(m)
        treedata = "treedata.txt"
        Phylo.write(njtree, treedata, "newick")
        f = open(treedata, "r")
        newicktr = f.read()
        tree = Tree(newicktr, format=1)
        return tree

    def unimog_createsmatrix(self):
        try:
            os.remove('dcj_matrix.txt')
        except:
            print("no file")
        finally:
            try:
                os.remove('UNIMOGfile.txt')
            except:
                print("no file")
            finally:
                finp = open('UNIMOGfile.txt', 'w')
                for m in range(len(self.labels)):
                    finp.write("a>" + self.labels[m] + "\n")
                    for COG in self.leaf_genomes[m]:
                        finp.write(str(COG) + "\n")
                finp.close()
                fmat = open('dcj_matrix.txt', 'w')
                fmat.close()
                os.system("java -jar UniMoG.jar -m=6 -d UNIMOGfile.txt >>dcj_matrix.txt")

                f = open("dcj_matrix.txt", "r")
                k = 0
                while k < 3:
                    line = f.readline()
                    lst = line.split()
                    if len(lst) == 1:
                        k += 1

                self.dm= []
                for m in range(self.leaves_number):
                    self.dm.append([])
                for m in range(1, self.leaves_number):
                    line = f.readline()
                    line = line[10:]
                    self.dm[m] = list(map(int, line.split()))
                for m in range(self.leaves_number):
                    self.dm[m].append(0)
        return

    def unimog(self):
        self.unimog_createsmatrix()
        ete_unimog_tree = self.dm_to_tree()
        d = self.test_tree.robinson_foulds(ete_unimog_tree, unrooted_trees=True)[0]
        self.rfd_unimog += d
        return

    def scg_createsmatrix(self):
        self.dm = []
        cogsets = []
        for m in range(self.leaves_number):
            cogsets.append([])
            self.dm.append([])
            cogsets[m] = set(self.leaf_genomes[m])

        for m in range(self.leaves_number - 1):
            le0 = len(cogsets[m])
            for n in range(m + 1, self.leaves_number):
                le1 = len(cogsets[n])
                mn_intersection = cogsets[m].intersection(cogsets[n])
                mn_intersection_size = len(mn_intersection)
                gene_content = -np.log(mn_intersection_size / min(le0, le1))
                self.dm[n].append(gene_content)
        for m in range(self.leaves_number):
            self.dm[m].append(0)
        return

    def sgc(self):
        self.scg_createsmatrix()
        ete_sgc_tree = self.dm_to_tree()
        d = self.test_tree.robinson_foulds(ete_sgc_tree, unrooted_trees=True)[0]
        self.rfd_sgc += d
        return

    def tddtr(self):

        self.tddtr_prog()
        d = self.test_tree.robinson_foulds(self.ete_tree_tddtr,unrooted_trees=True)[0]
        self.rfd_tddtr += d
        return
    def init_repetitions(self):

        self.genomes = [self.genome0]
        self.leaf_genomes = []
        self.random_tree = [[0, 0, 0, -1]]
        self.leaves = [0]
        self.labels = []
        self.test_tree = []
        self.gene_serial = self.n0
        self.node_serial_number = 1
        self.success = 1

        return
    def repetitions(self):
        
        for i in range(self.trials):
            self.init_repetitions()
            self.create_random_tree()

            self.newick()
            time_start = time.time()
            self.tddtr()
            time_end = time.time()
            self.time_tddtr += time_end - time_start
            if self.with_unimog:
                time_start = time.time()
                self.unimog()
                time_end = time.time()
                self.time_unimog += time_end - time_start
            time_start = time.time()
            self.sgc()
            time_end = time.time()
            self.time_sgc += time_end - time_start
            if self.with_unimog:
                print(str(i+1)+":","unimog rfd:", round(self.rfd_unimog / (i+1) / self.rfd_max,5), "scg rfd:",
                      round(self.rfd_sgc / (i+1) / self.rfd_max,5) , "tddtr rfd:",
                      round(self.rfd_tddtr  / (i+1) / self.rfd_max,5))
            else:
                print(str(i+1)+":",
                      round(self.rfd_sgc / (i+1) / self.rfd_max,5) , "tddtr rfd:",
                      round(self.rfd_tddtr  / (i+1) / self.rfd_max,5))                
                print("number of successes:", self.mean_success_rate)
                
        self.mean_success_rate = self.mean_success_rate / self.trials

        self.time_unimog =  self.time_unimog/ self.trials
        self.time_sgc =  self.time_sgc / self.trials
        self.time_tddtr =  self.time_tddtr / self.trials

        self.rfd_unimog =  self.rfd_unimog / self.rf_factor
        self.rfd_sgc = self.rfd_sgc / self.rf_factor
        self.rfd_tddtr = self.rfd_tddtr / self.rf_factor

        print("trials:",self.trials,"mean edge length:", self.mean_edge_length)
        print("leaves number", self.leaves_number, "genome length", self.n0 )
        if self.with_unimog:
            print("unimog rfd:",self.rfd_unimog,"scg rfd:",self.rfd_sgc,"tddtr rfd:",self.rfd_tddtr)
            print("unimog time:", self.time_unimog, "scg rfd:", self.time_sgc, "tddtr rfd:", self.time_tddtr)
        else:
            print("scg rfd:",self.rfd_sgc,"tddtr rfd:",self.rfd_tddtr)
            print("scg rfd:", self.time_sgc, "tddtr rfd:", self.time_tddtr)            
        print("mean success rate:", self.mean_success_rate)
        return
tddtr_obj = TDDTR()
tddtr_obj.repetitions()




