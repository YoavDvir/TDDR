import numpy as np
import pandas as pd
import math
import time
import random
from ete3 import Tree
from statistics import mean
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio.Phylo import draw_ascii
from itertools import groupby, product, filterfalse, starmap, islice, count
from scipy.optimize import fsolve
import argparse

def compare_tddr():
    tddtr_obj = TDDTR()

    parser = argparse.ArgumentParser(description ='run the TDDR on ATGCs real data from command line.')
    parser.add_argument('atgcnum',
                        metavar ='atgcnum',
                        type = int,
                        nargs = 1,
                        help ='serial number of the ATGC')
    parser.add_argument(dest ='tddtr',
                        action ='store_const',
                        const = tddtr_obj.tddtr_prog,
                        help ='compute the parameters theta and alpha')

    args = parser.parse_args()
    return args.tddtr(args.atgcnum[0])


def onlyfirst(l):
    return l[0]

def onlysecond(l):
    return l[1]

def cogonly(l):
    return list(map(onlyfirst, l))

class TDDTR:
    atgcnum = "9"               # The ATGC family number
    taxa_names = []             # Names of TAXAs
    cogs_lists = []             # List of COGS for each taxa. With repetition
    edgelength = []             # Hold a tree structure that includes edge directed distances and approximation of number                                # of COGs for internal nodes
    curnodelength = []          # list of number of distinct cogs for each node. Approximated for internal nodes.
    leavesgnm = []              # Set of COGs for each taxa. Without repetitions
    pairs = []                  # Set intersections betweens each pair of taxa
    log_leaves_size = []        # list of log of number of distinct COGS for each taxa
    nleaves = 0                 # number of taxas in the family
    nmax = 0                    # hold the number of leaves of the tree that we look for his cherry. Starts equal                                # to nleaves and then reduce by one for each repetition of the main loop.
    new_node_number = 0         # Each node is given a serial number. This is the serial number of the current new node.
    pair1 = 0                   # Leaf number of the cherry first leaf
    pair2 = 0                   # Leaf number of the cherry second leaf

    leaves_short_labels = []    # List of short labels for the taxa
    short_labels_tree = []      # Newick tree with short labels
    newick_tree = []            # Newick tree with full labels

    nodename = []               # list of node serial number of the current leaves
    intersect_res = []          # List of to matrices: lpairs, ltriplets. These matrices contain the log of the size
                                # of pairs intersections and triplets intersections
    leaf_out = []               # For every leaf a list for all triplet with that leaf that containes distance
                                # from the center and the other two leaves
    leaf_in = []                # For every leaf a list for all triplet with that leaf that containes distance
                                # to the center and the other two leaves
    mean_out = []               #  For each leaf estimator for the distance from the closest node to the leaf
    mean_in = []                #  For each leaf estimator for the distance to the closest node from the leaf

    pairs_dis = []              # Structure that containes all directional distances for each pair

    pairs_0_in = []             # For every pair of leaves m<n, all the triplets m,n,p the distance from m to the center
    pairs_0_out = []            # For every pair of leaves m<n, all the triplets m,n,p the distance from center to m
    pairs_1_in = []             # For every pair of leaves m<n, all the triplets m,n,p the distance from n to the center
    pairs_1_out = []            # For every pair of leaves m<n, all the triplets m,n,p the distance from center to n
    pairs_2_in = []             # For every pair of leaves m<n, all the triplets m,n,p the distance from p to the center
    pairs_2_out = []            # For every pair of leaves m<n, all the triplets m,n,p the distance from center to p
    pairs_2_id = []             # For every pair of leaves m<n, all the triplets m,n,p the node number of p


    err = []                    # For each taxa the number of distances of triplets with that leaf with negative value
    sumerr = []                 # For each taxa the sum of negative values in triplets that containes this leaf.

    def readatgc(self):
        #Reads the proper atgc file and fill the data to taxa_names, cogs_lists and nleaves
        atgcfile = "ATGC" + str(self.atgcnum) + "reduced.csv"
        df = pd.read_csv(atgcfile, header=None)
        atgcf = df.to_numpy()
        groups = []
        self.taxa_names = []

        atgcf = sorted(atgcf, key=onlysecond)

        for k, g in groupby(atgcf, onlysecond):
            groups.append(list(g))  # Store group iterator as a list
            self.taxa_names.append(k)
        self.cogs_lists = list(map(cogonly, groups))
        print('Taxa names:')
        print()
        for i in range(len(self.taxa_names)):
            print("L"+str(i) +": "+ self.taxa_names[i])
        self.nleaves = len(self.taxa_names)
        self.nmax = self.nleaves
        return

    def intersect(self):
        # Computes intersection sets for all pair of taxa (pairs) and compute the log of the intersection size (lpairs).
        # Compute the log of the intersection size of triplets of taxa (ltriplets).
        for m in range(self.nleaves):
            self.leavesgnm.append(set(self.cogs_lists[m]))
            self.log_leaves_size.append(np.log(len(self.leavesgnm[m])))
            self.curnodelength = list(self.log_leaves_size)

        for m in range(self.nleaves):
            self.pairs.append([])
            for n in range(self.nleaves):
                self.pairs[m].append([])
        lpairs = np.zeros([self.nleaves, self.nleaves])
        for m in range(self.nleaves - 1):
            for n in range(m + 1, self.nleaves):
                self.pairs[m][n] = self.leavesgnm[m] & self.leavesgnm[n]
                lpairs[m][n] = np.log(len(self.leavesgnm[m] & self.leavesgnm[n]))

        ltriplets = np.zeros([self.nleaves, self.nleaves, self.nleaves])
        for m in range(self.nleaves - 2):
            for n in range(m + 1, self.nleaves - 1):
                for p in range(n + 1, self.nleaves):
                    ltriplets[m][n][p] = np.log(len(self.pairs[m][n] & self.leavesgnm[p]))
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
        for m in range(self.nleaves):
            self.pairs_dis.append([])
            self.add_to_pairs_struct()
            self.leaf_in.append([])
            self.leaf_out.append([])
            self.err.append(0)
            self.sumerr.append(0)

            for n in range(self.nleaves):
                self.pairs_dis[m].append([])
                self.pairs_0_in[m].append([])
                self.pairs_0_out[m].append([])
                self.pairs_1_in[m].append([])
                self.pairs_1_out[m].append([])
                self.pairs_2_in[m].append([])
                self.pairs_2_out[m].append([])
                self.pairs_2_id[m].append([])
        return
    def preparing_loop(self):
        #Prepairing lists for the loop
        for m in range(self.nleaves):
            self.leaves_short_labels.append("L" + str(m))

        self.short_labels_tree = list(self.leaves_short_labels)
        self.newick_tree = list(self.taxa_names)
        self.nodename = list(range(self.nleaves))
        self.new_node_number = self.nleaves
        return
    def initial_values(self):
    # Filling directional distances triplet of taxa
        for m in range(self.nleaves - 2):
            for n in range(m + 1, self.nleaves - 1):
                for p in range(n + 1, self.nleaves):

                    dmout = self.intersect_res[0][n][p] - self.intersect_res[1][m][n][p]
                    dnout = self.intersect_res[0][m][p] - self.intersect_res[1][m][n][p]
                    dpout = self.intersect_res[0][m][n] - self.intersect_res[1][m][n][p]

                    dmin = (self.intersect_res[1][m][n][p] + self.log_leaves_size[m] - self.intersect_res[0][m][n]
                           - self.intersect_res[0][m][p])
                    dnin = (self.intersect_res[1][m][n][p] + self.log_leaves_size[n] - self.intersect_res[0][m][n]
                           - self.intersect_res[0][n][p])
                    dpin = (self.intersect_res[1][m][n][p] + self.log_leaves_size[p] - self.intersect_res[0][m][p]
                           - self.intersect_res[0][n][p])
                    if dmin < 0:
                        self.err[m] += 1
                        self.err[n] += 1
                        self.err[p] += 1
                        self.sumerr[m] += dmin
                        self.sumerr[n] += dmin
                        self.sumerr[p] += dmin
                        dmin = 0
                    if dnin < 0:
                        self.err[m] += 1
                        self.err[n] += 1
                        self.err[p] += 1
                        self.sumerr[m] += dnin
                        self.sumerr[n] += dnin
                        self.sumerr[p] += dnin
                        dnin = 0
                    if dpin < 0:
                        self.err[m] += 1
                        self.err[n] += 1
                        self.err[p] += 1
                        self.sumerr[m] += dpin
                        self.sumerr[n] += dpin
                        self.sumerr[p] += dpin
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
    def initial_direct_dist(self):
        #Initialize the direct distances structures
        self.initial_pairs_struct()
        self.initial_values()

        #print("sumerr", self.sumerr) #error reports
        #print("err", self.err)
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
            self.mean_out.append(result)

            self.leaf_in[m] = sorted(self.leaf_in[m])
            smeanin = []
            for n in range(self.nmax - 2):
                smeanin.append(self.leaf_in[m][n][0])
            result = mean(smeanin)
            self.mean_in.append(result)

        critlist = []
        for m in range(self.nmax - 1):
            for n in range(m + 1, self.nmax):
                pairslmeanin = mean(self.pairs_0_in[m][n]) - self.mean_in[m]
                pairslmeanout = mean(self.pairs_0_out[m][n]) - self.mean_out[m]
                pairsrmeanin = mean(self.pairs_1_in[m][n]) - self.mean_in[n]
                pairsrmeanout = mean(self.pairs_1_out[m][n]) - self.mean_out[n]
                crit = np.max([pairslmeanin, pairslmeanout, pairsrmeanin, pairsrmeanout])
                critlist.append([crit, m, n])

        critlist = sorted(critlist)
        self.pair1 = critlist[0][1]
        self.pair2 = critlist[0][2]

        return

    def updating_trees_data(self):
        # Updating trees structures
        self.nodename.append(self.new_node_number)

        self.curnodelength.append(( self.curnodelength[self.pair1] + self.mean_in[self.pair1] - self.mean_out[
            self.pair1] + self.curnodelength[self.pair2] + self.mean_in[self.pair2] - self.mean_out[self.pair2]) / 2)

        self.edgelength.append([
            "N" + str(self.new_node_number) + ".1", "dis->:"+ str(round(self.mean_in[self.pair1], 3))
            , "dis<-:"+str(round(self.mean_out[self.pair1], 3)),
            "size:"+str(round(np.exp(self.curnodelength[self.pair1]), 1)),
            self.short_labels_tree[self.pair1]])

        self.edgelength.append(["N" + str(self.new_node_number) + ".2","dis->:"+str(round(self.mean_in[self.pair2], 3)),
            "dis<-:"+str(round(self.mean_out[self.pair2], 3)),
            "size:"+str(round(np.exp(self.curnodelength[self.pair2]), 1)),
            self.short_labels_tree[self.pair2]])

        self.newick_tree.append("(" + self.newick_tree[self.pair1] + ":" + str(round(max([0, self.mean_in[self.pair1]
            + self.mean_out[self.pair1]]), 3)) + "," + self.newick_tree[ self.pair2] + ":"
            + str(round(max([0, self.mean_in[self.pair2] + self.mean_out[self.pair2]]), 3)) + ")")

        self.short_labels_tree.append("(" + self.short_labels_tree[self.pair1] + ":"
                + str(round(max([0, self.mean_in[self.pair1] + self.mean_out[self.pair1]]),3)) + ","
                + self.short_labels_tree[self.pair2] + ":" + str(round(max([0, self.mean_in[self.pair2]
                + self.mean_out[self.pair2]]), 3)) + ")" + "N" + str(self.new_node_number))

        self.short_labels_tree.pop(self.pair2)
        self.short_labels_tree.pop(self.pair1)
        self.newick_tree.pop(self.pair2)
        self.newick_tree.pop(self.pair1)
        self.curnodelength.pop(self.pair2)
        self.curnodelength.pop(self.pair1)
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

        self.curnodelength.append((self.curnodelength[0] + self.mean_in[0] - self.mean_out[0] + self.curnodelength[1]
                                   + self.mean_in[1] - self.mean_out[1] + self.curnodelength[2] + self.mean_in[2] -
                                   self.mean_out[2]) / 3)

        self.edgelength.append(["N" + str(self.new_node_number) + ".0", "dis->:"+str(round(self.mean_in[0], 3)),
            "dis<-:"+str(round(self.mean_out[0], 3)), "size:"+str(round(np.exp(self.curnodelength[0]), 1)),
                                self.short_labels_tree[0]])
        self.edgelength.append(["N" + str(self.new_node_number) + ".1", "dis->:"+str(round(self.mean_in[1], 3)),
            "dis<-:"+str(round(self.mean_out[1], 3)), "size:"+str(round(np.exp(self.curnodelength[1]), 1)),
                                self.short_labels_tree[1]])
        self.edgelength.append(["N" + str(self.new_node_number) + ".2", "dis->:"+str(round(self.mean_in[2], 3)),
            "dis<-:"+str(round(self.mean_out[2], 3)), "size:"+str(round(np.exp(self.curnodelength[2]), 1)),
                                self.short_labels_tree[2]])

        self.newick_tree.append("(" + self.newick_tree[0] + ":" + str(round(max([0, self.mean_in[0]
            + self.mean_out[0]]), 3)) + "," + self.newick_tree[1] + ":" + str(round(max([0, self.mean_in[1]
            + self.mean_out[1]]), 3)) + "," + self.newick_tree[2] + ":" + str(round(max([0, self.mean_in[2]
            + self.mean_out[2]]), 3)) + ")" + ";")

        self.short_labels_tree.append("(" + self.short_labels_tree[0] + ":" + str(round(max([0, self.mean_in[0]
            + self.mean_out[0]]), 3)) + "," + self.short_labels_tree[1] + ":" + str(
            round(max([0, self.mean_in[1] + self.mean_out[1]]), 3)) + "," + self.short_labels_tree[2] + ":" + str(
            round(max([0, self.mean_in[2] + self.mean_out[2]]), 3)) + ")" + ";")
        print()
        print("Newick of constructed tree:",self.short_labels_tree[3])
        print()
        print("Construction steps including: new node(mother) label, distance to sister, distance to mother, estimated size of mother,subtree of the sister:")
        print()
        for m in range(len(self.edgelength)):
            print(self.edgelength[m])
        tr = Tree(self.newick_tree[3], format=1)
        print()
        print("Constructed tree figure:")
        print()
        print(tr)

        unimog_treedata = "unimog" + str(self.atgcnum) + "_treedata.txt"
        f = open(unimog_treedata, "r")
        newicktr = f.read()
        unimog_tnj = Tree(newicktr, format=3)

        sgc_treedata = "sgc" + str(self.atgcnum) + "_treedata.txt"
        f = open(sgc_treedata, "r")
        newicktr = f.read()
        sgc_tnj = Tree(newicktr, format=3)
        trpng = "tddtr" + str(self.atgcnum) + ".png"
        tr.render(trpng, w=183, units="mm")

        treedata = "atgc" + str(self.atgcnum) + "_treedata.txt"
        f = open(treedata, "r")

        newicktr_atgc = f.read()
        tatgc = Tree(newicktr_atgc, format=3)

        maxdis = 2 * self.nleaves - 6
        d01 = tr.robinson_foulds(tatgc, unrooted_trees=True)[0] / maxdis
        d02 = tr.robinson_foulds(unimog_tnj, unrooted_trees=True)[0] / maxdis
        d03 = tr.robinson_foulds(sgc_tnj, unrooted_trees=True)[0] / maxdis
        d12 = tatgc.robinson_foulds(unimog_tnj, unrooted_trees=True)[0] / maxdis
        d13 = tatgc.robinson_foulds(sgc_tnj, unrooted_trees=True)[0] / maxdis
        d23 = unimog_tnj.robinson_foulds(sgc_tnj, unrooted_trees=True)[0] / maxdis

        atgcpng = "atgc" + str(self.atgcnum) + ".png"
        tatgc.render(atgcpng, w=183, units="mm")
        print()
        print("Tree figure downloaded from ATGCs site:")
        print(tatgc)
        print()
        print("Normalized Robinson-Foulds distances:")
        print()
        print("ATGC" +str(self.atgcnum), "TDDR-ATGC:",round(d01, 3), "TDDR-DJC:",round(d02, 3),
              "TDDR-SGC:",round(d03, 3), "ATGC-DCJ:",round(d12, 3), "ATGC-SGC:",round(d13, 3),
              "DCJ-SGC:",round(d23, 3), "max distance for normalization:", maxdis)

    def tddtr_prog(self,atgcnum):
        time0 = time.time()
        self.atgcnum = atgcnum
        self.readatgc()
        self.preparing_loop()
        self.intersect()
        self.initial_direct_dist()


        while self.nmax > 3:
            self.choosing_cherry()
            self.updating_data()
            self.nmax += -1

        self.program_end()
        time1 = time.time()
        return time1 - time0
print("Total time:",compare_tddr())

