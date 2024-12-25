# TDDR
You can find here two Phyton files:
1. trees simulations - sim_tr.py
   
  This program conducts simulations:
  
  It creates random trees and reconstructs them using TDDR, SGC+BioPhyton NJ, and DCJ + BioPhyton NJ.
  
  To run, Java 8 must be installed. You must also download UniMoG (BiBiServ2 - DCJ), which should be in the same folder as the code.

  command line: py sim_tr.py

  sim_tr.py arguments : 
  
  trials = 10       				 - Number of repetitions
  
  edge_ml = 0.05               - Mean edge length. Each directional edge length is from 
                                    		  exponential distribution with that mean.                                      

  leaves_n = 25                 - Number of leaves (taxa) in the random tree
  
  genes_n = 2000                - Number of genes in at the root genome

  command line: py sim_tr.py 10 0.05 25 2000
  
  Output: The results are printed on the screen.

2. real_data.py
   
   Test TDDR on ATDCs: ATDC005, ATDC007, ATDC008, ATDC009, ATDC032

   Compare the RF-resulting trees with trees constructed with the ATGC tree,  SGC+BioPhyton NJ tree, and DCJ + BioPhyton NJ tree. 

  Argument:
  
  atgcnum = "5"               The ATGC  number

  Input files: (all files included the respiratory and should be placed with the code).

  Files ATGCs genomes downloaded from ATGCs site and reduced to relevant columns:
  
  Columns: #cog, genomeassembly, protein.

  ATGC5reduced.csv
  
  ATGC7reduced.csv
  
  ATGC8reduced.csv
  
  ATGC9reduced.csv
  
  ATGC32reduced.csv

  Trees in Newick format: 

  Downloaded from ATGCs site -

  atgc5_treedata.txt
  
  atgc7_treedata.txt
  
  atgc8_treedata.txt
  
  atgc9_treedata.txt
  
  atgc32_treedata.txt

  Constructed by SGC+BioPhyton NJ -

  sgc5_treedata.txt
  
  sgc7_treedata.txt
  
  sgc8_treedata.txt
  
  sgc9_treedata.txt
  
  sgc32_treedata.txt

  Constructed by DCJ+BioPhyton NJ -

  dcj5_treedata.txt
  
  dcj7_treedata.txt
  
  dcj8_treedata.txt
  
  dcj9_treedata.txt
  
  dcj32_treedata.txt

  Output:

  Normalized RF distances are printed on the screen. 
  
  TDDR tree figure:
  
  one of:
  
  tddtr5.png, tddtr7.png, tddtr8.png, tddtr9.png, tddtr32.png 
  
  (according to the value of atgcnum).






