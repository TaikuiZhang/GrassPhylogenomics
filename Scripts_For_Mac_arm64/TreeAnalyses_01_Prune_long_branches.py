#!/usr/bin/env python3
import os,sys,re
import numpy as np
import collections
from Bio import Phylo
from Bio.Phylo.Consensus import *
from Bio.Phylo.Consensus import _BitString
from GRASSMODULE_arm64.Rooting import BranchAvg,Prune

def progress_bar(finished_number,tasks_numbers):
    percentage=round(int(finished_number)/int(tasks_numbers)*100)
    print("\rDone: {}%: ".format(percentage),"▋" * (percentage // 2), end="")
    sys.stdout.flush()

def handle():
    return None

def deleteDuplicatedElement(listA):
    return sorted(set(listA), key = listA.index)


if __name__ == "__main__":
    if len(sys.argv)<2:
        print ("Usage: python3 TreeAnalyses_01_Prune_long_branches.py RawTree.list")
        exit()
inList = sys.argv[1]
cmd="mkdir LB_removed"
os.system(cmd)
TreeNumbers=0
with open(inList,"r") as infile:
    for line in infile:
        TreeNumbers=TreeNumbers+1
TreeNumber=0
print("Processing ",TreeNumbers," gene trees")
with open(inList,"r") as infile:
    for line in infile:
        line=line.rstrip()
        spls = line.strip().split("\t")
        TreNum = spls[0]
        inTree = spls[1]
        TreeNumber=TreeNumber+1
        progress_bar(TreeNumber, TreeNumbers)
        tree = Phylo.read(inTree, "newick")
        AVGL=BranchAvg(tree)
        j=0
        MAX=len(tree.get_terminals())
        while(j<int(MAX)):
            j=j+1
            tree=Prune(tree)
            if len(tree.get_terminals())>6:
                AVGLj=BranchAvg(tree)
                if AVGLj==AVGL:
                    break
                else:
                    AVGL=AVGLj
            else:
                break
        if len(tree.get_terminals())>6:
            for Node in tree.get_nonterminals():
                Node.name=""
            outTree=str(inTree)+".LB.removed.tree1"
            Phylo.write(tree, outTree, "newick")
            outTrees=str(inTree)+".LB.removed.tree"
            fp = open(outTrees, "w")
            with open(outTree,"r") as infile:
                for line in infile:
                    tline=line.rstrip()
                    tail=':0.00000;'
                    tline=tline.replace(tail,";")
                    print(tline,file=fp)
            fp.close()
            cmd="rm -f "+str(outTree)
            os.system(cmd)
            cmd='mv '+str(outTrees)+' LB_removed/'
            os.system(cmd)
        
             