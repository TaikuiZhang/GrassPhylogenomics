#!/usr/bin/env python3
import os,sys,re
import numpy as np
import collections
from Bio import Phylo

from GRASSMODULE_arm64.Estimation import DupEstimation,CladeEstimationA,CladeEstimationB,CladeEstimationC,SplitTrees,Class_NR_SUBFAM

def handle():
    return None

def deleteDuplicatedElement(listA):
    return sorted(set(listA), key = listA.index)

def progress_bar(finished_number,tasks_numbers):
    percentage=round(int(finished_number)/int(tasks_numbers)*100)
    print("\rDone: {}%: ".format(percentage),"▋" * (percentage // 2), end="")
    sys.stdout.flush()

def retention_num(tree,reftree,dateAnno):
    DUP=DupEstimation(tree,reftree,dateAnno)
    KDUP=[ dup for dup in DUP if (dup.split(".")[0]=='Node2' or dup.split(".")[0]=='Node3' or dup.split(".")[0]=='Node4')]
    if len(KDUP)==0:
        Rentention=1
        Chloridoideae=CladeEstimationA('Node6',tree,reftree,dateAnno)
        Panicoideae=CladeEstimationA('Node8',tree,reftree,dateAnno)
        Pooideae=CladeEstimationA('Node14',tree,reftree,dateAnno)
        Bambusoideae=CladeEstimationA('Node18',tree,reftree,dateAnno)
        Oryzoideae=CladeEstimationA('Node21',tree,reftree,dateAnno)
        Pharoideae=CladeEstimationB('Pharus_latifolius',tree,dateAnno)
        Anomochlooideae=CladeEstimationB('Streptochaeta_angustifolia',tree,dateAnno)
        if ('(' not in str(Panicoideae.split("\t")[0])) and ('(' not in str(Bambusoideae.split("\t")[0])) and ('(' not in str(Chloridoideae.split("\t")[0])) and ('(' not in str(Pharoideae.split("\t")[0])) and ('(' not in str(Pooideae.split("\t")[0])) and ('(' not in str(Oryzoideae.split("\t")[0])) and ('(' not in str(Anomochlooideae.split("\t")[0])):
            Rentention='1C'
        else:
            if ('SD' in str(Panicoideae.split("\t")[0])) or ('SD' in str(Bambusoideae.split("\t")[0])) or ('SD' in str(Chloridoideae.split("\t")[0])) or ('SD' in str(Pharoideae.split("\t")[0])) or ('SD' in str(Pooideae.split("\t")[0])) or ('SD' in str(Oryzoideae.split("\t")[0])) or ('SD' in str(Anomochlooideae.split("\t")[0])):
                Rentention='1B'
            else:
                Rentention='1A'

        out=str(Rentention)+"\t"+str(Chloridoideae.split("\t")[0])
        out=out+"\t"+str(Panicoideae.split("\t")[0])+"\t"+str(Pooideae.split("\t")[0])
        out=out+"\t"+str(Bambusoideae.split("\t")[0])+"\t"+str(Oryzoideae.split("\t")[0])
        out=out+"\t"+str(Pharoideae.split("\t")[0])+"\t"+str(Anomochlooideae.split("\t")[0])
        out=out+"\t"+str(Chloridoideae.split("\t")[1])
        out=out+"\t"+str(Panicoideae.split("\t")[1])+"\t"+str(Pooideae.split("\t")[1])
        out=out+"\t"+str(Bambusoideae.split("\t")[1])+"\t"+str(Oryzoideae.split("\t")[1])
        out=out+"\t"+str(Pharoideae.split("\t")[1])+"\t"+str(Anomochlooideae.split("\t")[1])
        return out
    if len(KDUP)>0:
        Rentention=2
        Chloridoideae=CladeEstimationC('Node6',KDUP,tree,reftree,dateAnno)
        Panicoideae=CladeEstimationC('Node8',KDUP,tree,reftree,dateAnno)
        Pooideae=CladeEstimationC('Node14',KDUP,tree,reftree,dateAnno)
        Bambusoideae=CladeEstimationC('Node18',KDUP,tree,reftree,dateAnno)
        Oryzoideae=CladeEstimationC('Node21',KDUP,tree,reftree,dateAnno)
        Pharoideae=CladeEstimationC('Pharus_latifolius',KDUP,tree,reftree,dateAnno)
        Anomochlooideae=CladeEstimationC('Streptochaeta_angustifolia',KDUP,tree,reftree,dateAnno)
        if ('SubRL' in str(Panicoideae.split("\t")[0])) or ('SubRL' in str(Bambusoideae.split("\t")[0])) or ('SubRL' in str(Chloridoideae.split("\t")[0])) or ('SubRL' in str(Pharoideae.split("\t")[0])) or ('SubRL' in str(Pooideae.split("\t")[0])) or ('SubRL' in str(Oryzoideae.split("\t")[0])) or ('SubRL' in str(Anomochlooideae.split("\t")[0])):
            Rentention='2B'
        else:
            ChNum=0
            if '(' in str(Chloridoideae.split("\t")[0]):
                ChNum=int(Chloridoideae.split("(")[0])
            else:
                ChNum=int(Chloridoideae.split("\t")[0])
            
            PaNum=0
            if '(' in str(Panicoideae.split("\t")[0]):
                PaNum=int(Panicoideae.split("(")[0])
            else:
                PaNum=int(Panicoideae.split("\t")[0])
            
            PoNum=0
            if '(' in str(Pooideae.split("\t")[0]):
                PoNum=int(Pooideae.split("(")[0])
            else:
                PoNum=int(Pooideae.split("\t")[0])
            
            BaNum=0
            if '(' in str(Bambusoideae.split("\t")[0]):
                BaNum=int(Bambusoideae.split("(")[0])
            else:
                BaNum=int(Bambusoideae.split("\t")[0])
            
            OrNum=0
            if '(' in str(Oryzoideae.split("\t")[0]):
                OrNum=int(Oryzoideae.split("(")[0])
            else:
                OrNum=int(Oryzoideae.split("\t")[0])
            
            PhNum=0
            if '(' in str(Pharoideae.split("\t")[0]):
                PhNum=int(Pharoideae.split("(")[0])
            else:
                PhNum=int(Pharoideae.split("\t")[0])
            
            AnNum=0
            if '(' in str(Anomochlooideae.split("\t")[0]):
                AnNum=int(Anomochlooideae.split("(")[0])
            else:
                AnNum=int(Anomochlooideae.split("\t")[0])
            
            SUM=int(ChNum)+int(PaNum)+int(PoNum)+int(BaNum)+int(OrNum)+int(PhNum)+int(AnNum)
            
            if SUM==14:
                Rentention='2A'
            elif SUM==13:
                Rentention='2C'
            else:
                if (int(ChNum)!=1) and (int(PaNum)!=1) and (int(PoNum)!=1) and (int(BaNum)!=1) and (int(OrNum)!=1) and (int(PhNum)!=1) and (int(AnNum)!=1):
                    Rentention='2D'
                else:
                    Rentention='2E'
        out=str(Rentention)+"\t"+str(Chloridoideae.split("\t")[0])
        out=out+"\t"+str(Panicoideae.split("\t")[0])+"\t"+str(Pooideae.split("\t")[0])
        out=out+"\t"+str(Bambusoideae.split("\t")[0])+"\t"+str(Oryzoideae.split("\t")[0])
        out=out+"\t"+str(Pharoideae.split("\t")[0])+"\t"+str(Anomochlooideae.split("\t")[0])
        out=out+"\t"+str(Chloridoideae.split("\t")[1])
        out=out+"\t"+str(Panicoideae.split("\t")[1])+"\t"+str(Pooideae.split("\t")[1])
        out=out+"\t"+str(Bambusoideae.split("\t")[1])+"\t"+str(Oryzoideae.split("\t")[1])
        out=out+"\t"+str(Pharoideae.split("\t")[1])+"\t"+str(Anomochlooideae.split("\t")[1])
        return out

if __name__ == "__main__":
    print ("##Estimating rentention and loss events of the rho-derived duplicates in different grass subfamilies##")
    print ("##make by Taikui Zhang, PhD##")
    if len(sys.argv)<4:
        print ("Usage: python3 TreeAnalyses_03_estimate_retention.py Rho.tree_list ref.topology.tree Grass_tips Rho.tree_list.dup.txt")
        exit()
inList = sys.argv[1]
ref = sys.argv[2]
dateAnno = sys.argv[3]
OUT = sys.argv[4]
reftree = Phylo.read(ref, "newick")
fp = open(OUT, "w")

out="TreNum\tTree\tSubTrees\tRentention\tChloridoideae#\tPanicoideae#"
out=out+"\tPooideae#\tBambusoideae#\tOryzoideae#\tPharoideae#\tAnomochlooideae#"
out=out+"\tChloridoideae\tPanicoideae"
out=out+"\tPooideae\tBambusoideae\tOryzoideae\tPharoideae\tAnomochlooideae"
print(out,file=fp)

logfile="run_estimating.log"
fpl = open(logfile, "w")

TreeNumbers=0
with open(inList,"r") as infile:
    for line in infile:
        TreeNumbers=TreeNumbers+1

GD_ALL=[]
GDN=0
TreeNumber=0
H=0
print("Processing ",TreeNumbers," gene trees")
with open(inList,"r") as infile:
    for line in infile:
        line=line.rstrip()
        spls = line.strip().split("\t")
        TreNum = spls[0]
        inTree = spls[1]
        TreeNumber=TreeNumber+1
        progress_bar(TreeNumber, TreeNumbers)
        inTreeSize=os.path.getsize(inTree)
        if inTreeSize==0:
            print (TreNum,'size=',inTreeSize)
        else:
            tree = Phylo.read(inTree, "newick")
            jp1=len(tree.get_nonterminals())
            jp2=len(tree.get_terminals())
            if (jp1+1)!=jp2:
                print (TreNum,"polyphyly:",jp1,jp2)
            else:
                Ni=0
                for Node in tree.get_nonterminals():
                    Ni=Ni+1
                    Node.name='gn'+str(Ni)
                tree.root.confidence=100
                out=str(TreNum)+"\t"+str(inTree)
                DUP=DupEstimation(tree,reftree,dateAnno)
                ADUP=[ dup for dup in DUP if (dup.split(".")[0]=='Node1')]
                GDUP=[ dup for dup in DUP if (dup.split(".")[0]=='Node2')]
                PDUP=[ dup for dup in DUP if (dup.split(".")[0]=='Node3')]
                CDUP=[ dup for dup in DUP if (dup.split(".")[0]=='Node4')]
                logf=out+"\t"+str(len(ADUP))+"\t"+str(len(GDUP))+"\t"+str(len(PDUP))+"\t"+str(len(CDUP))
                print(logf,file=fpl)
                Subtrees=SplitTrees(tree,ADUP,GDUP,PDUP,CDUP)
                if len(Subtrees)>1:
                    snum=0  
                    for stree in Subtrees:
                        NSUBF=Class_NR_SUBFAM(stree.get_terminals(),dateAnno)
                        if len(stree.get_terminals())>5 and len(NSUBF)>3:
                            snum=snum+1
                            KDUP=retention_num(stree,reftree,dateAnno)
                            outA=out+'\t'+str(snum)+'\t'+KDUP
                            print(outA,file=fp)
                else:
                    for stree in Subtrees:
                        NSUBF=Class_NR_SUBFAM(stree.get_terminals(),dateAnno)
                        if len(stree.get_terminals())>5 and len(NSUBF)>3:
                            KDUP=retention_num(stree,reftree,dateAnno)
                            outA=out+'\t-\t'+KDUP
                            print(outA,file=fp)
       
fp.close()
fpl.close()