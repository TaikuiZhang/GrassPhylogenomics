#!/usr/bin/env python3
import os,sys,re
import numpy as np
import collections
from Bio import Phylo

from GRASSMODULE_arm64.Rooting import reroot,Tip_Prune,BLtree,SearchSisterLC,RootBranch,MAB,miniAncestralDup,SearchSister,SearchClade,SplitTips,RBL_MAD,ROOTatEvoTree
from GRASSMODULE_arm64.Estimation import NodeEstimation,TEXT_NR

def deleteDuplicatedElement(listA):
    return sorted(set(listA), key = listA.index)

def BranchAvg(tree):
    BranchLen=[]
    for clade in tree.get_terminals():
        if (clade != tree.root):
            CD=clade.branch_length
            CD=float(CD)
            BranchLen.append(CD)
    for clade in tree.get_nonterminals():
        CD=clade.branch_length
        if isinstance(CD,float)==True:
            CD=float(CD)
            BranchLen.append(CD)
    avg=np.mean(BranchLen)
    return avg

def progress_bar(finished_number,tasks_numbers):
    percentage=round(int(finished_number)/int(tasks_numbers)*100)
    print("\rDone "+format('['+str(finished_number)+'/'+str(tasks_numbers)+']')+" {}%: ".format(percentage),"▋" * (percentage // 2), end="")
    sys.stdout.flush()
    if int(finished_number)==int(tasks_numbers):
        print("\n")

if __name__ == "__main__":
    print("Root gene trees using EvoMADs")
    print("The initial run to root gene trees used the the minimal ancestor deviation(MAD)")
    print("approach (Tria et al., 2017 Nature Ecology & Evolution; see their mad.py script).")
    print("Additional analyses were performed to keep the minimal ancestral duplication(MAD).")
    print ("##make by Taikui Zhang, PhD##\n")
    if len(sys.argv)<3:
        print ("Usage: python3 TreeAnalyses_02_root_using_EvoMADs.py HOG.pineapple ref.topology.tree Grass_tips")
        exit()
inList = sys.argv[1]
ref = sys.argv[2]
dateAnno = sys.argv[3]
reftree = Phylo.read(ref, "newick")

cmd="mkdir RerootTree"
os.system(cmd)


logfile="run_rerooting.log"
fp = open(logfile, "w")

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
        HOG=spls[0]
        osp=spls[1]
        additionalBS=0
        TreeNumber=TreeNumber+1
        progress_bar(TreeNumber, TreeNumbers)
        inTree='LB_removed/'+HOG+".iqtree.treefile.LB.removed.tree"
        inTreeMAD='LB_removed/'+HOG+".iqtree.treefile.LB.removed.tree.MAD_rooted.tree"
        reTree = 'RerootTree/'+HOG+'.root.tree'
        if os.path.exists(inTree)==True:
            Atree = Phylo.read(inTree, "newick")
            if len(Atree.get_terminals())>5:
                cmd='python3 mad.py '+str(inTree)+' -n'
                os.system(cmd)
                MADtree = Phylo.read(inTreeMAD, "newick")
                tree = Phylo.read(inTree, "newick")
                CR=tree.root
                if len(CR.clades)>2:
                    new_A=CR.clades[-1]
                    ot=float(new_A.branch_length)/2
                    tree.root_with_outgroup(new_A,outgroup_branch_length=ot)
                tp=[]
                avgl=BranchAvg(MADtree)
                #print(tree)
                for Node in MADtree.get_nonterminals():
                    MA=Node.get_terminals()
                    if len(MA)<3:
                        SubAL=Node.branch_length
                        if SubAL > 10*avgl:
                            for tid in MA:
                                tp.append(tid.name)
                if len(tp)>0:
                    #print(tp)
                    tree=Tip_Prune(tree,tp)
                    MADtree=Tip_Prune(MADtree,tp)
                tp=[]
                #print(MADtree)
                CB=RootBranch(MADtree)
                node=MAB(tree,CB)
                Ntree=reroot(tree,node)
                TREEDB=[]
                againosp=0
                tree1=tree
                for gene in Ntree.get_terminals():
                    if gene.name==osp:
                        againosp=1
                        additionalBS=SearchSisterLC(gene,tree)
                        if additionalBS==None:
                            additionalBS=100
                for node in tree1.get_nonterminals():
                    if node.confidence==None:
                        node.confidence=100
                tree1.root.confidence=100
                tree1.ladderize(reverse=True)
                Phylo.write(tree1, reTree, "newick")
                
                DUP1=NodeEstimation(tree1,reftree,dateAnno)
                tree2=tree1
                SN=[]
                if againosp==1:
                    SN=SearchSister(tree1,osp)
                    tree2=reroot(tree1,osp)
                else:
                    additionalBS=100
                for node in tree2.get_nonterminals():
                    if node.confidence==None:
                        node.confidence=100
                tree2.root.confidence=additionalBS
                tree2.ladderize(reverse=True)
                DUP2=NodeEstimation(tree2,reftree,dateAnno)

                GDUP1=[ dup for dup in DUP1 if (dup.split(".")[0]=='Node2' or dup.split(".")[0]=='Node3' or dup.split(".")[0]=='Node4')]
                GDUP2=[ dup for dup in DUP2 if (dup.split(".")[0]=='Node2' or dup.split(".")[0]=='Node3' or dup.split(".")[0]=='Node4')]
                
                PDUP1=[ dup for dup in DUP1 if (dup.split(".")[0]=='Node1')]
                PDUP2=[ dup for dup in DUP2 if (dup.split(".")[0]=='Node1')]

                mDup12=miniAncestralDup(len(GDUP1),len(GDUP2),len(PDUP1),len(PDUP2))
                if againosp==0:
                    TPtree=Phylo.read(reTree, "newick")
                    TPD=[]
                    for node in tree2.get_nonterminals():
                        if node!=tree2.root:
                            GIDs=SearchClade(tree2,node)
                            tree4=reroot(TPtree,GIDs)
                            tree4.root.confidence=100
                            for nd in tree4.get_nonterminals():
                                if nd.confidence==None:
                                    nd.confidence=100
                            DUP4=NodeEstimation(tree4,reftree,dateAnno)
                            GDUP4=[ dup for dup in DUP4 if (dup.split(".")[0]=='Node2' or dup.split(".")[0]=='Node3' or dup.split(".")[0]=='Node4')]
                            if len(GDUP4)>0:
                                GDUP41=[ dup for dup in DUP4 if (dup.split(".")[0]=='Node2')]
                                GDUP42=[ dup for dup in DUP4 if (dup.split(".")[0]=='Node3')]
                                d41=0
                                e41=0
                                l41=len(GDUP41)
                                if len(GDUP41)>0:
                                    d41=GDUP41[0].split(".")[1]
                                    e41=GDUP41[0].split(".")[2]
                                    if l41==1:
                                        l41=-1
                                d42=0
                                e42=0
                                l42=len(GDUP42)
                                if len(GDUP42)>0:
                                    d42=GDUP42[0].split(".")[1]
                                    e42=GDUP42[0].split(".")[2]
                                    if l42==1:
                                        l42=-1
                                key=str(l41)+"\t"+str(d41)+"\t"+str(e41)+"\t"+str(l42)+"\t"+str(d42)+"\t"+str(e42)+"\t"+str(len(GDUP4))+"\t"+str(GIDs)
                                TPD.append(key)
                    
                    if len(TPD)>0:
                        TPD=sorted(TPD, key=lambda x:(float(x.split("\t")[1]),float(x.split("\t")[2])),reverse=True)
                        TPD=sorted(TPD, key=lambda x:(float(x.split("\t")[0])),reverse=False)
                        KEYt=TPD[0].split("\t")[0]
                        if int(KEYt)==0:
                            TPD=sorted(TPD, key=lambda x:(float(x.split("\t")[5]),float(x.split("\t")[4])),reverse=True)
                        #print(TPD)
                        KEY=TPD[0].split("\t")[-1]
                        KEYS=SplitTips(KEY)
                        #print(KEYS)
                        tree4=reroot(TPtree,KEYS)
                        tree4.root.confidence=100
                        DUP4=NodeEstimation(tree4,reftree,dateAnno)
                        GDUP4=[ dup for dup in DUP4 if (dup.split(".")[0]=='Node2' or dup.split(".")[0]=='Node3' or dup.split(".")[0]=='Node4')]
                        PDUP4=[ dup for dup in DUP4 if (dup.split(".")[0]=='Node1')]
                        mDup14=miniAncestralDup(len(GDUP1),len(GDUP4),len(PDUP1),len(PDUP4))
                        mDup24=miniAncestralDup(len(GDUP2),len(GDUP4),len(PDUP2),len(PDUP4))
                        if int(mDup12)==1:
                            if int(mDup14)==2:
                                tree4=RBL_MAD(tree4)
                                Phylo.write(tree4, reTree, "newick")
                            else:
                                tree1=RBL_MAD(TPtree)
                                Phylo.write(tree1, reTree, "newick")
                        else:
                            if int(mDup14)==2:
                                if int(mDup24)==2:
                                    tree4=RBL_MAD(tree4)
                                    Phylo.write(tree4, reTree, "newick")
                                else:
                                    tree2=RBL_MAD(TPtree)
                                    Phylo.write(tree2, reTree, "newick")
                            else:
                                tree2=RBL_MAD(TPtree)
                                Phylo.write(tree2, reTree, "newick")
                        logf=HOG+"\t2"
                        print(logf,file=fp)
                    else:
                        #TPtree=Phylo.read(reTree, "newick")
                        Outsp=ROOTatEvoTree(TPtree,reftree,dateAnno)
                        if len(Outsp)>0:
                            #print(Outsp)
                            tree1=reroot(TPtree,Outsp)
                            #print(tree1)
                            tree1=RBL_MAD(tree1)
                            Phylo.write(tree1, reTree, "newick")
                        else:
                            tree1=RBL_MAD(TPtree)
                            Phylo.write(tree1, reTree, "newick")
                        logf=HOG+"\t3"
                        print(logf,file=fp) 
                else:
                    TPtree=Phylo.read(reTree, "newick")
                    TPD=[]
                    for node in tree2.get_nonterminals():
                        if node!=tree2.root:
                            GIDs=SearchClade(tree2,node)
                            tree4=reroot(TPtree,GIDs)
                            tree4.root.confidence=100
                            for nd in tree4.get_nonterminals():
                                if nd.confidence==None:
                                    nd.confidence=100
                            DUP4=NodeEstimation(tree4,reftree,dateAnno)
                            GDUP4=[ dup for dup in DUP4 if (dup.split(".")[0]=='Node1' or dup.split(".")[0]=='Node2' or dup.split(".")[0]=='Node3' or dup.split(".")[0]=='Node4')]
                            GDUP41=[ dup for dup in DUP4 if (dup.split(".")[0]=='Node1')]
                            GDUP42=[ dup for dup in DUP4 if (dup.split(".")[0]=='Node2')]
                            GDUP43=[ dup for dup in DUP4 if (dup.split(".")[0]=='Node3')]
                            d41=0
                            e41=0
                            l41=len(GDUP41)
                            b41=0
                            if len(GDUP41)>0:
                                d41=GDUP41[0].split(".")[1]
                                e41=GDUP41[0].split(".")[2]
                                b41=GDUP41[0].split(".")[3]
                                if l41==1:
                                    l41=-1
                            l42=len(GDUP42)     
                            d42=0
                            e42=0
                            b42=0
                            if len(GDUP42)>0:
                                d42=GDUP42[0].split(".")[1]
                                e42=GDUP42[0].split(".")[2]
                                b42=GDUP42[0].split(".")[3]
                                if l42==1:
                                    l42=-1
                            l43=len(GDUP43)
                            d43=0
                            e43=0
                            b43=0
                            if len(GDUP43)>0:
                                d43=GDUP43[0].split(".")[1]
                                e43=GDUP43[0].split(".")[2]
                                b43=GDUP43[0].split(".")[3]
                                if l43==1:
                                    l43=-1 
                            if len(GDUP4)>0:
                                key=str(l41)+"\t"+str(e41)+"\t"+str(d41)+"\t"+str(b41)+"\t"+str(l42)+"\t"+str(e42)+"\t"+str(d42)+"\t"+str(b42)
                                key=key+"\t"+str(l43)+"\t"+str(e43)+"\t"+str(d43)+"\t"+str(b43)+"\t"+str(GIDs)
                                TPD.append(key)
                    if len(TPD)>0:
                        TPD=sorted(TPD, key=lambda x:(float(x.split("\t")[1]),float(x.split("\t")[2]),float(x.split("\t")[3])),reverse=True)
                        TPD=sorted(TPD, key=lambda x:(float(x.split("\t")[0])),reverse=False)
                        TPD=sorted(TPD, key=lambda x:(float(x.split("\t")[5]),float(x.split("\t")[6])),reverse=True)
                        TPD=sorted(TPD, key=lambda x:(float(x.split("\t")[4])),reverse=False)
                        KEYt=TPD[0].split("\t")[4]
                        if int(KEYt)==0:
                            TPD=sorted(TPD, key=lambda x:(float(x.split("\t")[9]),float(x.split("\t")[10])),reverse=True)
                        KEY=TPD[0].split("\t")[-1]
                        KEYS=SplitTips(KEY)
                        tree4=reroot(TPtree,KEYS)
                        tree4.root.confidence=100
                        DUP4=NodeEstimation(tree4,reftree,dateAnno)
                        GDUP4=[ dup for dup in DUP4 if (dup.split(".")[0]=='Node2' or dup.split(".")[0]=='Node3' or dup.split(".")[0]=='Node4')]
                        PDUP4=[ dup for dup in DUP4 if (dup.split(".")[0]=='Node1')]
                        mDup14=miniAncestralDup(len(GDUP1),len(GDUP4),len(PDUP1),len(PDUP4))
                        mDup24=miniAncestralDup(len(GDUP2),len(GDUP4),len(PDUP2),len(PDUP4))
                        if int(mDup12)==1:
                            if int(mDup14)==2:
                                tree4=RBL_MAD(tree4)
                                Phylo.write(tree4, reTree, "newick")    
                            else:
                                tree1=RBL_MAD(TPtree)
                                Phylo.write(tree1, reTree, "newick")
                        else:
                            if int(mDup14)==2:
                                if int(mDup24)==2:
                                    tree4=RBL_MAD(tree4)
                                    Phylo.write(tree4, reTree, "newick")
                                else:
                                    tree2=RBL_MAD(tree2)
                                    Phylo.write(tree2, reTree, "newick")
                            else:
                                if int(mDup24)==1:
                                    tree2=RBL_MAD(tree2)
                                    Phylo.write(tree2, reTree, "newick")
                                else:
                                    tree4=RBL_MAD(tree4)
                                    Phylo.write(tree4, reTree, "newick")
                        logf=HOG+"\t4"
                        print(logf,file=fp)
                    else:
                        logf=HOG+"\t5"
                        print(logf,file=fp)
                        if int(mDup12)>1:
                            tree2=RBL_MAD(tree2)
                            Phylo.write(tree2, reTree, "newick")
            else:
                logf=HOG+"\tNA\tNA\t4"
                print(logf,file=fp)
fp.close()
cmd='find RerootTree/ -name "*.tree" | cat -n > Rho.tree_list'
os.system(cmd)