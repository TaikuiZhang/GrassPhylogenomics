#!/usr/bin/env python3
import os,sys,re
import numpy as np
from itertools import islice
import collections

def handle():
    return None

def deleteDuplicatedElement(listA):
    return sorted(set(listA), key = listA.index)

def progress_bar(finished_number,tasks_numbers):
    percentage=round(int(finished_number)/int(tasks_numbers)*100)
    print("\rDone: {}%: ".format(percentage),"▋" * (percentage // 2), end="")
    sys.stdout.flush()

def LIST_NR(CLAS):
    Text="-"
    if len(CLAS)>0:
        Text='1: '+str(CLAS[0])
        for i in range(1,len(CLAS)):
            Text=Text+"; "+str(i+1)+": "+str(CLAS[i])
    return Text

def TEXT_NR(CLAS):
    Text="-"
    if len(CLAS)>0:
        Text=str(CLAS[0])
        for i in range(1,len(CLAS)):
            Text=Text+","+str(CLAS[i])
    return Text

def ANCHECK(CLUST):
    Y=0
    for i in range(0,len(CLUST)):
        for j in range(i+1, len(CLUST)):
            A=CLUST[i]
            B=CLUST[j]
            AB=A+"\t"+B
            if ANCIENT.get(AB)!=None:
                Y=1
    return Y

if __name__ == "__main__":
    print ("##Estimating rentention and loss events of the rho-derived duplicates in different grass subfamilies##")
    print ("##make by Taikui Zhang, PhD##")
inList = 'Rho.tree_list.dup.txt'
dateAnno = 'Verification/ge.id_list'
OUT = 'Rho.tree_list.dup.anno.txt'
fp = open(OUT, "w")
GH={}
with open(dateAnno,"r") as infile:
    for line in infile:
        line=line.rstrip()
        spls = line.strip().split("\t")
        A = spls[0]
        B = spls[1]
        GH[A]=B

GeAN="Verification/ancient_pairs_ABAB_earlier_than_rho"
ANCIENT={}
with open(GeAN,"r") as infile:
    for line in infile:
        line=line.rstrip()
        spls = line.strip().split("\t")
        A = spls[2]
        B = spls[3]
        AB=A+"\t"+B
        BA=B+"\t"+A
        ANCIENT[AB]=1
        ANCIENT[BA]=1

out="OG#\tSubTrees\tRentention\tDetail\tChloridoideae\tPanicoideae"
out=out+"\tPooideae\tBambusoideae\tOryzoideae\tPharoideae\tAnomochlooideae"
out=out+"\tChloridoideae\tPanicoideae"
out=out+"\tPooideae\tBambusoideae\tOryzoideae\tPharoideae\tAnomochlooideae"
print(out,file=fp)

TreeNumbers=0
with open(inList,"r") as infile:
    for line in islice(infile, 1, None):
        TreeNumbers=TreeNumbers+1

GD_ALL=[]
GDN=0
TreeNumber=0
H=0
I_DB=[]
II_DB=[]
III_DB=[]
IV_DB=[]
rls=0
print("Processing ",TreeNumbers," gene trees")
with open(inList,"r") as infile:
    for line in islice(infile, 1, None):
        line=line.rstrip()
        spls = line.strip().split("\t")
        TreeNumber=TreeNumber+1
        progress_bar(TreeNumber, TreeNumbers)
        HOG=spls[1].strip().split(".")[0]
        HOG=HOG.strip().split("/")[-1]
        SubTrees = spls[2]
        R="-"
        if spls[3]=="1A":
            R='I'
            I_DB.append(HOG)
        if spls[3]=="1B":
            R='I'
            I_DB.append(HOG)
        if spls[3]=="1C":
            R='II'
            II_DB.append(HOG)
        if spls[3]=="2A":
            R='III'
            III_DB.append(HOG)
        if spls[3]=="2B" or spls[3]=="2C" or spls[3]=="2D" or spls[3]=="2E":
            R='IV'
            IV_DB.append(HOG)
        out=str(HOG)+"\t"+str(SubTrees)+"\t"+str(R)+"\t"+str(spls[3])
        RLN=0
        for i in range(4,11):
            out=out+"\t"+str(spls[i])
            if str(spls[i])=='2(RL)':
                RLN=RLN+1
        MAPS=[]
        for i in range(11,18):
            if R=="I" or R=="II":
                x="-"
                if spls[i]!="-":
                    GD=[]
                    DB=spls[i].strip().split(",")
                    for gene in DB:
                        GID=gene.strip().split("|")[1]
                        gid=GH[GID]
                        MAPS.append(GID)
                        SP=gene.strip().split("|")[0]
                        t=str(SP)+"|"+str(gid)
                        GD.append(t)
                    x=TEXT_NR(GD)
                out=out+"\t"+str(x)
            if R=="III" or R=="IV":
                x="-"
                if spls[i]!="-":
                    GD=[]
                    DB=spls[i].strip().split(";")
                    for g in DB:
                        gd=g.strip().split(": ")[1]
                        gds=[]
                        for gene in gd.strip().split(","):
                            GID=gene.strip().split("|")[1]
                            gid=GH[GID]
                            MAPS.append(GID)
                            SP=gene.strip().split("|")[0]
                            t=str(SP)+"|"+str(gid)
                            gds.append(t)
                        y=TEXT_NR(gds)
                        GD.append(y)
                    x=LIST_NR(GD)
                out=out+"\t"+str(x)
        chy=ANCHECK(MAPS)
        if chy==0:
            print(out,file=fp)   
fp.close()