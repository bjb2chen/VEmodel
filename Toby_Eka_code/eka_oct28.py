# This code extracts normal mode eigen vectors
import os
import sys
import math
 

def extract_mode(fle_name):
    ##
    F=open(fle_name,"r")
    A=F.readlines()
    F.close()
    N=len(A)
    Lnum=[]
    for i in range(N):
        LA=A[i].split()
        if len(LA)>0 and LA[0]=="IR" and LA[1]=="INTENSITY:":
            Lnum.append([i,len(LA)])
        if len(LA)>0 and LA[0]=="TOTAL" and LA[1]=="NUMBER" and LA[2]=="OF" and LA[3]=="ATOMS":
            natom=int(LA[-1])
    ncoord=3*natom
    ndim=ncoord
    md=[]
    for i in range(ndim):
        md.append([])
    count=0
    for i in Lnum:
        ni=-1*i[1]+2
        for j in range(ncoord):
            for k in range(i[1]-2):
                md[k+count].append((A[i[0]+j+2].split()[ni:])[k])
        count += i[1]-2
    fout=open("Nmode.out","w")
    for i in range(ndim):
        for j in range(ncoord):
            fout.write("%12.8f "%float(md[j][i]))
        fout.write("\n")
    fout.close()
 

## Calling Main Function
if __name__=="__main__":
    extract_mode("nh3_ccd_mp2_c3v_gh.out")
