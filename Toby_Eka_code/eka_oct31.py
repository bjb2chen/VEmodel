##

# This code extract normal mode eigen vector

##

import os,sys,math,copy,itertools

 

def read_sample_input(fle_name,natm):

    F=open(fle_name,"r")

    AF=F.readlines()

    F.close()

    N=len(AF)

    for i in range(N):

        LA=AF[i].split()

        if len(LA)>0 and LA[0]=="$data":

            ia=i+3

    AM=[] # Atomic Mass

    for i in range(ia,ia+natm):

        AM.append(AF[i].split()[1])

    return AF[:ia],AF[ia+natm:],AM

def read_ref_geom(fle_name):

    F=open(fle_name,"r")

    A=F.readlines()

    F.close()

    N=len(A)

    Cord=[];Sym=[]

    for i in range(2,N): # Leave first two lines

        FA=A[i].split()

        Sym.append(FA[0]) # 0th column for Atom symbols

        for j in range(1,4): #Last Three Columns for x,y,z coordinates of individual atoms

            Cord.append(float(FA[j]))

    return Cord,Sym

def sumprod(A,B):

    sum0=0.0

    Na=len(A);Nb=len(B)

    for i in range(Na):

        sp=float(A[i])*float(B[i])

        sum0=sum0+sp

    return sum0

def write_new_geom(fle_name,Crd,A,Nl,iA,iB,AB):

    #

    natom=int(len(Crd)/3)

    ngrid=""

    Na=len(Crd)

    for i in range(Na):

        if Nl[i]>0:

            ngrid+="_q%i_pls_%s"%(i,Nl[i])

        if Nl[i]<0:

            ngrid+="_q%i_mns_%s"%(i,abs(Nl[i]))

 

    fout=open(fle_name+ngrid+".xyz","w")

    finp=open(fle_name+ngrid+".inp","w")

    for ai1 in iA:

        finp.write("%s"%ai1)

    fout.write("%i \n"%natom)

    fout.write("\n")

    for i in range(natom):

        fout.write("%s "%A[i])

        finp.write("%s %s "%(A[i],AB[i]))

        for j in range(3):

            fout.write("%12.8f"%Crd[i*3+j])

            finp.write("%12.8f"%Crd[i*3+j])

        fout.write("\n")

        finp.write("\n")

    for bi1 in iB:

        finp.write("%s"%bi1)

    fout.close()

    finp.close()

 

def distort_Nmode(fle_name,mis):

    ##

    natom=3

    F=open(fle_name,"r")

    A=F.readlines()

    F.close()

    Ncoord=len(A)

    Ndim=len(A[0].split())

    q0,ai=read_ref_geom("ref_geom.xyz") # Returns Reference Geometry, in a list (q0) and Atomic Symbons (ai)

    inpA,inpB,MA=read_sample_input("sample_gamessUS.inp",natom)

    qd=[];q1=[]

    for jmis in mis:

        for i in range(Ncoord):

            Ai=A[i].split()

            smpd=sumprod(jmis,Ai)

            qd.append(smpd)

            q1.append(q0[i]+qd[i])

        write_new_geom("H2S",q1,ai,jmis,inpA,inpB,MA)

#------------------------------

ndim=9;lst=[]

for i in range(ndim):

    lst.append([0])

iGrid=[[0,[1]],[7,[-0.1]]]

for ig in iGrid:

    lst[ig[0]]=ig[1]

lst1=list(itertools.product(*lst))

LS0=[]

for l in itertools.product(*lst):

    LS0.append(list(l))

if __name__=="__main__":

    distort_Nmode("Nmode.out",LS0)
