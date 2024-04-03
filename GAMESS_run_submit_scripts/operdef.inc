!    -*- fortran -*-
!***********************************************************************
!
!              MCTDH include file  operdef.inc
!
!***********************************************************************


      integer maxhtm,maxfac,maxmuld
      integer maxhpar,maxhop,maxpes,maxLMR
      parameter (maxhtm=200000) !PRCS:320000
      parameter (maxhop=200000) !PRCS:320000
      parameter (maxfac=800) !PRCS:800
      parameter (maxpes=800) !PRCS:800 ! maxpes.ge.maxrdf must hold!
      parameter (maxhpar=8) !PRCS:8 ! maxrdf is defined in runoper.inc
      parameter (maxmuld=64) !PRCS:64
      parameter (maxLMR=120) !PRCS:120

!-----------------------------------------------------------------------
! parameters defining operator modes (may be different from system)
!-----------------------------------------------------------------------
      integer nopmode,nopdof

      common /opmodes/ nopmode,nopdof

!-----------------------------------------------------------------------
! POINTERS FOR ARRAYS
!-----------------------------------------------------------------------
      integer zhop(maxhtm)

      common /operpoint/ zhop

!-----------------------------------------------------------------------
! the hterm and associated arrays
!-----------------------------------------------------------------------
      integer hterm(maxfac,maxhtm),htmdof(maxhtm),htmmode(maxhtm)
      integer htmmuld(maxhtm),htmfac(maxhtm),htmform(maxhtm)
      integer htmsym(maxhtm),htmadof(maxhtm),htmshift(maxhtm)
      integer htmndof(maxhtm)

      common /htmdef/ hterm,htmdof,htmmode,htmmuld,htmfac,htmform
      common /htmdef/ htmsym,htmadof,htmshift,htmndof

!-----------------------------------------------------------------------
! SOME MORE PARAMETERS NEEDED FOR HAMILTONIAN
!
! nopts(n)  : number of parameters storred in pesopts(n) for nth surface
!-----------------------------------------------------------------------
      integer hopfile(maxhop),hopdim,hopipar(maxhpar,maxhop)
      integer hopilab(maxhop),hopidef(maxhop),npesopt,nhoppar
      integer zpesopts(2*maxhpar,maxpes),nopts(maxpes)
      real*8 hoprpar(maxhpar,maxhop)
      real*8 hoppar(maxhpar,maxhop)
      character*(c2) hoplab(maxhop), hoplab0
      character*(c5) pesopts(maxpes)

      common /ihoplab/ hopfile,nhoppar,hopdim,hopilab,hopipar
      common /ihoplab/ npesopt,zpesopts,nopts,hopidef
      common /ahoplab/ hoplab0,hoplab,pesopts
      common /rhoplab/ hoprpar,hoppar

!-----------------------------------------------------------------------
! parameters used to cut high/low values of correlated PES surface,
! and parameters for natpotcut.
!-----------------------------------------------------------------------
      logical lcutpes,lsoftcutpes,lncutall
      integer ncutted(maxpes),nremoved,ncutn
      real*8  pesmin,pesmax,cutoff,natpotcut,npotcut(maxpes)
      real*8  ncutrlabel(maxpes),wpf(maxdim,maxpes)
      character*(c2) ncutlabel(maxpes)

      common /lcutpescb/ lcutpes,lsoftcutpes,lncutall
      common /icutpes / ncutted,nremoved,ncutn
      common /rcutpes/  pesmin,pesmax,cutoff,natpotcut,npotcut
      common /rcutpes/  ncutrlabel,wpf
      common /ccutpes/  ncutlabel

!-----------------------------------------------------------------------
! properties of the PES surface
!
! iqconst(f) : if 1, and dof is not in specified system, coordinate
!              is set to value in rqconst(f)
!              if 0, and dof is not in specified system, coordinate
!              is set to 0.0d0
!              if dof is in specified system, this flag is ignored.
!-----------------------------------------------------------------------
      integer iqconst(maxdim)
      real*8 vmin,vmax,rqconst(maxdim)

      common /ipesprop/ iqconst
      common /rpesprop/ vmin,vmax,rqconst

!-----------------------------------------------------------------------
! arrays storing information for multi-D operators
!
! mulddim: no. of multi-D sets of coordinates
!
! ldirmuld: if true, direct algorithm used for multi-dimensional
!           operators
!-----------------------------------------------------------------------
      integer mulddim,mulddof(maxdim,maxmuld),muldsize(maxmuld)
      integer nmulddof(maxmuld),vdmuld(maxmuld),ndmuld(maxmuld)
      integer mulddof2(maxdim,maxmuld)
      logical ldirmuld

      common /mulddef/ mulddim,mulddof,muldsize,nmulddof,vdmuld,ndmuld
      common /mulddef/ mulddof2
      common /mulddefl/ ldirmuld

!-----------------------------------------------------------------------
! PARAMETERS NEEDED FOR HAMILTONIAN
!
! hopsdim : length of array containing hops
!  lreal  : Hamiltonian and WF are real (improved relaxation).
!  lcap   : calculation employing CAPs (or complex potentials)
!  lcdvr  : any multidimensional Hamiltonian terms use the CDVR method
!  ltdvr  : any multidimensional Hamiltonian terms use the TDVR method
!  lanalpes : an analytic PES will be generated (e.g. for CDVR)
!  lwrnpot : Full information on natpots will be printed to op.log
!-----------------------------------------------------------------------
      integer      hopsdim,ktot
      logical      lreal,lrreal,lcap,lcdvr,ltdvr,lanalpes,lwrnpot

      common /iopinitg/ hopsdim,ktot
      common /lopinitg/ lcap,lcdvr,ltdvr,lanalpes,lreal,lrreal,lwrnpot

!-----------------------------------------------------------------------
! SOME MORE PARAMETERS NEEDED FOR HAMILTONIAN
!
! kzahl     : number of correlated Hamiltonian expansion coefficients
! koeff(k)  : expansion coefficients of correlated Hamiltonian
! ki(k),kf(k): defines states for Hamiltonian term <kf|H|ki>
! complete(m) : if true, spf basis is "complete" and (1-P)H phi =0
! diag(m,k)   : array of unitary spf hamiltonian matrices
! hsym(m,k)   : array of operator symmetries
! hlsym(m,k)  : array of hloch symmetries
! diagh(m)    : auxiliary diag array
! kinop(f,k)  : =1 if a correlated Hamiltonian term is non-diagonal,
!               =0 if diagonal,
!               <0 if side-diagonal.
! laddmode(m,nh)  : logical set to true if the mode m is to be included
!                   into the operator nh. (See keyword 'addmode').
! khzahl(m,s)  : number of uncorrelated coefficients for dimension m
!                and state s
! hamlab(nham) : string with name of nham-th Hamiltonian.
! hamilton(f,k) : integer with no. of kth Hamiltonian term for dof f.
! hamLMR(i,n)  : integer with no. of ith LMR term. n=1,2,3 -> R,M,L.
! zhun(m,s) : pointer to start of uncorrelated terms in hamilton array
! nhtmshift(m) : number of different htmshift values for mode m
! ihtmshift(i,m) : i-th htmshift value for mode m
!-----------------------------------------------------------------------
      integer maxham
      parameter (maxham=64) !PRCS:64

      integer kzahl(maxham),nmulpot(maxham),khtot(maxham)
      integer ki(maxkoe),kf(maxkoe),kinop(maxdim,maxkoe)
      integer nhtmshift(maxdim),ihtmshift(maxnhtmshift,maxdim)
      complex*16 koeff(maxkoe)
      logical complete(maxdim),diag(maxdim,maxkoe),diagh(maxdim)
      logical laddmode(maxdim,maxham),allcomplete

      integer khzahl(maxdim,maxsta,maxham),hamLMR(maxLMR,3),iLMR
      integer zhun(maxdim,maxsta,maxham),hamilton(maxdim,maxkoe)
      integer htmtype(maxhtm),htmdim,zham(maxham),hamdim,h00
      integer hamnum(maxkoe),hsym(maxdim,maxkoe),hlsym(maxdim,maxkoe)
      character*(c2) hamlab(maxham)

      common /iop/ kzahl,khzahl,hamilton,htmtype,zhun,khtot
      common /iop/ nmulpot,ki,kf,htmdim,zham,hamdim,hamnum,hsym
      common /iop/ hlsym,kinop,nhtmshift,ihtmshift,hamLMR,iLMR,h00
      common /cop/ koeff
      common /lop/ complete,diag,diagh,laddmode,allcomplete
      common /aop/ hamlab

!-----------------------------------------------------------------------
! flag to turn off use of diag flag and extraction of huncorr terms.
!-----------------------------------------------------------------------
      logical nodiag(maxham),usepack(maxham)
      common /hamflags/ nodiag, usepack

!-----------------------------------------------------------------------
! store the size of hteil for different hamiltonians
!-----------------------------------------------------------------------
      integer nhpmatdim(maxham)
      common /pmatsizes/ nhpmatdim

!-----------------------------------------------------------------------
!
! Flags for operator action

      integer   kfirst(maxdim,maxkoe)

      common /iop1/ kfirst

!-----------------------------------------------------------------------
! Variables needed for autocap option

! countpsi : number of output of psi in psi file
! autocap(k) : if > 0, k the term includes an autocap on dof f=autocap(k).
!              The cap is then activated when it is populated.
! indicap(f) : number of grid point where the CAP starts
! gpopcap(f) : grid population at the point where the CAP starts
! timecap(f) : time, after which the acap is allowed to be switched on.
!-----------------------------------------------------------------------
      integer countpsi,indicap(maxdim),autocap(maxkoe)
      real*8  gpopcap(maxdim),timecap(maxdim)

      common /iautcap/ countpsi,indicap,autocap
      common /rautcap/ gpopcap,timecap

!-----------------------------------------------------------------------
! Information about natpot terms:
! modc          - the contracted mode of the natpots in Potfit
! vpdim(m)      - Cumulative number of preceeding natural potentials for
!                 mode m in Potfit
! pdim(m)       - Number of natural potentials for mode m in Potfit
! npdim(m)      - Cumulative number of following natural potentials for
!                 mode m in Potfit
! maxsubpdim    - maximum size of a sub-particle grid (for diagonal
!                 terms acting on some but not all of the mode's DOFs)
!-----------------------------------------------------------------------

      integer  modc,vpdim(maxdim),pdim(maxdim),npdim(maxdim),maxsubpdim
      common /natdata/ modc,vpdim,pdim,npdim,maxsubpdim

!-----------------------------------------------------------------------
! Information need by "fast" feature
!
! lfast(n)      - true, if FAST algorithm for H(npot)*A is used
! lfastall      - true if FAST is used for all natpots
! nterm(n)      - the number of the first natpot term of the n-th natpot
! ntermlast(n)  - the number of the last natpot term of the n-th natpot
! natpmode(i,n) - mode for which natpot terms are pre-multiplied for the
!                 n-th natpot (i-th order "fast")
! mmodc(n)      - the contracted mode of the n-th natpot in MCTDH
! vmdim(m,n)    - Cumulative number of preceeding natural potentials for
!                 mode m for the n-th natpot in MCTDH
! mdim(m,n)     - Number of natural potentials for mode m for the n-th
!                 natpot in MCTDH
! nmdim(m,n)    - Cumulative number of following natural potentials for
!                 mode m for the n-th natpot in MCTDH
! fastorder(n)  - Order of the "fast" algorithm for the n-th natpot
! fastprem(i)   - Array of pre-multiplied terms
! fastnum(i)    - Array of terms to be taken with a given pre-multiplied
!                 term.
! fastdim(i)    - Number of all terms to be taken with a given
!                 pre-multiplied term.
! mdimf(m,n)    - Number of pre-multiplied terms for mode m and natpot n
! zfp(n)        - Pointer to fastprem array for given n
! zfn(n)        - Pointer to fastnum array for given n
! zfd(n)        - Pointer to fastdim array for given n
! fastn         -
! fastn0        -
! fastilabel    -
! fastlabel     - Labels of Natpots that use FAST
!-----------------------------------------------------------------------
      logical lfast(maxpes),lfastall
      common /fastl/ lfast,lfastall

      integer nterm(maxpes),ntermlast(maxpes),mmodc(maxpes)
      integer natpmode(maxdim,maxpes)
      integer vmdim(maxdim,maxpes),mdim(maxdim,maxpes)
      integer nmdim(maxdim,maxpes),fastorder(maxpes)
      integer fastprem(maxkoe),fastnum(maxkoe),fastdim(maxkoe)
      integer zfp(maxpes),zfn(maxpes),zfd(maxpes)
      integer mdimf(maxdim,maxpes),fastn,fastn0,fastilabel(maxpes)
      common /fasti/ nterm,ntermlast,mmodc,natpmode,vmdim,mdim,nmdim
      common /fasti/ fastorder,fastprem,fastdim,fastnum,mdimf,zfp,zfn
      common /fasti/ zfd,fastn,fastn0,fastilabel

      character*(c2) fastlabel(maxpes)
      common /fastc/ fastlabel

!-----------------------------------------------------------------------
! SOME MORE PARAMETERS NEEDED FOR EXPECTATION VALUES
!
! nexpect     : number of operators to derive the expectation-value
! expectlab() : labels for the expectation operators
! lreonly     : logical, only real parts will be output to expectation
! nexpect1    : number of operators for file expect1
! lreonly1    : logical, only real parts will be output to expect1
!-----------------------------------------------------------------------
      integer        nexpect, nexpect1
      logical        lreonly, lreonly1
      complex*16     expvalue(maxham)
      character*(c2) expectlab(maxham)

      common /cexpectation/ expectlab
      common /zexpectation/ expvalue
      common /iexpectation/ nexpect, nexpect1
      common /lexpectation/ lreonly, lreonly1
