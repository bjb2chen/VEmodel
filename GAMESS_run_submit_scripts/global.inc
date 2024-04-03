!   -*- fortran -*-
!***********************************************************************
!
!                         MCTDH include file GLOBAL
!
!  MCTDH include file defining various parameters and variables
!  global to the program.
!
!   V6.0 MB
!
!***********************************************************************

!-----------------------------------------------------------------------
! PARAMETERS FOR CHARACTER STRING LENGTHS
!
! c1, c2, c3, c4, c5, maxkeylen: various lengths
! If you change maxkeylen, also change it in:
! - einglib.f (runinpf)
! - vcglobal.inc
!-----------------------------------------------------------------------
      integer c1, c2, c3, c4, c5, maxkeylen,maxvec

      parameter (c1=16) !PRCS:16
      parameter (c2=32) !PRCS:32
      parameter (c3=64) !PRCS:64
      parameter (c4=80) !PRCS:80
      parameter (c5=240) !PRCS:240
      parameter (maxkeylen=200) !PRCS:200
      parameter (maxvec=32)    !PRCS:32
!-----------------------------------------------------------------------
! progver: code version number, defines output formats etc.
! tiny is added to all these numbers to ensure that version number is
! .ge. to that desired!
! Tiny is also used as a tolerance when real or complex numbers are
! compared.
! filever(i): file version number for file on channel i.
!             Defines output formats etc.
!  NB: filever must be dimensioned larger than largest channel no. used.
!-----------------------------------------------------------------------
      real*8  tiny
      integer maxfver
      parameter (tiny=1.0d-9) !PRCS:1.0d-9
      parameter (maxfver=800) !PRCS:800
      real*8 progver,filever(maxfver)
      common /version/ progver,filever


!-----------------------------------------------------------------------
! CHARACTER STRINGS TO PASS MESSAGES TO ERROR ROUTINE
!-----------------------------------------------------------------------
      character routine*(c3), message*(c5), inptit*(c5)

      common /errorhandling/ routine, message, inptit


!-----------------------------------------------------------------------
! chkerr: array used for checking system definition
!
! chkerr(1)   : Error reading file
! chkerr(2)   : Dofs differ from system
! chkerr(3)   : DVR bases differ from system
! chkerr(4)   : Mode combinations differ from system
! chkerr(5)   : Electronic basis type differs from system
! chkerr(6)   : Spf basis differs from system
! chkerr(7)   : Dofs not subset of system dofs
! chkerr(8)   : ? not for this system
! chkerr(9)   : Psi type differs from system
! chkerr(10)  : Density matrix type differs from system
!
!-----------------------------------------------------------------------
      integer nerr
      parameter (nerr=15)

      integer chkerr(nerr)

      common /rwerror/ chkerr


!-----------------------------------------------------------------------
! DIRECTORY DEFINITION FOR INPUT/OUTPUT FILES
!
! name          :  name of directory in which the output files are stored
! inputfile     : path of the inputfile (~.inp)
! options       : string containing the list of options
! laenge        : length of string name
! inplaenge     : length of string inputfile
! optlaenge     : length of string options
! rname         : directory from which to take restart file
! rlaenge       : length of string rname
! oname         : directory from which to take operator file
! olaenge       : length of string oname
! dname         : directory from which to take dvr file
! dlaenge       : length of string dname
! vname         : name of a vpot file to read
! vlength       : length of string vname
! wname         : name of a file with relevant region weights to be read
! wlength       : length of string wname
! operfile      : path and name of oper file
! opflaenge     : length of string operfile
! crossname     : directory from which to take restart file used for
!                 calculating the cross-correlation-function
! crosslng      : length of string crossname
! progname      : name of program ("$0")
! npotname      : name of natpot file to be read by cutnpot
! npotlaenge    : length of string npotname
!-----------------------------------------------------------------------
      character*(c5)  name,inputfile,rname,dname,oname,vname,operfile
      character*(c5)  options,runtitle,crossname,progname,wname
      character*(c5)  npotname,cstartup
      integer         laenge,inplaenge,optlaenge,rlaenge,dlaenge
      integer         olaenge,vlength,opflaenge,crosslng,lprogname
      integer         wlength, npotlaenge

      common /names/ name,inputfile,options,runtitle,rname,dname,oname
      common /names/ vname,operfile,crossname,progname,wname
      common /names/ npotname,cstartup
      common /lengths/ laenge,inplaenge,optlaenge,rlaenge,dlaenge
      common /lengths/ olaenge,vlength,opflaenge,crosslng,lprogname
      common /lengths/ wlength, npotlaenge

!-----------------------------------------------------------------------
! PARAMETERS FOR CALCULATIONS
!
! pi : pi is defined by 4d0*atan(1d0) in subroutine default.
! fs : time conversion from a.u. to fs (x[a.u.] = fs*x[fs]).
! ev : energy conversion from a.u. to ev (x[a.u.] = x[eV]/ev).
!   Note: If a.u. are output, fs=ev=1
! one: the variable 1
! oner: the variable 1.0d0
! True,False: the logicals .true. and .false., respectively
!-----------------------------------------------------------------------
      integer   one,zeroi,onei,twoi,threei,fouri,fivei,sixi,seveni
      integer   eighti,ninei,teni
      logical   true, false

      real*8    oner
      parameter (oner = 1.0d0)
      parameter (one=1,zeroi=0,onei=1,twoi=2,threei=3,fouri=4)
      parameter (fivei=5,sixi=6,seveni=7,eighti=8,ninei=9,teni=10)
      parameter (true = .true., false = .false.)

      real*8    pi, fs, ev
      common /constants/ pi, fs, ev

!-----------------------------------------------------------------------
! PARAMETERS read from input
!
! npar      : number of parameters read
! apar(n)   : string with name of nth parameter
! rpar(n)   : (real) value of nth parameter
!-----------------------------------------------------------------------
      integer maxpar
      parameter (maxpar=1000000) !PRCS:500000

      integer npar
      real*8 rpar(maxpar)
      character*(c2) apar(maxpar)

      common /iparam/ npar
      common /rparam/ rpar
      common /aparam/ apar

!-----------------------------------------------------------------------
!  POINTERS for the start positions of
!  the arrays stored in the memory arrays mc, mr, ms, mi and ml.
!
! see zeiger routines for definitions
!
!-----------------------------------------------------------------------

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !! If you add any further constant, change the mccount,..,mfcount   !!
! !! values in subroutine checkmemprop !!                             !!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer mcdim,mchpsi,mchunpsi,mcintpsi,mchin,mcrueck,mcexphin
      integer mcexprueck
      integer mchteil,mchloch,mcdicht1,mcdicht2,mcdicht4,mcdicht4o
      integer mchnatur,mcwmat,mcipbt,mcipft
      integer mcworkc,mccdvr2,mccdvr4,mcmfsum,mcoldpsi
      integer mcpredpsi,mcautopsi,mcautodtpsi,mchc,mcspf1d,mcphi
      integer mcpsiref,mcccdover,mcoctpsi
      integer mrdim,mrdicht3,mrort,mrdvrmat,mrtrafo,mrfftp
      integer mrworkr,mrgpop1,mrgpop2,mrcdvr3,mrhunew,mrhunev
      integer mrdicht3o,mrkinsph,mrhops,mrhr,mreigval,mreigerr
      integer mreigint,mrrk5wt,mrintr,ikrylov,ikryc
      integer midim,miworki,mijsph,mimsph,mifftfak,mihi,mijindx
      integer mldim,mlworkl,mlhl
      integer msdim,msworks
      integer mfdim,mfworkf

      common /mcpoint/ mcdim,mchpsi,mchunpsi,mcintpsi,mchin,mcrueck
      common /mcpoint/       mchteil,mchloch,mcdicht1,mcdicht2
      common /mcpoint/       mcdicht4,mchnatur
      common /mcpoint/       mcdicht4o,mcwmat,mcipbt,mcipft,mcworkc
      common /mcpoint/       mccdvr2,mccdvr4,mcmfsum,mcoldpsi
      common /mcpoint/       mcpredpsi,mcautopsi
      common /mcpoint/       mcautodtpsi,mchc,mcspf1d,mcphi
      common /mcpoint/       mcexphin,mcexprueck
      common /mcpoint/       mcpsiref,mcccdover,mcoctpsi

      common /mrpoint/ mrdim,mrdicht3,mrworkr,mrort,mrdvrmat
      common /mrpoint/ mrtrafo,mrfftp,mrdicht3o,mrgpop1
      common /mrpoint/ mrgpop2,mrcdvr3
      common /mrpoint/ mrhunew,mrhunev,mrkinsph,mrhops,mrhr
      common /mrpoint/ mreigval,mreigerr
      common /mrpoint/ mreigint,mrrk5wt,mrintr,ikrylov,ikryc

      common /mipoint/ midim,miworki,mijsph,mimsph,mifftfak
      common /mipoint/ mihi,mijindx

      common /mlpoint/ mldim,mlworkl,mlhl

      common /mspoint/ msdim,msworks

      common /mfpoint/ mfdim,mfworkf

!-----------------------------------------------------------------------
! WORK ARRAY SIZES
!
! workcdim : Length of the complex*16 scratch array
! workrdim : Size of real*8 scratch array
! workldim : Size of logical scratch array
! workidim : Size of integer scratch array
! worksdim : Size of complex*8 scratch array
! workfdim : Size of real*4 scratch array
!-----------------------------------------------------------------------
      integer workcdim,workrdim,workldim,workidim,worksdim,workfdim

      common /arraydim/ workcdim,workrdim,workldim,workidim
      common /arraydim/ worksdim,workfdim

!-----------------------------------------------------------------------
! LOGICALS TO CONTROL INPUT / OUTPUT OPTIONS
!
! If set to .true.,
!     neu      : calculation is a new run
!     cont     : calculation is a continuation run
!     restart  : calculation is a restart run
!     test     : runs up to propagation to test memory etc.
!     test0    : short test-run (via option -s)
!     lovwr    : files in directory are overwritten if neu
!     logisopen: the log file is open
!     lwarn    : is set to 'true' in subroutine 'warnmsg'
!     liopt    : is set to 'true' when input is read from name/input.
!     ltime    : controls timing output
!     lwrlog   : controls printing rst-file information (rstinfo) to log file.
!     lnpot    : flag if natpot files are used
!-----------------------------------------------------------------------
!     tm_zero  : UNIX time (in sec) of start-time of run. (set in mdefaults)
!-----------------------------------------------------------------------
      logical   neu,cont,restart,test,lovwr,lwarn,lnpot,ltime
      logical   test0,logisopen,liopt,lwrlog
      character sstop*8,scpu*10,swall*10
      integer   tm_zero

      common /strtopt/ neu,cont,restart,test,lovwr,lnpot,test0
      common /ioopt/   logisopen,ltime,lwarn,liopt,lwrlog
      common /stopstring/   sstop,scpu,swall
      common /stopunixtime/ tm_zero

!-----------------------------------------------------------------------
!
! FLAGS TO SPECIFY THE RUN-TYPE AND RW-FILE GENERATION
!
!  If set to .true.,
!     lrundvr  : A DVR file will be generated
!     lrunoper : An OPER file will be generated
!     lruninwf : A RESTART file will be generated
!     lrunpes  : A PES file will be generated
!     lrunprop : A propagation will be made
!     lrundiag : if true, the Hamiltonian is diagonalised using Lanczos
!     lrddvr   : A DVR file will be read
!     lrdoper  : An OPER file will be read
!     lrdinwf  : A RESTART file will be read
!     ldeldvr  : DVR file is deleted after use
!     ldelop   : OPER file is deleted after use
!     lexact   : numerically exact calculation
!     lrunml   : A multilayer propagation will be made
!     lreflex  : reflex calculation (malonaldehyde)
!     MGPF
!     ldvronly    : Generate DVR files only, then stop
!-----------------------------------------------------------------------
      logical lrundvr,lrunoper,lruninwf,lrunprop,lrddvr,lrdoper
      logical       lrdinwf,lrunpes,lrundiag,lexact,ldelop,ldeldvr
      logical       lrunml, lreflex, ldvronly

      common /runopts/ lrundvr,lrunoper,lruninwf,lrunprop,lrddvr
      common /runopts/ lrdoper,lrdinwf,lexact,lrundiag,lrunpes
      common /runopts/ ldelop,ldeldvr,lrunml,lreflex
      common /runopts/ ldvronly

!-----------------------------------------------------------------------
! ACCURACY VARIABLES
!
! macheps: machine precision, i. e. smallest positive number for which
!          1.0+macheps > 1.0  (macheps=2.2d-16)
!-----------------------------------------------------------------------
       real*8 macheps
       common /macheps/ macheps

!-----------------------------------------------------------------------
! VARIABLES USED TO READ INPUT FILE
! (defined here so that same size arrays are used in different routines)
!-----------------------------------------------------------------------
      integer maxkey
      parameter (maxkey=100)

      character*(maxkeylen) keyword(maxkey),keyorig(maxkey)
      integer        lc(maxkey)

      common /akey/ keyword,keyorig
      common /ikey/ lc


!-----------------------------------------------------------------------
! MULTILAYER VARIABLES
!-----------------------------------------------------------------------
!      integer maxtape
!      parameter (maxtape=1000)
!
!      integer   tape(maxtape)
!      integer   mladim,mlatopdim,mltotphidim
!
!      common /multilayer/ tape, mladim, mlatopdim, mltotphidim

!-----------------------------------------------------------------------
! array needed to specify which parts of read-write files are to be read
!-----------------------------------------------------------------------
      integer nrwdata
      parameter (nrwdata=15)

      logical dvrdata(nrwdata)

      common /rwlogicals/ dvrdata


!-----------------------------------------------------------------------
! integers for outputcontrol of pdensity:
! npd1...npd4 -- if set -- specify the modes for which the pdensity is
! written.
!-----------------------------------------------------------------------
      integer npd1,npd2,npd3,npd4
      common /npd/ npd1,npd2,npd3,npd4


!-----------------------------------------------------------------------
! DOFs used for gengmat
!-----------------------------------------------------------------------
      integer gmatf1,gmatf2
      common /blkgmat/ gmatf1,gmatf2

!-----------------------------------------------------------------------
! parameter maxnhtmshift:
! maximum number of different htmshift values for mode m.
! (this must be global as it is needed in runprop.inc and operdef.inc)
!-----------------------------------------------------------------------
      integer maxnhtmshift
      parameter (maxnhtmshift=100) !PRCS:100

!-----------------------------------------------------------------------
!       Transformation
!       variables for
!       Polyspherical coordinates to Cartesian coordinates
!       K1 added 2013.07.16
!
!       ltr to check if the surface demands the transformation (again in)
!       nbvec(maxvec) is the number of vectors read (may not needed)
!       nbdumm(maxvec) is the number of the dummy atom the vector passes
!       at1(maxvec) is the atom or dummy atom sitting on the initial
!               point of the vector.
!       at2(maxvec) is the atom or dummy atom sitting on the terminal
!               point of the vector.
!       vecwei(maxvec) defines the position of the dummy atom on the vector
!       icart defines the vector whose cartesian coordinates is needed
!       tord keeps the new order of coordinates (active one) and is consistant with Tnum
!       cart_in cartesian coordinates of the molecule on a gpoint
!       cord changes the order of cartesian coordinates of atoms
!-----------------------------------------------------------------------
      real*8 pfroz(maxvec*3)
      real*8 vecwei(maxvec)
      integer nbvec(maxvec), nbdumm(maxvec),tord(maxvec*3)
      integer cord(maxvec)
      integer icart(maxvec)
      integer frozen(maxvec*3)
      integer nbfr
      character*(1) at1(maxvec),at2(maxvec)
      logical lfr
      logical ltr

      common /trvar/   pfroz,vecwei,nbvec,nbdumm,tord,cord,icart,frozen
      common /trvar/   nbfr,at1,at2,ltr,lfr
