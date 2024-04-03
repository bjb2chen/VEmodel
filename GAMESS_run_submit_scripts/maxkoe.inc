!   -*- fortran -*-
!-----------------------------------------------------------------------
! PARAMETERS FOR ARRAY DIMENSIONS
!
!  maxkoe:   maximum number of hamiltonian terms.
!         If one enlarges maxkoe it will be in general neccessay to
!         enlarge maxhtm and maxhop on operdef.inc to a similar value.
!-----------------------------------------------------------------------
      integer       maxkoe

      parameter     (maxkoe=200000) !PRCS:320000

