      function fncpor(arg,ii,poros,stor,arginit)
        implicit none
!***********************************************************************
!     Routine:  fncpor.for                                             *
!     Function: Calculate porosity as a function of pr. head, using    *
!               storage coefficient to track changes in porosity.      *
!     Written:  February 2001, Summer 2005                             *
!     By:       M. Celia, M.J. Puma                                    *
!***********************************************************************
	  integer*4, parameter :: prec10 = selected_real_kind(10)
      integer*4, parameter :: nbdim  = 320
	  integer*4, parameter ::  matdim =  5

      real(kind = prec10) poros(matdim),stor(matdim),
     +		fncpor,arg,arginit
      integer*4 ii
      

      fncpor = poros(ii)+stor(ii)*(arg-arginit)

      return
      end
