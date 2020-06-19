      subroutine timchk(ntmax,itime,tmax,t,pr,nblocks)
        implicit none
!***********************************************************************
!     Routine:  timchk.for                                             *
!     Function: Check for time greater that tmax or number of time     *
!               steps greater than ntmax.                              *
!     Written:  January 2001, Summer 2005                              *
!     By:       M. Celia, M.J.Puma                                     *
!***********************************************************************

	  integer*4, parameter :: prec10 = selected_real_kind(10) 
      integer*4, parameter :: nbdim  = 320
      integer*4, parameter :: nuic   = 24
      
	  integer*4  ntmax, itime,i
      real(kind = prec10) tmax, t
      real(kind = prec10) pr(nbdim)   
     
      integer*4  nblocks 

!   Save last solution as possible initial cond. for next run
      if(itime.gt.ntmax.or.t.gt.tmax) then
         write(nuic,30) (pr(i), i=1,nblocks)
   30    format(2x,6e11.4)
         stop
      end if
!
      return
      end