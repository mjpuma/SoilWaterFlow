      subroutine pprint4(dz,sat,t,ztop,zroot,fluxin,fluxet,fluxev)
!***********************************************************************
!     Routine:  pprint4.for                                            *
!     Function: Write saturation profile to file                       *
!     Written:  April 2004, Summer 2005                                *
!     By:       M.J. Puma                                              *
!***********************************************************************

	  integer*4, parameter :: prec10 = selected_real_kind(10)
      integer*4, parameter :: nbdim  = 320
	  integer*4, parameter :: nupl5 =  36 
            
      real(kind = prec10) ztop,depth,zroot,fluxin,fluxet,fluxev
      real(kind = prec10) dz(nbdim),sat(nbdim),t
	  integer*4 i,numzroot

!
!   Output for printed file
      depth = ztop
      numzroot = zroot+10
      do i=1,numzroot
         depth = depth + dz(i)/2
         write(nupl5,100) t, depth, sat(i), fluxin, fluxet, fluxev 
         depth = depth + dz(i)/2
      end do

  100   format(e15.6,e15.4,4f9.4) 

      return
      end
