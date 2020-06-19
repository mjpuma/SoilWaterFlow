      subroutine pprint2(pr,dz,sat,t,ztop,nblocks,poros,prinit,
     + 					stor,imatl)
         implicit none
!***********************************************************************
!     Routine:  pprint2.for                                            *
!     Function: Output saturation profile (vertical)                   *
!     Written:  Summer 2005                      					   *
!     By:       M.J. Puma                                    		   *
!***********************************************************************

	  integer*4, parameter :: prec10 = selected_real_kind(10)
      integer*4, parameter ::  matdim =  5
      integer*4, parameter :: nbdim  = 320
      integer*4, parameter :: nupl2  = 23
                  
      real(kind = prec10) ztop,depth
      real(kind = prec10) fncpor,pr(nbdim),dz(nbdim),sat(nbdim),t
	  integer*4 nblocks,i,imatl(nbdim)
      
  	  real(kind = prec10) poros(matdim)      
      real(kind = prec10) prinit(nbdim),stor(matdim)

!   Output saturation profile
      depth = ztop
      do i=1,nblocks
         depth = depth + dz(i)/2
         write(nupl2,100) t, depth, fncpor(pr(i),imatl(i),
     +       			poros,stor,prinit(i)), sat(i)
         depth = depth + dz(i)/2
      end do

  100   format(e15.6,e15.4,2f9.4) 

      return
      end