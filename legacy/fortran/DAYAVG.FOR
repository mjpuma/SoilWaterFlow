      subroutine dayavg(qin,transum,evsum,qleak,dt,nblocks,pr,dz,sat,
     +                  poros,stor,prinit,imatl,daysatm,dqin,
     +					dtransum,devsum,dleak)
        implicit none
!***********************************************************************
!     Routine:  dayavg.for                                             *
!     Function: Sums various output variables to obtain daily level    *
!               values                                                 *
!     Written:  July 2003, Summer 2005                                 *
!     By:       M.J. Puma                                              *
!***********************************************************************
	  integer*4, parameter :: prec10 = selected_real_kind(10)
      integer*4, parameter :: nbdim  = 320
      integer*4, parameter ::  matdim =  5
            
      real(kind = prec10) dqin,dtransum,devsum,dleak
      real(kind = prec10) qin,transum,evsum,qleak
      real(kind = prec10) fncpor,dz(nbdim),pr(nbdim),sat(nbdim),dt
      real(kind = prec10) daysatm(3,nbdim),poros(matdim),stor(matdim),
     +			prinit(nbdim) 
      
	  integer*4 ndepth,nblocks,j,i,imatl(nbdim)


!   Initialize daily averaged vertically integrated pore space and 
!   pore water variables for different averaging depths
 
      do i = 1,nblocks
         daysatm(1,i) = 0.d0
         daysatm(2,i) = 0.d0
      enddo   

!   Vertically integrated pore space and water for various 
!   averaging depths and sum vertically integrated saturation over day
                        
      ndepth = 1 
      do i=1,nblocks
         do j = 1,ndepth  
            if (ndepth.gt.nblocks) EXIT
            daysatm(1,i) = daysatm(1,i)+fncpor(pr(j),imatl(j),poros,
     +                     stor,prinit(j))*dz(j)
            daysatm(2,i) = daysatm(2,i)+fncpor(pr(j),imatl(j),poros,
     + 					   stor,prinit(j))*sat(j)*dz(j)
         end do 
         if (ndepth.gt.nblocks) EXIT
         daysatm(3,i) = daysatm(3,i)+(daysatm(2,i)/daysatm(1,i))*dt  
         if(ndepth==1) then
            ndepth = ndepth +4
         else
            ndepth = ndepth + 5
         endif 
      end do   
      
!   Sum actual transpiration and stress over day 
      dqin = dqin + qin
      dtransum = dtransum + transum
      devsum = devsum + evsum
      dleak = dleak + qleak
!      dstatic = dstatic + static !* dt            

      return
      end
