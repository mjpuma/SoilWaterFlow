      subroutine pprint3(qin,transum,evsum,qleak,dt,nblocks,t,dz,sat,pr,
     +              poros,stor,prinit,imatl,daysatm,dqin,etmax,evmax,
     +				dtransum,devsum,dleak)
        implicit none
!***********************************************************************
!     Routine:  pprint3.for                                            *
!     Function: Write daily integrated output to file                  *
!     Written:  April 2003, Summer 2005                                *
!     By:       M.J. Puma                                              *
!***********************************************************************
	  integer*4, parameter :: prec10 = selected_real_kind(10)
      integer*4, parameter :: nbdim  = 320
      integer*4, parameter :: matdim =  5
      integer*4, parameter :: nupl4 =  27
            
      real(kind = prec10) dqin,dtransum,devsum,dleak
      real(kind = prec10) qin,transum,evsum,qleak
      real(kind = prec10) fncpor,pr(nbdim),dz(nbdim),sat(nbdim),dt,t
      real(kind = prec10) daysatm(3,nbdim),poros(matdim),stor(matdim),
     +			prinit(nbdim),etmax,evmax
      
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
            daysatm(1,i) = daysatm(1,i)+ fncpor(pr(j),imatl(j),poros,
     +                     stor,prinit(j))*dz(j)
            daysatm(2,i) = daysatm(2,i)+ fncpor(pr(j),imatl(j),poros,
     + 					   stor,prinit(j))*sat(j)*dz(j)
         end do
         if (ndepth.gt.nblocks) EXIT 
         daysatm(3,i) = daysatm(3,i)+(daysatm(2,i)/daysatm(1,i))*dt  
         if(ndepth.eq.1) then
            ndepth = ndepth +4
         else
            ndepth = ndepth + 5
         endif    
      end do  
                                              

!   Sum actual transpiration and stress over day 
      dqin = dqin + qin
      dtransum = dtransum + transum 
      devsum = devsum +evsum
      dleak = dleak + qleak
!      dstatic = dstatic + static !!* dt   
   
!   Output for daily printed file  
      write(nupl4,100) t, daysatm(3,9),daysatm(3,21),
     +                 dqin, dtransum, devsum, dleak,etmax,evmax
  100 format(e20.6,8f10.6)
        
!   Reset averaging variables
      dqin = 0.d0
      dtransum = 0.d0
      devsum = 0.d0
      dleak = 0.d0
!      dstatic = 0.d0
      
      do i = 1,nblocks
         daysatm(3,i) = 0.d0
      enddo

      return
      end
