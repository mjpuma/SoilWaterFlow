      subroutine pprint(transum,evsum,qleak,dt,nblocks,t,tday,pi,pr,
     +		etstart,fluxin,flxout,fluxet,fluxev,ratio,ratio1,etmax,       
     +		poros,prinit,stor,imatl,dz,sat,evmax)
        implicit none
!***********************************************************************
!     Routine:  pprint.for                                             *
!     Function: Write instantaneous output to print and plot files.    *
!     Written:  February 2001, April 2003, Summer 2005                 *
!     By:       A. Guswa, M.J. Puma                                    *
!***********************************************************************
	  
      integer*4, parameter :: prec10 = selected_real_kind(10)
      integer*4, parameter :: nbdim  = 320
      integer*4, parameter :: matdim =  5
      integer*4, parameter :: nupl1  = 22
      integer*4, parameter :: nupl3  = 26
            
      real(kind = prec10) transum,evsum,qleak,dt,t,tday
      real(kind = prec10) fncpor,pr(nbdim),dz(nbdim),sat(nbdim)
      real(kind = prec10) satinstm(3,nbdim) 
      
      real(kind = prec10) daytime,tpot,etmax(matdim),evmax(matdim),
     +                    pi,etstart,evpot
      real(kind = prec10) fluxin,flxout,fluxet,fluxev,ratio,ratio1 
      real(kind = prec10) errtpot,rateTR,rateEV,rateL,etevnorm

  	  real(kind = prec10) poros(matdim)      
      real(kind = prec10) prinit(nbdim),stor(matdim)
      
	  integer*4 ndepth,nblocks,j,i,imatl(nbdim) 
      
      do i = 1,nblocks
         satinstm(1,i) = 0.d0
         satinstm(2,i) = 0.d0
         satinstm(3,i) = 0.d0         
      enddo     
          
      ndepth = 1 
      do i=1,nblocks
         do j = 1,ndepth  
            if (ndepth.gt.nblocks) EXIT
            satinstm(1,i)=satinstm(1,i)+fncpor(pr(j),imatl(j),
     +       			poros,stor,prinit(j))*dz(j)
            satinstm(2,i)=satinstm(2,i)+fncpor(pr(j),imatl(j),
     +       			poros,stor,prinit(j))*sat(j)*dz(j)
         end do
         if (ndepth.gt.nblocks) EXIT 
         satinstm(3,i) = satinstm(2,i)/satinstm(1,i)  
         if(ndepth.eq.1) then
            ndepth = ndepth +4
         else
            ndepth = ndepth + 5
         endif
      enddo
      
      daytime = t-INT(t/tday)*tday
	  tpot = etmax(imatl(1)) - evmax(imatl(1))  !comment out for 12 hour
	  evpot = evmax(imatl(1))

!      tpot = (etmax(imatl(1)) - evmax(imatl(1)))*
!     +           (pi*sin(2*pi*(daytime-etstart))) 
!      evpot = evmax(imatl(1))*(pi*sin(2*pi*(daytime-etstart)))
      errtpot = 1.0e-8 
      
      rateTR = transum/dt
      rateEV = evsum/dt
      rateL = qleak/dt
      if (evpot.lt.errtpot) evpot = 0.d0
          
      if (tpot.lt.errtpot) then
         tpot = 0.d0  
         etevnorm = 0.d0
      else
         etevnorm = (rateTR + rateEV)/(tpot + evpot)                      
      endif 
      
!   Output for printed file                                           
      write(nupl1,100) t,etevnorm,tpot,rateTR,evpot,rateEV,rateL, 
     +                 satinstm(3,9),satinstm(3,21)                

      write(nupl3,120) t, fluxin, flxout, fluxet, fluxev, 
     +                   ratio, ratio1
  100 format(e20.6, 8f10.4)
  120 format(7e14.6) 
      
      rateTR = 0.d0
      rateEV = 0.d0
      rateL = 0.d0

      return
      end
