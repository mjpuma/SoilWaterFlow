      subroutine transpir(t,tday,etstart,etfin,pr,pi,rksat,
     + 		etmax,evmax,hwilt,csoil,croot,rtavgdens,dz,
     +		etmult,rweight,nroot,sres,expn,imatl,sat,
     +		trans, sumtrans) 
        implicit none
!***********************************************************************
!     Routine:  transpir.for                                           *
!     Function: Calculate the loss of water per unit depth due to      *
!               transpiration using a loss function similar to that    *
!               described by Cowan 1965                                *
!     Written:  May 2001, Feb. & Summer 2003, Jan. 2004, Summer 2005   *
!     By:       A. Guswa, M.J. Puma                                    *
!***********************************************************************

	  integer*4, parameter :: prec10 = selected_real_kind(10)
      
      integer*4, parameter :: nbdim  = 320
      integer*4, parameter ::  matdim =  5

      integer*4, parameter :: nuout  = 12

      real(kind = prec10)  t
      real(kind = prec10) tday, etstart, etfin
    
      real(kind = prec10) pr(nbdim),sat(nbdim), trans(nbdim)    

      real(kind = prec10) pi,rksat(matdim),sres(matdim),expn(matdim)
      real(kind = prec10) etmax(matdim), evmax(matdim),hwilt(matdim)                
                        
      real(kind = prec10) csoil, croot, rtavgdens(matdim)
      integer*4  imatl(nbdim)

      real(kind = prec10) dz(nbdim), etmult(nbdim)
      real(kind = prec10) rweight(nbdim)
      integer*4  nroot  

	  real(kind = prec10) dhp,f,f_lag,hp_lag,ttol,fnckrel,hp,
     + 					  sumtrans,tpot,errtpot,daytime
      integer*4 jj,j,i	

      daytime = t-INT(t/tday)*tday
      errtpot = 1.0e-8

!!      12-hour Diurnal Potential Transpiration Rate [cm/day]      
!!
!      if((daytime.lt.etstart).or.(daytime.gt.etfin))then
!           tpot = 0.d0
!      else
!           tpot = (etmax(imatl(1)) - evmax(imatl(1)))*
!     +           (pi*sin(2*pi*(daytime-etstart)))
     
!      if (tpot.lt.errtpot) then
         
!         sumtrans = 0.d0 
         
!         do i=1,nroot
!            trans(i) = 0.d0
!         enddo
      
!      else

!      24-hour Uniform Potential Transpiration Rate [cm/day]
         tpot = etmax(imatl(1)) - evmax(imatl(1))  !comment out for 12 hour

         sumtrans = 0.d0
         hp = hwilt(imatl(1))
         do j=1,nroot
              trans(j) = rweight(j) * dz(j)
     +                 * etmult(j) * ((pr(j)-hp) /
     +          ((csoil/(rtavgdens(1)*fnckrel(sat(j),imatl(j),
     +    		expn,sres)*rksat(1))) + (croot/rtavgdens(1))))
            if (trans(j).lt.0) then
               trans(j) = 0
            end if
            sumtrans = sumtrans + trans(j)
         end do
         if (sumtrans.le.tpot) then
            return
         else
            ttol = tpot/20000
            hp_lag = 0
            f_lag = -sumtrans
            f = sumtrans - tpot
            do j = 1,20
               dHp = (hp_lag-hp) * (f/(f-f_lag))
               hp_lag = hp
               f_lag = f
               hp = hp + dHp
               sumtrans = 0
               do jj = 1,nroot 
                  trans(jj) = rweight(jj) * dz(jj)
     +                 * etmult(jj) * ((pr(jj)-hp) /
     +          ((csoil/(rtavgdens(1)*fnckrel(sat(jj),imatl(jj),
     +			expn,sres)*rksat(1))) + (croot/rtavgdens(1))))
                  if (trans(jj).lt.0) then
                     trans(jj) = 0
                  end if
                  sumtrans = sumtrans + trans(jj)
               end do
               f = sumtrans - tpot 
               if (dabs(f).lt.ttol) then
                  return
               end if
            end do
            write(nuout,100)  
            end if
!         end if
!      end if 

  100 format(/'Maximum iterations exceeded in transpiration function')

      return
      end
