      function fncevap(arg,ii,time,swilt,sres,zevap,pi,
     +		evmax,etfin,etstart,tday)
        implicit none
!***********************************************************************
!     Routine:  fncevap.for                                            *
!     Function: Calculate the loss of water per unit depth due to      *
!               evaporation using the loss function similar to that    *
!               from Rodriguez-Iturbe et al.                           *
!     Written:  June 2001, February & June 2003                        *
!     By:       A. Guswa, M.J.Puma                                     *
!***********************************************************************

	  integer*4, parameter :: prec10 = selected_real_kind(10)
      integer*4, parameter ::  matdim =  5
      
      real(kind = prec10) fncevap,time,arg,swilt(matdim),sres(matdim),
     +		zevap,pi,evmax(matdim),evmaxd,etfin,etstart,tday,daytime

      integer*4 ii
      

      daytime = time-INT(time/tday)*tday   
!      if ((daytime.le.etstart).or.(daytime.ge.etfin)) then
!         fncevap = 0.d0                            
!      else
!         evmaxd = (evmax(ii)*pi*
!     +          sin(2*pi*(daytime-etstart)))/zevap  


!!      24-hour Uniform Potential Evaporation Rate [cm/day]
         evmaxd = evmax(ii)/zevap        !for 12 hour comment out

         if (arg.lt.sres(ii)) then
            fncevap = 0.d0
         elseif (arg.lt.swilt(ii)) then
            fncevap = evmaxd*((arg-sres(ii))/(swilt(ii)-sres(ii)))
         else
            fncevap = evmaxd
         end if
!      end if

      return
      end 
