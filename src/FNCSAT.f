      function fncsat(arg,ii,expn,sres,alpha,he,se)
        implicit none
!***********************************************************************
!     Routine:  fncsat.for                                             *
!     Function: Calculate saturation as a function of pr. head.        *
!     Written:  June 2001, Jan. 2004, Summer 2005                      *
!     By:       A. Guswa, M.J. Puma                                    *
!***********************************************************************

	  integer*4, parameter :: prec10 = selected_real_kind(10)
      integer*4, parameter :: nbdim  = 320
	  integer*4, parameter ::  matdim =  5

      real(kind = prec10)fncsat,arg,expn(matdim),sres(matdim),
     +       alpha(matdim),he(matdim),sattmp,se(matdim),bb,hh
      integer*4 ii


!   Brooks and Corey model (1964), which describes the portion of the 
!       curve only for pressure head values less than the bubbling 
!       pressure, or the pressure at which air will enter the soil
      if(arg.ge.0.d0) then
         sattmp = 1.d0
      elseif(arg.gt.he(ii)) then
           sattmp = 1.d0-(1.d0-Se(ii))*(arg/he(ii))
      else
         hh = -alpha(ii)/arg
         bb = 1/expn(ii) != lambda: pore size index
         sattmp = sres(ii)+(1-sres(ii))*hh**bb 
      end if
      fncsat = sattmp
      
!!   Campbell model, which describes the portion of the 
!!       curve only for pressure head values less than the bubbling 
!!       pressure, or the pressure at which air will enter the soil
!
!      if(arg.ge.0.d0) then
!         sattmp = 1.d0
!      elseif(arg.gt.he(ii)) then
!         sattmp = 1.d0-(1.d0-Se(ii))*(arg/he(ii))
!      else
!         hh = -alpha(ii)/arg
!         bb = 1/expn(ii) != lambda: pore size index
!         sattmp = hh**bb 
!      end if
!      fncsat = sattmp
           
      return
      end
