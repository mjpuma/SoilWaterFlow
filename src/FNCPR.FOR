      function fncpr(arg,ii,expn,sres,alpha,he,se)
        implicit none
!***********************************************************************
!     Routine:  fncpr.for                                              *
!     Function: Calculate pressure head as a function of saturation.   *
!     Written:  June 2001, January 2004, Summer 2005                   *
!     By:       A. Guswa, M.J. Puma                                     *
!***********************************************************************

	  integer*4, parameter :: prec10 = selected_real_kind(10)
      integer*4, parameter ::  matdim =  5

      real(kind = prec10) fncpr,arg,expn(matdim),sres(matdim),
     +       alpha(matdim),he(matdim),prtmp,se(matdim)
      integer*4 ii

!   Brooks and Corey model (1964), which describes the portion of the 
!       curve only for pressure head values less than the bubbling 
!       pressure, or the pressure at which air will enter the soil
      if(arg.gt.Se(ii)) then
         prtmp = he(ii)*(1-arg)/(1-Se(ii))
      else
         prtmp = -alpha(ii)*((arg-sres(ii))/(1-sres(ii)))**(-expn(ii))
      end if
      
      fncpr = prtmp

!!   Campbell model, which describes the portion of the 
!!       curve only for pressure head values less than the bubbling 
!!       pressure, or the pressure at which air will enter the soil
!!
!      if(arg.gt.Se(ii)) then
!         prtmp = he(ii)*(1-arg)/(1-Se(ii))
!      else
!         prtmp = -alpha(ii)*arg**(-expn(ii))
!      end if
      
!      fncpr = prtmp
        
      return
      end
