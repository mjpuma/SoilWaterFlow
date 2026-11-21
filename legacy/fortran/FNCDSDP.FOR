      function fncdsdp(arg,ii,alpha,sres,expn,se,he)
        implicit none
!***********************************************************************
!     Routine:  fncdsdp.for                                            *
!     Function: Calculate derivative of saturation with respect to     *
!               pressure (h in this case.)                             *
!     Written:  July 1990,  May 1992, Nov. 2000, Jan. 2001, Jan. 2004, *
!               Summer 2005  										   *
!     By:       M. Celia, A. Guswa, M.J. Puma                          *
!***********************************************************************

	  integer*4, parameter :: prec10 = selected_real_kind(10)
      integer*4, parameter ::  matdim =  5
      
      real(kind = prec10) fncdsdp,arg,expn(matdim),lambda,dsdptmp,
     +          alpha(matdim),sres(matdim),Se(matdim),he(matdim)
     
      integer*4 ii

!!   Brooks and Corey model (1964)
!!
      if(arg.ge.0.d0) then
         dsdptmp = 0.d0
      elseif(arg.gt.he(ii)) then
           dsdptmp = (1-Se(ii))/(-he(ii))
      else
         lambda= 1/expn(ii)
         dsdptmp = (1-sres(ii))*lambda/alpha(ii)*(arg/(-alpha(ii)))
     +          **(-1-lambda)   
      end if
      fncdsdp = dsdptmp

!!   Campbell
!!
!      if(arg.ge.0.d0) then
!         dsdptmp = 0.d0
!      elseif(arg.gt.he(ii)) then
!         dsdptmp = (1-Se(ii))/(-he(ii))
!      else
!         lambda=1/expn(ii)
!         dsdptmp=lambda/alpha(ii)*(arg/(-alpha(ii)))**(-1-lambda)
!	  end if
!
!      fncdsdp = dsdptmp


!   van Genuchten (1980) form (check if this is correct)
!      parama = 0.014
!      paramn= 1.5
!      paramm = 1 - (1/paramn)
!      htheta = (1/(1+(parama*abs(arg))**paramn))**paramm
!      if(arg.ge.0.d0) then
!         dsdptmp = 0.d0
!      else
!         dsdptmp = ((-parama*paramm*(1-sres(ii)))/(1-paramm))*
!     +           (htheta**(1/paramm))*((1-(htheta**(1/paramm)))
!     +                  **paramm)
!      endif
!      fncdsdp = -dsdptmp  

      return
      end
