      function fnckrel(arg,ii,expn,sres)
        implicit none
!***********************************************************************
!     Routine:  fnckrel.for                                            *
!     Function: Calculate relative permeability as a                   *
!               function of saturation.                                *
!     Written:  June 2001, January 2004, Summer 2005                   *
!     By:       A.Guswa, M.J. Puma                                     *
!***********************************************************************
	  integer*4, parameter :: prec10 = selected_real_kind(10)
      integer*4, parameter ::  matdim =  5

      real(kind = prec10) fnckrel,arg,expn(matdim),sres(matdim)
      integer*4 ii

!!   van Genuchten model (1980) 
!      paramn= 1.5
!      paramm = 1 - 1/paramn
!      effsat = (arg-sres(ii))/(1-sres(ii))
!      fnckrel = (effsat**0.5)*(1 - (1 - effsat**(1/paramm))**paramm)**2
!      if (fnckrel.lt.1.0e-11) then
!        fnckrel = 1.0e-11
!      end if


!!   Brooks and Corey model (1964)
!
        fnckrel = ((arg-sres(ii))/(1-sres(ii)))**(2.d0*expn(ii)+3)
        if (fnckrel.lt.1.0e-11) then
          fnckrel = 1.0e-11
        end if

!   Campbell model 
!
!        fnckrel = arg**(2.d0*expn(ii)+3)
!        if (fnckrel.lt.1.0e-11) then
!          fnckrel = 1.0e-11
!        end if

      return
      end
