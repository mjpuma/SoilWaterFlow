      subroutine matrix(dt,t,satold,trans,rksat,stor,dz,
     + 		evmult,eweight,nblocks,nbct,nbcb,bct,bcb,sat,
     +		pr,prinit,prold,expn,sres,alpha,he,se,zevap,pi,
     +	 	evmax,etfin,etstart,tday,swilt,imatl,poros,
     +		aa,rhs)
        implicit none
!***********************************************************************
!     Routine:  matrix.for                                             *
!     Function: Generate matrix entries for unsaturated flow equation  *
!               using block-centered finite difference formulation.    *
!  	  Notes:    Uses harmonic average for interblock saturated 		   *
!				conductivities, interblock relative permeabilities     *
!				are determined with an upwind scheme	 			   *
!     Written:  December 2000, February 2003,Summer 2005               *
!     By:       M. Celia, M.J. Puma                                    *
!***********************************************************************
	  integer*4, parameter :: prec10 = selected_real_kind(10)

      integer*4, parameter :: nbdim  = 320
      integer*4, parameter :: ntdim  = 14020
      integer*4, parameter ::  matdim =  5
     
      real(kind = prec10) bct, bcb
      real(kind = prec10) dt, t
      real(kind = prec10) aa(nbdim,3), rhs(nbdim)
      real(kind = prec10) pr(nbdim), prold(nbdim),
     +               sat(nbdim), satold(nbdim), trans(nbdim)    
      real(kind = prec10) rksat(matdim), stor(matdim)
      integer*4  imatl(nbdim)

      real(kind = prec10) dz(nbdim), evmult(nbdim),eweight(nbdim)
      integer*4  nblocks,nbct,nbcb
      
      real(kind = prec10) fncsat,sat0,pr0,fncevap,fncdsdp,diff2,
     +		xkrel2,delh2,dz2,xksat2,diff1,fnckrel,xkrel1,
     +		coefkr,delh1,dz1,xksat1,porold,por,fncpor
      integer*4 ne,i   

  	  real(kind = prec10) poros(matdim)      
      real(kind = prec10) prinit(nbdim),expn(matdim),sres(matdim),
     +       alpha(matdim),he(matdim),se(matdim),zevap,pi,
     +	 	evmax(matdim),etfin,etstart,tday,swilt(matdim) 

!   Loop over interior blocks
      do i=2,nblocks-1
         por = fncpor(pr(i),imatl(i),poros,stor,prinit(i))
         porold = fncpor(prold(i),imatl(i),poros,stor,prinit(i))
!         por = poros(imatl(i))
!         porold = por
         xksat1 = rksat(imatl(i-1))*rksat(imatl(i))*( dz(i-1)+dz(i) )
     +          / ( rksat(imatl(i))*dz(i-1)+rksat(imatl(i-1))*dz(i) )
           dz1 = (dz(i-1)+dz(i)) / 2.d0
           delh1 = (pr(i)-pr(i-1))/dz1 -1
           if (delh1.le.-1) then
             coefkr = 1
           elseif (delh1.ge.1) then
             coefkr = 0
           else
             coefkr = 0.5 - delh1/2.d0
           end if 
           xkrel1 = fnckrel(sat(i-1),imatl(i-1),expn,sres)*coefkr + 
     +            fnckrel(sat(i),imatl(i),expn,sres)*(1-coefkr)
!         xkrel1 = ( fnckrel(sat(i-1),imatl(i-1))*dz(i)
!     +              + fnckrel(sat(i),imatl(i))*dz(i-1) )
!     +            / ( dz(i-1)+dz(i) )
         diff1 = xksat1*xkrel1/dz1
         xksat2 = rksat(imatl(i))*rksat(imatl(i+1))*( dz(i)+dz(i+1) )
     +          / ( rksat(imatl(i+1))*dz(i)+rksat(imatl(i))*dz(i+1) )
           dz2 = (dz(i)+dz(i+1)) / 2.d0
           delh2 = (pr(i+1)-pr(i))/dz2 -1
           if (delh2.le.-1) then
             coefkr = 1
           elseif (delh2.ge.1) then
             coefkr = 0
           else
             coefkr = 0.5 - delh2/2.d0
           end if 
           xkrel2 = fnckrel(sat(i),imatl(i),expn,sres)*coefkr + 
     +            fnckrel(sat(i+1),imatl(i+1),expn,sres)*(1-coefkr)
!         xkrel2 = ( fnckrel(sat(i),imatl(i))*dz(i+1)
!     +              + fnckrel(sat(i+1),imatl(i+1))*dz(i) )
!     +            / ( dz(i)+dz(i+1) )
         diff2 = xksat2*xkrel2/dz2
         aa(i,1) = -diff1
         aa(i,2) = dz(i)/dt * ( por*fncdsdp(pr(i),imatl(i),alpha,sres,
     +				expn,se,he)+ stor(imatl(i))*sat(i) )
     +             + diff1 + diff2
         aa(i,3) = -diff2
!          revap = fncevap(sat(i),imatl(i),t)*evmult(i)*dz(i)
         rhs(i) = -dz(i)/dt * ( por*sat(i)-porold*satold(i) )
     +          - diff1*(pr(i)-pr(i-1)-dz1)
     +          + diff2*(pr(i+1)-pr(i)-dz2)
     +          - trans(i)
     +          - fncevap(sat(i),imatl(i),t,swilt,sres,zevap,pi,
     +		evmax,etfin,etstart,tday)*evmult(i)*eweight(i)*dz(i)
      end do


!   BOUNDARY CONDITIONS
!   1) Top block (i=1):
      por = fncpor(pr(1),imatl(1),poros,stor,prinit(1))
      porold = fncpor(prold(1),imatl(1),poros,stor,prinit(1))
!      por = poros(imatl(1))
!      porold = por
      xksat2 = rksat(imatl(1))*rksat(imatl(2))*( dz(1)+dz(2) )
     +       / ( rksat(imatl(2))*dz(1)+rksat(imatl(1))*dz(2) )
      dz2 = (dz(1)+dz(2)) / 2.d0
        delh2 = (pr(2)-pr(1))/dz2 -1
        if (delh2.le.-1) then
             coefkr = 1
        elseif (delh2.ge.1) then
             coefkr = 0
        else
             coefkr = 0.5 - delh2/2.d0
        end if 
        xkrel2 = fnckrel(sat(1),imatl(1),expn,sres)*coefkr + 
     +         fnckrel(sat(2),imatl(2),expn,sres)*(1-coefkr)
!      xkrel2 = ( fnckrel(sat(1),imatl(1))*dz(2)
!     +           + fnckrel(sat(2),imatl(2))*dz(1) )
!     +         / ( dz(1)+dz(2) )
      diff2 = xksat2*xkrel2/dz2
!
      if(nbct.eq.1) then
         pr0 = 2.d0*bct-pr(1)
!         sat0 = fncsat(pr0,imatl(1))
         sat0 = fncsat(bct,imatl(1),expn,sres,alpha,he,se)
         xksat1 = rksat(imatl(1))
           dz1 = dz(1)
           delh1 = (pr(1)-pr0)/dz1 -1
         if (delh1.le.-1) then
             coefkr = 1.d0
           elseif (delh1.ge.1) then
             coefkr = 0.d0
           else
             coefkr = 0.5 - delh1/2.d0
           end if 
           xkrel1 = fnckrel(sat0,imatl(1),expn,sres)*coefkr + 
     +            fnckrel(sat(1),imatl(1),expn,sres)*(1.0-coefkr)
!         xkrel1 = ( fnckrel(sat(1),imatl(1))
!     +              + fnckrel(sat0,imatl(1)) ) / 2.d0
         diff1 = xksat1*xkrel1/dz1
         aa(1,1) = 0.d0
         aa(1,2) = dz(1)/dt * ( por*fncdsdp(pr(1),imatl(1),alpha,
     +    		   sres,expn,se,he) + stor(imatl(1))*sat(1))
     +             + 2.d0*diff1 + diff2
         aa(1,3) = -diff2
         rhs(1) = -dz(1)/dt * ( por*sat(1)-porold*satold(1) )
     +          - diff1*(pr(1)-pr0-dz1)
     +          + diff2*(pr(2)-pr(1)-dz2)
     +          - trans(1)
     +          - fncevap(sat(1),imatl(1),t,swilt,sres,zevap,pi,
     +		evmax,etfin,etstart,tday)*evmult(1)*eweight(1)*dz(1)
      else if(nbct.eq.2) then
         aa(1,1) = 0.d0
         aa(1,2) = dz(1)/dt * ( por*fncdsdp(pr(1),imatl(1),alpha,
     +    		   sres,expn,se,he) + stor(imatl(1))*sat(1) )
     +             + diff2
         aa(1,3) = -diff2
         rhs(1) = -dz(1)/dt * ( por*sat(1)-porold*satold(1) )
     +          + bct
     +          + diff2*(pr(2)-pr(1)-dz2)
     +          - trans(1)
     +          - fncevap(sat(1),imatl(1),t,swilt,sres,zevap,pi,
     +		evmax,etfin,etstart,tday)*evmult(1)*eweight(1)*dz(1)
      end if


!    2) Bottom block (i=nblocks):
      ne = nblocks
      por = fncpor(pr(ne),imatl(ne),poros,stor,prinit(ne))
      porold = fncpor(prold(ne),imatl(ne),poros,stor,prinit(ne))
!      por = poros(imatl(ne))
!      porold = por
      xksat1 = rksat(imatl(ne-1))*rksat(imatl(ne))*( dz(ne-1)+dz(ne) )
     +       / ( rksat(imatl(ne))*dz(ne-1)+rksat(imatl(ne-1))*dz(ne) )
      dz1 = (dz(ne-1)+dz(ne)) / 2.d0
        delh1 = (pr(ne)-pr(ne-1))/dz1 -1
      if (delh1.le.-1) then
             coefkr = 1.d0
        elseif (delh1.ge.1) then
             coefkr = 0.d0
        else
             coefkr = 0.5 - delh1/2.d0
        end if 
        xkrel1 = fnckrel(sat(ne-1),imatl(ne-1),expn,sres)*coefkr + 
     +            fnckrel(sat(ne),imatl(ne),expn,sres)*(1.0-coefkr)
!      xkrel1 = ( fnckrel(sat(ne-1),imatl(ne-1))*dz(ne)
!     +           + fnckrel(sat(ne),imatl(ne))*dz(ne-1) )
!     +         / ( dz(ne-1)+dz(ne) )   
      diff1 = xksat1*xkrel1/dz1

      if(nbcb.eq.1) then
         pr0 = 2.d0*bcb-pr(ne)
         sat0 = fncsat(bcb,imatl(ne),expn,sres,alpha,he,se)
         xksat2 = rksat(imatl(ne))
         dz2 = dz(ne)
         delh2 = (pr0-pr(ne))/dz2 -1
           if (delh2.le.-1) then
             coefkr = 1.d0
           elseif (delh2.ge.1) then
             coefkr = 0.d0
           else
             coefkr = 0.5 - delh2/2.d0
           end if 
           xkrel2 = fnckrel(sat(ne),imatl(ne),expn,sres)*coefkr + 
     +         fnckrel(sat0,imatl(ne),expn,sres)*(1.0-coefkr)
!         xkrel2 = ( fnckrel(sat(ne),imatl(ne))
!     +              + fnckrel(sat0,imatl(ne)) ) / 2.d0
         diff2 = xksat2*xkrel2/dz2
         aa(ne,1) = -diff1
!         ftemp = fncdsdp(pr(ne),imatl(ne))
         aa(ne,2) = dz(ne)/dt * ( por*fncdsdp(pr(ne),imatl(ne),alpha,
     +    		    sres,expn,se,he) + stor(imatl(ne))*sat(ne) )
     +             + diff1 + 2.d0*diff2
         aa(ne,3) = 0.d0
         rhs(ne) = -dz(ne)/dt * ( por*sat(ne)-porold*satold(ne) )
     +      - diff1*(pr(ne)-pr(ne-1)-dz1)
     +      + diff2*(pr0-pr(ne)-dz2)
     +      - trans(ne)
     +      - fncevap(sat(ne),imatl(ne),t,swilt,sres,zevap,pi,
     +		evmax,etfin,etstart,tday)*evmult(ne)*eweight(ne)*dz(ne)
      else if(nbcb.eq.2) then
         aa(ne,1) = -diff1
         aa(ne,2) = dz(ne)/dt * ( por*fncdsdp(pr(ne),imatl(ne),alpha,
     +    		    sres,expn,se,he)+ stor(imatl(ne))*sat(ne) )
     +             + diff1
         aa(ne,3) = 0.d0
         rhs(ne) = -dz(ne)/dt * ( por*sat(ne)-porold*satold(ne) )
     +      - diff1*(pr(ne)-pr(ne-1)-dz1)
     +      - bcb
     +      - trans(ne)
     +      - fncevap(sat(ne),imatl(ne),t,swilt,sres,zevap,pi,
     +		evmax,etfin,etstart,tday)*evmult(ne)*eweight(ne)*dz(ne)
      end if

      return
      end
