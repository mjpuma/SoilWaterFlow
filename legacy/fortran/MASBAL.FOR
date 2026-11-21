      subroutine masbal(dt,trans,dz,evmult,etmult,eweight,nblocks,nroot,
     +		nbcb,nbct,rksat,poros,stor,prinit,expn,sres,imatl,
     +      alpha,he,se,zevap,pi,pr,prold,sat,sumtot,sumtotr,
     +	 	evmax,etfin,etstart,tday,bct,bcb,swilt,t,
     + 		ratio,ratio1,qin,qout,qleak,fluxin,flxout,
     +      fluxet,fluxev,transum,evsum)
        implicit none
!***********************************************************************
!     Routine:  masbal.for                                             *
!     Function: Compute mass balance for given time step.              *
!     Written:  July 1990, May  1992                                   *
!               January 2001, February 2003, Summer 2005               *
!     By:       M. Celia, A. Guswa, M.J.Puma                           *
!***********************************************************************

	  integer*4, parameter :: prec10 = selected_real_kind(10)
      
!      include 'dimen.par'
      integer*4, parameter :: nbdim  = 320
      integer*4, parameter :: ntdim  = 14020
      integer*4, parameter ::  matdim =  5

!      include 'units.par'
      integer*4, parameter :: nuscrn =  1
      integer*4, parameter :: nuin   = 11
      integer*4, parameter :: nuout  = 12
      integer*4, parameter :: nusave = 13
      integer*4, parameter :: nutemp = 15
      integer*4, parameter :: numb   = 16
      integer*4, parameter :: nusat  = 18
      integer*4, parameter :: nusatt = 19
      integer*4, parameter :: nucon  = 21
      integer*4, parameter :: nupl1  = 22
      integer*4, parameter :: nupl2  = 23
      integer*4, parameter :: nuic   = 24
      integer*4, parameter :: nuanm  = 25
      integer*4, parameter :: nupl3  = 26
      integer*4, parameter :: nupl4 =  27
      integer*4, parameter :: nulow =  28
      integer*4, parameter :: numid =  29
      integer*4, parameter :: nuhi =  30 
      integer*4, parameter :: nuvhi =  31 
      integer*4, parameter :: nupl5 =  36       

   
      real(kind = prec10) bct,bcb
      real(kind = prec10) dt,t
      real(kind = prec10)  fluxin, flxout, fluxet, fluxev, qin,
     +               sumtot, sumtotr,  ratio, ratio1, 
     +               transum, evsum, qleak 
      real(kind = prec10) pr(nbdim), prold(nbdim),  
     +               sat(nbdim), trans(nbdim)    

      real(kind = prec10) dz(nbdim),evmult(nbdim), etmult(nbdim)
      real(kind = prec10) eweight(nbdim)
      integer*4  nblocks, nroot,imatl(nbdim)    

      real(kind = prec10) denom2,denom1,fncevap,delh2,dz2,qout,
     +		fnckrel,xkrel,coefkr,delh1,dz1,xksat,sat0,
     +		pr0,addw,porold,fncpor,por,sold,fncsat,snew,
     +      sumroot,sumfull,rksat(matdim) 
      integer*4 ne,ii,i,nbcb,nbct
      
  	  real(kind = prec10) poros(matdim),stor(matdim)      
      real(kind = prec10) prinit(nbdim),expn(matdim),sres(matdim),
     +      alpha(matdim),he(matdim),se(matdim),swilt(matdim),
     +      zevap,pi,
     +	 	evmax(matdim),etfin,etstart,tday 

!   Calculate change in total mass
      sumfull = 0.d0
      sumroot = 0.d0
      do i=1,nblocks
         ii = imatl(i)
         snew = fncsat(pr(i),ii,expn,sres,alpha,he,se)
         sold = fncsat(prold(i),ii,expn,sres,alpha,he,se)
         por = fncpor(pr(i),ii,poros,stor,prinit(i))
         porold = fncpor(prold(i),ii,poros,stor,prinit(i))
!!!         st = stor(ii)*snew*(pr(i)-prold(i))
         addw = por*snew-porold*sold 
         sumfull = sumfull + addw*dz(i)
         sumroot = sumroot + addw*dz(i)*etmult(i)
      end do
!
      sumtot = sumtot+sumfull
      sumtotr = sumtotr+sumroot
!
!   Volumetric Flux at Top Boundary
      if(nbct.eq.2) then
         qin = bct
      elseif (nbct.eq.1) then
         pr0 = 2.d0*bct-pr(1)
         sat0 = fncsat(bct,imatl(1),expn,sres,alpha,he,se)
         xksat = rksat(imatl(1))
         dz1 = dz(1)
         delh1 = (pr(1)-pr0)/dz1 -1
         
         if (delh1.le.-1) then
             coefkr = 1
         elseif (delh1.ge.1) then
             coefkr = 0
         else
             coefkr = 0.5 - delh1/2.d0
         end if 
         
         xkrel = fnckrel(sat0,imatl(1),expn,sres)*coefkr + 
     +            fnckrel(sat(1),imatl(1),expn,sres)*(1-coefkr)
         qin = xksat*xkrel * ( (pr0-pr(1))/dz(1)+1.d0 )
      end if
!
!   Volumetrix Flux at Bottom Boundary
      ne = nblocks
      if(nbcb.eq.2) then
         qout = bcb
      elseif(nbcb.eq.1) then
         pr0 = 2.d0*bcb-pr(ne)
         sat0 = fncsat(bcb,imatl(ne),expn,sres,alpha,he,se)
         xksat = rksat(imatl(ne))
         dz2 = dz(ne)
         delh2 = (pr0-pr(ne))/dz2 -1
         
         if (delh2.le.-1) then
             coefkr = 1
         elseif (delh2.ge.1) then
             coefkr = 0
         else
             coefkr = 0.5 - delh2/2.d0
         end if 
         
         xkrel = fnckrel(sat(ne),imatl(ne),expn,sres)*coefkr + 
     +         fnckrel(sat0,imatl(ne),expn,sres)*(1-coefkr)
         qout = xksat*xkrel*((pr(ne)-pr0)/dz(ne) + 1.d0)
      end if
!
!   Losses due to transpiration
      transum = 0.d0
      do i=1,nblocks
        transum = transum + dt*trans(i)
      end do
!
!   Losses due to evaporation
      evsum = 0.d0
      do i=1,nblocks
        evsum = evsum + 
     +        dt*fncevap(sat(i),imatl(i),t,swilt,sres,zevap,pi,
     +		evmax,etfin,etstart,tday)*evmult(i)*eweight(i)*dz(i)
      end do  
      
!   Losses due to leakage out of the root zone      
      qleak = 0.d0
      xksat = rksat(imatl(nroot))*rksat(imatl(nroot+1))*
     +        (dz(nroot)+dz(nroot+1)) / ( rksat(imatl(nroot+1))
     +        *dz(nroot)+rksat(imatl(nroot))*dz(nroot+1) )
      dz2 = (dz(nroot)+dz(nroot+1)) / 2.d0
      delh2 = (pr(nroot+1)-pr(nroot))/dz2 -1
      if (delh2.le.-1) then
         coefkr = 1
      elseif (delh2.ge.1) then
         coefkr = 0
      else
         coefkr = 0.5 - delh2/2.d0
      end if 
         
      xkrel = fnckrel(sat(nroot),imatl(nroot),expn,sres)*coefkr + 
     +      fnckrel(sat(nroot+1),imatl(nroot+1),expn,sres)*(1-coefkr)
      qleak = xksat*xkrel*((pr(nroot)-pr(nroot+1))/dz(nroot) + 1.d0)    
           
!   Fluxes over time step dt and Total fluxes
      qin = qin*dt
      qout = qout*dt
      qleak = qleak*dt
      fluxin = fluxin+qin
      flxout = flxout+qout
      fluxet = fluxet+transum
      fluxev = fluxev+evsum
      denom1 = fluxin-flxout-fluxet-fluxev
      if(dabs(denom1).ge.1.d-8) then
         ratio = sumtot/denom1
      else
         ratio = 0.d0
      end if
      denom2 = qin-qout-transum-evsum
      if(dabs(denom2).ge.1.d-8) then
         ratio1 = sumfull/denom2
      else
         ratio1 = 0.d0
      end if
!
      return
      end