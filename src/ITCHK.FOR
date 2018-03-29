      subroutine itchk(bctop,bcbot,nbctop,nbcbot,iprad,dtprint,imax,
     +		itmax,iter,dtredu,dtmin,nprbad,nprtad,rmax,
     +		prold,satold,nblocks,t,dt,tprint,
     +		bct,bcb,nbct,nbcb,ipr,ipr2,ipr3,pr,sat,itest,
     + 		itime,nprtop,nprbot)
        implicit none
!***********************************************************************
!     Routine:  itchk.for                                              *
!     Function: Check for exceedence of max. number of iterations.     *
!     Written:  May 1992, January 2001, February 2003, Summer 2005     *
!     By:       M. Celia, M.J. Puma                                    *
!***********************************************************************

	  integer*4, parameter :: prec10 = selected_real_kind(10)

      integer*4, parameter :: nbdim  = 320
      integer*4, parameter :: ntdim  = 14020
      integer*4, parameter ::  matdim =  5
      integer*4, parameter :: nuout  = 12
          
      real(kind = prec10) bctop(ntdim), bcbot(ntdim),
     +              bct, bcb
      integer*4 nbctop(ntdim), nbcbot(ntdim), 
     +               nprtop, nprbot, nbct, nbcb

      integer*4  itime, ipr, ipr2, ipr3, iprad
      real(kind = prec10)  dt, t
      real(kind = prec10) tprint, dtprint
      integer*4  imax, itmax, iter, itest
      real(kind = prec10) dtredu, dtmin
      integer*4  nprbad, nprtad
      real(kind = prec10) rmax     

      real(kind = prec10) pr(nbdim), prold(nbdim),  
     +               sat(nbdim), satold(nbdim)    
      integer*4  nblocks, i      


      if(iter.gt.itmax) then
         t = t-dt
         dt = dt*dtredu
         if((dt.lt.dtmin).or.(dtredu.ge.1.d0)) then
            write(nuout,10) imax,itime,dt,dtmin,dtredu,rmax
   10       format(//
     +         5X,16(1H*),2X,'MAX ITERATIONS EXCEEDED',16(1H*)/
     +         26X,I5,' ITERATIONS'/
     +         25X,'TIME STEP',I5/
     +         'dt = ',e14.3,5x,'dtmin = ',e14.3,5x,'dtredu = ',e14.3/
     +         'rmax = ',e14.3/)
            stop
         else
!   Reset pressures, saturations and print, iteration, and 
!     boundary condition counters

            ipr = 0
            ipr2 = ipr2 - iprad
            ipr3 = ipr3 - iprad 
            tprint = tprint - iprad*dtprint
            itime = itime - 1
            itest = 1
            
            do i=1,nblocks
               pr(i) = prold(i)
               sat(i) = satold(i)
            end do
            
            nprtop = nprtop - nprtad
            nprbot = nprbot - nprbad 
            
            nbct = nbctop(nprtop)
            bct = bctop(nprtop)
            nbcb = nbcbot(nprbot)
            bcb = bcbot(nprbot)
            
         end if
      end if


      return
      end
