       subroutine bcchk(bctop,timet,bcbot,timeb,nbctop,nbcbot,t,dtmin,
     +      iqflag,ibcflag,nprtop,nprbot,bct,bcb,nbct,nbcb,
     +		nprtad,nprbad)
        implicit none
!***********************************************************************
!     Routine:  bcchk                                                  *
!     Function: Check boundary condition counters and update as        *
!               appropriate.                                           *
!     Written:  January 2001, Summer 2005                              *
!     By:       M. Celia, M.J. Puma                                    *
!***********************************************************************
!
	  integer*4, parameter :: prec10 = selected_real_kind(10) 

      integer*4, parameter :: nbdim  = 320
      integer*4, parameter :: ntdim  = 14020
      integer*4, parameter ::  matdim =  5
    
      real(kind = prec10) bctop(ntdim), timet(ntdim), bcbot(ntdim),
     +              timeb(ntdim), bct, bcb
      integer*4 nbctop(ntdim), nbcbot(ntdim), 
     +               nprtop, nprbot, nbct, nbcb, ibcflag

      real(kind = prec10) t, dtmin, tcheck
      integer*4  nprbad, nprtad, iqflag

      nprtad = 0
      nprbad = 0

      tcheck = abs(t - timet(nprtop))
      
      if((tcheck.lt.dtmin).and.(ibcflag.eq.0))then
        
         nprtad = 1
         nprtop = nprtop + nprtad
         ibcflag = 1
         
         if(iqflag.eq.0) then
            nbct = nbctop(nprtop-ibcflag)
            bct = bctop(nprtop-ibcflag)
         end if
         
         if(nbctop(nprtop).eq.1) then
            iqflag = 0
            nbct = nbctop(nprtop-ibcflag)
            bct = bctop(nprtop-ibcflag)
         end if
        
         
           
      elseif((ibcflag.eq.1).and.(tcheck.gt.dtmin))then
         
         nprtad = 0
         ibcflag = 0
       
         if(iqflag.eq.0) then
            nbct = nbctop(nprtop)
            bct = bctop(nprtop)
         end if
         
         if(nbctop(nprtop).eq.1) then
            iqflag = 0
            nbct = nbctop(nprtop)
            bct = bctop(nprtop)
         end if
          
      
      elseif((ibcflag.eq.1).and.(tcheck.lt.dtmin))then   
         
         ibcflag = 1
         nprtad = 1
         nprtop = nprtop + nprtad
         
         if(iqflag.eq.0) then
            nbct = nbctop(nprtop-ibcflag)
            bct = bctop(nprtop-ibcflag)
         end if
         
         if(nbctop(nprtop).eq.1) then
            iqflag = 0
            nbct = nbctop(nprtop-ibcflag)
            bct = bctop(nprtop-ibcflag)
         end if
      end if        

!     CHANGED nprbad from 1 to 0       
      if(t.gt.timeb(nprbot)) then
         nprbad = 0          
         nprbot = nprbot + nprbad
         nbcb = nbcbot(nprbot)
         bcb = bcbot(nprbot)
      end if
!
      return
      end
