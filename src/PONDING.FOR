      subroutine ponding(bctop,nbctop,nprtop,ipond,iscrn,dt,t,qin,pr,
     +		iqflag,bct,nbct,qbuild,qinflx) 
       implicit none
!***********************************************************************
!     Routine:  ponding.for                                            *
!     Function: Calculate appropriate parameters for ponding condition.*
!     Written:  July 1990, May  1992, January 2001, Summer 2005        *
!     By:       M. Celia, M.J. Puma                                    *
!***********************************************************************
	  integer*4, parameter :: prec10 = selected_real_kind(10)

      integer*4, parameter :: nbdim  = 320
      integer*4, parameter :: ntdim  = 14020
      integer*4, parameter :: matdim =  5
      integer*4, parameter :: nuout  = 12
      integer*4, parameter :: nuanm  = 25

   
      real(kind = prec10) bctop(ntdim), bct
      integer*4 nbctop(ntdim), nprtop, nbct
      integer*4 ipond,iscrn

      real(kind = prec10) dt,t,pr(nbdim)

      real(kind = prec10) qin,qbuild,qinflx
	  integer*4 iqflag

      if(iqflag.eq.0) then
         write(nuout,5) t
    5    format(//5x,'PONDING HAS BEGUN AT TIME T =',e15.6)
         if(iqflag.eq.0.and.iscrn.eq.1) then
            write(*,*) 'Ponding Underway'
         end if
         iqflag = 1
         nbct = 1
         bct = pr(1) 
!
!   Adjust Boundary Conditions if ponding is underway
!      if(iqflag.eq.1) then
        else
         qbuild = qbuild+bctop(nprtop)*dt-qin
         qinflx = qin/dt
!        if(ipr.eq.iprint) 
           write(2,10) t,qbuild,qinflx,bctop(nprtop)
   10    format(/5x,'PONDING INFORMATION (TIME =',e11.3,')'/
     +          10X,'BUILD-UP OF HEAD AT SURFACE',E19.3/
     +          10X,'CURRENT FLUX INTO DOMAIN AT SURFACE',E11.3/
     +          10X,'ORIGINAL FLUX SPECIFICATION',E19.3//) 
         
         write(nuanm,200) t,qin,dt,qbuild,qinflx,bctop(nprtop) 
  200    format(6f20.10)
     
         if(qbuild.gt.0.d0) then
            nbct = 1
            if(ipond.eq.1) then
               qbuild = 0.d0
            else if(ipond.eq.2) then
               bct = qbuild
            end if
         else
            nbct = nbctop(nprtop)
            bct = bctop(nprtop)
            iqflag = 0
            qbuild = 0.d0
         end if
      end if
!
!
      return
      end