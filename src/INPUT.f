      subroutine input(ztop,zbot,zevap,zroot,dz, 
     +               evmult,etmult,eweight,rweight,zbbot,zbtop,
     + 				 tstart,dt,tmax,tday,etstart,etfin,
     +               dtredu,dtincr,dtmin,dtmax,errpr,errres,
     +               pr,sat,poros,alpha,expn,rksat,stor,
     +               etmax,evmax,frac,sres,hwilt,hstar,swilt,
     +				 rtavgdens,Se,he,bctop,timet,bcbot,
     +				 timeb,nblocks,nevap,
     +				 nroot,ipond,mprint,iscrn,imax,
     +               itmax,itmin,nmatl,imatl,prinit,
     +               dtprint,nbctop,nbcbot,ntmax)    
        implicit none
!*********************************************************************
!   Routine:  input.f                                              *
!   Purpose:  Read input variables for block-centered plant1d        *
!              flow code.                                            *
!   Written:  December 2000, February 2003, Summer 2005              *
!   By:       M. Celia, M.J. Puma                                    *
!*********************************************************************
!
	  integer*4, parameter :: prec10 = selected_real_kind(10)      
      integer*4, parameter :: nuin = 11 
      integer*4, parameter :: nbdim  = 320
      integer*4, parameter :: matdim =  5	
      integer*4, parameter :: ntdim  = 14020  
	  
	  character*12 flname      
      character*75 title
      
	  real(kind = prec10) ztop,zbot,zevap,zroot,
     +      zbbot,zbtop,dz(nbdim), 
     + 		tstart,dt,tmax,tday,etstart,etfin,
     +      dtredu,dtincr,dtmin,dtmax,errpr,errres,
     +      pr(nbdim),sat(nbdim),poros(matdim),alpha(matdim),
     +      expn(matdim), rksat(matdim), stor(matdim),
     +      etmax(matdim), evmax(matdim),frac, 
     +      sres(matdim), hwilt(matdim), hstar(matdim),
     +      swilt(matdim),
     +		rtavgdens(matdim),Se(matdim), he(matdim), bctop(ntdim),
     +      timet(ntdim), bcbot(ntdim),timeb(ntdim),    
     +      eweight(nbdim),rweight(nbdim),
     + 		evmult(nbdim), etmult(nbdim),prinit(nbdim)
           
      integer*4 nblocks,nevap,nroot,ipond,mprint,iscrn,imax,
     +          itmax,itmin,nmatl,imatl(nbdim),ntmax
     
      real(kind = prec10) dtprint
      integer*4 nbctop(ntdim), nbcbot(ntdim)


!   Open Main Input File
      open(unit=nuin, file='plant1d.in')

!   Read Title for this Simulation
      read(nuin,10) title

!   Read Name of File for Geometry Input
      read(nuin,10) flname

!   Geometry Input
      call geometry(flname,ztop,zbot,zevap,zroot,dz, 
     +               evmult,etmult,eweight,rweight,nblocks,nevap, 
     +               nroot,zbbot,zbtop)

!   Read Name of File for Control Parameters
      read(nuin,10) flname

!   Control Parameter Input
      call control(flname,ipond,tstart,dt,tmax,tday,etstart, 
     +             etfin,mprint,iscrn,imax,itmax,itmin,dtprint,
     +             dtredu,dtincr,dtmin,dtmax,errpr,errres,ntmax)

!   Read Name of File for Material Properties
      read(nuin,10) flname

!   Materials Properties Input
      call propmat (flname,nblocks,poros,alpha,
     +                expn,rksat,stor,etmax,evmax,frac,sres,
     +                hwilt,hstar,swilt,rtavgdens,Se,he,nmatl,imatl)

!   Read Name of File for Initial Condition Input
      read(nuin,10) flname

!   Initial Condition Input
      call initial(flname,nblocks,imatl,ztop,dz,pr,sat,
     +		expn,sres,alpha,he,se,prinit)

!   Read Name of File for Boundary Condition Input
      read(nuin,10) flname

!   Boundary Condition Input
      call bcinput(flname,tmax,imatl,bctop,timet,bcbot,timeb,
     +            expn,sres,alpha,he,se,nbctop,nbcbot)

!   End of Input
      close(unit=nuin)

   10 format (a)

      return
      end
