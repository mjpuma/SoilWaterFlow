      subroutine geometry (flname,ztop,zbot,zevap,zroot,dz, 
     +               evmult,etmult,eweight,rweight,nblocks,nevap, 
     +               nroot,zbbot,zbtop)
        implicit none
!********************************************************************
!   Routine:  geometry.for
!   Purpose:  To read goemetric (grid and root zone) information
!              for cell-centered unsaturated flow code.
!   Written:  December 2000, February 2003, Summer 2005
!   By:       M. Celia, M.J. Puma                                   
!********************************************************************

	  integer*4, parameter :: prec10 = selected_real_kind(10)
      integer*4, parameter :: nbdim  = 320
      integer*4, parameter :: nutemp = 15      
      
	  character*12 flname

      real(kind = prec10) ztop, zbot, zevap, zroot, dz(nbdim), 
     +               evmult(nbdim), etmult(nbdim)
      real(kind = prec10) eweight(nbdim),rweight(nbdim)
      integer*4  nblocks, nevap, nroot 
      
	  real(kind = prec10) zbbot, zbtop,dztemp
	  integer*4 i,j,nfinish,nstart,nnotdz,dzconst,neweight,
     +         nrweight,ngrid

!   Open file
      open(unit=nutemp,file=flname)

!   Read Number of FD Blocks, Location of top boundary, Location of
!   bottom of evaporation and root zones:
      read(nutemp,*) nblocks, ztop, zevap, zroot

!   Read Block Sizes:
      read(nutemp,*) ngrid,nrweight,neweight
      if(ngrid.eq.1) then
         read(nutemp,*) dzconst, nnotdz
         do i=1,nblocks
            dz(i) = dzconst
         end do
         if(nnotdz.gt.0) then
            do i=1,nnotdz
               read(nutemp,*) nstart, nfinish, dztemp
               do j=nstart, nfinish
                  dz(j) = dztemp
               end do
            end do
         end if
      else
         read(nutemp,*) (dz(i), i=1,nblocks)
      end if

!   Determine location of bottom domain boundary:
      zbot = ztop
      do i=1,nblocks
         zbot = zbot + dz(i)
         rweight(i) = 0.d0
         eweight(i) = 0.d0
      end do

!   Determine blocks in which ET can occur:
      zbtop = ztop
      nroot = 0
      do i=1,nblocks
         zbbot = zbtop + dz(i)
         if(zroot.gt.zbbot) then
            etmult(i) = 1.d0
            nroot = i
         else if (zroot.le.zbtop) then
            etmult(i) = 0.d0
         else
            etmult(i) = (zroot-zbtop)/dz(i)
            nroot = i
         end if
         zbtop = zbbot
      end do

      zbtop = ztop
      do i=1,nblocks
         zbbot = zbtop + dz(i)
         if(zevap.gt.zbbot) then
            evmult(i) = 1.d0
            nevap = i
         else if (zevap.le.zbtop) then
            evmult(i) = 0.d0
         else
            evmult(i) = (zevap-zbtop)/dz(i)
            nevap = i
         end if
         zbtop = zbbot
      end do

!  Assign transpiration (root) weights as a function of depth
        if (nrweight.eq.1) then
           do i=1,nroot
              rweight(i) = 1.d0
           end do
        else
           do i=1,nroot
             read(nutemp,*) rweight(i)
           end do
        end if

!  Assign  evaporation weights as a function of depth
        if (neweight.eq.1) then
           do i=1,nevap
              eweight(i) = 1.d0
           end do
        else
           do i=1,nevap
              read(nutemp,*) eweight(i)
           end do
        end if
      
!  This ends the geometry input
        
        close(unit=nutemp)
        
      return
      end