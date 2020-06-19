      subroutine initial(flname,nblocks,imatl,ztop,dz,pr,sat,
     +		expn,sres,alpha,he,se,prinit)
        implicit none
!********************************************************************
!   Routine: initial.for
!   Purpose: To read initial conditions for FD unsat code
!   Written: January 2001, Summer 2005
!   By:      M. Celia, M.J. Puma
!********************************************************************

	  integer*4, parameter :: prec10 = selected_real_kind(10)

      character*12 flname

      integer*4, parameter :: nbdim  = 320
      integer*4, parameter :: nuout  = 12
      integer*4, parameter :: nutemp = 15
      integer*4, parameter ::  matdim =  5
      
 
      real(kind = prec10) pr(nbdim), sat(nbdim),ztop,dz(nbdim),zz,prref
      real(kind = prec10) sattemp,fncpr,prtemp,fncsat,value,zref 
      integer*4  imatl(nbdim),itemp,i,ndev,isat,iccode
      real(kind = prec10) prinit(nbdim),expn(matdim),sres(matdim),
     +       alpha(matdim),he(matdim),se(matdim)

      integer*4  nblocks 


!   Open Input File
      open (unit=nutemp,file=flname)

!   Read code for input structure
      read(nutemp,*) iccode, isat

!   Constant initial condition for either saturation
!       (isat = 1) or pressure head (isat (not =) 1)
      if(iccode.eq.1) then
         read(nutemp,*) ndev,value
         if(isat.ne.1) then
            do i=1,nblocks
               prinit(i) = value
               pr(i) = prinit(i)
               sat(i) = fncsat(value,imatl(i),expn,sres,alpha,
     +					he,se)
            end do
            if(ndev.gt.0) then
               do i=1,ndev
                  read(nutemp,*) itemp,prtemp
                  prinit(itemp) = prtemp
                  pr(itemp) = prinit(itemp)
                  sat(itemp) = fncsat(prtemp,imatl(itemp),
     +               expn,sres,alpha,he,se)
               end do
            end if
         else
            do i=1,nblocks
               sat(i) = value
               prinit(i) = fncpr(sat(i),imatl(i),expn,sres,alpha,he,se)
               pr(i) = prinit(i)
            end do
            if(ndev.gt.0) then
               do i=1,ndev
                  read(nutemp,*) itemp,sattemp
                  sat(itemp) = sattemp
                  prinit(itemp) = fncpr(sattemp,imatl(itemp),
     +             				expn,sres,alpha,he,se)
                  pr(itemp) = prinit(itemp)
               end do
            end if
         end if

!   Nonconstant I.C.:  read the initial values
      elseif(iccode.eq.2) then
         if(isat.eq.1) then
            read(nutemp,*) (sat(i), i=1,nblocks)
            do i=1,nblocks
               prinit(i) = fncpr(sat(i),imatl(i),expn,sres,alpha,he,se)
               pr(i) = prinit(i)
            end do
         else
            read(nutemp,*) (prinit(i), i=1,nblocks)
            do i=1,nblocks
               pr(i) = prinit(i)
               sat(i) = fncsat(pr(i),imatl(i),expn,sres,alpha,he,se)
            end do
         end if

!   Hydrostatic Initial Pressure
      elseif(iccode.eq.3) then
         read(nutemp,*) zref,prref
         zz = ztop
         do i=1,nblocks
            zz = zz + dz(i)/2.d0
            prinit(i) = prref+(zz-zref)
            pr(i) = prinit(i)
            sat(i) = fncsat(pr(i),imatl(i),expn,sres,alpha,he,se)
            zz = zz+dz(i)/2.d0
         end do
      else
         write(nuout,80)
   80    format(//68(1h*)/
     +    5x,'Improper Input for Initial Conditions; Check your input'
     +    /68(1h*)//)
      end if

      close(unit=nutemp)

      return
      end