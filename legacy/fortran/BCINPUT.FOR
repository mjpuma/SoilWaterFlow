      subroutine bcinput(flname,tmax,imatl,bctop,timet,bcbot,timeb,
     +            expn,sres,alpha,he,se,nbctop,nbcbot)
        implicit none
!********************************************************************
!   Routine:  bcinput.for
!   Purpose:  To read boundary conditions for FD Unsat Code.
!   Written:  January 2001, Summer 2005
!   By:       M. Celia, M. J. Puma
!********************************************************************

	  integer*4, parameter :: prec10 = selected_real_kind(10)

      character*12 flname

      integer*4, parameter :: nbdim  = 320
      integer*4, parameter :: ntdim  = 14020
      integer*4, parameter ::  matdim =  5


      integer*4, parameter :: nuout  = 12
      integer*4, parameter :: nutemp = 15
        
      real(kind = prec10) bctop(ntdim), timet(ntdim), bcbot(ntdim),
     +              timeb(ntdim),tmax,fncpr,expn(matdim),sres(matdim),
     +       alpha(matdim),he(matdim),se(matdim)
      integer*4 nbctop(ntdim), nbcbot(ntdim), ntopbc, nbotbc
      integer*4  imatl(nbdim),isat,i

 
!   Open file
      open(unit=nutemp,file=flname)
!
!   Read Number of Time Intervals used to define time-varying BC's
      read(nutemp,*) ntopbc,nbotbc
!
!   Top
      if(ntopbc.gt.ntdim) then
         write(nuout,70) ntopbc
   70    format(//1x,68(1h*)/5x,'WARNING FROM BOUNDARY COND INPUT:'/
     +           10X,'NUMBER OF B.C.CHANGES TOO LARGE'/
     +           10X,'MAX IS 20.  INPUT VALUE IS',I4,'.' /
     +           10X,'CHECK YOUR INPUT.  PROGRAM IS TERMINATING.',
     +            1x,68(1h*)//)
          stop
      else
         do i=1,ntopbc
            read(nutemp,*) nbctop(i),bctop(i),isat,timet(i)
            if(isat.eq.1.and.nbctop(i).eq.1) then
               bctop(i) = fncpr(bctop(i),imatl(i),expn,
     +                  sres,alpha,he,se)
            end if
         end do
         timet(ntopbc) = tmax
      end if
!
!     Bottom
      if(nbotbc.gt.ntdim) then
         write(nuout,70) nbotbc
         stop
      else
         do i=1,nbotbc
            read(nutemp,*) nbcbot(i),bcbot(i),isat,timeb(i)
            if(isat.eq.1.and.nbcbot(i).eq.1) then
               bcbot(i) = fncpr(bcbot(i),imatl(i),expn,
     +                    sres,alpha,he,se)
            end if
         end do
         timeb(nbotbc) = tmax
      end if
!
!
      return
      end

