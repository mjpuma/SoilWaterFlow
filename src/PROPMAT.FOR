      subroutine propmat (flname,nblocks,poros,alpha,
     +                expn,rksat,stor,etmax,evmax,frac,sres,
     +                hwilt,hstar,swilt,rtavgdens,Se,he,nmatl,imatl)
        implicit none
!********************************************************************
!   Routine:  propmat.for                                           *
!   Purpose:  To read material (plant and soil) properties for      *
!             the Richards model                                    *
!   Written:  December 2000, October 2003, Summer 2005              *
!   By:       M. Celia, M.J.Puma                                    *
!********************************************************************

	  integer*4, parameter :: prec10 = selected_real_kind(10)

      character*12 flname

      integer*4, parameter :: nbdim  = 320
      integer*4, parameter :: ntdim  = 14020
      integer*4, parameter ::  matdim =  5

      integer*4, parameter :: nuout  = 12
      integer*4, parameter :: nutemp = 15
      
      real(kind = prec10) poros(matdim), alpha(matdim),
     +                expn(matdim), rksat(matdim), stor(matdim)
      real(kind = prec10) etmax(matdim), evmax(matdim),frac, 
     +                sres(matdim), hwilt(matdim), hstar(matdim),
     +                swilt(matdim)                 
                        
      real(kind = prec10)  rtavgdens(matdim)
      real(kind = prec10) Se(matdim), he(matdim)
      integer*4  nmatl, imatl(nbdim),i,j,nblocks

	  integer*4 mattemp, nfinish,nstart,nnotmat,matconst,nparam
      real(kind = prec10) fncpr,fncsat
      
!   Open file
      open(unit=nutemp,file=flname)

!   Material (Soil and Plant) Properties 
      read(nutemp,*) nmatl
      do i=1,nmatl
         read(nutemp,*) j,poros(j),alpha(j),expn(j),rksat(j),
     +                  stor(j),sres(j),Se(j)
         read(nutemp,*) frac,rtavgdens(j),etmax(j),evmax(j),
     +                  hwilt(j),hstar(j)
         if(j.ne.i) then
            write(nuout,10) i,j
   10       format(1x,70(1h*)/5x,'WARNING FROM INPUT READ:'/
     +          10X,'CHECK INPUT FOR PROPERTIES OF MATERIALS NUMBERED',
     +          I3,' AND',I3/1X,70(1H*))
         end if
!   Heads for wilting and stomatal closure are stored as negative values
           if (hwilt(j).gt.0) then
            hwilt(j) = -hwilt(j)
         end if
         if (hstar(j).gt.0) then
            hstar(j) = -hstar(j)
         end if
         he(j) = 0.d0
         he(j) = fncpr(Se(j),j,expn,sres,alpha,he,se)
         he(j) = -alpha(j)
         swilt(j) = fncsat(hwilt(j),j,expn,sres,alpha,he,se)
      end do

!   Read Material Properties for each Block:
      read(nutemp,*) nparam
      if(nparam.eq.1) then
         read(nutemp,*) matconst, nnotmat
         do i=1,nblocks
            imatl(i) = matconst
         end do
         if(nnotmat.gt.0) then
            do i=1,nnotmat
               read(nutemp,*) nstart, nfinish, mattemp
               do j=nstart, nfinish
                  imatl(j) = mattemp
               end do
            end do
         end if
      else
         read(nutemp,*) (imatl(i), i=1,nblocks)
      end if

!   This ends material properties input.
      close(unit=nutemp)

      return
      end

