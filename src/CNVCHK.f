      subroutine cnvchk (errpr,errres,rmax,princr,nblocks,
     + 			pr,sat,imatl,ntest,expn,sres,alpha,he,se)
        implicit none
!***********************************************************************
!     Routine:  cnvchk.for                                             *
!     Function: Check iteration for convergence.                       *
!     Written:  May 1992, January 2001, Summer 2005                    *
!     By:       M. Celia, M.J. Puma                                    *
!***********************************************************************

	  integer*4, parameter :: prec10 = selected_real_kind(10)
      integer*4, parameter :: prec10int = selected_int_kind(10)
      integer*4, parameter :: nbdim  = 320
      integer*4, parameter :: ntdim  = 14020
      integer*4, parameter ::  matdim =  5


      real(kind = prec10) errpr, errres, rmax    
      real(kind = prec10) pr(nbdim),princr(nbdim),sat(nbdim),se(matdim),
     +		expn(matdim),sres(matdim),alpha(matdim),he(matdim)
      integer*4  imatl(nbdim),nblocks    

      real(kind = prec10) fncsat,prmax
      integer*4 ntest,i,nres,npr
	  integer(kind = prec10int) ptest
      
      npr = 0
      nres = 0
      prmax = 0.d0
      do i=1,nblocks
         pr(i) = pr(i)+princr(i)
         sat(i) = fncsat(pr(i),imatl(i),expn,sres,alpha,he,se)
         ptest = dabs(princr(i))
         if(ptest.gt.prmax) then
            prmax = ptest
         end if
      end do

      if(prmax.gt.errpr) then
         npr = 1
      end if
      if(rmax.gt.errres) then
         nres = 1
      end if
      
      ntest = npr+nres

      return
      end
