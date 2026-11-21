      subroutine thomas(a,d,u,n,ndim)
!**********************************************************************
!     Thomas algorithm for the solution of tri-diagonal matrices.
!     Matrix a is assumed to be in banded form, with size nx3.
!     The first column of a contains entries in the left diagonal,
!     the second column contains entries from the main diagonal,
!     and the third column contains entries from the right diagonal.
!     The vector d is the right side vector.
!     The vector u is the unknown to be solved for.
!     n is the actual number of equations.
!     ndim is the dimension of the arrays.
!     n must be not greater than ndim.
!**********************************************************************
      implicit real*8(a-h,o-z)

      dimension a(ndim,3),d(ndim),u(ndim)

!   Define new entries for matrix elements by Gauss elimination
      a(1,1) = 0.d0
      a(1,3) = a(1,3)/a(1,2)
      d(1) = d(1)/a(1,2)
      a(1,2) = 1.d0
      if(n.eq.1) goto 15
      do 10 i=2,n
      q = a(i,2)-a(i,1)*a(i-1,3)
      d(i) = (d(i)-a(i,1)*d(i-1))/q
      a(i,1) = 0.d0
      a(i,2) = 1.d0
      a(i,3) = a(i,3)/q
   10 continue

!     solve
   15 u(n) = d(n)
      if(n.eq.1) return
      nn1 = n-1
      do 20 i=1,nn1
      jj = n-i
      u(jj) = d(jj)-a(jj,3)*u(jj+1)
   20 continue
!
      return
      end