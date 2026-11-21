      subroutine control(flname,ipond,tstart,dt,tmax,tday,etstart, 
     +             etfin,mprint,iscrn,imax,itmax,itmin,dtprint,
     +             dtredu,dtincr,dtmin,dtmax,errpr,errres,ntmax)
        implicit none
!********************************************************************
!   Routine:  control.for                                           *
!   Purpose:  To read control parameters for the unsat flow code.   *
!   Written:  December 2000, Summer 2005                            *
!   By:       M. Celia, M.J. Puma                                   *
!********************************************************************
 
	  integer*4, parameter :: prec10 = selected_real_kind(10)
      integer*4, parameter :: nbdim  = 320
      integer*4, parameter :: nutemp = 15  

      character*12 flname

      integer*4 ipond
      real(kind = prec10) tstart, dt, tmax
      real(kind = prec10) tday, etstart, etfin
      integer*4  mprint, iscrn
      integer*4  imax, itmax, itmin,ntmax 
      real(kind = prec10) dtredu, dtincr, dtmin, dtmax
      real(kind = prec10) errpr, errres
		
      real(kind = prec10) dtprint     
  

!   Open file
      open(unit=nutemp,file=flname)

!   Read Starting Time
      read(nutemp,*) tstart

!   Read Initial Time Step Size
      read(nutemp,*) dt

!   Read Max. simulation time
      read(nutemp,*) tmax

!   Read Max. Number of Time Steps
      read(nutemp,*) ntmax

!   Read Length of Time corresponding to One Day (for ET)
      read(nutemp,*) tday

!   Read Time of Day to Start ET
      read(nutemp,*) etstart

!   Read Time of Day to End ET
      read(nutemp,*) etfin

!   Read Integer Flag for Ponding Option
      read(nutemp,*) ipond

!   Read angle 1-D domain makes with the vertical direction
!     read(nutemp,*) phi

!   Read Print Interval for vertically-integrated saturation
      read(nutemp,*) dtprint

!   Read multiple for printing saturation profile
      read(nutemp,*) mprint

!   Read flag to control printing to the screen
      read(nutemp,*) iscrn
   
!   Read Max. Number of Nonlinear Iterations
      read(nutemp,*) imax

!   Read min and max number of interations before a change in time step
      read(nutemp,*) itmin
      read(nutemp,*) itmax

!   Read time step multipliers for decreasing or increasing the time step
      read(nutemp,*) dtredu
      read(nutemp,*) dtincr

!   Read minimum and maximum time step sizes
      read(nutemp,*) dtmin
      read(nutemp,*) dtmax

!   Read Error Tolerances for Nonlinear Convergence
      read(nutemp,*) errpr
      read(nutemp,*) errres
 
!   This ends control parameter input.

      close(unit=nutemp)

      return
      end

