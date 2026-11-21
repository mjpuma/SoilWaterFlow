      Program plant1d
        implicit none
!***********************************************************************
!     Routine:  plant1d.for                                            *
!     Function: Main routine for 1-D block-centered finite difference  *
!               code to solve unsaturated flow equation with soil-	   *
!               water loss due to transpiration and evaporation        *
!     Written:  Summer 2005          								   *
!     By:       M.J. Puma                        			  		   *
!***********************************************************************
    
	  integer*4, parameter :: prec10 = selected_real_kind(10)

!      include 'string.chr'       
!      character*75 title 
!      include 'string.inc'
            
!      include 'dimen.par'
      integer*4, parameter :: nbdim  = 320
      integer*4, parameter :: ntdim  = 14020
      integer*4, parameter ::  matdim =  5

!      include 'units.par'
      integer*4, parameter :: nuscrn =  1
      integer*4, parameter :: nuin   = 11
      integer*4, parameter :: nuout  = 12
      integer*4, parameter :: nusave = 13
      integer*4, parameter :: nutemp = 15
      integer*4, parameter :: numb   = 16
      integer*4, parameter :: nusat  = 18
      integer*4, parameter :: nusatt = 19
      integer*4, parameter :: nucon  = 21
      integer*4, parameter :: nupl1  = 22
      integer*4, parameter :: nupl2  = 23
      integer*4, parameter :: nuic   = 24
      integer*4, parameter :: nuanm  = 25
      integer*4, parameter :: nupl3  = 26
      integer*4, parameter :: nupl4 =  27
      integer*4, parameter :: nulow =  28
      integer*4, parameter :: numid =  29
      integer*4, parameter :: nuhi =  30 
      integer*4, parameter :: nuvhi =  31 
      integer*4, parameter :: nupl5 =  36       

!      include 'bcpr.inc'     
      real(kind = prec10) bctop(ntdim), timet(ntdim), bcbot(ntdim),
     +      timeb(ntdim), bct, bcb
      integer*4 nbctop(ntdim), nbcbot(ntdim),nprtop,nprbot, 
     +		nbct, nbcb, ibcflag

!	  include 'control.inc'
      integer*4 ipond
      integer*4  ntmax, itime, ipr, ipr2, ipr3, iprad
      real(kind = prec10) tstart, dt, dtold2, tmax, t
      real(kind = prec10) tday, etstart, etfin
      real(kind = prec10) tprint, dtprint
      integer*4  mprint, iscrn
      integer*4  imax, itmax, itmin, iter, itest
      real(kind = prec10) dtredu, dtincr, dtmin, dtmax
      integer*4  nprbad, nprtad
      real(kind = prec10) errpr, errres, rmax, errtrans      

!	  include 'masbal.inc'
      real(kind = prec10) tmass, fluxin, flxout, fluxet, fluxev, qin,
     +               sumtot, sumtotr, ratio, ratio1, 
     +               transum, evsum, qleak 

!	  include 'matrix.inc'
      real(kind = prec10) aa(nbdim,3), rhs(nbdim)
      real(kind = prec10) pr(nbdim), prold(nbdim), princr(nbdim), 
     +               sat(nbdim), satold(nbdim), trans(nbdim),
     +               daysatm(3,nbdim), satinstm(3,nbdim)     

!	  include 'param.inc'
      real(kind = prec10) pi, poros(matdim), alpha(matdim),
     +                expn(matdim), rksat(matdim), stor(matdim)
      real(kind = prec10) etmax(matdim), evmax(matdim),frac, 
     +                sres(matdim), hwilt(matdim), hstar(matdim),
     +                swilt(matdim)                 
                        
      real(kind = prec10) csoil, croot, rtavgdens(matdim)
      real(kind = prec10) Se(matdim), he(matdim)
      integer*4  nmatl, imatl(nbdim)

     
!	  include 'pond.inc'
      real(kind = prec10) qbuild
      integer*4  iqflag
     
!	  include 'space.inc'
      real(kind = prec10) ztop, zbot, zevap, zroot, dz(nbdim), 
     +               evmult(nbdim), etmult(nbdim)
      real(kind = prec10) eweight(nbdim),rweight(nbdim)
      integer*4  nblocks, nevap, nroot      


      real(kind = prec10) qinflx,qout,sumtrans,dtold,durvhi,
     +		durhi,durmid,durlow,static,dqin,dleak,devsum,
     +		dtransum,zbtop,zbbot
      integer*4 i,ntest,numvhi,numhi,nummid,numlow,
     +		ivhicros,ihicros,imidcros,ilowcros
      real(kind = prec10) prinit(nbdim)
      
	  real(kind = prec10) sumlow,summid,sumhi,sumvhi 

!   Open read file and files to save specific pieces of output
      open(nuin,file='plant1d.in')
      open(nuout,file='plant1d.out')
      open(numb,file='plant1d.mb')
      open(nupl1,file='plant1d.pl1')
      open(nupl2,file='plant1d.pl2')
      open(nupl3,file='plant1d.pl3')
      open(nupl4,file='plant1d.pl4')
      open(nupl5,file='plant1d.pl5')
      open(nuanm,file='plant1d.anm')
      open(nucon,file='plant1d.con')
      open(nuic,file='plant1d.ic')  
      open(nulow,file='plant1d.low')
      open(numid,file='plant1d.md')
      open(nuhi,file='plant1d.hi')
      open(nuvhi,file='plant1d.vhi')

!   Input data
      call input(ztop,zbot,zevap,zroot,dz, 
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
     
!   Initialize variables
      call initmvs(bctop,bcbot,nbctop,nbcbot,tstart,dt,dtprint,
     + 	bct,bcb,nprtop,nprbot,nbct,nbcb,ibcflag,itime,ipr,sat,rksat,
     + 	etmax,evmax,frac,hwilt,rtavgdens,zroot,dz,etmult,nblocks,
     +  poros,stor,prinit,pr,imatl,hstar,expn,sres,alpha,he,se,
     +	ipr2,ipr3,dtold2,t,tprint,errtrans,tmass,fluxin,flxout,
     +	fluxet,fluxev,sumtot,sumtotr,prold,princr,satold,trans,daysatm,
     +	satinstm,croot,qbuild,iqflag,dtransum,devsum,dleak,dqin,static,
     + 	ilowcros,imidcros,ihicros,ivhicros,numlow,nummid,numhi,numvhi,
     + 	durlow,durmid,durhi,durvhi,pi,csoil,sumlow,
     +  summid,sumhi,sumvhi)
     
      iprad = 1  

!*************************
!     Begin Time Loop
!*************************

!     New Time Step 
      
    1 if((t+dt).ge.(tprint-dt/10)) then
         dtold = dt
         dt = tprint-t
         t = tprint
         ipr = 1
         iprad = 1 
         ipr2 = ipr2+iprad
         ipr3 = ipr3+iprad 
         tprint = t+dtprint     
      
      elseif((t+dt).ge.(timet(nprtop)-dt/10))then
         dtold = dt
         dt =  timet(nprtop)-t 
         if (dt.lt.dtmin)then
            dt = dtold
         endif
         t =  timet(nprtop) 
         iprad = 0
      
      else
         t = t+dt
         iprad = 0
      end if
      
      itime = itime+1

!   Boundary condition counters
      call bcchk(bctop,timet,bcbot,timeb,nbctop,nbcbot,t,dtmin,
     +      iqflag,ibcflag,nprtop,nprbot,bct,bcb,nbct,nbcb,
     +		nprtad,nprbad)

!   Check for excedence of max time or max number of time steps
      call timchk(ntmax,itime,tmax,t,pr,nblocks)

!   Initialize iteration parameters
      call initer(iter,itest)

!   Iteration Loop
    2 iter = iter+1
      rmax = 0.d0

!   Check for excedence of max number of iterations
      call itchk(bctop,bcbot,nbctop,nbcbot,iprad,dtprint,imax,
     +		itmax,iter,dtredu,dtmin,nprbad,nprtad,rmax,
     +		prold,satold,nblocks,t,dt,tprint,
     +		bct,bcb,nbct,nbcb,ipr,ipr2,ipr3,pr,sat,itest,
     + 		itime,nprtop,nprbot)
        
      if(itest.eq.1) goto 1

!   Subroutine to calculate local transpiration rates, explicit in time
      call transpir(t,tday,etstart,etfin,pr,pi,rksat,
     + 		etmax,evmax,hwilt,csoil,croot,rtavgdens,dz,
     +		etmult,rweight,nroot,sres,expn,imatl,sat,
     +		trans, sumtrans)  

!   Generate Matrix Equation
      call matrix(dt,t,satold,trans,rksat,stor,dz,
     + 		evmult,eweight,nblocks,nbct,nbcb,bct,bcb,sat,
     +		pr,prinit,prold,expn,sres,alpha,he,se,zevap,pi,
     +	 	evmax,etfin,etstart,tday,swilt,imatl,poros, 
     +		aa,rhs)

!   Solve Equations
      call thomas(aa,rhs,princr,nblocks,nbdim)

!   Update iterated solution and check for convergence
      call cnvchk (errpr,errres,rmax,princr,nblocks,
     + 		pr,sat,imatl,ntest,expn,sres,alpha,he,se)

!   Goto next iteration if solution has not converged
      if(ntest.ge.1) goto 2

!   If solution has converged, then continue

!   Mass Balance
      call masbal(dt,trans,dz,evmult,etmult,eweight,nblocks,nroot,
     +		nbcb,nbct,rksat,poros,stor,prinit,expn,sres,imatl,
     +      alpha,he,se,zevap,pi,pr,prold,sat,sumtot,sumtotr,
     +	 	evmax,etfin,etstart,tday,bct,bcb,swilt,t,
     + 		ratio,ratio1,qin,qout,qleak,fluxin,flxout,
     +      fluxet,fluxev,transum,evsum)

!   Subroutine to calculate static stress explicitly in time
      call stressone(t,tday,pi,etstart,tmax,transum,dt,
     +		etmax,imatl,evmax, 
     +		ilowcros,imidcros,ihicros,ivhicros,
     +		numlow,nummid,numhi,numvhi,
     + 		durlow,durmid,durhi,durvhi,
     +		sumlow,summid,sumhi,sumvhi)

!   Update Solution Vector
      do i=1,nblocks
         prold(i) = pr(i)
         satold(i) = sat(i)
      end do

!   Write information to screen
      if(iscrn.eq.1) then
         write(*,20) itime,t,dt,iter
   20    format(1x,'time step',i8,',  time =',f10.3,
     +             ',  dt =',e10.4,',  number of iter',i4) 
      end if

!   Check for possible ponding
      if((nbct.eq.2.and.pr(1).gt.0.d0.and.ipond.gt.0).or.
     +   iqflag.eq.1) then
         call ponding(bctop,nbctop,nprtop,ipond,iscrn,dt,t,qin,pr,
     +		iqflag,bct,nbct,qbuild,qinflx)  
      end if

!   Compute and Print Daily-Integrated, Vertically Integrated Saturation
!         and Transpiration, and Daily Integrated Stress
      if(ipr3.eq.mprint) then
         ipr3 = 0
         call pprint3(qin,transum,evsum,qleak,dt,nblocks,t,dz,sat,pr,
     +              poros,stor,prinit,imatl,daysatm,dqin,etmax,evmax,
     +				dtransum,devsum,dleak)
      else 
         call dayavg(qin,transum,evsum,qleak,dt,nblocks,pr,dz,sat,
     +                  poros,stor,prinit,imatl,daysatm,dqin,
     +					dtransum,devsum,dleak)
      end if 
        
!   Print Saturation Profile
      if(ipr2.eq.mprint) then
         ipr2 = 0
         call pprint2(pr,dz,sat,t,ztop,nblocks,poros,prinit,
     + 					stor,imatl)
      end if 
    
!   Print Vertically-Integrated (over root zone) Saturation
      if(ipr.eq.1) then
         ipr = 0
         call pprint(transum,evsum,qleak,dt,nblocks,t,tday,pi,pr,
     +		etstart,fluxin,flxout,fluxet,fluxev,ratio,ratio1,etmax,       
     +		poros,prinit,stor,imatl,dz,sat,evmax)
         !if ((t.gt.100).and.(t.le.200))then
         !  call pprint4(dz,sat,t,ztop,zroot,fluxin,fluxet,fluxev)
         !endif
         dt = dtold
      end if

!   Adjust Time Step
      if(iter.lt.itmin) then
         dt = dt*dtincr
         if(dt.gt.dtmax) then
            dt = dtmax
         end if
      end if

      if(ipr.gt.-1) then
         goto 1
      end if

      stop
      end
