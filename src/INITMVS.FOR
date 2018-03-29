      subroutine initmvs(bctop,bcbot,nbctop,nbcbot,tstart,dt,dtprint,
     + 	bct,bcb,nprtop,nprbot,nbct,nbcb,ibcflag,itime,ipr,sat,rksat,
     + 	etmax,evmax,frac,hwilt,rtavgdens,zroot,dz,etmult,nblocks,
     +  poros,stor,prinit,pr,imatl,hstar,expn,sres,alpha,he,se,
     +	ipr2,ipr3,dtold2,t,tprint,errtrans,tmass,fluxin,flxout,
     +	fluxet,fluxev,sumtot,sumtotr,prold,princr,satold,trans,daysatm,
     +	satinstm,croot,qbuild,iqflag,dtransum,devsum,dleak,dqin,static,
     + 	ilowcros,imidcros,ihicros,ivhicros,numlow,nummid,numhi,numvhi,
     + 	durlow,durmid,durhi,durvhi,pi,csoil,sumlow,
     +  summid,sumhi,sumvhi)
        implicit none
!***********************************************************************
!     Routine:  initmvs.for                                            *
!     Function: Initialize certain Matrics, Vectors, and Scalars.      *
!     Written:  January 2001, February 2003                            *
!     By:       M. Celia, M.J. Puma                                    *
!***********************************************************************
	  integer*4, parameter :: prec10 = selected_real_kind(10)

      integer*4, parameter :: nbdim  = 320
      integer*4, parameter :: ntdim  = 14020
      integer*4, parameter ::  matdim =  5

      integer*4, parameter :: nucon  = 21
      integer*4, parameter :: nupl1  = 22
      integer*4, parameter :: nupl2  = 23

      real(kind = prec10) bctop(ntdim), bcbot(ntdim),
     +               bct, bcb
      integer*4 nbctop(ntdim), nbcbot(ntdim), 
     +               nprtop, nprbot, nbct, nbcb, ibcflag     


      integer*4  itime, ipr, ipr2, ipr3
      real(kind = prec10) tstart, dt, dtold2, t
      real(kind = prec10) tprint, dtprint
      real(kind = prec10) errtrans 
      
      real(kind = prec10) tmass, fluxin, flxout, fluxet, fluxev, 
     +               sumtot, sumtotr 
      
      real(kind = prec10) pr(nbdim), prold(nbdim), princr(nbdim), 
     +               sat(nbdim), satold(nbdim), trans(nbdim),
     +               daysatm(3,nbdim), satinstm(3,nbdim) 
      
      real(kind = prec10)  rksat(matdim)
      real(kind = prec10) etmax(matdim), evmax(matdim),frac, 
     +                hwilt(matdim), hstar(matdim),tpotnoon            
                        
      real(kind = prec10) csoil, croot, rtavgdens(matdim)
      integer*4  imatl(nbdim)
      
      real(kind = prec10) qbuild
      integer*4  iqflag
      
      real(kind = prec10)  zroot, zrn, dz(nbdim), etmult(nbdim)
      integer*4  nblocks,i
     
      real(kind = prec10)  dtransum, devsum, dleak, dqin

      real(kind = prec10) static
      integer*4 ilowcros,imidcros,ihicros, ivhicros
      integer*4 numlow, nummid, numhi, numvhi
      real(kind = prec10) durlow, durmid, durhi, durvhi, sumlow,
     + 		summid,sumhi,sumvhi 

 	  real(kind = prec10) pi,fncsat,fnckrel,fncpor
      
      real(kind = prec10) poros(matdim),stor(matdim)

      real(kind = prec10) prinit(nbdim),expn(matdim),sres(matdim),
     +       alpha(matdim),he(matdim),se(matdim)

!   Initialize
      zrn = 0.d0
      pi = 4.d0*atan(1.d0)
!      ne = nblocks
!      hwilt = fncpr(swilt(imatl(1)),imatl(1))
!      hstar = fncpr(sstar(imatl(1)),imatl(1))
!      tsstar = 0.d0
      do i=1,nblocks
         prold(i) = pr(i)
         princr(i) = 0.d0
         satold(i) = sat(i)
         trans(i) = 0.d0
         zrn = zrn + fncpor(pr(i),imatl(i),poros,stor,prinit(i))
     +    			*dz(i)*etmult(i)
      end do
      ! Diurnal and uniform potential
!      tpotnoon = pi*(etmax(imatl(1))-evmax(imatl(1)))
      tpotnoon = etmax(imatl(1))-evmax(imatl(1))
      csoil = (zroot/tpotnoon)*((fnckrel(fncsat(hstar,imatl(1),
     +			expn,sres,alpha,he,se),imatl(1),expn,sres)
     +         *rksat(1)*rtavgdens(1))/
     +         (1.d0 - fnckrel(fncsat(hstar,imatl(1),     
     +			expn,sres,alpha,he,se),imatl(1),expn,sres)))*
     +         (hstar(imatl(1))-hwilt(imatl(1))+
     +         hwilt(imatl(1))*frac)
      
      croot = (zroot/tpotnoon)*(-hwilt(imatl(1))*frac*rtavgdens(1))
     +          - (csoil/rksat(1))
                                                
!      ctrans = zroot*( (hstar(imatl(1))-hwilt(imatl(1))) /
!     +         tpotnoon+hwilt(imatl(1))*frac/tpotnoon)*
!     +         (fnckrel(fncsat(hstar,imatl(1)),imatl(1)) /
!     +         (1.d0 - fnckrel(fncsat(hstar,imatl(1)),imatl(1))))
!      rplant = -hwilt(imatl(1))*zroot*frac/tpotnoon - ctrans

      rewind(nupl1)
      rewind(nupl2)
      rewind(nucon)

      itime = 0
      ipr = 0
      ipr2 = 0 
      ipr3 = 0
      t = tstart
      tprint = tstart+dtprint
      dtold2 = dt
      
      nprtop = 1
      nprbot = 1
      nbct = nbctop(1)
      bct = bctop(1)
      nbcb = nbcbot(1)
      bcb = bcbot(1)
      ibcflag = 0
      
      tmass = 0.d0
      fluxin = 0.d0
      flxout = 0.d0
      fluxet = 0.d0
      fluxev  =0.d0
      sumtot = 0.d0
      sumtotr = 0.d0
      qbuild = 0.d0
      iqflag = 0  
      errtrans = 0

!  Initialize Stress Variables   
      static = 0.d0
        
      ilowcros = 0
      imidcros = 0
      ihicros = 0
      ivhicros = 0
      
      durlow = 0.d0
      durmid = 0.d0
      durhi = 0.d0
      durvhi = 0.d0
     
      numlow = 0 
      nummid = 0      
      numhi = 0
      numvhi = 0     

	  sumlow = 0.d0
      summid=0.d0
      sumhi=0.d0
      sumvhi=0.d0
      
!  Initialize Daily Average Variables         
      dtransum = 0.d0
      devsum = 0.d0
      dleak = 0.d0
      dqin = 0.d0
      !dstatic = 0.d0
 
      do i = 1,nblocks
         satinstm(1,i) = 0.d0
         satinstm(2,i) = 0.d0
         satinstm(3,i) = 0.d0      
         daysatm(1,i) = 0.d0
         daysatm(2,i) = 0.d0
         daysatm(3,i) = 0.d0
      enddo       



      return
      end
