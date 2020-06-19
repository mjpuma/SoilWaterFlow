      subroutine stressone(t,tday,pi,etstart,tmax,transum,dt,
     +		etmax,imatl,evmax, 
     +		ilowcros,imidcros,ihicros,ivhicros,
     +		numlow,nummid,numhi,numvhi,
     + 		durlow,durmid,durhi,durvhi,
     +		sumlow,summid,sumhi,sumvhi)
        implicit none
!***********************************************************************
!     Routine:  stressone.for                                          *
!     Function: Static stress, stress duration, and # of Upcrossings   *
!               for each time step                                     *
!     Written:  Feb. - April 2003                                      *
!     By:       M.J.Puma                                               *
!***********************************************************************
	  integer*4, parameter :: prec10 = selected_real_kind(10)
      integer*4, parameter :: matdim =  5
      integer*4, parameter :: nbdim  = 320      
      
	  integer*4, parameter :: nulow= 28
      integer*4, parameter :: numid= 29
      integer*4, parameter :: nuhi= 30 
      integer*4, parameter :: nuvhi= 31 
      
	  real(kind = prec10) static
      integer*4 ilowcros,imidcros,ihicros, ivhicros
      integer*4 numlow, nummid, numhi, numvhi
      real(kind = prec10) durlow, durmid, durhi, durvhi 
      
	  real(kind = prec10) daytime,t,tday,tpot,etmax(matdim),
     +  evmax(matdim),pi,
     + 	etstart,errtpot,dt,transum,tmax,avgvhi,tvhiend,tvhibegin,
     +	sumvhi,avghi,thiend,thibegin,sumhi,avgmid,tmidend,tmidbegin,
     +	summid,avglow,tlowend,tlowbgin,sumlow,transdiff,errtime,
     +  timeend,stressq
      integer*4 imatl(nbdim) 

!   For diurnal ET case      
      daytime = t-INT(t/tday)*tday
!      tpot = (etmax(imatl(1)) - evmax(imatl(1)))*
!     +           (pi*sin(2*pi*(daytime-etstart)))
!      errtpot = 1.0e-8 
!      if (tpot.lt.errtpot) then
!         static = 0.d0                         
!      else   
!!	For temporally constant ET case        
         tpot = etmax(imatl(1)) - evmax(imatl(1)) 
      
!   Stress Parameters
         stressq = 1 
      
!   Maximum Time Check      
         timeend = tmax - t
         errtime = 1.0e-8  
              

!  Calculate the static stress based on transpiration     
         transdiff = tpot-transum
         if (transdiff.lt.0.d0) then
            static = 0
         else   
            static = ((tpot-transum)/tpot)**stressq
         endif       
  
!  Calculate duration of stress excursions between 0 and 0.2    
         if (static.lt.0.2) then
         
            if (ilowcros.eq.0) then  
               sumlow = static * dt
               ilowcros = 1
               tlowbgin = t
               numlow = numlow + 1      
            elseif(timeend.lt.errtime)then
               ilowcros = 0 
               tlowbgin  = t
               tlowend = tmax     
               durlow = tlowend - tlowbgin  
               avglow = sumlow/durlow
               write(nulow,100) numlow, tlowbgin, tlowend,
     +                       durlow, sumlow, avglow
               return  
            else
               sumlow = sumlow + static*dt
            endif    
      
         elseif (ilowcros.eq.1 )then
            ilowcros = 0 
            tlowend = t    
            durlow = tlowend - tlowbgin  
            avglow = sumlow/durlow
            write(nulow,100) numlow, tlowbgin, tlowend,
     +                       durlow, sumlow, avglow
         endif 
  
!  Calculate duration of stress excursions between 0.2 and 0.5
             
         if ((static.gt.0.2).and.(static.lt.0.5)) then
         
            if (imidcros.eq.0) then  
               summid = static * dt
               imidcros = 1
               tmidbegin = t
               nummid = nummid + 1 
            elseif(timeend.lt.errtime) then     
               imidcros = 0 
               tmidbegin  = t
               tmidend = tmax      
               durmid = tmidend - tmidbegin  
               avgmid = summid/durmid
               write(numid,120) nummid, tmidbegin, tmidend,
     +                       durmid, summid, avgmid
               return
            else
               summid = summid + static*dt
            endif    
      
         elseif  (imidcros.eq.1 )then
            imidcros = 0 
            tmidend = t    
            durmid = tmidend - tmidbegin  
            avgmid = summid/durmid
            write(numid,120) nummid, tmidbegin, tmidend,
     +                       durmid, summid, avgmid
         endif          
  
!  Calculate duration of stress excursions between 0.5 and 0.8
             
         if ((static.gt.0.5).and.(static.lt.0.8)) then
         
            if (ihicros.eq.0) then  
               sumhi = static * dt
               ihicros = 1
               thibegin = t
               numhi = numhi + 1      
            elseif(timeend.lt.errtime) then        
               ihicros = 0 
               thibegin  = t
               thiend = tmax     
               durhi = thiend - thibegin  
               avghi = sumhi/durhi
               write(nuhi,140) numhi, thibegin, thiend,
     +                       durhi, sumhi, avghi
               return
            else
               sumhi = sumhi + static*dt
            endif    
      
         elseif (ihicros.eq.1 )then
            ihicros = 0 
            thiend = t    
            durhi = thiend - thibegin  
            avghi = sumhi/durhi
            write(nuhi,140) numhi, thibegin, thiend,
     +                       durhi, sumhi, avghi
         endif         
  
!  Calculate duration of stress excursions greater than 0.8             
         if (static.gt.0.8) then
         
            if (ivhicros.eq.0) then  
               sumvhi = static * dt
               ivhicros = 1
               tvhibegin = t
               numvhi = numvhi + 1 
            elseif(timeend.lt.errtime) then
               ivhicros = 0 
               tvhibegin  = t
               tvhiend = tmax   
               durvhi = tvhiend - tvhibegin  
               avgvhi = sumvhi/durvhi
               write(nuvhi,160) numvhi, tvhibegin, tvhiend,
     +                       durvhi, sumvhi, avgvhi
               return         
            else
               sumvhi = sumvhi + static*dt
            endif    
      
         elseif (ivhicros.eq.1 )then
            ivhicros = 0 
            tvhiend = t    
            durvhi = tvhiend - tvhibegin  
            avgvhi = sumvhi/durvhi
            write(nuvhi,160) numvhi, tvhibegin, tvhiend,
     +                       durvhi, sumvhi, avgvhi
         endif         
              
  100    format(i10, 5f16.6)
  120    format(i10, 5f16.6)  
  140    format(i10, 5f16.6)
  160    format(i10, 5f16.6)         
  
!      endif
      return
      end         