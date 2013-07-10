        program getwave
*
*       compute wavelet transform and significances.
*       based on wavetest.f by Christopher Torrence and Gilbert P. Compo
*       for more info see http://paos.colorado.edu/research/wavelets/
*
        implicit none
#include "param.inc"
#include "getopts.inc"
#include "recfac.h"
        integer jtotmax,nvarmax,nsig
        parameter(jtotmax=88,nvarmax=4,nsig=100)
***        parameter(npadmax=2**18)
        integer i,j,k,k1,n,nperyear,ninter,i0,i1,subscale,jtot,j1,j2
     +        ,ibase2,npad,mother,m,yr,mo,npadmax
        real eta,p
        real,allocatable :: data(:,:,:)
        complex c
*       I am too lazy to convert to single precision, which should be
*       enough
        double precision,allocatable :: y(:),coi(:)
        double precision pi,dt,s0,dj,param
        double precision scale(jtotmax),period(jtotmax)
        double complex,allocatable :: wave(:,:)
        double precision clag0,clag1,siglvl(nsig),dof(jtotmax)
        double precision fft_theor(jtotmax),signif(jtotmax,nsig),ymean
     +        ,variance
        double precision Cdelta,psi0
*       real again
        integer nvars,ivars(2,nvarmax),iens,mens1,mens,iu
        real xx(1),yy(1),zz(jtotmax),ss(nsig),sl(nsig),power,signi,dy
        real,allocatable :: coifield(:,:),powfield(:,:),sigfield(:,:)
        character file*1024,ctlfile*1024,datfile*1024,title*132
        character vars(nvarmax)*10,lvars(nvarmax)*40,units(nvarmax)*20
        integer iargc
        real factln,gammln
        external factln,gammln
        pi = 4*atan(1.d0)
*
*       process arguments
*
        n = iargc()
        if ( n.lt.4 ) then
            print *,'usage: wave infile type param options outfile.ctl'
            stop
        endif
        call getarg(2,file)
        if ( file(1:6).eq.'morlet' ) then
            print '(a)','# using Morlet wavelet'
            mother = 0
        elseif ( file(1:4).eq.'paul' ) then
            print '(a)','# using Paul wavelet'
            mother = 1
        elseif ( file(1:3).eq.'dog' ) then
            print '(a)','# using DOG wavelet'
            mother = 2
        else
            write(0,*) 'error: please specify one of morlet|paul|dog'
            call abort
        endif
        call getarg(3,file)
        read(file,*,err=900) param
        print '(a,f6.2,a)','# using ',param,' wiggles/pattern'
        call getarg(1,file)
        allocate(data(npermax,yrbeg:yrend,0:nensmax))
        call readensseries(file,data,npermax,yrbeg,yrend,nensmax
     +       ,nperyear,mens1,mens,vars(1),units(1),lstandardunits
     +       ,lwrite)
        call getopts(4,n-1,nperyear,yrbeg,yrend,.true.,mens1,mens)
        if ( mens.gt.0 ) then
            write(0,'(a,i4,a,i4,a)') 'Averaging ensemble members ',nens1
     +            ,' to ',nens2,'<br>'
            write(*,'(a,i4,a,i4)') '# Averaging ensemble members ',nens1
     +            ,' to ',nens2
        endif
*
*       sum series
*
        if ( lsum.gt.1 ) then
	    print '(a,i4,a)','# summing over ',lsum,' periods'
            do iens=nens1,nens2
                call sumit(data(1,yrbeg,iens),npermax,nperyear
     +                ,yrbeg,yrend,lsum,oper)
            enddo
        endif
*       
*       log, sqrt
*
        if ( logscale ) then
	    print '(a)','# taking logarithm'
            do iens=nens1,nens2
                call takelog(data(1,yrbeg,iens),npermax,nperyear
     +                ,yrbeg,yrend)
            enddo
        endif
        if ( sqrtscale ) then
	    print '(a)','# taking sqrt'
            do iens=nens1,nens2
                call takesqrt(data(1,yrbeg,iens),npermax,nperyear
     +                ,yrbeg,yrend)
            enddo
        endif
*
*       differentiate data
*
        if ( ndiff.ne.0 ) then
            print '(a,i4)','# taking differences/averaging ',ndiff
            do iens=nens1,nens2
                call diffit(data(1,yrbeg,iens),npermax,nperyear
     +                ,yrbeg,yrend,ndiff)
            enddo
        endif
*
*       anomalies - necessary if we consider more than one month
*
        if ( lsel.gt.1 .and. ndiff.le.0 .or. anom ) then
            print '(a)','# taking anomalies'
            do iens=nens1,nens2
                call anomal(data(1,yrbeg,iens),npermax,nperyear,yrbeg
     +                ,yrend,yr1,yr2)
            enddo
        endif
*
        if ( m1.eq.0 ) then
            j1 = 1
            j2 = nperyear
            dt = 1./nperyear
        else
            j1 = m1
            j2 = m2
            call month2period(j1,nperyear,1)
            call month2period(j2,nperyear,0)
            if ( j2.lt.j1 ) j2 = j1 + nperyear
            dt = 1
        endif
*
*       find first, last valid data
*       
        do yr=yr1-1,yr2
            do iens=nens1,nens2
                do mo=j1,j2
                    j = mo
                    call normon(j,yr,i,nperyear)
                    if ( i.lt.yr1 .or. i.gt.yr2 ) cycle
                    if ( data(j,i,iens).lt.1e33 ) goto 100
                enddo
            enddo
        enddo
        goto 110
  100   continue
        yr1 = i
  110   continue
        do yr=yr2,yr1-1,-1
            do iens=nens1,nens2
                do mo=j1,j2
                    j = mo
                    call normon(j,yr,i,nperyear)
                    if ( i.lt.yr1 .or. i.gt.yr2 ) cycle
                    if ( data(j,i,iens).lt.1e33 ) goto 200
                enddo
            enddo
        enddo
        goto 210
  200   continue
        yr2 = i
  210   continue
        if ( lwrite ) print *,'found valid data in ',yr1,yr2
*       
*       loop over ensemble members
*
        npadmax = 4*(yr2-yr1+1)*nperyear
        allocate(y(npadmax))
        do iens=nens1,nens2
*
*           make a vector of values
*
            n = 0
            do i=yr1,yr2
                if ( j1.ne.j2 ) then
                    do j=1,nperyear
                        n = n + 1
                        if ( n.gt.npadmax ) then
                            do iu=0,6,6
                                write(iu,*) 'wave: error: linear array '
     +                               //'too small: ',npadmax
                            enddo
                            call abort
                        endif
                        if ( j.ge.j1 .and. j.le.j2 .or.
     +                       j.ge.j1-nperyear .and. j.le.j2-nperyear )
     +                       then
                            y(n) = data(j,i,iens)
                        else
                            y(n) = 3e33
                        endif
                        if ( lwrite) print *,'n,y(n) = ',n,y(n)
                    enddo
                else
                    n = n + 1
                    if ( n.gt.npadmax ) then
                        do iu=0,6,6
                            write(iu,*) 'wave: error: linear array '
     +                           //'too small: ',npadmax
                        enddo
                        call abort
                    endif
                    y(n) = data(j,i,iens)
                    if ( lwrite) print *,'n,y(n) = ',n,y(n)
                endif
            enddo
*
*           replace missing data by a linear interpolation
*
            ninter = 0
            i0 = -1
            i1 = -1
            do i=1,n
                if ( y(i).gt.1e33 ) then
*                   search next defined point
                    if ( i1.eq.-1 ) then
                        do j=i+1,n
                            if ( y(j).lt.1e33 ) then
                                i1 = j
                                goto 300
                            endif
                        enddo
*                       no more valid values - extrapolating is different
                        goto 400
  300                   continue
                    endif
                    if ( i0.gt.0 ) then
*                       interpolate
                        ninter = ninter + 1
                        y(i) = ((i-i0)*y(i1) + (i1-i)*y(i0))/(i1-i0)
                    endif
                else
*                   set previous defined point
                    i0 = i
                    i1 = -1
                endif
            enddo
  400       continue
*
*           determine total lenghth of vector.  Note that I add a longer
*           stretch of padding, at least equal to the data.
*
            ibase2 = int(log(dble(n))/log(2.d0)) + 2
            npad = nint(2.d0**ibase2)
            write(*,'(a,i5,a,i6,a)')
     +            '# The time series has been extended from ',n,' to ',
     +            npad,' elements with a straight line'
            if ( npad.gt.npadmax ) then
                do iu=0,6,6
                    write(iu,*)
     +                   'wave: error: array y too small, extend from '
     +                   ,npadmax,' to ',npad
                enddo
                call abort
            endif
*       
*       fill padding area with data interpolated between last and first
*       value.  i0 still has the position of the last valid datum, find
*       first valid one and interpolate
*       
            do i1=1,n
                if ( y(i1).lt.1e33 ) goto 410
            enddo
  410       continue
            do i=1,i1-1
                y(i) = ((i+npad-i0)*y(i1) + (i1-i)*y(i0))/(i1-i0+npad)
            enddo
            do i=i0+1,npad
                y(i) = ((i-i0)*y(i1) + (i1+npad-i)*y(i0))/(i1-i0+npad)
            enddo
*
*           Wavelet transform
*
*           print numbers for illustration plot
            if ( iens.eq.nens1 ) then
            if ( mother.eq.0 ) then
                do i=-400,400
                    eta = i/100.
                    print *,eta,
     +                    cos(param*eta)*exp(-eta**2/2)/sqrt(sqrt(pi)),
     +                    sin(param*eta)*exp(-eta**2/2)/sqrt(sqrt(pi))
                enddo
            elseif ( mother.eq.1 ) then
                m = nint(param)
                if (abs(param-m).gt.0.01 ) then
                    write(0,*)
     +                    'error: only works for integer parameter, not'
     +                    ,param
                    call abort
                endif
                do i=-400,400
                    eta = i/100.
                    c = CMPLX(0.,2.)**m * exp(factln(m) - factln(2*m)/2)
     +                    /sqrt(pi) * CMPLX(1.,-eta)**(-m-1)
                    print *,eta,REAL(c),AIMAG(c)
                enddo
            elseif ( mother.eq.2 ) then
                m = nint(param)
                if (abs(param-m).gt.0.01 ) then
                    write(0,*)
     +                    'error: only works for integer parameter, not'
     +                    ,param
                    call abort
                endif
                do i=-400,400
                    eta = i/100.
*                   hand-derived, I hope there are not too many errors
                    if ( m.eq.0 ) then
                        p = 1
                    elseif ( m.eq.1 ) then
                        p = eta
                    elseif ( m.eq.2 ) then
                        p = eta**2 - 1
                    elseif ( m.eq.3 ) then
                        p = eta*(eta**2 - 3)
                    elseif ( m.eq.4 ) then
                        p = eta**2*(eta**2 - 6) + 3
                    elseif ( m.eq.5 ) then
                        p = eta*(eta**2*(eta**2 - 10) + 15)
                    elseif ( m.eq.6 ) then
                        p = eta**2*(eta**2*(eta**2 - 15) + 45) - 15
                    elseif ( m.eq.7 ) then
                        p = eta*(eta**2*(eta**2*(eta**2 - 21) + 105) -
     +                        105)
                    elseif ( m.eq.8 ) then
                        p = eta**2*(eta**2*(eta**2*(eta**2 - 28) + 210)
     +                        -520) + 105
                    else
                        write(0,*) 'Sorry, DOG plot for m = ',m
     +                        ,' not yet ready'
                        goto 500
                    endif
                    print *,eta,p*exp(-eta**2/2)*exp(-gammln(m+0.5)/2),0
     +                    .
                enddo
  500           continue
            else
                write(0,*) 'Other wavelets not yet supported',mother
            endif
            endif
*       
*           a few more parameters; I follow the comments in wavetest.f
*
            if ( mother.eq.0 ) then
                s0 = dt
            elseif ( mother.eq.1 ) then
                s0 = dt/4
            else
                s0 = 2*dt
            endif
            subscale = 8
            dj=1./subscale
            jtot=11*subscale
            if ( jtot.gt.jtotmax ) then
                write(0,*) 'wave: error: increase jtotmax to ',jtot
                call abort
            endif
*       
*           get the wavelet transform
*       
            print '(a)','# computing wavelet transform...'
            if ( lwrite ) then
                print *,'calling wave with'
                print *,'jtot,jtotmax = ',jtot,jtotmax
                do i=1,npad
                    print *,i,y(i)
                enddo
            endif
            if ( iens.eq.nens1 ) then
                allocate(wave(npadmax,jtotmax))
                allocate(coi(npadmax))
            end if
            call wavelet(n,npadmax,y,dt,mother,param,s0,dj,jtot,npad,
     &            wave,scale,period,coi)
            if ( lwrite ) then
                print *,'wave returned'
                do i=1,npad
                    print *,i,(wave(i,jtot),j=1,jtot)
                enddo
            endif
*           convert coi to a field
            if ( iens.eq.nens1 ) then
                allocate(coifield(npadmax,jtotmax))
                allocate(powfield(npadmax,jtotmax))
                allocate(sigfield(npadmax,jtotmax))
                do i=1,n
                    do j=1,jtot
                        if ( period(j).gt.coi(i) ) then
                            coifield(i,j) = 0
                        else
                            coifield(i,j) = 1
                        endif
                    enddo
                enddo
            endif
*       
*           local significance tests
*       
            ymean = 0
            do i=1,n
                ymean = ymean + y(i)
            enddo
            ymean = ymean/n
            clag0 = 0
            clag1 = 0
            do i=1,n-1
                clag0 = clag0 + (y(i)-ymean)**2
                clag1 = clag1 + (y(i)-ymean)*(y(i+1)-ymean)
            enddo
            clag1 = clag1/clag0
            if ( lwrite ) print *,'clag1 = ',clag1
            do k=1,nsig
                siglvl(k) = 10.d0**(-5.d0*k/nsig)
                if ( lwrite ) print *,'calling wave_signif',k,siglvl(k)
                call wave_signif(0,n,y,dt,mother,param,dj,jtot,
     +                scale,period,clag1,siglvl(k),dof,fft_theor
     +                ,signif(1,k),ymean,variance,Cdelta,psi0,lwrite)
                if ( lwrite ) print *,'back from wave_signif'
            enddo
*       
*           collect information from ensemble members
*       
            do j=1,jtot
                do k=1,nsig
                    ss(k) = signif(j,k)
                    sl(k) = siglvl(k)
                enddo
                do i=1,n
                    power = abs(wave(i,j))**2
                    call hunt(ss,nsig,power,k1)
                    k1 = min(max(1,k1-2),n-3)
                    call polint(ss(k1),sl(k1),4,power,signi,dy)
                    if ( iens.eq.nens1 ) then
                        powfield(i,j) = power
                        sigfield(i,j) = signi
                    else
                        powfield(i,j) = powfield(i,j) + power
                        sigfield(i,j) = 3e33 ! I have to think about this
                    endif
                enddo
            enddo
        enddo
        if ( nens2.gt.nens1 ) then
            do j=1,jtot
                do i=1,n
                    powfield(i,j) = powfield(i,j)/(nens2-nens1)
                    sigfield(i,j) = -3e33
                enddo
            enddo
        endif
        print '(a)','# writing output...'
        call getarg(iargc(),ctlfile)
        datfile = ctlfile(1:index(ctlfile,'.ctl '))//'grd'
*       convert to single precision
        xx(1) = 0
        yy(1) = 0
        do j=1,jtot
            zz(j) = period(j)
        enddo
        nvars = 3
        vars(1) = 'power'
        units(1) = '('//trim(units(1))//')^2'
        lvars(1) = 'power of wavelet transform'
        ivars(1,1) = jtot
        ivars(2,1) = 99
        vars(2) = 'coi'
        lvars(2) = 'cone of influence'
        units(2) = '1'
        ivars(1,2) = jtot
        ivars(2,2) = 99
        vars(3) = 'prob'
        lvars(3) = 'p=value'
        units(3) = '1'
        ivars(1,3) = jtot
        ivars(2,3) = 99
        do i=len(file),1,-1
            if ( file(i:i).eq.'/' ) goto 777
        enddo
  777   continue
        call args2title(title)
        call writectl(ctlfile,datfile,1,xx,1,yy,jtot,zz,n,nint(1/dt),
     +        yr1,j1,3e33,title,nvars,vars,ivars,lvars,units)
        ctlfile=ctlfile(1:index(ctlfile,'.ctl')-1)//'f.ctl'
        do j=1,jtot
            zz(j) = 1/zz(j)
        enddo
        call writectl(ctlfile,datfile,1,xx,1,yy,jtot,zz,n,nint(1/dt),
     +        yr1,j1,3e33,title,nvars,vars,ivars,lvars,units)
        open(1,file=datfile,access='direct',recl=jtot*recfa4,
     +        status='old',err=800)
        close(1,status='delete')
  800   continue
        open(1,file=datfile,access='direct',recl=3*jtot*recfa4,
     +        status='new')
        do i=1,n
            write(1,rec=i) (powfield(i,j),j=1,jtot),
     +            (coifield(i,j),j=1,jtot),
     +            (sigfield(i,j),j=1,jtot)
        enddo
        close(1)
        call savestartstop(yr1,yr2)
        goto 999
  900   write(0,*) 'Expected real param, but found ',file
        call abort
  999   continue
        end
*  #[ rest:
#if 0
        subroutine restwavetest
C these parameters depend on the particular time series
      PARAMETER (n=504,dt=0.25D0,s0=dt)
      PARAMETER (subscale=4)
      PARAMETER (dj=1.D0/subscale,jtot=11*subscale)
C Note: for accurate reconstruction and wavelet-derived variance
C     do not pad with zeroes, set s0=dt (for Paul set s0=dt/4), and use
C     a large "jtot" (even though the extra scales will be within
C     the cone of influence).
C     For plotting purposes, it is only necessary to use
C     s0=2dt (for Morlet) and "jtot" from Eqn(10) Torrence&Compo(1998).

      INTEGER mother,ibase2,npad
      DOUBLE PRECISION sst(n),recon_sst(n),param,pi
      DOUBLE PRECISION scale(jtot),period(jtot),coi(n)
      DOUBLE COMPLEX wave(n,jtot)

      INTEGER i,j,isigtest,javg1,javg2
      DOUBLE PRECISION lag1,siglvl,dof(jtot)
      DOUBLE PRECISION fft_theor(jtot),signif(jtot),ymea,nvariance
      DOUBLE PRECISION recon_mean,recon_vari
      DOUBLE PRECISION Cdelta,psi0
      DOUBLE PRECISION global_ws(jtot),global_signif(jtot)
      DOUBLE PRECISION savg_dof(jtot),savg_signif(jtot),sstENSO(n)
        
C*************************************************** Significance testing

C** local significance test
      isigtest = 0
      lag1 = 0.72D0
      siglvl = 0.05D0
      CALL WAVE_SIGNIF (isigtest,n,sst,dt,mother,param,dj,jtot,
     &	     scale,period,lag1,siglvl,dof,fft_theor,signif,
     &	     ymean,variance,Cdelta,psi0)


C** global wavelet spectrum & significance test
      isigtest = 1
      lag1 = 0.72D0
      siglvl = 0.05D0
      DO 10 j=1,jtot
	DO 20 i=1,n
	  global_ws(j) = global_ws(j) + ABS(wave(i,j))**2
20	CONTINUE
	global_ws(j) = global_ws(j)/n
	dof(j) = n - scale(j)
10    CONTINUE

      CALL WAVE_SIGNIF (isigtest,n,sst,dt,mother,param,dj,jtot,
     &	     scale,period,lag1,siglvl,dof,fft_theor,global_signif,
     &	     ymean,variance,Cdelta,psi0)


C** scale-average time series & significance test
      isigtest = 2
      lag1 = 0.72D0
      siglvl = 0.05D0
C    scale average between 2 and 7.9 years
      savg_dof(1) = 2.0D0
      savg_dof(2) = 7.9D0
C    find the "j"-values that correspond to savg_dof(1) & savg_dof(2)
      javg1 = 0
      javg2 = 0
      DO 30 j=1,jtot
	IF ((scale(j).GE.savg_dof(1)).AND.(javg1.EQ.0)) javg1 = j
	IF (scale(j).LE.savg_dof(2)) javg2 = j
30    CONTINUE
C   call wave_signif first, to get the value of "Cdelta"
      CALL WAVE_SIGNIF (isigtest,n,sst,dt,mother,param,dj,jtot,
     &	   scale,period,lag1,siglvl,savg_dof,fft_theor,savg_signif,
     &	   ymean,variance,Cdelta,psi0)
C   construct the scale-averaged time series [Eqn(24)]
      DO 50 i=1,n
	sstENSO(i) = 0.D0
	DO 60 j=javg1,javg2
	  sstENSO(i) = sstENSO(i) + (ABS(wave(i,j))**2)/scale(j)
60	CONTINUE
	sstENSO(i) = dj*dt*sstENSO(i)/Cdelta
50    CONTINUE


C************************************************************* print results
      PRINT*,' n=',n
      PRINT*,' dt=',dt
      PRINT*,' mother=',mother
      PRINT*,' param=',param
      PRINT*,' s0=',s0
      PRINT*,' dj=',dj
      PRINT*,' jtot=',jtot
      PRINT*,' npad=',npad
      PRINT'(/,"Let w = wave(n/2,j)",/)'
      PRINT'(A4,7A10)',"j","Scale","Period","ABS(w)^2","phase(w)",
     &  "5%signif","Global","GWS5%sig"
      PRINT'(I4,7F10.3)',(j,scale(j),period(j),
     &   ABS(wave(n/2,j))**2,
     &   ATAN2(DIMAG(wave(n/2,j)),DBLE(wave(n/2,j)))*180.D0/pi,
     &   signif(j),global_ws(j),global_signif(j),j=1,jtot)
      PRINT'(/,A,F10.3)',
     &    ' Scale-average degrees of freedom = ',savg_dof(1)
      PRINT'(A,F10.3,/)',
     &    ' Scale-avg 5% significance level  = ',savg_signif(1)


C************************************************************ Reconstruction

C** construct the wavelet derived variance (Parseval's theorem)  [Eqn(14)]
C   Cdelta & psi0 are returned from WAVE_SIGNIF
      recon_vari = 0.D0
      DO 900 i=1,n
	DO 1000 j=1,jtot
	  recon_vari = recon_vari + (ABS(wave(i,j))**2)/scale(j)
1000	CONTINUE
900   CONTINUE
      recon_vari = dj*dt*recon_vari/(Cdelta*n)
      PRINT'(A,F14.5)',' Reconstructed variance=',recon_vari
      PRINT'(A,F14.5)',' Original variance   =',variance
      PRINT'(A,F14.5,A,/)',' Ratio = ',recon_vari/variance,
     &     ' (this is low due to padding with zeroes)'

C** reconstruct the time series [Eqn(11)]
C   check mean and RMS difference of reconstructed time series
      recon_mean=0.D0
      recon_vari = 0.D0
      DO 1100 i=1,n
	recon_sst(i)=0.D0
	DO 1200 j=1,jtot
	  recon_sst(i) = recon_sst(i)+(DBLE(wave(i,j)))/SQRT(scale(j))
1200	CONTINUE
	recon_sst(i) = dj*SQRT(dt)*recon_sst(i)/(Cdelta*psi0)
	recon_vari = recon_vari+(sst(i)-ymean-recon_sst(i))**2
	recon_mean = recon_mean + recon_sst(i)
1100  CONTINUE
      recon_mean = recon_mean/n
      recon_vari = SQRT(recon_vari/n)

      PRINT'(A,F14.6)',' Reconstructed mean=',recon_mean
      PRINT'(A,F14.6)',' Original mean   =',ymean
      PRINT'(A,F14.6,/)',' Root-mean-square difference of time series=',
     &      recon_vari

      END
#endif
*  #] rest:
