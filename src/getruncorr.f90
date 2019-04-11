subroutine getruncorr(dindx,ddata,lfirst,dddata,ndata,j1,j2 &
    ,lag,k,month,nperyear,imens,indxmx,indx &
    ,data,npermax,yrbeg,yrend,nensmax,n0 &
    ,string,lboot,lprint,rmin,rmax,zdif)

!       compute running correlations

    implicit none
    include 'getopts.inc'

    integer :: lperyear,lyr,lensmax
    parameter(lperyear=366,lyr=200,lensmax=31)
    integer :: ndata,j1,j2,lag,k,month,nperyear,indxmx,npermax,yrbeg,yrend,nensmax,n0
    integer :: imens(0:indxmx)
    real :: dindx(ndata),ddata(ndata),dddata(ndata)
    logical :: lfirst(ndata)
    real :: data(npermax,yrbeg:yrend,0:nensmax),indx(npermax,yrbeg:yrend,0:nensmax,indxmx)
    logical :: lboot,lprint
    real :: rmin,rmax,zdif,savedata(lperyear,lyr,0:lensmax), &
        saveindx(lperyear,lyr,0:lensmax),sig(1),a,b,siga,sigb,chi2,q
    character*(*) string

    integer :: yr1s,yr2s,yr1d,yr2d,yr1i,yr2i,n,yr,yrp,m,iiens,iens,jens,i
    integer :: yrmo(3,1)       ! dummy argument, not used
    real :: result,dresult(-2:2),prob,zmin,zmax,z
    logical :: lblank,llwrite
    integer :: init,iu
    save init
    data init /0/

    lblank = lprint
    if ( lprint ) then
        write(14,'(2a)') '# ',string
        write(14,'(a)') '#year  #pts   corr        prob '// &
            '       2.5%         16%         50%'// &
            '         84%       97.5%'
    end if
    if ( .false. ) then
        call rsunit(iu)
        write(string,'(a,i3.3,a)') 'check',init,'.txt'
        open(iu,file=string)
        init = init + 1
        write(iu,'(a,i3)') 'at beginning of getruncorr',k
        do iiens=nens1,nens2
            if ( imens(0) > 0 ) then
                iens = iiens
            else
                iens = 0
            end if
            if ( imens(k) > 0 ) then
                jens = iiens
            else
                jens = 0
            end if
            do yr=yr1,yr2
                do m=1,nperyear
                    if ( data(m,yr,iens) < 3e33 .or. &
                    indx(m,yr,jens,k) < 3e33 ) then
                        write(iu,'(i3,i5,i3,2g20.12)') iiens,yr,m, &
                            data(m,yr,iens),indx(m,yr,jens,k)
                    end if
                end do
            end do
            close(iu)
        end do
    end if
    yr1s = yr1
    yr2s = yr2
    rmin = +3e33
    rmax = -3e33
    zmin = +3e33
    zmax = -3e33
    do yr=yr1s,yr2s
        if ( irunvar > 0 ) then
!           moving window
            yr1 = max(yrbeg,yr - nyrwindow/2 + 1)
            yr2 = min(yrend,yr1 + nyrwindow - 1)
        else
!           running start date, fixed end date
            yr1 = yr
            yr2 = yr2s
            if ( yr2-yr1 < 4 ) cycle
        end if
        if ( ldetrend ) then
!           rescue data
            if ( nperyear > lperyear .or. yr2-yr1+1 > lyr ) then
                write(0,*) 'getruncorr: error: fixed array '// &
                    'too small; ',nperyear,lperyear,yr2-yr1+1,lyr
                call exit(-1)
            end if
            if ( lag == 0 ) then
                yr1i = yr1
                yr2i = yr2
                yr1d = yr1
                yr2d = yr2
            else if ( fix2 ) then
                yr1i = yr1
                yr2i = yr2
                if ( lag > 0 ) then
                    yr1d = yr1d + (lag-1)/nperyear + 1
                    yr2d = yr2d + (lag-1)/nperyear + 1
                else
                    yr1d = yr1d + (lag+1)/nperyear - 1
                    yr2d = yr2d + (lag+1)/nperyear - 1
                end if
            else
                yr1d = yr1
                yr2d = yr2
                if ( lag > 0 ) then
                    yr1i = yr1i - (lag-1)/nperyear - 1
                    yr2i = yr2i - (lag-1)/nperyear - 1
                else
                    yr1i = yr1i - (lag+1)/nperyear + 1
                    yr2i = yr2i - (lag+1)/nperyear + 1
                end if
            end if
            if ( max(yr2d-yr1d+1,yr2i-yr1i+1) > lyr ) then
                write(0,*) 'getruncorr: error: not enough room', &
                    ' to rescue detrended data, increase lyr to ', &
                    max(yr2d-yr1d+1,yr2i-yr1i+1)
                call exit(-1)
            end if
            do iiens=nens1,nens2
                if ( imens(0) > 0 ) then
                    iens = iiens
                else
                    iens = 0
                end if
                if ( imens(k) > 0 ) then
                    jens = iiens
                else
                    jens = 0
                end if
                if ( max(iens,jens) > lensmax ) then
                    write(0,*) 'getruncorr: error: increase lensmax' &
                        ,iens,jens,lensmax
                    call exit(-1)
                end if
                do yrp=yr1d,yr2d
                    do m=1,nperyear
                        savedata(m,yrp-yr1d+1,iens) = data(m,yrp,iens)
                    end do
                end do
                do yrp=yr1i,yr2i
                    do m=1,nperyear
                        saveindx(m,yrp-yr1i+1,jens) = indx(m,yrp,jens,k)
                    end do
                end do
            end do
            call detrend(data(1,yr1d,iens),npermax,nperyear,yr1d &
                ,yr2d,yr1d,yr2d,m1,m2,lsel)
            call detrend(indx(1,yr1i,jens,k),npermax,nperyear, &
                yr1i,yr2i,yr1i,yr2i,m1,m2,lsel)
        end if
        n = 0
        llwrite = lwrite
        lwrite = .false. 
        call filllinarray(dindx,ddata,lfirst,dddata,ndata,n,j1,j2 &
            ,lag,k,nperyear,imens,indxmx,indx,data,npermax,yrbeg &
            ,yrend,nensmax,-999,-999,yrmo)
        if ( ldetrend ) then
!           restore data
            do iiens=nens1,nens2
                if ( imens(0) > 0 ) then
                    iens = iiens
                else
                    iens = 0
                end if
                if ( imens(k) > 0 ) then
                    jens = iiens
                else
                    jens = 0
                end if
                do yrp=yr1d,yr2d
                    do m=1,nperyear
                        data(m,yrp,iens) = savedata(m,yrp-yr1d+1,iens)
                    end do
                end do
                do yrp=yr1i,yr2i
                    do m=1,nperyear
                        indx(m,yrp,jens,k) = saveindx(m,yrp-yr1i+1,jens)
                    end do
                end do
            end do
        end if
        lwrite = llwrite
        if ( minnum <= 0 ) then
            if ( month == 0 ) then
                m = 12*(nyrwindow-1)
            else
                m = (nyrwindow-1)
            end if
        else
            m = minnum
        end if
        if ( n < (1+nens2-nens1)*m ) then
            if ( lwrite ) print *,'getruncorr: not enough points:', &
                n,(1+nens2-nens1)*m
            if ( lblank ) write(14,'(a)')
            lblank = .false. 
            goto 700
        else if ( lwrite ) then
            print *,'getruncorr: enough points:',n,max(m,minnum)
        end if
        lblank = lprint
        if ( abs(irunvar) == 1 ) then
!           running correlatons
            call printcorr(dindx,ddata,lfirst,dddata,yrmo,n,n0,j1,j2 &
                ,month,nperyear,lag,string,lboot, .false. ,result &
                ,dresult,prob)
            if ( lprint .and. result < 1e33 ) then
                write(14,'(i4,i6,99g12.4)') yr,n,result,prob,dresult
            end if
            if ( lwrite ) then
                write(*,'(i4,i6,99g12.4)') yr,n,result,prob,dresult
            end if
            if ( abs(result) < 1e30 ) then
                rmin = min(rmin,result)
                rmax = max(rmax,result)
                if ( abs(result) < 1 ) then
                    z = (1+result)/(1-result)
                    z = 0.5*log(z)
                else if ( result > 0 ) then
                    z = +1e10
                else
                    z = -1e10
                end if
                zmin = min(zmin,z)
                zmax = max(zmax,z)
            end if
        else if ( abs(irunvar) == 2 ) then
!           running regressions, assume normality for the time being
            call fit(dindx,ddata,n,sig,0,a,b,siga,sigb,chi2,q)
            if ( lprint .and. b < 1e33 ) then
                write(14,'(i4,i6,99g12.4)') yr,n,b,q,b-2*sigb,b-sigb &
                    ,b,b+sigb,b+2*sigb
            end if
            if ( lwrite ) then
                do i=1,n
                    print '(i4,2g16.8)',i,dindx(i),ddata(i)
                end do
                write(*,'(i4,i6,99g12.4)') yr,n,b,q,b-2*sigb,b-sigb &
                    ,b,b+sigb,b+2*sigb
            end if
            rmin = min(rmin,b)
            rmax = max(rmax,b)
            zmin = rmin
            zmax = rmax
        else
            write(0,*) 'getruncorr: error: irunvar = ',irunvar
            write(*,*) 'getruncorr: error: irunvar = ',irunvar
            call exit(-1)
        end if
        if ( lprint ) call keepalive(yr-yr1s,yr2s-yr1s)
700     continue
    end do
    if ( rmax < -1e33 ) then
        rmax = 3e33
        zdif = 3e33
    else
        zdif = zmax - zmin
    end if
    if ( lprint ) write(14,'(a)')
    yr1 = yr1s
    yr2 = yr2s
    if ( .false. ) then
        call rsunit(iu)
        write(string,'(a,i3.3,a)') 'check',init,'.txt'
        open(iu,file=string)
        init = init + 1
        write(iu,'(a,i3)') 'at end of getruncorr',k
        do iiens=nens1,nens2
            if ( imens(0) > 0 ) then
                iens = iiens
            else
                iens = 0
            end if
            if ( imens(k) > 0 ) then
                jens = iiens
            else
                jens = 0
            end if
            do yr=yr1,yr2
                do m=1,nperyear
                    if ( data(m,yr,iens) < 3e33 .or. &
                    indx(m,yr,jens,k) < 3e33 ) then
                        write(iu,'(i3,i5,i3,2g20.12)') iiens,yr,m, &
                        data(m,yr,iens),indx(m,yr,jens,k)
                    end if
                end do
            end do
        end do
        close(iu)
    end if
end subroutine getruncorr
