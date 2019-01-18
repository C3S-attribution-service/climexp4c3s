subroutine project(var,nperyear,yr1,yr2,nens1,nens2,xx,nx,yy,ny, &
    field,pattern,nxf,nyf,npermax,firstyr,lastyr,month,minfac, &
    lnomissing,lonlymonth,anom,lwrite)

!   project a field on a pattern to obtain a time series

    implicit none
    integer :: nperyear,yr1,yr2,nens1,nens2,nx,ny,nxf,nyf,npermax,firstyr,lastyr,month,iens
    real :: var(nperyear,yr1:yr2,0:nens2),xx(nx),yy(ny), &
        field(nxf,nyf,npermax,firstyr:lastyr,0:nens2), &
        pattern(nxf,nyf,npermax,0:1),minfac
    logical :: lnomissing,lonlymonth,anom,lwrite
    real :: zz(1)

    zz(1) = 0
    call  project3(var,nperyear,yr1,yr2,nens1,nens2,xx,nx,yy,ny &
        ,zz,1,field,pattern,nxf,nyf,1,npermax,firstyr,lastyr &
        ,month,minfac,lnomissing,lonlymonth,anom,lwrite)
end subroutine project

subroutine project3(var,nperyear,yr1,yr2,nens1,nens2,xx,nx,yy,ny &
    ,zz,nz,field,pattern,nxf,nyf,nzf,npermax,firstyr,lastyr &
    ,month,minfac,lnomissing,lonlymonth,anom,lwrite)

!   project a field on a pattern to obtain a time series
!   3D version of proejct

    implicit none
    integer :: nperyear,yr1,yr2,nens1,nens2,nx,ny,nz,nxf,nyf,nzf &
        ,npermax,firstyr,lastyr,month,iens
    real :: var(nperyear,yr1:yr2,0:nens2),xx(nx),yy(ny),zz(nz), &
        field(nxf,nyf,nzf,npermax,firstyr:lastyr,0:nens2), &
        pattern(nxf,nyf,nzf,npermax,0:1),minfac
    logical :: lnomissing,lonlymonth,anom,lwrite

    integer :: yr01,mo,m,y,yr,i,j,k,m1,m2
    integer,allocatable :: nn(:,:,:,:)
    real :: s,w,wt
    real,allocatable :: wx(:),wy(:),wz(:),mean(:,:,:,:),avemean(:)
    logical :: xrev,yrev,xwrap
    logical,allocatable :: pointhasdata(:,:,:)

!       get weights and indices

    if ( lwrite ) then
        print *,'project: input'
        print *,'nperyear,yr1,yr2,nens1,nens2 = ',nperyear,yr1,yr2,nens1,nens2
        print *,'xx,nx,yy,ny = ',(xx(i),i=1,min(3,nx)),'...',xx(nx),nx, &
            (yy(i),i=1,min(3,ny)),'...',yy(ny),ny
        print *,'nxf,nyf,month,minfac = ',nxf,nyf,month,minfac
        print *,'lnomissing = ',lnomissing
    endif
    allocate(wx(nx))
    allocate(wy(ny))
    allocate(wz(ny))
    call getxyprop(xx,nx,yy,ny,xrev,yrev,xwrap)
    call getweights('x',xx,wx,nx,xwrap,lwrite)
    call getweights('y',yy,wy,ny, .false. ,lwrite)
    call getweights('z',zz,wz,nz, .false. ,lwrite)
    yr01 = 1
    if ( month == 0 ) then
        yr01 = 0
        month = nperyear
    endif
    if ( lonlymonth ) then
        call getj1j2(m1,m2,month,nperyear, .false. )
    else
        m1 = 1
        m2 = nperyear
    endif
    if ( lwrite ) then
        print *,'pattern = ',m1
        do j=ny,1,-1
            print '(i2,15f6.1)',j,(pattern(i,j,1,1+mod(m1-1,nperyear),yr01),i=1,min(nx,15))
        end do
        print *,'field = ',(yr1+yr2)/2,m1,nens1
        do j=ny,1,-1
            print '(i2,15f6.1)',j,(field(i,j,1,1+mod(m1-1,nperyear), &
                (yr1+yr2)/2,nens1),i=1,min(nx,15))
        end do
    end if

!   get mean values and sum these

    allocate(mean(nxf,nyf,nzf,nperyear))
    allocate(nn(nxf,nyf,nzf,nperyear))
    allocate(avemean(nperyear))
    allocate(pointhasdata(nxf,nyf,nzf))
    if ( .not. lnomissing ) then
        pointhasdata = .false. 
        call getensmean3(mean,nn,nxf,nyf,nzf,nperyear,field,nxf,nyf &
            ,nzf,nperyear,yr1,yr2,nens1,nens2,nx,ny,nz,yr1,1 &
            ,nperyear*(yr2-yr1+1), .false. )
        do m=m1,m2
            mo = m
            if ( mo > nperyear ) mo = mo - nperyear
            s = 0
            w = 0
            wt = 0
            do k=1,nz
                do j=1,ny
                    do i=1,nx
                        if ( pattern(i,j,k,month,yr01) < 1e33 ) then
                            wt = wt + wx(i)*wy(j)*wz(k)
                            if ( mean(i,j,k,mo) < 1e33 ) then
                                pointhasdata(i,j,k) = .true. 
                                s = s + wx(i)*wy(j)*wz(k)*mean(i,j,k,mo)* &
                                    pattern(i,j,k,month,yr01)
                                w = w + wx(i)*wy(j)*wz(k)
                            end if
                        end if
                    end do
                end do
            end do
            if ( w > minfac*wt ) then
                avemean(mo) = s/w
                if ( lwrite ) print *,'avemean(',mo,') = ',avemean(mo)
            else
                if ( lwrite ) print *,'w,minfac*wt = ',w,minfac*wt
                avemean(mo) = 3e33
            end if
        end do
    else
        avemean = 0
        pointhasdata = .true. 
    end if

!   sum

    do iens=nens1,nens2
        do y=yr1-1,yr2
            call keepalive(y-yr1+1,yr2-yr1+1)
            do m=m1,m2
                mo = m
                call normon(mo,y,yr,nperyear)
                if ( yr < yr1 .or. yr > yr2 ) cycle
                s = 0
                w = 0
                wt = 0
                do k=1,nz
                    do j=1,ny
                        do i=1,nx
                            if ( pointhasdata(i,j,k) .and. &
                                 pattern(i,j,k,month,yr01) < 1e33 ) then
                                wt = wt + wx(i)*wy(j)*wz(k)
                                if ( field(i,j,k,mo,yr,iens) < 1e33 ) then
                                    s = s+wx(i)*wy(j)*wz(k)* &
                                        (field(i,j,k,mo,yr,iens) - mean(i,j,k,mo)) &
                                        *pattern(i,j,k,month,yr01)
                                    w = w + wx(i)*wy(j)*wz(k)
                                end if
                            end if
                        end do
                    end do
                end do
                if ( w > minfac*wt .and. avemean(mo) < 1e33 ) then
                    var(mo,yr,iens) = s/w + avemean(mo)
                    if ( lwrite ) then
                        print *,'var(',mo,yr,iens,') = ',var(mo,yr,iens)
                        print *,'s,w,wt,avemean(mo) = ',s,w,wt,avemean(mo)
                    endif
                else
                    if ( lwrite ) then
                        print *,'avemean(',mo,') = ',avemean(mo)
                        print *,'w,minfac,wt = ',w,minfac,wt
                    endif
                    var(mo,yr,iens) = 3e33
                end if
            end do
        enddo
    enddo
    deallocate(wx)
    deallocate(wy)
    deallocate(mean)
    deallocate(nn)
    deallocate(avemean)
    deallocate(pointhasdata)
    if ( .false. ) then
        print *,'var = '
        do iens=nens1,nens2
            print *,'iens = ',iens
            do yr=yr1,yr2
                print *,yr,(var(mo,yr,iens),mo=1,nperyear)
            end do
        end do
    end if
end subroutine project3
