subroutine getmean(mean,nn,nxmax,nymax,npermax,field,nxf,nyf &
    ,nperyear,yrbeg,yrend,nx,ny,yrbegin,mobegin,nt,lwrite)
    implicit none
    integer,intent(in) :: nxmax,nymax,npermax,nxf,nyf,nperyear,yrbeg,yrend,nx,ny,yrbegin,mobegin,nt
    integer,intent(out) :: nn(nxmax,nymax,npermax)
    real,intent(in) :: field(nxf,nyf,nperyear,yrbeg:yrend)
    real,intent(out) :: mean(nxmax,nymax,npermax)
    logical,intent(in) :: lwrite
    call getensmean(mean,nn,nxmax,nymax,npermax,field,nxf,nyf &
        ,nperyear,yrbeg,yrend,0,0,nx,ny,yrbegin,mobegin,nt,lwrite)
end subroutine getmean

subroutine getensmean(mean,nn,nxmax,nymax,npermax,field,nxf,nyf &
    ,nperyear,yrbeg,yrend,nens1,nens2,nx,ny,yrbegin,mobegin,nt,lwrite)
    implicit none
    integer,intent(in) :: nxmax,nymax,npermax,nxf,nyf,nperyear,yrbeg,yrend,nens1 &
        ,nens2,nx,ny,yrbegin,mobegin,nt
    integer,intent(out) :: nn(nxmax,nymax,npermax)
    real,intent(in) :: field(nxf,nyf,nperyear,yrbeg:yrend,0:nens2)
    real,intent(out) :: mean(nxmax,nymax,npermax)
    logical,intent(in) :: lwrite
    call getensmean3(mean,nn,nxmax,nymax,1,npermax,field,nxf,nyf,1 &
        ,nperyear,yrbeg,yrend,nens1,nens2,nx,ny,1,yrbegin,mobegin &
        ,nt,lwrite)
end subroutine getensmean

subroutine getensmean3(mean,nn,nxmax,nymax,nzmax,npermax,field &
    ,nxf,nyf,nzf,nperyear,yrbeg,yrend,nens1,nens2,nx,ny,nz &
    ,yrbegin,mobegin,nt,lwrite)

!   compute means

    implicit none
    integer,intent(in) :: nxmax,nymax,nzmax,npermax,nxf,nyf,nzf,nperyear,yrbeg &
        ,yrend,nens1,nens2,nx,ny,nz,yrbegin,mobegin,nt
    integer,intent(out) :: nn(nxmax,nymax,nzmax,npermax)
    real,intent(in) :: field(nxf,nyf,nzf,nperyear,yrbeg:yrend,0:nens2)
    real,intent(out) :: mean(nxmax,nymax,nzmax,npermax)
    logical,intent(in) :: lwrite
    integer :: month,year,i,j,k,m,iens

    if ( lwrite ) then
        print *,'getensmean3: called with'
        print *,'  nx,nxf,nxmax     = ',nx,nxf,nxmax
        print *,'  ny,nyf,nymax     = ',ny,nyf,nymax
        print *,'  nz,nzf,nzmax     = ',nz,nzf,nzmax
        print *,'  nperyear,npermax = ',nperyear,npermax
        print *,'  yrbeg,yrend      = ',yrbeg,yrend
        print *,'  nens1,nens2      = ',nens1,nens2
        print *,'  yrbegin,mobegin,nt ',yrbegin,mobegin,nt
    end if
    do month=1,nperyear
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    mean(i,j,k,month) = 0
                end do
            end do
        end do
    end do
    do month=1,nperyear
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    nn(i,j,k,month) = 0
                end do
            end do
        end do
    end do
    do iens=nens1,nens2
        month = mobegin-1
        year = yrbegin
        do m=1,nt
            month = month + 1
            if ( month > nperyear ) then
                month = month - nperyear
                year = year + 1
                call keepalive1('Computing mean of year',year,yrbegin+nt/nperyear)
            endif
            do k=1,nz
                do j=1,ny
                    do i=1,nx
                        if ( field(i,j,k,month,year,iens) < 1e33 ) then
                            mean(i,j,k,month) = mean(i,j,k,month) + field(i,j,k,month,year,iens)
                            nn(i,j,k,month) = nn(i,j,k,month) + 1
                        end if
                    end do
                end do
            end do
        end do
    end do
    do month=1,nperyear
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    if ( nn(i,j,k,month) > 0 ) then
                        mean(i,j,k,month) = mean(i,j,k,month)/nn(i,j,k,month)
                    else
                        mean(i,j,k,month) = 3e33
                    end if
                end do
            end do
        end do
    end do
    if ( lwrite ) then
        do month=1,nperyear
            print *,'month = ',month
            do k=1,nz
                print '(a,i4)','level ',k
                do j=ny,1,-1
                    print '(i3,1000f8.2)',j,(mean(i,j,k,month),i=1,nx)
                    print '(i3,1000i8)',j,(nn(i,j,k,month),i=1,nx)
                end do
            end do
        end do
    end if
end subroutine getensmean3
