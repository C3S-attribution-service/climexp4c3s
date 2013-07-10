        subroutine getmean(mean,nn,nxmax,nymax,npermax,field,nxf,nyf
     +       ,nperyear,yrbeg,yrend,nx,ny,yrbegin,mobegin,nt,lwrite)
        implicit none
        integer nxmax,nymax,npermax,nxf,nyf,nperyear,yrbeg,yrend,nx,ny
     +       ,yrbegin,mobegin,nt
        integer nn(nxmax,nymax,npermax)
        real mean(nxmax,nymax,npermax),
     +        field(nxf,nyf,nperyear,yrbeg:yrend)
        logical lwrite
        call getensmean(mean,nn,nxmax,nymax,npermax,field,nxf,nyf
     +       ,nperyear,yrbeg,yrend,0,0,nx,ny,yrbegin,mobegin,nt,lwrite)
        end
        subroutine getensmean(mean,nn,nxmax,nymax,npermax,field,nxf,nyf
     +       ,nperyear,yrbeg,yrend,nens1,nens2,nx,ny,yrbegin,mobegin,nt
     +       ,lwrite)
*
*       compute means
*       
        implicit none
        integer nxmax,nymax,npermax,nxf,nyf,nperyear,yrbeg,yrend,nens1
     +       ,nens2,nx,ny,yrbegin,mobegin,nt
        integer nn(nxmax,nymax,npermax)
        real mean(nxmax,nymax,npermax),
     +        field(nxf,nyf,nperyear,yrbeg:yrend,0:nens2)
        logical lwrite
        integer month,year,i,j,k,iens
*
        if ( lwrite ) then
            print *,'getensmean: called with'
            print *,'  nx,nxf,nxmax     = ',nx,nxf,nxmax
            print *,'  ny,nyf,nymax     = ',ny,nyf,nymax
            print *,'  nperyear,npermax = ',nperyear,npermax
            print *,'  yrbeg,yrend      = ',yrbeg,yrend
            print *,'  nens1,nens2      = ',nens1,nens2
            print *,'  yrbegin,mobegin,nt ',yrbegin,mobegin,nt
        endif
        do month=1,nperyear
            do j=1,ny
                do i=1,nx
                    mean(i,j,month) = 0
                enddo
            enddo
        enddo
        do month=1,nperyear
            do j=1,ny
                do i=1,nx
                    nn(i,j,month) = 0
                enddo
            enddo
        enddo
        do iens=nens1,nens2
            month = mobegin-1
            year = yrbegin
            do k=1,nt
                month = month + 1
                if ( month.gt.nperyear ) then
                    month = month - nperyear
                    year = year + 1
                    call keepalive(year,yrbegin+nt/nperyear)
                endif
                do j=1,ny
                    do i=1,nx
                        if ( field(i,j,month,year,iens).lt.1e33 ) then
                            mean(i,j,month) = mean(i,j,month) + 
     +                           field(i,j,month,year,iens)
                            nn(i,j,month) = nn(i,j,month) + 1
                        endif
                    enddo
                enddo
            enddo
        enddo
        do month=1,nperyear
            do j=1,ny
                do i=1,nx
                    if ( nn(i,j,month).gt.0 ) then
                        mean(i,j,month) = mean(i,j,month)/nn(i,j,month)
                    else
                        mean(i,j,month) = 3e33
                    endif
                enddo
            enddo
        enddo
        if ( lwrite ) then
            do month=1,nperyear
                print *,'month = ',month
                do j=ny,1,-1
                    print '(i3,1000f8.2)',j,(mean(i,j,month),i=1,nx)
                    print '(i3,1000i8)',j,(nn(i,j,month),i=1,nx)
                enddo
            enddo
        endif
        end
