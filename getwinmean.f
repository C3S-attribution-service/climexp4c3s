        subroutine getwinmean(mean,nn,nxmax,nymax,npermax,field,nxf,nyf
     +        ,nperyear,yrbeg,yrend,nx,ny,yrbegin,mobegin,nt,x1,x2,y1,y2
     +        ,lwrite)
*
*       compute means
*       
        implicit none
        integer nxmax,nymax,npermax,nxf,nyf,nperyear,yrbeg,yrend,nx,ny
     +        ,yrbegin,mobegin,nt,x1,x2,y1,y2
        integer nn(nxmax,nymax,npermax)
        real mean(nxmax,nymax,npermax),
     +        field(nxf,nyf,nperyear,yrbeg:yrend)
        logical lwrite
        integer month,year,i,j,k,ii
        integer normx
        external normx
*
        if ( lwrite ) then
            print *,'getwinmean: called with'
            print *,'  nx,nxf,nxmax     = ',nx,nxf,nxmax
            print *,'  ny,nyf,nymax     = ',ny,nyf,nymax
            print *,'  x1,x2,y1,y2      = ',x1,x2,y1,y2
            print *,'  nperyear,npermax = ',nperyear,npermax
            print *,'  yrbeg,yrend      = ',yrbeg,yrend
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
        month = mobegin-1
        year = yrbegin
        do k=1,nt
            month = month + 1
            if ( month.gt.nperyear ) then
                month = month - nperyear
                year = year + 1
                call keepalive1('Computing mean of ',
     +               year,yrbegin+nt/nperyear)
            endif
            if ( lwrite ) print *,'month,year = ',month,year
            do j=y1,y2
                do i=x1,x2
                    ii = normx(i,nx)
                    if ( field(ii,j,month,year).lt.1e33 ) then
                        mean(ii,j,month) = mean(ii,j,month) + 
     +                        field(ii,j,month,year)
                        nn(ii,j,month) = nn(ii,j,month) + 1
                    endif
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
        if ( .false. .and. lwrite ) then
            do month=1,nperyear
                print *,'month = ',month
                do j=1,ny
                    print '(i3,1000f10.2)',j,(mean(i,j,month),i=1,nx)
                    print '(i3,1000i10)',j,(nn(i,j,month),i=1,nx)
                enddo
            enddo
        endif
        end


        subroutine getenswinmean(mean,nn,nxmax,nymax,nzmax,npermax,field
     +       ,nxf,nyf,nzf,nperyear,yrbeg,yrend,nens1,nens2,nx,ny,nz
     +       ,yrbegin,mobegin,nt,x1,x2,y1,y2,lwrite)
*
*       compute mean and variance of a field
*
        implicit none
        integer nxmax,nymax,nzmax,npermax,nxf,nyf,nzf,nperyear,yrbeg
     +       ,yrend,nens1,nens2,nx,ny,nz,yrbegin,mobegin,nt,x1,x2,y1,y2
        integer nn(nxmax,nymax,nzmax,npermax)
        real mean(nxmax,nymax,nzmax,npermax,2),
     +       sd(nxmax,nymax,nzmax,npermax),
     +       field(nxf,nyf,nzf,nperyear,yrbeg:yrend,nens1:nens2)
        logical lwrite
        integer month,year,i,j,k,it,ii,iens,moment
        integer normx
        external normx
*
        if ( lwrite ) then
            print *,'getenswinmean: called with'
            print *,'  nx,nxf,nxmax     = ',nx,nxf,nxmax
            print *,'  ny,nyf,nymax     = ',ny,nyf,nymax
            print *,'  nz,nzf,nzmax     = ',nz,nzf,nzmax
            print *,'  x1,x2,y1,y2      = ',x1,x2,y1,y2
            print *,'  nperyear,npermax = ',nperyear,npermax
            print *,'  yrbeg,yrend      = ',yrbeg,yrend
            print *,'  nens1,nens2      = ',nens1,nens2
            print *,'  yrbegin,mobegin,nt ',yrbegin,mobegin,nt
        endif
        do moment=1,2
            do month=1,nperyear
                do k=1,nz
                    do j=1,ny
                        do i=1,nx
                            mean(i,j,k,month,moment) = 0
                        end do
                    end do
                end do
            end do
        end do
!
!       mean
!
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
            call keepalive(iens-nens1,2*(nens2-nens1))
            month = mobegin-1
            year = yrbegin
            do it=1,nt
                month = month + 1
                if ( month.gt.nperyear ) then
                    month = month - nperyear
                    year = year + 1
                endif
                if ( lwrite ) print *,'month,year,iens = ',month,year
     +               ,iens
                do k=1,nz
                    do j=y1,y2
                        do i=x1,x2
                            ii = normx(i,nx)
                            if ( field(ii,j,k,month,year,iens).lt.1e33 )
     +                           then
                                mean(ii,j,k,month,1) =
     +                               mean(ii,j,k,month,1) + 
     +                               field(ii,j,k,month,year,iens)
                                nn(ii,j,k,month) = nn(ii,j,k,month) + 1
                            end if
                        end do
                    end do
                end do
            enddo
        enddo
        do month=1,nperyear
            do k=1,nz
                do j=1,ny
                    do i=1,nx
                        if ( nn(i,j,k,month).gt.0 ) then
                            mean(i,j,k,month,1) =
     +                           mean(i,j,k,month,1)/nn(i,j,k,month)
                        else
                            mean(i,j,k,month,1) = 3e33
                        end if
                    end do
                end do
            end do
        end do
!
!       variance
!
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
            call keepalive(iens-nens1+(nens2-nens1),2*(nens2-nens1))
            month = mobegin-1
            year = yrbegin
            do it=1,nt
                month = month + 1
                if ( month.gt.nperyear ) then
                    month = month - nperyear
                    year = year + 1
                endif
                if ( lwrite ) print *,'month,year,iens = ',month,year
     +               ,iens
                do k=1,nz
                    do j=y1,y2
                        do i=x1,x2
                            ii = normx(i,nx)
                            if ( field(ii,j,k,month,year,iens).lt.1e33 )
     +                           then
                                mean(ii,j,k,month,2) =
     +                               mean(ii,j,k,month,2)
     +                               + (field(ii,j,k,month,year,iens)-
     +                               mean(ii,j,k,month,1))**2
                                nn(ii,j,k,month) = nn(ii,j,k,month) + 1
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
                        if ( nn(i,j,k,month).gt.0 ) then
                            mean(i,j,k,month,2) =
     +                           mean(i,j,k,month,2)/nn(i,j,k,month)
                        else
                            mean(i,j,k,month,2) = 3e33
                        end if
                    end do
                end do
            end do
        end do
!
!       debug info
!
        if ( .false. .and. lwrite ) then
            do month=1,nperyear
                print *,'month = ',month
                do k=1,nz
                    do j=1,ny
                        print '(2i3,1000f10.2)',j,k,(mean(i,j,k,month,1)
     +                       ,i=1,nx)
                        print '(2i3,1000f10.2)',j,k,(mean(i,j,k,month,2)
     +                       ,i=1,nx)
                        print '(2i3,1000i10)',j,k,(nn(i,j,k,month),i=1
     +                       ,nx)
                    end do
                end do
            end do
        end if
        end

        integer function normx(i,nx)
        implicit none
        integer i,nx
        if ( i.le.0 ) then
            normx = i + (-i/nx+1)*nx
        elseif ( i.gt.nx ) then
            normx = 1 + mod(i-1,nx)
        else
            normx = i
        endif
        end
