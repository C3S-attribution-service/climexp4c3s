subroutine filllinarray(dindx,ddata,lfirst,dddata,ndata,n,j1,j2 &
    ,lag,k,nperyear,imens,indxmx,indx,data,npermax,yrbeg,yrend &
    ,nensmax,yrstart,yrstop,yrmo)

!   fill linear arrays without absent values

    implicit none
    include 'getopts.inc'
    integer :: ndata,n,j1,j2,lag,k,nperyear,npermax,yrbeg,yrend,nensmax,indxmx, &
        yrstart,yrstop,yrmo(3,ndata)
    integer :: imens(0:indxmx)
    real :: dindx(ndata),ddata(ndata),dddata(ndata)
    logical :: lfirst(ndata)
    real :: indx(npermax,yrbeg:yrend,0:nensmax,indxmx),data(npermax,yrbeg:yrend,0:nensmax)
    integer :: yr,jj,j,i,m,ii,iiens,iens,jens,if,jm,jp,im,ip
    logical :: lastvalid

    if ( lwrite ) then
        print *,'filllinarray: nens1,nens2,imens(0),imens(k) = ', &
            nens1,nens2,imens(0),imens(k)
        print *,'              yr1,yr2 = ',yr1,yr2
        print *,'              j1,j2 = ',j1,j2
        if ( .false. ) then
            do iens=0,0
                print *,'data = ',iens
                call printdatfile(6,data(1,yrbeg,iens),npermax,nperyear,yrbeg,yrend)
            end do
            do iens=0,0
                print *,'indx = ',iens,k
                call printdatfile(6,indx(1,yrbeg,iens,k),npermax,nperyear,yrbeg,yrend)
            end do
        end if
    end if
    if ( nens2 > max(imens(0),imens(k)) ) then
        write(0,*) 'filllinarray: error: nens2 > mens: ',nens2,imens(0),imens(k)
        write(*,*) 'filllinarray: error: nens2 > mens: ',nens2,imens(0),imens(k)
        call exit(-1)
    end if

    n = 0
    do iiens=nens1,nens2
        lastvalid = .false. 
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
        if ( nfittime > 0 ) then
            call derivative(2*nfittime+1,data(1,yrbeg,iens),indx(1,yrbeg,iens,indxmx), &
                npermax,yrbeg,yrend,nperyear,minfac,lwrite)
        end if
        do yr=yr1-1,yr2
            if ( j1 /= j2 .and. j2-j1+1 /= nperyear ) then
                lastvalid = .false. 
            end if
            do jj=j1,j2
                if ( fix2 ) then
                    j = jj+lag
                else
                    j = jj
                end if
                call normon(j,yr,i,nperyear)
                if ( i < yr1 .or. i > yr2 ) goto 710
                m = j-lag
                call normon(m,i,ii,nperyear)
                if ( ii < yr1 .or. ii > yr2 ) goto 710
                if ( .false. ) then
                    print *,'data(',j,i,iens,')    = ',data(j,i,iens)
                    print *,'indx(',m,ii,jens,k,') = ',indx(m,ii,jens,k)
                end if
                if (  data(j,i,iens) < 1e33 .and. indx(m,ii,jens,k) < 1e33 .and. &
                    ( lconting .or. ( &
                    (data(j,i,iens) <= maxdata) .eqv. &
                    (data(j,i,iens) >= mindata) .eqv. &
                    (maxdata >= mindata) ) .and. ( &
                    (indx(m,ii,jens,k) <= maxindx) &
                    .eqv. &
                    (indx(m,ii,jens,k) >= minindx) &
                    .eqv.(maxindx >= minindx) ) ) ) &
                        then
                    if ( n == 0 .and. lwrite ) then
                        print *,'filllinarray: first valid point '
                    end if
                    if ( lwrite ) then
                        print '(i3,i5,i3,g12.4,i3,i5,i3,g12.4,i3)',j,i,iens,data(j,i,iens), &
                            m,ii,jens,indx(j,i,iens,k),k
                    end if
                    n = n+1
                    if ( n > ndata ) then
                        write(0,*) 'filllinarray: error: n>ndata ',ndata
                        write(*,*) 'filllinarray: error: n>ndata ',ndata
                        call exit(-1)
                    end if
                    ddata(n) = data(j,i,iens)
                    dindx(n) = indx(m,ii,jens,k)
                    lfirst(n) = .not. lastvalid
                    lastvalid = .true. 
                    if ( yrstop /= -999 ) then
                        yrmo(1,n) = i
                        yrmo(2,n) = j
                        yrmo(3,n) = iens
                        yrstart = min(yrstart,i,ii)
                        yrstop  = max(yrstop,i,ii)
                    end if
                    if ( lwrite .and. lfirst(n) ) print *,'boundary at ',n,j,i,iens
                    if ( nfittime > 0 ) then
                        if ( indx(j,i,iens,indxmx) > 1e33 ) then
                            n = n - 1
                            goto 710
                        end if
                        dddata(n) = indx(j,i,iens,indxmx)
                    end if   ! nfittime
                else
                    goto 710
                end if       ! valid data point
                goto 720
            710 continue
!               invalid data point
                lastvalid = .false. 
            720 continue
            end do           ! month jj
        end do               ! year yr
    end do                   ! iens
end subroutine filllinarray
