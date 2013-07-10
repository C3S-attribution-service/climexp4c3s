        subroutine getspectrum(data,npermax,nperyear,yrbeg,yrend,yr1,yr2
     +       ,mens1,mens,nens1,nens2,j1,j2,epx,epy,ndata,nout,ofac
     +       ,yrstart,yrstop)
*
*       compute a spectrum from the time series in data
*
        implicit none
        integer npermax,nperyear,yrbeg,yrend,yr1,yr2,mens1,mens,nens1
     +       ,nens2,ndata,nout,yrstart,yrstop
        real data(npermax,yrbeg:yrend,0:nens2)
        real epx(ndata),epy(ndata),ofac
        integer iens,n,n1,n2,yr,mo,i,j
        integer,allocatable :: nepx(:)
        real hifac,sx,sy
        real,allocatable :: x(:),y(:),px(:),py(:)
        allocate(x(ndata))
        allocate(y(ndata))
        allocate(px(4*ndata))
        allocate(px(4*ndata))
        allocate(nepx(4*ndata))
*
*       loop over ensemble members
*       
        do iens=nens1,nens2
*
*           fill arrays
*       
            n = 0
            n1 = 0
            do yr=yr1-1,yr2
                do mo=j1,j2
                    j = mo
                    call normon(j,yr,i,nperyear)
                    if ( yr.lt.yr1 .or. yr.gt.yr2 ) cycle
                    if ( abs(data(j,i,iens)).lt.1e33 ) then
                        n = n + 1
                        yrstart = min(yrstart,i)
                        yrstop  = max(yrstop,i)
                        x(n) = i + (j-0.5)/nperyear
                        if ( nperyear.gt.12 ) then
                            x(n) = x(n)*365.24
                        endif
                        y(n) = data(j,i,iens)
                        if ( n1.eq.0 ) then
                            n1 = nperyear*i+j
                        endif
                        n2 = nperyear*i+j
                    endif
                enddo
            enddo
*       
*           compute/estimate other parameters
*       
            if ( m1.eq.0 ) then
                hifac = real(n)/real(n2-n1+1)
            else
                hifac = real(n)/real((n2-n1)/nperyear+1)
            endif
            if ( hifac.eq.1 ) then
                ofac = 1
            else
                ofac = 4        ! see how it works
            endif
*       
*           call period (Numerical recipes p 572)
*           take care of dependent data!
*
            call period(x,y,n,ofac,hifac,px,py,4*ndata,nout,jmax,prob)
*
*           average
*       
            if ( avex.gt.1 ) then
                do i=nint(ofac),nout-avex,avex
                    sx = px(i)
                    sy = py(i)
                    do j=1,avex-1
                        sx = sx + px(i+j)
                        sy = sy + py(i+j)
                    enddo
                    n = 1+i/avex
                    px(n) = sx/avex
                    py(n) = sy/avex
                enddo
                ofac = 1+nint(ofac)/avex
                nout = n
            endif
*
*           collect ensemble informnation
*
            if ( iens.eq.nens1 ) then
                nn = nout
                do i=1,nn
                    nepx(i) = 1
                    epx(i) = px(i)
                    epy(i) = py(i)
                    if ( lwrite ) print *,i,epx(i),epy(i),nepx(i)
                enddo
            else
                do i=1,nout
                    if ( epx(i).eq.px(i) ) then
                        nepx(i) = nepx(i) + 1
                        epy(i) = epy(i) + py(i)
                    else        ! unequal array sizes - choose nearest
                        do j=1,nn-1
                            if ( px(j).lt.(epx(j)+epx(j+1))/2 ) then
                                goto 100
                            endif
                        enddo
  100                   continue
                        epx(j) = (nepx(j)*epx(j) + px(i))/(nepx(j) + 1)
                        nepx(j) = nepx(j) + 1
                        epy(j) = epy(j) + py(i)
                    endif
                enddo
            endif
        enddo
        do i=1,nn
            epy(i) = epy(i)/nepx(i)
        enddo
        end
