        subroutine myloess(data,npermax,nperyear,yrbeg,yrend,nmonth
     +           ,minfac,filtertype,hilo,yearmonth,family,lwrite)
!
!       my own LOESS routine; I cannot get the horrible netlib ones to work
!       input:
!       data(1:nperyear,yrbeg:yrend)    contains the data with first dimension days, months, ...
!       nmonth                          use data from point-nmonth to point+nmongth (ie 2*nmonth+1 length)
!       minfac                          only consider intervals with minfac% valid data
!       filtertype                      loess1: 1st degree, loess2: 2nd degree
!       hilo                            hi: high-pass filtering, lo: low-pass filtering
!       yearmonth                       year: filter in the year direction for each month (day)
!                                       month: filter consecutive months (days)
!       family                          should be 'gaus*'
!       lwrite                          debug flag
!
!       output:
!       data(1:nperyear,yrbeg:yrend)    contains the data with first dimension days, months, ...
!
        implicit none
        integer npermax,nperyear,yrbeg,yrend,nmonth
        real data(npermax,yrbeg:yrend),minfac
        character filtertype*(*),hilo*(*),yearmonth*(*),family*(*)
        logical lwrite
        integer i,j,n,d,degree,yr,mo,yy,mm
        real,allocatable :: x(:),y(:),sig(:),filtered(:,:)

        if ( lwrite ) then
            print *,'myloess: nmonth,minfac = ',nmonth,minfac
            print *,'         filtertype,hilo = ',filtertype,hilo
            print *,'         yearmonth,family = ',yearmonth,family
            print *,'         npermax,nperyear = ',npermax,nperyear
        end if

        if ( filtertype.eq.'loess1' ) then
            degree = 1
        else if ( filtertype.eq.'loess2' ) then
            degree = 2
        else
            write(0,*) 'myloess: error: unknown filtertype '
     +           ,trim(filtertype)
            call abort
        end if
        if ( family(1:4).ne.'gaus' ) then
            write(0,*) 'myloess: error: family ',trim(family),
     +           ' not yet ready'
            call abort
        end if

        allocate(x(2*nmonth+1),y(2*nmonth+1),sig(2*nmonth+1))
        allocate(filtered(nperyear,yrbeg:yrend))

        if ( yearmonth.eq.'month' ) then
            i = 0
            do yr=yrbeg,yrend
                do mo=1,nperyear
                    n = 0
                    do j=-nmonth,nmonth
                        mm = mo + j
                        call normon(mm,yr,yy,nperyear)
                        if ( yy.lt.yrbeg .or. yy.gt.yrend ) cycle
                        if ( data(mm,yy).lt.1e33 ) then
                            n = n + 1
                            x(n) = j
                            y(n) = data(mm,yy)
                            sig(n) = 1/sqrt(
     +                           (1-(real(abs(j))/real(nmonth+1))**3)**3
     +                           )
                        end if
                    end do
                    call loessfit(filtered(mo,yr),x,y,sig,n,minfac
     +                   ,nmonth,degree,lwrite)
                end do
            end do
        else if ( yearmonth.eq.'year' ) then
            do mo=1,nperyear
                i = 0
                do yr=yrbeg,yrend
                    n = 0
                    do j=-nmonth,nmonth
                        yy = yr + j
                        if ( yy.lt.yrbeg .or. yy.gt.yrend ) cycle
                        if ( data(mo,yy).lt.1e33 ) then
                            n = n + 1
                            x(n) = j
                            y(n) = data(mo,yy)
                            sig(n) = 1/sqrt(
     +                           (1-(real(abs(j))/real(nmonth+1))**3)**3
     +                           )
                        end if
                    end do
                    call loessfit(filtered(mo,yr),x,y,sig,n,minfac
     +                   ,nmonth,degree,lwrite)
                end do
             end do
        else
            write(0,*) 'myloess: error: unknown value fopr yearmonth '
     +           ,trim(yearmonth)
            call abort
        end if

        if ( hilo(1:2).eq.'lo' ) then
            do mo=1,nperyear
                do yr=yrbeg,yrend
                    data(mo,yr) = filtered(mo,yr)
                end do
            end do
        else if ( hilo(1:2).eq.'hi' ) then
            do mo=1,nperyear
                do yr=yrbeg,yrend
                    if ( data(mo,yr).lt.1e33 .and. 
     +                   filtered(mo,yr).lt.1e33 ) then
                        data(mo,yr) = data(mo,yr) - filtered(mo,yr)
                    else
                        data(mo,yr) = 3e33
                    end if
                end do
            end do
        else
            write(0,*) 'myloess: error: unknown value for hilo: '
     +           ,trim(hilo)
            call abort
        end if

        deallocate(x,y,sig)
        deallocate(filtered)

        end subroutine

        subroutine loessfit(filtered,x,y,sig,n,minfac,nmonth,degree
     +       ,lwrite)
        implicit none
        integer n,nmonth,degree
        real filtered,x(n),y(n),sig(n),minfac
        logical lwrite
        integer i
        real a,b,c,siga,sigb,sigc,chi2,q,r,prob,df

        if ( n.lt.minfac*(2*nmonth+1) .or. n.le.nmonth .or.
     +       n.lt.2+degree ) then
            filtered = 3e33
        else if ( degree.eq.1 ) then
            if ( lwrite ) then
                print *,'myloess: x,y,sig = ',n
                do i=1,n
                    print *,x(i),y(i),sig(i)
                end do
            end if
            call fit(x,y,n,sig,1,a,b,siga,sigb,chi2,q)
            filtered = a
        else if ( degree.eq.2 ) then
            df = n - 3          ! not needed
            call fit2(x,y,n,a,b,c,siga,sigb,sigc,r,prob,df,0,lwrite)
            filtered = a            
        else
            write(0,*) 'myloess: error: degree ',degree,
     +           ' not yet ready'
            call abort
        end if
        end subroutine
