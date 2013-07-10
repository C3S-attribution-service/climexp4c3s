        program variations
*
*       Get the standard deviation of an index, smoothened 
*       over n-year intervals
*
*       GJvO hack, dec-1997
*
        implicit none
        integer yrbeg,yrend
        parameter (yrbeg=1740,yrend=2020)
        integer i,j,k,l,nn(-1:12),nyear,year,month,j1,j2
        real data(12,yrbeg:yrend),s1(-1:12),s2(-1:12),adata,a
        logical lwrite,txtdat
        character line*128
        integer iargc
        lwrite = .FALSE.
*
*       check arguments
*
        if ( iargc().ne.2 ) then
            print *,'usage: variate file n'
            stop
        endif
        call getarg(2,line)
        read(line,*) nyear
        print '(a,i3,a)','Averaging over ',nyear,' years'
        call getarg(1,line)
*
*       init
*
        do j=yrbeg,yrend
            do i=1,12
                data(i,j) = 3e33
            enddo
        enddo
*
*       read data from station file downloaded from
*       http://www.ncdc.noaa.gov/ghcn/ghcnV1.CLIMVIS.html
*
        call rdtxtdat(line,data,yrbeg,yrend)
*
*       compute summed variance
*
        print '(a)','              mean     s.d.'
        print '(a)',
     +    '             12-mnth ave'//
     +      '     All year        Jan             Feb'//
     +      '             Mar             Apr             May'//
     +      '             Jun             Jul             Aug'//
     +      '             Sep             Oct             Nov'//
     +      '             Dec'
        do year=yrbeg,yrend-nyear+1
            do month=-1,12
                if ( month.le.0 ) then
                    j1 = 1
                    j2 = 12
                else
                    if ( nyear.lt.1 ) stop
                    j1 = month
                    j2 = month
                endif
                nn(month) = 0
                s1(month) = 0
                s2(month) = 0
                do k=0,nyear-1
                    do j=j1,j2
                        if ( month.ge.0 ) then
                            adata = data(j,year+k)
                        else
                            adata = 0
                            do l=-6,6
                                if ( j+l.lt.0 ) then
                                    if ( year+k.gt.yrbeg ) then
                                        if ( abs(l).eq.6 ) then
                                            a = data(j+l+12,year-1+k)/2
                                        else
                                            a = data(j+l+12,year-1+k)
                                        endif
                                    else
                                        a = 3e33
                                    endif
                                elseif ( j+l.lt.12 ) then
                                    if ( abs(l).eq.6 ) then
                                        a = data(j+l,year+k)/2
                                    else
                                        a = data(j+l,year+k)
                                    endif
                                else
                                    if ( year+k.lt.yrend ) then
                                        if ( abs(l).eq.6 ) then
                                            a = data(j+l-12,year+1+k)/2
                                        else
                                            a = data(j+l-12,year+1+k)
                                        endif
                                    else
                                        a = 3e33
                                    endif
                                endif
                                if ( a.lt.1e33 .and. adata.lt.1e33 )
     +                              then
                                    adata = adata + a
                                else
                                    adata = 3e33
                                endif
                            enddo
                            if ( adata.lt.1e33 ) then
                                adata = adata/12
                            endif
                        endif
                        if ( adata.lt.1e33 ) then
                            nn(month) = nn(month) + 1
                            s1(month) = s1(month) + adata
                            s2(month) = s2(month) + adata**2
                        endif
                    enddo
                enddo
                if ( nn(month).gt.0 ) then
                    s1(month) = s1(month)/nn(month)
                    s2(month) = s2(month)/nn(month)
                    if ( s2(month).lt.0 ) then
                        print *,'Error: var<0 ',s2(month)
                        s2(month) = 0
                    endif
                    s2(month) = sqrt(s2(month))
                else
                    s1(month) = -999.99
                    s2(month) = -999.99
                endif
            enddo
            if ( nn(0).gt.11*nyear ) then
                print '(2i5,28f8.2)',year,year+nyear-1,
     +              (s1(j),s2(j),j=-1,12)
            endif
        enddo
*
*       error messages
*
        stop
  900   print *,'error reading station data ',line
        if ( txtdat ) then
            print*,'expecting: <lf>yyyy<lf>text<lf>12*(month value<lf>)'
        else
            print *,'expecting: yyyy val1 ...val12'
        endif
        stop
        end
