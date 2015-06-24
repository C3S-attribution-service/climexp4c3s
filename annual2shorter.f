        subroutine annual2shorter(data,npermax,yrbeg,yrend,nperyear,
     +       newdata,npermaxnew,yrbegnew,yrendnew,npernew,m1,n,nfac,
     +       lwrite)
!
!       distribute (nfac>1) or set equal (nfac=1) the values in the lower-frequency series
!       data into the higher-frequency series newdata in n months/days starting at month/day m1
!
        implicit none
        integer npermax,yrbeg,yrend,nperyear,
     +       npermaxnew,yrbegnew,yrendnew,npernew,
     +       m1,n,nfac
        real data(npermax,yrbeg:yrend),
     +       newdata(npermaxnew,yrbegnew:yrendnew)
        logical lwrite
        integer yr,mo,dy,i,k,yr1,yr2,dpm(12,2)
        integer,external :: leap
        data dpm /31,28,31,30,31,30,31,31,30,31,30,31,
     +            31,29,31,30,31,30,31,31,30,31,30,31/

        yr1 = max(yrbeg,yrbegnew)
        yr2 = min(yrend,yrendnew)
        newdata = 3e33
        if ( nperyear.eq.1 ) then
            if ( lwrite ) print *,'annual2shorter: annual data'
            do yr=yr1,yr2
                do k=1,n
                    mo = m1 + k - 1
                    call normon(mo,yr,i,npernew)
                    if ( i.ge.yrbegnew .and. i.le.yrendnew ) then
                        if ( data(1,yr).lt.1e33 ) then
                            newdata(mo,i) = data(1,yr)/nfac
                        endif
                    endif
                enddo
            enddo
        elseif ( mod(npernew,nperyear).eq.0 ) then
            if ( nfac.ne.1 ) nfac = npernew/nperyear
            do yr=yr1,yr2
                do mo=1,npernew
                    k = 1 + (mo-1)*nperyear/npernew
                    if ( data(k,yr).lt.1e33 ) then
                        newdata(mo,yr) = data(k,yr)/nfac
                    endif
                enddo
            enddo
        elseif ( npernew.eq.366 .and. nperyear.eq.12 ) then
            do yr=yr1,yr2
                k = 0
                do mo=1,12
                    if ( nfac.ne.1 ) nfac = dpm(mo,leap(yr))
                    do dy=1,dpm(mo,2)
                        k = k + 1
                        if ( k.eq.60 .and. leap(yr).eq.1 ) then
                            newdata(k,yr) = 3e33
                        elseif ( data(mo,yr).lt.1e33 ) then
                            newdata(k,yr) = data(mo,yr)/nfac
                        endif
                    enddo
                enddo
            enddo
        else
            write(0,*) 'annual2longer: error: cannet handle ',
     +           'conversion from ',nperyear,' to ',npernew ,
     +           ' points per year yet'
            write(*,*) 'annual2longer: error: cannet handle ',
     +           'conversion from ',nperyear,' to ',npernew ,
     +           ' points per year yet'
        endif
        end
