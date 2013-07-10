        subroutine sumit(data,npermax,nperyear,yrbeg,yrend,lsum,oper)
*
*       sum lsum consecutive months into the first one
*
        implicit none
        integer npermax,nperyear,yrbeg,yrend,lsum
        real data(npermax,yrbeg:yrend)
        character*1 oper
        integer i,j,j0,j1,mo,m,n,ksum
        logical lwrite
        parameter (lwrite=.false.)
        integer leap
        integer,save :: dpm(12,2)
        data dpm 
     +       /31,28,31,30,31,30,31,31,30,31,30,31
     +       ,31,29,31,30,31,30,31,31,30,31,30,31/
*
        if ( lsum.lt.1 ) return
        if ( lwrite ) then
            print *,'sumit: npermax,nperyear = ',npermax,nperyear
            print *,'       yrbeg,yrend      = ',yrbeg,yrend
            print *,'       lsum,oper        = ',lsum,oper
            print *,'       data(1990)       = ',
     +           (data(i,1990),i=1,nperyear)
            print *,'       data(1991)       = ',
     +           (data(i,1991),i=1,nperyear)
        endif
        if ( nperyear.le.350 ) then
            do i=yrbeg,yrend
                do j=1,nperyear
                    call sumone(data,npermax,nperyear,yrbeg,yrend,i,j
     +                   ,lsum,oper)
                enddo
            enddo
        else
*           daily data
            do i=yrbeg,yrend
                j0 = 1
                do mo=1,12
                    j1 = j0 - 1 + dpm(mo,leap(i))
                    if ( lsum.lt.30 ) then
*                       compute lsum-day averages (...) every lsum days
                        do j=j0,j1,lsum
                            call sumone(data,npermax,nperyear,yrbeg
     +                           ,yrend,i,j,lsum,oper)
*                           put the intermediate values to undefined
                            do m=j+1,min(j1,j+lsum-1)
                                data(j,i) = 3e33
                            enddo
                        enddo
                    elseif (lsum.lt.360 ) then
*                       sum months or longer
*                       I interprete 30 as monthly, 360-366 as annual
                        ksum = 0
                        do n=1,lsum/30
                            m = 1 + mod(n,12)
                            ksum = ksum + dpm(m,leap(i))
                        enddo
                        call sumone(data,npermax,nperyear,yrbeg
     +                       ,yrend,i,j,ksum,oper)                        
                    else
*                       sum years or longer
                        ksum = 0
                        do n=1,nint(lsum/365.24)
                            ksum = ksum + 364 + leap(i-1+n)
                        enddo
                        call sumone(data,npermax,nperyear,yrbeg
     +                       ,yrend,i,j,ksum,oper)
                    endif
                enddo
            enddo
        endif
        if ( lwrite ) then
            print *,'sumit: after averaging'
            print *,'       data(1990)       = ',
     +           (data(i,1990),i=1,nperyear)
            print *,'       data(1991)       = ',
     +           (data(i,1991),i=1,nperyear)
        endif
*
        end

        subroutine sumone(data,npermax,nperyear,yrbeg,yrend,i,j,lsum
     +       ,oper)
        implicit none
        integer npermax,nperyear,yrbeg,yrend,i,j,lsum
        character oper*1
        real data(npermax,yrbeg:yrend)
        integer n,m,ii
        real s
        s = data(j,i)
        if ( data(j,i).lt.1e33 ) then
            do n=1,lsum-1
                m = j+n
                call normon(m,i,ii,nperyear)
                if ( ii.gt.yrend ) then
                    s = 3e33
                elseif ( s.lt.1e33 .and. data(m,ii).lt.1e33 ) then
                    if ( oper.eq.'+' .or. oper.eq.'v' ) then
                        s = s + data(m,ii)
                    elseif ( oper.eq.'a' ) then
                        s = max(s,data(m,ii))
                    elseif ( oper.eq.'i' ) then
                        s = min(s,data(m,ii))
                    else
                        print *,'sumit: error: cannot handle oper = '
     +                       ,oper,', only [+vai]'
                        call abort
                    endif
                else
                    s = 3e33
                endif
            enddo
            if ( oper.eq.'v' .and. s.lt.1e33 )
     +           then
                s = s/lsum
            endif
            data(j,i) = s
        endif
        end
