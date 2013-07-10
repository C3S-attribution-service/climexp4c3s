        subroutine sumit(data,npermax,nperyear,yrbeg,yrend,lsum,oper)
*
*       sum lsum consecutive months into the first one
*
        implicit none
        integer npermax,nperyear,yrbeg,yrend,lsum
        real data(npermax,yrbeg:yrend)
        character*1 oper
        integer i,j,m,n,ii
        real s,absent
        parameter (absent=3e33)
        logical lwrite
        parameter (lwrite=.false.)
        integer,external :: leap
*
        if ( lsum.gt.1 ) then
            if ( lwrite ) then
                print *,'sumit: npermax,nperyear = ',npermax,nperyear
                print *,'       yrbeg,yrend      = ',yrbeg,yrend
                print *,'       lsum,oper        = ',lsum,oper
                print *,'       data(1990)       = ',
     +                (data(i,1990),i=1,nperyear)
                print *,'       data(1991)       = ',
     +                (data(i,1991),i=1,nperyear)
            endif
            do i=yrbeg,yrend
                do j=1,nperyear
                    s = data(j,i)
                    if ( data(j,i).lt.0.9*absent ) then
                        do n=1,lsum-1
                            m = j+n
                            if ( nperyear.eq.366 .and. leap(i).eq.1
     +                           .and. j.lt.60 .and. m.ge.60 ) m = m + 1
                            call normon(m,i,ii,nperyear)
                            if ( ii.gt.yrend ) then
                                s = absent
                            elseif ( s.lt.0.9*absent .and. 
     +                                data(m,ii).lt.0.9*absent ) then
                                if ( oper.eq.'+' .or. oper.eq.'v' ) then
                                    s = s + data(m,ii)
                                elseif ( oper.eq.'a' ) then
                                    s = max(s,data(m,ii))
                                elseif ( oper.eq.'i' ) then
                                    s = min(s,data(m,ii))
                                else
                                    print *,'sumit: error: cannot '
     +                                    ,'handle oper = ',oper
     +                                    ,', only [+vai]'
                                    call abort
                                endif
                            else
                                s = absent
                            endif
                        enddo
                        if ( oper.eq.'v' .and. s.lt.0.9*absent ) then
                            s = s/lsum
                        endif
                        data(j,i) = s
                    endif
                enddo
            enddo
            if ( lwrite ) then
                print *,'sumit: after averaging'
                print *,'       data(1990)       = ',
     +                (data(i,1990),i=1,nperyear)
                print *,'       data(1991)       = ',
     +                (data(i,1991),i=1,nperyear)
            endif
        endif
*
        end
