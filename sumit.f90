subroutine sumit(data,npermax,nperyear,yrbeg,yrend,lsum,oper)

!   sum lsum consecutive months/days/... into the first one

    implicit none
    integer :: npermax,nperyear,yrbeg,yrend,lsum
    real :: data(npermax,yrbeg:yrend)
    character(1) :: oper
    integer :: i,j,m,n,ii
    real :: s
    real,parameter :: absent=3e33
    logical,parameter :: lwrite=.FALSE.
    integer,external :: leap

    if ( lsum > 1 ) then
        if ( lwrite ) then
            print *,'sumit: npermax,nperyear = ',npermax,nperyear
            print *,'       yrbeg,yrend      = ',yrbeg,yrend
            print *,'       lsum,oper        = ',lsum,oper
            print *,'       data(1990)       = ',(data(i,1990),i=1,nperyear)
            print *,'       data(1991)       = ',(data(i,1991),i=1,nperyear)
        endif
        do i=yrbeg,yrend
            do j=1,nperyear
                s = data(j,i)
                if ( data(j,i) < 0.9*absent ) then
                    do n=1,lsum-1
                        m = j+n
                        if ( nperyear == 366 .AND. leap(i) == 1 .AND. j < 60 .AND. m >= 60 ) m = m + 1
                        call normon(m,i,ii,nperyear)
                        if ( ii > yrend ) then
                            s = absent
                        elseif ( s < 0.9*absent .AND. data(m,ii) < 0.9*absent ) then
                            if ( oper == '+' .OR. oper == 'v' ) then
                                s = s + data(m,ii)
                            elseif ( oper == 'a' ) then
                                s = max(s,data(m,ii))
                            elseif ( oper == 'i' ) then
                                s = min(s,data(m,ii))
                            else
                                print *,'sumit: error: cannot handle oper = ',oper,', only [+vai]'
                                call abort
                            endif
                        else
                            s = absent
                        endif
                    enddo
                    if ( oper == 'v' .AND. s < 0.9*absent ) then
                        s = s/lsum
                    endif
                    data(j,i) = s
                endif
            enddo
        enddo
        if ( lwrite ) then
            print *,'sumit: after averaging'
            print *,'       data(1990)       = ',(data(i,1990),i=1,nperyear)
            print *,'       data(1991)       = ',(data(i,1991),i=1,nperyear)
        endif
    endif

end subroutine sumit
