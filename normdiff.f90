program normdiff
!
!       compute the normalized diference of two time series, 
!       like the NAO or SOI
!
!       X'(m,y) = [X(m,y)-Xave(m)]/\sigma(X)   with s.d. yearly or monthly
!       C(m,y) = A'(m,y) - B'(m,y)
!
    implicit none
    include 'param.inc'
    integer :: i,j,nperyear,my1,my2,n(npermax),iadd,nmetadata1
    real :: data1(npermax,yrbeg:yrend),data2(npermax,yrbeg:yrend)
    character :: file1*1023,file2*1023,var1*40,var2*40,units1*20,units2*20,string*3
    character :: lvar1*120,lvar2*120,svar1*120,svar2*120,history1*50000,history2*50000
    character :: metadata1(2,100)*2000,metadata2(2,100)*2000,title*250
    logical :: lwrite
    integer :: iargc,llen

    lwrite = .false.
    if ( iargc() < 3 ) then
        print *,'usage: normdiff series1 series2 '// &
 &           'none|monthly|yearly|full none|monthly|yearly|full'// &
 &           ' [add|ave]'
        stop
    endif
    
    call getarg(1,file1)
    call readseriesmeta(file1,data1,npermax,yrbeg,yrend,nperyear,var1 &
 &       ,units1,lvar1,svar1,history1,metadata1,.false.,lwrite)
    call getarg(2,file2)
    if ( file2 == 'null' .or. file2 == 'nothing' ) then
        data2 = 0
        i = nperyear
        var2 = 'null'
        units2 = ' '
        lvar2 = ' '
        svar2 = ' '
        history2 = ' '
        metadata2 = ' '
    else
        call readseriesmeta(file2,data2,npermax,yrbeg,yrend,i,var2 &
 &           ,units2,lvar2,svar2,history2,metadata2,.false.,lwrite)
    endif
    write(title,'(4a)') 'Difference between ',trim(file1),' and ',trim(file2)
    print '(2a)','# ',trim(title)
    if ( i /= nperyear ) then
        write(0,*) 'normdiff: cannot interpolate in time (yet)',i,nperyear
        write(*,*) 'normdiff: cannot interpolate in time (yet)',i,nperyear
        call exit(-1)
    endif
    call getmy(3,my1,nperyear)
    call getmy(4,my2,nperyear)
    call getarg(5,string)
    if ( string == 'add' ) then
        iadd = 1
    else if ( string == 'ave' ) then
        iadd = 2
    else
        iadd = 0
    end if
    if ( my1 >= 0 ) then
        call anomal(data1,npermax,nperyear,yrbeg,yrend,yrbeg,yrend)
        call normalize(my1,data1,nperyear)
    end if
    if ( file2 /= 'null' .and. file2 /= 'nothing' ) then
        if ( my1 >= 0 ) then
            call anomal(data2,npermax,nperyear,yrbeg,yrend,yrbeg &
 &               ,yrend)
            call normalize(my1,data2,nperyear)
        end if
    endif
    do i=yrbeg,yrend
        do j=1,nperyear
            if ( data1(j,i) < 1e33 .and. data2(j,i) < 1e33 ) then
                if ( iadd == 0 ) then
                    data1(j,i) = data1(j,i) - data2(j,i)
                else if ( iadd == 1 ) then
                    data1(j,i) = data1(j,i) + data2(j,i)
                else if ( iadd == 2 ) then
                    data1(j,i) = (data1(j,i) + data2(j,i))/2
                end if
            else
                data1(j,i) = 3e33
            endif
        enddo
    enddo
    if ( my2 >= 0 ) then
        call normalize(my2,data1,nperyear)
    end if
    if ( my1 <= 0 .and. my2 <= 0 ) then
        if ( iadd == 0 ) then
            write(*,'(8a)') '# diff [',trim(units1),'] difference ', &
 &               'of ',trim(var1),' and ',trim(var2)
        else if ( iadd == 1 ) then
            write(*,'(8a)') '# sum [',trim(units1),'] sum ', &
 &               'of ',trim(var1),' and ',trim(var2)
        else if ( iadd == 2 ) then
            write(*,'(7a)') '# ave [',trim(units1),'] average ', &
 &               'of ',trim(var1),' and ',trim(var2)
        end if
    else
        if ( iadd == 0 ) then
            write(*,'(7a)') '# diff [1] normalised difference of ' &
 &               ,trim(var1),' and ',trim(var2)
        else if ( iadd == 1 ) then
            write(*,'(7a)') '# sum [1] normalised sum of ' &
 &               ,trim(var1),' and ',trim(var2)
        else
            write(*,'(7a)') '# ave [1] normalised average of ' &
 &               ,trim(var1),' and ',trim(var2)
        end if
    end if
    call merge_metadata(metadata1,nmetadata1,metadata2,' ',history2,'series2_')
    if ( file2 /= 'null' .and. file2 /= 'nothing' ) then
        nmetadata1 = nmetadata1 + 1
        metadata1(1,nmetadata1) = 'file2'
        metadata1(2,nmetadata1) = file2    
    end if
    call printmetadata(6,file1,' ',title,history1,metadata1)
    call printdatfile(6,data1,npermax,nperyear,yrbeg,yrend)
end program

subroutine getmy(i,my,nperyear)
    implicit none
    integer i,my,nperyear
    character*1 chr
    call getarg(i,chr)
    if ( chr == 'n' ) then
        my = 0
    elseif ( chr == 'm' ) then
        my = nperyear
        print '(a)','# Timeseries are normalized per month'
    elseif ( chr == 'y' ) then
        my = 1
        print '(a)','# Timeseries are normalized per year'
    elseif ( chr == 'f' ) then
        my = -1
        print '(a)','# Full timeseries are not normalized'
    else
        write(0,*) 'normdiff: expecting ''m'' or ''y'', not ',chr
        write(*,*) 'normdiff: expecting ''m'' or ''y'', not ',chr
        call exit(-1)
    endif
end subroutine

subroutine normalize(my,data,nperyear)
    implicit none
    include 'param.inc'
    integer my,nperyear
    real data(npermax,yrbeg:yrend)
    integer i,j,jj,n,nn(npermax)
    real s1(npermax),s2(npermax)

    n = min(nperyear,my)
    if ( n == 0 ) return
    do jj=1,n
        nn(jj) = 0
        s1(jj) = 0
        s2(jj) = 0
    enddo
    do i=yrbeg,yrend
        do j=1,nperyear
            if ( data(j,i) < 1e33 ) then
                jj = min(j,my)
                nn(jj) = nn(jj) + 1
                s1(jj) = s1(jj) + data(j,i)
                s2(jj) = s2(jj) + data(j,i)**2
            endif
        enddo
    enddo
    do jj=1,n
        if ( nn(jj) > 1 ) then
            s1(jj) = s1(jj)/nn(jj)
            s2(jj) = sqrt(s2(jj)/nn(jj) - s1(jj)**2)
        else
            s1(jj) = 3e33
            s2(jj) = 3e33
        endif
    enddo
    do i=yrbeg,yrend
        do j=1,nperyear
            if ( data(j,i) < 1e33 ) then
                jj = min(j,my)
                if ( s1(jj) < 1e33 ) then
                    data(j,i) = (data(j,i)-s1(jj))/s2(jj)
                else
                    data(j,i) = 3e33
                endif
            endif
        enddo
    enddo
end subroutine
