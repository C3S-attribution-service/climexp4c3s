        program normdiff
*
*       compute the normalized diference of two time series, 
*       like the NAO or SOI
*
*       X'(m,y) = [X(m,y)-Xave(m)]/\sigma(X)   with s.d. yearly or monthly
*       C(m,y) = A'(m,y) - B'(m,y)
*
        implicit none
#include "param.inc"
        integer i,j,nperyear,my1,my2,n(npermax),iadd
        real  data1(npermax,yrbeg:yrend),data2(npermax,yrbeg:yrend)
        character*255 file,var1*40,var2*40,units1*20,units2*20,string*3
        logical lwrite
        integer iargc,llen

        lwrite = .false.
        if ( iargc().lt.3 ) then
            print *,'usage: normdiff series1 series2 '//
     +           'none|monthly|yearly|full none|monthly|yearly|full'//
     +           ' [add|ave]'
            stop
        endif
        
        print '(a)','# Difference between'
        call getarg(1,file)
        print '(2a)','# ',file(1:llen(file))
        call readseries(file,data1,npermax,yrbeg,yrend,nperyear,var1
     +       ,units1,.false.,lwrite)
        call getarg(2,file)
        print '(2a)','# ',file(1:llen(file))
        if ( file.eq.'null' .or. file.eq.'nothing' ) then
            data2 = 0
            i = nperyear
            units2 = ' '
            var2 = 'null'
        else
            call readseries(file,data2,npermax,yrbeg,yrend,i,var2
     +           ,units2,.false.,lwrite)
        endif
        if ( i.ne.nperyear ) then
            write(0,*) 'normdiff: cannot interpolate in time (yet)'
            write(*,*) 'normdiff: cannot interpolate in time (yet)'
            call abort
        endif
        call getmy(3,my1,nperyear)
        call getmy(4,my2,nperyear)
        call getarg(5,string)
        if ( string.eq.'add' ) then
            iadd = 1
        else if ( string.eq.'ave' ) then
            iadd = 2
        else
            iadd = 0
        end if
        if ( my1.ge.0 ) then
            call anomal(data1,npermax,nperyear,yrbeg,yrend,yrbeg,yrend)
            call normalize(my1,data1,nperyear)
        end if
        if ( file.ne.'null' .and. file.ne.'nothing' ) then
            if ( my1.ge.0 ) then
                call anomal(data2,npermax,nperyear,yrbeg,yrend,yrbeg
     +               ,yrend)
                call normalize(my1,data2,nperyear)
            end if
        endif
        do i=yrbeg,yrend
            do j=1,nperyear
                if ( data1(j,i).lt.1e33 .and. data2(j,i).lt.1e33 ) then
                    if ( iadd.eq.0 ) then
                        data1(j,i) = data1(j,i) - data2(j,i)
                    else if ( iadd.eq.1 ) then
                        data1(j,i) = data1(j,i) + data2(j,i)
                    else if ( iadd.eq.2 ) then
                        data1(j,i) = (data1(j,i) + data2(j,i))/2
                    end if
                else
                    data1(j,i) = 3e33
                endif
            enddo
        enddo
        if ( my2.ge.0 ) then
            call normalize(my2,data1,nperyear)
        end if
        if ( my1.le.0 .and. my2.le.0 ) then
            if ( iadd.eq.0 ) then
                write(*,'(8a)') '# diff [',trim(units1),'] difference ',
     +               'of ',trim(var1),' and ',trim(var2)
            else if ( iadd.eq.1 ) then
                write(*,'(8a)') '# sum [',trim(units1),'] sum ',
     +               'of ',trim(var1),' and ',trim(var2)
            else if ( iadd.eq.2 ) then
                write(*,'(7a)') '# ave [',trim(units1),'] average ',
     +               'of ',trim(var1),' and ',trim(var2)
            end if
        else
            if ( iadd.eq.0 ) then
                write(*,'(7a)') '# diff [1] normalised difference of '
     +               ,trim(var1),' and ',trim(var2)
            else if ( iadd.eq.1 ) then
                write(*,'(7a)') '# sum [1] normalised sum of '
     +               ,trim(var1),' and ',trim(var2)
            else
                write(*,'(7a)') '# ave [1] normalised average of '
     +               ,trim(var1),' and ',trim(var2)
            end if
        end if
        call printdatfile(6,data1,npermax,nperyear,yrbeg,yrend)
        end

        subroutine getmy(i,my,nperyear)
        implicit none
        integer i,my,nperyear
        character*1 chr
        call getarg(i,chr)
        if ( chr.eq.'n' ) then
            my = 0
            print '(a)','# Timeseries are not normalized'
        elseif ( chr.eq.'m' ) then
            my = nperyear
            print '(a)','# Timeseries are normalized per month'
        elseif ( chr.eq.'y' ) then
            my = 1
            print '(a)','# Timeseries are normalized per year'
        elseif ( chr.eq.'f' ) then
            my = -1
            print '(a)','# Full timeseries are not normalized'
        else
            write(0,*) 'normdiff: expecting ''m'' or ''y'', not ',chr
            write(*,*) 'normdiff: expecting ''m'' or ''y'', not ',chr
            call abort
        endif
        end

        subroutine normalize(my,data,nperyear)
        implicit none
#include "param.inc"
        integer my,nperyear
        real data(npermax,yrbeg:yrend)
        integer i,j,jj,n,nn(npermax)
        real s1(npermax),s2(npermax)

        n = min(nperyear,my)
        if ( n.eq.0 ) return
        do jj=1,n
            nn(jj) = 0
            s1(jj) = 0
            s2(jj) = 0
        enddo
        do i=yrbeg,yrend
            do j=1,nperyear
                if ( data(j,i).lt.1e33 ) then
                    jj = min(j,my)
                    nn(jj) = nn(jj) + 1
                    s1(jj) = s1(jj) + data(j,i)
                    s2(jj) = s2(jj) + data(j,i)**2
                endif
            enddo
        enddo
        do jj=1,n
            if ( nn(jj).gt.1 ) then
                s1(jj) = s1(jj)/nn(jj)
                s2(jj) = sqrt(s2(jj)/nn(jj) - s1(jj)**2)
            else
                s1(jj) = 3e33
                s2(jj) = 3e33
            endif
        enddo
        do i=yrbeg,yrend
            do j=1,nperyear
                if ( data(j,i).lt.1e33 ) then
                    jj = min(j,my)
                    if ( s1(jj).lt.1e33 ) then
                        data(j,i) = (data(j,i)-s1(jj))/s2(jj)
                    else
                        data(j,i) = 3e33
                    endif
                endif
            enddo
        enddo
        end
