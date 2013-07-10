        program dat2grads
*
*       convert a dat file from correlate into a grads ctl/dat file combo
*       with 1 data point.  Useful to check field versions of programs against 
*       point versions :-(
*
        implicit none
#include "recfac.h"
#include "param.inc"
        integer i,nperyear,yr,mo,yr1,yr2
        real data(npermax,yrbeg:yrend),lon,lat,alt
        character infile*256,file*256,line*128,unit*3,var*40,units*20
        logical lwrite
        integer iargc,llen,get_endian

        lwrite = .false.
        if ( iargc().ne.2 ) then
            print *,'usage: dat2grads infile.dat outfile.ctl'
            stop
        endif
        call getarg(1,infile)
        lon = 0
        lat = 0
        alt = 0
        open(1,file=infile,status='old')
        do
            read(1,'(a)') line
            if ( line(1:1).ne.'#' ) exit
            i = index(line,'N, ')
            if ( i.ne.0 ) call getcoordinate(line,i,lat)
            i = index(line,'E, ')
            if ( i.ne.0 ) call getcoordinate(line,i,lon)
            i = index(line,'m       ')
            if ( i.ne.0 ) call getcoordinate(line,i,alt)
        end do
        close(1)
        
        call readseries(infile,data,npermax,yrbeg,yrend,nperyear,var
     +       ,units,.false.,lwrite)
        yr1 = yrend
        yr2 = yrbeg
        do yr=yrbeg,yrend
            do mo=1,nperyear
                if ( data(mo,yr).lt.1e33 ) then
                    yr1 = min(yr1,yr)
                    yr2 = max(yr2,yr)
                endif
            enddo
        enddo
        call getarg(2,file)
        open(2,file=file,status='unknown')
        
        i = index(file,'.ctl')
        if ( i.eq.0 ) then
            write(0,*) 'error: cannot find ''.ctl'' in ',file
            call abort
        endif
        file(i:) = '.grd'
        open(3,file=file,status='unknown',form='unformatted',
     +       access='direct',recl=recfa4)
        i = 0
        do yr=yr1,yr2
            do mo=1,nperyear
                i = i + 1
                write(3,rec=i) data(mo,yr)
            enddo
        enddo
        close(3)

        write(2,'(2a)') 'DSET ^',file(1:llen(file))
        write(2,'(2a)') 'TITLE converted with dat2grads from ',
     +       infile(1:llen(infile))
        write(2,'(a)') 'UNDEF 3e33'
        if ( get_endian().eq.-1 ) then
            write(2,'(a)') 'OPTIONS LITTLE_ENDIAN'
        else
            write(2,'(a)') 'OPTIONS BIG_ENDIAN'
        endif
        write(2,'(a,f8.3,a)') 'XDEF 1 LINEAR ',lon,' 1'
        write(2,'(a,f8.3,a)') 'YDEF 1 LINEAR ',lat,' 1'
        write(2,'(a,f8.1,a)') 'ZDEF 1 LINEAR ',alt,' 1'
        if ( nperyear.eq.1 ) then
            unit = '1YR'
        elseif ( nperyear.eq.12 ) then
            unit = '1MO'
        elseif ( nperyear.eq.73 ) then
            unit = '5DY'
        elseif ( nperyear.eq.366 ) then
            unit = '1DY'
        else
            write(0,*) 'error: cannot handle nperyear = ',nperyear,
     +           ' yet'
            call abort
        endif
        write(2,'(a,i6,a,i4.4,2a)') 'TDEF ',nperyear*(yr2-yr1+1),
     +       ' LINEAR 1JAN',yr1,' ',unit
        write(2,'(a)') 'VARS 1'
        write(2,'(6a)') trim(var),' 0 99 ',trim(var),' [',trim(units)
     +       ,']'
        write(2,'(a)') 'ENDVARS'
        close(2)
        end

        subroutine getcoordinate(line,i,lat)
        implicit none
        integer i
        real lat
        character line*(*)
        integer j
        if ( i.ne.0 ) then
            j = i-1
            do while ( line(j:j).eq.'.' .or. line(j:j).eq.'-' .or.
     +           (ichar(line(j:j)).ge.ichar('0') .and.
     +           ichar(line(j:j)).le.ichar('9') ) )
                j = j - 1
            end do
            j = j + 1
            if ( j.eq.i ) return
            print *,'getcoordinate: coordinate = ',line(j:i-1)
            read(line(j:i-1),*) lat
            print *,'getcoordinate: coordinate = ',lat
        endif
        end
