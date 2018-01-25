program selectyear

!   trivial little program to select a given range of years from a
!   standard station data file

    implicit none
    include 'param.inc'
    integer :: yr,mo,i,nperyear,yr1,yr2
    real :: data(npermax,yrbeg:yrend)
    character :: file*1023,line*128,var*40,units*80,lvar*120,svar*120,history*50000, &
        metadata(2,100)*2000,title*500
    logical :: lvalid,lwrite
    integer :: iargc,llen

    lwrite = .false. 
    if ( iargc() < 3 .or. iargc() > 4 ) then
        print *,'usage: selectyear yrbegin yrend datfile [dummy]'
        call exit(-1)
    endif
    call getarg(1,line)
    read(line,*,err=901) yr1
    call getarg(2,line)
    read(line,*,err=901) yr2
    call getarg(3,file)
    call readseriesmeta(file,data,npermax,yrbeg,yrend,nperyear,var,units &
        ,lvar,svar,history,metadata,.false.,lwrite)
    title = ' '
    call copyheadermeta(file,6,title,history,metadata)
    write(6,'(8a)') '# ',trim(var),' [',trim(units),'] ',trim(lvar)
    call printdatfile(6,data(1,yr1),npermax,nperyear,yr1,yr2)
    goto 999
901 write(0,*) 'selectyear: error reading year from ',line
    call exit(-1)
999 continue
end program selectyear
