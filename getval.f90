program getval

!   retrieve a value from a time series,
!   accompanied by mean and standard deviation

    implicit none
    include 'param.inc'
    include 'getopts.inc'
    integer :: unit,yr,mo,n,nperyear,i,j,m
    real :: data(npermax,yrbeg:yrend),xx(yrend-yrbeg+1),x1,s1
    logical :: lopen
    character file*256,var*20,units*20
    integer :: iargc,llen

    lwrite = .false. 
    if ( iargc() < 3 ) then
        print *,'usage: getval file end2 yr [mon mo [name name]]'
        stop
    endif
    call getarg(1,file)
    call readseries(file,data,npermax,yrbeg,yrend,nperyear,var,units,lstandardunits,lwrite)
    if ( nperyear == 0 ) goto 999
    call getopts(2,iargc(),nperyear,yrbeg,yrend,.false.,0,0)
    if ( yr2a > yrend .or. yr2a < yrbeg ) then
        if ( yr2b > yrend .or. yr2b < yrbeg ) then
            write(*,*) 'getval: error: year undefined',yr2a,yr2b
            write(0,*) 'getval: error: year undefined',yr2a,yr2b
            call exit(-1)
        else
            yr2a = yr2b ! hack
        end if
    endif
    if ( nperyear == 1 ) m1 = 1
    if ( m1 == 0 ) then
        write(*,*) 'getval: error: month undefined'
        write(0,*) 'getval: error: month undefined'
        call exit(-1)
    endif
    call invgetdymo(day0,m1,m,nperyear)
    if ( lwrite ) then
        print *,'getval: dy,mo  ',day0,m1
        print *,'getval: data(',m,yr2a,') ',data(m,yr2a)
    end if
    if ( lsum > 1 ) then
        call sumit(data,npermax,nperyear,yrbeg,yrend,lsum,oper)
    endif
    n = 0
    do yr=yrbeg,yrend
        if ( data(m,yr) < 1e33 ) then
            n = n + 1
            xx(n) = data(m,yr)
        endif
    enddo
    call getmoment(1,xx,n,x1)
    call getmoment(2,xx,n,s1)
    inquire(unit=11,opened=lopen)
    if ( .not. lopen ) then
        unit = 6
    else
        unit = 11
    endif
    i = index(file,'/', .true. ) + 1
    j = index(file,'.') - 1
    if ( data(m,yr2a) < 1e30 ) then
        write(*,'(3a,g16.4,3a)') &
            '<table class=realtable width="100%"><tr><td width=300>' &
            ,trim(namestring),'<td>',data(m,yr2a),' ',trim(units),'</table>'
    endif
    write(unit,'(2i5,g12.4,f6.1,i2,2g12.4,3f6.1)') &
    yr2a,m1,data(m,yr2a),-1.,0,x1,s1,0.,0.,0.
999 continue
end program getval
