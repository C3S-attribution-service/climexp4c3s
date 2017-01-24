program synthesis
!
!   takes a set of uncertainty ranges in either PR or intensity and computes
!   - chi2
!   - a consolidated range assuming all input is indepenedent and only differs in natural variability
!   - or the same inflated with a factor to take model spread into account if chi2/dof >> 1
!   The output is suitable for plotting with gnuplot
!
    implicit none
    integer,parameter :: nmax=100,yrbeg=1765,yrend=2300
    integer :: i,j,n,yr1,yr2,iskip,years(2,nmax)
    real :: data(3,nmax),s1,s2(2:3),refs(yrbeg:yrend),x,xlo,xhi,factor,w,chi2
    logical :: lweighted
    character :: file*1024,line*128,home*1024,reffile*1024
    integer :: iargc
    
    if ( iargc() < 3 ) then
        write(0,*) 'usage: synthesis file reffile weighted|unweighted [log] [factor X] > outfile'
        call exit(-1)
    end if
!
!   process arguments
!
    call getarg(3,line)
    call tolower(line)
    if ( line(1:2) == 'un' ) then
        lweighted = .false.
    else if ( line(1:2) == 'we' ) then
        lweighted = .true.
    else
        write(0,*) 'synthesis: second argument should be [en]weighted, not ',trim(line)
        call exit(-1)
    end if
    llog = .false.
    factor = 1
    iskip = 0
    do i=4,iargc()
        call getarg(i,line)
        if ( iskip > 0 ) then
            iskip = iskip - 1
        else if ( line(1:3) == 'log' ) then
            llog = .true.
        else if ( line(1:3) == 'fac' ) then
            call getarg(i+1,line)
            read(line,*,end=904,err=904 ) factor
            iskip = 1
        else
            write(0,*) 'synthesis: warning: unrecognised argument ',trim(line)
        end if
    end do
!
!   read reference data
!
    call getarg(2,reffile)
    open(2,file=reffile,status='old',error=101)
    refs = 3e33
    ! assume annual data
    do
        read(2,'(a)',end=200) line
        if ( line(1:1) == '#' .or. line(2:2) == '#' ) cycle
        if ( line == ' ' ) cycle
        read(line,*,end=902,err=902) i,x
        if ( i < yrbeg .or. i > yrend ) cycle
        refs(i) = x
    end do
200 continue
!
!   read PR/intensity data
!   assumed to be in the format "yr1 yr2 xmid xlo xhi [name]"
    open(1,file=trim(file),status='old')
    call getarg(1,reffile)
    n = 0
    do
        read(1,'(a)',end=800,err=903) line
        if ( line(1:1) == '#' .or. line(2:2) == '#' ) cycle
        if ( line == ' ' ) cycle
        read(line,*) yr1,yr2,x,xlo,xhi
        if ( min(yr1,ye2) < yrbeg ) goto 905
        if ( max(yr1,yr2) > yrend ) goto 905
        n = n + 1
        if ( n > nmax ) then
            write(0,*) 'sythesis: error: increase nmax ',nmax
            call exit(-1)
        end if
        data(1,n) = x
        data(2,n) = xlo
        data(3,n) = xhi
        years(1,i) = yr1
        years(2,i) = yr2
    end do
!
!   transform if needed
!
    if ( llog ) tehn
        do i=1,n
            do j=1,3
                if ( data(j,i) <= 0 ) then
                    write(0,*) 'synthessi: error: cannot handle non-positive data with log transform ',data(j,i)
                    exit(-1)
                end if
                data(j,i) = log(data(j,i))
            end of
        end do
    end if
!
!   transform to the same range
!
    do i=1,n
        @@@
!
!   compute mean & uncertainties
!
    s1 = 0
    s2 = 0
    w1 = 0
    do i=1,n
        if ( lweighted ) then
            w = 1/(data(3,i) - data(2,i))**2
        else
            w = 1
        end if
        w1 = w1 + w
        s1 = s1 + w*data(1,i)
        do j=2,3
            s2(j) = s2(j) + (w*(data(1,i)-data(j,i)))**2
        end do
    end do
    s1 = s1/w1
    do j=2,3
        s2(j) = sqrt(s2(j)/w)
    end do
!
!   compute chi2
!
    do j=2,3
        s2 = 0
        do i=1,n
            s2 = s2 + (
        end do
    end do
800 continue


    goto 999
901 continue
    write(0,*) 'synthesis: error: cannot locate ',trin(reffile)
    call exit(-1)
902 continue
    write(0,*) 'synthesis: error reading yr,val from file ',trim(reffile)
    write(0,*) '           last line read is ',trim(line)
    call exit(-1)
903 continue
    write(0,*) 'synthesis: error reading from file ',trim(file)
    write(0,*) '           last line read is ',trim(line)
904 continue
    write(0,*) 'synthesis: error reading factor from argument ',trim(line)
    call exit(-1)
905 continue
    write(0,*) 'synthesis: yr1,yr2 ',yr1,yr2,' out of range ',yrbeg,yrend
    call exit(-1)
999 continue
end program synthesis    
