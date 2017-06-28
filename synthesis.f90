program synthesis
!
!   takes a set of uncertainty ranges in either PR or intensity and computes
!   - chi2
!   - a consolidated range assuming all input is independent and only differs in natural variability
!   - or the same inflated with a factor to take model spread into account if chi2/dof >> 1
!   - or an unweighted average
!   - or no average
!   Input is a file in the format 
!   # title
!   yr0 yr1 bestfit lowbound highbound conf name
!   ...
!   A blank line separates observations from models
!   Conf is 95 for 95% CI and so on.
!   Arguments:
!   - Infile: see above
!   - Reffile to transform begin dates to the earliest one, can be abbreviated to GMST, CO2, time or none
!   - [un]weighted: see above.
!   - [log] Apply logarithm before assuming normal distributions (differing for lower and higher half always)
!   - [perc] The input is a percentage, will be converted to a fraction 1+perc/100
!   - [flipsign] Indicates it is a reverse attribution, ie, the numbers should be normalised to the highest year.
!   - [factor X] auto|number inflate error bars to make chi1/dof one, or by hand
!   The output is suitable for plotting with gnuplot
!
    implicit none
    integer,parameter :: nmax=100,yrbeg=1765,yrend=2300
    integer :: i,j,n,yr1,yr2,iskip,years(2,nmax),nblank(nmax),oldyears(2,nmax),m,ncolor,colors(3)
    real :: data(3,nmax),s1,s2,w1,ss2(2:3),refs(yrbeg:yrend),x,xlo,xhi,factor,w,chi2,f,ci(nmax),cc
    logical :: lnoave,lweighted,llog,lperc,lflip,printit
    character :: file*1024,line*128,home*1024,reffile*1024,names(nmax)*40
    integer :: iargc
    
    if ( iargc() < 3 ) then
        write(0,*) 'usage: synthesis file reffile weighted|unweighted|noave [log] [perc] [flipsign] [factor X] > outfile'
        call exit(-1)
    end if
!
!   process arguments
!
    call getarg(3,line)
    call tolower(line)
    lweighted = .false.
    lnoave = .false.
    if ( line(1:2) == 'un' ) then
        lweighted = .false.
    else if ( line(1:2) == 'we' ) then
        lweighted = .true.
    else if ( line(1:2) == 'no' ) then
        lnoave = .true.
    else
        write(0,*) 'synthesis: third argument should be [un]weighted, not ',trim(line)
        call exit(-1)
    end if
    llog = .false.
    lperc = .false.
    lflip = .false.
    factor = 1
    iskip = 0
    do i=4,iargc()
        call getarg(i,line)
        if ( iskip > 0 ) then
            iskip = iskip - 1
        else if ( line(1:3) == 'log' ) then
            llog = .true.
        else if ( line(1:3) == 'per' ) then
            lperc = .true.
        else if ( line(1:3) == 'fli' ) then
            lflip = .true.
        else if ( line(1:3) == 'fac' ) then
            call getarg(i+1,line)
            if ( line(1:4) == 'auto' ) then
                factor = -1
            else
                read(line,*,end=904,err=904 ) factor
            end if
            iskip = 1
        else
            write(0,*) 'synthesis: warning: unrecognised argument ',trim(line)
        end if
    end do
!
!   read reference data
!
    call getarg(2,reffile)
    if ( reffile == 'GMST' ) reffile = 'NASAData/giss_al_gl_a_4yrlo.dat'
    if ( reffile == 'CO2' ) reffile = 'CDIACData/co2_annual.dat'
    if ( reffile == 'time' ) then
        do i=yrbeg,yrend
            refs(i) = i
        end do
    else if ( reffile == 'none' ) then
        refs = 0
    else
        open(2,file=reffile,status='old',err=901)
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
200     continue
    end if
!
!   read PR/intensity data
!   assumed to be in the format "yr1 yr2 xmid xlo xhi ci name"
!
    call getarg(1,file)
    open(1,file=trim(file),status='old')
    n = 0
    nblank = 0
    do
        read(1,'(a)',end=800,err=903) line
        if ( line(1:1) == '#' .or. line(2:2) == '#' ) then
            print '(a)',trim(line)
            cycle
        end if
        if ( line == ' ' ) then
            if ( n > 0 ) nblank(n) = nblank(n) + 1
            cycle
        end if
        read(line,*) yr1,yr2,x,xlo,xhi,cc
        if ( min(yr1,yr2) < yrbeg ) goto 905
        if ( max(yr1,yr2) > yrend ) goto 905
        n = n + 1
        if ( n > nmax ) then
            write(0,*) 'sythesis: error: increase nmax ',nmax
            call exit(-1)
        end if
        years(1,n) = yr1
        years(2,n) = yr2
        data(1,n) = x
        data(2,n) = xlo
        data(3,n) = xhi
        ci(n) = cc
        do j=1,len(line)
            if ( ( ichar(line(j:j)).ge.ichar('a') .and. ichar(line(j:j)).le.ichar('z') &
              .or. ichar(line(j:j)).ge.ichar('A') .and. ichar(line(j:j)).le.ichar('Z') ) & 
                .and. line(j:j+1) /= 'E-' .and. line(1:2) /= 'E+' ) then
                names(n) = line(j:)
                exit
            end if
        end do
    end do
800 continue
!
!   sanity check
!
    if ( .not. lnoave .and. lweighted .and. factor == -1 .and. n <= 3 ) then
        write(0,*) 'Autoscaling does not work for such a small number of results.'
        write(0,*) 'Specifiy the scaling by hand or use an unweighted average'
        write(0,*) 'I''ll do the latter now.'
        lweighted = .false.
    end if
!
!   transform to logarithm if needed
!
    if ( llog ) then
        do i=1,n
            do j=1,3
                if ( data(j,i) <= 0 ) then
                    write(0,*) 'synthesis: error: cannot handle non-positive data with log transform ',data(j,i)
                    call exit(-1)
                end if
                data(j,i) = log(data(j,i))
            end do
        end do
    end if
    if ( lperc ) then
        do i=1,n
            do j=1,3
                if ( data(j,i) <= -100 ) then
                    write(0,*) 'synthesis: error: cannot handle non-positive data with log transform ',data(j,i)
                    call exit(-1)
                end if
                data(j,i) = log(1+data(j,i)/100)
            end do
        end do
    end if
!
!   transform to the same (maximum) range
!
    if ( reffile /= 'none' ) then
        yr1 = years(1,1)
        s1 = refs(years(1,1))
        do i=2,n
            if ( s1 > 1e33 .or. refs(years(1,i)) > 1e33 ) then
                if ( s1 < 1e33 ) yr1 = years(1,i)
                write(0,*) 'synthesis: error: reference series is undefined at ',yr1
                call exit(-1)
            end if
            if ( refs(years(1,i)) < s1 .neqv. lflip ) then
                yr1 = years(1,i)
                s1 = refs(years(1,i))
            end if
        end do
        yr2 = years(2,1)
        s2 = refs(years(2,1))
        do i=2,n
            if ( yr2 /= years(2,i) .neqv. lflip ) then
                yr2 = years(2,i)
                s2 = refs(years(2,i))
            end if
        end do
        m = 1
        oldyears(1,1) = yr1
        oldyears(2,1) = yr2
        do i=1,n
            f = (s2-s1)/(refs(years(2,i))-refs(years(1,i)))
            do j=1,3
                data(j,i) = f*data(j,i)
            end do
            printit = .true.
            do j=1,m
                if ( years(1,i) == oldyears(1,j) .and. &
                     years(2,i) == oldyears(2,j) ) printit = .false.
            end do
            if ( printit ) then
                print '(a,i4,a,i4,a,i4,a,i4,a,f6.2)', '# transformed from ', &
                    years(1,i),'-',years(2,i),' to ',yr1,'-',yr2,' with factor ',f
                m = m + 1
                oldyears(1,m) = years(1,i)
                oldyears(2,m) = years(2,i)
            end if
        end do
        print '(a,i4,a,i4)','# using common interval ',yr1,'-',yr2
        print '(a,2g14.4)','# with reference values ',s1,s2 
    end if
!
!   compute mean & uncertainties
!
    s1 = 0
    ss2 = 0
    w1 = 0
    do i=1,n
        if ( lweighted ) then
            if ( abs(data(3,i) - data(2,i)) < 1e-6 ) then
                write(0,*) 'synthesis: error: upper and lower bound are equal: ',i,data(2,i),data(3,i)
                call exit(-1)
            end if
            w = 1/(data(3,i) - data(2,i))**2
        else
            w = 1
        end if
        w1 = w1 + w
        s1 = s1 + w*data(1,i)
        do j=2,3
            ss2(j) = ss2(j) + (w*(data(1,i)-data(j,i)))**2
        end do
    end do
    s1 = s1/w1
    do j=2,3
        ss2(j) = sqrt(ss2(j)/w1)
    end do
    nblank(n) = 1
    n = n + 1
    data(1,n) = s1
    data(2,n) = s1 - ss2(2)
    data(3,n) = s1 + ss2(3)
    years(1,n) = yr1
    years(2,n) = yr2
    if ( lweighted ) then
        names(n) = 'Weighted average'
    else
        names(n) = 'Average'
    end if
!
!   compute chi2
!
    chi2 = 0
    do i=1,n-1
        if ( data(1,i) > s1 ) then
            chi2 = chi2 + ((data(1,i)-s1)/(data(1,i)-data(2,i)))**2
        else
            chi2 = chi2 + ((s1-data(1,i))/(data(3,i)-data(1,i)))**2
        end if
        if ( ci(i) /= 95 ) then
            write(0,*) 'Confidence intervals other than 95% not yet supported',i,ci(i),' ',names(i)
            call exit(-1)
        end if
    end do
    chi2 = chi2*4 ! transformed to use 95%~2sigma CIs
    print '(a,g14.4)','# chi2/dof = ',chi2/(n-2)
    if ( .not.lnoave .and. chi2/(n-2) > 1 .and. factor == -1 ) then
        factor = sqrt(chi2/(n-2)) ! n-1 but n has been increased by 1 to accomodate the mean
    end if
!
!   scale error bounds
!
    if ( .not.lnoave .and. lweighted .and. factor /= 1 ) then
        data(2,n) = s1 - factor*ss2(2)
        data(3,n) = s1 + factor*ss2(3)
        write(names(n)(len_trim(names(n))+1:),'(a,f4.2)') ' inflated ',factor
    end if
!
!   transform back
!
    if ( llog ) then
        do i=1,n
            do j=1,3
                data(j,i) = exp(data(j,i))
            end do
        end do
    end if
    if ( lperc ) then
        do i=1,n
            do j=1,3
                data(j,i) = 100*(-1+exp(data(j,i)))
            end do
        end do
    end if
!
!   output
!
    colors(1) = 3 ! blue in classic gnuplot colours
    colors(2) = 1 ! red
    colors(3) = 4 ! purple
    ncolor = 1
    if ( lnoave ) n = n-1
    do i=1,n
        print '(2i5,3f12.3,i3,4a)',yr1,yr2,(data(j,i),j=1,3),colors(ncolor),' "',trim(names(i)),'"'
        do j=1,nblank(i)
            ncolor = ncolor + 1
        end do
    end do
!
!   error messages, ye olde Fortran way
!
    goto 999
901 continue
    write(0,*) 'synthesis: error: cannot locate ',trim(reffile)
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
