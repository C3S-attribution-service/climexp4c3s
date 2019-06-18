program synthesis
!
!   takes a set of uncertainty ranges in either PR or intensity and computes
!   - chi2
!   - a consolidated range assuming all input is independent and only differs in natural variability
!   - or the same with a constant model error into account if chi2/dof >> 1
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
!   - [sigma_mod X] auto|number add model error to make chi1/dof one, or by hand
!   The output is suitable for plotting with gnuplot
!
    use ZbrentToGSL
    implicit none
    integer,parameter :: nmax=100,yrbeg=1765,yrend=2300
    integer :: i,j,n,nobs,yr1,yr2,iskip,years(2,nmax),oldyears(2,nmax),m,ncolor,colors(3), &
        icolor
    real :: data(5,nmax),obs(5),mod(5),syn(5),s,s1,s2,ss2(2:3),refs(yrbeg:yrend),x,xlo,xhi, &
        sig_obs,sig_mod,chi2,f,ci(nmax),cc,w1,w2
    logical :: lnoave,lweighted,llog,lperc,lflip,printit,lwrite
    character :: file*1024,line*128,home*1024,reffile*1024,names(nmax)*40
    real,external :: syn_func

    integer :: syn_n,syn_nobs
    real :: syn_data(5,1000)
    common /syn_dat/ syn_n,syn_nobs,syn_data

    if ( command_argument_count() < 3 ) then
        write(0,*) 'usage: synthesis file reffile weighted|unweighted|noave [log] [perc] [flipsign] [sig_mod X] > outfile'
        call exit(-1)
    end if
!
!   process arguments
!
    call get_command_argument(3,line)
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
    lwrite = .false.
    llog = .false.
    lperc = .false.
    lflip = .false.
    sig_mod = -1
    iskip = 0
    do i=4,command_argument_count()
        call get_command_argument(i,line)
        if ( iskip > 0 ) then
            iskip = iskip - 1
        else if ( line(1:3) == 'deb' .or. line(1:3) == 'lwr' ) then
            lwrite = .true.
        else if ( line(1:3) == 'log' ) then
            llog = .true.
        else if ( line(1:3) == 'per' ) then
            lperc = .true.
        else if ( line(1:3) == 'fli' ) then
            lflip = .true.
        else if ( line(1:3) == 'sig' ) then
            call get_command_argument(i+1,line)
            if ( line(1:4) == 'auto' ) then
                sig_mod = -1
            else
                read(line,*,end=904,err=904 ) sig_mod
            end if
            iskip = 1
        else
            write(0,*) 'synthesis: warning: unrecognised argument ',trim(line)
        end if
    end do
    if ( lperc ) llog = .false.
!
!   read reference data
!
    call get_command_argument(2,reffile)
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
    call get_command_argument(1,file)
    open(1,file=trim(file),status='old')
    n = 0
    do
        read(1,'(a)',end=800,err=903) line
        if ( line(1:1) == '#' .or. line(2:2) == '#' ) then
            print '(a)',trim(line)
            cycle
        end if
        if ( line == ' ' ) then
            if ( n > 0 ) then
                nobs = n
                cycle
            end if
        end if
        read(line,*,err=906) yr1,yr2,x,xlo,xhi,cc
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
        data(4,n) = 3e33
        data(5,n) = 3e33
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
    if ( lwrite ) then
        print *,'read ',nobs,' observations and ',n-nobs,' models'
    end if
!
!   sanity check
!
    if ( .not. lnoave .and. lweighted .and. sig_mod == -1 .and. n <= 3 ) then
        write(0,*) 'Autoscaling does not work for such a small number of results.'
        write(0,*) 'Specifiy the scaling by hand or use an unweighted average'
        write(0,*) 'I''ll do the latter now.'
        lweighted = .false.
    end if
!
!   transform to logarithm if needed
!
    if ( llog ) then
        if ( lwrite ) print *,'taking log'
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
        if ( lwrite ) print *,'converting percentages'
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
    do i=1,n-1
        if ( ci(i) /= 95 ) then
            write(0,*) 'Confidence intervals other than 95% not yet supported',i,ci(i),' ',names(i)
            call exit(-1)
        end if
    end do
!
!   compute mean & uncertainties of observations
!
!   natural variability: assume 100% correlated
    do j=1,3
        s1 = 0
        do i=1,nobs
            s1 = s1 + data(j,i)
        end do
        obs(j) = s1/nobs
    end do

!   representation error from scatter of mean
    if ( nobs > 1 .and. lweighted ) then
        s2 = 0
        do i=1,nobs
            s2 = s2 + (data(1,i)-obs(1))**2
        end do
        sig_obs = sqrt(s2/(nobs-1))
        obs(2) = obs(1) - sqrt( (obs(1)-obs(2))**2 + sig_obs**2 )
        obs(3) = obs(1) + sqrt( (obs(3)-obs(1))**2 + sig_obs**2 )
    else
        sig_obs = 0 ! cannot estimate it from one point...
    end if
    do i=1,nobs
        data(4,i) = data(1,i) - sqrt( (data(1,i)-data(2,i))**2 + sig_obs**2 )
        data(5,i) = data(1,i) + sqrt( (data(3,i)-data(1,i))**2 + sig_obs**2 )
    end do
    obs(4) = obs(2)
    obs(5) = obs(3)
    if ( lwrite ) print *,'found combined observational estimate ',obs
    if ( sig_obs > 0 ) then
        if ( llog ) then
            print '(a,g12.3)','# representation uncertainty (2&sigma;) factor ',exp(sig_obs)
        else if ( lperc ) then
            print '(a,g12.3,a)','# representation uncertainty (2&sigma;) ',100*(-1+exp(sig_obs)),'%'
        else
            print '(a,g12.3)','# representation uncertainty (2&sigma;) ',sig_obs
        end if
    end if
 !
!   compute mean & uncertainties of models
!
    call getsynmean(lweighted,data,n,nobs,sig_mod)
    if ( lwrite ) print *,'first guess of model mean is ',data(:,n+1)
    call getsynchi2(data,n,nobs,sig_mod,chi2)
    if ( lweighted ) then
        print '(a,g14.4)','# model chi2/dof = ',chi2/(n-1)
    end if
    if ( lweighted .and. .not.lnoave .and. chi2/(n-1) > 1 .and. sig_mod < 0 ) then
        ! compute sig_mod to make chi2/dof = 1
        syn_n = n ! copy to common as I am a F77 programmer and do not feel comfortable with f90 global variables
        syn_nobs = nobs
        syn_data(:,1:n) = data(:,1:n)
        s1 = 0
        s2 = data(3,n+1) - data(2,n+1)
        if ( syn_func(s1) <= 0 ) then
            write(0,*) 'synthesis: internal error: syn_func(0) should be > 0, not ',syn_func(s1)
            call exit(-1)
        end if
        s = syn_func(s2)
        do while ( s2 < 1e10 )
            s2 = 2*s2
            if ( s2 > 1e10 ) then
                write(0,*) 'synthesis: error: cannot find zero: ',s2
                call exit(-1)
            end if
            s = syn_func(s2)
            if ( s < 0 ) exit
        end do
        if ( lwrite ) print *,'bracketed sig_mod by ',s1,s2,', calling Brent'
        sig_mod = zbrent(syn_func,s1,s2,1e-4)
        !!!write(names(n)(len_trim(names(n))+1:),'(a,f4.2)') ' added ',sig_mod
        ! copy results back
        mod(1:3) = syn_data(1:3,n+1)
    else
        sig_mod = 0
        mod(1:3) = data(1:3,n+1)
    end if
    mod(4) = mod(2)
    mod(5) = mod(3)
    do i=nobs+1,n
        data(4,i) = data(1,i) - sqrt( (data(1,i)-data(2,i))**2 + sig_mod**2 )
        data(5,i) = data(1,i) + sqrt( (data(3,i)-data(1,i))**2 + sig_mod**2 )
    end do
    if ( lwrite ) print *,'found combined model estimate ',mod
!
!   compute synthesis: weighted average of obervations and models
!
    if ( obs(2) >= obs(3) ) then
        write(0,*) 'synthesis: internal error: obs = ',obs
        call exit(-1)
    end if
    if ( mod(2) >= mod(3) ) then
        write(0,*) 'synthesis: internal error: mod = ',mod
        call exit(-1)
    end if
    if ( lweighted ) then
        ! weighted mean coloured
        w1 = 1/(obs(3)-obs(2))**2
        w2 = 1/(mod(3)-mod(2))**2
        syn(1) = (w1*obs(1) + w2*mod(1))/(w1+w2)
        syn(2) = syn(1) - sqrt( (w1*(obs(1)-obs(2)))**2 + (w2*(mod(1)-mod(2)))**2 )/(w1+w2)
        syn(3) = syn(1) + sqrt( (w1*(obs(3)-obs(1)))**2 + (w2*(mod(3)-mod(1)))**2 )/(w1+w2)
        ! unweighted mean of observations and models box
        s1 = (obs(1) + mod(1))/2
        syn(4) = s1 - sqrt( (obs(1)-obs(2))**2 + (mod(1)-mod(2))**2 )/2
        syn(5) = s1 + sqrt( (obs(3)-obs(1))**2 + (mod(3)-mod(1))**2 )/2
    else
        call getsynmean(lweighted,data,n,0,sig_mod)
        syn(1:3) = data(1:3,n+1)
        syn(4) = syn(2)
        syn(5) = syn(3)
    end if
    if ( lwrite ) print *,'found synthesised estimate ',syn
!
!   transform back
!
    if ( llog ) then
        do i=1,n
            do j=1,5
                data(j,i) = exp(data(j,i))
            end do
        end do
        do j=1,5
            obs(j) = exp(obs(j))
            mod(j) = exp(mod(j))
            syn(j) = exp(syn(j))
        end do
    end if
    if ( lperc ) then
        do i=1,n
            do j=1,5
                data(j,i) = 100*(-1+exp(data(j,i)))
            end do
        end do
        do j=1,5
            obs(j) = 100*(-1+exp(obs(j)))
            mod(j) = 100*(-1+exp(mod(j)))
            syn(j) = 100*(-1+exp(syn(j)))
        end do
    end if
!
!   output
!
   if ( sig_mod > 0 ) then
        if ( llog ) then
            print '(a,g12.3)','# model uncertainty (2&sigma;) factor ',exp(sig_mod)
        else if ( lperc ) then
            print '(a,g12.3,a)','# model uncertainty (2&sigma;) ',100*(-1+exp(sig_mod)),'%'
        else
            print '(a,g12.3)','# model uncertainty (2&sigma;) ',sig_mod
        end if
    else if ( lweighted ) then
        print '(a,g12.3)','# model uncertainty is neglected'
    end if
    if ( nobs == 1 ) then
        icolor = 2
    else
        icolor = 1
    end if
    do i=1,nobs
        call printsynline(yr1,yr2,data(1,i),icolor,names(i))
    end do
    if ( .not. lnoave .and. nobs > 1 ) call printsynline(yr1,yr2,obs,2,'observations')
    do i=nobs+1,n
        call printsynline(yr1,yr2,data(1,i),3,names(i))
    end do
    if ( .not. lnoave ) then
        call printsynline(yr1,yr2,mod,4,'models')
        if ( lweighted ) then
            call printsynline(yr1,yr2,syn,5,'synthesis')
            !!!call printsynline(yr1,yr2,syn,4,'{/:Bold synthesis}') gs chokes on his...
        else
            call printsynline(yr1,yr2,syn,5,'average')
        end if
    end if
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
    write(0,*) 'synthesis: error reading sig_mod from argument ',trim(line)
    call exit(-1)
905 continue
    write(0,*) 'synthesis: yr1,yr2 ',yr1,yr2,' out of range ',yrbeg,yrend
    call exit(-1)
906 continue
    write(0,*) 'synthesis: error reading yr1,yr2,x,xlo,xhi,ci from ',trim(line)
    call exit(-1)
999 continue
end program synthesis    

subroutine getsynmean(lweighted,data,n,nobs,sig_mod)
!
!   compute the mean of model data, either weighted with model error or unweighted
!   n is the number of input points, ..(n) denotes the average.
!
    implicit none
    logical,intent(in) :: lweighted
    integer,intent(in) :: n,nobs
    real,intent(in) :: sig_mod
    real,intent(inout) :: data(5,n+1)
    integer :: i,j
    real s1,ss2(3),w,w1
!
    s1 = 0
    ss2 = 0
    w1 = 0
    do i=nobs+1,n
        if ( lweighted ) then
            if ( abs(data(3,i) - data(2,i)) < 1e-6 ) then
                write(0,*) 'synthesis: error: upper and lower bound are equal: ',i,data(2,i),data(3,i)
                call exit(-1)
            end if
            if ( sig_mod > 0 ) then
                w = 1/((data(3,i) - data(2,i))**2 + sig_mod**2)
            else
                w = 1/(data(3,i) - data(2,i))**2
            end if
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
        ss2(j) = sqrt(ss2(j))/w1
        if ( sig_mod > 0 ) then
            ss2(j) = sqrt(ss2(j)**2 + sig_mod**2)
        end if
    end do
    data(1,n+1) = s1
    data(2,n+1) = s1 - ss2(2)
    data(3,n+1) = s1 + ss2(3)
end subroutine getsynmean

subroutine getsynchi2(data,n,nobs,sig_mod,chi2)
!
!   compute chi2
!
    implicit none
    integer,intent(in) :: n,nobs
    real,intent(in) :: data(5,n+1)
    real,intent(out) :: sig_mod,chi2
    integer :: i
    real :: s1

    chi2 = 0
    s1 = data(1,n+1)
    do i=nobs+1,n
        if ( sig_mod > 0 ) then
            if ( data(1,i) > s1 ) then
                chi2 = chi2 + (data(1,i)-s1)**2/((data(1,i)-data(2,i))**2 + sig_mod**2)
            else
                chi2 = chi2 + (s1-data(1,i))**2/((data(3,i)-data(1,i))**2 + sig_mod**2)
            end if
        else
            if ( data(1,i) > s1 ) then
                chi2 = chi2 + ((data(1,i)-s1)/(data(1,i)-data(2,i)))**2
            else
                chi2 = chi2 + ((s1-data(1,i))/(data(3,i)-data(1,i)))**2
            end if
        end if
    end do
    chi2 = chi2*4 ! transformed to use 95%~2sigma CIs
end subroutine getsynchi2

real function syn_func(x)
!
!   function to be minimised
!
    implicit none
    real,intent(in) :: x
    real :: sig_mod,chi2
    integer :: syn_n,syn_nobs
    real :: syn_data(5,1000)
    common /syn_dat/ syn_n,syn_nobs,syn_data
    
    sig_mod = x
    call getsynmean(.true.,syn_data,syn_n,syn_nobs,sig_mod)
    call getsynchi2(syn_data,syn_n,syn_nobs,sig_mod,chi2)
    syn_func = chi2 - (syn_n-1)
    !!!print *,'@@@ sig_mod,mean,bounds,chi2 = ',sig_mod,syn_data(1:3,syn_n),chi2,syn_n,syn_nobs
end function syn_func

subroutine printsynline(yr1,yr2,data,colour,name)
    implicit none
    integer,intent(in) :: yr1,yr2,colour
    real,intent(in) :: data(5)
    character,intent(in) :: name*(*)
    integer :: j
    print '(2i5,5g12.3,i3,4a)',yr1,yr2,(data(j),j=1,5),colour,' "',trim(name),'"'
end subroutine