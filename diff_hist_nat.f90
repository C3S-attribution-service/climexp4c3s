program dff_hist_nat
!
!   subtracts two asymmetric uncertainty ranges, eg from  historical and natural ensembles of a model
!
!   options: log: assume PDFs look reasonable after transforming to logarithm
!            perc: input is a percentage, use 1+p/100.
!   input:   file or stdin: X Xlo [...] Xhi CI. At the moment only CI=95 is supported.
!
!   At the momengt just assume two half-gaussians with differing widths.
!   Will be improved with better algorithms (thanks for the hints, David S).
!
    implicit none
    integer   :: i,j
    real      :: xx(3,2),diff(3),ci(2)
    logical   :: logscale,lperc,lexist
    character :: file*1024,string*128
    integer   :: iargc

    logscale = .false.
    lperc = .false.
    if ( iargc() > 0 ) then
        call getarg(1,file)
        inquire(file=trim(file),exist=lexist)
    else
        lexist = .false.
    end if
    if ( .not.lexist .or. file == '-' ) then
        print *,'Hist X Xlo Xhi CI? '
        read '(a)',string
        call rm_ellipsis(string)
        read(string,*) (xx(j,1),j=1,3),ci(1)
        print *,'Nat  X Xlo Xhi CI? '
        read '(a)',string
        call rm_ellipsis(string)
        read(string,*) (xx(j,2),j=1,3),ci(2)
    else
        open(1,file=trim(file),status='old')
        i = 1
10      continue
        read(1,'(a)',end=20) string
        if ( string == ' ' .or. index(string(1:2),'#') /= 0 ) goto 10
        call rm_ellipsis(string)
        read(string,*) (xx(j,i),j=1,3),ci(i)
        i = i+1
        if ( i > 2 ) goto 20
        goto 10
20      continue
    end if
    do i=1,iargc()
        if ( lexist .and. i == 1 ) cycle
        call getarg(i,string)
        if ( string(1:3) == 'log' ) logscale = .true.
        if ( string(1:4) == 'perc' ) lperc = .true.
    end do
  
    if ( ci(1) /= ci(2) ) then
        write(0,*) 'dff_hist_nat: error: can only handle equal CI , not ',ci
        call exit(-1)
    end if
    if ( lperc ) then
        xx = 1 + xx/100
    end if
    if ( logscale ) then
        do i=1,2
            do j=1,3
                if ( xx(j,i) > 0 ) then
                    xx(j,i) = log(xx(j,i))
                else
                    write(0,*) 'dff_hist_nat: cannpot use logscale for zero or negative data ',xx(j,i)
                    call exit(-1)
                end if
            end do
        end do
    end if
   
    diff(1) = xx(1,1) - xx(1,2)
    !!!print *,'@@@ ',xx(1,1),' - ',xx(1,2),' = ',diff(1),exp(diff(1))
    diff(2) = diff(1) - sqrt( (xx(1,1)-xx(2,1))**2 + (xx(1,2)-xx(3,2))**2 ) ! lower bound of hist + upper one of nat
    !!!print *,'@@@ ',diff(1),' - ',diff(1)-diff(2),' = ',diff(2),exp(diff(2))
    diff(3) = diff(1) + sqrt( (xx(1,1)-xx(3,1))**2 + (xx(1,2)-xx(2,2))**2 ) ! upper bound of hist + lower one of nat
    !!!print *,'@@@ ',diff(1),' + ',diff(3)-diff(1),' = ',diff(3),exp(diff(3))

    if ( logscale ) then
        diff = exp(diff)
    end if
    if ( lperc ) then
        diff = 100*(diff-1)
    end if
    
    print *,diff,ci(1)

end program

subroutine rm_ellipsis(string)
    implicit none
    character :: string*(*)
    integer   :: i

!   there probably is a standard library routine for this now
    do i=1,len(string)-3
        if ( string(i:i+2) == '...' ) string(i:i+2) = '   '
    end do
end subroutine