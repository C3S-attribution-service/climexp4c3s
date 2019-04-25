subroutine moment(data,n,ave,adev,sdev,var,skew,curt)
!
!   replacement for Numerical recipes routine
!   very inefficient but computers are fast these days...
!
    use fgsl
    implicit none
    integer :: n
    real :: data(n),ave,adev,sdev,var,skew,curt
    integer(fgsl_size_t) :: dim
    real(fgsl_double) :: dave,dvar,dsdev,dadev,dskew,dcurt
    real(fgsl_double), allocatable :: ddata(:)
    logical,parameter :: lwrite=.false.

    if ( lwrite ) then
        print *,'entering moment, n = ',n
        print *,'data = ',data(1:n)
    endif
    ave = 3e33
    adev = 3e33
    sdev = 3e33
    var = 3e33
    skew = 3e33
    curt = 3e33
 !
!   call FGSL routines
!    
    allocate(ddata(n))
    ddata = data ! C double
    dim = n ! correct data type I think
    if ( n < 1 ) return
    dave = fgsl_stats_mean(ddata,1_fgsl_size_t,dim)
    ave = dave
    if ( n < 2 ) return    
    dvar = fgsl_stats_variance_m(ddata,1_fgsl_size_t,dim,dave)
    var = dvar
    dsdev = sqrt(dvar)
    sdev = dsdev
    dadev = fgsl_stats_absdev_m(ddata,1_fgsl_size_t,dim,dave)
    adev = dadev
    if ( n < 3 .or. dsdev == 0 ) return
    dskew = fgsl_stats_skew_m_sd(ddata,1_fgsl_size_t,dim,dave,dsdev)
    skew = dskew
    if ( n < 4 ) return
    dcurt = fgsl_stats_kurtosis_m_sd(ddata,1_fgsl_size_t,dim,dave,dsdev)
    curt = dcurt
    deallocate(ddata)!

end subroutine moment