subroutine getreftime(ncid,varid,tt,nt,firstmo,firstyr &
    ,nperyear,iperyear,lwrite)

!   extract the reference times from the netcdf file

    implicit none
    include 'netcdf.inc'
    integer,intent(in) :: ncid,varid,nt,firstmo,firstyr,nperyear,iperyear
    real*8,intent(out) :: tt(nt)
    logical,intent(in) :: lwrite
    integer :: it,n,status
    real*8,allocatable :: ttu(:)
    logical,allocatable :: tdefined(:)
    character :: units*(nf_max_name),ltime*40

    if ( lwrite ) print *,'retrieving reftime',ncid,varid,nt
!   read data
    allocate(ttu(nt))
    status = nf_get_var_double(ncid,varid,tt)
    if ( status /= nf_noerr ) call handle_err(status,'nf_get_var_int reftimes')
    if ( lwrite ) then
        do it=1,nt
            print *,it,tt(it)
        enddo
    endif
!   get unique time steps
    n = 1
    ttu(1) = tt(1)
    if ( lwrite ) print *,'found new ttu(',1,') = ',ttu(1)
    do it=2,nt
        if ( tt(it) > ttu(n) ) then
            n = n + 1
            ttu(n) = tt(it)
            if ( lwrite ) print *,'found new ttu(',n,') = ',ttu(n)
        endif
    enddo
!   convert to characters - only the month, I am
!   not interested in the years.
    allocate(tdefined(n))
    call getperyear(ncid,varid,ttu,n,firstmo,firstyr,nperyear,iperyear,ltime,tdefined,n,lwrite)
    deallocate(ttu)
    deallocate(tdefined)
    if ( lwrite ) then
        print *,'getreftime: firstmo,firstyr = ',firstmo,firstyr
        print *,'getreftime: nperyear,iperyear = ',nperyear,iperyear
    end if
end subroutine getreftime

subroutine getleadtime(ncid,varid,tt,nt,leadunits,lwrite)

!   extract the lead times from the netcdf file

    implicit none
    include 'netcdf.inc'
    integer,intent(in) :: ncid,varid,nt
    real*8,intent(out) :: tt(nt)
    character,intent(out) :: leadunits*(*)
    logical,intent(in) :: lwrite
    integer :: it,n,status

    if ( lwrite ) print *,'retrieving leadtime',ncid,varid
!   get units attribute
    call gettextatt(ncid,varid,'units',leadunits,lwrite)
    call tolower(leadunits)
!   read data
    status = nf_get_var_double(ncid,varid,tt)
    if ( status /= nf_noerr ) call handle_err(status,'nf_get_var_int leadtimes')
    if ( lwrite ) then
        do it=1,nt
            print *,it,tt(it)
        enddo
    endif
end subroutine getleadtime
