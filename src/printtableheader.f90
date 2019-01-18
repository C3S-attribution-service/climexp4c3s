subroutine printtableheader(nperyear,nens1,nens2)
    implicit none
    integer :: nperyear,nens1,nens2
    integer :: iens
    if ( nperyear == 12 ) then
        write(10,'(a,100(a,i2.2))') '  yr month      obs', &
        ('         ens',iens,iens=nens1,nens2)
    elseif ( nperyear >= 360 .and. nperyear <= 366 ) then
        write(10,'(a,100(a,i2.2))') '  yr  day       obs', &
        ('         ens',iens,iens=nens1,nens2)
    else
        write(10,'(a,100(a,i2.2))') '  yr period     obs', &
        ('         ens',iens,iens=nens1,nens2)
    endif
end subroutine printtableheader
