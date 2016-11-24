subroutine getmaskbox(mask,nx,ny,xwrap,cutoff,x1,x2,y1,y2)
!
!   compute a box encompassing all 1s in a mask
!   in: mask(nx,ny)     real 0=no data, 1 =data, fractional values possible
!       xwrap           if .true. the field wraps around the earth
!       cutoff          include points with mask>=cutoff
!   out: x1,x2,y1,y1    corners of box (possible the whole field)
!
    implicit none
    integer :: nx,ny,x1,x2,y1,y2
    real    :: mask(nx,ny),cutoff
    logical :: xwrap
    integer :: ix,ixx,iy
    
    ! if a mask, find encompassing box
    ! latitude is easy
    y1 = 0
    do iy=1,ny
        do ix=1,nx
            if ( lsmask(ix,iy) > 0.5 .and. y1 == 0 ) y1 = iy
        end do
    end do
    y2 = ny
    do iy=ny,1,-1
        do ix=1,nx
            if ( lsmask(ix,iy) > 0.5 .and. y2 == ny ) y2 = iy
        end do
    end do
    if ( .not.xwrap ) then
        ! longitude also if it does not wrap
        x1 = 0
        do ix=1,nx
            do iy=1,ny
                if ( lsmask(ix,iy) > 0.5 .and. x1 == 0 ) x1 = ix
            end do
        end do
        x2 = 0
        do ix=nx,1,-1
            do iy=1,ny
                if ( lsmask(ix,iy) > 0.5 .and. x2 == 0 ) x2 = ix
            end do
        end do
        if ( x1 == 0 .or. x2 == 0 ) then
            x1 = 0
            x2 = 0
            y1 = 0
            y2 = 0
            return
        end if
    else ! lwrap
        ! first search for outside the mask, next inside.
        x1 = -999
        do ix=1,nx
            if ( x1 == -999 ) then
                allzero = .true.
                do iy=1,ny
                    if ( lsmask(ix,iy) > 0.5 ) allzero = .false.
                end do
                if ( allzero ) x1 = 0
            end if
            if ( ix == 0 ) then
                do iy=1,ny
                    if ( lsmask(ix,iy) > 0.5 ) x1 = ix
                end do
            end if
        end do
        if ( x1 == -999 ) x1 = 1
        if ( x1 == 0 ) then
            x1 = 0
            x2 = 0
            y1 = 0
            y2 = 0
            return
        end if
        x2 = -999
        do ix=nx,1,-1
            if ( x2 == -999 ) then
                allzero = .true.
                do iy=1,ny
                    if ( lsmask(ix,iy) > 0.5 ) allzero = .false.
                end do
                if ( allzero ) x2 = 0
            end if
            if ( x2 == 0 ) then
                do iy=1,ny
                    if ( lsmask(ix,iy) > 0.5 ) x2 = ix
                end do
            end if
        end do            
        if ( x2 == -999 ) x2 = nx
        if ( x2 == 0 ) then
            x1 = 0
            x2 = 0
            y1 = 0
            y2 = 0
            return
        end if
        if ( x2 < x1 ) x2 = x2 + nx            
    end if ! lwrap
end subroutine