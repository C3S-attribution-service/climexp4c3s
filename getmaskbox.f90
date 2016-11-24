subroutine getmaskbox(mask,nx,ny,xwrap,cutoff,x1,x2,y1,y2,lwrite)
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
    logical :: xwrap,lwrite
    integer,allocatable :: xwithdata(:)
    integer :: ix,ixx,iy,ix1,xx1,xx2
    logical :: allzero,founddata
    
    ! if a mask, find encompassing box
    ! latitude is easy
    y1 = 0
    do iy=1,ny
        do ix=1,nx
            if ( mask(ix,iy) > cutoff .and. y1 == 0 ) y1 = iy
        end do
    end do
    y2 = ny
    do iy=ny,1,-1
        do ix=1,nx
            if ( mask(ix,iy) > cutoff .and. y2 == ny ) y2 = iy
        end do
    end do
    if ( y1 > ny .or. y2 < 1 ) then
        x1 = 0
        x2 = 0
        y1 = 0
        y2 = 0
        return
    end if
    if ( lwrite ) print *,'getmaskbox: found latitudes ',y1,y2
    if ( .not.xwrap ) then
        ! longitude also if it does not wrap
        x1 = 0
        do ix=1,nx
            do iy=1,ny
                if ( mask(ix,iy) > cutoff .and. x1 == 0 ) x1 = ix
            end do
        end do
        x2 = 0
        do ix=nx,1,-1
            do iy=1,ny
                if ( mask(ix,iy) > cutoff .and. x2 == 0 ) x2 = ix
            end do
        end do
        if ( x1 == 0 .or. x2 == 0 ) then
            if ( lwrite ) print *,'getmaskbox: cannot find longitudes ',x1,x2
            x1 = 0
            x2 = 0
            y1 = 0
            y2 = 0
            return
        end if
        if ( lwrite ) print *,'getmaskbox: found longitudes ',x1,x2
    else ! lwrap
        ! find all the longitudes without data
        allocate(xwithdata(nx))
        founddata = .false.
        do ix=1,nx
            allzero = .true.
            do iy=1,ny
                if ( mask(ix,iy) > cutoff ) allzero = .false.
            end do
            if ( allzero ) then
                xwithdata(ix) = 0
            else
                xwithdata(ix) = 1
                founddata = .true.
            end if
        end do
        if ( .not.founddata ) then
            x1 = 0
            x2 = 0
            y1 = 0
            y2 = 0
            return
        end if
        ! find largest hole
        x1 = 0
        x2 = -1
        xx1 = 0
        xx2 = 0
        do ixx=1,2*nx-1
            ix = ixx
            if ( ix > nx ) ix = ix - nx
            ix1 =ix + 1
            if ( ix1 > nx ) ix1 = ix1 - nx
            if ( xwithdata(ix) == 1 .and. xwithdata(ix1) == 0 ) then
                xx1 = ixx+1
            else if ( xwithdata(ix) == 0 .and. xwithdata(ix1) == 1 ) then
                xx2 = ixx
                if ( xx2-xx1 .gt. x2-x1 ) then
                    if ( lwrite ) print *,'found larger hole,: ',xx1,xx2
                    x1 = xx1
                    x2 = xx2
                end if
            end if
        end do
        ! we now have the hole, invert to get the box around the mask
        xx1 = x1
        xx2 = x2
        x1 = xx2+1
        x2 = xx1-1
        if ( x1 > x2 ) then
            if ( x1 > nx ) then
                x1 = x1 - nx
            else
                x2 = x2 + nx
            end if
        end if
        if ( x1 > nx ) then
            x1 = x1 - nx
            x2 = x2 - nx
        end if
        if ( lwrite ) print *,'getmaskbox: found wrapping longitudes ',x1,x2
    end if ! lwrap
end subroutine