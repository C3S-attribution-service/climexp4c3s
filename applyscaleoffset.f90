subroutine applyscaleoffset(ncid,jvars1,field,nxf,nyf,nzf &
    ,nperyear,yrbeg,yrend,nx,ny,nz,yr1,yr2,lwrite)
    implicit none
    integer :: ncid,jvars1,nxf,nyf,nzf,nperyear,yrbeg,yrend
    integer :: nx,ny,nz,yr1,yr2
    real :: field(nxf,nyf,nzf,nperyear,yrbeg:yrend)
    logical :: lwrite
    integer :: i,j,jx,jy,jz
    real :: scale,offset
    call applyscaleoffsetens(ncid,jvars1,field,nxf,nyf,nzf &
        ,nperyear,yrbeg,yrend,0,nx,ny,nz,yr1,yr2,0,0,lwrite)

end subroutine applyscaleoffset

subroutine applyscaleoffsetens(ncid,jvars1,field,nxf,nyf,nzf &
    ,nperyear,yrbeg,yrend,nensmax,nx,ny,nz,yr1,yr2,mens1,mens,lwrite)

!   there may be a scale and/or offset in the file!

    implicit none
    integer :: ncid,jvars1,nxf,nyf,nzf,nperyear,yrbeg,yrend,nensmax,mens1,mens
    integer :: nx,ny,nz,yr1,yr2
    real :: field(nxf,nyf,nzf,nperyear,yrbeg:yrend,0:nensmax)
    logical :: lwrite
    integer :: i,j,jx,jy,jz,iens
    real :: scale,offset

    call getrealattopt(ncid,jvars1,'scale_factor',scale,lwrite)
    call getrealattopt(ncid,jvars1,'add_offset',offset,lwrite)
!       returns 3e33 if not found...
    if ( offset > 1e33 ) offset = 0
    if ( scale < 1e33 .and. scale /= 1 ) then
        if ( offset /= 0 ) then
!           going through all the data takes an awful lot of time,
!           so do not do it twice...
            if ( lwrite ) print *,'readncfile: scaling and offsetting with ',scale,offset
            call keepalive1('Applying scale and offset ',0,2)
            do iens=mens1,mens
                do i=yr1,yr2
                    do j=1,nperyear
                        do jz=1,max(nz,1)
                            do jy=1,max(ny,1)
                                do jx=1,max(nx,1)
                                    if ( field(jx,jy,jz,j,i,iens) < 2e33 ) then
                                        field(jx,jy,jz,j,i,iens) = &
                                            field(jx,jy,jz,j,i,iens)*scale + offset
                                    end if
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        else
            if ( lwrite ) print *,'readncfile: scaling with ',scale
            call keepalive1('Applying scale ',1,2)
            do iens=mens1,mens
                do i=yr1,yr2
                    do j=1,nperyear
                        do jz=1,max(nz,1)
                            do jy=1,max(ny,1)
                                do jx=1,max(nx,1)
                                    if ( field(jx,jy,jz,j,i,iens) < 2e33 ) then
                                        field(jx,jy,jz,j,i,iens) = &
                                            field(jx,jy,jz,j,i,iens)*scale
                                    end if
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        endif
    else if ( offset /= 0 ) then
        if ( lwrite ) print *,'readncfile: offsetting ',offset
        call keepalive1('Applying offset ',2,2)
        do iens=mens1,mens
            do i=yr1,yr2
                do j=1,nperyear
                    do jz=1,max(nz,1)
                        do jy=1,max(ny,1)
                            do jx=1,max(nx,1)
                                if ( field(jx,jy,jz,j,i,iens) < 2e33 ) then
                                    field(jx,jy,jz,j,i,iens) = &
                                        field(jx,jy,jz,j,i,iens) + offset
                                end if
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end if
end subroutine applyscaleoffsetens

