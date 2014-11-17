subroutine checkgridequal(nx,ny,xx,yy,nx1,ny1,xx1,yy1)
    integer nx,ny,nx1,ny1
    real xx(nx),yy(ny),xx1(nx1),yy1(ny1)
    integer nz,nz1
    real zz(1),zz1(1)
    nz = 1
    zz(1) = 0
    nz1 = 1
    zz1(1) = 0
    call checkgridequal3d(nx,ny,nz,xx,yy,zz,nx1,ny1,nz1,xx1,yy1,zz1)
end subroutine
subroutine checkgridequal3d(nx,ny,nz,xx,yy,zz,nx1,ny1,nz1,xx1,yy1,zz1)
    implicit none
    integer nx,ny,nz,nx1,ny1,nz1
    real xx(nx),yy(ny),zz(nz),xx1(nx1),yy1(ny1),zz1(nz1)
    logical lequal
    call aregridsequal3d(nx,ny,nz,xx,yy,zz,nx1,ny1,nz1,xx1,yy1,zz1,lequal,.false.)
    if ( .not.lequal ) then
        call aregridsequal3d(nx,ny,nz,xx,yy,zz,nx1,ny1,nz1,xx1,yy1,zz1,lequal,.true.)
        write(0,*) 'checkgridequal3d: error: unequal grids '
        write(*,*) 'checkgridequal3d: error: unequal grids '
        call exit(-1)
    endif
end subroutine
subroutine aregridsequal(nx,ny,xx,yy,nx1,ny1,xx1,yy1,lequal,lwrite)
!
!   check that the 2D grids are equal
!
    implicit none
    integer nx,ny,nx1,ny1
    real xx(nx),yy(ny),xx1(nx1),yy1(ny1)
    logical lequal,lwrite
    integer nz,nz1
    real zz(1),zz1(1)
    nz = 1
    zz(1) = 0
    nz1 = 1
    zz1(1) = 0
    call aregridsequal3d(nx,ny,nz,xx,yy,zz,nx1,ny1,nz1,xx1,yy1,zz1,lequal,lwrite)
end subroutine
subroutine aregridsequal3d(nx,ny,nz,xx,yy,zz,nx1,ny1,nz1,xx1,yy1,zz1,lequal,lwrite)
!
!   check that the 3D grids are equal
!
    implicit none
    integer nx,ny,nz,nx1,ny1,nz1
    real xx(nx),yy(ny),zz(nz),xx1(nx1),yy1(ny1),zz1(nz1)
    logical lequal,lwrite
    integer i,j
!
    lequal = .true.
    if ( nx1.ne.nx ) then
        lequal = .false.
        if ( lwrite ) print *,'aregridsequal3d: nx1 != nx: ',nx1,nx
        return
    endif
    if ( ny1.ne.ny ) then
        lequal = .false.
        if ( lwrite ) print *,'aregridsequal3d: ny1 != ny: ',ny1,ny
        return
    endif
    if ( nz1.ne.nz ) then
        lequal = .false.
        if ( lwrite ) print *,'aregridsequal3d: nz1 != nz: ',nz1,nz
        return
    endif
    do i=1,nx
        if ( abs(xx(i)-xx1(i)).gt.0.01 ) then
            lequal = .false.
            if ( lwrite ) print *,'aregridsequal: xx != xx1: ',i,xx(i),xx1(i)
            return
        endif
    enddo
    do i=1,ny
        if ( abs(yy(i)-yy1(i)).gt.0.01 ) then
            lequal = .false.
            if ( lwrite ) print *,'aregridsequal: yy != yy1: ',i,yy(i),yy1(i)
            return
        endif
    enddo
    do i=1,nz
        if ( abs(zz(i)-zz1(i)).gt.0.01 ) then
            lequal = .false.
            if ( lwrite ) print *,'aregridsequal: zz != zz1: ',i,zz(i),zz1(i)
            return
        endif
    enddo
end subroutine
