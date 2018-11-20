subroutine applylsmask(field,lsmask,nx,ny,nz,nperyear,yr1,yr2,nens1,nens2,lsmasktype,lwrite)
  !
  !       apply the land sea mask according to the masking in lsmasktype
  !       (indicates what should be defined, the rest is undefined)
  !
  implicit none
  integer nx,ny,nz,nperyear,yr1,yr2,nens1,nens2
  real field(nx,ny,nz,nperyear,yr1:yr2,nens1:nens2),lsmask(nx,ny)
  character lsmasktype*(*)
  logical lwrite
  integer i,j,k,mo,yr,iens
!
  if ( lwrite ) print *,'applylsmask: lsmasktype = ',lsmasktype
  if ( lsmasktype.eq.'all' ) return
  if ( lsmasktype.eq.'sea' .or. lsmasktype.eq.'nots' ) then
     do iens=nens1,nens2
        do yr=yr1,yr2
           do mo=1,nperyear
              do k=1,nz
                 do j=1,ny
                    do i=1,nx
                       if ( abs(lsmask(i,j)).gt.1e-4 .eqv. lsmasktype.eq.'sea' ) then
                          field(i,j,k,mo,yr,iens) = 3e33
                       end if
                    end do
                 end do
              end do
           end do
        end do
     enddo
  else if ( lsmasktype.eq.'land' .or. lsmasktype.eq.'notl' ) then
     do iens=nens1,nens2
        do yr=yr1,yr2
           do mo=1,nperyear
              do k=1,nz
                 do j=1,ny
                    do i=1,nx
                       if ( abs(lsmask(i,j)-1).gt.1e-4 .eqv. lsmasktype.eq.'land' ) then
                          field(i,j,k,mo,yr,iens) = 3e33
                       end if
                    end do
                 end do
              end do
           end do
        end do
     enddo
  else if ( lsmasktype.eq.'5lan' .or. lsmasktype.eq.'5sea' ) then
     do iens=nens1,nens2
        do yr=yr1,yr2
           do mo=1,nperyear
              do k=1,nz
                 do j=1,ny
                    do i=1,nx
                       if ( lsmask(i,j).lt.0.5 .eqv. lsmasktype.eq.'5lan' ) then
                          field(i,j,k,mo,yr,iens) = 3e33
                       end if
                    end do
                 end do
              end do
           end do
        end do
     enddo
  else
     write(0,*) 'applylsmask: error: donot know what to do for lsmasktype = ',lsmasktype
     call exit(-1)
  end if
end subroutine applylsmask
