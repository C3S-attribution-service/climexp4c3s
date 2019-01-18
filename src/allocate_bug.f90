program fail
  implicit none
  integer*8 :: ix,iy,dy,yr
  integer*8 :: status
  integer*8 :: nx,ny,nperyear,firstyr,lastyr
  real*4, allocatable :: field(:,:,:,:)
 
  nx = 464
  ny = 201
  nperyear = 366
  firstyr= 1950
  lastyr = 2012

  print *,'n = ',nx*ny*nperyear*(lastyr-firstyr+1)
 
  allocate(field(nx,ny,nperyear,firstyr:lastyr),&
           stat=status)
  print *,"status = ",status
  if (status .ne. 0) then
     print *,"Allocation failed ..."
     stop 1
  end if

  do yr=firstyr,lastyr
     print *,"yr = ",yr
     do dy=1,nperyear
        do iy=1,ny
           do ix=1,nx
              field(ix,iy,dy,yr) = 1000*yr + dy + iy/1000. + ix/1000000.
           end do
        end do
     end do
  end do
  do yr=firstyr,lastyr
     do dy=1,nperyear
        do iy=1,ny
           do ix=1,nx
              print *,'field(',ix,iy,dy,yr,') = ',field(ix,iy,dy,yr)
           end do
        end do
     end do
  end do

end program fail
