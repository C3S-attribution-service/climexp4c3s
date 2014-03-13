subroutine readonencfield(ncid,jvars,lsmask,nxls,nyls,lwrite)
  !
  ! read a single X-Y field from a netcdf file
  !
  implicit none
  include 'netcdf.inc'
  integer ncid,jvars(6,1),nxls,nyls
  real lsmask(nxls,nyls)
  logical lwrite
  integer status,jx,jy
  real undef,x,scale,offset
  
  status = nf_get_var_real(ncid,jvars(1,1),lsmask)
  if ( status.ne.nf_noerr ) call handle_err(status,' reading one netcdf field')
  !       
  ! change undefs to our convention (3e33)!
  !
  call getrealattopt(ncid,jvars(1,1),'_FillValue',undef,lwrite)
  call getrealattopt(ncid,jvars(1,1),'missing_value',x,lwrite)
  if ( undef.eq.3e33 .and. x.ne.3e33 ) undef = x
  if ( lwrite ) print *,'readonencfield: found undef = ',undef
  if ( undef.ne.3e33 ) then
     if ( lwrite ) print *,'changing undef from ',undef,' to 3e33'
     if ( abs(undef).gt.1000 ) then
        ! relative diff
        do jy=1,max(nyls,1)
           do jx=1,max(nxls,1)
              if ( abs(lsmask(jx,jy)-undef).lt.1e-6*abs(undef) ) then
                 lsmask(jx,jy) = 3e33
              end if
           end do
        end do
     else
        ! absolute diff
        do jy=1,max(nyls,1)
           do jx=1,max(nxls,1)
              if ( abs(lsmask(jx,jy)-undef).lt.1e-5 ) then
                 lsmask(jx,jy) = 3e33
              end if
           end do
        end do
     end if
  end if
  !
  ! there may be a scale and/or offset in the file!
  !
  call getrealattopt(ncid,jvars(1,1),'scale_factor',scale,lwrite)
  call getrealattopt(ncid,jvars(1,1),'add_offset',offset,lwrite)
  ! returns 3e33 if not found...
  if ( offset.gt.1e33 ) offset = 0
  if ( scale.gt.1e33 ) scale = 1
  if ( offset.ne.0  .or. scale.ne.1 ) then
     if ( lwrite ) print *,'readncfile: scaling and offsetting with ',scale,offset
     call keepalive(0,0)
     do jy=1,max(nyls,1)
        do jx=1,max(nxls,1)
           if ( lsmask(jx,jy).lt.2e33 ) then
              lsmask(jx,jy) = lsmask(jx,jy)*scale + offset
           end if
        end do
     end do
  end if
end subroutine readonencfield
