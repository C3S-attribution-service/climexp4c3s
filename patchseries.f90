program patchseries
!
!       fill in missing values in a time series or extend it in time using
!       data froma  second time series. A linear regression on the overlap is 
!       used to adjust the second series to the first one.
!
    implicit none
    include 'param.inc'
    integer :: i,j,dy,mo,yr,nperyear,nperyear2,n,nperday,dpm(12),nmetadata
    real,allocatable :: maindata(:,:),auxdata(:,:)
    real :: xx(30*(yrend-yrbeg+1)),yy(30*(yrend-yrbeg+1))
    real :: scale(12),offset(12),a,b,siga,sigb,chi2,q,sig(1),sx(12),sy(12)
    character :: var*100,units*120,var2*100,units2*20,mainfile*1023,auxfile*1023
    character :: method*4,lvar*200,svar*100,lvar2*200,svar2*100, &
        metadata(2,100)*2000,metadata2(2,100)*2000,history*10000,history2*10000
    logical :: lwrite,lreversefit,lstandardunits
    data dpm /31,29,31,20,31,30,31,31,30,31,30,31/

    lwrite = .false.
    lreversefit = .false.
    lstandardunits = .true.
    sx = 0
    sy = 0
    if ( command_argument_count() < 2 ) then
        write(0,*) 'usage: patchseries mainseries auxseries [regr|bias|none]'
        write(0,*) 'patches holes in mainseries using data from '// &
             'auxseries linearly regressed on mainseries | biascorrected to mainseries '// &
             '| without corrections (default)'
        stop
    endif
    allocate(maindata(npermax,yrbeg:yrend),auxdata(npermax,yrbeg:yrend))
    call get_command_argument(1,mainfile)
    if ( lwrite ) print *,'reading ',trim(mainfile)
    call readseriesmeta(mainfile,maindata,npermax,yrbeg,yrend,nperyear &
         ,var,units,lvar,svar,history,metadata,lstandardunits,lwrite)
    call get_command_argument(2,auxfile)
    if ( lwrite ) print *,'reading ',trim(auxfile)
    call readseriesmeta(auxfile,auxdata,npermax,yrbeg,yrend,nperyear2  &
         ,var2,units2,lvar2,svar2,history2,metadata2,lstandardunits,lwrite)
    if ( nperyear.ne.nperyear2 ) then
        write(0,*) 'patchseries: error: can only handle series '// &
             'with the same time resolution, not ',nperyear,nperyear2
        write(*,*) 'patchseries: error: can only handle series '// &
             'with the same time resolution, not ',nperyear,nperyear2
        call exit(-1)
    end if
    if ( units /= units2 ) then
        write(0,*) 'patchseries: warning: unequal units ',trim(units),' vs ',trim(units2)
    end if
    if ( command_argument_count() >= 3 ) then
        call get_command_argument(3,method)
        ! "noscale" was recognised by a version that I mistakenly skipped
        if ( method == 'nosc' ) method = 'bias'
    else
        method = 'none'
    end if
!
!   determine regression coefficients per month or less
!
    if ( nperyear.le.12 ) then
        sx = 0
        sy = 0
        do mo=1,nperyear
            n = 0
            do yr=yrbeg,yrend
                if ( maindata(mo,yr).lt.1e33 .and. &
 &                   auxdata(mo,yr).lt.1e33 ) then
                    n = n + 1
                    xx(n) = auxdata(mo,yr)
                    sx(mo) = sx(mo) + xx(n)
                    yy(n) = maindata(mo,yr)
                    sy(mo) = sy(mo) + yy(n)
                end if
            end do
            if ( n.gt.0 ) then
                sx(mo) = sx(mo)/n
                sy(mo) = sy(mo)/n
            end if
            if ( method == 'none' ) then
                scale(mo) = 1
                offset(mo) = 0
            else if ( method == 'bias' ) then
                if ( n.lt.7 ) then
                    scale(mo) = 3e33
                    offset(mo) = 3e33
                else
                    scale(mo) = 1
                    offset(mo) = sy(mo) - sx(mo)
                end if
            else if ( method == 'regr' ) then
                if ( n.lt.7 ) then ! arbitrary
                    scale(mo) = 3e33
                    offset(mo) = 3e33
                else if ( lreversefit ) then
                    call fit(xx,yy,n,sig,0,offset(mo),scale(mo),siga,sigb,chi2,q)
                else
                    call fit(yy,xx,n,sig,0,a,b,siga,sigb,chi2,q)
                    scale(mo) = 1/b
                    offset(mo) = -a/b
                end if
            else
                write(0,*) 'patchseries: error: unknown method ',method
                call exit(-1)
            end if
        end do
        do yr=yrbeg,yrend
            do mo=1,nperyear
                if ( maindata(mo,yr).gt.1e33 .and. &
 &                   auxdata(mo,yr).lt.1e33 ) then
                    maindata(mo,yr) = offset(mo) + scale(mo)*auxdata(mo,yr)
                end if
            end do
        end do
    else if ( nperyear.ge.260 ) then ! daily or 6-hourly frequency
        if ( nperyear.eq.360 .or. nperyear.eq.4*360 ) then
            dpm = 30
            nperday = nperyear/360
        else if ( nperyear.eq.365 .or. nperyear.eq.4*365 ) then
            dpm(2) = 28
            nperday = nperyear/365
        else
            nperday = nint(nperyear/366.)
        end if
        sx = 0
        sy = 0
        do mo=1,12
            n = 0
            do yr=yrbeg,yrend
                do dy=1,nperday*dpm(mo)
                    call invgetdymo(dy,mo,j,nperyear)
                    if ( maindata(j,yr).lt.1e33 .and. &
 &                       auxdata(j,yr).lt.1e33 ) then
                        n = n + 1
                        xx(n) = auxdata(j,yr)
                        sx(mo) = sx(mo) + xx(n)
                        yy(n) = maindata(j,yr)
                        sy(mo) = sy(mo) + yy(n)
                    end if
                end do
            end do
            if ( n.gt.0 ) then
                sx(mo) = sx(mo)/n
                sy(mo) = sy(mo)/n
            end if
            if ( method == 'none' ) then
                scale(mo) = 1
                offset(mo) = 0
            else if ( method == 'bias' ) then
                if ( n.lt.7 ) then
                    scale(mo) = 3e33
                    offset(mo) = 3e33
                else
                    scale(mo) = 1
                    offset(mo) = sy(mo) - sx(mo)
                end if
            else if ( method == 'regr' ) then
                if ( n.lt.7 ) then ! arbitrary
                    scale(mo) = 3e33
                    offset(mo) = 3e33
                else if ( lreversefit ) then
                    call fit(xx,yy,n,sig,0,offset(mo),scale(mo),siga,sigb,chi2,q)
                else
                    call fit(yy,xx,n,sig,0,a,b,siga,sigb,chi2,q)
                    scale(mo) = 1/b
                    offset(mo) = -a/b
                end if
            else
                write(0,*) 'patchseries: error: unknown method ',method
                call exit(-1)
            end if
        end do
        do yr=yrbeg,yrend
            do mo=1,12
                do dy=1,nperday*dpm(mo)
                    call invgetdymo(dy,mo,j,nperyear)
                    if ( maindata(j,yr).gt.1e33 .and. &
 &                       auxdata(j,yr).lt.1e33 ) then
                        maindata(j,yr) = offset(mo) + scale(mo)*auxdata(j,yr)
                    end if
                end do
            end do
        end do
    else
       write(0,*) 'merging pentad or decadal time series not ready'
       call exit(-1) 
    end if
!
!   output
!
    call printvar(6,var,units,lvar)
    call merge_metadata(metadata,nmetadata,metadata2,' ',history2,'aux_')
    call copyheadermeta(mainfile,6,trim(mainfile)//' extended with '//trim(auxfile),history,metadata)
    if ( method == 'bias' ) then
        write(6,'(a,12f8.3)') '# bias_correction :: ',(offset(mo),mo=1,min(12,nperyear))
    else if ( method == 'regr' ) then
        write(6,'(a,12f8.3)') '# patch_scale :: ',(scale(mo),mo=1,min(12,nperyear))
        write(6,'(a,12f8.3)') '# bias_correction :: ',(offset(mo),mo=1,min(12,nperyear))
    end if
    call printdatfile(6,maindata,npermax,nperyear,yrbeg,yrend)
end program
