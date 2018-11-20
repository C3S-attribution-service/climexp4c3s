subroutine printtable(fcst,obs,npermax,yrbeg,yrend,nmodel1 &
    ,nmodel2,nens1,nens2,yr1,yr2,m1,m2,nperyear,n &
    ,yrstart,yrstop)
    implicit none
    integer :: npermax,yrbeg,yrend,nmodel1,nmodel2,nens1,nens2,yr1,yr2 &
        ,m1,m2,nperyear,n,yrstart,yrstop
    real :: obs(npermax,yrbeg:yrend),fcst(npermax,yrbeg:yrend,0:nens2)
    integer :: nmodel,yr,month,mo1,mo,iens,j1,j2
    logical :: lvalid
    character :: line*80

    integer,save :: init=0

    if ( nmodel1 == -1 ) then
    
!       obs + ensemble
    
        if ( init == 0 ) then
            init = 1
            call printtableheader(nperyear,nens1,nens2)
        endif
        do yr=yr1,yr2
            do month=m1,m2
                call getj1j2(j1,j2,month,nperyear, .false. )
                do mo1=j1,j2
                    mo = mo1
                    if ( mo > 12 ) mo = mo - 12
                    if ( .true. ) then
!                       demand at least one valid forecast
                        lvalid = .false. 
                        do iens=nens1,nens2
                            if ( fcst(mo,yr,iens) < 1e33 ) lvalid = .true. 
                        enddo
                    else
                        lvalid = .true. 
                    endif
                    do iens=nens1,nens2
                        if ( fcst(mo,yr,iens) > 1e33 ) then
                            fcst(mo,yr,iens) = -999.9
                        endif
                    enddo
!                   and a valid observation
                    if ( lvalid .and. obs(mo,yr) < 1e30 ) then
                        n = n + 1
                        yrstart = min(yrstart,yr)
                        yrstop  = max(yrstop,yr)
                        write(10,'(i4,i5,100g14.6)') yr,mo,obs(mo,yr),(fcst(mo,yr,iens),iens=nens1,nens2)
                    endif
                enddo
            enddo
        enddo
    else
    
!       perfect model
    
        if ( init == 0 ) then
            init = 1
            call printtableheader(nperyear,nens1,nens2-1)
        endif
        n = 0
        do nmodel=nmodel1,nmodel2
            do yr=yr1,yr2
                do month=m1,m2
                    call getj1j2(j1,j2,month,nperyear, .false. )
                    do mo1=j1,j2
                        mo = mo1
                        if ( mo > 12 ) mo = mo - 12
!                       demand at least one valid forecast
                        lvalid = .false. 
                        do iens=nens1,nmodel-1
                            if ( fcst(mo,yr,iens) < 1e33 ) lvalid = .true. 
                        enddo
                        do iens=nmodel+1,nens2
                            if ( fcst(mo,yr,iens) < 1e33 ) lvalid = .true. 
                        enddo
!                       and a valid perfect model
                        if ( lvalid .and. fcst(mo,yr,nmodel) < 1e30 ) then
                            n = n + 1
                            write(10,'(i4,i5,100g14.6)') yr,mo &
                                ,fcst(mo,yr,nmodel) &
                                ,(fcst(mo,yr,iens),iens=nens1,nmodel-1) &
                                ,(fcst(mo,yr,iens),iens=nmodel+1,nens2)
                        endif
                    enddo
                enddo
            enddo
        enddo
    endif
end subroutine printtable
