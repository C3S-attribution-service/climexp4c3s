subroutine getcorr(series1,nmax1,nper1,fy1,ly1,series2,nmax2,nper2,fy2,ly2,j1,j2,corr,lwrite)
!
!   Compute the correlation between the two series.
!   Only the correlation, no significances or anything.
!
    implicit none
    integer nmax1,nper1,fy1,ly1,nmax2,nper2,fy2,ly2,j1,j2
    real series1(nmax1,fy1:ly1),series2(nmax2,fy2:ly2),corr
    logical lwrite
    integer yyr,yr,mm,mo,ntot,nmax,fyr,lyr,nperyear
    real prob,z,ax,sxx,ay,syy,sxy,df
    real,allocatable :: xx(:),yy(:)
    
    fyr = max(fy1,fy2)
    lyr = min(ly1,ly2)
    if ( lyr.lt.fyr .or. j2.lt.j1 ) then
        corr = 3e33
        return
    end if
    if ( nper1 /= nper2 ) then
        write(0,*) 'getcorr: error: cannot handle unequal time scales yet',nper1,nper2
        call exit(-1)
    end if
    nperyear = min(nper1,nper2)
    nmax = nperyear*(j2-j1+1)*(lyr-fyr+1)
    if ( lwrite ) then
        print *,'nmax = ',nmax,nperyear,j1,j2,fyr,lyr
        if ( .true. ) then
            do yr=fyr,lyr
                if ( series1(j1,yr).lt.1e33 .and. &
     &               series2(j1,yr).lt.1e33 ) then
                    print *,yr,series1(j1,yr),series2(j1,yr)
                end if
            end do
        end if
    end if
    allocate(xx(nmax),yy(nmax))
    ntot = 0
    do yyr=fyr,lyr
        do mm=j1,j2
            mo = mm
            call normon(mo,yyr,yr,nperyear)
            if ( yr.ge.fyr .and. yr.le.lyr ) then
                if ( series1(mo,yr).lt.1e33 .and. series2(mo,yr).lt.1e33 ) then
                    ntot = ntot + 1
                    if ( ntot > nmax ) then
                        write(0,*) 'getcorr: internal error: ntot > nmax ',ntot,nmax
                        call exit(-1)
                    end if
                    xx(ntot) = series1(mo,yr)
                    yy(ntot) = series2(mo,yr)
                end if
            end if
        end do
    end do
    df = ntot - 2 ! not used
    call pearsnxx(xx,yy,ntot,corr,prob,z,ax,sxx,ay,syy,sxy,df)
end subroutine