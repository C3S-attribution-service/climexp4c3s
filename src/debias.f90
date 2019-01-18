subroutine debiasmean(obs,fcst,npermax,nperyear,yrbeg,yrend &
    ,yr1,yr2,nens1,nens2,var,lwrite)
    implicit none
    integer,intent(in) :: npermax,nperyear,yrbeg,yrend,yr1,yr2,nens1,nens2
    real,intent(inout) :: obs(npermax,yrbeg:yrend),fcst(npermax,yrbeg:yrend,0:nens2)
    logical,intent(in) :: lwrite
    character,intent(in) :: var*(*)
    integer :: yr,mo,n1,n2,nn2,iens,i,j
    real :: s1,s2,ss1,ss2,zero
    integer,allocatable :: yrzero(:,:)
    real,allocatable :: newfcst(:,:,:),savobs(:,:)
    logical :: isvarprcp
    integer,save :: init=0

    if ( lwrite ) then
        print *,'debiasmean: input  ',var
        print *,'nperyear    = ',nperyear
        print *,'yr1,yr2     = ',yr1,yr2
        print *,'nens1,nens2 = ',nens1,nens2
    endif

    allocate(newfcst(nperyear,yr1:yr2,nens1:nens2))
    if ( isvarprcp(var) ) then
        if ( init == 0 ) then
            init = 1
            print '(a)','# using multiplicative bias correction'
        endif
        allocate(savobs(nperyear,yr1:yr2))
        allocate(yrzero(nperyear,nens1:nens2))
        do yr=yr1,yr2
            do mo=1,nperyear
                savobs(mo,yr) = obs(mo,yr)
            enddo
        enddo
        if ( lwrite ) print *,'debiasmean: taking log10(obs)'
        call takelog(obs(1,yr1),npermax,nperyear,yr1,yr2)
        do iens=nens1,nens2
            do mo=1,nperyear
                yrzero(mo,iens) = -999
                do yr=yr1,yr2
                    if ( fcst(mo,yr,iens) <= 0 ) then
                        yrzero(mo,iens) = yr
                        exit
                    endif
                enddo
            enddo
            if ( lwrite ) print *,'debiasmean: taking log10(fcst) ',iens
            call takelog(fcst(1,yr1,iens),npermax,nperyear,yr1,yr2)
        enddo
    endif
    do mo=1,nperyear
        s1 = 0
        n1 = 0
        do yr=yr1,yr2
            if ( obs(mo,yr) < 1e33 ) then
                n1 = n1 + 1
                s1 = s1 + obs(mo,yr)
            endif
        enddo
        s2 = 0
        n2 = 0
        do iens=nens1,nens2
            do yr=yr1,yr2
                if ( fcst(mo,yr,iens) < 1e33 ) then
                    n2 = n2 + 1
                    s2 = s2 + fcst(mo,yr,iens)
                endif
            enddo
        enddo
        if ( n1 == 0 .or. n2 == 0 ) then
            if ( lwrite ) print *,'debiasmean: no valid data for ',mo
            do iens=nens1,nens2
                do yr=yr1,yr2
                    newfcst(mo,yr,iens) = 3e33
                enddo
            enddo
            goto 800
        endif
        if ( lwrite ) then
            print *,'debiasmean: mean obs = ',mo,s1/n1
            print *,'            mean fcst= ',mo,s2/n2
        endif
    
!       jackknife
    
        do yr=yr1,yr2
            if ( obs(mo,yr) < 1e33 ) then
                if ( n1 == 1 ) then
                    ss1 = 3e33
                else
                    ss1 = (s1 - obs(mo,yr))/(n1-1)
                endif
            else
                ss1 = s1/n1
            endif
            ss2 = s2
            nn2 = n2
            do iens=nens1,nens2
                if ( fcst(mo,yr,iens) < 1e33 ) then
                    ss2 = ss2 - fcst(mo,yr,iens)
                    nn2 = nn2 - 1
                endif
            enddo
            if ( nn2 <= 0 ) then
                ss2 = 3e33
            else
                ss2 = ss2/nn2
            endif
            if ( lwrite ) print *,'debiasmean: bias = ',mo,yr,ss1,ss2,ss1-ss2
            if ( ss1 < 1e33 .and. ss2 < 1e33 ) then
                do iens=nens1,nens2
                    newfcst(mo,yr,iens) = fcst(mo,yr,iens) + ss1 - ss2
                enddo
            else
                do iens=nens1,nens2
                    newfcst(mo,yr,iens) = 3e33
                enddo
            endif
        enddo
    800 continue
    enddo

!   copy back to original arrays, setting the rest to undefined
!   for security reasons

    do iens=nens1,nens2
        do yr=yrbeg,yr1-1
            do mo=1,nperyear
                fcst(mo,yr,iens) = 3e33
            enddo
        enddo
        do yr=yr1,yr2
            do mo=1,nperyear
                fcst(mo,yr,iens) = newfcst(mo,yr,iens)
            enddo
        enddo
        do yr=yr2+1,yrend
            do mo=1,nperyear
                fcst(mo,yr,iens) = 3e33
            enddo
        enddo
    enddo
    if ( isvarprcp(var) ) then
        do yr=yr1,yr2
            do mo=1,nperyear
                obs(mo,yr) = savobs(mo,yr)
            enddo
        enddo
        deallocate(savobs)
        do iens=nens1,nens2
            if ( lwrite ) then
                print *,'debiasmean: taking 10^fcst ',iens
                print *,'yr1,yr2 = ',yr1,yr2
                do i=yr1,yr2
                    print '(i4,12g12.4)',i,(fcst(j,i,iens),j=1,nperyear)
                enddo
            endif
            call takeexp(fcst(1,yr1,iens),npermax,nperyear,yr1,yr2)
            do mo=1,nperyear
                if ( yrzero(mo,iens) /= -999 ) then
                    zero = fcst(mo,yrzero(mo,iens),iens)
                    do yr=yr1,yr2
                        if ( fcst(mo,yr,iens) == zero ) fcst(mo,yr,iens) = 0
                    enddo
                endif
            enddo
        enddo
        deallocate(yrzero)
    endif
    deallocate(newfcst)

end subroutine debiasmean

subroutine debiasvar(obs,fcst,npermax,nperyear,yrbeg,yrend &
    ,yr1,yr2,nens1,nens2,var,lwrite)
    implicit none
    integer,intent(in) :: npermax,nperyear,yrbeg,yrend,yr1,yr2,nens1,nens2
    real,intent(inout) :: obs(npermax,yrbeg:yrend),fcst(npermax,yrbeg:yrend,0:nens2)
    logical,intent(in) :: lwrite
    character,intent(in) :: var*(*)
    integer :: yr,mo,n1,n2,nn2,iens,i,j
    real :: s1,s2,s1m,s2m,t1,t2,ss1,ss2,tt1,tt2,zero
    integer,allocatable :: yrzero(:,:)
    real,allocatable :: newfcst(:,:,:),savobs(:,:)
    logical :: isvarprcp
    integer,save :: init=0,nunlikely=0

    if ( lwrite ) then
        print *,'debiasvar: input  ',var
        print *,'nperyear    = ',nperyear
        print *,'yr1,yr2     = ',yr1,yr2
        print *,'nens1,nens2 = ',nens1,nens2
    endif

    allocate(newfcst(nperyear,yr1:yr2,nens1:nens2))
    if ( isvarprcp(var) ) then
        if ( init == 0 ) then
            init = 1
            print '(a)','# using multiplicative bias correction'
        endif
        allocate(savobs(nperyear,yr1:yr2))
        allocate(yrzero(nperyear,nens1:nens2))
        do yr=yr1,yr2
            do mo=1,nperyear
                savobs(mo,yr) = obs(mo,yr)
            enddo
        enddo
        if ( lwrite ) print *,'debiasmean: taking log10(obs)'
        call takelog(obs(1,yr1),npermax,nperyear,yr1,yr2)
        do iens=nens1,nens2
            do mo=1,nperyear
                yrzero(mo,iens) = -999
                do yr=yr1,yr2
                    if ( fcst(mo,yr,iens) <= 0 ) then
                        yrzero(mo,iens) = yr
                        exit
                    endif
                enddo
            enddo
            if ( lwrite ) print *,'debiasvar: taking log10(fcst) ',iens
            call takelog(fcst(1,yr1,iens),npermax,nperyear,yr1,yr2)
        enddo
    endif
    do mo=1,nperyear
        s1 = 0
        n1 = 0
        do yr=yr1,yr2
            if ( obs(mo,yr) < 1e33 ) then
                n1 = n1 + 1
                s1 = s1 + obs(mo,yr)
            endif
        enddo
        s2 = 0
        n2 = 0
        do iens=nens1,nens2
            do yr=yr1,yr2
                if ( fcst(mo,yr,iens) < 1e33 ) then
                    n2 = n2 + 1
                    s2 = s2 + fcst(mo,yr,iens)
                endif
            enddo
        enddo
        if ( n1 <= 2 .or. n2 <= 2 ) then
            if ( lwrite ) print *,'debiasmean: not enough valid data for ',mo
            do iens=nens1,nens2
                do yr=yr1,yr2
                    newfcst(mo,yr,iens) = 3e33
                enddo
            enddo
            goto 800
        endif
        s1m = s1/n1
        s2m = s2/n2
        if ( lwrite ) then
            print *,'debiasvar: mean obs = ',mo,s1m
            print *,'           mean fcst= ',mo,s2m
        endif
        t1 = 0
        do yr=yr1,yr2
            if ( obs(mo,yr) < 1e33 ) then
                t1 = t1 + (obs(mo,yr)-s1m)**2
            endif
        enddo
        t2 = 0
        do iens=nens1,nens2
            do yr=yr1,yr2
                if ( fcst(mo,yr,iens) < 1e33 .and. s2m < 1e33 ) then
                    t2 = t2 + (fcst(mo,yr,iens)-s2m)**2
                endif
            enddo
        enddo
        if ( lwrite ) then
            print *,'debiasvar: s.d. obs = ',mo,sqrt(t1/(n1-1))
            print *,'           s.d. fcst= ',mo,sqrt(t2/(n2-1))
        endif
    
!       jackknife
    
        do yr=yr1,yr2
            if ( obs(mo,yr) < 1e33 ) then
                ss1 = (s1 - obs(mo,yr))/(n1-1)
                tt1 = (t1 - (obs(mo,yr)-s1m)**2 &
                    - 2*obs(mo,yr)/n1*(obs(mo,yr) -s1m)/n1 &
                    + (n1-1)/real(n1**2)*obs(mo,yr)**2)/(n1-2)
            else
                ss1 = s1/n1
                tt1 = tt1/(n1-1)
            endif
            ss2 = s2
            nn2 = n2
            tt2 = t2
            do iens=nens1,nens2
                if ( fcst(mo,yr,iens) < 1e33 ) then
                    ss2 = ss2 - fcst(mo,yr,iens)
                    tt2 = tt2 - (fcst(mo,yr,iens)-s2m)**2 &
                        - 2*fcst(mo,yr,iens)/n2*(fcst(mo,yr,iens) - s2m)/n2 &
                        + (n2-1)/real(n2**2)*fcst(mo,yr,iens)**2
                    nn2 = nn2 - 1
                endif
            enddo
            if ( nn2 <= 0 ) then
                ss2 = ss1  ! no bias correction
            else
                ss2 = ss2/nn2
            endif
            if ( nn2 <= 1 ) then
                tt2 = tt1  ! no bias correction
            else
                tt2 = tt2/(nn2-1)
            endif
            if ( lwrite ) then
                print *,'bias mean = ',mo,yr,ss1,ss2,ss1-ss2
                print *,'bias var  = ',mo,yr,tt1,tt2,tt1/tt2
            endif
            if ( tt2 < 0.1*tt1 ) then
                tt2 = 0.1*tt1
                if ( init <= 1 ) then
                    init = 2
                    print '(a)','# refusing to inflate variance by more than a factor 10'
                endif
            endif
            do iens=nens1,nens2
                if ( fcst(mo,yr,iens) < 1e33 ) then
                    newfcst(mo,yr,iens) = (fcst(mo,yr,iens) - ss2)*sqrt(tt1/tt2) + ss1
                    if ( newfcst(mo,yr,iens) > 1e20 ) then
                        if ( nunlikely < 20 ) then
                            nunlikely = nunlikely + 1
                            print *,'debiasvar: unlikely value for debiassed forecast '
                            print *,'           s1,s2   = ',s1,s2
                            print *,'           ss1,ss2 = ',ss1,ss2
                            print *,'           t1,t2   = ',t1,t2
                            print *,'           tt1,tt2 = ',tt1,tt2
                            print *,'           fcst(',mo,yr,iens,') = ',fcst(mo,yr,iens)
                            print *,'           newfcst(',mo,yr,iens,') = ',newfcst(mo,yr,iens)
                        endif
                    endif
                else
                    newfcst(mo,yr,iens) = 3e33
                endif
            enddo
        enddo
    800 continue
    enddo

!   copy back to original arrays, setting the rest to undefined
!   for security reasons

    do iens=nens1,nens2
        do yr=yrbeg,yr1-1
            do mo=1,nperyear
                fcst(mo,yr,iens) = 3e33
            enddo
        enddo
        do yr=yr1,yr2
            do mo=1,nperyear
                fcst(mo,yr,iens) = newfcst(mo,yr,iens)
            enddo
        enddo
        do yr=yr2+1,yrend
            do mo=1,nperyear
                fcst(mo,yr,iens) = 3e33
            enddo
        enddo
    enddo
    if ( isvarprcp(var) ) then
        if ( lwrite ) print *,'debiasmean: copying back saved obs'
        do yr=yr1,yr2
            do mo=1,nperyear
                obs(mo,yr) = savobs(mo,yr)
            enddo
        enddo
        deallocate(savobs)
        do iens=nens1,nens2
            if ( lwrite ) print *,'debiasmean: taking 10^fcst ',iens
            call takeexp(fcst(1,yr1,iens),npermax,nperyear,yr1,yr2)
            do mo=1,nperyear
                if ( yrzero(mo,iens) /= -999 ) then
                    zero = fcst(mo,yrzero(mo,iens),iens)
                    do yr=yr1,yr2
                        if ( fcst(mo,yr,iens) == zero ) fcst(mo,yr,iens) = 0
                    enddo
                endif
            enddo
        enddo
        deallocate(yrzero)
    endif
    deallocate(newfcst)

end subroutine debiasvar

subroutine debiasall(obs,fcst,npermax,nperyear,yrbeg,yrend &
    ,yr1,yr2,nens1,nens2,var,lwrite)

!   adjust the array fcst so that for each year the PDF of all
!   other years is equal to the PDF of all other observations.

    implicit none
    integer,intent(in) :: npermax,nperyear,yrbeg,yrend,yr1,yr2,nens1,nens2
    real,intent(inout) :: obs(npermax,yrbeg:yrend),fcst(npermax,yrbeg:yrend,0:nens2)
    character,intent(in) :: var*(*)
    logical,intent(in) :: lwrite
    integer :: yr,mo,iens,i,j,nyr,nobs1,nfcst1
    real :: perc
    real,allocatable :: obs1(:),fcst1(:),newfcst(:,:)

!   no jack-knife possible with only one year
    if ( yrend == yrbeg ) return

!   allocate working arrays
    nyr = yrend-yrbeg
    allocate(obs1(nyr))
    allocate(fcst1(nyr*(nens2-nens1+1)))
    allocate(newfcst(yrbeg:yrend,nens1:nens2))

!   loop over months - all independent
    do mo=1,nperyear
    
!       jack-knife
        do yr=yrbeg,yrend
        
!           copy to working arrays
            if ( yr > yrbeg ) obs1(1:yr-yrbeg) = obs(mo,yrbeg:yr-1)
            if ( yr < yrend ) obs1(yr-yrbeg+1:nyr) = obs(mo,yr+1:yrend)
            call nrsort(nyr,obs1)
            do iens=nens1,nens2
                do i=yrbeg,yr-1
                    fcst1(1+i-yrbeg+nyr*(iens-nens1)) = fcst(mo,i,iens)
                enddo
                do i=yr+1,yrend
                    fcst1(i-yrbeg+nyr*(iens-nens1)) = fcst(mo,i,iens)
                enddo
            enddo
            call nrsort(nyr*(nens2-nens1+1),fcst1)
        
!           count the number of valid obs1,fcst1 elements
!           (there must be an f90 intrinsic for this)
            do nobs1=nyr,1,-1
                if ( obs1(nobs1) < 1e33 ) exit
            enddo
            if ( nobs1 == 0 ) then
                newfcst = 3e33
                exit
            endif
            do nfcst1=1,nyr*(nens2-nens1+1)
                if ( fcst1(nfcst1) > 1e33 ) exit
            enddo
            nfcst1 = nfcst1 - 1
            if ( nfcst1 == 0 ) then
                newfcst = 3e33
                exit        ! if there is only one year it's not enough
            endif
        
!           compute the percentile of the current year forecast in
!           the PDF of all other years
            do iens = nens1,nens2
                call val2frac(fcst1,nfcst1,fcst(mo,yr,iens),perc)
                call frac2val(perc,obs1,nobs1,newfcst(yr,iens))
            enddo
        enddo               ! yr
    
!       copy back
        fcst(mo,:,:) = newfcst(:,:)
    enddo                   ! mo

!   deallocate arrays
    deallocate(obs1)
    deallocate(fcst1)
    deallocate(newfcst)

end subroutine debiasall

logical function isvarprcp(var)
    implicit none
    character,intent(in) :: var*(*)
    if ( var == 'prec' .or. var == 'tp' .or. var == 'prcp' .or. &
         var == 'prate' .or. var == 'rr' .or. var == 'pr' ) then
        isvarprcp = .true. 
    else
        isvarprcp = .false. 
    endif
end function isvarprcp
