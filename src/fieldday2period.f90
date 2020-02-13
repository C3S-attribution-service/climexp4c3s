subroutine fieldallday2period( &
    oldfield,nperyear,lvalid, &
    newfield,nperyearnew, &
    nx,ny,yrbeg,yrend,oper,lgt,cut,minfac,add_option, &
    itype,var,units,lwrite)

    implicit none
    integer :: nperyear,nx,ny,nperyearnew,yrbeg,yrend,itype,add_option
    real :: oldfield(nx,ny,nperyear,yrbeg:yrend), &
        newfield(nx,ny,abs(nperyearnew),yrbeg:yrend), &
        cut(nx,ny,nperyear),minfac
    character :: oper*3,lgt*1,var*(*),units*(*)
    logical :: lvalid(nx,ny,nperyear,yrbeg:yrend),lwrite
    integer :: j,jj,j1,j2,jold,mo,n,dpm(12),dpm365(12)
    data dpm   / 31,29,31,30,31,30,31,31,30,31,30,31/
    data dpm365/ 31,28,31,30,31,30,31,31,30,31,30,31/

    if ( lwrite ) then
        print *,'fieldallday2period: nperyear,nperyearnew = ',nperyear,nperyearnew
        print *,'                    nx,ny,yrbeg,yrend = ',nx,ny,yrbeg,yrend
        print *,'                    oper,lgt,cut,itype= ',oper,lgt,cut(1,1,1),itype
        print *,'                    var,units         = ',var,units
        if ( yrend > yrbeg ) then
            print *,'oldfield = ',oldfield(1,1,1,yrbeg),oldfield(1,1 &
            ,1,yrbeg+1),'...,',oldfield(1,1,1,yrend)
        else
            print *,'oldfield = ',oldfield(1,1,1,yrbeg)
        endif
    endif

    if ( nperyearnew == -1 ) then
        if ( nperyear < 365 ) then
            j1 = 1 + nperyear/2
        else if ( nperyear == 365 ) then
            j1 = 182        ! 1 July
        else if ( nperyear == 366 ) then
            j1 = 183        ! 1 July including Feb 29
        else
            write(0,*) 'allday2period: error: cannot handle '// &
                'shifted annual data for sub-daily data yet'
            call exit(-1)
        end if
        j2 = j1 + nperyear - 1
        call fieldday2period( &
            oldfield,nperyear,lvalid, &
            newfield,abs(nperyearnew), &
            nx,ny,yrbeg,yrend,j1,j2,1,oper,lgt,cut, &
            minfac,add_option,itype,lwrite)
    else if ( nperyearnew == 1 ) then
        call fieldday2period( &
            oldfield,nperyear,lvalid, &
            newfield,nperyearnew, &
            nx,ny,yrbeg,yrend,1,nperyear,1,oper,lgt,cut, &
            minfac,add_option,itype,lwrite)
    elseif ( nperyearnew == 4 ) then
        if ( nperyear /= 360 ) then
            j = -31
        else
            j = -30
        endif
        jold = nint(j/(366./nperyear))
        do mo=0,11
            if ( nperyear == 360 ) then
                j = j + 30
                jj = j
            elseif ( nperyear == 365 ) then
                if ( mo > 0 ) then
                    j = j + dpm365(mo)
                else
                    j = j + dpm365(mo+12)
                endif
                jj = nint(j/(365./nperyear))
            else
                if ( mo > 0 ) then
                    j = j + dpm(mo)
                else
                    j = j + dpm(mo+12)
                endif
                jj = nint(j/(366./nperyear))
            endif
            if ( mo == 2 .or. mo == 5 .or. mo == 8 .or. mo == 11 ) &
            then
                call keepalive2('Season',1+mo/4,4, .true. )
                call fieldday2period( &
                    oldfield,nperyear,lvalid, &
                    newfield,nperyearnew, &
                    nx,ny,yrbeg,yrend,jold+1,jj,(mo+1)/3,oper,lgt &
                    ,cut,minfac,add_option,itype,lwrite)
                jold = jj
            endif
        enddo
    elseif ( nperyearnew == 12 ) then
        j = 0
        jold = 0
        do mo=1,12
            if ( nperyear == 360 ) then
                j = j + 30
                jj = j
            elseif ( nperyear == 365 ) then
                j = j + dpm365(mo)
                jj = nint(j/(365./nperyear))
            else
                j = j + dpm(mo)
                jj = nint(j/(366./nperyear))
            endif
            call keepalive2('Month',mo,12, .true. )
            call fieldday2period( &
                oldfield,nperyear,lvalid, &
                newfield,nperyearnew, &
                nx,ny,yrbeg,yrend,jold+1,jj,mo,oper,lgt,cut, &
                minfac,add_option,itype,lwrite)
            jold = jj
        enddo
    elseif ( nperyearnew == 36 ) then
        j = 0
        jold = 0
        do mo=1,12
            call keepalive2('Month',mo,12, .true. )
            j = j + 10
            if ( nperyear == 360 .or. nperyear == 366 ) then
                jj = j
            else
                jj = nint(j/(365./nperyear))
            endif
            call fieldday2period( &
                oldfield,nperyear,lvalid, &
                newfield,nperyearnew, &
                nx,ny,yrbeg,yrend,jold+1,jj,3*mo-2,oper,lgt,cut, &
                minfac,add_option,itype,lwrite)
            jold = jj
            j = j + 10
            if ( nperyear == 360 .or. nperyear == 366 ) then
                jj = j
            else
                jj = nint(j/(365./nperyear))
            endif
            call fieldday2period( &
                oldfield,nperyear,lvalid, &
                newfield,nperyearnew, &
                nx,ny,yrbeg,yrend,jold+1,jj,3*mo-1,oper,lgt,cut, &
                minfac,add_option,itype,lwrite)
            jold = jj
            if ( nperyear == 360 ) then
                j = j + 10
                jj = j
            elseif ( nperyear == 365 ) then
                j = j + dpm365(mo) - 20
                jj = j
            elseif ( nperyear == 366 ) then
                j = j + dpm(mo) - 20
                jj = j
            else
                j = j + dpm(mo) - 20
                jj = nint(j/(365./nperyear))
            endif
            call fieldday2period( &
                oldfield,nperyear,lvalid, &
                newfield,nperyearnew, &
                nx,ny,yrbeg,yrend,jold+1,jj,3*mo,oper,lgt,cut, &
                minfac,add_option,itype,lwrite)
            jold = jj
        enddo
    else
        n = nint(real(nperyear)/real(nperyearnew))
        do j=1,nperyearnew
            call keepalive2('Period',j,nperyearnew, .true. )
            call fieldday2period( &
                oldfield,nperyear,lvalid, &
                newfield,nperyearnew, &
                nx,ny,yrbeg,yrend,n*(j-1)+1,n*j,j,oper,lgt,cut, &
                minfac,add_option,itype,lwrite)
        enddo
    endif

    call adjustunits(oper,nperyear,nperyearnew,units,lwrite)

end subroutine fieldallday2period

subroutine fieldday2period( &
    oldfield,nperyear,lvalid, &
    newfield,nperyearnew, &
    nx,ny,yrbeg,yrend,j1,j2,jnew,oper,lgt,cut, &
    minfac,add_option,itype,lwrite)

!       operates on oldfield (j1:j2,) to make newfield(jnew,)
!       oper = mean|sd|var|min|max|num|sum
!       lgt  = ' '|<|>
!       cut  = cut-off
!       itype=0 scalar data
!       itype=1 direction data

!       15-jun-2004 also compute mean when there are missing data
!       2-nov-2005 adapted for fields

    implicit none
    integer :: nx,ny,nperyear,nperyearnew,yrbeg,yrend,j1,j2,jnew,itype,add_option
    real :: oldfield(nx,ny,nperyear,yrbeg:yrend), &
        newfield(nx,ny,nperyearnew,yrbeg:yrend), &
        cut(nx,ny,nperyear),minfac
    character oper*3,lgt*1
    logical lvalid(nx,ny,nperyear,yrbeg:yrend),lwrite
    integer :: yr,mo,i,j,n,jx,jy,ntot,offset,lfirst,ncount
    real :: s,s1,s2,sx,sy,mean,st

    if ( lwrite ) then
        print *,'fieldday2period: parameters are '
        print *,'nperyear,nperyearnew= ',nperyear,nperyearnew
        print *,'yrbeg,yrend         = ',yrbeg,yrend
        print *,'j1,j2,jnew          = ',j1,j2,jnew
        print *,'oper,lgt            = ',oper,lgt
        print *,'itype               = ',itype
        if ( lgt /= ' ' ) then
            print *,'cut(1,1)            = '
            print *,(cut(1,1,j),j=1,nperyear)
        endif
        print *,'oldfield(1,1)        = '
        do yr=yrbeg,yrend,max(1,yrend-yrbeg)
            do mo=j1,j2
                j = mo
                call normon(j,yr,i,nperyear)
                if ( yrbeg == 0 .and. yrend == 0 ) i = yr
                if ( i < yrbeg .or. i > yrend ) cycle
                if ( oldfield(1,1,j,i) < 1e33 ) then
                    print *,yr,j,oldfield(1,1,j,i)
                endif
            enddo
        enddo
    endif
    if ( itype == 1 .and. (oper == 'min' .or. oper == 'max' .or. &
         oper == 'nti' .or. oper == 'xti' .or. oper == 'fti' .or. &
         oper == 'lti' .or.  oper == 'sd' .or. oper == 'con') ) then
        write(0,*) 'fieldday2period: cannot take min,max,sd of vector data'
        write(*,*) 'fieldday2period: cannot take min,max,sd of vector data'
        call exit(-1)
    endif

    ncount = 0
    do yr=yrbeg,yrend
        ncount = ncount + 1
        call keepalive2('Year',ncount,yrend-yrbeg+1, .true. )
        do jy=1,ny
            do jx=1,nx
                call keepalive1('Grid point',jx+nx*(jy-1),nx*ny)
                if ( (oper == 'sd' .or. oper == 'var' ) .and. itype == 0 ) then
!                   precompute the mean for numerical stability
                    n = 0
                    s = 0
                    do mo=j1,j2
                        j = mo
                        call normon(j,yr,i,nperyear)
                        if ( yrbeg == 0 .and. yrend == 0 ) i = yr
                        if ( i < yrbeg ) cycle
!                           skip Feb 29
                        if ( nperyear == 366 .and. j == 60 .and. &
                        oldfield(jx,jy,j,i) > 1e33 ) goto 40
!                           but generate an "invalid" for any other
!                           invalid data [this may be relaxed later...]
                        if ( oldfield(jx,jy,j,i) > 1e33 ) goto 40
                        n = n + 1
                        s = s + oldfield(jx,jy,j,i)
                     40 continue
                    enddo
                    if ( n > 1 ) then
                        mean = s/n
                    else
                        mean = 0
                    endif
                endif
                if ( oper == 'mea' .or. oper == 'sd ' .or. &
                     oper == 'var' .or. &
                     oper == 'sum' .or. oper == 'bel' .or. &
                     oper == 'abo' .or. oper == 'con' ) then
                    if ( itype == 0 ) then
                        s = 0
                    elseif ( itype == 1 ) then
                        sx = 0
                        sy = 0
                    else
                        write(0,*) 'fieldday2period: error: unknown itype ',itype
                        write(*,*) 'fieldday2period: error: unknown itype ',itype
                        call exit(-1)
                    endif
                    if ( oper == 'sd ' .or. oper == 'var' ) s2 = 0
                elseif ( oper == 'min' .or. oper == 'nti' ) then
                    if ( itype == 0 ) then
                        s = 3e33
                    elseif ( itype == 1 ) then
                        sx = 3e33
                        sy = 3e33
                    endif
                elseif ( oper == 'max' .or. oper == 'xti' ) then
                    if ( itype == 0 ) then
                        s = -3e33
                    elseif ( itype == 1 ) then
                        sx = -3e33
                        sy = -3e33
                    endif
                endif
                st = 3e33   ! time of max/min
                s1 = 0      ! for consecutive days
                n = 0
                ntot = 0
                lfirst = 9999
                offset = 0
                do mo=j1,j2
                    j = mo
                    call normon(j,yr,i,nperyear)
                    if ( yrbeg == 0 .and. yrend == 0 ) i = yr
                    if ( i < yrbeg .or. i > yrend ) cycle
!                   skip Feb 29
                    if ( nperyear == 366 .and. j == 60 .and. oldfield(jx,jy,j,i) > 1e33 ) then
                        offset = -1
                        cycle
                    end if
!                   but generate an "invalid" for any other invalid
!                   data [this may be relaxed later...]
                    if ( oldfield(jx,jy,j,i) > 1e33 ) then
!                       it has been relaxed
!**                     newfield(jnew,i) = 3e33
                        cycle
                    endif
                    if ( lfirst == 9999 ) lfirst = mo-j1
!                   this test should be exactly the same as in daily2longerfield
                    if ( ( oper == 'mea' .or. oper == 'sum' ) .and. &
                            lgt == ' ' .and. add_option > 0 ) then
                        if ( lvalid(jx,jy,j,i) ) ntot = ntot + 1
                    else
                        ntot = ntot + 1
                    end if
                    if ( lgt == ' ' .or. &
                         lgt == '<' .and. oldfield(jx,jy,j,i) < cut(jx,jy,j) .or. &
                         lgt == '>' .and. oldfield(jx,jy,j,i) > cut(jx,jy,j) ) then
                        n = n + 1
                        if ( oper == 'fti' ) then
                            if ( st == 3e33 ) st = mo - 0.5 + offset
                        else if ( oper == 'lti' ) then
                            st = mo - 0.5 + offset
                        else if ( oper == 'mea' .or. oper == 'sum' &
                                  .or. &
                                  oper == 'sd ' .or. oper == 'var' ) &
                                then
                            if ( itype == 0 ) then
                                s = s + oldfield(jx,jy,j,i)
                            elseif ( itype == 1 ) then
                                sx = sx + cos(atan(1.)*oldfield(jx,jy,j,i)/45)
                                sy = sy + sin(atan(1.)*oldfield(jx,jy,j,i)/45)
                            endif
                            if ( oper == 'sd ' .or. oper == 'var' ) &
                            then
                                s2 = s2 &
                                + (oldfield(jx,jy,j,i)-mean)**2
                            endif
                        elseif ( oper == 'abo' ) then
                            s = s + oldfield(jx,jy,j,i)-cut(jx,jy,j)
                        elseif ( oper == 'bel' ) then
                            s = s + cut(jx,jy,j)-oldfield(jx,jy,j,i)
                        elseif ( oper == 'min' .or. oper == 'nti' ) then
                            if ( oldfield(jx,jy,j,i) < s ) then
                                s = oldfield(jx,jy,j,i)
                                st = mo - 0.5 + offset
                            end if
                        elseif ( oper == 'max' .or. oper == 'xti') then
                            if ( oldfield(jx,jy,j,i) > s ) then
                                s = oldfield(jx,jy,j,i)
                                st = mo - 0.5 + offset
                            end if
                        else if ( oper == 'con' ) then
                            s1 = s1 + 1
                            s = max(s,s1)
                        endif
                        if ( lwrite .and. &
                        jx == (nx+1)/2 .and. jy == (ny+1)/2 ) &
                        print *,'partial sum: ',j,n,s
                    else
                        if ( oper == 'con' ) then
                            s1 = 0
                        end if
                        if ( lwrite .and. jx == (nx+1)/2 .and. jy == (ny+1)/2 ) &
                            print *,j,oldfield(jx,jy,j,i),lgt,cut(jx,jy,j)
                    endif
                enddo       ! j
                if ( .false. .and. lwrite .and. lfirst < 9999 ) then
                    write(*,*) yr,'lfirst,minfac*(j2-j1+1) = ',lfirst,minfac*(j2-j1+1)
                    write(*,*) yr,'ntot,minfac*nperyear/nperyearnew = ', &
                                   ntot,minfac*nperyear/nperyearnew
                end if
                if ( lfirst > minfac*(j2-j1+1) .or. &
                        ntot < minfac*nperyear/nperyearnew ) then
                    newfield(jx,jy,jnew,yr) = 3e33
                elseif ( n == 0 .and. &
                    (oper == 'mea' .or. oper == 'min' .or. &
                     oper == 'max' .or. oper == 'xti' .or. &
                     oper == 'nti' .or. oper == 'fti' .or. &
                     oper == 'lti' ) ) then
                    newfield(jx,jy,jnew,yr) = 3e33
                elseif ( oper == 'mea' ) then
                    if ( itype == 0 ) then
                        if ( lwrite .and. &
                        jx == (nx+1)/2 .and. jy == (ny+1)/2 ) &
                            print *,'s,n,s/n = ',s,n,s/n
                        newfield(jx,jy,jnew,yr) = s/n
                    elseif ( itype == 1 ) then
                        if ( lwrite .and. &
                        jx == (nx+1)/2 .and. jy == (ny+1)/2 ) then
                            print *,'sx,n,sx/n = ',sx,n,sx/n
                            print *,'sy,n,sy/n = ',sy,n,sy/n
                            print *,'atan2(sy,sx)',atan2(sy,sx)
                        endif
                        newfield(jx,jy,jnew,yr) = 45*atan2(sy,sx)/atan(1.)
                        if ( newfield(jx,jy,jnew,yr) < 0 ) &
                            newfield(jx,jy,jnew,yr) = newfield(jx,jy,jnew,yr) + 360
                    endif
                elseif ( oper == 'sd ' .or. oper == 'var' ) then
                    newfield(jx,jy,jnew,yr) = s2/n - (s/n-mean)**2
                    if ( oper == 'sd' ) then
                        if ( newfield(jx,jy,jnew,yr) >= 0 ) then
                            newfield(jx,jy,jnew,yr) = &
                            sqrt(newfield(jx,jy,jnew,yr))
                        else
                            newfield(jx,jy,jnew,yr) = 3e33
                        end if
                    end if
                elseif ( oper == 'sum' .or. oper == 'bel' .or. &
                         oper == 'abo' .or. oper == 'min' .or. &
                         oper == 'max' ) then
                    if ( itype == 0 ) then
                        newfield(jx,jy,jnew,yr) = s
                    elseif ( itype == 1 ) then
                        newfield(jx,jy,jnew,yr) = 45*atan2(sy,sx)/atan(1.)
                        if ( newfield(jx,jy,jnew,yr) < 0 ) &
                            newfield(jx,jy,jnew,yr) = newfield(jx,jy,jnew,yr) + 360
                    endif
                elseif ( oper == 'num' ) then
                    newfield(jx,jy,jnew,yr) = n
                elseif ( oper == 'nti' .or. oper == 'xti' .or. &
                    oper == 'fti' .or. oper == 'lti' ) then
                    newfield(jx,jy,jnew,yr) = st
                elseif ( oper == 'con' ) then
                    newfield(jx,jy,jnew,yr) = s
                else
                    write(0,*) 'fieldday2period: error: unknown operation ',oper
                    call exit(-1)
                endif
                if ( lwrite .and. jx == (nx+1)/2 .and. jy == (ny+1)/2 ) then
                    print *,'fieldday2period: computed newfield'
                    print *,yr,jnew,newfield(jx,jy,jnew,yr),s,st,itype
                endif
            100 continue
            enddo           ! jx
        enddo               ! jy
    enddo                   ! yr
end subroutine fieldday2period

subroutine adjustnames(oper,nperyear,nperyearnew,lgt,pcut,punits,lvars,cell_methods,lsum)

!   adjust names to reflect the operation just performed

    implicit none
    integer :: nperyear,nperyearnew,lsum
    real :: pcut
    character oper*3,lgt*1,punits*(*),lvars*(*),cell_methods*(*)
    integer :: i,j
    character prefix*255,postfix*255
            
    call nperyear2string(nperyearnew,prefix)
    if ( oper == 'mea' ) then
        prefix = trim(prefix)//' mean'
    elseif ( oper == 'sd ' ) then
        prefix = trim(prefix)//' standard deviation'
    elseif ( oper == 'var' ) then
        prefix = trim(prefix)//' variance'
    elseif ( oper == 'min' ) then
        prefix = trim(prefix)//' minimum'
    elseif ( oper == 'max' ) then
        prefix = trim(prefix)//' maximum'
    elseif ( oper == 'nti' ) then
        prefix = 'time of '//trim(prefix)//' minimum'
    elseif ( oper == 'xti' ) then
        prefix = 'time of '//trim(prefix)//' maximum'
    elseif ( oper == 'fti' ) then
        prefix = 'time of '//trim(prefix)//' first'
    elseif ( oper == 'lti' ) then
        prefix = 'time of '//trim(prefix)//' last'
    elseif ( oper == 'num' ) then
        prefix = trim(prefix)//' number'
    elseif ( oper == 'con' ) then
        prefix = trim(prefix)//' max number of consecutive number'
    elseif ( oper == 'bel' .or. oper == 'abo' .or. &
        oper == 'sum' ) then
        prefix = trim(prefix)//' sum '
    else
        write(0,*) 'adjustnames: unknown operation: ',oper
        write(*,*) 'adjustnames: unknown operation: ',oper
        call exit(-1)
    endif
    prefix = trim(prefix)//' of'
    i = len_trim(prefix)
    if ( lsum > 1 ) then
        write(prefix(i+2:),'(i3,a)') lsum,'-'
    end if
    i = len_trim(prefix)
    if ( prefix(i:i) == '-' ) then
        j = 1
    else
        j = 2
    end if
    call nperyear2string(nperyear,prefix(i+j:))

    if ( lgt == ' ' ) then
        postfix = ' '
    else
        if ( lgt == '<' ) then
            postfix = 'below'
        elseif ( lgt == '>' ) then
            postfix = 'above'
        endif
        write(postfix(7:),'(g20.6,2a)') pcut,' ',punits
    endif
    lvars = trim(prefix)//' '//trim(lvars)//' '//trim(postfix)
    call deletedoublespaces(lvars)
    if ( cell_methods == ' ' ) then
        cell_methods = trim(prefix)//' data '//trim(postfix)
    else
        cell_methods = trim(prefix)//' '//trim(postfix)//' '// &
        trim(cell_methods)
    endif
    call deletedoublespaces(cell_methods)
end subroutine 

subroutine deletedoublespaces(string)
    implicit none
    character*(*) string
    integer :: i,j,k,state
    if ( string == ' ' ) return
    if ( string(1:1) == ' ' ) then
        state = 0
    else
        state = 1
    endif
    i = 1
    do
        if ( state == 0 ) then
!       search for more whitespace
            j = i+1
            do while ( string(j:j) == ' ' )
                j = j + 1
                if ( j > len(string) ) return
            end do
            if ( j > i+1 ) then
            !                   multiple whitespace => shorten
                string(i+1:) = string(j:)
            endif
            i = i+1
            state = 1
        else
            j = i+1
            do while ( string(j:j) /= ' ' )
                j = j + 1
                if ( j > len(string) ) return
            end do
            i = j
            state = 0
        endif
    end do
end subroutine deletedoublespaces
