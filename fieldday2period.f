        subroutine fieldallday2period(
     +       oldfield,nperyear,lvalid,
     +       newfield,nperyearnew,
     +       nx,ny,yrbeg,yrend,oper,lgt,cut,minfac,itype,var,units,
     +       lwrite)
        implicit none
        integer nperyear,nx,ny,nperyearnew,yrbeg,yrend,itype
        real oldfield(nx,ny,nperyear,yrbeg:yrend),
     +       newfield(nx,ny,abs(nperyearnew),yrbeg:yrend),
     +       cut(nx,ny,nperyear),minfac
        character oper*3,lgt*1,var*(*),units*(*)
        logical lvalid(nx,ny,nperyear,yrbeg:yrend),lwrite
        integer j,jj,j1,jold,mo,n,dpm(12),dpm365(12)
        data dpm   / 31,29,31,30,31,30,31,31,30,31,30,31/
        data dpm365/ 31,28,31,30,31,30,31,31,30,31,30,31/
*
        if ( lwrite ) then
            print *,'fieldallday2period: nperyear,nperyearnew = '
     +           ,nperyear,nperyearnew
            print *,'                    nx,ny,yrbeg,yrend = ',nx,ny
     +           ,yrbeg,yrend
            print *,'                    oper,lgt,cut,itype= ',oper,lgt
     +           ,cut(1,1,1),itype
            print *,'                    var,units         = ',var,units
            if ( yrend.gt.yrbeg ) then
                print *,'oldfield = ',oldfield(1,1,1,yrbeg),oldfield(1,1
     +               ,1,yrbeg+1),'...,',oldfield(1,1,1,yrend)
            else
                print *,'oldfield = ',oldfield(1,1,1,yrbeg)
            endif
        endif

        if ( nperyearnew.eq.-1 ) then
            if ( nperyear.lt.365 ) then
                j1 = 1 + nperyear/2
            else if ( nperyear.eq.365 ) then
                j1 = 182        ! 1 July
            else if ( nperyear.eq.366 ) then
                j1 = 183        ! 1 July including Feb 29
            else
                write(0,*) 'allday2period: error: cannot handle '//
     +               'shifted annual data for sub-daily data yet'
                call abort
            end if
            call fieldday2period(
     +           oldfield,nperyear,lvalid,
     +           newfield,abs(nperyearnew),
     +           nx,ny,yrbeg,yrend,1,nperyear,1,oper,lgt,cut,
     +           minfac,itype,lwrite)
        else if ( nperyearnew.eq.1 ) then
            call fieldday2period(
     +           oldfield,nperyear,lvalid,
     +           newfield,nperyearnew,
     +           nx,ny,yrbeg,yrend,1,nperyear,1,oper,lgt,cut,
     +           minfac,itype,lwrite)
        elseif ( nperyearnew.eq.4 ) then
            if ( nperyear.ne.360 ) then
                j = -31
            else
                j = -30
            endif
            jold = nint(j/(366./nperyear))
            do mo=0,11
                if ( nperyear.eq.360 ) then
                    j = j + 30
                    jj = j
                elseif ( nperyear.eq.365 ) then
                    if ( mo.gt.0 ) then
                        j = j + dpm365(mo)
                    else
                        j = j + dpm365(mo+12)
                    endif
                    jj = nint(j/(365./nperyear))
                else
                    if ( mo.gt.0 ) then
                        j = j + dpm(mo)
                    else
                        j = j + dpm(mo+12)
                    endif
                    jj = nint(j/(366./nperyear))
                endif
                if ( mo.eq.2 .or. mo.eq.5 .or. mo.eq.8 .or. mo.eq.11 )
     +               then
                    call fieldday2period(
     +                   oldfield,nperyear,lvalid,
     +                   newfield,nperyearnew,
     +                   nx,ny,yrbeg,yrend,jold+1,jj,(mo+1)/3,oper,lgt
     +                   ,cut,minfac,itype,lwrite)
                    jold = jj
                endif
            enddo
        elseif ( nperyearnew.eq.12 ) then
            j = 0
            jold = 0
            do mo=1,12
                if ( nperyear.eq.360 ) then
                    j = j + 30
                    jj = j
                elseif ( nperyear.eq.365 ) then
                    j = j + dpm365(mo)
                    jj = nint(j/(365./nperyear))
                else
                    j = j + dpm(mo)
                    jj = nint(j/(366./nperyear))
                endif
                call fieldday2period(
     +               oldfield,nperyear,lvalid,
     +               newfield,nperyearnew,
     +               nx,ny,yrbeg,yrend,jold+1,jj,mo,oper,lgt,cut,
     +               minfac,itype,lwrite)
                jold = jj
            enddo
        elseif ( nperyearnew.eq.36 ) then
            j = 0
            jold = 0
            do mo=1,12
                j = j + 10
                if ( nperyear.eq.360 .or. nperyear.eq.366 ) then
                    jj = j
                else
                    jj = nint(j/(365./nperyear))
                endif
                call fieldday2period(
     +               oldfield,nperyear,lvalid,
     +               newfield,nperyearnew,
     +               nx,ny,yrbeg,yrend,jold+1,jj,3*mo-2,oper,lgt,cut,
     +               minfac,itype,lwrite)
                jold = jj
                j = j + 10
                if ( nperyear.eq.360 .or. nperyear.eq.366 ) then
                    jj = j
                else
                    jj = nint(j/(365./nperyear))
                endif
                call fieldday2period(
     +               oldfield,nperyear,lvalid,
     +               newfield,nperyearnew,
     +               nx,ny,yrbeg,yrend,jold+1,jj,3*mo-1,oper,lgt,cut,
     +               minfac,itype,lwrite)
                jold = jj
                if ( nperyear.eq.360 ) then
                    j = j + 10
                    jj = j
                elseif ( nperyear.eq.365 ) then
                    j = j + dpm365(mo) - 20
                    jj = j
                elseif ( nperyear.eq.366 ) then
                    j = j + dpm(mo) - 20
                    jj = j
                else
                    j = j + dpm(mo) - 20
                    jj = nint(j/(365./nperyear))
                endif
                call fieldday2period(
     +               oldfield,nperyear,lvalid,
     +               newfield,nperyearnew,
     +               nx,ny,yrbeg,yrend,jold+1,jj,3*mo,oper,lgt,cut,
     +               minfac,itype,lwrite)
                jold = jj
            enddo
        else
            n = nint(real(nperyear)/real(nperyearnew))
            do j=1,nperyearnew
                call fieldday2period(
     +               oldfield,nperyear,lvalid,
     +               newfield,nperyearnew,
     +               nx,ny,yrbeg,yrend,n*(j-1)+1,n*j,j,oper,lgt,cut,
     +               minfac,itype,lwrite)
            enddo
        endif
*
        call adjustunits(oper,nperyear,nperyearnew,units,lwrite)
*
        end

        subroutine fieldday2period(
     +       oldfield,nperyear,lvalid,
     +       newfield,nperyearnew,
     +       nx,ny,yrbeg,yrend,j1,j2,jnew,oper,lgt,cut,
     +       minfac,itype,lwrite)
*
*       operates on oldfield (j1:j2,) to make newfield(jnew,)
*       oper = mean|sd|var|min|max|num|sum
*       lgt  = ' '|<|>
*       cut  = cut-off
*       itype=0 scalar data
*       itype=1 direction data
*
*       15-jun-2004 also compute mean when there are missing data
*       2-nov-2005 adapted for fields
*
        implicit none
        integer nx,ny,nperyear,nperyearnew,yrbeg,yrend,j1,j2,jnew,itype
        real oldfield(nx,ny,nperyear,yrbeg:yrend),
     +       newfield(nx,ny,nperyearnew,yrbeg:yrend),
     +       cut(nx,ny,nperyear),minfac
        character oper*3,lgt*1
        logical lvalid(nx,ny,nperyear,yrbeg:yrend),lwrite
        integer yr,mo,i,j,n,jx,jy,ntot,offset,lfirst
        real s,s1,s2,sx,sy,mean,st
*
        if ( lwrite ) then
            print *,'fieldday2period: parameters are '
            print *,'nperyear,nperyearnew= ',nperyear,nperyearnew
            print *,'yrbeg,yrend         = ',yrbeg,yrend
            print *,'j1,j2,jnew          = ',j1,j2,jnew
            print *,'oper,lgt            = ',oper,lgt
            print *,'itype               = ',itype
            if ( lgt.ne.' ' ) then
                print *,'cut(1,1)            = '
                print *,(cut(1,1,j),j=1,nperyear)
            endif
            print *,'oldfield(1,1)        = '
            do yr=yrbeg,yrend,max(1,yrend-yrbeg)
                do mo=j1,j2
                    j = mo
                    call normon(j,yr,i,nperyear)
                    if ( yrbeg.eq.0 .and. yrend.eq.0 ) i = yr
                    if ( i.lt.yrbeg ) cycle
                    if ( oldfield(1,1,j,i).lt.1e33 ) then
                        print *,yr,j,oldfield(1,1,j,i)
                    endif
                enddo
            enddo
        endif
        if ( itype.eq.1 .and. (oper.eq.'min' .or. oper.eq.'max' .or.
     +       oper.eq.'nti' .or. oper.eq.'xti' .or. oper.eq.'fti' .or.
     +       oper.eq.'lti' .or.  oper.eq.'sd' .or. oper.eq.'con') ) then
            write(0,*) 'fieldday2period: cannot take min,max,sd of '//
     +           'vector data'
            write(*,*) 'fieldday2period: cannot take min,max,sd of '//
     +           'vector data'
            call abort
        endif
*
        do yr=yrbeg,yrend
            call keepalive1('Year ',yr-yrbeg+1,yrend-yrbeg+1)
            do jy=1,ny
                do jx=1,nx
                    if ( (oper.eq.'sd' .or. oper.eq.'var' ) .and. 
     +                   itype.eq.0 ) then
*                       precompute the mean for numerical stability
                        n = 0
                        s = 0
                        do mo=j1,j2
                            j = mo
                            call normon(j,yr,i,nperyear)
                            if ( yrbeg.eq.0 .and. yrend.eq.0 ) i = yr
                            if ( i.lt.yrbeg ) cycle
*                           skip Feb 29
                            if ( nperyear.eq.366 .and. j.eq.60 .and. 
     +                           oldfield(jx,jy,j,i).gt.1e33 ) goto 40
*                           but generate an "invalid" for any other
C                           invalid data [this may be relaxed later...]
                            if ( oldfield(jx,jy,j,i).gt.1e33 ) goto 40
                            n = n + 1
                            s = s + oldfield(jx,jy,j,i)
 40                         continue
                        enddo
                        if ( n.gt.1 ) then
                            mean = s/n
                        else
                            mean = 0
                        endif
                    endif
                    if ( oper.eq.'mea' .or. oper.eq.'sd ' .or.
     +                   oper.eq.'var' .or.
     +                   oper.eq.'sum' .or. oper.eq.'bel' .or. 
     +                   oper.eq.'abo' .or. oper.eq.'con' ) then
                        if ( itype.eq.0 ) then
                            s = 0
                        elseif ( itype.eq.1 ) then
                            sx = 0
                            sy = 0
                        else
                            write(0,*) 'fieldday2period: error: unknown
     +                           itype ',itype
                            write(*,*) 'fieldday2period: error: unknown
     +                           itype ',itype
                            call abort
                        endif
                        if ( oper.eq.'sd ' .or. oper.eq.'var' ) s2 = 0
                    elseif ( oper.eq.'min' .or. oper.eq.'nti' ) then
                        if ( itype.eq.0 ) then
                            s = 3e33
                        elseif ( itype.eq.1 ) then
                            sx = 3e33
                            sy = 3e33
                        endif
                    elseif ( oper.eq.'max' .or. oper.eq.'xti' ) then
                        if ( itype.eq.0 ) then
                            s = -3e33
                        elseif ( itype.eq.1 ) then
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
                        if ( yrbeg.eq.0 .and. yrend.eq.0 ) i = yr
                        if ( i.lt.yrbeg .or. i.gt.yrend ) cycle
*                       skip Feb 29
                        if ( nperyear.eq.366 .and. j.eq.60 .and. 
     +                       oldfield(jx,jy,j,i).gt.1e33 ) then
                            offset = -1
                            cycle
                        end if
*                       but generate an "invalid" for any other invalid
*                       data [this may be relaxed later...]
                        if ( oldfield(jx,jy,j,i).gt.1e33 ) then
***                         newfield(jnew,i) = 3e33
                            cycle
                        endif
                        if ( lfirst.eq.9999 ) lfirst = mo-j1
                        if ( lvalid(jx,jy,j,i) ) ntot = ntot + 1
                        if ( lgt.eq.' ' .or.
     +                       lgt.eq.'<' .and. 
     +                       oldfield(jx,jy,j,i).lt.cut(jx,jy,j) .or.
     +                       lgt.eq.'>' .and. 
     +                       oldfield(jx,jy,j,i).gt.cut(jx,jy,j) ) then
                            n = n + 1
                            if ( oper.eq.'fti' ) then
                                if ( st.eq.3e33 ) st = mo - 0.5 + offset
                            else if ( oper.eq.'lti' ) then
                                st = mo - 0.5 + offset
                            else if ( oper.eq.'mea' .or. oper.eq.'sum'
     +                               .or.
     +                               oper.eq.'sd ' .or. oper.eq.'var' )
     +                               then
                                if ( itype.eq.0 ) then
                                    s = s + oldfield(jx,jy,j,i)
                                elseif ( itype.eq.1 ) then
                                    sx = sx + cos(atan(1.)
     +                                   *oldfield(jx,jy,j,i)/45)
                                    sy = sy + sin(atan(1.)
     +                                   *oldfield(jx,jy,j,i)/45)
                                endif
                                if ( oper.eq.'sd ' .or. oper.eq.'var' )
     +                               then
                                    s2 = s2
     +                                   + (oldfield(jx,jy,j,i)-mean)**2
                                endif
                            elseif ( oper.eq.'abo' ) then
                                s = s +oldfield(jx,jy,j,i)-cut(jx,jy,j)
                            elseif ( oper.eq.'bel' ) then
                                s = s +cut(jx,jy,j)-oldfield(jx,jy,j,i)
                            elseif ( oper.eq.'min' .or. oper.eq.'nti' )
     +                               then
                                if ( oldfield(jx,jy,j,i).lt.s ) then
                                    s = oldfield(jx,jy,j,i)
                                    st = mo - 0.5 + offset
                                end if
                            elseif ( oper.eq.'max' .or. oper.eq.'xti')
     +                               then
                                if ( oldfield(jx,jy,j,i).gt.s ) then
                                    s = oldfield(jx,jy,j,i)
                                    st = mo - 0.5 + offset
                                end if
                            else if ( oper.eq.'con' ) then
                                s1 = s1 + 1
                                s = max(s,s1)
                            endif
                            if ( lwrite .and.
     +                           jx.eq.(nx+1)/2 .and. jy.eq.(ny+1)/2 )
     +                           print *,'partial sum: ',j,n,s
                        else
                            if ( oper.eq.'con' ) then
                                s1 = 0
                            end if
                            if ( lwrite .and.
     +                           jx.eq.(nx+1)/2 .and. jy.eq.(ny+1)/2 )
     +                           print *,j,oldfield(jx,jy,j,i),lgt
     +                           ,cut(jx,jy,j)
                        endif
                    enddo       ! j
                    if ( lfirst.gt.minfac*(j2-j1+1) .or. 
     +                   ntot.lt.minfac*nperyear/nperyearnew ) then
                        newfield(jx,jy,jnew,i) = 3e33
                    elseif ( n.eq.0 .and. 
     +                   (oper.eq.'mea' .or. oper.eq.'min' .or. 
     +                       oper.eq.'max' .or. oper.eq.'xti' .or. 
     +                       oper.eq.'nti' .or. oper.eq.'fti' .or. 
     +                       oper.eq.'lti' ) ) then
                        newfield(jx,jy,jnew,i) = 3e33
                    elseif ( oper.eq.'mea' ) then
                        if ( itype.eq.0 ) then
                            if ( lwrite .and.
     +                           jx.eq.(nx+1)/2 .and. jy.eq.(ny+1)/2 )
     +                           print *,'s,n,s/n = ',s,n,s/n
                            newfield(jx,jy,jnew,i) = s/n
                        elseif ( itype.eq.1 ) then
                            if ( lwrite .and.
     +                           jx.eq.(nx+1)/2 .and. jy.eq.(ny+1)/2 )
     +                           then
                                print *,'sx,n,sx/n = ',sx,n,sx/n
                                print *,'sy,n,sy/n = ',sy,n,sy/n
                                print *,'atan2(sy,sx)',atan2(sy,sx)
                            endif
                            newfield(jx,jy,jnew,i) = 45*atan2(sy,sx)
     +                           /atan(1.)
                            if ( newfield(jx,jy,jnew,i).lt.0 )
     +                           newfield(jx,jy,jnew,i) =
     +                           newfield(jx,jy,jnew,i) + 360
                        endif
                    elseif ( oper.eq.'sd ' .or. oper.eq.'var' ) then
                        newfield(jx,jy,jnew,i) = s2/n - (s/n-mean)**2
                        if ( oper.eq.'sd' ) then
                            if ( newfield(jx,jy,jnew,i).ge.0 ) then
                                newfield(jx,jy,jnew,i) = 
     +                               sqrt(newfield(jx,jy,jnew,i))
                            else
                                newfield(jx,jy,jnew,i) = 3e33
                            end if
                        end if
                    elseif ( oper.eq.'sum' .or. oper.eq.'bel' .or. 
     +                       oper.eq.'abo' .or. oper.eq.'min' .or. 
     +                       oper.eq.'max' ) then
                        if ( itype.eq.0 ) then
                            newfield(jx,jy,jnew,i) = s
                        elseif ( itype.eq.1 ) then
                            newfield(jx,jy,jnew,i) = 45*atan2(sy,sx)
     +                           /atan(1.)
                            if ( newfield(jx,jy,jnew,i).lt.0 )
     +                           newfield(jx,jy,jnew,i) =
     +                           newfield(jx,jy,jnew,i) + 360
                        endif
                    elseif ( oper.eq.'num' ) then
                        newfield(jx,jy,jnew,i) = n
                    elseif ( oper.eq.'nti' .or. oper.eq.'xti' .or.
     +                       oper.eq.'fti' .or. oper.eq.'lti' ) then
                        newfield(jx,jy,jnew,i) = st
                    elseif ( oper.eq.'con' ) then
                        newfield(jx,jy,jnew,i) = s
                    else
                        write(0,*) 'fieldday2period: error: '//
     +                       'unknown operation ',oper
                        call abort
                    endif
                    if ( lwrite .and.
     +                   jx.eq.(nx+1)/2 .and. jy.eq.(ny+1)/2 ) then
                        print *,'fieldday2period: computed newfield'
***                     if ( newfield(jx,jy,i,i).lt.1e33 ) then
                        print *,i,jnew,newfield(jx,jy,jnew,i),s,st,itype
***                     endif
                    endif
 100                continue
                enddo           ! jx
            enddo               ! jy
            !!!call keepalive(yr-yrbeg+1,yrend-yrbeg+1)            
        enddo                   ! yr
        end

        subroutine  adjustnames(oper,nperyear,nperyearnew,lgt,pcut
     +       ,punits,lvars,cell_methods)
!
!       adjust names to reflect the operation just performed
!
        implicit none
        integer nperyear,nperyearnew
        real pcut
        character oper*3,lgt*1,punits*(*),lvars*(*),cell_methods*(*)
        integer i
        character prefix*255,postfix*255
        
        call nperyear2string(nperyearnew,prefix)        
        if ( oper.eq.'mea' ) then
            prefix = trim(prefix)//' mean'
        elseif ( oper.eq.'sd ' ) then
            prefix = trim(prefix)//' standard deviation'
        elseif ( oper.eq.'var' ) then
            prefix = trim(prefix)//' variance'
        elseif ( oper.eq.'min' ) then
            prefix = trim(prefix)//' minimum'
        elseif ( oper.eq.'max' ) then
            prefix = trim(prefix)//' maximum'
        elseif ( oper.eq.'nti' ) then
            prefix = 'time of '//trim(prefix)//' minimum'
        elseif ( oper.eq.'xti' ) then
            prefix = 'time of '//trim(prefix)//' maximum'
        elseif ( oper.eq.'fti' ) then
            prefix = 'time of '//trim(prefix)//' first'
        elseif ( oper.eq.'lti' ) then
            prefix = 'time of '//trim(prefix)//' last'
        elseif ( oper.eq.'num' ) then
            prefix = trim(prefix)//' number'
        elseif ( oper.eq.'con' ) then
            prefix = trim(prefix)//' max number of consecutive number'
        elseif ( oper.eq.'bel' .or. oper.eq.'abo' .or. 
     +           oper.eq.'sum' ) then
            prefix = trim(prefix)//' sum '
        else
            write(0,*) 'adjustnames: unknown operation: ',oper
            write(*,*) 'adjustnames: unknown operation: ',oper
            call abort
        endif
        prefix = trim(prefix)//' of'
        i = len_trim(prefix)
        call nperyear2string(nperyear,prefix(i+2:))

        if ( lgt.eq.' ' ) then
            postfix = ' '
        else
            if ( lgt.eq.'<' ) then
                postfix = 'below'
            elseif ( lgt.eq.'>' ) then
                postfix = 'above'
            endif
            write(postfix(7:),'(g20.6,2a)') pcut,' ',punits
        endif
        lvars = trim(prefix)//' '//trim(lvars)//' '//trim(postfix)
        call deletedoublespaces(lvars)
        if ( cell_methods.eq.' ' ) then
            cell_methods = trim(prefix)//' data '//trim(postfix)
        else
            cell_methods = trim(prefix)//' '//trim(postfix)//' '//
     +           trim(cell_methods)
        endif
        call deletedoublespaces(cell_methods)
        end

        subroutine nperyear2string(nperyear,string)
        implicit none
        integer nperyear
        character string*(*)
        if ( nperyear.gt.366 ) then
            write(string,'(i1,a)') nint(24*365.25/nperyear),'-hr'
        elseif ( nperyear.ge.360 ) then
            string = 'daily'
        elseif ( nperyear.gt.12 ) theN
            write(string,'(i3,a)') nint(365.25/nperyear),'-dy'
        elseif ( nperyear.eq.12 ) then
            string = 'monthly'
        elseif ( nperyear.eq.4 ) then
            string = 'seasonal'
        elseif ( nperyear.eq.1 ) then
            string = 'annual'
        else
            write(0,*) 'nperyear2string: error: cannot handle ',nperyear
     +           ,' yet'
            write(*,*) 'nperyear2string: error: cannot handle ',nperyear
     +           ,' yet'
            call abort
        endif
        end

        subroutine deletedoublespaces(string)
        implicit none
        character*(*) string
        integer i,j,k,state
        if ( string.eq.' ' ) return
        if ( string(1:1).eq.' ' ) then
            state = 0
        else
            state = 1
        endif
        i = 1
        do
            if ( state.eq.0 ) then
!               search for more whitespace
                j = i+1
                do while ( string(j:j).eq.' ' )
                    j = j + 1
                    if ( j.gt.len(string) ) return
                end do
                if ( j.gt.i+1 ) then
!                   multiple whitespace => shorten
                    string(i+1:) = string(j:)
                endif
                i = i+1
                state = 1
            else
                j = i+1
                do while ( string(j:j).ne.' ' )
                    j = j + 1
                    if ( j.gt.len(string) ) return
                end do
                i = j
                state = 0
            endif
        end do
        end
