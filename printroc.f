        subroutine probroc(obs,fcst,ntime,nens,thobs,pthobs,thmod,pthmod
     +       ,lsame,lprint,lwrite,area)
*
*       vary probability (number of ensemble members above
*       threshold) to construct the ROC curve
*
        implicit none
        integer ntime,nens
        real obs(ntime),fcst(nens,ntime),thobs,pthobs,thmod,pthmod,area
        logical lsame,lprint,lwrite
        integer itime,n,iens
        real,allocatable :: prob(:)
*
        if ( lsame ) then
            call getcutoff(thmod,pthobs,obs,1,1,1,ntime,1,ntime
     +           ,1,1,0)
        elseif ( pthmod.gt.0 .and. pthmod.lt.100 ) then
            call getenscutoff(thmod,pthmod,fcst,1,1,1,ntime
     +           ,nens-1,0,nens-1,1,ntime,1,1,0)
            if ( lprint ) then
                write(0,'(a,f6.2,a,g14.6,a)')
     +               'Threshold at ',pthmod,'% corresponds to model '
     +               ,thmod,'<br>'
            endif
        endif
        allocate(prob(ntime))
        do itime=1,ntime
            prob(itime) = 0
            n = 0
            do iens=1,nens
                if ( fcst(iens,itime).lt.1e30 .and.
     +               fcst(iens,itime).ne.-999.9 ) then
                    n = n + 1
                    if ( fcst(iens,itime).gt.thmod ) then
                        prob(itime) = prob(itime) + 1
                    endif
                endif
            enddo
            if ( n.gt.0 ) then
                prob(itime) = prob(itime)/n
            else
                prob(itime) = 3e33
            endif
        enddo
        call printroc(obs,prob,ntime,1,thobs,pthobs,lprint,lwrite,area)
        deallocate(prob)
        end

        subroutine printroc(obs,fcst,ntime,nens,thobs,pthobs
     +       ,lprint,lwrite,area)
*
*       compute the ROC curve and area under it
*       home-grown, I hope fast, algorithm without any binning
*
        implicit none
        integer ntime,nens
        real obs(ntime),fcst(nens,ntime),thobs,pthobs,area
        logical lprint,lwrite
        integer i,itime,iens,iabove,nabove,ibelow,nbelow,nhit,nmis
     +       ,oldhit,oldmis
        real old
        real,allocatable :: above(:),below(:)
        character line*80
        logical lweb
        integer,save :: noutside,nprint
        data noutside,nprint /0,2/
*
*       if needed, convert percentage thobs into aboslute one
*
        if ( pthobs.gt.0 .and. pthobs.lt.100 ) then
            call getcutoff(thobs,pthobs,obs,1,1,1,ntime,1,ntime
     +           ,1,1,0)
            if ( lprint ) write(0,'(a,f6.2,a,g14.6,a)') 'Thobs at '
     +           ,pthobs,'% corresponds to observations ',thobs,
     +           '<br>'
        endif
        if ( thobs.gt.1e33 ) then
            if ( lwrite ) print *,'printroc: thrsehold undefined'
            area = 3e33
            return
        endif
*
*       first make two lists of forecasts for observations 
*       that were above the threshold and below the thobs
*
        allocate(above(nens*ntime+1))
        allocate(below(nens*ntime+1))
        nbelow = 0
        nabove = 0
        do itime=1,ntime
            if ( obs(itime).lt.1e30  .and. 
     +           abs(obs(itime)+999.9).gt.0.01 ) then
                if ( obs(itime).lt.thobs ) then
                    do iens=1,nens
                        if ( fcst(iens,itime).lt.1e30 .and.
     +                       abs(fcst(iens,itime)+999.9).gt.0.01 ) then
                            nbelow = nbelow + 1
                            below(nbelow) = fcst(iens,itime)
                        endif
                    enddo
                    if ( lwrite ) print *,'below: ',(fcst(iens,itime),
     +                   iens=1,min(nens,4))
                else
                    do iens=1,nens
                        if ( fcst(iens,itime).lt.1e30 .and.
     +                       abs(fcst(iens,itime)+999.9).gt.0.01 ) then
                            nabove = nabove + 1
                            above(nabove) = fcst(iens,itime)
                        endif
                    enddo
                    if ( lwrite ) print *,'above: ',(fcst(iens,itime),
     +                   iens=1,min(nens,4))
                endif
            endif
        enddo
        if ( lwrite ) then
            print *,'nbelow = ',nbelow
            print *,'nabove = ',nabove
        endif
        if ( nbelow.eq.0 .or. nabove.eq.0 ) then
            noutside = noutside + 1
            if ( noutside.eq.1 ) then
                call getenv('REMOTE_ADDR',line)
                if ( line.ne.' ' ) then
                    lweb = .true.
                else
                    lweb = .false.
                endif
                write(0,*) 'printroc: cut-off (defined over all ',
     +               'data) ',thobs,' is outside the range ',
     +               'observed during the forecast period<br>'
                if ( .false. ) then
                    if ( lweb ) then
                        write(0,*) '<table>'
                        do itime=1,ntime
                            write(0,*) '<tr><td>',obs(itime)
                            do iens=1,nens
                                write(0,*) '<td>',fcst(iens,itime)
                            end do
                        enddo
                        write(0,*) '</table>'
                    else
                        do itime=1,ntime
                            write(0,*) obs(itime),
     +                           (fcst(iens,itime),iens=1,nens)
                        enddo
                    endif
                endif
            elseif ( noutside.ge.nprint ) then
                nprint = 2*nprint
                write(0,*) 'printroc: cut-off outside range at '
     +               ,noutside,' points<br>'
            endif
            area = 3e33
            return
        endif
*
*       sort the lists from small to large
*
        call nrsort(nbelow,below)
        below(nbelow+1) = 2*abs(below(nbelow)) ! something bigger
        call nrsort(nabove,above)
        above(nabove+1) = 2*abs(above(nabove))
*
*       and peel them off
*
        ibelow = 1
        iabove = 1
!       upper right-hand corner
        nhit = nabove
        nmis = nbelow
        oldhit = nabove
        oldmis = nbelow
        old = 3e33
        area = 0
 100    continue
 1000   format(2f8.5,2g14.6)
        if ( lwrite ) print *,iabove,above(iabove),ibelow,below(ibelow)
        if ( above(iabove).ne.old .and. below(ibelow).ne.old ) then
            if ( lprint ) print 1000,nmis/real(nbelow),nhit/real(nabove)
     +           ,old,min(above(iabove),below(ibelow))
            old = min(above(iabove),below(ibelow))
            area = area + ((nbelow-oldmis)+(0.5*(oldmis-nmis)))
     +           *(oldhit-nhit)
            if ( lwrite ) print *,'area = ',area,nmis,oldmis,nhit,oldhit
            oldhit = nhit
            oldmis = nmis
        endif
!       switch the forecasts one by one from above to below
        if ( iabove.lt.nabove. and. 
     +       ( ibelow.gt.nbelow .or. above(iabove).lt.below(ibelow) ) )
     +       then
!           predicted above: one hit fewer
            iabove = iabove + 1
            nhit = nhit - 1
            goto 100
        elseif ( iabove.lt.nabove ) then
!           predicted below: one false alarm fewer
            ibelow = ibelow + 1
            nmis = nmis - 1
            goto 100
        endif
!       last point
        if ( lprint ) print 1000,0.,0.,old,1.
        area = area + ((nbelow-oldmis)+(0.5*(oldmis)))*(oldhit)
        if ( lwrite ) print *,'area = ',area,0,oldmis,0,oldhit
        area = area/(real(nabove)*nbelow)
        if ( lprint ) print '(a,f5.3)','# area = ',area
        end

        subroutine readthreshold(string,th,pth)
        implicit none
        character string*(*)
        real th,pth
        integer i
        i = index(string,'%')
        if ( i.ne.0 ) then
            th = 3e33
            read(string(:i-1),*) pth
        else
            pth = 3e33
            read(string,*) th
        endif
        end
