        subroutine latlonint(temp,nt,ntimes,nlat,nlon,lats,lons,iiwrite
     +        ,tmin,tmax,dt)
*
*       interpolate latitude and longitude
*       latitude: if there is a measurement within maxlat degree
*       longitude: if there is a measurement within maxlon degrees
*       or if the temperature difference is less than dt.
*
	implicit none
	integer ntimes,nt,nlat,nlon,iiwrite
	real lats(nlat),lons(nlon)
	real temp(nlon,nlat,nt),tmin,tmax,dt
        real maxtemp
	integer i,j,k,l,n,statelat,statelon,lastlat,nextlat
     +	      ,lastlon,nextlon,ii(2),jj(2),i1,j1,iwrite,iwritesav
        character idigit*1
	integer maxlat,maxlon
        logical again
	data maxlat,maxlon /2,10/
*
*       init
*
        iwrite = iiwrite        ! because we will change it
        if ( tmin.ne.0 ) then
            i = log10(abs(tmin))
        else
            i = 0
        endif
        if ( tmin.lt.0 ) i = i + 1
        if ( tmax.ne.0 ) then
            j = log10(abs(tmax))
        else
            j = 0
        endif
        write(idigit,'(i1)') 5+max(i,j)
*       
*       loop over time, depths
*
        do k=1,ntimes
*       
*           any valid values?
            n = 0
            do i=1,nlon
                do j=1,nlat
                    if ( temp(i,j,k).lt.1e33 ) then
                        n = n + 1
                        if ( n.le.2 ) then
                            ii(n) = i
                            jj(n) = j
                        endif
                    endif
                enddo
            enddo
            if ( iwrite.ge.2 ) print '(a,i8,a,i4,i3,a,i3)'
     +            ,'latlonint: found ',n,' valid points at time ',k
*           only one point
            if ( n.le.1 ) goto 800
*	    adjacent points
            if ( n.eq.2 .and.
     +            jj(1).eq.jj(2) .and. abs(ii(1)-ii(2)).eq.1 .or. 
     +            ii(1).eq.ii(2) .and. abs(jj(1)-jj(2)).eq.1 ) goto 800
*
*	    debug output
            if ( iwrite.ge.2 ) then
                do j=nlat,1,-1
                    print '(i4,100f'//idigit//'.2)',
     +                   nint(lats(j)),(temp(i,j,k),i=1,nlon)
                enddo
                print '(4x,100i'//idigit//')',(nint(lons(i)),i=1,nlon)
            endif
*
*	    there are more missing longitudes
            n = 0
            do j=1,nlat
*               to guard against errors
                lastlon = 1
                nextlon = nlon
                statelon = 0
                do i=1,nlon
*                   to guard against errors
                    lastlat = 1
                    nextlat = nlat
*		    search for interval in longitude
                    if ( temp(i,j,k).lt.1e33 ) then
*			valid point
                        statelon = 1
                        lastlon = i
                        nextlon = i
                        if ( iwrite.ge.4 ) print *,'valid point'
                    elseif ( statelon.eq.0 ) then
*			no interpolation possible
                        statelon = 0
                        if ( iwrite.ge.4 ) print *
     +                        ,'no lon interpolation'
                    elseif ( statelon.eq.1 ) then
*			invalid point, search for interval
                        do nextlon=i+1,nlon
                            if ( temp(nextlon,j,k).lt.1e33 ) then
                                statelon = 2
                                if ( iwrite.ge.4 ) print *
     +                                ,'lon interval ',lastlon,nextlon
                                goto 110
                            endif
                        enddo
*                       to avoid errors later on
                        nextlon = nlon
                        statelon = 0
                        if ( iwrite.ge.4 ) print *
     +                        ,'no lon interpolation'
  110                   continue
                    endif
                    statelat = 0
                    if ( statelon.eq.0 .or. statelon.eq.2 ) then
*			search for interval in latitude
*			(for each point in longitude anew)
                        do lastlat=j-1,1,-1
                            if ( temp(i,lastlat,k).lt.1e33 ) then
*				valid point
                                statelat = 1
                                if ( iwrite.ge.4 ) print *
     +                                ,'found lastlat ',lastlat
                                goto 210
                            endif
                        enddo
                        lastlat = nlat
                        statelat = 0
                        if ( iwrite.ge.4 ) print *,'no lastlat'
  210                   continue
                        if ( statelat.ne.0 ) then
                            do nextlat=j+1,nlat,1
                                if ( temp(i,nextlat,k).lt.1e33 ) then
*				    valid point
                                    statelat = 2
                                    if ( iwrite.ge.4 ) print *
     +                                    ,'found nextlat ',nextlat
                                    goto 220
                                endif
                            enddo
                            nextlat = nlat
                            statelat = 0
                            if ( iwrite.ge.4 ) print *,'no nextlat'
  220                       continue
                        endif   ! found firstlat
                    endif       ! invalid point
*
*		    interpolate!
*
                    again = .FALSE.
  300               continue    ! comefrom: recomputation with debugging
                    if ( statelon.eq.2 .and. (
     +                    lons(i)-lons(lastlon).le.maxlon .or.
     +                    lons(nextlon)-lons(i).le.maxlon .or.
     +                    abs(temp(lastlon,j,k) - temp(nextlon,j,k))
     +                    .le.dt )) then
                        if ( statelat.eq.2 .and. (
     +                        lats(j)-lats(lastlat).le.maxlat .or.
     +                        lats(nextlat)-lats(j).le.maxlat .or.
     +                        abs(temp(i,lastlat,k) - temp(i,nextlat,k))
     +                        .le.dt ) ) then
*
*			    interpolate both latitude and longitude
*                           weighted by maxlon:maxlat
                            n = n + 1
                            temp(i,j,k) = ( (
     +                            (lats(nextlat) - lats(j))
     +                            *temp(i,lastlat,k) + 
     +                            (lats(j) - lats(lastlat))
     +                            *temp(i,nextlat,k) )/maxlat + (
     +                            (lons(nextlon) - lons(i))
     +                            *temp(lastlon,j,k) + 
     +                            (lons(i) - lons(lastlon))
     +                            *temp(nextlon,j,k) )/maxlon )/(
     +                            (lons(nextlon)-lons(lastlon))
     +                            /real(maxlon) + 
     +                            (lats(nextlat)-lats(lastlat))
     +                            /real(maxlat) )
                            if ( iwrite.ge.3 ) then
                                print *,'latlonint: interpolating both '
     +                                ,i,j
                                print '(a,i3,i4,f8.2)',
     +                               'last ',lastlat,nint(lats(lastlat))
     +                               ,temp(i,lastlat,k),
     +                               ' new ',j,nint(lats(j)),temp(i,j,k)
     +                               ,
     +                               'next ',nextlat,nint(lats(nextlat))
     +                               ,temp(i,nextlat,k),
     +                               'last ',lastlon,nint(lons(lastlon))
     +                               ,temp(lastlon,j,k),
     +                               ' new ',i,nint(lons(i)),temp(i,j,k)
     +                               ,
     +                               'next ',nextlon,nint(lons(nextlon))
     +                               ,temp(nextlon,j,k)
                            endif
                        else    ! latitude interpolatable?
*
*			    interpolate only longitude
                            n = n + 1
                            temp(i,j,k) = ( 
     +                            (lons(nextlon) - lons(i))
     +                            *temp(lastlon,j,k) + 
     +                            (lons(i) - lons(lastlon))
     +                            *temp(nextlon,j,k) )/
     +                            (lons(nextlon)-lons(lastlon))
                            if ( iwrite.ge.3 ) then
                                print *,'latlonint: '//
     +                                'interpolating longitude',i,j
                                print '(a,i3,i4,f8.2)',
     +                               'last ',lastlon,nint(lons(lastlon))
     +                                ,temp(lastlon,j,k),
     +                               ' new ',i,nint(lons(i)),
     +                               temp(i,j,k),
     +                               'next ',nextlon,nint(lons(nextlon))
     +                               ,temp(nextlon,j,k)
                            endif
                        endif   ! latitude interpolatable?
                    else        ! longitude interplatable?
                        if ( statelat.eq.2 .and. (
     +                        lats(j)-lats(lastlat).le.maxlat .or.
     +                        lats(nextlat)-lats(j).le.maxlat .or.
     +                        abs(temp(i,lastlat,k) - temp(i,nextlat,k))
     +                        .le.dt ) ) then
*
*			    interpolate only latitude
                            n = n + 1
                            temp(i,j,k) = ( 
     +                            (lats(nextlat) - lats(j))
     +                            *temp(i,lastlat,k) + 
     +                            (lats(j) - lats(lastlat))
     +                            *temp(i,nextlat,k) )/
     +                            (lats(nextlat)-lats(lastlat))
                            if ( iwrite.ge.3 ) then
                                print *,'latlonint: '//
     +                                'interpolating latitude',i,j
                                print '(a,i3,i4,f8.2)',
     +                               'last ',lastlat,nint(lats(lastlat))
     +                               ,temp(i,lastlat,k),
     +                               ' new ',j,nint(lats(j)),
     +                               temp(i,j,k),
     +                               'next ',nextlat,nint(lats(nextlat))
     +                               ,temp(i,nextlat,k)
                            endif
                        endif   ! latitude interpolatable
                    endif       ! longitude interpolatable
                    if (  temp(i,j,k).ne.3e33 .and.(
     +                    temp(i,j,k).lt.tmin .or.
     +                    temp(i,j,k).gt.tmax ) ) then
                        if ( again ) then
                            print *,'==========='
                            iwrite = iwritesav
                            again = .FALSE.
                        else    ! again?
                            print *,'latlonint: error: T=',temp(i,j,k)
     +                            ,' recomputing with debugging'
                            iwritesav = iwrite
                            iwrite = 4
                            again = .TRUE.
                            goto 300
                        endif   ! again?
                    endif       ! funny temperature
                enddo           ! longitudes
            enddo               ! latitudes
*
*	    debug output
            if ( n.gt.0 .and. iwrite.ge.2 ) then
                do j=nlat,1,-1
                    print '(i4,100f'//idigit//'.2)',
     +                    nint(lats(j)),(temp(i,j,k),i=1,nlon)
                enddo
                print '(4x,100i'//idigit//')',(nint(lons(i)),i=1,nlon)
            endif
*
*		    comefrom: any valid values?
  800       continue
        enddo                   ! loop over fields
	end

