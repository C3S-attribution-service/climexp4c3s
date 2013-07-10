*  #[ normsvd1:
        subroutine normsvd1(svdxy,nxf,nyf,nzf,nx,ny,nz,n,normalization)
*
*       normalize the SVDs according to normalization
*       1: largest point value is one (this makes it easiest to plot)
*       2: spatial variance is one
*
        implicit none
        integer nxf,nyf,nzf,nx,ny,nz,n,normalization
        real svdxy(nxf,nyf,nzf,n)
        integer i,j,k,ii,n1
        real xmax,xmin,scale,s1,s2,mean
*       
        do ii=1,n
            if ( normalization.eq.1 ) then
                xmax = -1e30
                xmin = +1e30
                do j=1,ny
                    do i=1,nx
                        do k=1,nz
                            if ( svdxy(i,j,k,ii).lt.1e30 ) then
                                xmin = min(xmin,svdxy(i,j,k,ii))
                                xmax = max(xmax,svdxy(i,j,k,ii))
                            endif
                        end do
                    end do
                end do
                if ( -xmin.gt.xmax  ) then
                    scale = xmin
                else
                    scale = xmax
                end if
            else if ( normalization.eq.2 ) then
                s1 = 0
                n1 = 0
                do j=1,ny
                    do i=1,nx
                        do k=1,nz
                            if ( svdxy(i,j,k,ii).lt.1e30 ) then
                                n1 = n1 + 1
                                s1 = s1 + svdxy(i,j,k,ii)
                            end if
                        end do
                    end do
                end do
                mean = s1/n1
                s1 = 0
                s2 = 0
                do j=1,ny
                    do i=1,nx
                        do k=1,nz
                            if ( svdxy(i,j,k,ii).lt.1e30 ) then
                                s1 = s1 + (svdxy(i,j,k,ii) - mean)
                                s2 = s2 + (svdxy(i,j,k,ii) - mean)**2
                            end if
                        end do
                    end do
                end do
                s2 = s2/n1 - (s1/n1)**2
                if ( s2.gt.0 ) then
                    scale = sqrt(s2)
                else
                    write(0,*) 'normsvd1: error: negative variance',s2
                    call abort
                endif
            else
                write(0,*) 'normsvd1: error: normalization '
     +               ,normalization,' is not yet supported'
                write(*,*) 'normsvd1: error: normalization '
     +               ,normalization,' is not yet supported'
                call abort
            endif
            do j=1,ny
                do i=1,nx
                    do k=1,nz
                        if ( svdxy(i,j,k,ii).lt.1e30 ) then
                            svdxy(i,j,k,ii) = svdxy(i,j,k,ii)/scale
                        end if
                    end do
                end do
            end do
        end do
        end
*  #] normsvd1:
*  #[ normsvd2:
        subroutine normsvd2(pc,svdxy,nperyear,yr1,yr2,j1,j2,
     +       nens1,nens2,neigen,nxf,nyf,nzf,nx,ny,nz,normalization,
     +       lnomissing,lwrite)
*
*       normalize the SVDs according to normalization
*       -2: time series variance is one
*
        implicit none
        integer nperyear,yr1,yr2,j1,j2,nens1,nens2,neigen,
     +       nxf,nyf,nzf,nx,ny,nz,normalization
        real pc(nperyear,yr1:yr2,0:nens2,neigen),
     +       svdxy(nxf,nyf,nzf,neigen)
        logical lnomissing,lwrite
        integer i,j,k,ii,n1,iens,mo,yr,m,y
        real xmin,xmax,scale,s1,s2,mean
*
        if ( lwrite ) then
            print *,'normsvd2: input'
            print *,'nperyear,yr1,yr2,j1,j2 = ',nperyear,yr1,yr2,j1,j2
            print *,'nens1,nens2,neigen     = ',nens1,nens2,neigen
            print *,'nxf,nyf,nzf,nx,ny,nz   = ',nxf,nyf,nzf,nx,ny,nz
            print *,'normalization          = ',normalization
        end if
*
        do ii=1,neigen
            if ( normalization.eq.-1 ) then
                xmax = -1e30
                xmin = +1e30
                do iens=nens1,nens2
                    do y=yr1-1,yr2
                        do m=j1,j2
                            mo=m
                            call normon(mo,y,yr,nperyear)
                            if ( yr.lt.yr1 .or. yr.gt.yr2 ) cycle
                            if ( pc(mo,yr,iens,ii).lt.1e30 ) then
                                xmin=min(xmin,pc(mo,yr,iens,ii))
                                xmax=max(xmax,pc(mo,yr,iens,ii))
                            end if
                        end do
                    end do
                end do
                if ( -xmin.gt.xmax  ) then
                    scale = xmin
                else
                    scale = xmax
                endif
            else if ( normalization.eq.-2 ) then
                s1 = 0
                n1 = 0
                do iens=nens1,nens2
                    do y=yr1-1,yr2
                        do m=j1,j2
                            mo=m
                            call normon(mo,y,yr,nperyear)
                            if ( yr.lt.yr1 .or. yr.gt.yr2 ) cycle
                            if ( pc(mo,yr,iens,ii).lt.1e30 ) then
                                n1 = n1 + 1
                                s1 = s1 + pc(mo,yr,iens,ii)
                            end if
                        end do
                    end do
                end do
                mean = s1/n1
                s1 = 0
                s2 = 0
                do iens=nens1,nens2
                    do y=yr1-1,yr2
                        do m=j1,j2
                            mo=m
                            call normon(mo,y,yr,nperyear)
                            if ( yr.lt.yr1 .or. yr.gt.yr2 ) cycle
                            if ( pc(mo,yr,iens,ii).lt.1e30 ) then
                                s1 = s1 + (pc(mo,yr,iens,ii) -
     +                               mean)
                                s2 = s2 + (pc(mo,yr,iens,ii) -
     +                               mean)**2
                            end if
                        end do
                    end do
                end do
                s2 = s2/n1 - (s1/n1)**2
                if ( s2.gt.0 ) then
                    scale = sqrt(s2)
                else
                    write(0,*) 'normsvd2: error: negative variance',s2
                    call abort
                end if
            else
                write(0,*) 'normsvd2: error: normalization '
     +               ,normalization,' is not yet supported'
                write(*,*) 'normsvd2: error: normalization '
     +               ,normalization,' is not yet supported'
                call abort
            end if
            do iens=nens1,nens2
                do y=yr1-1,yr2
                    do m=j1,j2
                        mo=m
                        call normon(mo,y,yr,nperyear)
                        if ( yr.lt.yr1 .or. yr.gt.yr2 ) cycle
                        if ( pc(mo,yr,iens,ii).lt.1e30 ) then
                            pc(mo,yr,iens,ii) = 
     +                           pc(mo,yr,iens,ii)/scale
                        end if
                    end do
                end do
            end do
            do k=1,nz
                do j=1,ny
                    do i=1,nx
                        if ( svdxy(i,j,k,ii).lt.1e30 ) then
                            svdxy(i,j,k,ii) = svdxy(i,j,k,ii)*scale
                        end if
                    end do
                end do
            end do
        end do
        end
*  #] normsvd2:
*  #[ normeof1:

        subroutine normeof1(eofxy,nxf,nyf,nx,ny,n,normalization)
*
*       normalize the EOFs according to normalization
*       1: largest point value is one (this makes it easiest to plot)
*       2: spatial variance is one
*
        implicit none
        integer nxf,nyf,nx,ny,n,normalization
        real eofxy(nxf,nyf,n)
        integer i,j,ii,n1
        real xmax,xmin,scale,s1,s2,mean
*       
        do ii=1,n
            if ( normalization.eq.1 ) then
                xmax = -1e30
                xmin = +1e30
                do j=1,ny
                    do i=1,nx
                        if ( eofxy(i,j,ii).lt.1e30 ) then
                            xmin = min(xmin,eofxy(i,j,ii))
                            xmax = max(xmax,eofxy(i,j,ii))
                        endif
                    enddo
                enddo
                if ( -xmin.gt.xmax  ) then
                    scale = xmin
                else
                    scale = xmax
                endif
            elseif ( normalization.eq.2 ) then
                s1 = 0
                n1 = 0
                do j=1,ny
                    do i=1,nx
                        if ( eofxy(i,j,ii).lt.1e30 ) then
                            n1 = n1 + 1
                            s1 = s1 + eofxy(i,j,ii)
                        endif
                    enddo
                enddo
                mean = s1/n1
                s1 = 0
                s2 = 0
                do j=1,ny
                    do i=1,nx
                        if ( eofxy(i,j,ii).lt.1e30 ) then
                            s1 = s1 + (eofxy(i,j,ii) - mean)
                            s2 = s2 + (eofxy(i,j,ii) - mean)**2
                        endif
                    enddo
                enddo
                s2 = s2/n1 - (s1/n1)**2
                if ( s2.gt.0 ) then
                    scale = sqrt(s2)
                else
                    write(0,*) 'normeof1: error: negative variance',s2
                    call abort
                endif
            else
                write(0,*) 'normeof1: error: normalization '
     +               ,normalization,' is not yet supported'
                write(*,*) 'normeof1: error: normalization '
     +               ,normalization,' is not yet supported'
                call abort
            endif
            do j=1,ny
                do i=1,nx
                    if ( eofxy(i,j,ii).lt.1e30 ) then
                        eofxy(i,j,ii) = eofxy(i,j,ii)/scale
                    endif
                enddo
            enddo
        enddo
        end
*  #] normeof1:
*  #[ normeof2:
        subroutine normeof2(pc,eofxy,nperyear,yr1,yr2,j1,j2,
     +       nens1,nens2,neigen,nxf,nyf,nx,ny,normalization,
     +       lnomissing,lwrite)
*
*       normalize the EOFs according to normalization
*       -2: time series variance is one
*
        implicit none
        integer nperyear,yr1,yr2,j1,j2,nens1,nens2,neigen,nxf,nyf,nx,ny
     +       ,normalization
        real pc(nperyear,yr1:yr2,0:nens2,neigen),
     +       eofxy(nxf,nyf,neigen)
        logical lnomissing,lwrite
        integer i,j,ii,n1,iens,mo,yr,m,y
        real xmin,xmax,scale,s1,s2,mean
*
        if ( lwrite ) then
            print *,'normeof2: input'
            print *,'nperyear,yr1,yr2,j1,j2 = ',nperyear,yr1,yr2,j1,j2
            print *,'nens1,nens2,neigen     = ',nens1,nens2,neigen
            print *,'nxf,nyf,nx,ny          = ',nxf,nyf,nx,ny
            print *,'normalization          = ',normalization
        endif
*       
        do ii=1,neigen
            if ( normalization.eq.-1 ) then
                xmax = -1e30
                xmin = +1e30
                do iens=nens1,nens2
                    do y=yr1-1,yr2
                        do m=j1,j2
                            mo=m
                            call normon(mo,y,yr,nperyear)
                            if ( yr.lt.yr1 .or. yr.gt.yr2 ) cycle
                            if ( pc(mo,yr,iens,ii).lt.1e30 ) then
                                xmin=min(xmin,pc(mo,yr,iens,ii))
                                xmax=max(xmax,pc(mo,yr,iens,ii))
                            endif
                        enddo
                    enddo
                enddo
                if ( -xmin.gt.xmax  ) then
                    scale = xmin
                else
                    scale = xmax
                endif
            elseif ( normalization.eq.-2 ) then
                s1 = 0
                n1 = 0
                do iens=nens1,nens2
                    do y=yr1-1,yr2
                        do m=j1,j2
                            mo=m
                            call normon(mo,y,yr,nperyear)
                            if ( yr.lt.yr1 .or. yr.gt.yr2 ) cycle
                            if ( pc(mo,yr,iens,ii).lt.1e30 ) then
                                n1 = n1 + 1
                                s1 = s1 + pc(mo,yr,iens,ii)
                            endif
                        enddo
                    enddo
                enddo
                mean = s1/n1
                s1 = 0
                s2 = 0
                do iens=nens1,nens2
                    do y=yr1-1,yr2
                        do m=j1,j2
                            mo=m
                            call normon(mo,y,yr,nperyear)
                            if ( yr.lt.yr1 .or. yr.gt.yr2 ) cycle
                            if ( pc(mo,yr,iens,ii).lt.1e30 ) then
                                s1 = s1 + (pc(mo,yr,iens,ii) -
     +                               mean)
                                s2 = s2 + (pc(mo,yr,iens,ii) -
     +                               mean)**2
                            endif
                        enddo
                    enddo
                enddo
                s2 = s2/n1 - (s1/n1)**2
                if ( s2.gt.0 ) then
                    scale = sqrt(s2)
                else
                    write(0,*) 'normeof2: error: negative variance',s2
                    call abort
                endif
            else
                write(0,*) 'normeof2: error: normalization '
     +               ,normalization,' is not yet supported'
                write(*,*) 'normeof2: error: normalization '
     +               ,normalization,' is not yet supported'
                call abort
            endif
            do iens=nens1,nens2
                do y=yr1-1,yr2
                    do m=j1,j2
                        mo=m
                        call normon(mo,y,yr,nperyear)
                        if ( yr.lt.yr1 .or. yr.gt.yr2 ) cycle
                        if ( pc(mo,yr,iens,ii).lt.1e30 ) then
                            pc(mo,yr,iens,ii) = 
     +                           pc(mo,yr,iens,ii)/scale
                        endif
                    enddo
                enddo
            enddo
            do j=1,ny
                do i=1,nx
                    if ( eofxy(i,j,ii).lt.1e30 ) then
                        eofxy(i,j,ii) = eofxy(i,j,ii)*scale
                    endif
                enddo
            enddo
        enddo
        end
*  #] normeof2:
