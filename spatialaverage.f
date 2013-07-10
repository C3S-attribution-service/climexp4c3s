        subroutine spatialaverage(field,xx,yy,nx,ny,nperyear,firstyr
     +       ,lastyr,avex,avey,lwrite)
!
!       take (in-place) spatial average of field
!
        implicit none
        integer nx,ny,nperyear,firstyr,lastyr,avex,avey
        real field(nx,ny,nperyear,firstyr:lastyr),xx(nx),yy(ny)
        logical lwrite
        integer yr,mo,i,j,ii,jj,iii,jjj,n
        real s,ss,minfac
        minfac = 0.45
!
!       nothing to do?
!
        if ( avex.eq.1 .and. avey.eq.1 ) return
!
!       prevent crashes later on...
!
        if ( avex.gt.nx ) avex = nx
        if ( avey.gt.ny ) avey = ny
!
!       average field
!
        do yr=firstyr,lastyr
            do mo=1,nperyear
                if ( lwrite ) print *,'spatialverage: averaging yr,mo '
     +               ,yr,mo
                do i=1,nx/avex
                    do j=1,ny/avey                        
                        n = 0
                        s = 0
                        do ii=1,avex
                            do jj=1,avey                                
                                iii = (i-1)*avex+ii
                                jjj = (j-1)*avey+jj
                                if ( iii.gt.nx ) cycle
                                if ( jjj.gt.ny ) cycle
                                ss = field(iii,jjj,mo,yr)
                                if ( ss.lt.1e33 ) then
                                    n = n + 1
                                    s = s + ss
                                endif
                            enddo
                        enddo
                        if ( n.gt.minfac*avex*avey ) then
                            field(i,j,mo,yr) = s/n
                        else
                            field(i,j,mo,yr) = 3e33
                        endif
                    enddo
                enddo
            enddo
        enddo
*
*       new axes
*
        do i=1,nx/avex
            s = 0
            do ii=1,avex
                s = s + xx((i-1)*avex+ii)
            enddo
            xx(i) = s/avex
        enddo
        nx = nx/avex
*
        do j=1,ny/avey
            s = 0
            do jj=1,avey
                s = s + yy((j-1)*avey+jj)
            enddo
            yy(j) = s/avey
        enddo
        ny = ny/avey
        if ( lwrite ) then
            print *,'new x-axis ',nx,xx(1),xx(min(2,nx)),'...',xx(nx)
            print *,'new y-axis ',ny,yy(1),yy(min(2,ny)),'...',yy(ny)
        end if
!
        end
