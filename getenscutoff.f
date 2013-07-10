        subroutine getenscutoff(cut,pcut,dat,npermax,nperyear,yrbeg
     +        ,yrend,nensmax,nens1,nens2,yr1,yr2,j1,j2,lag)
*
*       compute the cut-off in absolute units given a cut-off as a
*       percentage (pcut), data and the range over which the data 
*       should be considered
*       14-jan-2002
*       If the data is daily, a 5-day window is considered.
*
        implicit none
        integer nmax
        integer npermax,nperyear,yrbeg,yrend,yr1,yr2,nensmax,nens1,nens2
     +        ,j1,j2,lag
        real cut,pcut,dat(npermax,yrbeg:yrend,0:nensmax)
*
        integer yr,jj,m,j,i,ii,n,k,kdif,iens
        real,allocatable :: a(:)
        logical lwrite
        parameter (lwrite=.false.)
*
*       for daily data, smooth by considering two days on either side
        if ( nperyear.eq.366 ) then
            kdif = 2
        else
            kdif = 0
        endif
*
*       get linear array with data
        nmax = 12000
 10     continue
        allocate(a(nmax))
        n = 0
        do iens=nens1,nens2
            do i=yr1,yr2
                do j=j1,j2
                    do k=j-kdif,j+kdif
                        m = k-lag
                        call normon(m,i,ii,nperyear)
                        if ( ii.lt.yrbeg .or.ii.gt.yrend ) goto 710
                        if ( dat(m,ii,iens).lt.1e33 ) then
                            n = n+1
                            if ( n.gt.nmax ) then
                                deallocate(a)
                                nmax = nmax*2
                                goto 10
                            endif
                            a(n) = dat(m,ii,iens)
                        endif
  710                   continue
                    enddo
                enddo
            enddo
        enddo
        call getcut(cut,pcut,n,a)
        if ( lwrite ) then
            print *,'getcutoff: pcut          = ',pcut
            print *,'           yr1,yr2,j1,j2 = ',yr1,yr2,j1,j2
            print *,'           i,n           = ',i,n
            print *,'           cut           = ',cut
        endif
        deallocate(a)
        end
