        subroutine filllinarray(dindx,ddata,lfirst,dddata,ndata,n,j1,j2
     +       ,lag,k,nperyear,imens,indxmx,indx,data,npermax,yrbeg,yrend
     +       ,nensmax,filter,yrstart,yrstop,yrmo)
*
*       fill linear arrays without absent values
*
        implicit none
#include "getopts.inc"
        integer ndata,n,j1,j2,lag,k,nperyear,npermax,yrbeg,yrend,nensmax
     +        ,indxmx,yrstart,yrstop,yrmo(2,ndata)
        integer imens(0:indxmx)
        real dindx(ndata),ddata(ndata),dddata(ndata)
        logical lfirst(ndata)
        real indx(npermax,yrbeg:yrend,0:nensmax,indxmx),
     +        data(npermax,yrbeg:yrend,0:nensmax)
        real filter(100)
        integer yr,jj,j,i,m,ii,iiens,iens,jens,if,jm,jp,im,ip
        logical lastvalid
*
        if ( lwrite ) then
            print *,'filllinarray: nens1,nens2,imens(0),imens(k) = ',
     +           nens1,nens2,imens(0),imens(k)
            print *,'              yr1,yr2 = ',yr1,yr2
            print *,'              j1,j2 = ',j1,j2
            if ( .false. ) then
                do iens=0,0
                    print *,'data = ',iens
                    call printdatfile(6,data(1,yrbeg,iens),npermax
     +                   ,nperyear,yrbeg,yrend)
                end do
                do iens=0,0
                    print *,'indx = ',iens,k
                    call printdatfile(6,indx(1,yrbeg,iens,k),npermax
     +                   ,nperyear,yrbeg,yrend)
                end do                
            end if
        endif
        if ( nens2.gt.max(imens(0),imens(k)) ) then
            write(0,*) 'filllinarray: error: nens2 > mens: ',nens2
     +           ,imens(0),imens(k)
            write(*,*) 'filllinarray: error: nens2 > mens: ',nens2
     +           ,imens(0),imens(k)
            call abort
        endif
*
        n = 0
        do iiens=nens1,nens2
            lastvalid = .false.
            if ( imens(0).gt.0 ) then
                iens = iiens
            else
                iens = 0
            endif
            if ( imens(k).gt.0 ) then
                jens = iiens
            else
                jens = 0
            endif
            do yr=yr1-1,yr2
                if ( j1.ne.j2 .and. j2-j1+1.ne.nperyear ) then
                    lastvalid = .false.
                endif
                do jj=j1,j2
                    if ( fix2 ) then
                        j = jj+lag
                    else
                        j = jj
                    endif
                    call normon(j,yr,i,nperyear)
                    if ( i.lt.yr1 .or.i.gt.yr2 ) goto 710
                    m = j-lag
                    call normon(m,i,ii,nperyear)
                    if ( ii.lt.yr1 .or.ii.gt.yr2 ) goto 710
                    if ( .false. ) then
                        print *,'data(',j,i,iens,')    = ',
     +                       data(j,i,iens)
                        print *,'indx(',m,ii,jens,k,') = ',
     +                       indx(m,ii,jens,k)
                    end if
                    if (  data(j,i,iens).lt.1e33 .and.
     +                    indx(m,ii,jens,k).lt.1e33 .and.
     +                    ( lconting .or. (
     +                    (data(j,i,iens).le.maxdata) .eqv.
     +                    (data(j,i,iens).ge.mindata) .eqv.
     +                    (maxdata.ge.mindata) ) .and. (
     +                    (indx(m,ii,jens,k).le.maxindx)
     +                    .eqv.
     +                    (indx(m,ii,jens,k).ge.minindx)
     +                    .eqv.(maxindx.ge.minindx) ) ) )
     +                    then
                        if ( n.eq.0 .and. lwrite ) then
                            print *,'filllinarray: first valid point '
                        endif
                        if ( lwrite ) then
                            print '(i3,i5,i3,g12.4,i3,i5,i3,g12.4,i3)',
     +                            j,i,iens,data(j,i,iens),
     +                            m,ii,jens,indx(j,i,iens,k),k
                        endif
                        n = n+1
                        if ( n.gt.ndata ) then
                            write(0,*) 'filllinarray: error: n>ndata '
     +                            ,ndata
                            write(*,*) 'filllinarray: error: n>ndata '
     +                            ,ndata
                            call abort
                        endif
                        ddata(n) = data(j,i,iens)
                        dindx(n) = indx(m,ii,jens,k)
                        lfirst(n) = .not.lastvalid
                        lastvalid = .true.
                        if ( yrstop.ne.-999 ) then
                            yrmo(1,n) = i
                            yrmo(2,n) = j
                            yrstart = min(yrstart,i,ii)
                            yrstop  = max(yrstop,i,ii)
                        endif
                        if ( lwrite .and. lfirst(n) ) print *
     +                       ,'boundary at ',n,j,i,iens
                        if ( nfittime.gt.0 ) then
                            dddata(n) = 0
                            do if=1,nfittime
                                jm = j-if
                                call normon(jm,i,im,nperyear)
                                jp = j+if
                                call normon(jp,i,ip,nperyear)
                                if (  ip.gt.yrend .or. im.lt.yrbeg )
     +                                then
                                    n = n - 1
                                    goto 710
                                endif
                                if (  data(jp,ip,iens).gt.1e33 .or.
     +                                data(jm,im,iens).gt.1e33 ) then
                                    n = n - 1
                                    goto 710
                                endif
                                dddata(n) = dddata(n) + filter(1+if)*(
     +                                data(jm,im,iens) - data(jp,ip,iens
     +                                ))
                            enddo ! ip=1,nfittime
                            indx(j,i,iens,indxmx) = dddata(n)
                        endif   ! nfittime
                    else
                        goto 710
                    endif       ! valid data point
                    goto 720
 710                continue
*                   invalid data point
                    lastvalid = .false.
 720                continue
                enddo           ! month jj
            enddo               ! year yr
        enddo                   ! iens
        end
