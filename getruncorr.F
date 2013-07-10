        subroutine getruncorr(dindx,ddata,lfirst,dddata,ndata,j1,j2
     +        ,lag,k,month,nperyear,imens,indxmx,indx
     +        ,data,npermax,yrbeg,yrend,nensmax,n0
     +        ,filter,string,lboot,lprint,rmin,rmax,zdif)
*
*       compute running correlations
*
        implicit none
        integer lperyear,lyr,lensmax
        parameter(lperyear=366,lyr=200,lensmax=31)
        integer ndata,j1,j2,lag,k,month,nperyear,indxmx,npermax
     +        ,yrbeg,yrend,nensmax,n0
        integer imens(0:indxmx)
        real dindx(ndata),ddata(ndata),dddata(ndata)
        logical lfirst(ndata)
        real data(npermax,yrbeg:yrend,0:nensmax),
     +        indx(npermax,yrbeg:yrend,0:nensmax,indxmx)
        real filter(100)
        logical lboot,lprint
        real rmin,rmax,zdif,savedata(lperyear,lyr,0:lensmax),
     +       saveindx(lperyear,lyr,0:lensmax),sig(1),a,b,siga,sigb,chi2,
     +       q
        character*(*) string
#include "getopts.inc"
        integer yr1s,yr2s,yr1d,yr2d,yr1i,yr2i,n,yr,yrp,m,iiens,iens,jens
     +       ,i
        integer yrmo(2,1)       ! dummy argument, not used
        real result,dresult(-2:2),prob,zmin,zmax,z
        logical lblank,llwrite
        integer init,iu
        save init
        data init /0/
*
        lblank = lprint
        if ( lprint ) then
            write(14,'(2a)') '# ',string
            write(14,'(a)') '#year  #pts   corr        prob '//
     +           '       2.5%         16%         50%'//
     +           '         84%       97.5%'
        endif
        if ( .false. ) then
            call rsunit(iu)
            write(string,'(a,i3.3,a)') 'check',init,'.txt'
            open(iu,file=string)
            init = init + 1
            write(iu,'(a,i3)') 'at beginning of getruncorr',k
            do iiens=nens1,nens2
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
                do yr=yr1,yr2
                    do m=1,nperyear
                        if ( data(m,yr,iens).lt.3e33 .or. 
     +                       indx(m,yr,jens,k).lt.3e33 ) then
                            write(iu,'(i3,i5,i3,2g20.12)') iiens,yr,m,
     +                           data(m,yr,iens),indx(m,yr,jens,k)
                        endif
                    enddo
                enddo
                close(iu)
            enddo
        endif
        yr1s = yr1
        yr2s = yr2
        rmin = +3e33
        rmax = -3e33
        zmin = +3e33
        zmax = -3e33
        do yr=yr1s,yr2s
            yr1 = max(yrbeg,yr - nyrwindow/2 + 1)
            yr2 = min(yrend,yr1 + nyrwindow - 1)
            if ( ldetrend ) then
*               rescue data
                if ( nperyear.gt.lperyear .or. yr2-yr1+1.gt.lyr ) then
                    write(0,*) 'getruncorr: error: fixed array '//
     +                   'too small; ',nperyear,lperyear,yr2-yr1+1,lyr
                    call abort
                endif
                if ( lag.eq.0 ) then
                    yr1i = yr1
                    yr2i = yr2
                    yr1d = yr1
                    yr2d = yr2
                elseif ( fix2 ) then
                    yr1i = yr1
                    yr2i = yr2
                    if ( lag.gt.0 ) then
                        yr1d = yr1d + (lag-1)/nperyear + 1
                        yr2d = yr2d + (lag-1)/nperyear + 1
                    else
                        yr1d = yr1d + (lag+1)/nperyear - 1
                        yr2d = yr2d + (lag+1)/nperyear - 1
                    endif
                else
                    yr1d = yr1
                    yr2d = yr2
                    if ( lag.gt.0 ) then
                        yr1i = yr1i - (lag-1)/nperyear - 1
                        yr2i = yr2i - (lag-1)/nperyear - 1
                    else
                        yr1i = yr1i - (lag+1)/nperyear + 1
                        yr2i = yr2i - (lag+1)/nperyear + 1
                    endif
                endif
                if ( max(yr2d-yr1d+1,yr2i-yr1i+1).gt.lyr ) then
                    write(0,*) 'getruncorr: error: not enough room',
     +                   ' to rescue detrended data, increase lyr to ',
     +                   max(yr2d-yr1d+1,yr2i-yr1i+1)
                    call abort
                endif
                do iiens=nens1,nens2
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
                    if ( max(iens,jens).gt.lensmax ) then
                        write(0,*) 'getruncorr: error: increase lensmax'
     +                       ,iens,jens,lensmax
                        call abort
                    endif
                    do yrp=yr1d,yr2d
                        do m=1,nperyear
                            savedata(m,yrp-yr1d+1,iens) = 
     +                           data(m,yrp,iens)
                        enddo
                    enddo
                    do yrp=yr1i,yr2i
                        do m=1,nperyear
                            saveindx(m,yrp-yr1i+1,jens) = 
     +                           indx(m,yrp,jens,k)
                        enddo
                    enddo
                enddo
                call detrend(data(1,yr1d,iens),npermax,nperyear,yr1d
     +               ,yr2d,yr1d,yr2d,m1,m2,lsel)
                call detrend(indx(1,yr1i,jens,k),npermax,nperyear,
     +               yr1i,yr2i,yr1i,yr2i,m1,m2,lsel)
            endif
            n = 0
            llwrite = lwrite
            lwrite = .false.
            call filllinarray(dindx,ddata,lfirst,dddata,ndata,n,j1,j2
     +           ,lag,k,nperyear,imens,indxmx,indx,data,npermax,yrbeg
     +           ,yrend,nensmax,filter,-999,-999,yrmo)
            if ( ldetrend ) then
*               restore data
                do iiens=nens1,nens2
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
                    do yrp=yr1d,yr2d
                        do m=1,nperyear
                            data(m,yrp,iens) =
     +                           savedata(m,yrp-yr1d+1,iens)
                        enddo
                    enddo
                    do yrp=yr1i,yr2i
                        do m=1,nperyear
                            indx(m,yrp,jens,k) = 
     +                           saveindx(m,yrp-yr1i+1,jens)
                        enddo
                    enddo
                enddo
            endif
            lwrite = llwrite
            if ( minnum.le.0 ) then
                if ( month.eq.0 ) then
                    m = 12*(nyrwindow-1)
                else
                    m = (nyrwindow-1)
                endif
            else
                m = minnum
            endif
            if ( n.lt.(1+nens2-nens1)*m ) then
                if ( lwrite ) print *,'getruncorr: not enough points:',
     +               n,(1+nens2-nens1)*m
                if ( lblank ) write(14,'(a)')
                lblank = .FALSE.
                goto 700
            elseif ( lwrite ) then
                print *,'getruncorr: enough points:',n,max(m,minnum)
            endif
            lblank = lprint
            if ( irunvar.eq.1 ) then
*               running correlatons
                call printcorr(dindx,ddata,lfirst,dddata,yrmo,n,n0,j1,j2
     +               ,month,nperyear,lag,string,lboot,.false.,result
     +               ,dresult,prob)
                if ( lprint .and. result.lt.1e33 ) then
                    write(14,'(i4,i6,99g12.4)') yr,n,result,prob,dresult
                endif
                if ( lwrite ) then
                    write(*,'(i4,i6,99g12.4)') yr,n,result,prob,dresult
                endif
                if ( abs(result).lt.1e30 ) then
                    rmin = min(rmin,result)
                    rmax = max(rmax,result)
                    if ( abs(result).lt.1 ) then
                        z = (1+result)/(1-result)
                        z = 0.5*log(z)
                    elseif ( result.gt.0 ) then
                        z = +1e10
                    else
                        z = -1e10
                    endif
                    zmin = min(zmin,z)
                    zmax = max(zmax,z)
                endif
            elseif ( irunvar.eq.2 ) then
*               running regressions, assume normality for the time being
                call fit(dindx,ddata,n,sig,0,a,b,siga,sigb,chi2,q)
                if ( lprint .and. b.lt.1e33 ) then
                    write(14,'(i4,i6,99g12.4)') yr,n,b,q,b-2*sigb,b-sigb
     +                   ,b,b+sigb,b+2*sigb
                endif
                if ( lwrite ) then
                    do i=1,n
                        print '(i4,2g16.8)',i,dindx(i),ddata(i)
                    end do
                    write(*,'(i4,i6,99g12.4)') yr,n,b,q,b-2*sigb,b-sigb
     +                   ,b,b+sigb,b+2*sigb
                endif
                rmin = min(rmin,b)
                rmax = max(rmax,b)
                zmin = rmin
                zmax = rmax
            else
                write(0,*) 'getruncorr: error: irunvar = ',irunvar
                write(*,*) 'getruncorr: error: irunvar = ',irunvar
                call abort
            endif
            if ( lprint ) call keepalive(yr-yr1s,yr2s-yr1s)
  700       continue
        enddo
        if ( rmax.lt.-1e33 ) then
            rmax = 3e33
            zdif = 3e33
        else
            zdif = zmax - zmin
        endif
        if ( lprint ) write(14,'(a)')
        yr1 = yr1s
        yr2 = yr2s
        if ( .false. ) then
            call rsunit(iu)
            write(string,'(a,i3.3,a)') 'check',init,'.txt'
            open(iu,file=string)
            init = init + 1
            write(iu,'(a,i3)') 'at end of getruncorr',k
            do iiens=nens1,nens2
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
                do yr=yr1,yr2
                    do m=1,nperyear
                        if ( data(m,yr,iens).lt.3e33 .or. 
     +                       indx(m,yr,jens,k).lt.3e33 ) then
                            write(iu,'(i3,i5,i3,2g20.12)') iiens,yr,m,
     +                           data(m,yr,iens),indx(m,yr,jens,k)
                        endif
                    enddo
                enddo
            enddo
            close(iu)
        endif
        end
