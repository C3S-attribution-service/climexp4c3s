        subroutine getruncorr(dindx,ddata,lfirst,dddata,ndata,j1,j2
     +        ,lag,k,month,nperyear,imens,indxmx,indx
     +        ,data,npermax,yrbeg,yrend,nensmax,n0
     +        ,filter,string,lboot,lprint,rmin,rmax,zdif)
*
*       compute running correlations
*
        implicit none
        integer ndata,j1,j2,lag,k,month,nperyear,indxmx,npermax
     +        ,yrbeg,yrend,nensmax,n0
        integer imens(indxmx)
        real dindx(ndata),ddata(ndata),dddata(ndata)
        logical lfirst(ndata)
        real data(npermax,yrbeg:yrend,0:nensmax),
     +        indx(npermax,yrbeg:yrend,0:nensmax,indxmx)
        real filter(100)
        logical lboot,lprint
        real rmin,rmax,zdif
        character*(*) string
#include "getopts.inc"
        integer yr1s,yr2s,n,yr,m
        real result,dresult(-2:2),prob,zmin,zmax,z
        logical lblank,llwrite
*
        lblank = lprint
        if ( lprint ) then
            write(14,'(a)') '#year  #pts  prob 2.5% 16% 50% 84% 97.5%'
            write(14,'(2a)') '# ',string
        endif
        yr1s = yr1
        yr2s = yr2
        rmin = +3e33
        rmax = -3e33
        zmin = +3e33
        zmax = -3e33
        do yr=yr1s,yr2s
            yr1 = yr - nyrwindow/2
            yr2 = yr1 + nyrwindow - 1
            n = 0
            llwrite = lwrite
            lwrite = .false.
            call filllinarray(dindx,ddata,lfirst,dddata,ndata,n,j1,j2
     +           ,lag,k,nperyear,imens,indxmx,indx,data,npermax,yrbeg
     +           ,yrend,nensmax,filter)
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
            if ( n.lt.m ) then
                if ( lwrite ) print *
     +               ,'getruncorr: not enough points:',n,m
                if ( lblank ) write(14,'(a)')
                lblank = .FALSE.
                goto 700
            elseif ( lwrite ) then
                print *,'getruncorr: enough points:',n,max(m,minnum)
            endif
            lblank = lprint
            call printcorr(dindx,ddata,lfirst,dddata,n,n0,j1,j2,month
     +           ,nperyear,lag,string,lboot,.false.,result,dresult,prob)
            if ( lprint .and. result.lt.1e33 ) then
                write(14,'(i4,i6,99g11.2)') yr,n,result,prob,dresult
            endif
            if ( lwrite ) then
                write(*,'(i4,i6,99g11.2)') yr,n,result,prob,dresult
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
        end
