        program plotdatrel
*
*       compute relative anomalies
*
        implicit none
#include "param.inc"
        integer i,j,k,l,n,year,loop,type,mean,nperyear,yr1,yr2,nens1
     +        ,nens2,iens,nperyearnew
        real data(npermax,yrbeg:yrend,0:nensmax),s(2),nextx,lasty
        logical lwrite,lanomal,lskip,lexist,ensemble
        character line*255,ensfile*255,var*20,units*20
	integer iargc,llen
        lwrite = .FALSE.
	lskip = .FALSE.
        call getenv('LWRITE_PLOTDAT',line)
        if ( line(1:1).eq.'t' .or. line(1:1).eq.'T' ) then
            lwrite = .true.
        endif
*
*       init
*
        if ( iargc().lt.2 ) then
            write(0,*)
     +           'usage: plotdatrel nperyear [anom [yr1 yr2]] datafile'
            stop
        endif
*
        call getarg(1,line)
        read(line,*,err=900) nperyearnew
        if ( iargc().gt.2 ) then
            call getarg(2,line)
            if ( line(1:4).ne.'anom' ) then
                write(0,*) 'plotdat: error: expecting ''anom'''
                call abort
            endif
            lanomal = .TRUE.
            if ( iargc().eq.3 ) then
                yr1 = yrbeg
                yr2 = yrend
                call getarg(3,line)
            elseif ( iargc().eq.4 ) then
                call getarg(3,line)
                read(line,*,err=901) yr1
                call getarg(4,line)
                read(line,*,err=901) yr2
                call getarg(5,line)
            else
                write(0,*) 'plotdat: error: expecting 5 args, not '
     +                ,iargc()
                call abort
            endif
        else
            lanomal = .FALSE.
            call getarg(2,line)
        endif
        nens1 = 0
        if ( index(line,'%%').eq.0 .and. index(line,'++').eq.0 ) then
            ensemble = .FALSE.
            nens2 = 0
        else
            ensemble = .TRUE.
            do nens2=0,nensmax
                ensfile=line
                call filloutens(ensfile,nens2)
                if ( lwrite ) write(0,*) 'looking for file '
     +                ,trim(ensfile)
                inquire(file=ensfile,exist=lexist)
                if ( .not.lexist ) then
                    if ( nens2.eq.0 ) then
                        nens1=1
                        cycle
                    else
                        exit
                    endif
                endif
            enddo
            nens2 = nens2 - 1
        endif
        if ( lwrite ) write(0,*) 'nens1,nens2 = ',nens1,nens2
        do iens=nens1,nens2
            ensfile = line
            if ( ensemble ) call filloutens(ensfile,iens)
            if ( lwrite ) write(0,*) 'reading file '
     +            ,ensfile(1:llen(ensfile))
            call readseries(ensfile,data(1,yrbeg,iens),npermax,yrbeg
     +           ,yrend,nperyear,var,units,.false.,lwrite)
            if ( lanomal ) call anomalrel(data(1,yrbeg,iens),npermax
     +            ,nperyear,yrbeg,yrend,yr1,yr2,nperyearnew)
        enddo
*
*       print
*       
        do iens=nens1,nens2
            if ( ensemble ) then
                print '(a,i4)','# ensemble member ',iens
                print '(a)'
            endif
            do i=yrbeg,yrend
                do j=1,nperyear
                    k = j
                    n = nperyear
                    if ( nperyear.eq.366 .and. (mod(i,4).ne.0 .or. (
     +                    mod(i,100).eq.0 .and. mod(i,400).ne.0)) ) then
                        n = 365
                        if ( j.gt.60 ) k = j-1
                        if ( j.eq.60 .and. data(j,i,iens).lt.1e33 ) then
                            print '(a,i4,f14.6)'
     +                            ,'# skipping valid data on 29feb',i,
     +                            data(j,i,iens)
                        endif
                    endif
                    if ( abs(data(j,i,iens)).lt.1e33 ) then
                        print '(f10.4,g14.6)',i+(k-1.)/n,data(j,i,iens)
                        nextx = i+real(k)/n
                        lasty = data(j,i,iens)
                        lskip = .TRUE.
                    elseif ( lskip .and. .not.(n.eq.366 .and. j.eq.60) )
     +                        then
*                       to counter peculiarity in gnuplot steps plotting
                        print '(f10.4,g14.6,a)',nextx,lasty,
     +                        '# repeat last y to get nice gnuplot plot'
                        print '(a)'
                        lskip = .FALSE.
                    endif
                enddo
            enddo
        enddo
*
        goto 999
 900    write(0,*) 'plotdat: expecting nperyear, found ',trim(line)
        call abort
 901    write(0,*) 'plotdat: expecting year, found ',trim(line)
        call abort
 999    continue
        end
