        program plotdat
*
*       convert 14-column data into a format suitable for gnuplot,
*       taking anomalies if requested
*
        implicit none
#include "param.inc"
#include "getopts.inc"
        integer i,j,k,l,n,year,loop,type,nperyear,iens,iarg,mens
     +       ,mens1,ayr1,ayr2
        real data(npermax,yrbeg:yrend,0:nensmax),mean(npermax)
     +       ,s,nextx,lasty
        logical lanomal,lskip,lexist
        character line*255,ensfile*255,var*40,units*40,title*1000
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
        if ( iargc().eq.0 ) then
            print *,'usage: plotdat [anom [yr1 yr2]] '//
     +           '[normsd [ave N]] datafile'
            stop
        endif
*
        if ( iargc().gt.1 ) then
            call getarg(1,line)
            if ( line(1:4).eq.'anom' ) then
                lanomal = .TRUE.
                if ( iargc().lt.3 ) then
                    ayr1 = yrbeg
                    ayr2 = yrend
                    iarg = 2
                else
                    call getarg(2,line)
                    read(line,*,err=901) ayr1
                    ayr1 = max(yrbeg,ayr1)
                    call getarg(3,line)
                    read(line,*,err=901) ayr2
                    ayr2 = min(yrend,ayr2)
                    iarg = 4
                endif
            else
                lanomal = .FALSE.
                iarg = 1
            end if
        endif
  
        title = '#'
        do i=0,iargc()
            call getarg(i,line)
            title(len_trim(title)+2:) = line
        end do
        print '(a)',trim(title)
        call getarg(iargc(),line)
        lstandardunits = .false. ! because we cannot call getopts before we know mens1,mens...
        call readensseries(line,data,npermax,yrbeg,yrend,nensmax
     +        ,nperyear,mens1,mens,var,units,lstandardunits,lwrite)
        call getopts(iarg,iargc()-1,nperyear,yrbeg,yrend,.true.,mens1
     +       ,mens)
        if ( lnormsd ) print '(a)','# taking relative anomalies'
        if ( lstandardunits ) then
            do iens=mens1,mens
                call makestandardseries(data(1,yrbeg,iens),
     +               npermax,yrbeg,yrend,nperyear,var,units,lwrite)
            end do
        endif

        if ( lanomal ) then
            if ( lensanom ) then
                call ensanomalclim(data,npermax,nperyear,yrbeg,yrend,
     +               nens1,nens2,ayr1,ayr2,mean)
                if ( lnormsd ) then
                    call takerelanom(data,mean,npermax,yrbeg,yrend,
     +                   nens1,nens2,nperyear,lsum,lwrite)
                end if
            else
                do iens=nens1,nens2
                    call ensanomalclim(data,npermax,nperyear,yrbeg,yrend
     +                   ,iens,iens,ayr1,ayr2,mean)
                    if ( lnormsd ) then
                        call takerelanom(data,mean,npermax,yrbeg,yrend,
     +                       iens,iens,nperyear,lsum,lwrite)
                    end if
                end do
            end if
        end if                  ! lanomal
*
*       print
*       
        do iens=nens1,nens2
            if ( nens2.gt.nens1 ) then
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
  901   write(0,*) 'plotdat: expecting year, found ',line
        call abort
  999   continue
        end

