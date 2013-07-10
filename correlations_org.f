        program correlate
*
*       correlate station parameters from www.ncdc.noaa.gov or other
*       sources with all kinds of indices.  Currently included are:
*           SOI SLP index, from 1866 to 1997, from Phil Jones via Marleen
*           4 NINO SST indices, from 1950 to present, from 
*               http://nic.fb4.noaa.gov:80/data/cddb/
*           NAO SLP index, from 1865 to 1996
*           SIDC sunspot index, from Zuerich,  http://www.astroinfo.ch/sunspot
*           SIDC sunspot cycle length, from Zuerich
*
*       GJvO hack, nov-1997, revised dec-1997, jan-1998
*
        implicit none
        integer yrbeg,yrend,indxmx
        parameter (yrbeg=1600,yrend=2020,indxmx=8)
        integer i,ii,j,k,l,m,n,j1,j2,lag,lsum,year,month,
     +        iorder(12,yrbeg:yrend)
        logical lwrite,logscale,txtdat,dump,lincl(indxmx),anom,lrank
        real data(12,yrbeg:yrend),indx(12,yrbeg:yrend,indxmx),
     +      anino(8),val12(12),adata,aindx,sxx,syy,sxy,a,yrmin(25)
     +      ,yrmax(25),xmin,xmax,rise,fall,slength(25),yr,r,t,df,prob
     +      ,absent
        parameter (absent=3e33)
        character line*128,strindx(indxmx)*8,oper*1
        real betai
        integer iargc
        data strindx /
     +      'SOI     ',
     +      'NINO12  ','NINO3   ','NINO4   ','NINO3.4 ',
     +      'NAO     ',
     +      'sunspot ','  length'/
        data lincl /8*.FALSE./
        lwrite = .FALSE.
*
*       check arguments
*
        n = iargc()
        if ( n.lt.1 ) then
            print *,'usage: correlate datafile [lag n] [sum|max|min m] '
     +            //' [log|rank]'
            stop
        endif
        logscale = .FALSE.
        lrank = .FALSE.
        lag  = 0
        lsum = 1
        dump = .FALSE.
        anom = .FALSE.
        do i=2,n
            call getarg(i,line)
            if ( line(1:3).eq.'log' ) then
                logscale = .TRUE.
                print *,'Set logscale on'
            endif
            if ( line(1:4).eq.'rank' ) then
                lrank = .TRUE.
                print *,'Set rank correllations on'
            endif
            if ( line(1:3).eq.'lag' ) then
                call getarg(i+1,line)
                read(line,*,err=905) lag
                print *,'Using lag ',lag
            endif
            if ( line(1:3).eq.'sum' ) then
                call getarg(i+1,line)
                read(line,*,err=906) lsum
                oper = '+'
            endif
            if ( line(1:3).eq.'max' ) then
                call getarg(i+1,line)
                read(line,*,err=906) lsum
                oper = 'a'
            endif
            if ( line(1:3).eq.'min' ) then
                call getarg(i+1,line)
                read(line,*,err=906) lsum
                oper = 'i'
            endif
            if ( line(1:4).eq.'dump' ) then
                dump = .TRUE.
                print *,'Dumping on dump.dat'
            endif
            if ( line(1:4).eq.'anom' ) then
                anom = .TRUE.
                print '(a)','# taking anomalies wrt climatology data'
            endif
            if ( line(1:3).eq.'soi' ) then
                lincl(1) = .TRUE.
            endif
            if ( line(1:6).eq.'nino12' ) then
                lincl(2) = .TRUE.
            endif
            if ( line(1:6).eq.'nino3 ' ) then
                lincl(3) = .TRUE.
            endif
            if ( line(1:6).eq.'nino4 ' ) then
                lincl(4) = .TRUE.
            endif
            if ( line(1:7).eq.'nino3.4' ) then
                lincl(5) = .TRUE.
            endif
            if ( line(1:3).eq.'nao' ) then
                lincl(6) = .TRUE.
            endif
            if ( line(1:7).eq.'sunspot' ) then
                lincl(7) = .TRUE.
            endif
            if ( line(1:9).eq.'sunlength' ) then
                lincl(8) = .TRUE.
            endif
        enddo
        do i=1,indxmx
            if ( lincl(i) ) goto 100
        enddo
*       default set-up: SOI, NINO3, NAO
        lincl(1) = .TRUE.
        lincl(3) = .TRUE.
        lincl(6) = .TRUE.
  100   continue
        call getarg(1,line)
*
*       init
*
        if ( dump ) then
            open(10,file='dump.dat',status='unknown')
            write(10,'(2a)') '# ',line(1:index(line,' '))
        endif
        call makeabsent(data,12,yrbeg,yrend)
        do k=1,indxmx
            call makeabsent(indx(1,yrbeg,k),12,yrbeg,yrend)
        enddo
*
*       read data from station file downloaded from
*       http://www.ncdc.noaa.gov/ghcn/ghcnV1.CLIMVIS.html
*
        call rdtxtdat(line,data,yrbeg,yrend)
        if ( logscale ) then
            do i=yrbeg,yrend
                do j=1,12
                    if ( data(j,i).lt.0.9*absent ) then
                        if ( data(j,i).eq.0 ) data(j,i) = 0.1
                        data(j,i) = log(data(j,i))
                    endif
                enddo
            enddo
        endif
*
*       get the SOI data
*
        if ( lincl(1) ) then
            call readdat(indx(1,yrbeg,1),yrbeg,yrend,'soi.knmi')
        endif
*
*       get old NINO indices from Alexey Kaplan
*
        if ( lincl(2).or.lincl(3).or.lincl(4).or.lincl(5) ) then
            if ( lwrite ) print *,'Opening file sstkap.indices'
            open(1,file='sstkap.indices',status='old')
*           skip header
            do i=1,1
                read(1,'(a)') line
            enddo
*           data is in the format 
*           yyyy mm NINO12 dNINO12 NINO3 dNINO3 NINO4 dNINO4 NINO3.4 dNINO3.4
            do i=1,10000
                read(1,*,err=902,end=150) year,month,anino
                if ( year.lt.200 ) year = year + 1900
                do j=1,4
                    indx(month,year,j+1) = anino(2*j)
                    if ( indx(month,year,j+1).gt.10 .or.
     +                   indx(month,year,j+1).lt.-10 ) then
                        indx(month,year,j+1) = absent
                    endif
                enddo
                if ( lwrite ) print *,year,month,(anino(2*j),j=1,4)
            enddo
  150       continue
            close(1)
        endif
*
*       and the NINO indices from http://nic.fb4.noaa.gov:80/data/cddb/
*
        if ( lincl(2).or.lincl(3).or.lincl(4).or.lincl(5) ) then
            if ( lwrite ) print *,'Opening file sstoi.indices'
            open(1,file='sstoi.indices',status='old')
*           skip header
            do i=1,1
                read(1,'(a)') line
            enddo
*           data is in the format 
*           yy mm NINO12 dNINO12 NINO3 dNINO3 NINO4 dNINO4 NINO3.4 dNINO3.4
            do i=1,10000
                read(1,*,err=902,end=200) year,month,anino
                if ( year.lt.200 ) year = year + 1900
                do j=1,4
                    indx(month,year,j+1) = anino(2*j)
                    if ( indx(month,year,j+1).gt.10 .or.
     +                   indx(month,year,j+1).lt.-10 ) then
                        indx(month,year,j+1) = absent
                    endif
                enddo
                if ( lwrite ) print *,year,month,(anino(2*j),j=1,4)
            enddo
  200       continue
            close(1)
        endif
*
*       and the NAO data
*
        if ( lincl(6) ) then
            call readdat(indx(1,yrbeg,6),yrbeg,yrend,'nao')
        endif
*
*       and the sunspot data, 1749-1991
*
        if ( lincl(7) ) then
            if ( lwrite ) print *,'Opening file sunhist_1.html'
            open(1,file='sunhist_1.html',status='old')
*           data is in the format
*           yyyy indx(1) ... indx(12)
            do i=1,10000
                read(1,'(a)',err=904,end=400) line
                if ( line(1:1).eq.'1' ) then
  390               continue
                    j = index(line,'*')
                    if ( j.ne.0 ) then
*                   get rid of footnotes
                        line(j:j) = ' '
                        goto 390
                    endif
                    read(line,*,err=904) year,val12
                    do j=1,12
                        if ( val12(j).gt.500 .or. val12(j).lt.0 ) then
                            val12(j) = absent
                        endif
                        indx(j,year,7) = val12(j)
                        if ( lwrite ) print *,year,j,val12(j)
                    enddo
                endif
            enddo
  400       continue
            close(1)
        endif
        if ( lincl(8) ) then
            if ( lwrite ) print *,'Opening file sunhist_2.html'
            open(1,file='sunhist_2.html',status='old')
*           data is in the format
*           no yrmin valmin yrmax valmax rise fall length
            do i=1,10000
                read(1,'(a)',err=904,end=500) line
                if ( line(1:2).eq.'  ' .and.
     +                (line(3:3).ge.'0' .and. line(3:3).le.'9' .or. 
     +                 line(3:3).eq.' ') .and.
     +                (line(4:4).ge.'0' .and. line(4:4).le.'9') ) then
                    read(line,*,err=904) j
                    if ( j.lt.22 ) then
                        read(line,*,err=904) j,yrmin(j),xmin,yrmax(j)
     +                        ,xmax,rise,fall,slength(j)
                    else
                        read(line,*,err=904) j,yrmin(j),xmin,yrmax(j)
     +                        ,xmax
                    endif
                    if ( lwrite ) print *,'read ',j,yrmin(j),xmin
     +                    ,yrmax(j),xmax,rise,fall,slength(j)
                endif
            enddo
  500       continue
            j = 1
            do year=yrbeg,yrend
                do month=1,12
                    yr = year + (month-0.5)/12
                    if ( j.le.22 ) then
                        if ( yr.gt.yrmin(j) ) then
                            j = j + 1
                        endif
                        if ( j.gt.1 .and. j.le.22 ) then
                            indx(month,year,8) = slength(j-1)
                            if ( lwrite ) print *,year,month,slength(j-1
     +                            )
                        endif
                    endif
                enddo
            enddo
            close(1)
        endif
*
*       anomalies
*
        if ( anom ) then
            call anomal(data,yrbeg,yrend)
        endif
*
*       sum
*
        if ( lsum.gt.1 ) then
            call sumit(data,12,yrbeg,yrend,lsum,oper)
            do k=1,indxmx
                call sumit(indx(1,yrbeg,k),12,yrbeg,yrend,lsum,'+')
            enddo
        endif
*
*       correlate!
*
        print '(a)',
     +      '========================================================='
        print '(a)',
     +      ' index    corr  sign.  no        data           index' 
        print '(a)',
     +      '  name             % pnts    mean    s.d.    mean    s.d.'
        print '(a)',
     +      '========================================================='
        do month=0,12
            if ( month.eq.0 ) then
                j1 = 1
                j2 = 12
                print *,'All year:'
                if ( dump ) write(10,'(a)') '# All year'
            else
                j1 = month
                j2 = month
                print '(a,12i3)','Month: ',(1+mod(month+l-1,12),
     +                l=0,lsum-1)
                if ( dump ) write(10,'(a,12i3)') '# Month: ',(1
     +              +mod(month+l-1,12),l=0,lsum-1)
            endif
            do k=1,indxmx
                if ( .not.lincl(k) ) goto 800
                n = 0
                adata = 0
                aindx = 0
                if ( dump ) write(10,'(2a)') '# ',strindx(k)
                do i=yrbeg,yrend
                    do j=j1,j2
                        m = j-lag
                        call normon(m,i,ii)
                        if ( ii.lt.yrbeg .or.ii.gt.yrend ) goto 710
                        if ( data(j,i).lt.1e33 .and. 
     +                       indx(m,ii,k).lt.1e33 ) then
                            n = n+1
                            adata = adata + data(j,i)
                            aindx = aindx + indx(m,ii,k)
                        endif
  710                   continue
                    enddo
                enddo
                if ( n.gt.0 ) then
                    adata = adata/n
                    aindx = aindx/n
                    sxx = 0
                    syy = 0
                    sxy = 0
                    do i=yrbeg,yrend
                        do j=j1,j2
                            m = j-lag
                            call normon(m,i,ii)
                            if ( ii.lt.yrbeg .or.ii.gt.yrend ) goto 720
                            if ( data(j,i).lt.1e33 .and. 
     +                          indx(m,ii,k).lt.1e33 ) then
                                a = data(j,i)
                                sxx = sxx + (a - adata)**2
                                syy = syy + (indx(m,ii,k) - aindx)**2
                                sxy = sxy + (a - adata)*(indx(m,ii,k)
     +                              -aindx)
                                if ( dump ) write(10,'(2f12.4,i5,i3)')
     +                                indx(m,ii,k)/lsum,a,i,j
  720                           continue
                            endif
                        enddo
                    enddo
                    r = sxy/sqrt(sxx*syy)
                    df = real(n-2)
                    if ( 1-r**2.gt.0 ) then
                        t = r*sqrt(df/(1-r**2))
                    else
                        t = 10000
                    endif
                    prob = betai(0.5*df,0.5,df/(df+t**2))
                    if ( sxx.eq.0 ) sxx = 999.99
                    if ( syy.eq.0 ) syy = 999.99
                    print 1000,strindx(k),r,100*(1-prob),n,adata
     +                  ,sqrt(sxx/n),aindx,sqrt(syy/n)
                    if ( dump ) write(10,1000)'# '//strindx(k),r,100*(1
     +                  -prob),n,adata,sqrt(sxx/n),aindx,sqrt(syy/n) 
                endif
  800           continue
            enddo
        enddo
*       
 1000   format(a,f6.2,f6.1,i5,4f8.2)
*
*       error messages
*
        stop
  902   print *,'error reading NINO data'
        stop
  904   print *,'error reading sunspot data'
        stop
  905   print *,'error reading lag value ',line
        stop
  906   print *,'error reading sum value ',line
        end
