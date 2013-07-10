        program sstoi2dat
*
*       convert the file sstoi.indices to files nino*.dat
*       readable by the correlation programs.
*
        implicit none
        integer yrbeg,yrend
        parameter (yrbeg=1856,yrend=2020)
        integer i,j,k,year,month,ii(8)
        real anino(8),indx(12,yrbeg:yrend,2:5)
        character line*80,file*9
*
*       init
*
        do k=2,5
            do i=yrbeg,yrend
                do j=1,12
                    indx(j,i,k) = -999.9
                enddo
            enddo
        enddo
*
*       old data
*
        if ( .true. ) then
        open(1,file='sstkap.indices',status='old')
*       skip header
        do i=1,1
            read(1,'(a)') line
        enddo
*       data is in the format 
*       yy mm NINO12 dNINO12 NINO3 dNINO3 NINO4 dNINO4 NINO3.4 dNINO3.4
        do i=1,10000
            read(1,*,err=902,end=100) year,month,anino
            if ( year.lt.200 ) year = year + 1900
            do j=1,4
                if ( indx(month,year,j+1).gt.10 .or.
     +               indx(month,year,j+1).lt.-10 ) then
                    indx(month,year,j+1) = -999.9
                endif
                indx(month,year,j+1) = anino(2*j)
            enddo
        enddo
  100   continue
        close(1)
        end if
*
*       new data
*
        open(1,file='sstoi.indices',status='old')
*       skip header
        do i=1,1
            read(1,'(a)') line
        enddo
*       data is in the format 
*       yy mm NINO12 dNINO12 NINO3 dNINO3 NINO4 dNINO4 NINO3.4 dNINO3.4
        do i=1,10000
            read(1,*,err=902,end=200) year,month,anino
            if ( year.lt.200 ) year = year + 1900
***            if ( year.ge.1975 ) then
                do j=1,4
                    if ( indx(month,year,j+1).gt.10 .or.
     +                    indx(month,year,j+1).lt.-10 ) then
                        indx(month,year,j+1) = -999.9
                    endif
                    indx(month,year,j+1) = anino(2*j)
                enddo
***            endif
        enddo
  200   continue
        close(1)
*
*       write data
*
        do j=2,5
            write(file,'(a,i1,a)') 'nino',j,'.dat'
            open(1,file=file,status='unknown')
            if ( j.eq.2 ) then
                write(1,'(2a)') '# NINO12 index, SST anomalies in ',
     +               '10S-EQ, 80W-90W'
            elseif ( j.eq.3 ) then
                write(1,'(2a)') '# NINO3 index, SST anomalies in ',
     +               '5S-5N, 90W-150W'
            elseif ( j.eq.4 ) then
                write(1,'(2a)') '# NINO4 index, SST anomalies in ',
     +               '5S-5N, 160E-150W'
            elseif ( j.eq.5 ) then
                write(1,'(2a)') '# NINO3.4 index, SST anomalies in ',
     +               '5S-5N, 120W-170W'
            endif
            write(1,'(4a)') '# 1856-1949: <a href=',
     +           'http://ingrid.ldgo.columbia.edu/SOURCES/',
     +           '.Indices/.nino/.KAPLAN/',
     +           '>Kaplan reconstruction</a>'
            write(1,'(4a)') '# 1950-now: <a href=',
     +           'http://www.cpc.noaa.gov/data/indices/',
     +           '>CPC (Reynolds OI SST)</a>'
            call date_and_time(values=ii)
            write(1,'(a,i4,a,i2.2,a,i2.2)') '# last updated: ',ii(1),'-'
     +           ,ii(2),'-',ii(3)
            write(1,'(a)') '# SSTA [C]'
            do i=yrbeg,yrend
                do k=1,12
                    if ( indx(k,i,j).gt.-900 ) goto 210
                enddo
                goto 220
  210           continue
                write(1,'(i5,12f8.2)') i,(indx(k,i,j),k=1,12)
  220           continue
            enddo
            close(1)
        enddo
        stop
*
*       error messages
*
  902   print *,'error reading NINO data'
        stop
        end
