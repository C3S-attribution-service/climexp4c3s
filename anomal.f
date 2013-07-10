        subroutine anomal(data,npermax,nperyear,yrbeg,yrend,yr1,yr2)
        implicit none
        integer npermax,nperyear,yrbeg,yrend,yr1,yr2
        real data(npermax,yrbeg:yrend)
        call ensanomal(data,npermax,nperyear,yrbeg,yrend,0,0,yr1,yr2)
        end

        subroutine ensanomal(data,npermax,nperyear,yrbeg,yrend,nens1
     +       ,nens2,yr1,yr2)
*
*       subtract the climatology determined from the data itself over yr1:yr2
*
        implicit none
        integer npermax,nperyear,yrbeg,yrend,yr1,yr2,nens1,nens2
        real data(npermax,yrbeg:yrend,0:nens2)
        real,allocatable :: mean(:)
        allocate(mean(nperyear))
        call ensanomalclim(data,npermax,nperyear,yrbeg,yrend,nens1
     +       ,nens2,yr1,yr2,mean)
        deallocate(mean)
        end

        subroutine anomalclim(data,npermax,nperyear,yrbeg,yrend,yr1,yr2
     +       ,mean)
        integer npermax,nperyear,yrbeg,yrend,yr1,yr2
        real data(npermax,yrbeg:yrend),mean(nperyear)
        call ensanomalclim(data,npermax,nperyear,yrbeg,yrend,0,0,yr1,yr2
     +       ,mean)
        end

        subroutine ensanomalclim(data,npermax,nperyear,yrbeg,yrend,nens1
     +       ,nens2,yr1,yr2,mean)
*
*       subtract the climatology determined from the data itself over yr1:yr2
*       and also return the climatology
*
        integer npermax,nperyear,yrbeg,yrend,yr1,yr2,nens1,nens2
        real data(npermax,yrbeg:yrend,0:nens2),mean(nperyear)
        integer i,j,n,nn,firstyr,lastyr,iens
        real s,absent,minfac
        parameter (absent=3e33)
        logical lwrite
        lwrite = .false.
*
        if ( lwrite ) call printdat('anomal: before',
     +		data,npermax,nperyear,yrbeg,yrend)
        do j=1,nperyear
            do iens=nens1,nens2
                do firstyr=yr1,yr2
                    if ( data(j,firstyr,iens).lt.1e33 ) goto 10
                enddo
            enddo
            if ( lwrite ) write(*,*) 'no valid data for j=',j
            goto 800
   10       continue
            do iens=nens1,nens2
                do lastyr=yr2,yr1,-1
                    if ( data(j,lastyr,iens).lt.1e33 ) goto 20
                enddo
            enddo
            if ( lwrite ) write(*,*) 'should never come here'
            goto 800
   20       continue
            if ( lwrite ) write(*,*) 'firstyr,lastyr = ',firstyr,lastyr
            n = 0
            nn = 0
            s = 0
            do iens=nens1,nens2
                do i=firstyr,lastyr
                    nn = nn + 1
                    if ( data(j,i,iens).lt.0.9*absent ) then
                        n = n + 1
                        s = s + data(j,i,iens)
                    endif
                enddo
            enddo
***            minfac = min(0.9,max(0.1,1.5-log(1+nperyear*(lastyr-firstyr
***     +            +1.)*n/nn/nperyear)/4*n/nn))
****           otherwise lots of data gets too thin...
***            minfac = minfac/2
***            if ( lwrite ) write(*,*) 'minfac,n,nn = ',minfac,n,nn
            minfac = 0   ! let's take anomalies whenever there are 2 datapoints
            if ( n.gt.minfac*nn .and. n.gt.1 ) then
                s = s/n
                do iens=nens1,nens2
                    do i=yrbeg,yrend
                        data(j,i,iens) = data(j,i,iens) - s
                    enddo
                enddo
                mean(j) = s
            else
                do iens=nens2,nens2
                    do i=yrbeg,yrend
                        data(j,i,iens) = 3e33
                    enddo
                enddo
                mean(j) = 3e33
            endif
  800       continue
        enddo
        if ( lwrite ) call printdat('anomal: after',
     +		data,npermax,nperyear,yrbeg,yrend)
        end
*
        subroutine printdat(string,data,npermax,nperyear,yrbeg,yrend)
        implicit none
        integer npermax,nperyear,yrbeg,yrend
        real data(npermax,yrbeg:yrend)
        integer i,j
        character string*(*)
        logical lvalid
*
        print *,string,npermax,nperyear,yrbeg,yrend
        do i=yrbeg,yrend
            lvalid = .FALSE.
            do j=1,nperyear
                if ( data(j,i).lt.1e33 ) lvalid = .TRUE.
            enddo
            if ( lvalid ) print '(i5,366f8.1)',
     +              i,(data(j,i),j=1,nperyear)
        enddo
	end
