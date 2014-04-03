c--------------------------------------------------------------
c
c provisional CHANGE LOG
c
c  original:     GvdS 
c  modification: GvdS 20100706 
c                routine write2file: data is to file from 1900 to 2020
c                the start in 1900 (regardless when the first non-absent data
c                is found) was already implemented, the end in 2020 is new. The
c                last year written to file was the last year for which there was data
c  modification: GvdS 20120905
c                routine adjustdate: this routine is rewritten to handle both days shift
c                forward in time and backwards in time. There is no limit the the length
c                of the shift
c  modification: GvdS 20120919
c                routine addtospell compressed and made more efficient
c  modification: GvdS 20130606
c                added new routines which used to be part of indexSPI
c
c--------------------------------------------------------------
      subroutine readData(infile,fyear,lyear,stagrp,a,qc)
      implicit none
      include 'comgeneral.h'

      real*8        a(yrbeg:yrend,12,31)
      integer       qc(yrbeg:yrend,12,31)
      real*8        absenr,absentr
      integer       i,j,k,n,absen,fyear,lyear,stagrp
      integer       year,month,day,value,qcval,date
      character*128 infile
      character     c1*1
      logical       cold

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

c initialize array
      do i=yrbeg,yrend
        do j=1,12
          do k=1,31
            a(i,j,k) = absentr
            qc(i,j,k) = 9
          enddo
        enddo
      enddo

      cold = .true.

      open(1,file=infile,status='old')

c skip headers
      i = 1
 100  read(1,*,err=921,end=200) stagrp,date,value,qcval

        year = date/10000
        month = (date-year*10000)/100
        day = date-year*10000-month*100

c       write(6,1000) year,month,day,value,qcval
c1000   format(I4,I2,I2,I6,I3)

        if(year+1.lt.yrbeg) call error('decrease yrbeg')

        if(value.gt.absen) then
          if(cold) then
            fyear = year
            cold=.false.
          endif
          a(year,month,day) = dble(value)
          qc(year,month,day) = qcval
          lyear = year
        else
          a(year,month,day) = absentr
          qc(year,month,day) = 9
        endif
      goto 100

  200 close(1)
      goto 900

 921  write(6,*) 'error reading file ',infile

 900  continue

      if((fyear.eq.0).or.(lyear.eq.0)) 
     +   call error('readData: file empty?')

      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine convert2char(i,word)
      implicit none

      integer         i
      character       word*6,one*1,two*2,three*3,four*4

      if(i.lt.10) then
        write (one, '(I1)') i
        word = '00000'//one
      elseif(i.lt.100) then
        write (two, '(I2)') i
        word = '0000'//two
      elseif(i.lt.1000) then
        write (three, '(I3)') i
        word = '000'//three
      elseif(i.lt.10000) then
        write (four, '(I4)') i
        word = '00'//four
      else
        write(6,*) 'file number out of range ',i
      endif

      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine setcalendar(cal)
      !
      ! set calendar type, default gregorian
      !
      character cal*(*)
      character calendar*20
      common /lengthofmonth_calendar/ calendar
      save /lengthofmonth_calendar/
      calendar = cal
      if ( calendar.eq.' ' ) calendar = 'gregorian'
      end
c
c---------------------------------------------------------------------------
c
      subroutine lengthofmonth(iyear,imonth,length)
      implicit none

      integer         iyear,imonth,length
      logical         leap
      character calendar*20
      common /lengthofmonth_calendar/ calendar
      save /lengthofmonth_calendar/
	  !
      !	see http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.6/cf-conventions.html#calendar
      !
      if ( len_trim(calendar).lt.4 .or. calendar.eq.' ' ) then
        calendar = 'gregorian'
      end if
	  if ( calendar.eq.'gregorian' .or. calendar.eq.'standard' ) then
	    if ( iyear.le.1582 ) then
	      write(0,*) 'lengthofmonth: error: cannot handle change '//
     +      'from Julian to Gregorian calendar in 1582 yet'
	      call abort
	    end if
	  end if
	  if ( calendar.eq.'julian' ) then
	    write(0,*) 'lengthofmonth: error: cannot handle julian '//
     +    'calendar yet'
	    call abort
	  end if

      if ( calendar.eq.'360_day' ) then
        length = 30
      else
        ! original code
      if(imonth.le.6) then
        if(imonth.eq.1) then
          length = 31
        elseif(imonth.eq.2) then
          if ( calendar.eq.'noleap' .or. calendar.eq.'365_day' ) then
            length = 28
          else if ( calendar.eq.'all_leap' .or. calendar.eq.'366_day' ) 
     +      then
            length = 29
          else if( calendar.eq.'gregorian' .or. calendar.eq.'standard'
     +        .or. calendar.eq.'proleptic_gregorian' ) then
            call leapyr(iyear,leap)
            if(leap) then
              length = 29
            else
              length = 28
            endif
          else
            write(0,*) 'lengthofmonth: error: cannot handle calendar ',
     +        trim(calendar),' yet'
          end if
        elseif(imonth.eq.3) then
          length = 31
        elseif(imonth.eq.4) then
          length = 30
        elseif(imonth.eq.5) then
          length = 31
        else
          length = 30
        endif
      else
        if(imonth.eq.7) then
          length = 31
        elseif(imonth.eq.8) then
          length = 31
        elseif(imonth.eq.9) then
          length = 30
        elseif(imonth.eq.10) then
          length = 31
        elseif(imonth.eq.11) then
          length = 30
        else
          length = 31
        endif
      endif
      end if

      !!!print *,'lengthofmonth: returning ',length
      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine leapyr(iyear,leap)
      implicit none

      integer         iyear
      logical         leap

      if(mod(iyear,4).eq.0) then
        leap = .true.
c this takes care of years like 1700-1800-1900 no leap years
c and 1600 and 2000: leap years
        if(mod(iyear,100).eq.0) leap = .false.
        if(mod(iyear,400).eq.0) leap = .true.
      else
        leap = .false.
      endif

      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine bounds(season,firstmonth,lastmonth)
      implicit none

      integer         season,firstmonth,lastmonth

      if(season.le.4) then
        if(season.eq.1) then
c whole year
          firstmonth = 1
          lastmonth = 12
        elseif(season.eq.3) then
c warm half year
          firstmonth = 4
          lastmonth = 9
        else
          call error('invalid input')
        endif
      else
        if(season.eq.5) then
c spring - MAM
          firstmonth = 3
          lastmonth = 5
        elseif(season.eq.6) then
c summer - JJA
          firstmonth = 6
          lastmonth = 8
        elseif(season.eq.7) then
c autumn - SON
          firstmonth = 9
          lastmonth = 11
        elseif(season.ge.8) then
c monthly "seasons"
          firstmonth = season - 7
          lastmonth = firstmonth
        endif
      endif

      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine calcSeason(fyear,lyear,c,b)
c this routine sums the info in array c over the appropriate seasons and outputs
c the various time-averages in array b
c
      implicit none
      include 'comgeneral.h'

      real*8        b(yrbeg:yrend,nseason)
      integer       c(yrbeg:yrend,12,2)
      integer       i,j,k,length,absen,nm,xm,fyear,fm,lm,lyear
      real*8        absenr,absentr
      logical       pass

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

      do i=fyear,lyear
        do j=1,nseason
          pass = .false.
          xm = 0
          nm = 0
          if(j.eq.1) then
            do k=1,12
              if(c(i,k,1).gt.absen) then
                xm = xm + c(i,k,1)
                nm = nm + c(i,k,2)
              endif
            enddo
            if(nm.ge.pyear) pass = .true.

          elseif(j.eq.2) then
            if(i.gt.yrbeg)then
              do k=10,12
                if(c(i-1,k,1).gt.absen) then
                  xm = xm + c(i-1,k,1)
                  nm = nm + c(i-1,k,2)
                endif
              enddo
            endif
            do k=1,3
              if(c(i,k,1).gt.absen) then
                xm = xm + c(i,k,1)
                nm = nm + c(i,k,2)
              endif
            enddo
            if(nm.ge.p6month) pass = .true.

          elseif(j.eq.3) then
            do k=4,9
              if(c(i,k,1).gt.absen) then
                xm = xm + c(i,k,1)
                nm = nm + c(i,k,2)
              endif
            enddo
            if(nm.ge.p6month) pass = .true.

          elseif(j.eq.4) then
            if(i.gt.yrbeg)then
              if(c(i-1,12,1).gt.absen) then
                xm = c(i-1,12,1)
                nm = c(i-1,12,2)
              endif
            endif
            do k=1,2
              if(c(i,k,1).gt.absen) then
                xm = xm + c(i,k,1)
                nm = nm + c(i,k,2)
              endif
            enddo
            if(nm.ge.pseason) pass = .true.

          elseif(j.le.7)then
            call bounds(j,fm,lm)
            do k=fm,lm
              if(c(i,k,1).gt.absen) then
                xm = xm + c(i,k,1)
                nm = nm + c(i,k,2)
              endif
            enddo
            if(nm.ge.pseason) pass = .true.

          else
            k = j - 7
            if(c(i,k,1).gt.absen) then
              xm = c(i,k,1)
              nm = c(i,k,2)
            endif
            if(nm.ge.pmonth) pass = .true.

          endif

          if(pass) then
            b(i,j) = dble(xm)
          else
            b(i,j) = absentr
          endif

        enddo
      enddo

      return
      end
c
c---------------------------------------------------------------------
c
      subroutine calcSeason_GJ(fyear,lyear,nswitch,c,b)
c this routine sums the info in array c over the appropriate seasons and outputs
c the various time-averages in array b
c
c GvdS 11/07/2013 added nswitch for output of
c    annual values only (nswitch = 1)
c    halfyear values only (nswitch = 2)
c    seasonal values only (nswitch = 3)
c    monthly values only (nswitch = 4)
c    when e.g. annual values are requested, then array b is filled with missing for all seasons 
c
      implicit none
      include 'comgeneral.h'

      real*8        b(yrbeg:yrend,nseason)
      integer       c(yrbeg:yrend,12,2)
      integer       i,j,k,length,absen,nm,xm,fyear,fm,lm,lyear,nswitch
      real*8        absenr,absentr

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

c initialize
      do i=yrbeg,yrend
        do j=1,nseason
           b(i,j) = absentr
        enddo
      enddo

      if(nswitch.eq.1) then
         j=1
         do i=fyear,lyear
            xm = 0
            nm = 0
            do k=1,12
              if(c(i,k,1).gt.absen) then
                xm = xm + c(i,k,1)
                nm = nm + c(i,k,2)
              endif
            enddo
            if(nm.ge.pyear) b(i,j) = dble(xm)
         enddo
      elseif(nswitch.eq.2) then
         j=2
         xm = 0
         nm = 0
         if(i.gt.yrbeg)then
           do k=10,12
             if(c(i-1,k,1).gt.absen) then
               xm = xm + c(i-1,k,1)
               nm = nm + c(i-1,k,2)
             endif
           enddo
         endif
         do k=1,3
           if(c(i,k,1).gt.absen) then
             xm = xm + c(i,k,1)
             nm = nm + c(i,k,2)
           endif
         enddo
         if(nm.ge.p6month) b(i,j) = dble(xm)
  
         j=3
         xm = 0
         nm = 0
         do k=4,9
           if(c(i,k,1).gt.absen) then
             xm = xm + c(i,k,1)
             nm = nm + c(i,k,2)
           endif
         enddo
         if(nm.ge.p6month) b(i,j) = dble(xm)

      elseif(nswitch.eq.3)then
         j=4
         xm = 0
         nm = 0
         if(i.gt.yrbeg)then
           if(c(i-1,12,1).gt.absen) then
             xm = c(i-1,12,1)
             nm = c(i-1,12,2)
           endif
         endif
         do k=1,2
           if(c(i,k,1).gt.absen) then
             xm = xm + c(i,k,1)
             nm = nm + c(i,k,2)
           endif
         enddo
         if(nm.ge.pseason) b(i,j) = dble(xm)

         do j=5,7
           xm = 0
           nm = 0
           call bounds(j,fm,lm)
           do k=fm,lm
             if(c(i,k,1).gt.absen) then
               xm = xm + c(i,k,1)
               nm = nm + c(i,k,2)
             endif
           enddo
           if(nm.ge.pseason) b(i,j) = dble(xm)
         enddo

      else
         do j=8,19
           xm = 0
           nm = 0
           if(c(i,j,1).gt.absen) then
             xm = c(i,j,1)
             nm = c(i,j,2)
           endif
           if(nm.ge.pmonth) b(i,j) = dble(xm)
         enddo
      endif

      return
      end
c
c---------------------------------------------------------------------
c
      subroutine calcSeasonR(fyear,lyear,c,b)
c this routine is similiar to calcSeason, but here the input are real*8 arrays
c
      implicit none
      include 'comgeneral.h'

      real*8        b(yrbeg:yrend,nseason)
      real*8        c(yrbeg:yrend,12,2)
      real*8        nm,xm
      integer       i,j,k,length,absen,fyear,fm,lm,lyear
      real*8        absenr,absentr
      logical       pass

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

      do i=fyear,lyear
        do j=1,nseason
          pass= .false.
          xm = 0.0d0
          nm = 0.0d0
          if(j.eq.1) then
            do k=1,12
              if(c(i,k,1).gt.absenr) then
                xm = xm + c(i,k,1)
                nm = nm + c(i,k,2)
              endif
            enddo
            if(nint(nm).ge.pyear) pass = .true.

          elseif(j.eq.2) then
            if(i.gt.yrbeg)then
              do k=10,12
                if(c(i-1,k,1).gt.absenr) then
                  xm = xm + c(i-1,k,1)
                  nm = nm + c(i-1,k,2)
                endif
              enddo
            endif
            do k=1,3
              if(c(i,k,1).gt.absenr) then
                xm = xm + c(i,k,1)
                nm = nm + c(i,k,2)
              endif
            enddo
            if(nint(nm).ge.p6month) pass = .true.

          elseif(j.eq.3) then
            do k=4,9
              if(c(i,k,1).gt.absenr) then
                xm = xm + c(i,k,1)
                nm = nm + c(i,k,2)
              endif
            enddo
            if(nint(nm).ge.p6month) pass = .true.

          elseif(j.eq.4) then
            if(i.gt.yrbeg)then
              if(c(i-1,12,1).gt.absenr) then
                xm = c(i-1,12,1)
                nm = c(i-1,12,2)
              endif
            endif
            do k=1,2
              if(c(i,k,1).gt.absenr) then
                xm = xm + c(i,k,1)
                nm = nm + c(i,k,2)
              endif
            enddo
            if(nint(nm).ge.pseason) pass = .true.

         elseif(j.le.7) then
            call bounds(j,fm,lm)
            do k=fm,lm
              if(c(i,k,1).gt.absenr) then
                xm = xm + c(i,k,1)
                nm = nm + c(i,k,2)
              endif
            enddo
            if(nint(nm).ge.pseason) pass = .true.

          else
            k = j - 7
            if(c(i,k,1).gt.absen) then
              xm = c(i,k,1)
              nm = c(i,k,2)
            endif
            if(nm.ge.pmonth) pass = .true.

          endif

          if(pass) then
            b(i,j) = xm
          else
            b(i,j) = absentr
          endif

        enddo
      enddo

c test
c     do i=yrbeg,yrend
c       write(6,123) i,(b(i,j), j=1,nseason)
c     enddo
c123  format(I6,7e12.2)

      return
      end
c
c---------------------------------------------------------------------
c
      subroutine calcAverage(fyear,lyear,c,b)
c this routine is similiar to calcSeasonR, but here the input 
c is averaged rather than summed
c
c this routine gives the average *per day*
      implicit none
      include 'comgeneral.h'

      real*8        b(yrbeg:yrend,nseason)
      real*8        c(yrbeg:yrend,12,2)
      real*8        nm,xm
      integer       i,j,k,length,absen,fyear,fm,lm,lyear
      real*8        absenr,absentr
      logical       pass

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

      do i=fyear,lyear
        do j=1,nseason
          pass= .false.
          xm = 0.0d0
          nm = 0.0d0
          if(j.eq.1) then
            do k=1,12
              if(c(i,k,1).gt.absenr) then
                xm = xm + c(i,k,1)
                nm = nm + c(i,k,2)
              endif
            enddo
            if(nint(nm).ge.pyear) pass = .true.

          elseif(j.eq.2) then
            if(i.gt.yrbeg)then
              do k=10,12
                if(c(i-1,k,1).gt.absenr) then
                  xm = xm + c(i-1,k,1)
                  nm = nm + c(i-1,k,2)
                endif
              enddo
            endif
            do k=1,3
              if(c(i,k,1).gt.absenr) then
                xm = xm + c(i,k,1)
                nm = nm + c(i,k,2)
              endif
            enddo
            if(nint(nm).ge.p6month) pass = .true.

          elseif(j.eq.3) then
            do k=4,9
              if(c(i,k,1).gt.absenr) then
                xm = xm + c(i,k,1)
                nm = nm + c(i,k,2)
              endif
            enddo
            if(nint(nm).ge.p6month) pass = .true.

          elseif(j.eq.4) then
            if(i.gt.yrbeg)then
              if(c(i-1,12,1).gt.absenr) then
                xm = c(i-1,12,1)
                nm = c(i-1,12,2)
              endif
            endif
            do k=1,2
              if(c(i,k,1).gt.absenr) then
                xm = xm + c(i,k,1)
                nm = nm + c(i,k,2)
              endif
            enddo
            if(nint(nm).ge.pseason) pass = .true.

         elseif(j.le.7) then
            call bounds(j,fm,lm)
            do k=fm,lm
              if(c(i,k,1).gt.absenr) then
                xm = xm + c(i,k,1)
                nm = nm + c(i,k,2)
              endif
            enddo
            if(nint(nm).ge.pseason) pass = .true.

          else
            k = j - 7
            if(c(i,k,1).gt.absen) then
              xm = c(i,k,1)
              nm = c(i,k,2)
            endif
            if(nint(nm).ge.pmonth) pass = .true.

          endif

          if(pass) then
            b(i,j) = xm/dble(nm)
          else
            b(i,j) = absentr
          endif

        enddo
      enddo

      return
      end
c
c---------------------------------------------------------------------
c
      subroutine calcAverageSPI(fyear,lyear,c,b)
c this routine is similiar to calcSeasonR, but here the input 
c is averaged rather than summed
c
      implicit none
      include 'comgeneral.h'

      real*8        b(yrbeg:yrend,nseason)
      real*8        c(yrbeg:yrend,12,2)
      real*8        nm,xm
      integer       sl(nseason)
      integer       i,j,k,length,absen,fyear,fm,lm,lyear
      real*8        absenr,absentr
      logical       pass
      data          (sl(j), j=1,nseason) /
     +        12,6,6,3,3,3,3,1,1,1,1,1,1,1,1,1,1,1,1/

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

      do i=fyear,lyear
        do j=1,nseason
          pass= .false.
          xm = 0.0d0
          nm = 0.0d0
          if(j.eq.1) then
            do k=1,12
              if(c(i,k,1).gt.absenr) then
                xm = xm + c(i,k,1)
                nm = nm + c(i,k,2)
              endif
            enddo
            if(nint(nm).ge.pyear) pass = .true.

          elseif(j.eq.2) then
            if(i.gt.yrbeg)then
              do k=10,12
                if(c(i-1,k,1).gt.absenr) then
                  xm = xm + c(i-1,k,1)
                  nm = nm + c(i-1,k,2)
                endif
              enddo
            endif
            do k=1,3
              if(c(i,k,1).gt.absenr) then
                xm = xm + c(i,k,1)
                nm = nm + c(i,k,2)
              endif
            enddo
            if(nint(nm).ge.p6month) pass = .true.

          elseif(j.eq.3) then
            do k=4,9
              if(c(i,k,1).gt.absenr) then
                xm = xm + c(i,k,1)
                nm = nm + c(i,k,2)
              endif
            enddo
            if(nint(nm).ge.p6month) pass = .true.

          elseif(j.eq.4) then
            if(i.gt.yrbeg)then
              if(c(i-1,12,1).gt.absenr) then
                xm = c(i-1,12,1)
                nm = c(i-1,12,2)
              endif
            endif
            do k=1,2
              if(c(i,k,1).gt.absenr) then
                xm = xm + c(i,k,1)
                nm = nm + c(i,k,2)
              endif
            enddo
            if(nint(nm).ge.pseason) pass = .true.

         elseif(j.le.7) then
            call bounds(j,fm,lm)
            do k=fm,lm
              if(c(i,k,1).gt.absenr) then
                xm = xm + c(i,k,1)
                nm = nm + c(i,k,2)
              endif
            enddo
            if(nint(nm).ge.pseason) pass = .true.

          else
            k = j - 7
            if(c(i,k,1).gt.absen) then
              xm = c(i,k,1)
              nm = c(i,k,2)
            endif
            if(nint(nm).ge.pmonth) pass = .true.

          endif

          if(pass) then
            b(i,j) = xm/dble(sl(j))
          else
            b(i,j) = absentr
          endif

        enddo
      enddo

      return
      end
c
c---------------------------------------------------------------------
c
      subroutine calcSeasonR_GJ(fyear,lyear,nswitch,c,b)
c this routine is similiar to calcSeason, but here the input are real*8 arrays
c
c GvdS 11/07/2013 added nswitch for output of
c    annual values only (nswitch = 1)
c    halfyear values only (nswitch = 2)
c    seasonal values only (nswitch = 3)
c    monthly values only (nswitch = 4)
c    when e.g. annual values are requested, then array b is filled with missing for all seasons 
c
      implicit none
      include 'comgeneral.h'

      real*8        b(yrbeg:yrend,nseason)
      real*8        c(yrbeg:yrend,12,2)
      integer       i,j,k,length,absen,fyear,fm,lm,lyear,nswitch
      real*8        nm,xm
      real*8        absenr,absentr

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

c initialize
      do i=yrbeg,yrend
        do j=1,nseason
           b(i,j) = absentr
        enddo
      enddo

      if(nswitch.eq.1) then
         j=1
         do i=fyear,lyear
            xm = 0.0
            nm = 0.0
            do k=1,12
              if(c(i,k,1).gt.absenr) then
                xm = xm + c(i,k,1)
                nm = nm + c(i,k,2)
              endif
            enddo
            if(nint(nm).ge.pyear) b(i,j) = xm
         enddo
      elseif(nswitch.eq.2) then
         j=2
         xm = 0.0
         nm = 0.0
         if(i.gt.yrbeg)then
           do k=10,12
             if(c(i-1,k,1).gt.absenr) then
               xm = xm + c(i-1,k,1)
               nm = nm + c(i-1,k,2)
             endif
           enddo
         endif
         do k=1,3
           if(c(i,k,1).gt.absenr) then
             xm = xm + c(i,k,1)
             nm = nm + c(i,k,2)
           endif
         enddo
         if(nint(nm).ge.p6month) b(i,j) = xm
  
         j=3
         xm = 0.0
         nm = 0.0
         do k=4,9
           if(c(i,k,1).gt.absenr) then
             xm = xm + c(i,k,1)
             nm = nm + c(i,k,2)
           endif
         enddo
         if(nint(nm).ge.p6month) b(i,j) = xm

      elseif(nswitch.eq.3)then
         j=4
         xm = 0.0
         nm = 0.0
         if(i.gt.yrbeg)then
           if(c(i-1,12,1).gt.absenr) then
             xm = c(i-1,12,1)
             nm = c(i-1,12,2)
           endif
         endif
         do k=1,2
           if(c(i,k,1).gt.absenr) then
             xm = xm + c(i,k,1)
             nm = nm + c(i,k,2)
           endif
         enddo
         if(nint(nm).ge.pseason) b(i,j) = xm

         do j=5,7
           xm = 0.0
           nm = 0.0
           call bounds(j,fm,lm)
           do k=fm,lm
             if(c(i,k,1).gt.absenr) then
               xm = xm + c(i,k,1)
               nm = nm + c(i,k,2)
             endif
           enddo
           if(nint(nm).ge.pseason) b(i,j) = xm
         enddo

      else
         do j=8,19
           xm = 0.0
           nm = 0.0
           if(c(i,j,1).gt.absenr) then
             xm = c(i,j,1)
             nm = c(i,j,2)
           endif
           if(nint(nm).ge.pmonth) b(i,j) = xm
         enddo
      endif

      return
      end
c
c---------------------------------------------------------------------
c
      subroutine calcAverage_GJ(fyear,lyear,nswitch,c,b)
c this routine is similiar to calcSeason, but here the input are real*8 arrays
c
c GvdS 11/07/2013 added nswitch for output of
c    annual values only (nswitch = 1)
c    halfyear values only (nswitch = 2)
c    seasonal values only (nswitch = 3)
c    monthly values only (nswitch = 4)
c    when e.g. annual values are requested, then array b is filled with missing for all seasons 
c
      implicit none
      include 'comgeneral.h'

      real*8        b(yrbeg:yrend,nseason)
      real*8        c(yrbeg:yrend,12,2)
      integer       i,j,k,length,absen,fyear,fm,lm,lyear,nswitch
      real*8        nm,xm
      real*8        absenr,absentr

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

c initialize
      do i=yrbeg,yrend
        do j=1,nseason
           b(i,j) = absentr
        enddo
      enddo

      if(nswitch.eq.1) then
         j=1
         do i=fyear,lyear
            xm = 0.0
            nm = 0.0
            do k=1,12
              if(c(i,k,1).gt.absenr) then
                xm = xm + c(i,k,1)
                nm = nm + c(i,k,2)
              endif
            enddo
            if(nint(nm).ge.pyear) b(i,j) = xm
         enddo
      elseif(nswitch.eq.2) then
         j=2
         xm = 0.0
         nm = 0.0
         if(i.gt.yrbeg)then
           do k=10,12
             if(c(i-1,k,1).gt.absenr) then
               xm = xm + c(i-1,k,1)
               nm = nm + c(i-1,k,2)
             endif
           enddo
         endif
         do k=1,3
           if(c(i,k,1).gt.absenr) then
             xm = xm + c(i,k,1)
             nm = nm + c(i,k,2)
           endif
         enddo
         if(nint(nm).ge.p6month) b(i,j) = xm
  
         j=3
         xm = 0.0
         nm = 0.0
         do k=4,9
           if(c(i,k,1).gt.absenr) then
             xm = xm + c(i,k,1)
             nm = nm + c(i,k,2)
           endif
         enddo
         if(nint(nm).ge.p6month) b(i,j) = xm

      elseif(nswitch.eq.3)then
         j=4
         xm = 0.0
         nm = 0.0
         if(i.gt.yrbeg)then
           if(c(i-1,12,1).gt.absenr) then
             xm = c(i-1,12,1)
             nm = c(i-1,12,2)
           endif
         endif
         do k=1,2
           if(c(i,k,1).gt.absenr) then
             xm = xm + c(i,k,1)
             nm = nm + c(i,k,2)
           endif
         enddo
         if(nint(nm).ge.pseason) b(i,j) = xm

         do j=5,7
           xm = 0.0
           nm = 0.0
           call bounds(j,fm,lm)
           do k=fm,lm
             if(c(i,k,1).gt.absenr) then
               xm = xm + c(i,k,1)
               nm = nm + c(i,k,2)
             endif
           enddo
           if(nint(nm).ge.pseason) b(i,j) = xm
         enddo

      else
         do j=8,19
           xm = 0.0
           nm = 0.0
           if(c(i,j,1).gt.absenr) then
             xm = c(i,j,1)
             nm = c(i,j,2)
           endif
           if(nint(nm).ge.pmonth) b(i,j) = xm
         enddo
      endif

      return
      end
c
c---------------------------------------------------------------------
c
      subroutine calcAverageSPI_GJ(fyear,lyear,nswitch,c,b)
c this routine is similiar to calcSeasonR, but here the input 
c is averaged rather than summed
c
c this routine gives the average *per day*
c
c GvdS 11/07/2013 added nswitch for output of
c    annual values only (nswitch = 1)
c    halfyear values only (nswitch = 2)
c    seasonal values only (nswitch = 3)
c    monthly values only (nswitch = 4)
c    when e.g. annual values are requested, then array b is filled with missing for all seasons 
c
      implicit none
      include 'comgeneral.h'

      real*8        b(yrbeg:yrend,nseason)
      real*8        c(yrbeg:yrend,12,2)
      integer       i,j,k,length,absen,fyear,fm,lm,lyear,nswitch
      real*8        nm,xm
      real*8        absenr,absentr
      integer       sl(nseason)
      data          (sl(j), j=1,nseason) /
     +        12,6,6,3,3,3,3,1,1,1,1,1,1,1,1,1,1,1,1/

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

c initialize
      do i=yrbeg,yrend
        do j=1,nseason
           b(i,j) = absentr
        enddo
      enddo

      if(nswitch.eq.1) then
         j=1
         do i=fyear,lyear
            xm = 0.0
            nm = 0.0
            do k=1,12
              if(c(i,k,1).gt.absenr) then
                xm = xm + c(i,k,1)
                nm = nm + c(i,k,2)
              endif
            enddo
            if(nint(nm).ge.pyear) b(i,j) = xm/dble(sl(j))
         enddo
      elseif(nswitch.eq.2) then
         j=2
         xm = 0.0
         nm = 0.0
         if(i.gt.yrbeg)then
           do k=10,12
             if(c(i-1,k,1).gt.absenr) then
               xm = xm + c(i-1,k,1)
               nm = nm + c(i-1,k,2)
             endif
           enddo
         endif
         do k=1,3
           if(c(i,k,1).gt.absenr) then
             xm = xm + c(i,k,1)
             nm = nm + c(i,k,2)
           endif
         enddo
         if(nint(nm).ge.p6month) b(i,j) = xm/dble(sl(j))
  
         j=3
         xm = 0.0
         nm = 0.0
         do k=4,9
           if(c(i,k,1).gt.absenr) then
             xm = xm + c(i,k,1)
             nm = nm + c(i,k,2)
           endif
         enddo
         if(nint(nm).ge.p6month) b(i,j) = xm/dble(sl(j))

      elseif(nswitch.eq.3)then
         j=4
         xm = 0.0
         nm = 0.0
         if(i.gt.yrbeg)then
           if(c(i-1,12,1).gt.absenr) then
             xm = c(i-1,12,1)
             nm = c(i-1,12,2)
           endif
         endif
         do k=1,2
           if(c(i,k,1).gt.absenr) then
             xm = xm + c(i,k,1)
             nm = nm + c(i,k,2)
           endif
         enddo
         if(nint(nm).ge.pseason) b(i,j) = xm/dble(sl(j))

         do j=5,7
           xm = 0.0
           nm = 0.0
           call bounds(j,fm,lm)
           do k=fm,lm
             if(c(i,k,1).gt.absenr) then
               xm = xm + c(i,k,1)
               nm = nm + c(i,k,2)
             endif
           enddo
           if(nint(nm).ge.pseason) b(i,j) = xm/dble(sl(j))
         enddo

      else
         do j=8,19
           xm = 0.0
           nm = 0.0
           if(c(i,j,1).gt.absenr) then
             xm = c(i,j,1)
             nm = c(i,j,2)
           endif
           if(nint(nm).ge.pmonth) b(i,j) = xm
         enddo
      endif

      return
      end
c
c---------------------------------------------------------
c 
      subroutine addtospell(val,qc,pval,pqc,gt,thresh,nlength)
c this routine determines if the value "val" contributes to a spell
c
c if gt = .true., then the inequality (val.ge.thresh) is tested
c if gt = .false., then the inequality (val.lt.thresh) is tested
      implicit none
      include 'comgeneral.h'

      integer        nlength,absen,qc,pqc
      real*8         val,pval,absentr,absenr,thresh
      logical        gt

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

c if val is absent, then it does not contribute to the spell 
      if(qc.eq.0) then

        if(     ((gt).and.(val.gt.thresh)).or. 
     +     ((.not.gt).and.(val.lt.thresh))) then

c if the previous val is absent, then it does not contribute to the spell
          if(pqc.eq.0) then

            if(     ((gt).and.(pval.gt.thresh)).or. 
     +         ((.not.gt).and.(pval.lt.thresh))) then
c a(i,j,k) contributes to an established spell
              nlength = nlength + 1
            else
c a(i,j,k) starts a new spell
              nlength = 1
            endif
          endif
        else
c a(i,j,k) ends an established spell since it fails to meet the inequality
          nlength = 0
        endif
      else
c a(i,j,k) is absent of suspect and no spell is started or continued
        nlength = 0
      endif

      return
      end
c
c---------------------------------------------------------------------
c
c     subroutine addtospell(val,qc,pval,pqc,gt,thresh,nlength)
c this routine determines if the value "val" contributes to a spell
c
c if gt = .true., then the inequality (val.ge.thresh) is tested
c if gt = .false., then the inequality (val.lt.thresh) is tested
c     implicit none
c     include 'comgeneral.h'

c     integer        nlength,absen,qc,pqc
c     real*8         val,pval,absentr,absenr,thresh
c     logical        now,prev,gt

c     absen = absent + 1
c     absenr = dble(absen)
c     absentr = dble(absent)

c if val is absent, then it does not contribute to the spell 
c     now = .false.
c     if(qc.eq.0) then

c       if(     ((gt).and.(val.gt.thresh)).or. 
c    +     ((.not.gt).and.(val.lt.thresh))) now = .true.

c     endif

c it only makes sense to continue if now = .true.
c     if (now) then

c if the previous val is absent, then it does not contribute to the spell
c       prev = .false.
c       if(pqc.eq.0) then

c         if(     ((gt).and.(pval.gt.thresh)).or. 
c    +       ((.not.gt).and.(pval.lt.thresh))) prev = .true.

c       endif

c       if(prev) then
c a(i,j,k) contributes to an established spell: now = .true. and prev = .true.
c         nlength = nlength + 1
c       else
c a(i,j,k) starts a new spell: now = .true. and prev = .false.
c         nlength = 1
c       endif

c     else
c either a(i,j,k) ends an established spell or
c both a(i,j,k) and its previous value are either absent or fail to meet the threshold:
c no spell is started or continued
c       nlength = 0
c     endif

c     return
c     end
c
c---------------------------------------------------------------------------
c
      subroutine adjustdate(io,jo,ko,no,ii,jj,kk)
c given the date (yr, mo, dy) = (i,j,k), this routine returns the date
c of the day n days further or n days backwards in (ii,jj,kk)
      implicit none

      integer       io,jo,ko,no
      integer       i,j,k,n,ii,jj,kk,length

c store the original input (i,j,k) in (io,j0,k0)
      i = io
      j = jo
      k = ko
      n = no
        !!!print *,'adjustdate: input i,j,k,n = ',i,j,k,n

c if n=0: then we have the simplest case
      if(n.eq.0) then
        ii = i
        jj = j
        kk = k
      elseif(n.gt.0)then
 500    call lengthofmonth(i,j,length)
        if(k+n.le.length)then
          ii = i
          jj = j
          kk = k + n
          goto 501
        else
          j=j + 1
          if(j.eq.13) then
            j=1
            i=i + 1
          endif
          n = n - (length - k) - 1
          k = 1
        endif
        goto 500
      elseif(n.lt.0)then
 502    if(k+n.ge.1)then
          ii = i
          jj = j
          kk = k + n
          goto 501
        else
          j = j - 1
          if(j.eq.0)then
            j=12
            i=i-1
          endif
          n = n + k
          call lengthofmonth(i,j,length)
          k = length
          goto 502
        endif
      endif

 501  continue

      return
      end
c
c---------------------------------------------------------------------------
c
c     subroutine adjustdate(i,j,k,n,ii,jj,kk)
c given the date (yr, mo, dy) = (i,j,k), this routine returns the date
c of the day n days further in (ii,jj,kk)
c     implicit none

c     integer       i,j,k,n,ii,jj,kk,length

c if n=0: then we have the simplest case
c     if(n.eq.0) then
c       ii = i
c       jj = j
c       kk = k
c     else
c       if(k.lt.15) then
c         if(k+n.lt.1) then
c we are at the beginning of the month
c           if(j.eq.1) then
c             ii = i - 1
c             jj = 12
c             kk = k + n + 31
c           else
c             call lengthofmonth(i,j-1,length)
c             ii = i
c             jj = j - 1
c             kk = k + n + length
c           endif
c         else
c           ii = i
c           jj = j
c           kk = k + n
c         endif
c       else
c         call lengthofmonth(i,j,length)
c         if(k+n.gt.length) then
c we are at the end of the month
c           if(j.eq.12) then
c             ii = i + 1
c             jj = 1
c             kk = k + n - 31
c           else
c             ii = i
c             jj = j + 1
c             kk = k + n - length
c           endif
c         else
c           ii = i
c           jj = j
c           kk = k + n
c         endif
c       endif
c     endif

c     return
c     end
c
c---------------------------------------------------------------------------
c
      subroutine readIndexData(infile,b,ex)
      implicit none
      include 'comgeneral.h'

      real*8        b(yrbeg:yrend,nseason)
      real*8        absentr
      integer       ndum(nseason)
      integer       i,j,absen
      character*128 infile
      logical       ex

      absentr = -9999.99

      inquire(file=infile,exist=ex)

      if(.not.ex) goto 900

c initialize array
      do i=yrbeg,yrend
        do j=1,nseason
          b(i,j) = absentr
        enddo
      enddo

      open(1,file=infile,status='old')

c read file 
 100    read(1,*,end=200) i,(ndum(j), j=1,nseason)

        do j=1,nseason
          b(i,j) = dble(ndum(j))/100.0d0
        enddo
      goto 100

  200 close(1)

 900  continue

      return
      end
c
c---------------------------------------------------------------------------
c
      integer function nzhang(i)
c this function is used in the Zhang et al. approach of determining percentiles
c it returns the argument when it is \in [calyrbeg,calyrend]. Otherwise it returns
c the value calyrend+1
c
c when the zhang et al. approach is disabled (zhang = .false.), then it returns
c the value calyrend+1
      implicit none
      include 'comgeneral.h'

      integer      i

      nzhang = calyrend+1

      if(zhang) then
        if((i.ge.calyrbeg).and.(i.le.calyrend)) nzhang = i
      endif

      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine smoothPercentiles(perc)
c this routine smoothes percentiles with a simple weighted running mean of length
c npoint
      implicit none
      include 'comgeneral.h'

      integer       mxlength
      parameter     (mxlength=(calyrend-calyrbeg+1)*nwl + 10)
      real*8        perc(calyrbeg:calyrend+1,12,31)
      real*8        pXX(12,31)

      real*8        dum(3*365)
      integer       i,j,k,l,ii,jj,kk,length,ndata,m
      integer       absen
      integer       yearday2date(365,2)
      real*8        absenr,absentr

c if zhang = .false., then all years in perc are equal. it saves time to do
c the smoothng for one year and then replace the other years
c
c if zhang = .true., then every year needs to be done separately

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

      if (zhang) then
        do m=calyrbeg,calyrend+1
c put the data in a long 1-D array
          i=0
          do j=1,12
            call lengthofmonth(1973,j,length)
            do k=1,length
              i = i + 1

c the array yearday2date is needed for the inverse transformation
              yearday2date(i,1) = j
              yearday2date(i,2) = k

              dum(i) = perc(m,j,k)
              dum(i+365) = dum(i)
              dum(i+730) = dum(i)
            enddo
          enddo

c smooth the ipercentiles with a weighted running mean
          call smooth(dum,yearday2date,pXX,npoint)

          do j=1,12
            do k=1,31
              perc(m,j,k) = pXX(j,k)
            enddo
          enddo
        enddo
      else
        m=calyrend+1
c put the data in a long 1-D array
        i=0
        do j=1,12
          call lengthofmonth(1973,j,length)
          do k=1,length
            i = i + 1

c the array yearday2date is needed for the inverse transformation
            yearday2date(i,1) = j
            yearday2date(i,2) = k

            dum(i) = perc(m,j,k)
            dum(i+365) = dum(i)
            dum(i+730) = dum(i)
          enddo
        enddo

c smooth the ipercentiles with a weighted running mean
        call smooth(dum,yearday2date,pXX,npoint)

        do m=calyrbeg,calyrend+1
          do j=1,12
            do k=1,31
              perc(m,j,k) = pXX(j,k)
            enddo
          enddo
        enddo
      endif

      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine smooth(dum,yearday2date,anncycle,np)
      implicit none
      include 'comgeneral.h'

      real*8        anncycle(12,31),dum(3*365),w(maxstep)
      integer       i,j,k,length,np,np2
      integer       yearday2date(365,2)
      real*8        xd

c determine the weights in the weighted running mean
      if(np.gt.maxstep) call error('# points in running mean too large')

      call detweights(w,np)

c smooth the data in array dum
      np2 = (np-1)/2

      do i=1,365
        xd = 0.0D0
        do k=1,np
          xd = xd + w(k)*dum(365+i-np2-1+k)
        enddo
c       write(6,*) i,anncycle(yearday2date(i,1),yearday2date(i,2)),xd
        anncycle(yearday2date(i,1),yearday2date(i,2)) = xd
      enddo

      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine detweights(w,nmean)
      implicit none
      include 'comgeneral.h'

      integer      nmean,i
      real*8       w(maxstep),x,s,dx

      if(nmean.gt.1) then
        dx=4d0/(real(nmean-1))
        s=0d0
        do i=1,nmean
          x=2d0-(nmean-i)*dx
          w(i)=exp(-x*x)
          s=s+w(i)
        enddo
        do i=1,nmean
          w(i)=w(i)/s
        enddo
      else
        w(1) = 1.0d0
      endif

      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine write2file(newdata,indid,sid,fyear,lyear,b,c)
      implicit none
      include 'comgeneral.h'

      real*8        b(yrbeg:yrend,nseason)
      real*8        c(yrbeg:yrend,nseason)
      integer       i,j,sid,fyear,lyear,indid,absen
      character*128 newdata
      real*8        absenr,absentr

      absen = absent + 1
      absenr = dble(absen)

c write to file
c note that the season numbering starts with 0 (=year)
      open(3,file=newdata,status='new')

c fill the part from 1900-fyear-1 with absent data
      do i=1900,fyear-1
        do j=1,nseason
          write(3,1003) indid,sid,j-1,i,-999999,-999999
        enddo
      enddo

      do i=fyear,lyear
        do j=1,nseason
          if(b(i,j).gt.absenr) then
            write(3,1003) indid,sid,j-1,i,
     +        nint(100*b(i,j)),nint(100*c(i,j))
          else
            write(3,1003) indid,sid,j-1,i,-999999,nint(100*c(i,j))
          endif
        enddo
      enddo

c fill the part from lyear+1 to 2020 with absent data
      do i=lyear+1,2020
        do j=1,nseason
          write(3,1003) indid,sid,j-1,i,-999999,-999999
        enddo
      enddo

      close(3)

1003  format(4(I6,','),I10,',',I10)

      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine xscale(fyear,lyear,xfac,a)
c
c this subroutine scales the array a with factor xfac, starting from year fyear
c
      implicit none
      include 'comgeneral.h'

      real*8        a(yrbeg:yrend,12,31)
      real*8        xfac
      integer       fyear,absen,i,j,k,length,lyear
      real*8        absenr,absentr

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

      do i=max(fyear,yrbeg),min(lyear,yrend)
        do j=1,12
          call lengthofmonth(i,j,length)
          do k=1,length
            if(a(i,j,k).gt.absenr) a(i,j,k)=xfac*a(i,j,k)
          enddo
        enddo
      enddo

      return
      end
c
c---------------------------------------------------------------------------
c
      integer function findindexid(cindex)
c
c this integer function relates the name of the index to the ind_id
c
      implicit none
      include 'comgeneral.h'

      integer           i
      character*10      cindex

      i=0
 10   i=i+1
      if(indexnames(i).eq.cindex) then
        findindexid = indexids(i)
        goto 30
      else
        if(i.eq.mxindices) then
          write(6,*) cindex
          goto 20
        else
          goto 10
        endif
      endif
       
 20   call error('index id not found')

 30   return

      end
c
c-----------------------------------------------------------------------
c
      subroutine lowessfilter(fyear,lyear,b,c)
c-----------------------------------------------------------------------
c     This module low-pass filters the index records using a lowess filter.
c     The code is based on routines provided by:
* wsc@research.bell-labs.com Mon Dec 30 16:55 EST 1985
* W. S. Cleveland
* Bell Laboratories
* Murray Hill NJ 07974
c
c-----------------------------------------------------------------------
      implicit none
      include 'comgeneral.h'

      real*8        b(yrbeg:yrend,nseason)
      real*8        c(yrbeg:yrend,nseason)
      integer       fyear,lyear,i,j,mt,nval,absen,niter,nrec
      parameter     (mt=yrend-yrbeg+1)
      real*8        x(mt),y(mt),ys(mt),rw(mt),res(mt)
      real*8        absenr,absentr,f,delta
      logical       prsnt(yrbeg:yrend)

*             The following are data and output from LOWESS that  can
*        be  used  to check your implementation of the routines.  The
*        notation (10)v means 10 values of v.
*
*
*
*
*        X values:
*          1  2  3  4  5  (10)6  8  10  12  14  50
*
*        Y values:
*           18  2  15  6  10  4  16  11  7  3  14  17  20  12  9  13  1  8  5  19
*
*
*        YS values with F = .25, NSTEPS = 0, DELTA = 0.0
*         13.659  11.145  8.701  9.722  10.000  (10)11.300  13.000  6.440  5.596
*           5.456  18.998
*
*        YS values with F = .25, NSTEPS = 0 ,  DELTA = 3.0
*          13.659  12.347  11.034  9.722  10.511  (10)11.300  13.000  6.440  5.596
*            5.456  18.998
*
*        YS values with F = .25, NSTEPS = 2, DELTA = 0.0
*          14.811  12.115  8.984  9.676  10.000  (10)11.346  13.000  6.734  5.744
*            5.415  18.998

c  test driver for lowess: checked all three tests on June 12th 2009
c  for expected output, see introduction
c     x(1) = 1
c     x(2) = 2
c     x(3) = 3
c     x(4) = 4
c     x(5) = 5
c     x(6) = 6
c     x(7) = 6
c     x(8) = 6
c     x(9) = 6
c     do i=10,15
c       x(i) = 6
c     enddo
c     x(16) = 8
c     x(17) = 10
c     x(18) = 12
c     x(19) = 14
c     x(20) =50 
c     y(1) = 18.0
c     y(2) = 2.0
c     y(3) = 15.0
c     y(4) = 6.0
c     y(5) = 10.0
c     y(6) = 4.0
c     y(7) = 16.0
c     y(8) = 11.0
c     y(9) = 7.0
c     y(10) = 3.0
c     y(11) = 14.0
c     y(12) = 17.0
c     y(13) = 20.0
c     y(14) = 12.0
c     y(15) = 9.0
c     y(16) = 13.0
c     y(17) = 1.0
c     y(18) = 8.0
c     y(19) = 5.0
c     y(20) = 19.0
c     call lowess(x,y,20,.25,0,0.,ys,rw,res)
c     do i=1,20
c       write(6,*) i,ys(i)
c     enddo
c     write(6,*) 
c     call lowess(x,y,20,.25,0,3.,ys,rw,res)
c     do i=1,20
c       write(6,*) i,ys(i)
c     enddo
c     write(6,*) 
c     call lowess(x,y,20,.25,2,0.,ys,rw,res)
c     do i=1,20
c       write(6,*) i,ys(i)
c     enddo
c     stop

      absen = absent + 1
      absenr = dble(absen)
c this gives nint(100*absentr) = -999999
      absentr = dble(absent)-0.99

c the R-manual suggests that delta = 1/100 of the range of x is an appropriate default value
      delta = 0.01d0*(lyear-fyear+1)
c the values below are from CalcIndex_index.R
      f = 30.0/dble(lyear-fyear+1)
      niter=3

      do j=1,nseason
        nval = 0
        do i=fyear,lyear
          prsnt(i)=.false.
          if(b(i,j).gt.absenr) then
            nval = nval + 1
            x(nval) = dble(i-fyear+1)
            y(nval) = b(i,j)
            prsnt(i)=.true.
          endif
        enddo

c if the window if 25 yrs or longer: lowess filter makes sense (AKT, from CalcIndex_index.R)
        if(nval.gt.24) then
           call lowess(x,y,nval,f,niter,delta,ys,rw,res)

           nrec = 0
           do i=fyear,lyear
             if(prsnt(i)) then
               nrec=nrec + 1
               c(i,j)= ys(nrec)
             else
               c(i,j) = absentr
             endif
           enddo
           if(nrec.ne.nval) call error('strange problem with lowess')

        else
           do i=fyear,lyear
             c(i,j) = absentr
           enddo
        endif

      enddo

      return
      end
c
c-----------------------------------------------------------------------------
c
      subroutine lowess(x, y, n, f, nsteps, delta, ys, rw, res)
      integer n
      integer nsteps
      real*8 x(n), y(n), f, delta, ys(n), rw(n)
      real*8 res(n)
      integer nright, min0, max0, i, j, ifix
      integer iter, last, m1, m2, ns, nleft
      real*8 abs, cut, cmad, r, d1, d2
      real*8 c1, c9, alpha, denom, float
      logical ok
      if (n .ge. 2) goto 1
         ys(1) = y(1)
         return
c at least two, at most n points
   1  ns = max0(min0(int(f*float(n)), n), 2)
      iter = 1
         goto  3
   2     iter = iter+1
   3     if (iter .gt. nsteps+1) goto  22
c robustness iterations
         nleft = 1
         nright = ns
c index of prev estimated point
         last = 0
c index of current point
         i = 1
   4        if (nright .ge. n) goto  5
c move nleft, nright to right if radius decreases
               d1 = x(i)-x(nleft)
c if d1<=d2 with x(nright+1)==x(nright), lowest fixes
               d2 = x(nright+1)-x(i)
               if (d1 .le. d2) goto  5
c radius will not decrease by move right
               nleft = nleft+1
               nright = nright+1
               goto  4
c fitted value at x(i)
   5        call lowest(x, y, n, x(i), ys(i), nleft, nright, res, iter
     +     .gt. 1, rw, ok)
            if (.not. ok) ys(i) = y(i)
c all weights zero - copy over value (all rw==0)
            if (last .ge. i-1) goto 9
               denom = x(i)-x(last)
c skipped points -- interpolate
c non-zero - proof?
               j = last+1
                  goto  7
   6              j = j+1
   7              if (j .ge. i) goto  8
                  alpha = (x(j)-x(last))/denom
                  ys(j) = alpha*ys(i)+(1.0-alpha)*ys(last)
                  goto  6
   8           continue
c last point actually estimated
   9        last = i
c x coord of close points
            cut = x(last)+delta
            i = last+1
               goto  11
  10           i = i+1
  11           if (i .gt. n) goto  13
c find close points
               if (x(i) .gt. cut) goto  13
c i one beyond last pt within cut
               if (x(i) .ne. x(last)) goto 12
                  ys(i) = ys(last)
c exact match in x
                  last = i
  12           continue
               goto  10
c back 1 point so interpolation within delta, but always go forward
  13        i = max0(last+1, i-1)
  14        if (last .lt. n) goto  4
c residuals
         do  15 i = 1, n
            res(i) = y(i)-ys(i)
  15        continue
         if (iter .gt. nsteps) goto  22
c compute robustness weights except last time
         do  16 i = 1, n
            rw(i) = abs(res(i))
  16        continue
         !!!call sort(rw, n)
         call nrsort(n,rw) ! use Numerical Recipes
         m1 = n/2+1
         m2 = n-m1+1
c 6 median abs resid
         cmad = 3.0*(rw(m1)+rw(m2))
         c9 = .999*cmad
         c1 = .001*cmad
         do  21 i = 1, n
            r = abs(res(i))
            if (r .gt. c1) goto 17
               rw(i) = 1.
c near 0, avoid underflow
               goto  20
  17           if (r .le. c9) goto 18
                  rw(i) = 0.
c near 1, avoid underflow
                  goto  19
  18              rw(i) = (1.0-(r/cmad)**2)**2
  19        continue
  20        continue
  21        continue
         goto  2
  22  return
      end
c
c
c
      subroutine lowest(x, y, n, xs, ys, nleft, nright, w, userw
     +, rw, ok)
      integer n
      integer nleft, nright
      real*8 x(n), y(n), xs, ys, w(n), rw(n)
      logical userw, ok
      integer nrt, j
      real*8 abs, a, b, c, h, r
      real*8 h1, sqrt, h9, amax1, range
      range = x(n)-x(1)
      h = max(xs-x(nleft), x(nright)-xs)
      h9 = .999*h
      h1 = .001*h
c sum of weights
      a = 0.0
      j = nleft
         goto  2
   1     j = j+1
   2     if (j .gt. n) goto  7
c compute weights (pick up all ties on right)
         w(j) = 0.
         r = abs(x(j)-xs)
         if (r .gt. h9) goto 5
            if (r .le. h1) goto 3
               w(j) = (1.0-(r/h)**3)**3
c small enough for non-zero weight
               goto  4
   3           w(j) = 1.
   4        if (userw) w(j) = rw(j)*w(j)
            a = a+w(j)
            goto  6
   5        if (x(j) .gt. xs) goto  7
c get out at first zero wt on right
   6     continue
         goto  1
c rightmost pt (may be greater than nright because of ties)
   7  nrt = j-1
      if (a .gt. 0.0) goto 8
         ok = .false.
         goto  16
   8     ok = .true.
c weighted least squares
         do  9 j = nleft, nrt
c make sum of w(j) == 1
            w(j) = w(j)/a
   9        continue
         if (h .le. 0.) goto 14
            a = 0.0
c use linear fit
            do  10 j = nleft, nrt
c weighted center of x values
               a = a+w(j)*x(j)
  10           continue
            b = xs-a
            c = 0.0
            do  11 j = nleft, nrt
               c = c+w(j)*(x(j)-a)**2
  11           continue
            if (sqrt(c) .le. .001*range) goto 13
               b = b/c
c points are spread out enough to compute slope
               do  12 j = nleft, nrt
                  w(j) = w(j)*(b*(x(j)-a)+1.0)
  12              continue
  13        continue
  14     ys = 0.0
         do  15 j = nleft, nrt
            ys = ys+w(j)*y(j)
  15        continue
  16  return
      end
c
c---------------------------------------------------------------------------
c
      subroutine smoothcycle(xcycle)
c this routine smooth the annual cycle by simple weighted running mean
c to the annual cycle
c
      implicit none
      include 'comgeneral.h'

      integer       mmax
      parameter     (mmax=12*31*3)
      real*8        xcycle(12,31)
      real*8        xdum(mmax),ydum(mmax)
      real*8        a(mmax)
      integer       i,j,k,l,absen,nrec,np,length
      real*8        absenr,absentr,x

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

      nrec = 0
      np = 0

c make an array of three whole years
      do i=1,3
        do j=1,12
          call lengthofmonth(1973,j,length)
          do k=1,length
            nrec = nrec + 1
            if(xcycle(j,k).gt.absenr) then
              np=np+1
              xdum(np) = dble(nrec)
              ydum(np) = xcycle(j,k)
            endif
          enddo
        enddo
      enddo

c calculate a simple 5-day running mean
      nrec=365
      do j=1,12
        call lengthofmonth(1973,j,length)
        do k=1,length
          nrec = nrec + 1

          x=0.0d0
          np = 0
          do l=-2,2
            if(ydum(nrec+l).gt.absenr) then
              np = np + 1
              x = x + ydum(nrec+l)
            endif
          enddo
          if(np.gt.0) then
            xcycle(j,k) = x/dble(np)
          else
            xcycle(j,k) = absent
          endif

        enddo
      enddo

c add leap day
      if((xcycle(2,28).gt.absenr).and.(xcycle(3,1).gt.absenr)) then
        xcycle(2,29) = (xcycle(2,28) + xcycle(3,1))/2.0d0
      else
        xcycle(2,29) = absentr
      endif

      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine smoothcycle_poly(xcycle)
c this routine smooth the annual cycle by fitting a polynomial approximation
c to the annual cycle
c
      implicit none
      include 'comgeneral.h'

      integer       ifail,kplus1,kp1max,nrows,mmax
      parameter     (kp1max=25,nrows=kp1max,mmax=12*31*3)
      real*8        xcycle(12,31)
      real*8        xdum(mmax),ydum(mmax)
      real*8        a(nrows,kp1max),ak(kp1max),s(kp1max),w(mmax)
      real*8        work1(3*mmax),work2(2*kp1max)
      integer       i,j,k,absen,nrec,np,length,nrecmax
      real*8        absenr,absentr,x,fit
      external      e02adf,e02aef

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

      nrec = 0
      np = 0

c make an array of 14 months
      do k=1,31
        nrec = nrec + 1
        if(xcycle(12,k).gt.absenr) then
          np=np+1
          xdum(np) = dble(nrec)
          ydum(np) = xcycle(12,k)
          w(np) = 1.0
        endif
      enddo
      do j=1,12
        call lengthofmonth(1973,j,length)
        do k=1,length
          nrec = nrec + 1
          if(xcycle(j,k).gt.absenr) then
            np=np+1
            xdum(np) = dble(nrec)
            ydum(np) = xcycle(j,k)
            w(np) = 1.0
          endif
        enddo
      enddo
      do k=1,31
        nrec = nrec + 1
        if(xcycle(1,k).gt.absenr) then
          np=np+1
          xdum(np) = dble(nrec)
          ydum(np) = xcycle(1,k)
          w(np) = 1.0
        endif
      enddo

      nrecmax = nrec

      if(np.ge.25) then
c this makes the fit with Chebyshev polynomials
        kplus1=25
        ifail = 0
        call e02adf(np,kplus1,nrows,xdum,ydum,w,work1,work2,a,s,ifail)

        do j=1,kplus1
          ak(j) = a(kplus1,j)
        enddo

c this calculates the smoothed cycled by evaluating the fit
        nrec = 31
        do j=1,12
          call lengthofmonth(1973,j,length)
          do k=1,length
            nrec = nrec + 1
            x = dble((nrec-1)-(nrecmax-nrec))/dble(nrecmax-1)
            ifail = 0
            call e02aef(kplus1,ak,x,fit,ifail)
            xcycle(j,k) = fit
          enddo
        enddo

c add leap year
        xcycle(2,29) = (xcycle(2,28) + xcycle(3,1))/2.0d0
      endif

      return
      end
c
c---------------------------------------------------------------------------
c
        SUBROUTINE SPIPE3(nrun,pp,par1,par2,par3,pzero,spi,probne,
     1  pcpacc,xmom1,xmom2,xmom3,amssng,num,numpos,maxyrs,index,y,
     2  x,temparr)
c    
c       INPUT:  pp                one dimensional array of input data
c                                 where input is monthly precip data 
c                                 beginning with a January value
c               nrun              length of moving totals
c               amssng            value for missing data
c               maxyrs            maximum possible number of years of
c                                 data
c
c       OUTPUT: par1,par2,par3    parameters (mean, s.d., skewness) of
c                                 Pearson Type III distribution used to
c                                 fit the input data > 0
c               xmom1,xmom2,xmom3 l-moments of the input data > 0
c               pzero             relative  frequency of input data = 0
c               spi               standardized precip index
c               probne            probability of observed data <= x
c               pcpacc            data accumulated over length nrun
c               num               number of non-missing observations
c               numpos            number of non-zero observations
c
c       This program used the program of John Kleist, Colorado State
c       Univ., that computed spi values from a 2 parameter ML gamma
c       probability distribution, as a basis for the moving totals, etc.
c       code.  I added the PE3 l-moment computations.  The subroutines
c       and functions in upper case came from Jon Hosking, IBM.
c               Ned Guttman
c               Nov. 1997
c

        real*8    pp(maxyrs*12),probne(maxyrs*12),pcpacc(maxyrs*12),
     1            pzero(12),spi(maxyrs*12)
        integer   num(12),numpos(12)
        real*8    par1(12),par2(12),par3(12),index(maxyrs*12),
     1            xmom1(12),xmom2(12),xmom3(12),xmom(3),para(3),
     2            y(maxyrs*12),x(maxyrs),temparr(maxyrs+1),
     3            cdfpe3,derf,dlgamma,gamind,quastn,amssng,zero
c
        amssn = amssng + 1.0d0
c
c
c       save par1,par2,par3 to 8 decimal places if you want to use
c       them later for any additional calculations...not keeping
c       enough decimal places leads to errors in probabilities!     
c
c
c       the first nrun-1 index values will be missing
c
        do j=1,nrun-1   
          index(j)=amssng
          probne(j)=amssng
          pcpacc(j)=amssng
        enddo
c
c       sum nrun precip values; store them in appropriate index location
c
c       if any value is missing, set the sum to missing
c
        do 30 j=nrun,maxyrs*12
          index(j)=0.0
          do 20 i=0,nrun-1
            if(pp(j-i).gt.amssn)then
              index(j)=index(j)+pp(j-i)
              pcpacc(j)=index(j)
            else
              index(j)=amssng
              probne(j)=amssng
              pcpacc(j)=amssng
              goto 30
            endif
 20       continue
 30     continue
       
c
c       to maintain seasonality, do everything by month
c
        do 50 i=0,11
           n=0
           nz=0
           np=0
           if (nrun.le.12) then
             do 40 j=nrun+i,maxyrs*12,12
               if(index(j).gt.amssn) then
c
c          this routine calculates lmoments and parameters of the PE3
c          distribution for all input data.  if you want to compute the
c          parameters for a specific period of record, then the limits
c          of j must be set appropriately for temparr. it is recommended
c          that for both monitoring and historical perspective that
c          all data be used for computing the probabilities, and that
c          the historical time series be the "latest" and most complete
c          output.
c
c          n-count for all non-missing data, including zeroes
c
                 n=n+1
                 if(index(j).gt.0) then
c
c          n-count for all non-missing, non-zero data; temparr is for
c          non-zero, non-missing precipitation
c
                   np=np+1
                   temparr(np)=index(j)
                 else 
c
c          n-count for all non-missing, zero data
c
                   nz=nz+1
                 endif
               endif
 40          continue
        elseif (nrun.gt.12) then
          j=nrun+i 
c
c         step the sequence of data by the length of nrun to get
c         independent samples 
c
          if (nrun.gt.12.and.nrun.le.24) nstep=24
          if (nrun.gt.24.and.nrun.le.36) nstep=36
          if (nrun.gt.36.and.nrun.le.48) nstep=48
          if (nrun.gt.48.and.nrun.le.60) nstep=60
          if (nrun.gt.60.and.nrun.le.72) nstep=72
 41       if (j.gt.maxyrs*12) GOTO 42
c
c          this routine calculates lmoments and parameters of the PE3
c          distribution for all input data.  if you want to compute the
c          parameters for a specific period of record, then the limits
c          of j must be set appropriately for temparr
c
             if(index(j).gt.amssn) then
               n=n+1
               if(index(j).gt.0) then
                 np=np+1
                 temparr(np)=index(j)
               else 
                 nz=nz+1
               endif
               j=j+nstep
               goto 41
             else
c
c          look for next non-missing nrun>12 value
c
               j=j+12
               goto 41
             endif
           endif
 42     im=mod(nrun+i-1,12)+1
        pzero(im)=dble(nz)/dble(n)
        num(im)=n
        numpos(im)=np
c
c       order the data
c
        call sort_hosking(temparr,x,np,maxyrs)
c   
c       fit l-moments for non-zero precip; if 2nd moment=0, routine
c       fails and all 3 moments set to zero (condition indicates all
c       data values are equal).  also, if not enough non-zero data
c       (less than 3), routine fails and all moments set to zero.
c
        zero=0.0D0
        call samlmr(x,np,xmom,3,zero,zero,ifail) 
        xmom1(im)=xmom(1)
        xmom2(im)=xmom(2)
        xmom3(im)=xmom(3)
c       
c       compute parameters of Pearson type III (3-parameter gamma); the
c       three parameters are the mean, s.d., skewness of the NON-ZERO data.
c       The mean for all the data is (1-pzero)*par1 or (1-pzero)*xmom1
c       since par1=xmom1
c
        call pelpe3(xmom,para,ifail)
c
c       if ifail not = 0, then all 3 parameters are set to zero
c       if routine fails, output is sent to unit 6 by pelpe3 routine
c
        par1(im)=para(1)
        par2(im)=para(2)
        par3(im)=para(3)
 50     continue
c
c       compute the probability (cdfpe3), take into account the 
c       mixed distribution if pzero not equal to zero, truncate
c       the probability from .001 to .999, transform the probability to
c       spi (quastn). ifail=1 indicates invalid parameter 2 (negative).
c       if routines fail, output is sent to unit 6 by called routines.
c
        do 60 j=nrun,maxyrs*12
          im=mod(j-1,12)+1
          para(1)=par1(im)
          para(2)=par2(im)
          para(3)=par3(im)
c
c         set missing values to amssng
c
 61       if(index(j).le.amssn) then
             probne(j)=amssng
             spi(j)=amssng
             goto 60
          else
             index(j)=cdfpe3(index(j),para,ifail,amssng)
c
c         if cdf routine fails, set probne and spi = amssng
c
             if(ifail.ne.0) then
               index(j)=amssng
               goto 61
             endif
             y(j)=pzero(im)+(1-pzero(im))*index(j)
c
c            force the probabilities and therefore the spi to be bounded
c            by +/- 3.09
c
c            write(6,*) j,y(j),index(j)
             if(y(j).gt..999)y(j)=.999
             if(y(j).lt..001)y(J)=.001
             probne(j)=y(j)
             index(j)=quastn(y(j))
             spi(j)=index(j)
          endif
 60    continue
       return
       end
c
c
      DOUBLE PRECISION FUNCTION CDFPE3(X,PARA,ifail,amssng)
C***********************************************************************
C*                                                                     *
C*  FORTRAN CODE WRITTEN FOR INCLUSION IN IBM RESEARCH REPORT RC20525, *
C*  'FORTRAN ROUTINES FOR USE WITH THE METHOD OF L-MOMENTS, VERSION 3' *
C*                                                                     *
C*  J. R. M. HOSKING                                                   *
C*  IBM RESEARCH DIVISION                                              *
C*  T. J. WATSON RESEARCH CENTER                                       *
C*  YORKTOWN HEIGHTS                                                   *
C*  NEW YORK 10598, U.S.A.                                             *
C*                                                                     *
C*  VERSION 3     AUGUST 1996                                          *
C*                                                                     *
C***********************************************************************
C
C  DISTRIBUTION FUNCTION OF THE PEARSON TYPE 3 DISTRIBUTION
C
C  OTHER ROUTINES USED: DERF,DLGAMA,GAMIND
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION PARA(3)
      DATA ZERO/0D0/,HALF/0.5D0/,ONE/1D0/,TWO/2D0/,FOUR/4D0/
      DATA RTHALF/0.70710 67811 86547 524D0/
C
C         SMALL IS USED TO TEST WHETHER SKEWNESS IS EFFECTIVELY ZERO
C
      DATA SMALL/1D-6/
C
      ifail=0
      CDFPE3=ZERO
      IF(PARA(2).LE.ZERO)GOTO 1000
      GAMMA=PARA(3)
      IF(DABS(GAMMA).LE.SMALL)GOTO 10
      ALPHA=FOUR/(GAMMA*GAMMA)
      Z=TWO*(X-PARA(1))/(PARA(2)*GAMMA)+ALPHA
      IF(Z.GT.ZERO)CDFPE3=GAMIND(Z,ALPHA,DLGAMA(ALPHA))
      IF(GAMMA.LT.ZERO)CDFPE3=ONE-CDFPE3
      RETURN
C
C         ZERO SKEWNESS
C
   10 Z=(X-PARA(1))/PARA(2)
      CDFPE3=HALF+HALF*DERF(Z*RTHALF)
      RETURN
C
c1000 WRITE(6,7000)
1000  ifail=1
      cdfpe3=amssng
      RETURN
C
 7000 FORMAT(' *** ERROR *** ROUTINE CDFPE3 : PARAMETERS INVALID')
      END
c
c
      SUBROUTINE PELPE3(XMOM,PARA,ifail)
C***********************************************************************
C*                                                                     *
C*  FORTRAN CODE WRITTEN FOR INCLUSION IN IBM RESEARCH REPORT RC20525, *
C*  'FORTRAN ROUTINES FOR USE WITH THE METHOD OF L-MOMENTS, VERSION 3' *
C*                                                                     *
C*  J. R. M. HOSKING                                                   *
C*  IBM RESEARCH DIVISION                                              *
C*  T. J. WATSON RESEARCH CENTER                                       *
C*  YORKTOWN HEIGHTS                                                   *
C*  NEW YORK 10598, U.S.A.                                             *
C*                                                                     *
C*  VERSION 3     AUGUST 1996                                          *
C*                                                                     *
C***********************************************************************
C
C  PARAMETER ESTIMATION VIA L-MOMENTS FOR THE PEARSON TYPE 3 DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  XMOM   * INPUT* ARRAY OF LENGTH 3. CONTAINS THE L-MOMENTS LAMBDA-1,
C                  LAMBDA-2 AND TAU-3.
C  PARA   *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS THE PARAMETERS
C                  IN THE ORDER MU, SIGMA, GAMMA (MEAN, S.D., SKEWNESS).
C
C  OTHER ROUTINES USED: DLGAMA
C
C  METHOD: RATIONAL APPROXIMATION IS USED TO EXPRESS ALPHA, THE SHAPE
C  PARAMETER OF THE GAMMA DISTRIBUTION, AS A FUNCTION OF TAU-3.
C  RELATIVE ACCURACY OF THE APPROXIMATION IS BETTER THAN 3E-5.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION XMOM(3),PARA(3)
      DATA ZERO/0D0/,THIRD/0.33333333D0/,HALF/0.5D0/,ONE/1D0/,TWO/2D0/
C
C         SMALL IS USED TO TEST WHETHER SKEWNESS IS EFFECTIVELY ZERO
C
      DATA SMALL/1D-6/
C
C         CONSTANTS USED IN MINIMAX APPROXIMATIONS
C
      DATA C1,C2,C3/ 0.2906D0,  0.1882D0,  0.0442D0/
      DATA D1,D2,D3/ 0.36067D0,-0.59567D0, 0.25361D0/
      DATA D4,D5,D6/-2.78861D0, 2.56096D0,-0.77045D0/
      DATA PI3,ROOTPI/9.4247780D0,1.7724539D0/
C
      ifail=0
      T3=DABS(XMOM(3))
      IF(XMOM(2).LE.ZERO.OR.T3.GE.ONE)GOTO 1000
      IF(T3.LE.SMALL)GOTO 100
      IF(T3.GE.THIRD)GOTO 10
      T=PI3*T3*T3
      ALPHA=(ONE+C1*T)/(T*(ONE+T*(C2+T*C3)))
      GOTO 20
   10 CONTINUE
      T=ONE-T3
      ALPHA=T*(D1+T*(D2+T*D3))/(ONE+T*(D4+T*(D5+T*D6)))
   20 CONTINUE
      RTALPH=DSQRT(ALPHA)
      BETA=ROOTPI*XMOM(2)*DEXP(DLGAMA(ALPHA)-DLGAMA(ALPHA+HALF))
      PARA(1)=XMOM(1)
      PARA(2)=BETA*RTALPH
      PARA(3)=TWO/RTALPH
      IF(XMOM(3).LT.ZERO)PARA(3)=-PARA(3)
      RETURN
C
C         ZERO SKEWNESS
C
  100 CONTINUE
      PARA(1)=XMOM(1)
      PARA(2)=XMOM(2)*ROOTPI
      PARA(3)=ZERO
      RETURN
C
c1000 WRITE(6,7000)
 1000 DO 1010 I=1,3
 1010 PARA(I)=ZERO
      ifail=1
      RETURN
C
 7000 FORMAT(' *** ERROR *** ROUTINE PELPE3 : L-MOMENTS INVALID')
      END
c
c
c
      SUBROUTINE SAMLMR(X,N,XMOM,NMOM,A,B,ifail)
C***********************************************************************
C*                                                                     *
C*  FORTRAN CODE WRITTEN FOR INCLUSION IN IBM RESEARCH REPORT RC20525, *
C*  'FORTRAN ROUTINES FOR USE WITH THE METHOD OF L-MOMENTS, VERSION 3' *
C*                                                                     *
C*  J. R. M. HOSKING                                                   *
C*  IBM RESEARCH DIVISION                                              *
C*  T. J. WATSON RESEARCH CENTER                                       *
C*  YORKTOWN HEIGHTS                                                   *
C*  NEW YORK 10598, U.S.A.                                             *
C*                                                                     *
C*  VERSION 3     AUGUST 1996                                          *
C*                                                                     *
C***********************************************************************
C
C  SAMPLE L-MOMENTS OF A DATA ARRAY
C
C  PARAMETERS OF ROUTINE:
C  X      * INPUT* ARRAY OF LENGTH N. CONTAINS THE DATA, IN ASCENDING
C                  ORDER.
C  N      * INPUT* NUMBER OF DATA VALUES
C  XMOM   *OUTPUT* ARRAY OF LENGTH NMOM. ON EXIT, CONTAINS THE SAMPLE
C                  L-MOMENTS L-1, L-2, T-3, T-4, ... .
C  NMOM   * INPUT* NUMBER OF L-MOMENTS TO BE FOUND. AT MOST MAX(N,20).
C  A      * INPUT* ) PARAMETERS OF PLOTTING
C  B      * INPUT* ) POSITION (SEE BELOW)
C
C  FOR UNBIASED ESTIMATES (OF THE LAMBDA'S) SET A=B=ZERO. OTHERWISE,
C  PLOTTING-POSITION ESTIMATORS ARE USED, BASED ON THE PLOTTING POSITION
C  (J+A)/(N+B)  FOR THE J'TH SMALLEST OF N OBSERVATIONS. FOR EXAMPLE,
C  A=-0.35D0 AND B=0.0D0 YIELDS THE ESTIMATORS RECOMMENDED BY
C  HOSKING ET AL. (1985, TECHNOMETRICS) FOR THE GEV DISTRIBUTION.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(N),XMOM(NMOM),SUM(20)
      DATA ZERO/0D0/,ONE/1D0/
      ifail=0
      IF(NMOM.GT.20.OR.NMOM.GT.N)GOTO 1000
      DO 10 J=1,NMOM
   10 SUM(J)=ZERO
      IF(A.EQ.ZERO.AND.B.EQ.ZERO)GOTO 50
      IF(A.LE.-ONE.OR.A.GE.B)GOTO 1010
C
C         PLOTTING-POSITION ESTIMATES OF PWM'S
C
      DO 30 I=1,N
      PPOS=(I+A)/(N+B)
      TERM=X(I)
      SUM(1)=SUM(1)+TERM
      DO 20 J=2,NMOM
      TERM=TERM*PPOS
   20 SUM(J)=SUM(J)+TERM
   30 CONTINUE
      DO 40 J=1,NMOM
   40 SUM(J)=SUM(J)/N
      GOTO 100
C
C         UNBIASED ESTIMATES OF PWM'S
C
   50 DO 70 I=1,N
      Z=I
      TERM=X(I)
      SUM(1)=SUM(1)+TERM
      DO 60 J=2,NMOM
      Z=Z-ONE
      TERM=TERM*Z
   60 SUM(J)=SUM(J)+TERM
   70 CONTINUE
      Y=N
      Z=N
      SUM(1)=SUM(1)/Z
      DO 80 J=2,NMOM
      Y=Y-ONE
      Z=Z*Y
   80 SUM(J)=SUM(J)/Z
C
C         L-MOMENTS
C
  100 K=NMOM
      P0=ONE
      IF(NMOM-NMOM/2*2.EQ.1)P0=-ONE
      DO 120 KK=2,NMOM
      AK=K
      P0=-P0
      P=P0
      TEMP=P*SUM(1)
      DO 110 I=1,K-1
      AI=I
      P=-P*(AK+AI-ONE)*(AK-AI)/(AI*AI)
  110 TEMP=TEMP+P*SUM(I+1)
      SUM(K)=TEMP
  120 K=K-1
      XMOM(1)=SUM(1)
      IF(NMOM.EQ.1)RETURN
      XMOM(2)=SUM(2)
      IF(SUM(2).EQ.ZERO)GOTO 1020
      IF(NMOM.EQ.2)RETURN
      DO 130 K=3,NMOM
  130 XMOM(K)=SUM(K)/SUM(2)
      RETURN
C
 1000 ifail=1
      do i=1,3
        xmom(i)=0
      enddo
c     WRITE(6,7000)
      RETURN
c1010 WRITE(6,7010)
 1010 continue
      RETURN
 1020 ifail=1
      do i=1,3
        xmom(i)=0
      enddo
c     WRITE(6,7020)
      RETURN
C
 7000 FORMAT(' *** ERROR *** ROUTINE SAMLMR : PARAMETER NMOM INVALID')
 7010 FORMAT(' *** ERROR *** ROUTINE SAMLMR :',
     *  ' PLOTTING-POSITION PARAMETERS INVALID')
 7020 FORMAT(' *** ERROR *** ROUTINE SAMLMR : ALL DATA VALUES EQUAL')
      END
c
c
c
      DOUBLE PRECISION FUNCTION DERF(X)
C***********************************************************************
C*                                                                     *
C*  FORTRAN CODE WRITTEN FOR INCLUSION IN IBM RESEARCH REPORT RC20525, *
C*  'FORTRAN ROUTINES FOR USE WITH THE METHOD OF L-MOMENTS, VERSION 3' *
C*                                                                     *
C*  J. R. M. HOSKING                                                   *
C*  IBM RESEARCH DIVISION                                              *
C*  T. J. WATSON RESEARCH CENTER                                       *
C*  YORKTOWN HEIGHTS                                                   *
C*  NEW YORK 10598, U.S.A.                                             *
C*                                                                     *
C*  VERSION 3     AUGUST 1996                                          *
C*                                                                     *
C***********************************************************************
C
C  ERROR FUNCTION
C
C  BASED ON ALGORITHM 5666, J.F.HART ET AL. (1968) 'COMPUTER
C  APPROXIMATIONS'
C
C  ACCURATE TO 15 DECIMAL PLACES
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DATA ZERO/0D0/,ONE/1D0/,TWO/2D0/,THREE/3D0/,FOUR/4D0/,P65/0.65D0/
C
C         COEFFICIENTS OF RATIONAL-FUNCTION APPROXIMATION
C
      DATA P0,P1,P2,P3,P4,P5,P6/
     *  0.22020 68679 12376 1D3,    0.22121 35961 69931 1D3,
     *  0.11207 92914 97870 9D3,    0.33912 86607 83830 0D2,
     *  0.63739 62203 53165 0D1,    0.70038 30644 43688 1D0,
     *  0.35262 49659 98910 9D-1/
      DATA Q0,Q1,Q2,Q3,Q4,Q5,Q6,Q7/
     *  0.44041 37358 24752 2D3,   0.79382 65125 19948 4D3,
     *  0.63733 36333 78831 1D3,   0.29656 42487 79673 7D3,
     *  0.86780 73220 29460 8D2,   0.16064 17757 92069 5D2,
     *  0.17556 67163 18264 2D1,   0.88388 34764 83184 4D-1/
C
C         C1 IS SQRT(2), C2 IS SQRT(2/PI)
C         BIG IS THE POINT AT WHICH DERF=1 TO MACHINE PRECISION
C
      DATA C1/1.4142 13562 37309 5D0/
      DATA C2/7.9788 45608 02865 4D-1/
      DATA BIG/6.25D0/,CRIT/5D0/
C
      DERF=ZERO
      IF(X.EQ.ZERO)RETURN
      XX=DABS(X)
      IF(XX.GT.BIG)GOTO 20
      EXPNTL=DEXP(-X*X)
      ZZ=DABS(X*C1)
      IF(XX.GT.CRIT)GOTO 10
      DERF=EXPNTL*((((((P6*ZZ+P5)*ZZ+P4)*ZZ+P3)*ZZ+P2)*ZZ+P1)*ZZ+P0)/
     *  (((((((Q7*ZZ+Q6)*ZZ+Q5)*ZZ+Q4)*ZZ+Q3)*ZZ+Q2)*ZZ+Q1)*ZZ+Q0)
      IF(X.GT.ZERO)DERF=ONE-TWO*DERF
      IF(X.LT.ZERO)DERF=TWO*DERF-ONE
      RETURN
C
   10 DERF=EXPNTL*C2/(ZZ+ONE/(ZZ+TWO/(ZZ+THREE/(ZZ+FOUR/(ZZ+P65)))))
      IF(X.GT.ZERO)DERF=ONE-DERF
      IF(X.LT.ZERO)DERF=DERF-ONE
      RETURN
C
   20 DERF=ONE
      IF(X.LT.ZERO)DERF=-ONE
      RETURN
      END
c
c
c
      DOUBLE PRECISION FUNCTION DLGAMA(X)
C***********************************************************************
C*                                                                     *
C*  FORTRAN CODE WRITTEN FOR INCLUSION IN IBM RESEARCH REPORT RC20525, *
C*  'FORTRAN ROUTINES FOR USE WITH THE METHOD OF L-MOMENTS, VERSION 3' *
C*                                                                     *
C*  J. R. M. HOSKING                                                   *
C*  IBM RESEARCH DIVISION                                              *
C*  T. J. WATSON RESEARCH CENTER                                       *
C*  YORKTOWN HEIGHTS                                                   *
C*  NEW YORK 10598, U.S.A.                                             *
C*                                                                     *
C*  VERSION 3     AUGUST 1996                                          *
C*                                                                     *
C***********************************************************************
C
C  LOGARITHM OF GAMMA FUNCTION
C
C  BASED ON ALGORITHM ACM291, COMMUN. ASSOC. COMPUT. MACH. (1966)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA SMALL,CRIT,BIG,TOOBIG/1D-7,13D0,1D9,2D36/
C
C         C0 IS 0.5*LOG(2*PI)
C         C1...C7 ARE THE COEFFTS OF THE ASYMPTOTIC EXPANSION OF DLGAMA
C
      DATA C0,C1,C2,C3,C4,C5,C6,C7/
     *   0.91893 85332 04672 742D 0,  0.83333 33333 33333 333D-1,
     *  -0.27777 77777 77777 778D-2,  0.79365 07936 50793 651D-3,
     *  -0.59523 80952 38095 238D-3,  0.84175 08417 50841 751D-3,
     *  -0.19175 26917 52691 753D-2,  0.64102 56410 25641 026D-2/
C
C         S1 IS -(EULER'S CONSTANT), S2 IS PI**2/12
C
      DATA S1/-0.57721 56649 01532 861D 0/
      DATA S2/ 0.82246 70334 24113 218D 0/
C
      DATA ZERO/0D0/,HALF/0.5D0/,ONE/1D0/,TWO/2D0/
      DLGAMA=ZERO
      IF(X.LE.ZERO)GOTO 1000
      IF(X.GT.TOOBIG)GOTO 1000
C
C         USE SMALL-X APPROXIMATION IF X IS NEAR 0, 1 OR 2
C
      IF(DABS(X-TWO).GT.SMALL)GOTO 10
      DLGAMA=DLOG(X-ONE)
      XX=X-TWO
      GOTO 20
   10 IF(DABS(X-ONE).GT.SMALL)GOTO 30
      XX=X-ONE
   20 DLGAMA=DLGAMA+XX*(S1+XX*S2)
      RETURN
   30 IF(X.GT.SMALL)GOTO 40
      DLGAMA=-DLOG(X)+S1*X
      RETURN
C
C         REDUCE TO DLGAMA(X+N) WHERE X+N.GE.CRIT
C
   40 SUM1=ZERO
      Y=X
      IF(Y.GE.CRIT)GOTO 60
      Z=ONE
   50 Z=Z*Y
      Y=Y+ONE
      IF(Y.LT.CRIT)GOTO 50
      SUM1=SUM1-DLOG(Z)
C
C         USE ASYMPTOTIC EXPANSION IF Y.GE.CRIT
C
   60 SUM1=SUM1+(Y-HALF)*DLOG(Y)-Y+C0
      SUM2=ZERO
      IF(Y.GE.BIG)GOTO 70
      Z=ONE/(Y*Y)
      SUM2=((((((C7*Z+C6)*Z+C5)*Z+C4)*Z+C3)*Z+C2)*Z+C1)/Y
   70 DLGAMA=SUM1+SUM2
      RETURN
C
 1000 WRITE(6,7000)X
      RETURN
C
 7000 FORMAT(' *** ERROR *** ROUTINE DLGAMA :',
     *  ' ARGUMENT OUT OF RANGE :',D24.16)
      END
c
c
c
      DOUBLE PRECISION FUNCTION GAMIND(X,ALPHA,G)
C***********************************************************************
C*                                                                     *
C*  FORTRAN CODE WRITTEN FOR INCLUSION IN IBM RESEARCH REPORT RC20525, *
C*  'FORTRAN ROUTINES FOR USE WITH THE METHOD OF L-MOMENTS, VERSION 3' *
C*                                                                     *
C*  J. R. M. HOSKING                                                   *
C*  IBM RESEARCH DIVISION                                              *
C*  T. J. WATSON RESEARCH CENTER                                       *
C*  YORKTOWN HEIGHTS                                                   *
C*  NEW YORK 10598, U.S.A.                                             *
C*                                                                     *
C*  VERSION 3     AUGUST 1996                                          *
C*                                                                     *
C***********************************************************************
C
C  THE INCOMPLETE GAMMA INTEGRAL
C
C  BASED ON ALGORITHM AS239, APPL. STATIST. (1988) VOL.37 NO.3
C
C  PARAMETERS OF ROUTINE:
C  X      * INPUT* ARGUMENT OF FUNCTION (UPPER LIMIT OF INTEGRATION)
C  ALPHA  * INPUT* SHAPE PARAMETER
C  G      * INPUT* LOG(GAMMA(ALPHA)). MUST BE SUPPLIED BY THE PROGRAM,
C                  E.G. AS DLGAMA(ALPHA).
C
C  OTHER ROUTINES USED: DERF
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA ZERO/0D0/,HALF/0.5D0/,ONE/1D0/,TWO/2D0/,THREE/3D0/,X13/13D0/,
     *  X36/36D0/,X42/42D0/,X119/119D0/,X1620/1620D0/,X38880/38880D0/,
     *  RTHALF/0.70710 67811 86547 524D0/
C
C         EPS,MAXIT CONTROL THE TEST FOR CONVERGENCE OF THE SERIES AND
C           CONTINUED-FRACTION EXPANSIONS.
C         OFL IS A LARGE NUMBER, USED TO RESCALE THE CONTINUED FRACTION.
C         UFL IS SUCH THAT EXP(UFL) IS JUST .GT. ZERO.
C         AHILL CONTROLS THE SWITCH TO HILL'S APPROXIMATION.
C
      DATA EPS/1D-12/,MAXIT/100000/,OFL/1D30/,UFL/-180D0/,AHILL/1D4/
      GAMIND=ZERO
      IF(ALPHA.LE.ZERO)GOTO 1000
      IF(X.LT.ZERO)GOTO 1010
      IF(X.EQ.ZERO)RETURN
C
      IF(ALPHA.GT.AHILL)GOTO 100
      IF(X.GT.ONE.AND.X.GE.ALPHA)GOTO 50
C
C         SERIES EXPANSION
C
      SUM=ONE
      TERM=ONE
      A=ALPHA
      DO 10 IT=1,MAXIT
      A=A+ONE
      TERM=TERM*X/A
      SUM=SUM+TERM
      IF(TERM.LE.EPS)GOTO 20
   10 CONTINUE
      WRITE(6,7020)
   20 ARG=ALPHA*DLOG(X)-X-G+DLOG(SUM/ALPHA)
      GAMIND=ZERO
      IF(ARG.GE.UFL)GAMIND=DEXP(ARG)
      RETURN
C
C         CONTINUED-FRACTION EXPANSION
C
   50 CONTINUE
      A=ONE-ALPHA
      B=A+X+ONE
      TERM=ZERO
      PN1=ONE
      PN2=X
      PN3=X+ONE
      PN4=X*B
      RATIO=PN3/PN4
      DO 70 IT=1,MAXIT
      A=A+ONE
      B=B+TWO
      TERM=TERM+ONE
      AN=A*TERM
      PN5=B*PN3-AN*PN1
      PN6=B*PN4-AN*PN2
      IF(PN6.EQ.ZERO)GOTO 60
      RN=PN5/PN6
      DIFF=DABS(RATIO-RN)
      IF(DIFF.LE.EPS.AND.DIFF.LE.EPS*RN)GOTO 80
      RATIO=RN
   60 PN1=PN3
      PN2=PN4
      PN3=PN5
      PN4=PN6
      IF(DABS(PN5).LT.OFL)GOTO 70
      PN1=PN1/OFL
      PN2=PN2/OFL
      PN3=PN3/OFL
      PN4=PN4/OFL
   70 CONTINUE
      WRITE(6,7020)
   80 ARG=ALPHA*DLOG(X)-X-G+DLOG(RATIO)
      GAMIND=ONE
      IF(ARG.GE.UFL)GAMIND=ONE-DEXP(ARG)
      RETURN
C
C         ALPHA IS LARGE: USE HILL'S APPROXIMATION (N.L. JOHNSON AND
C         S. KOTZ, 1970, 'CONTINUOUS UNIVARIATE DISTRIBUTIONS 1', P.180)
C
C         THE 'DO 110' LOOP CALCULATES 2*(X-ALPHA-ALPHA*DLOG(X/ALPHA)),
C         USING POWER-SERIES EXPANSION TO AVOID ROUNDING ERROR
C
  100 CONTINUE
      R=ONE/DSQRT(ALPHA)
      Z=(X-ALPHA)*R
      TERM=Z*Z
      SUM=HALF*TERM
      DO 110 I=1,12
      TERM=-TERM*Z*R
      SUM=SUM+TERM/(I+TWO)
      IF(DABS(TERM).LT.EPS)GOTO 120
  110 CONTINUE
  120 WW=TWO*SUM
      W=DSQRT(WW)
      IF(X.LT.ALPHA)W=-W
      H1=ONE/THREE
      H2=-W/X36
      H3=(-WW+X13)/X1620
      H4=(X42*WW+X119)*W/X38880
      Z=(((H4*R+H3)*R+H2)*R+H1)*R+W
      GAMIND=HALF+HALF*DERF(Z*RTHALF)
      RETURN
C
 1000 WRITE(6,7000)ALPHA
      RETURN
 1010 WRITE(6,7010)X
      RETURN
C
 7000 FORMAT(' *** ERROR *** ROUTINE GAMIND :',
     *  ' SHAPE PARAMETER OUT OF RANGE :',D16.8)
 7010 FORMAT(' *** ERROR *** ROUTINE GAMIND :',
     *  ' ARGUMENT OF FUNCTION OUT OF RANGE :',D16.8)
 7020 FORMAT(' ** WARNING ** ROUTINE GAMIND :',
     *  ' ITERATION HAS NOT CONVERGED. RESULT MAY BE UNRELIABLE.')
      END
c
c
c
      DOUBLE PRECISION FUNCTION QUASTN(F)
C***********************************************************************
C*                                                                     *
C*  FORTRAN CODE WRITTEN FOR INCLUSION IN IBM RESEARCH REPORT RC20525, *
C*  'FORTRAN ROUTINES FOR USE WITH THE METHOD OF L-MOMENTS, VERSION 3' *
C*                                                                     *
C*  J. R. M. HOSKING                                                   *
C*  IBM RESEARCH DIVISION                                              *
C*  T. J. WATSON RESEARCH CENTER                                       *
C*  YORKTOWN HEIGHTS                                                   *
C*  NEW YORK 10598, U.S.A.                                             *
C*                                                                     *
C*  VERSION 3     AUGUST 1996                                          *
C*                                                                     *
C***********************************************************************
C
C  QUANTILE FUNCTION OF THE STANDARD NORMAL DISTRIBUTION
C
C  BASED ON ALGORITHM AS241, APPL. STATIST. (1988) VOL.37 NO.3
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA ZERO/0D0/,HALF/0.5D0/,ONE/1D0/
      DATA SPLIT1/0.425D0/,SPLIT2/5D0/,CONST1/0.180625D0/,CONST2/1.6D0/
C
C         COEFFICIENTS OF RATIONAL-FUNCTION APPROXIMATIONS
C
      DATA A0,A1,A2,A3,A4,A5,A6,A7,B1,B2,B3,B4,B5,B6,B7/
     *                                0.33871 32872 79636 661D  1,
     *  0.13314 16678 91784 377D  3,  0.19715 90950 30655 144D  4,
     *  0.13731 69376 55094 611D  5,  0.45921 95393 15498 715D  5,
     *  0.67265 77092 70087 009D  5,  0.33430 57558 35881 281D  5,
     *  0.25090 80928 73012 267D  4,  0.42313 33070 16009 113D  2,
     *  0.68718 70074 92057 908D  3,  0.53941 96021 42475 111D  4,
     *  0.21213 79430 15865 959D  5,  0.39307 89580 00927 106D  5,
     *  0.28729 08573 57219 427D  5,  0.52264 95278 85285 456D  4/
      DATA C0,C1,C2,C3,C4,C5,C6,C7,D1,D2,D3,D4,D5,D6,D7/
     *                                0.14234 37110 74968 358D  1,
     *  0.46303 37846 15654 530D  1,  0.57694 97221 46069 141D  1,
     *  0.36478 48324 76320 461D  1,  0.12704 58252 45236 838D  1,
     *  0.24178 07251 77450 612D  0,  0.22723 84498 92691 846D -1,
     *  0.77454 50142 78341 408D -3,  0.20531 91626 63775 882D  1,
     *  0.16763 84830 18380 385D  1,  0.68976 73349 85100 005D  0,
     *  0.14810 39764 27480 075D  0,  0.15198 66656 36164 572D -1,
     *  0.54759 38084 99534 495D -3,  0.10507 50071 64441 684D -8/
      DATA E0,E1,E2,E3,E4,E5,E6,E7,F1,F2,F3,F4,F5,F6,F7/
     *                                0.66579 04643 50110 378D  1,
     *  0.54637 84911 16411 437D  1,  0.17848 26539 91729 133D  1,
     *  0.29656 05718 28504 891D  0,  0.26532 18952 65761 230D -1,
     *  0.12426 60947 38807 844D -2,  0.27115 55568 74348 758D -4,
     *  0.20103 34399 29228 813D -6,  0.59983 22065 55887 938D  0,
     *  0.13692 98809 22735 805D  0,  0.14875 36129 08506 149D -1,
     *  0.78686 91311 45613 259D -3,  0.18463 18317 51005 468D -4,
     *  0.14215 11758 31644 589D -6,  0.20442 63103 38993 979D-14/
C
      Q=F-HALF
      IF(DABS(Q).GT.SPLIT1)GOTO 10
      R=CONST1-Q*Q
      QUASTN=Q*(((((((A7*R+A6)*R+A5)*R+A4)*R+A3)*R+A2)*R+A1)*R+A0)
     *        /(((((((B7*R+B6)*R+B5)*R+B4)*R+B3)*R+B2)*R+B1)*R+ONE)
      RETURN
   10 R=F
      IF(Q.GE.ZERO)R=ONE-F
      IF(R.LE.ZERO)GOTO 1000
      R=DSQRT(-DLOG(R))
      IF(R.GT.SPLIT2)GOTO 20
      R=R-CONST2
      QUASTN=(((((((C7*R+C6)*R+C5)*R+C4)*R+C3)*R+C2)*R+C1)*R+C0)
     *      /(((((((D7*R+D6)*R+D5)*R+D4)*R+D3)*R+D2)*R+D1)*R+ONE)
      GOTO 30
   20 R=R-SPLIT2
      QUASTN=(((((((E7*R+E6)*R+E5)*R+E4)*R+E3)*R+E2)*R+E1)*R+E0)
     *      /(((((((F7*R+F6)*R+F5)*R+F4)*R+F3)*R+F2)*R+F1)*R+ONE)
   30 IF(Q.LT.ZERO)QUASTN=-QUASTN
      RETURN
C
 1000 WRITE(6,7000)F
      QUASTN=ZERO
      RETURN
C
 7000 FORMAT(' *** ERROR *** ROUTINE QUASTN :',
     *  ' ARGUMENT OF FUNCTION INVALID')
      END
c
c
c
      SUBROUTINE SORT_HOSKING(temparr,X,N,maxyrs)
C***********************************************************************
C*                                                                     *
C*  FORTRAN CODE WRITTEN FOR INCLUSION IN IBM RESEARCH REPORT RC20525, *
C*  'FORTRAN ROUTINES FOR USE WITH THE METHOD OF L-MOMENTS, VERSION 3' *
C*                                                                     *
C*  J. R. M. HOSKING                                                   *
C*  IBM RESEARCH DIVISION                                              *
C*  T. J. WATSON RESEARCH CENTER                                       *
C*  YORKTOWN HEIGHTS                                                   *
C*  NEW YORK 10598, U.S.A.                                             *
C*                                                                     *
C*  VERSION 3     AUGUST 1996                                          *
C*                                                                     *
C***********************************************************************
C
C  SORTS THE ARRAY X INTO ASCENDING ORDER
C
C  PARAMETERS OF ROUTINE:
C  X      *IN/OUT* ARRAY OF LENGTH N. CONTAINS THE NUMBERS TO BE SORTED.
C                  ON EXIT, CONTAINS THE SORTED NUMBERS.
C  N      * INPUT* NUMBER OF ELEMENTS TO BE SORTED
C
C  METHOD USED IS SHELL SORT WITH SEQUENCE OF INCREMENTS AS IN
C  D.F.KNUTH (1969) 'THE ART OF COMPUTER PROGRAMMING', VOL.3, P.95
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(N),temparr(maxyrs+1)
      IF(N.LE.1)RETURN
      do i=1,n
        x(i)=temparr(i)
      enddo
      J=4
      DO 10 I=1,100
      J=3*J+1
      IF(J.GE.N)GOTO 20
   10 CONTINUE
   20 CONTINUE
      M=(J/3)
      DO 60 MM=1,100
      M=M/3
      IF(M.EQ.0)RETURN
      DO 50 I=M+1,N
      TEST=X(I)
      J=I
      DO 30 JJ=1,100
      J=J-M
      IF(J.LE.0)GOTO 40
      IF(TEST.GE.X(J))GOTO 40
   30 X(J+M)=X(J)
   40 CONTINUE
   50 X(J+M)=TEST
   60 CONTINUE
      END
c
c---------------------------------------------------------------------
c
      subroutine calcMaxSeason(fyear,lyear,c,b)
c this routine determines the max. value in the appropriate time span of each season 
c
      implicit none
      include 'comgeneral.h'

      real*8        b(yrbeg:yrend,nseason)
      real*8        c(yrbeg:yrend,12)
      integer       i,j,k,length,absen,nm,xm,fyear,fm,lm,lyear
      real*8        absenr,absentr,xmax

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

      do i=fyear,lyear
        do j=1,nseason
          xmax = -1.0D6
          if(j.eq.1) then
            do k=1,12
              if((c(i,k).gt.absenr).and.(c(i,k).gt.xmax)) xmax = c(i,k)
            enddo

          elseif(j.eq.2) then
            if(i.gt.yrbeg)then
              do k=10,12
                if((c(i-1,k).gt.absenr).and.(c(i-1,k).gt.xmax)) 
     +              xmax = c(i-1,k)
              enddo
            endif
            do k=1,3
              if((c(i,k).gt.absenr).and.(c(i,k).gt.xmax)) xmax = c(i,k)
            enddo

          elseif(j.eq.3) then
            do k=4,9
              if((c(i,k).gt.absenr).and.(c(i,k).gt.xmax)) xmax = c(i,k)
            enddo

          elseif(j.eq.4) then
            if(i.gt.yrbeg)then
              k=12
              if((c(i-1,k).gt.absenr).and.(c(i-1,k).gt.xmax)) 
     +             xmax = c(i-1,k)
            endif
            do k=1,2
              if((c(i,k).gt.absenr).and.(c(i,k).gt.xmax)) xmax = c(i,k)
            enddo

          elseif(j.le.7)then
            call bounds(j,fm,lm)
            do k=fm,lm
              if((c(i,k).gt.absenr).and.(c(i,k).gt.xmax)) xmax = c(i,k)
            enddo

          else
            if((c(i,j-7).gt.absenr).and.(c(i,j-7).gt.xmax)) 
     +              xmax = c(i,j-7)
          endif

          b(i,j) = xmax

        enddo
      enddo

      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine calcMinSeason(fyear,lyear,c,b)
c this routine determines the min. value in the appropriate time span of each season 
c
      implicit none
      include 'comgeneral.h'

      real*8        b(yrbeg:yrend,nseason)
      real*8        c(yrbeg:yrend,12)
      integer       i,j,k,length,absen,nm,xm,fyear,fm,lm,lyear
      real*8        absenr,absentr,xmin

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

      do i=fyear,lyear
        do j=1,nseason
          xmin = 1.0D6
          if(j.eq.1) then
            do k=1,12
              if((c(i,k).gt.absenr).and.(c(i,k).lt.xmin)) xmin = c(i,k)
            enddo

          elseif(j.eq.2) then
            if(fyear.gt.yrbeg)then
            do k=10,12
              if((c(i-1,k).gt.absenr).and.(c(i-1,k).lt.xmin)) 
     +            xmin = c(i-1,k)
            enddo
            endif
            do k=1,3
              if((c(i,k).gt.absenr).and.(c(i,k).lt.xmin)) xmin = c(i,k)
            enddo

          elseif(j.eq.3) then
            do k=4,9
              if((c(i,k).gt.absenr).and.(c(i,k).lt.xmin)) xmin = c(i,k)
            enddo

          elseif(j.eq.4) then
            if(fyear.gt.yrbeg)then
              k=12
              if((c(i-1,k).gt.absenr).and.(c(i-1,k).lt.xmin)) 
     +             xmin = c(i-1,k)
            endif
            do k=1,2
              if((c(i,k).gt.absenr).and.(c(i,k).lt.xmin)) xmin = c(i,k)
            enddo

          elseif(j.le.7)then
            call bounds(j,fm,lm)
            do k=fm,lm
              if((c(i,k).gt.absenr).and.(c(i,k).lt.xmin)) xmin = c(i,k)
            enddo

          else
            if((c(i,j-7).gt.absenr).and.(c(i,j-7).lt.xmin)) 
     +              xmin = c(i,j-7)
          endif

          b(i,j) = xmin

        enddo
      enddo

      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine previousday(yr,mo,dy,yrp,mop,dyp)
c this subroutine gives the year,month and day of the previous day, given the
c input year, month, day
c
      implicit none
      include 'comgeneral.h'

      integer     yr,mo,dy,yrp,mop,dyp

      if(dy.eq.1) then
        if(mo.eq.1) then
          yrp = yr-1
          mop = 12
          call lengthofmonth(yrp,mop,dyp)
        else
          yrp = yr
          mop = mo - 1
          call lengthofmonth(yr,mop,dyp)
        endif
      else
        yrp = yr
        mop = mo
        dyp = dy - 1
      endif

      return
      end


	subroutine error(string)
	character *(*) string
	write(0,*) trim(string)
	write(*,*) trim(string)
	call abort
	end

