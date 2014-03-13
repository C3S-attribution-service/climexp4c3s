      subroutine calcRXXp(fyear,lyear,a,qc,p,rsum,b)
c this routine calculates the number of days with R > XX percentile
c
c the procedure is to calculate the number of cold days per month first,
c and then to calculate the values for each of the seasons
      implicit none
      include 'comgeneral.h'

      real*8        a(yrbeg:yrend,12,31),p(calyrbeg:calyrend+1,nseason)
      integer       qc(yrbeg:yrend,12,31)
      real*8        b(yrbeg:yrend,nseason),rsum(yrbeg:yrend,nseason)
      integer       c(yrbeg:yrend,nseason,2),nmonths
      integer       fyear,lyear,i,j,k,l,length,absen
      integer       nsum,npre,nzhang,ii,jj,nm
      real*8        absenr,absentr,thresh,xsum
      logical       pass

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

c theshold for wet days
      thresh = 1.0d0

c initialize array
      do i=fyear,lyear
        do l=1,nseason
          c(i,l,1) = 0
          c(i,l,2) = 0
          rsum(i,l) = absentr
          b(i,l) = absentr
        enddo
      enddo

      do i=fyear,yrend
        do l=1,nseason
          xsum = 0.0d0
          nsum = 0
          npre = 0
          do j=1,nmonths(l)
            call selectmonth(l,i,j,ii,jj)
            call lengthofmonth(ii,jj,length)
            do k=1,length
              if((qc(ii,jj,k).eq.0).and.
     +           (p(nzhang(ii),l).gt.absen)) then
                npre = npre + 1
                if((a(ii,jj,k).ge.thresh).and.
     +             (a(ii,jj,k).gt.p(nzhang(ii),l))) then
                  nsum = nsum + 1
                  xsum = xsum + a(ii,jj,k)
                endif
              endif
            enddo
          enddo
          c(i,l,1) = nsum
          c(i,l,2) = npre
          rsum(i,l) = xsum
        enddo
      enddo

c check if the number of entries per season meet the threshold
      do i=fyear,lyear
        do j=1,nseason
          nm = c(i,j,2)
          pass = .false.
          if(j.eq.1) then
            if(nm.ge.pyear) pass = .true.
          elseif(j.le.3) then
            if(nm.ge.p6month) pass = .true.
          elseif(j.le.7)then
            if(nm.ge.pseason) pass = .true.
          else
            if(nm.ge.pmonth) pass = .true.
          endif

          if(pass) then
            b(i,j) = dble(c(i,j,1))
          else
            b(i,j) = absentr
            rsum(i,j) = absentr
          endif

        enddo
      enddo

      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine selectmonth(season,i,j,year,month)
c this function is used in the percentile calculations for RR
c given the season and the number of the month of that season, it returns the 
c adjusted year and the calender month
c
      implicit none
      include 'comgeneral.h'

      integer       season,i,j,year,month,ndum

      if(season.le.4) then
        if(season.eq.1) then
          month = j
          year = i
        elseif(season.eq.2)then
          ndum = mod(j+9,12)
          if(ndum.eq.0) ndum = 12
          month = ndum
          if(j.le.3) then
            year = i-1
          else
            year = i
          endif
        elseif(season.eq.3)then
          month = j+3
          year = i
        else
          ndum = mod(j+11,12)
          if(ndum.eq.0) ndum = 12
          month = ndum
          if(j.eq.1) then
            year = i-1
          else
            year = i
          endif
        endif
      else
        if(season.eq.5)then
          month = j+2
          year = i
        elseif(season.eq.6)then
          month = j+5
          year = i
        elseif(season.eq.7)then
          month = j+8
          year = i
        else
          month = season - 7
          year = i
        endif
      endif

      return
      end
c
c---------------------------------------------------------------------------
c
