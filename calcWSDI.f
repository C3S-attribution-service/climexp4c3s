      subroutine calcWSDI(fyear,lyear,a,qc,p90,b)
c this routine calculates the number warm spell duration index
c
c the procedure is:
c 1) to check if T > T_90p  for each day
c    if so, an integer dummy array is filled with "1" for that day, 0 otherwise
c    the length of this dummy array is 365 days, plus 3 at either end of the array
c 2) data is cumulated if the entry is "1"
c ---- steps 1 and 2 are combined ---
c 3) a moving 6-day window is then applied and data in the dummy array is set
c    to "true" if a day is part of a 6-day consecutive period with value 1
c 3) monthly means of the index are then calculated
c 4) seasonal means are then calculated
c
      implicit none
      include 'comgeneral.h'

      real*8        a(yrbeg:yrend,12,31)
      integer       qc(yrbeg:yrend,12,31)
      real*8        b(yrbeg:yrend,nseason)
      real*8        p90(calyrbeg:calyrend+1,12,31)
      integer       ndum(-5:375)
      logical       ldum(-5:375)
      integer       c(yrbeg:yrend,12,2)
      integer       fyear,lyear,i,j,k,length
      integer       absen,nsum,npre,nday,nprev,maxday
      integer       nzhang
      real*8        absenr,absentr
      logical       leap

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

c initialize array
      do i=fyear,lyear
        do j=1,12
          do k=1,2
            c(i,j,k) = 0
          enddo
        enddo

        do j=1,nseason
          b(i,j) = absentr
        enddo
      enddo

c check if a(i,j,k) passes the threshold and make the array cumulative
      nprev = 0
      do i=fyear,lyear

        do nday=-5,375
          ndum(nday) = 0
          ldum(nday) = .false.
        enddo

c add the last 5 days of Dec of year i-1 to the array ndum
        nday = -6
        do k=27,31
          nday = nday + 1
          if((qc(i-1,12,k).eq.0).and.
     +       (p90(nzhang(i-1),12,k).gt.absen)) then
            if(a(i-1,12,k).gt.p90(nzhang(i-1),12,k))
     +        ndum(nday) = nprev + 1
          endif
          nprev = ndum(nday)
        enddo

        do j=1,12
          call lengthofmonth(i,j,length)
          npre = 0
          do k=1,length
            nday = nday + 1
            if((qc(i,j,k).eq.0).and.(p90(nzhang(i),j,k).gt.absen)) then
              npre = npre + 1
              if(a(i,j,k).gt.p90(nzhang(i),j,k)) 
     +          ndum(nday) = nprev + 1
            endif
            nprev = ndum(nday)
          enddo
          c(i,j,2) = npre
        enddo

c add the first 5 days of Jan of year i+1 to the array ndum
        do k=1,6
          nday = nday + 1
          if((qc(i+1,1,k).eq.0).and.
     +       (p90(nzhang(i+1),1,k).gt.absen)) then
            if(a(i+1,1,k).gt.p90(nzhang(i+1),1,k)) 
     +          ndum(nday) = nprev + 1
          endif
          nprev = ndum(nday)
        enddo

c use a 6-day moving window to assess which days are part of 6-day consecutive periods
        call leapyr(i,leap)
        if(leap) then
          maxday = 366
        else
          maxday = 365
        endif

c initialize array
        do nday=-5,maxday+5
          ldum(nday) = .false.
        enddo

c go through the array: if value .eq. 6 -> set value plus 5 preceeding values to "on"
c                       if value .ge. 7 -> set value to "on"

        do nday=1,maxday
          if(ndum(nday).eq.6) then
            do k=1,6
              ldum(nday+1-k) = .true.
            enddo
          endif
          if(ndum(nday).ge.7) ldum(nday) = .true.
        enddo

c calculate monthly means of WSDI
        nday = 0
        do j=1,12
          call lengthofmonth(i,j,length)
          nsum = 0
          do k=1,length
            nday = nday + 1
            if(ldum(nday)) nsum = nsum + 1
          enddo
          c(i,j,1) = nsum
        enddo

      enddo

c average the monthly indices over the appropriate months
      call calcSeason(fyear,lyear,c,b)

      return
      end
c
c---------------------------------------------------------------------------
c
