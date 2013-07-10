      subroutine calcRXday(fyear,lyear,a,qc,b)
c this routine calculates the maximum amount of daily rain 
c
c the procedure is to calculate the max. amount of daily rain per month first,
c and then to calculate the values for each of the seasons
c
c array cp is an array with info if the number of "present" values pass the threshold
      implicit none
      include 'comgeneral.h'

      real*8        a(yrbeg:yrend,12,31)
      integer       qc(yrbeg:yrend,12,31)
      real*8        b(yrbeg:yrend,nseason)
      real*8        cx(yrbeg:yrend,12)
      integer       c(yrbeg:yrend,12)
      integer       fyear,lyear,i,j,k,length,absen,nsum,npre
      integer       nm,fm,lm
      real*8        absenr,absentr,xmax

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

c initialize array
      do i=fyear,lyear
        do j=1,12
          c(i,j) = 0
          cx(i,j) = absentr
        enddo

        do j=1,nseason
          b(i,j) = absentr
        enddo
      enddo

      do i=fyear,lyear
        do j=1,12
          call lengthofmonth(i,j,length)
          xmax = -1.0
          npre = 0
          do k=1,length
            if(qc(i,j,k).eq.0) then
              npre = npre + 1
              if(a(i,j,k).gt.xmax) xmax = a(i,j,k)
            endif
          enddo
          cx(i,j) = xmax
          c(i,j) = npre
        enddo
      enddo

c calculate maximum values for the seasons based on monthly values
      call calcMaxSeason(fyear,lyear,cx,b)

c fill array cp: it contains the information whether the number of "present" values in each season
c below of above the threshold
c initialize array
      do i=fyear,lyear
        do j=1,12
          call lengthofmonth(i,j,length)
          npre = 0.0d0
          do k=1,length
            if(qc(i,j,k).eq.0) npre = npre + 1
          enddo
          c(i,j) = npre
        enddo
      enddo

      do i=fyear,lyear
        do j=1,nseason
          nm = 0
          if(j.eq.1) then
            do k=1,12
              nm = nm + c(i,k)
            enddo
            if(nm.lt.pyear) b(i,j) = absentr

          elseif(j.eq.2) then
            if ( i.gt.yrbeg ) then
              do k=10,12
                nm = nm + c(i-1,k)
              enddo
            end if
            do k=1,3
              nm = nm + c(i,k)
            enddo
            if(nm.lt.p6month) b(i,j) = absentr

          elseif(j.eq.3) then
            do k=4,9
              nm = nm + c(i,k)
            enddo
            if(nm.lt.p6month) b(i,j) = absentr

          elseif(j.eq.4) then
            if ( i.gt.yrbeg ) then
              nm = c(i-1,12)
            end if
            do k=1,2
              nm = nm + c(i,k)
            enddo
            if(nm.lt.pseason) b(i,j) = absentr

          elseif(j.le.7) then
            call bounds(j,fm,lm)
            do k=fm,lm
              nm = nm + c(i,k)
            enddo
            if(nm.lt.pseason) b(i,j) = absentr
          else
            call bounds(j,fm,lm)
            do k=fm,lm
              nm = nm + c(i,k)
            enddo
            if(nm.lt.pmonth) b(i,j) = absentr
          endif

        enddo
      enddo

      return
      end
c
c---------------------------------------------------------------------------
c      
