      subroutine calcR10mm(fyear,lyear,a,qc,b)
c this routine calculates the number of days with heavy precip (RR => thresh mm)
c
c the procedure is to calculate the number of heavy precip days per month first,
c and then to calculate the values for each of the seasons
      implicit none
      include 'comgeneral.h'

      real*8        a(yrbeg:yrend,12,31)
      integer       qc(yrbeg:yrend,12,31)
      real*8        b(yrbeg:yrend,nseason)
      integer       c(yrbeg:yrend,12,2)
      integer       fyear,lyear,i,j,k,length,absen,nsum,npre
      real*8        absenr,absentr,thresh
      parameter     (thresh=10.0d0)

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

      do i=fyear,yrend
        do j=1,12
          call lengthofmonth(i,j,length)
          nsum = 0
          npre = 0
          do k=1,length
            if(qc(i,j,k).eq.0) then
              npre = npre + 1
              if(a(i,j,k).ge.thresh) nsum = nsum + 1
            endif
          enddo
          c(i,j,1) = nsum
          c(i,j,2) = npre
        enddo
      enddo

c average the monthly indices over the appropriate months
      call calcSeason(fyear,lyear,c,b)

      return
      end
c
c---------------------------------------------------------------------------
c
