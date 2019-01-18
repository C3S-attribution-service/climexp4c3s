      subroutine calcRXXptot(fyear,lyear,a,qc,rsum,b)
c this routine calculates the precipitation fraction due to (extremely) wet days
c
c input:
c   fyear : first year with non-missing data
c   lyear : last year with non-missing data
c   a     : daily data array
c   rsum  : sum of precip. in excess of some percentile
c output:
c   b     : array with precip fraction due to moderate wet days
c   
      implicit none
      include 'comgeneral.h'

      real*8        a(yrbeg:yrend,12,31)
      integer       qc(yrbeg:yrend,12,31)
      real*8        b(yrbeg:yrend,nseason),rsum(yrbeg:yrend,nseason)
      real*8        btot(yrbeg:yrend,nseason)
      real*8        c(yrbeg:yrend,12,2)
      integer       fyear,lyear,i,j,k,l,length,absen,npre,nzhang
      real*8        absenr,absentr,xsum

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

c initialize array
      do i=fyear,lyear
        do j=1,12
          c(i,j,1) = 0.0d0
          c(i,j,2) = 0.0d0
        enddo
        do l=1,nseason
          b(i,l) = absentr
          btot(i,l) = absentr
        enddo
      enddo

      do i=fyear,lyear
        do j=1,12
          call lengthofmonth(i,j,length)
          xsum = 0.0d0
          npre = 0
          do k=1,length
            if(qc(i,j,k).eq.0) then
              npre = npre + 1
              xsum = xsum + a(i,j,k)
            endif
          enddo
          c(i,j,1) = xsum
          c(i,j,2) = dble(npre)
        enddo
      enddo

c average the monthly indices over the appropriate months
      call calcSeasonR(fyear,lyear,c,btot)

c calculate the fraction
      do i=fyear,lyear
        do j=1,nseason
          if((rsum(i,j).gt.absenr).and.(btot(i,j).gt.absenr)) then
            if(btot(i,j).gt.0.0d0) then
              b(i,j) = 1.0d2*rsum(i,j)/btot(i,j)
            else
              b(i,j) = 0.0d0
            endif
          else
            b(i,j) = absentr
          endif
        enddo
      enddo

      return
      end
c
c---------------------------------------------------------------------------
c
