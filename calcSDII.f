      subroutine calcSDII(fyear,lyear,a,qc,b)
c this routine calculates a simple daily intensity index
c
c the procedure is to calculate the number of wet days per month first,
c and then to calculate the values for each of the seasons
      implicit none
      include 'comgeneral.h'

      real*8        a(yrbeg:yrend,12,31)
      integer       qc(yrbeg:yrend,12,31)
      real*8        b(yrbeg:yrend,nseason)
      real*8        cr(yrbeg:yrend,12,2)
      integer       fyear,lyear,i,j,k,length,absen,nsum,npre
      real*8        absenr,absentr,xsum

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

c initialize array
      do i=fyear,lyear
        do j=1,nseason
          b(i,j) = absentr
        enddo
      enddo

      do i=fyear,lyear
        do j=1,12
          call lengthofmonth(i,j,length)
          xsum = 0.0d0
          nsum = 0
          npre = 0
          do k=1,length
            if(qc(i,j,k).eq.0) then
              npre = npre + 1
              if(a(i,j,k).ge.1.0d0) then
                nsum = nsum + 1
                xsum = xsum + a(i,j,k)
              endif
            endif
          enddo
          if(nsum.gt.0) then
            cr(i,j,1) = xsum
          else
            cr(i,j,1) = 0.0d0
          endif
          cr(i,j,2) = dble(npre)
        enddo
      enddo

c average SDII
      call calcAverage(fyear,lyear,cr,b)

      return
      end
c
c---------------------------------------------------------------------------
c
