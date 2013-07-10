      subroutine calcdPabs24(fyear,lyear,a,qc,b)
c this routine calculates the monthly means of the storminess index dp(abs)24
c
c dp(abs)24 is the absolute 24-hourly atmospheric surface pressure variation
c
c the procedure is to calculate the monthly means first,
c and then to calculate the values for each of the seasons
      implicit none
      include 'comgeneral.h'

      real*8        a(yrbeg:yrend,12,31),anew(yrbeg:yrend,12,31)
      integer       qc(yrbeg:yrend,12,31)
      real*8        b(yrbeg:yrend,nseason)
      real*8        c(yrbeg:yrend,12,2)
      integer       i,j,k,length,absen,fyear,lyear,ip,jp,kp
      real*8        xsum,npre
      real*8        absenr,absentr

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

c calculate pressure difference day_i - day_(i-1)
c tis is the basis of the index
      do i=fyear,lyear
        do j=1,12
          call lengthofmonth(i,j,length)
          do k=1,length
            call previousday(i,j,k,ip,jp,kp)
            if((qc(i,j,k).eq.0).and.(qc(ip,jp,kp).eq.0)) then
              anew(i,j,k) = abs( a(i,j,k) - a(ip,jp,kp) )
            else
              anew(i,j,k) = absentr
            endif
          enddo
        enddo
      enddo

      do i=fyear,lyear
        do j=1,12
          call lengthofmonth(i,j,length)
          xsum = 0.0d0
          npre = 0.0d0
          do k=1,length
            if(anew(i,j,k).gt.absenr) then
              npre = npre + 1.0d0
              xsum = xsum + anew(i,j,k)
            endif
          enddo
          c(i,j,1) = xsum
          c(i,j,2) = npre
        enddo
      enddo

c calculate average 
      call calcAverage(fyear,lyear,c,b)

      return
      end
c
c---------------------------------------------------------------------------
c
