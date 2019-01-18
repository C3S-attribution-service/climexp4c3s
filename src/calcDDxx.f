      subroutine calcDDxx(fyear,lyear,lthresh,uthresh,a,qc,b)
c this routine calculates the number of days with lthresh < a <= uthresh
c
c for winds from the north, this gives a problem: transform winds with dd \in [0, 45]
c to [360, 405] before going through the procedure.
c
c the procedure is to calculate the number of days per month first,
c and then to calculate the values for each of the seasons
      implicit none
      include 'comgeneral.h'

      real*8        a(yrbeg:yrend,12,31)
      integer       qc(yrbeg:yrend,12,31)
      real*8        b(yrbeg:yrend,nseason)
      integer       c(yrbeg:yrend,12,2)
      integer       fyear,lyear,i,j,k,length,absen,nsum,npre
      real*8        absenr,absentr,lthresh,uthresh

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

      do i=fyear,lyear
        do j=1,12
          call lengthofmonth(i,j,length)
          nsum = 0
          npre = 0
          do k=1,length
            if(qc(i,j,k).eq.0)then
              npre = npre + 1
c this is a fix to be able to count the days with northerly winds
              if(a(i,j,k).le.45.0d0) a(i,j,k)=a(i,j,k)+360.0

              if((a(i,j,k).gt.lthresh).and.(a(i,j,k).le.uthresh)) 
     +          nsum = nsum + 1
            endif
          enddo
          c(i,j,1) = nsum
          c(i,j,2) = npre
c         write(6,*) i,nsum,npre
        enddo
      enddo

c accumulate the monthly indices over the appropriate months
      call calcSeason(fyear,lyear,c,b)

      return
      end
c
c---------------------------------------------------------------------------
c
