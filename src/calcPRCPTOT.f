      subroutine calcPRCPTOT(fyear,lyear,a,qc,b)
c this routine calculates precipitation amounts due to wet days
c
c the procedure is to calculate accumulated precipitation per month first,
c and then to calculate the values for each of the seasons
      implicit none
      include 'comgeneral.h'

      real*8        a(yrbeg:yrend,12,31)
      integer       qc(yrbeg:yrend,12,31)
      real*8        b(yrbeg:yrend,nseason)
      real*8        c(yrbeg:yrend,12,2),c1(yrbeg:yrend,12)
      integer       fyear,lyear,i,j,k,length,absen,npre,nm,fm,lm
      real*8        absenr,absentr,xsum
      logical       cp(yrbeg:yrend,nseason)

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

c initialize array
      do i=fyear,lyear
        do j=1,12
          do k=1,2
            c(i,j,k) = 0.0d0
          enddo
        enddo

        do j=1,nseason
          b(i,j) = absentr
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
                if(a(i,j,k).ge.1.0d0) then
                   xsum = xsum + a(i,j,k)
                endif
             endif
          enddo
          c(i,j,1) = xsum
          c(i,j,2) = dble(npre)
        enddo
      enddo

c average the monthly indices over the appropriate months
      call calcSeasonR(fyear,lyear,c,b)

      return
      end
c
c---------------------------------------------------------------------------
c
