      subroutine calcRX5day(fyear,lyear,a,qc,b)
c this routine calculates the maximum amount of daily rain accumulated over 5 day periods
c
c the procedure is to calculate the max. amount of daily rain per month first,
c and then to calculate the values for each of the seasons
c
c array cp is an array with info if the number of "present" values pass the threshold
c the array cp is input to this routine
      implicit none
      include 'comgeneral.h'

      real*8        a(yrbeg:yrend,12,31),adum(yrbeg:yrend,12,31)
      integer       qc(yrbeg:yrend,12,31),qcdum(yrbeg:yrend,12,31)
      real*8        b(yrbeg:yrend,nseason)
      integer       fyear,lyear,i,j,k,l,ii,jj,kk,length,absen,nsum,npre
      integer       nm,fm,lm,N
      parameter     (N=5)
      real*8        absenr,absentr,xsum
      logical       cp(yrbeg:yrend,nseason)

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

c initialize array
      do i=fyear,lyear
        do j=1,12
          do k=1,31
            adum(i,j,k) = absentr
          enddo
        enddo

        do j=1,nseason
          b(i,j) = absentr
        enddo
      enddo

c accumulate the rain over N day period preceeding day (i,j,k) and put it in array adum
      do i=fyear,lyear
        do j=1,12
          call lengthofmonth(i,j,length)
          do k=1,length
            xsum = 0.0d0
            npre = 0
            !!!write(6,*) i,j,k
            do l=0,N-1
              call adjustdate(i,j,k,-l,ii,jj,kk)
              if ( ii.ge.yrbeg ) then
                if(qc(ii,jj,kk).eq.0) then
                  xsum = xsum + a(ii,jj,kk)
                  npre = npre + 1
                endif
              end if
            enddo
            if(npre.eq.N) then
              adum(i,j,k) = xsum
              qcdum(i,j,k)=0
            else
              qcdum(i,j,k)=1
            endif
          enddo
        enddo
      enddo

      call calcRXday(fyear,lyear,adum,qcdum,b)

      return
      end
