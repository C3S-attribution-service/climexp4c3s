      subroutine calcMIN(fyear,lyear,a,qc,b)
c this routine calculates the minimum value of the input array a
c
c the procedure is to calculate the min. value of a per month first,
c and then to calculate the values for each of the seasons
c
      implicit none
      include 'comgeneral.h'

      real*8        a(yrbeg:yrend,12,31)
      integer       qc(yrbeg:yrend,12,31)
      real*8        b(yrbeg:yrend,nseason)
      real*8        c(yrbeg:yrend,12)
      integer       cm(yrbeg:yrend,12)
      integer       fyear,lyear,i,j,k,length,absen,nsum,npre
      integer       nm,fm,lm
      real*8        absenr,absentr,xmax,xmin
      logical       cp(yrbeg:yrend,nseason)

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

c initialize array
      do i=fyear,lyear
        do j=1,12
          cm(i,j) = 0
          c(i,j) = absentr
        enddo

        do j=1,nseason
          b(i,j) = absentr
          cp(i,j) = .false.
        enddo
      enddo

      do i=fyear,lyear
        do j=1,12
          call lengthofmonth(i,j,length)
          xmin = 1.0D6
          npre = 0
          do k=1,length
            if(qc(i,j,k).eq.0) then
              npre = npre + 1
              if(a(i,j,k).lt.xmin) xmin = a(i,j,k)
            endif
          enddo
          c(i,j) = xmin
          cm(i,j) = npre
        enddo
      enddo

c fill array cp: it contains the information whether the number of "present" values in each season
c below of above the threshold

      do i=fyear,lyear
        do j=1,nseason
          nm = 0
          if(j.eq.1) then
            do k=1,12
              nm = nm + cm(i,k)
            enddo
            if(nm.ge.pyear) cp(i,j) = .true.
  
          elseif(j.eq.2) then
            if(fyear.gt.yrbeg)then
              do k=10,12
                nm = nm + cm(i-1,k)
              enddo
            endif
            do k=1,3
              nm = nm + cm(i,k)
            enddo
            if(nm.ge.p6month) cp(i,j) = .true.
  
          elseif(j.eq.3) then
            do k=4,9
              nm = nm + cm(i,k)
            enddo
            if(nm.ge.p6month) cp(i,j) = .true.
  
          elseif(j.eq.4) then
            if(fyear.gt.yrbeg) nm = cm(i-1,12)
            do k=1,2
              nm = nm + cm(i,k)
            enddo
            if(nm.ge.pseason) cp(i,j) = .true.
  
          elseif(j.le.7) then
            call bounds(j,fm,lm)
            do k=fm,lm
              nm = nm + cm(i,k)
            enddo
            if(nm.ge.pseason) cp(i,j) = .true.

          else
            nm = cm(i,j-7)
            if(nm.ge.pmonth) cp(i,j) = .true.
          endif

        enddo
      enddo

c calculate maximum values for the seasons based on monthly values
      call calcMinSeason(fyear,lyear,c,b)

c mask array b with the info in cp
      do i=fyear,lyear
        do j=1,nseason
          if(.not.cp(i,j)) then
            b(i,j) = absentr
          endif
        enddo
      enddo

      return
      end
c
c---------------------------------------------------------------------------
c      
