      subroutine calcCFD(fyear,lyear,a,qc,b)
c this routine determines the maximum number of consecutive frost days (index CFD)
c
c the procedure is to go through the data array, detect sequences of frost days and
c characterize these sequences by their first date and last date (and/or length)
c
c Then determine for each sequence to which season it contributes
      implicit none
      include 'comgeneral.h'

      real*8        a(yrbeg:yrend,12,31)
      integer       qc(yrbeg:yrend,12,31)
      real*8        b(yrbeg:yrend,nseason)
      integer       c(yrbeg:yrend,12)
      integer       fyear,lyear,i,j,k,length,absen
      integer       fy,fm,fd,ly,lm,ld,nlength,qcp
      integer       season,spell,nmax,npre,nm
      real*8        absenr,absentr,ap
      logical       gt

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

c we test for frost days: T < 0
      gt = .false.

c initialize array
      do i=fyear,lyear
        do j=1,nseason
          b(i,j) = absentr
        enddo
      enddo

      do season=1,nseason
        if((season.ne.2).and.(season.ne.4))then
          call bounds(season,fm,lm)
          ap = 1.0d0
          qcp = 0
          do i=fyear,lyear
            nmax = -1
            spell = 0
            do j=fm,lm
              call lengthofmonth(i,j,length)
              do k=1,length

                call addtospell(a(i,j,k),qc(i,j,k),
     +                          ap,      qcp,gt,0.0d0,spell)

                if(spell.gt.nmax) nmax = spell

                ap = a(i,j,k)
                qcp = qc(i,j,k)
              enddo
            enddo
            b(i,season) = dble(max(nmax,0))
          enddo
        elseif(season.eq.2)then
          ap = 1.0d0
          qcp = 0
          do i=fyear+1,lyear
            spell = 0
            nmax = -1
            do j=10,12
              call lengthofmonth(i-1,j,length)
              do k=1,length
  
                call addtospell(a(i-1,j,k),qc(i-1,j,k),
     +                          ap,        qcp,gt,0.0d0,spell)
  
                if(spell.gt.nmax) nmax = spell

                ap = a(i,j,k)
                qcp = qc(i,j,k)
              enddo
            enddo
            do j=1,3
              call lengthofmonth(i,j,length)
              do k=1,length

                call addtospell(a(i,j,k),qc(i,j,k),
     +                          ap,      qcp,gt,0.0d0,spell)

                if(spell.gt.nmax) nmax = spell

                ap = a(i,j,k)
                qcp = qc(i,j,k)
              enddo
            enddo
            b(i,season) = dble(max(nmax,0))
          enddo
        else
          ap = 1.0d0
          qcp = 0
          do i=fyear+1,yrend
            spell = 0
            nmax = -1
            j=12
            call lengthofmonth(i-1,j,length)
            do k=1,length

              call addtospell(a(i-1,j,k),qc(i-1,j,k),
     +                        ap,        qcp,gt,0.0d0,spell)

              if(spell.gt.nmax) nmax = spell

              ap = a(i,j,k)
              qcp = qc(i,j,k)
            enddo
            do j=1,2
              call lengthofmonth(i,j,length)
              do k=1,length

                call addtospell(a(i,j,k),qc(i,j,k),
     +                          ap,      qcp,gt,0.0d0,spell)

                if(spell.gt.nmax) nmax = spell

                ap = a(i,j,k)
                qcp = qc(i,j,k)
              enddo
            enddo
            b(i,season) = dble(max(nmax,0))
          enddo
        endif
      enddo

c fill array c: it contains the information whether the number of "present" values in each season
c below of above the threshold
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

