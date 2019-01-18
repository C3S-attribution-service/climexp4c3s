      subroutine calcCWD(fyear,lyear,a,qc,b)
c this routine determines the maximum number of consecutive wet days (index CWD)
c
c the procedure is to go through the data array, detect sequences of wet days and
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
      integer       season,spell,npre,nm
      real*8        absenr,absentr,ap
      logical       gt

c test for dry days: RR => 1 mm
      gt = .true.

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)
      fy = fyear

c initialize array
      do i=fyear,lyear
        do j=1,nseason
          b(i,j) = absentr
        enddo
      enddo

      do season=1,nseason
        if((season.ne.2).and.(season.ne.4))then
          call bounds(season,fm,lm)
c set previous value
          ap = absentr
          qcp = 9
          do i=fyear,lyear
            spell = 0
            do j=fm,lm
              call lengthofmonth(i,j,length)
              do k=1,length
                call addtospell(a(i,j,k),qc(i,j,k),
     +                          ap,      qcp,gt,1.0d0,spell)

                if(spell.gt.nint(b(i,season))) then
                  b(i,season)=dble(spell)
                end if
                ap = a(i,j,k)
                qcp = qc(i,j,k)
              enddo
            enddo
          enddo
        elseif(season.eq.2)then
c make sure that fyear > yrbeg
          !!!fyear = max(fyear,yrbeg+1) TOTALLY WRONG, ALREADY CHECKED, MESSES UP THE MINIMUM FRACTION CODE
          ap = 1.0d0
          qcp = 0
          do i=fyear+1,yrend
            spell = 0
            do j=10,12
              call lengthofmonth(i-1,j,length)
              do k=1,length
  
                call addtospell(a(i-1,j,k),qc(i-1,j,k),
     +                          ap,        qcp,gt,1.0d0,spell)
  
                if(spell.gt.nint(b(i,season))) b(i,season)=dble(spell)

                ap = a(i,j,k)
                qcp = qc(i,j,k)
              enddo
            enddo
            do j=1,3
              call lengthofmonth(i,j,length)
              do k=1,length

                call addtospell(a(i,j,k),qc(i,j,k),
     +                          ap,      qcp,gt,1.0d0,spell)

                if(spell.gt.nint(b(i,season))) b(i,season)=dble(spell)

                ap = a(i,j,k)
                qcp = qc(i,j,k)
              enddo
            enddo
          enddo
        else
c make sure that fyear > yrbeg
          !!!fyear = max(fyear,yrbeg+1) TOTALLY WRONG, ALREADY CHECKED, MESSES UP THE MINIMUM FRACTION CODE
          ap = 1.0d0
          qcp = 0
          do i=fyear+1,yrend
            spell = 0
            j=12
            call lengthofmonth(i-1,j,length)
            do k=1,length

              call addtospell(a(i-1,j,k),qc(i-1,j,k),
     +                        ap,        qcp,gt,1.0d0,spell)

              if(spell.gt.nint(b(i,season))) b(i,season)=dble(spell)

              ap = a(i,j,k)
              qcp = qc(i,j,k)
            enddo
            do j=1,2
              call lengthofmonth(i,j,length)
              do k=1,length

                call addtospell(a(i,j,k),qc(i,j,k),
     +                          ap,      qcp,gt,1.0d0,spell)

                if(spell.gt.nint(b(i,season))) b(i,season)=dble(spell)

                ap = a(i,j,k)
                qcp = qc(i,j,k)
              enddo
            enddo
          enddo
        endif
      enddo

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
