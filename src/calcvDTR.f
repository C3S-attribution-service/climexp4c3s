c
c---------------------------------------------------------------------------
c
      subroutine calcvDTR(fyear,lyear,a,qc,b)
c this subroutine calculates the input to calculate index vDTR. Here we simply
c calculate for day 1,365: ( (tx_i - tn_i) - (tx_{i-1} - tn_{i-1} )
c input a is already tx_i - tn_i
c 
      implicit none
      include 'comgeneral.h'

      real*8        a(yrbeg:yrend,12,31)
      integer       qc(yrbeg:yrend,12,31)
      real*8        anew(yrbeg:yrend,12,31)
      real*8        b(yrbeg:yrend,nseason)
      real*8        absentr,absenr,xm
      integer       i,j,k,l,fyear,lyear,ip,jp,kp,length,absen,nm,fm,lm
      logical       pass

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

c initialize b
      do i=fyear,lyear
        do j=1,nseason
          b(i,j) = absentr
        enddo
      enddo

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

c----
      do i=fyear,lyear
        do j=1,nseason
          pass= .false.
          xm = 0.0d0
          nm = 0
          if(j.eq.1) then
            do k=1,12
              call lengthofmonth(i,k,length)
              do l=1,length
                if((anew(i,k,l).gt.absenr).and.
     +            ((k.ne.1).or.(l.ne.1))) then
                  xm = xm + anew(i,k,l)
                  nm = nm + 1
                endif
              enddo
            enddo
            if(nm.ge.pyear) pass = .true.

          elseif(j.eq.2) then
            if(i.gt.yrbeg)then
              do k=10,12
                call lengthofmonth(i,k,length)
                do l=1,length
                  if((anew(i,k,l).gt.absenr).and.
     +              ((k.ne.10).or.(l.ne.1))) then
                    xm = xm + anew(i,k,l)
                    nm = nm + 1
C                    write(6,*) i,k,l,anew(i,k,l)
                  endif
                enddo
              enddo
            endif
            do k=1,3
              call lengthofmonth(i,k,length)
              do l=1,length
                if(anew(i,k,l).gt.absenr) then
                  xm = xm + anew(i,k,l)
                  nm = nm + 1
C                  write(6,*) i,k,l,anew(i,k,l)
                endif
              enddo
            enddo
            if(nm.ge.p6month) pass = .true.

          elseif(j.eq.3) then
            do k=4,9
              call lengthofmonth(i,k,length)
              do l=1,length
                if((anew(i,k,l).gt.absenr).and.
     +            ((k.ne.4).or.(l.ne.1))) then
                  xm = xm + anew(i,k,l)
                  nm = nm + 1
                endif
              enddo
            enddo
            if(nm.ge.p6month) pass = .true.

          elseif(j.eq.4) then
            if(i.gt.yrbeg)then
              k=12
              call lengthofmonth(i,k,length)
              do l=1,length
                if((anew(i,k,l).gt.absenr).and.(l.ne.1)) then
                  xm = xm + anew(i,k,l)
                  nm = nm + 1
                endif
              enddo
            endif
            do k=1,2
              call lengthofmonth(i,k,length)
              do l=1,length
                if(anew(i,k,l).gt.absenr) then
                  xm = xm + anew(i,k,l)
                  nm = nm + 1
                endif
              enddo
            enddo
            if(nm.ge.pseason) pass = .true.

         elseif(j.le.7) then
            call bounds(j,fm,lm)
            do k=fm,lm
              call lengthofmonth(i,k,length)
              do l=1,length
                if((anew(i,k,l).gt.absenr).and.
     +            ((k.ne.fm).or.(l.ne.1))) then
                  xm = xm + anew(i,k,l)
                  nm = nm + 1
                endif
              enddo
            enddo
            if(nm.ge.pseason) pass = .true.
         else
            k = j - 7
              call lengthofmonth(i,k,length)
              do l=1,length
                if((anew(i,k,l).gt.absenr).and.(l.ne.1)) then
                  xm = xm + anew(i,k,l)
                  nm = nm + 1
                endif
              enddo
            if(nm.ge.pmonth) pass = .true.

          endif

          if(pass) then
            b(i,j) = xm/dble(nm)
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
