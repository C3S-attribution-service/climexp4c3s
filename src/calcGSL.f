      subroutine calcGSL(fyear,lyear,a,qc,b)
c this routine determines the growing season length (index GSL)
c GSL is defined for annual values
c
c the procedure is to go through the data array, detect the first 6-day sequences with TG > 5,
c detect the first 6-day sequences with TG < 5 and
c calculate the number of days between the sequences
c
      implicit none
      include 'comgeneral.h'

      real*8        a(yrbeg:yrend,12,31)
      integer       qc(yrbeg:yrend,12,31)
      real*8        b(yrbeg:yrend,nseason)
      integer       fyear,i,j,k,length,absen,lyear
      integer       fy,fm,fd,ly,lm,ld,nlength
      integer       season,spellf,spelll,npres
      real*8        absenr,absentr,ptemp
      logical       start,halt

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

c initialize array
      do i=fyear,lyear
        do j=1,nseason
          b(i,j) = absentr
        enddo
      enddo

      ptemp = 0.0d0
      do i=fyear,lyear
        start = .false.
        halt = .false.
        spellf = 0
        spelll = 0
        npres = 0
        do j=1,12
          call lengthofmonth(i,j,length)
          do k=1,length

c check for no. of present values
          if(qc(i,j,k).eq.0) npres = npres + 1

          if(.not.start) then
            if(qc(i,j,k).eq.0) then
              if(a(i,j,k).gt.5.0d0) then
                if((ptemp.gt.5.0d0).and.(ptemp.gt.absenr)) then
c a(i,j,k) contributes to the sequence of >5C days
                  spellf = spellf + 1
                else
c a(i,j,k) starts a new cold spell
                  spellf = 1
                endif
              else
c a(i,j,k) ends the sequece of warm days
                spellf = 0
              endif
            else
c a(i,j,k) ends the sequece of warm days since it is absent
              spellf = 0
            endif
          endif

 123      format(3I5,2f8.2,I4)

          if((spellf.eq.6).and.(.not.start)) then
            start = .true.
c set nlength to 5; in the next step it is updated to length 6 (as it should)
            nlength = 5
          endif

c sum the number of days 
          if((start).and.(.not.halt).and.(qc(i,j,k).eq.0))
     +          nlength = nlength + 1

c look for the first sequence of 6 days after July 1st with temp < 5C
          if((start).and.(.not.halt).and.(j.ge.7)) then
            if(qc(i,j,k).eq.0) then
              if(a(i,j,k).lt.5.0d0) then
                if((ptemp.lt.5.0d0).and.(ptemp.gt.absenr)) then
c a(i,j,k) contributes to the sequence of <5C days
                  spelll = spelll + 1
                else
c a(i,j,k) starts a new cold spell
                  spelll = 1
                endif
              else
c a(i,j,k) ends the sequece of cold days
                spelll = 0
              endif
            else
c a(i,j,k) ends the sequece of cold days since it is absent
              spelll = 0
            endif
          endif

          if(spelll.eq.6) halt = .true.

          ptemp = a(i,j,k)

          enddo
        enddo

c nlength now also includes the sequence of 6 days with days < 5C: subtract 6 days
        if(npres.gt.pyear) b(i,1) = nlength - 6

      enddo

      return
      end
c
c---------------------------------------------------------------------------
c
