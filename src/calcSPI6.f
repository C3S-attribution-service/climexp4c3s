      subroutine calcSPI6(fyear,lyear,a,qc,b)
c this routine calculates the standardized precipitation index
c
c     CAVIAT
c        SPI values in central part of the precip distribution can
c        be considered to be accurate and reasonable if the number
c        of independent non-zero accumulated precipitation amounts
c        is greater than about 40; values in the tails of the
c        distribution can be considered accurate and reasonable
c        if the number of non-zero accumulated amounts is greater
c        than about 60.
c
c     written by Ned Guttman in January 1998; some of the software
c     and logic based on software provided by John Kleist, Colorado
c     State Univ.
c
      implicit none
      include 'comgeneral.h'

      real*8        a(yrbeg:yrend,12,31)
      integer       qc(yrbeg:yrend,12,31)
      real*8        b(yrbeg:yrend,nseason)
      integer       fyear,lyear,i,j,k,length,absen,npre,nm,fm,lm,nrec
      real*8        absenr,absentr,xsum

      integer       nlen,maxyrs
      parameter     (nlen=7,maxyrs=yrend-yrbeg+1)

      real*8        pp(maxyrs*12),probne(maxyrs*12),
     1              pcpacc(maxyrs*12),pzero(12),spi(maxyrs*12),
     2              spiout(12),prob(12),pacc(12)
      real*8        spical(yrbeg:yrend,12,2)
      integer       num(12),numpos(12),nrun(nlen)
      real*8        par1(12),par2(12),par3(12),xmom1(12),
     1              xmom2(12),xmom3(12),x(maxyrs),y(maxyrs*12),
     2              index(maxyrs*12),temparr(maxyrs+1)
       data nrun /1,2,3,6,9,12,24/

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

c initialize array
      do i=fyear,lyear
        do j=1,nseason
          b(i,j) = absentr
        enddo
      enddo

c the input to Ned Guttman's SPI calculation are monthly precip. sums
c in a one dimensional array, starting with January
      nrec = 0
      do i=fyear,lyear
        do j=1,12
          nrec = nrec + 1
          call lengthofmonth(i,j,length)
          xsum = 0.0d0
          npre = 0
          do k=1,length
            if(qc(i,j,k).eq.0) then
              npre = npre + 1
              xsum = xsum + a(i,j,k)
            endif
          enddo
c require at least "percentage" of the data to be available
          if(dble(npre)/dble(length).ge.percentage) then
            pp(nrec) = xsum
          else
            pp(nrec) = absentr
          endif
c         write(6,1004) i,j,pp(nrec)
        enddo
      enddo


c----------------------------------
c
c Subroutine spipe3 (and everything it calls) is part of Ned's code
c
c----------------------------------

c calculating only i=3 relates to a 3-month window for calculating SPI
      i=4
c
c     compute SPI values, probabilities, accumulated precip
c
      call spipe3 (nrun(i),pp,par1,par2,par3,pzero,spi,probne,
     1     pcpacc,xmom1,xmom2,xmom3,absentr,num,numpos,lyear-fyear+1,
     2     index,y,x,temparr)
c
c
c      put array spi in the usual form
c
      do j=fyear,lyear
        do k=1,12
c         write(6,*) j,k,spi((j-fyear)*12+k)
          spical(j,k,1) = spi((j-fyear)*12+k)
c the second place contains the # entries of the month: put it at 31
c the SPI algorithm can't handle absent data and has filled this in already
          spical(j,k,2) = 31.0d0
        enddo
      enddo

c average
      call calcAverageSPI(fyear,lyear,spical,b)

      return
      end
c
c---------------------------------------------------------------------------
c
