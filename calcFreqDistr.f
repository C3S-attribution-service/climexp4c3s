      subroutine calcFreqDistr(a,qc,p10,p90)
c this routine calculates the 10th and 90th percentiles for a nwl length window 
c centered on each calender day in the period calyrbeg, calyrend
      implicit none
      include 'comgeneral.h'

      integer       mxlength,mxperc
      parameter     (mxlength=(calyrend-calyrbeg+1)*nwl + 10,mxperc=5)
      real*8        a(yrbeg:yrend,12,31)
      integer       qc(yrbeg:yrend,12,31)
      real*8        p10(calyrbeg:calyrend+1,12,31)
      real*8        p90(calyrbeg:calyrend+1,12,31),dum(3*365)
      real*8        xdata(mxlength)
      integer       i,j,k,l,ii,jj,kk,m,n,length,ndata,npres
      integer       nwl2,absen,nperc
      integer       nthresh(mxperc)
      real*8        absenr,absentr,percentiles(mxperc)

      common/ndatablock/ndata
      common/xdatablock/xdata

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

      nwl2 = (nwl - 1)/2

      nperc = 2
      nthresh(1) = 10
      nthresh(2) = 90

      do j=1,12
        call lengthofmonth(i,j,length)
        do k=1,length

          if (zhang) then
c do the Zhang et al. approach: leave one year out
            do m=calyrbeg,calyrend

              npres = 0
              ndata = 0
              do i=calyrbeg,calyrend
                if(m.eq.i) goto 120

c count how many years are present
                if(qc(i,j,k).eq.0) npres=npres+1

c pool data from an nwl-length window
                do l=-nwl2,nwl2
                  call adjustdate(i,j,k,l,ii,jj,kk)
                  if(qc(ii,jj,kk).eq.0) then
                    ndata = ndata + 1
                    xdata(ndata) = a(ii,jj,kk)
                    !!!print *,'adding1 ',ndata,i,j,k,l,ii,jj,kk
                  endif
                enddo
 120            continue
              enddo
c add the extra year
              n = m+1
              if(n.gt.calyrend) n = calyrbeg
              if(qc(n,j,k).eq.0) npres=npres+1
c pool data from an nwl-length window
              do l=-nwl2,nwl2
                call adjustdate(n,j,k,l,ii,jj,kk)
                if(qc(ii,jj,kk).eq.0) then
                  ndata = ndata + 1
                  xdata(ndata) = a(ii,jj,kk)
                  !!!print *,'adding2 ',ndata,n,j,k,l,ii,jj,kk
                endif
              enddo

c if the amount of data does not equal or exceed the availability threshold: set percentiles to absent
              if(dble(npres)/dble(calyrend-calyrbeg+1).
     +                                 ge.percentage) then

                if(.not.fancy) then
c use the "poor men's" percentile estimate: sort data and take the 10th & 90th percentile
                  call poormenperc(nperc,nthresh,percentiles)
                else
c use the "fancy men's" percentile estimate: fit a three-parameter Gamma function
                  call fancymenperc(nperc,nthresh,percentiles)
                endif 

                p10(m,j,k) = percentiles(1)
                p90(m,j,k) = percentiles(2)
              else 
                p10(m,j,k) = absentr
                p90(m,j,k) = absentr
              endif

            enddo
          endif
c now calculate percentiles over the complete calibration period: no years left out
          m = calyrend+1
          ndata = 0
          do i=calyrbeg,calyrend
c pool data from an nwl-length window
            do l=-nwl2,nwl2
              call adjustdate(i,j,k,l,ii,jj,kk)
              if(qc(ii,jj,kk).eq.0) then
                ndata = ndata + 1
                xdata(ndata) = a(ii,jj,kk)
              endif
            enddo
          enddo

c if the amount of data does not equal or exceed the availability threshold: set percentiles to absent
          if(dble(npres)/dble(calyrend-calyrbeg+1).ge.percentage) then

            if(.not.fancy) then
c use the "poor men's" percentile estimate: sort data and take the 10th & 90th percentile
              call poormenperc(nperc,nthresh,percentiles)
            else
c use the "fancy men's" percentile estimate: fit a three-parameter Gamma function
              call fancymenperc(nperc,nthresh,percentiles)
            endif

            p10(m,j,k) = percentiles(1)
            p90(m,j,k) = percentiles(2)
          else
            p10(m,j,k) = absentr
            p90(m,j,k) = absentr
          endif

        enddo
      enddo

c there is a chance that not all calenderdays are filled with a value if the zhang approach
c is used. Fill these missing calender days with a polynomial approximation
c     if(zhang) then
c       call fillValue(p10)
c       call fillValue(p90)
c     endif

c if the Zhang et al. approach is NOT followed: then replace p10(calyrbeg:calyrend)
c with the value in p10(calyrend+1)
      if(.not.zhang) then
        do m=calyrbeg,calyrend
          do j=1,12
            call lengthofmonth(m,j,length)
            do k=1,length
              p10(m,j,k) = p10(calyrend+1,j,k)
              p90(m,j,k) = p90(calyrend+1,j,k)
            enddo
          enddo
        enddo
      endif

c one may want to smooth the percentiles
c     call smoothPercentiles(p10)
c     call smoothPercentiles(p90)

c Feb. 29th is always a problem: substitute the average between Feb 28th and March 1st
      do i=calyrbeg,calyrend+1
        if((p10(i,2,28).gt.absenr).and.(p10(i,3,1).gt.absenr))then
          p10(i,2,29) = (p10(i,2,28) + p10(i,3,1))/2.0d0
        else
          p10(i,2,29) = absentr
        endif

        if((p90(i,2,28).gt.absenr).and.(p90(i,3,1).gt.absenr))then
          p90(i,2,29) = (p90(i,2,28) + p90(i,3,1))/2.0d0
        else
          p90(i,2,29) = absentr
        endif
      enddo

 123  format(e14.4)

c test
c     do m=calyrbeg,calyrend+1
c       do j=1,12
c         call lengthofmonth(i,j,length)
c         do k=1,length
c           write(6,234) m,j,k, p90(m,j,k)
c         enddo
c       enddo
c     enddo

 234  format(I5,2I3,f10.2)

      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine poormenperc(nperc,nthresh,percentiles)
c this routine calculates the 10th and 90th percentile of the distribution based
c on the data in the array xdum
c
c here the simple approach is used
      implicit none
      include 'comgeneral.h'

      integer       mxlength
      parameter     (mxlength=(calyrend-calyrbeg+1)*nwl + 10)
      integer       N,nperc
      real*8        xdata(mxlength)
      integer       i,j,k,nthresh(nperc),nthreshd(nperc)
      real*8        percentiles(nperc)
c NAG things
      integer       ifail,irank(mxlength),absen
      real*8        absenr,absentr
      external      m01daf

      common/ndatablock/N
      common/xdatablock/xdata

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

      if(N.ge.pdayspresent) then
c sort array xdum into ascending order
        ifail = 0
        call m01daf(xdata,1,N,'Ascending',irank,ifail)

        do j=1,nperc
          nthreshd(j) = nint(dble(nthresh(j)*N)/100.0d0)

          do i=1,N
            if(irank(i).eq.nthreshd(j)) percentiles(j) = xdata(i)
          enddo
        enddo
      else
        do j=1,nperc
          percentiles(j) = absentr
        enddo
      endif

      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine fancymenperc(nperc,nthresh,percentiles)
c this routine calculates the 10th and 90th percentile of the distribution based
c on the data in the array xdum
c
c here the fancy approach is used
      implicit none
      include 'comgeneral.h'

      integer       mxlength
      parameter     (mxlength=(calyrend-calyrbeg+1)*nwl + 10)
      integer       N,nperc,absen
      real*8        absenr,absentr
      real*8        xdata(mxlength),array(mxlength)
      integer       i,j,k,nthresh(nperc)
      real*8        percentiles(nperc)
c NAG things
      integer       ifail,nclass
      external      g01aaf,g08cgf,g01aef,e04abf,g01fff
      real*8        g01fff,xmx,xmn,xthresh

      real*8        dint
      parameter     (nclass=40,dint=1.0)
      real*8        xmin,xmax,s2,s3,s4,wtsum,wt(mxlength),chisq
      real*8        eps,t,a,b
      real*8        cint(nclass)
      integer       ifreq(nclass)
      integer       maxcal
      integer       length,ndum,mdum,probdist(200)
      real*8        xd,x,alpha,beta,D,offset,xmean
      external      m01daf,func

      common/ndatablock/N
      common/xdatablock/xdata
      common/xfreq/cint
      common/shape/xmean,alpha,beta

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

      if(N.ge.pdayspresent) then

c fit a Gamma distribution to the data using maximum-likelihood
c theory is in Wilks p89
c
c if the data is negatively skewed, then it should be multiplied with -1 to make
c it positively skewed.
c the gamma distribution assumes a positively skewed distribution
        ifail = 0
        call g01aaf(N,xdata,0,wt,xmean,s2,s3,s4,xmin,xmax,wtsum,ifail)

c flip the distribution if neccessary
        if(s3.lt.0.0) then
          call flip(N,xdata)
          xd = -xmin
          xmin = -xmax
          xmax = xd
          xmean = -xmean
        endif

c the NAG routine g08cgf requires the distribution to be on the positive real axis
c shift the distribution with a distance -xmin + 1
        call shift(N,xmin-1.0,xdata,array)
        do i=1,N
          xdata(i) = array(i)
        enddo
        xmean = xmean - (xmin-1.0)

c construct the class-boundaries
        do i=1,nclass
          cint(i) = (i-1)*dint
        enddo

c construct frequency distribution of array
c       ifail = 0

c       call g01aef(N,nclass,xdata,1,cint,ifreq,xmn,xmx,ifail)

c       do i=1,nclass
c         write(6,*) i,-(cint(i)+(xmin-1.0)),ifreq(i)
c       enddo

        a = -1.0
        b = 0.99
        eps = 1.0D-6
        t = 1.0D-4
        x = 0.0
        maxcal = 50
        ifail = 0

        call e04abf(func,eps,t,a,b,maxcal,x,chisq)

c this calculates the percentiles
        do i=1,nperc
          xthresh = dble(nthresh(i))/100.0d0
          ifail = 0
          percentiles(i) = g01fff(xthresh,alpha,beta,0.0,ifail)
        enddo
  
c transform the percentiles back to the original distribution
c if s3: flip it!
        if(s3.lt.0.0) then
          if(nperc.gt.2) call error('can only do this trick with two')
          b = percentiles(1)
          percentiles(1) = percentiles(2)
          percentiles(2) = b
        endif

        do i=1,nperc
          percentiles(i) = percentiles(i) + (xmin - 1.0D0) + x

          if(s3.lt.0.0) then
            percentiles(i) = -percentiles(i)
          endif
        enddo
      else
        do i=1,nperc
          percentiles(i) = absentr
        enddo
      endif

c     write(6,59) vp10,vp90

 59   format('P10, P90: ',3f8.2)
 60   format('Skewness coefficient: ',e14.4)
 62   format('Chi squared: ',I5,e14.4)
 63   format('optimal shift: ',e14.4)
 65   format('no. iterations: ',I4)

      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine func(offset,chisq)
      implicit none
      include 'comgeneral.h'

      integer       mxlength
      parameter     (mxlength=(calyrend-calyrbeg+1)*nwl + 10)
      integer       nclass
      parameter     (nclass=40)
      integer       N
      real*8        xdata(mxlength),array(mxlength)
      external      g01aaf
      integer       i,ifail,ndf
      integer       ifreq(nclass)
      real*8        xd,D,alpha,beta,chisq,p
      real*8        cint(nclass),par(2),prob(nclass),eval(nclass)
      real*8        offset,chisqi(nclass)
      real*8        wt(mxlength),xmean,s2,s3,s4,xmin,xmax,wtsum
      external      g08cgf,g01aef

      common/ndatablock/N
      common/xdatablock/xdata
      common/xfreq/cint
      common/shape/xmean,alpha,beta

c shift the distribution
      call shift(N,offset,xdata,array)

c construct frequency distribution of shifted array
      ifail = 0

      call g01aef(N,nclass,array,1,cint,ifreq,xmin,xmax,ifail)

c     write(6,*) 'shift: ',offset
c     do i=1,nclass
c       write(6,*) i,cint(i),ifreq(i)
c     enddo

      xd = 0.0D0
      do i=1,N
        xd = xd + log(xdata(i) - offset)
      enddo
      D = log(xmean - offset) - xd/dble(N)
      alpha = (1.0 + sqrt(1.0 + 4*D/3.0))/(4.0*D)
      beta = (xmean - offset)/alpha

c     write(6,*) 'alpha & beta: ',alpha, beta

      par(1) = alpha
      par(2) = beta

c do the goodness-of-fit test
      ifail = 1
      call g08cgf(nclass,ifreq,cint,'G',par,0,prob,chisq,p,ndf,
     +            eval,chisqi,ifail)

c if ifail = 10, then we can go on
      if((ifail.ne.0).and.(ifail.ne.10)) then
        write(6,63) ifail
        stop
      endif

c     do i=1,nclass
c       write(6,61) i,cint(i),ifreq(i),eval(i),chisqi(i)
c     enddo
c     write(6,*)

c     write(6,62) chisq

 61   format(I4,e14.4,I4,2e14.4)
 62   format('Chi squared: ',e14.4)
 63   format('ifail is: ',I4)

      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine flip(n,x)
      implicit none

      integer       i,n
      real*8        x(n)

      do i=1,n
        x(i) = -1.0*x(i)
      enddo

      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine shift(n,offset,x,y)
      implicit none

      integer       i,n
      real*8        offset,x(n),y(n)

      do i=1,n
        y(i) = x(i) - offset
      enddo

      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine fillValue_polynomial(perc)
c this routine fills any absent value in perc with value derived from
c neighboring values.
c
c the procedure is:
c 1) check if there are absent values at all
c 2) if there is a year in perc with missing values:
c    3) fit a polynomial through the year with the missing data
c    4) evaluate the polynomial on the calender date with the missing data
c
      implicit none
      include 'comgeneral.h'

      integer       norder,M
      parameter     (norder=6,M=366)
      integer       absen,i,j,k,length,np,nq,ifail
      integer       npmin,npmax
      real*8        perc(calyrbeg:calyrend+1,12,31)
      real*8        absentr,absenr,xarg,fit
      logical       missing
      real*8        x(M),y(M),w(M),s(norder+1),a(norder+1,norder+1)
      real*8        ak(norder+1),work1(3*M),work2(2*(norder+1))
      external      e02adf,e02aef

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)


c check if there is a missing value
      do i=calyrbeg,calyrend+1
        missing  = .false.
        do j=1,12
          call lengthofmonth(i,j,length)
          do k=1,length
            if(perc(i,j,k).lt.absenr) missing = .true.
          enddo
        enddo
c
c if missing = .true., fit polynomial
        if(missing)then
          np = 0
          nq = 0
          do j=1,12
            call lengthofmonth(i,j,length)
            do k=1,length
              np = np + 1
              if(perc(i,j,k).gt.absenr) then
                nq = nq + 1
                x(nq) = dble(np)
                y(nq) = perc(i,j,k)
                w(nq) = 1.0d0
              endif
            enddo
          enddo

c fit polynomial
          ifail = 0
          write(6,*) 'nq: ',nq
          call e02adf(nq,norder+1,norder+1,x,y,w,work1,work2,a,s,ifail)

c fill perc with the evaluation of the polynomial fit
          npmax = np
          npmin = 1
          np = 0
          do j=1,12
            call lengthofmonth(i,j,length)
            do k=1,length
              np = np + 1
              if(perc(i,j,k).lt.absenr) then
                ifail = 0
                xarg = dble(np - npmin)/dble(npmax - npmin)
                call e02aef(norder+1,ak,xarg,fit,ifail)
                perc(i,j,k) = fit
              endif
            enddo
          enddo
        endif
      enddo

      return
      end
