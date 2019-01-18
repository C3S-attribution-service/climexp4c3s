      subroutine calcFreqDistrRR(a,qc,p75,p95,p99)
c this routine calculates the 75th, 95th and 99th percentiles, averaged over each season,
c in the period calyrbeg, calyrend
c
c only data from wet days (RR > 1.0mm) are entered
c
c here a simple two-parameter Gamma fit is used (rather than the three-parameter or 
c Pearson type III fit)
      implicit none
      include 'comgeneral.h'

      integer       mxlength,mxperc
      parameter     (mxlength=(calyrend-calyrbeg+1)*365,mxperc=5)
      real*8        a(yrbeg:yrend,12,31)
      integer       qc(yrbeg:yrend,12,31)
      real*8        p75(calyrbeg:calyrend+1,nseason)
      real*8        p95(calyrbeg:calyrend+1,nseason)
      real*8        p99(calyrbeg:calyrend+1,nseason)
      real*8        dum(3*365)
      real*8        xdata(mxlength)
      integer       i,j,k,l,ii,jj,kk,length,ndata,nperc,m,nzhang,n
      integer       nmonths,month
      integer       nwl2,absen,nthresh(mxperc)
      real*8        absenr,absentr,percentiles(mxperc),thresh

      common/datablock/xdata
      common/intblock/ndata

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

      thresh = 1.0d0

      nperc = 3
      nthresh(1) = 75
      nthresh(2) = 95
      nthresh(3) = 99

      do l=1,nseason
        if (zhang) then
c do the Zhang et al. approach: leave one year out
          do m=calyrbeg,calyrend
            ndata = 0
            do i=calyrbeg,calyrend
              if(m.eq.i) goto 120
              do j=1,nmonths(l)
                call lengthofmonth(1973,month(l,j),length)
                do k=1,length
                  if((qc(i,month(l,j),k).eq.0).and.
     +               (a(i,month(l,j),k).gt.thresh)) then
                    ndata = ndata + 1
                    xdata(ndata) = a(i,month(l,j),k)
                  endif
                enddo
              enddo
 120          continue

              if(.not.fancy) then
c use the "poor men's" percentile estimate: sort data and take the relevant percentiles
                call poormenpercRR(nperc,nthresh,percentiles)
              else
c use the "fancy men's" percentile estimate: fit a two-parameter Gamma function
                call gammafitperc(nperc,nthresh,percentiles)
              endif

            enddo
            p75(m,l) = percentiles(1)
            p95(m,l) = percentiles(2)
            p99(m,l) = percentiles(3)
          enddo
        endif

c now calculate percentiles over the complete calibration period: no years left out
        m = calyrend+1
        ndata = 0
        do i=calyrbeg,calyrend
          do j=1,nmonths(l)
            call lengthofmonth(1973,month(l,j),length)
            do k=1,length
              if((qc(i,month(l,j),k).eq.0).and.
     +           (a(i,month(l,j),k).gt.thresh)) then
                ndata = ndata + 1
                xdata(ndata) = a(i,month(l,j),k)
              endif
            enddo
          enddo
        enddo

        if(.not.fancy) then
c use the "poor men's" percentile estimate: sort data and take the percentiles
          call poormenpercRR(nperc,nthresh,percentiles)
        else
c use the "fancy men's" percentile estimate: fit a two-parameter Gamma function
          call gammafitperc(nperc,nthresh,percentiles)
        endif

        p75(m,l) = percentiles(1)
        p95(m,l) = percentiles(2)
        p99(m,l) = percentiles(3)

      enddo

c if the Zhang et al. approach is NOT followed: then replace p10(calyrbeg:calyrend)
c with the value in p10(calyrend+1)
      if(.not.zhang) then
        do m=calyrbeg,calyrend
          do l=1,nseason
            p75(m,l) = p75(calyrend+1,l)
            p95(m,l) = p95(calyrend+1,l)
            p99(m,l) = p99(calyrend+1,l)
          enddo
        enddo
      endif

c     write(6,123) (p75(calyrend+1,l), l=1,nseason)
c     write(6,123) (p95(calyrend+1,l), l=1,nseason)
c     write(6,123) (p99(calyrend+1,l), l=1,nseason)
 123  format(7e14.4)

      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine poormenpercRR(nperc,nthresh,percentiles)
c this routine calculates the percentile in nthresh of the distribution based
c on the data in the array xdum
c
c here the simple approach is used
      implicit none
      include 'comgeneral.h'

      integer       mxlength,mxperc
      parameter     (mxlength=(calyrend-calyrbeg+1)*365,mxperc=5)
      integer       N,nperc,absen
      real*8        absentr,absenr
      real*8        xdata(mxlength)
      integer       i,j,k,nthresh(mxperc),nthreshd(mxperc)
      real*8        percentiles(nperc)
c NAG things
      integer       ifail,irank(mxlength)
      external      m01daf

      common/datablock/xdata
      common/intblock/N

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
c         write(6,*) j,nthreshd(j),percentiles(j)
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
      subroutine gammafitperc(nperc,nthresh,percentiles)
c this routine calculates the percentile of the distribution based
c on the data in the array xdum
c
c here we fit a Gamma distribution to the data using maximum-likelihood
c theory is in Wilks p89
      implicit none
      include 'comgeneral.h'

      integer       mxlength
      parameter     (mxlength=(calyrend-calyrbeg+1)*365)
      integer       N,nperc,absen
      real*8        absentr,absenr
      real*8        xdata(mxlength),array(mxlength)
      integer       i,j,k,nthresh(nperc)
      real*8        percentiles(nperc)
c NAG things
      integer       ifail
      external      g01aaf,g08cgf,g01aef,g01fff
      real*8        g01fff,xmx,xmn

      integer       nclass
      parameter     (nclass=100)
      real*8        xmin,xmax,s2,s3,s4,wtsum,wt(mxlength),chisq
      real*8        eps,t,a,b
      real*8        cint(nclass),par(2),prob(nclass),eval(nclass)
      real*8        chisqi(nclass)
      integer       ifreq(nclass)
      integer       maxcal,ndf
      real*8        xd,x,alpha,beta,D,offset,xmean,p,xthresh
      logical       test

      common/datablock/xdata
      common/intblock/N

      absen = absent + 1
      absenr = dble(absen)
      absentr = dble(absent)

      if(N.ge.pdayspresent) then
        test = .false.

c calculate the mean of the distribution and test for positive skewness
        ifail = 0
        call g01aaf(N,xdata,0,wt,xmean,s2,s3,s4,xmin,xmax,wtsum,ifail)

c add a skewness test
c     if(s3.lt.0.0) call error('data should be positively skewed')

        xd = 0.0D0
        do i=1,N
          xd = xd + log(xdata(i))
        enddo
        D = log(xmean) - xd/dble(N)
        alpha = (1.0 + sqrt(1.0 + 4*D/3.0))/(4.0*D)
        beta = xmean/alpha

c this offers the possibility of making a goodness-of-fit test
        if(test) then
          par(1) = alpha
          par(2) = beta

c construct the clas boundaries
          do i=1,nclass
            cint(i) = (i-1)*1.0d0
          enddo

c do the goodness-of-fit test
          ifail = 1
          call g08cgf(nclass,ifreq,cint,'G',par,0,prob,chisq,p,ndf,
     +            eval,chisqi,ifail)

          write(6,62) chisqi

c if ifail = 10, then we can go on
          if((ifail.ne.0).and.(ifail.ne.10)) then
            write(6,63) ifail
            stop
          endif
        endif

c this calculates the percentiles
        do i=1,nperc
          xthresh = dble(nthresh(i))/100.0d0
          ifail = 0
          percentiles(i) = g01fff(xthresh,alpha,beta,0.0,ifail)
        enddo
      else
        do i=1,nperc
          percentiles(i) = absentr
        enddo
      endif

 62   format('Chi squared: ',e14.4)
 63   format('error in g08cgf- ifail is: ',I4)

      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine fillValue(perc)
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
c
c---------------------------------------------------------------------------
c
      integer function nmonths(i)
c this function is used in the percentile calculations for RR
c it returns the number of months in season i (i=1,..,nseason)
c
      implicit none
      include 'comgeneral.h'

      integer       i,j
      integer       sl(nseason)
      data          (sl(j), j=1,nseason) / 
     +         12,6,6,3,3,3,3,1,1,1,1,1,1,1,1,1,1,1,1/

      if(i.gt.nseason) call error('out of bounds: function nmonths')

      nmonths = sl(i)

      return
      end
c
c---------------------------------------------------------------------------
c
      integer function month(season,j)
c this function is used in the percentile calculations for RR
c given the season, it returns the month number
c
      implicit none
      include 'comgeneral.h'

      integer       j,season,ndum

      if(season.le.4) then
        if(season.eq.1) then
          month = j
        elseif(season.eq.2)then
          ndum = mod(j+9,12)
          if(ndum.eq.0) ndum = 12
          month = ndum
        elseif(season.eq.3)then
          month = j+3
        else
          ndum = mod(j+11,12)
          if(ndum.eq.0) ndum = 12
          month = ndum
        endif
      else
        if(season.eq.5)then
          month = j+2
        elseif(season.eq.6)then
          month = j+5
        elseif(season.eq.7)then
          month = j+8
        else
          month = season - 7
        endif
      endif

      return
      end 
c
c---------------------------------------------------------------------------
c
