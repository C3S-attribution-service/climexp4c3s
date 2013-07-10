      FUNCTION dgammln(xx)
*       hacked up gammln (NUmerical recipes) to return dlog(Gamma(x))/dx
        implicit none
      REAL dgammln,xx
      INTEGER j
      DOUBLE PRECISION ser,dser,stp,tmp,x,y,cof(6)
      double precision xgammln
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=log(tmp)+(x+.5d0)/tmp-1
      ser=1.000000000190015d0
      dser = 0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
        dser = dser - cof(j)/y**2
11    continue
      dgammln=tmp+dser/ser-1/x
***        print *,'dgammln  ',dgammln
***        print *,'     cmp ',real((xgammln(x+1d-7)-xgammln(x-1d-7))/2d-7)
      return
      END

      double precision FUNCTION xgammln(x)
        implicit none
      DOUBLE PRECISION x
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      xgammln=tmp+log(stp*ser/x)
      return
      END
