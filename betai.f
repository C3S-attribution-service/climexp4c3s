C  (C) Copr. 1986-92 Numerical Recipes Software +.-).
      FUNCTION betai(a,b,x)
      REAL betai,a,b,x
CU    USES betacf,gammln
      REAL bt,betacf,gammln
      if(.not.(x.ge.0..and.x.le.1.))then
        print *,'betai: input a,b,x = ',a,b,x
        print *,'bad argument x in betai'
***        betai = 0./0.
        return
      endif
      if(x.eq.0..or.x.eq.1.)then
        bt=0.
      else
        bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.-x))
      endif
      if(x.lt.(a+1.)/(a+b+2.))then
        betai=bt*betacf(a,b,x)/a
        return
      else
        betai=1.-bt*betacf(b,a,1.-x)/b
        return
      endif
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +.-).
      FUNCTION betacf(a,b,x)
      INTEGER MAXIT
      REAL betacf,a,b,x,EPS,FPMIN
      PARAMETER (MAXIT=200,EPS=3.e-7,FPMIN=1.e-30)
      INTEGER m,m2
      REAL aa,c,d,del,h,qab,qam,qap
      qab=a+b
      qap=a+1.
      qam=a-1.
      c=1.
      d=1.-qab*x/qap
      if(abs(d).lt.FPMIN)d=FPMIN
      d=1./d
      h=d
      do 11 m=1,MAXIT
        m2=2*m
        aa=m*(b-m)*x/((qam+m2)*(a+m2))
        d=1.+aa*d
        if(abs(d).lt.FPMIN)d=FPMIN
        c=1.+aa/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        h=h*d*c
        aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
        d=1.+aa*d
        if(abs(d).lt.FPMIN)d=FPMIN
        c=1.+aa/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1
11    continue
        print *,'betacf: called with a,b,x = ',a,b,x
        print *,'        abs(del-1) = ',abs(del-1),eps
        print *,'a or b too big, or MAXIT too small in betacf'
        betacf = 3e33
        return
1     betacf=h
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +.-).
