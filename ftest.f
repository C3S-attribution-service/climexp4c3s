C  (C) Copr. 1986-92 Numerical Recipes Software +.-).
      SUBROUTINE ftestx(data1,n1,a1,data2,n2,a2,f,prob)
        ! GJvO: added lag-1 autocorrelations a1,a2 to decrease the number of degrees of freedom
        ! in case of significant autocorrelations.
        implicit none
      INTEGER n1,n2
      REAL f,prob,data1(n1),data2(n2),a1,a2
CU    USES avevar,betai
      REAL ave1,ave2,df1,df2,var1,var2,betai,s
      call avevar(data1,n1,ave1,var1)
      call avevar(data2,n2,ave2,var2)
        if ( a1.lt.1/sqrt(real(n1)) .and. a1.gt.exp(-1.) ) then
            df1=n1-1
        else
            df1 = (n1-1)*(-log(a1))
        end if
        if ( a2.lt.1/sqrt(real(n2)) .and. a2.gt.exp(-1.) ) then
            df2=n2-1
        else
            df2 = (n2-1)*(-log(a2))
        end if
      if(var1.gt.var2)then
        f=var1/var2
      else
        f=var2/var1
        s = df1
        df1=df2
        df2=s
      endif
      prob=2.*betai(0.5*df2,0.5*df1,df2/(df2+df1*f))
      if(prob.gt.1.)prob=2.-prob
      return
      END
