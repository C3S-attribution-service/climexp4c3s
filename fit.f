C  (C) Copr. 1986-92 Numerical Recipes Software +.-).
      SUBROUTINE fit(x,y,ndata,sig,mwt,a,b,siga,sigb,chi2,q)
*     fits a+b*x to (x,y)
      INTEGER mwt,ndata
      REAL a,b,chi2,q,siga,sigb,sig(ndata),x(ndata),y(ndata)
CU    USES gammq
      INTEGER i
      REAL sigdat,ss,st,st2,sx,sxoss,sy,syoss,t,wt,gammq
*       special case causes crashes later on...
        if ( ndata.eq.1 ) then
            a = y(1)
            b = 0
            siga = 3e33
            sigb = 3e33
            chi2 = 0
            q = 3e33
        endif
      sx=0.
      sy=0.
      st=0.
      st2=0.
      b=0.
      if(mwt.ne.0) then
        ss=0.
        do 11 i=1,ndata
          wt=1./(sig(i)**2)
          ss=ss+wt
          sx=sx+x(i)*wt
          sy=sy+y(i)*wt
11      continue
      else
        do 12 i=1,ndata
          sx=sx+x(i)
          sy=sy+y(i)
12      continue
        ss=float(ndata)
      endif
      sxoss=sx/ss
      syoss=sy/ss
      if(mwt.ne.0) then
        do 13 i=1,ndata
          t=(x(i)-sxoss)/sig(i)
          st=st+t
          st2=st2+t*t
          b=b+t*y(i)/sig(i)
13      continue
      else
        do 14 i=1,ndata
          t=x(i)-sxoss
          st=st+t
          st2=st2+t*t
          b=b+t*y(i)
14      continue
      endif
!GJvO
      if ( st2.eq.0 ) then
        b = 0
        sigb = 0
        a = sy/ss
        siga = 1/sqrt(ss)
      else
!GJvO
!!!          print *,b,syoss*st,st
      b=b/st2
      a=(sy-sx*b)/ss
      siga=sqrt((1.+sx*sx/(ss*st2))/ss)
      sigb=sqrt(1./st2)
      endif
      chi2=0.
        if ( ndata.eq.2 ) then
            siga = 3e33
            sigb = 3e33
            q = 1
        else
            if(mwt.eq.0) then
                do 15 i=1,ndata
                    chi2=chi2+(y(i)-a-b*x(i))**2
   15           continue
                q=1.
                sigdat=sqrt(chi2/(ndata-2))
                siga=siga*sigdat
                sigb=sigb*sigdat
            else
                do 16 i=1,ndata
                    chi2=chi2+((y(i)-a-b*x(i))/sig(i))**2
   16           continue
                q=gammq(0.5*(ndata-2),0.5*chi2)
            endif
        endif
      return
      END
