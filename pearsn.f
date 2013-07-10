C  (C) Copr. 1986-92 Numerical Recipes Software +.-).
      SUBROUTINE pearsnxx(x,y,n,r,prob,z,ax,sxx,ay,syy,sxy,df)
      implicit none
      INTEGER n
      REAL prob,r,z,x(n),y(n),TINY
      PARAMETER (TINY=1.e-20)
CU    USES betai
      INTEGER j
      REAL ax,ay,df,sxx,sxy,syy,t,xt,yt,betai,arg3
        integer ierr,nerr
        save ierr,nerr
        data ierr /0/
        data nerr /1/
      ax=0.
      ay=0.
      do 11 j=1,n
        ax=ax+x(j)
        ay=ay+y(j)
11    continue
      ax=ax/n
      ay=ay/n
      sxx=0.
      syy=0.
      sxy=0.
      do 12 j=1,n
        xt=x(j)-ax
        if ( abs(xt).lt.1e-6*abs(ax) ) xt = 0
        yt=y(j)-ay
        if ( abs(yt).lt.1e-6*abs(ay) ) yt = 0
        sxx=sxx+xt**2
        syy=syy+yt**2
        sxy=sxy+xt*yt
12    continue
      if ( sxx*syy.eq.0 ) then
	r=0
      else
        r=sxy/sqrt(sxx*syy)
***        if ( abs(r).lt.1e-7 ) then
***            print *,'Found incredibly small correlation: ',r
***            print *,'sxx,syy,sxy = ',sxx,syy,sxy
***            print *,'ax,ay = ',ax,ay
***            do j=1,n
***                print *,j,x(j),x(j)-ax,y(j),y(j)-ay
***            enddo
***        endif
      endif
      if ( r.gt.1 ) r = 1
      if ( r.lt.-1 ) r = -1
      if ( abs(r-1.).lt.1e-6 ) then
        z = 3e33
        prob = 0
        return
      endif
      z = ((1.+r)+TINY)/((1.-r)+TINY)
      if ( z.gt.0 ) then
        z=0.5*log(z)
      else
        write(0,*) 'pearsnx: error: arg log<=0: ',z
        write(0,*) '                r = ',r
        z = 3e33
      endif
      t = (((1.-r)+TINY)*((1.+r)+TINY))
      if ( .not. df.gt.0 ) then
          ierr = ierr + 1
          if ( ierr.ge.nerr ) then
              nerr = 2*nerr
              write(*,*) 'pearsn: df<=0: ',df
              write(0,*) 'pearsn: df<=0: ',df
          end if
          prob = 3e33
        return
      endif
      if ( t.gt.0 ) then
        t=r*sqrt(df/t)
      else
        write(0,*) 'pearsnx: error: arg sqrt<=0: ',t
        write(0,*) '                r = ',r
        t = 3e33
      endif
      arg3 = df/(df+t**2)
      if ( .not. (arg3.ge.0..and.arg3.le.1.) ) then
        print *,'pearsnx: error: arg betai = ',df/(df+t**2)
	print *,'                df = ',df
	print *,'                r  = ',r
	print *,'                t  = ',t
	do j=1,n
	  print *,j,x(j),x(j)-ax,y(j),y(j)-ay
	enddo
	return
      endif
      prob=betai(0.5*df,0.5,arg3)
C     prob=erfcc(abs(z*sqrt(n-1.))/1.4142136)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +.-).
      SUBROUTINE apearsnxx(x,y,n,r,prob,z,ax,sxx,ay,syy,sxy,df)
*	no mean subtracted
      implicit none
      INTEGER n
      REAL prob,r,z,x(n),y(n),TINY
      PARAMETER (TINY=1.e-20)
CU    USES betai
      INTEGER j
      REAL ax,ay,df,sxx,sxy,syy,t,xt,yt,betai,arg3
      if ( .not. df.gt.0 ) then
        write(0,*) 'pearsn: error: df<=0: ',df
        call abort
      endif
      ax=0.
      ay=0.
      sxx=0.
      syy=0.
      sxy=0.
      do 12 j=1,n
        xt=x(j)-ax
        yt=y(j)-ay
        sxx=sxx+xt**2
        syy=syy+yt**2
        sxy=sxy+xt*yt
12    continue
      if ( sxx*syy.eq.0 ) then
	r=0
      else
        r=sxy/sqrt(sxx*syy)
      endif
      if ( abs(r-1.).lt.1e-6 ) then
        z = 3e33
        prob = 0
        return
      endif
      z=0.5*log(((1.+r)+TINY)/((1.-r)+TINY))
      t=r*sqrt(df/(((1.-r)+TINY)*((1.+r)+TINY)))
      arg3 = df/(df+t**2)
      if ( .not. (arg3.ge.0..and.arg3.le.1.) ) then
        print *,'pearsnx: error: arg betai = ',df/(df+t**2)
	print *,'                df = ',df
	print *,'                r  = ',r
	print *,'                t  = ',t
	do j=1,n
	  print *,j,x(j),x(j)-ax,y(j),y(j)-ay
	enddo
	return
      endif
      prob=betai(0.5*df,0.5,arg3)
C     prob=erfcc(abs(z*sqrt(n-1.))/1.4142136)
      return
      END

