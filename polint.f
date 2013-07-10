C  (C) Copr. 1986-92 Numerical Recipes Software +.-).
      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      integer nerror
      save nerror
      data nerror /0/
      ns=1
      dif=abs(x-xa(1))
        if ( n.gt.NMAX ) then
            do iu=0,6,6
                write(iu,*) 'polint: increse NMAX from ',NMAX,' to ',n
            enddo
            call abort
        endif
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)then
              if ( nerror.eq.0 ) then
                  do iu=0,6,6
                      write(iu,*) 'failure in polint, den = ',den
                      write(iu,*)' x,y,dy = ',x,y,dy
                      do k=1,n
                          write(iu,*) k,xa(k),ya(k)
                      enddo
                  enddo
              endif
              nerror = nerror + 1
              y = 3e33
              dy= 3e33
              return
          endif
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
