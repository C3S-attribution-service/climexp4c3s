C  (C) Copr. 1986-92 Numerical Recipes Software +.-).
        SUBROUTINE moment(data,n,ave,adev,sdev,var,skew,curt)
        implicit none
        INTEGER n
        REAL adev,ave,curt,sdev,skew,var,data(n)
        INTEGER j
        REAL*8 p,s,ep,rvar,rskew,rcurt
        logical lwrite
        parameter (lwrite=.false.)
        if ( lwrite ) then
            print *,'entering moment, n = ',n
            print *,'data = ',(data(j),j=1,n)
        endif
        if (n.le.1)then
            write(0,*) 'n must be at least 2 in moment'
            ave = 3e33
            adev = 3e33
            sdev = 3e33
            var = 3e33
            skew = 3e33
            curt = 3e33
            return
        endif
        s=0.
        do 11 j=1,n
            s=s+data(j)
   11   continue
        ave=s/n
        adev=0.
        rvar=0.
        rskew=0.
        rcurt=0.
        ep=0.
        do 12 j=1,n
            s=data(j)-ave
            ep=ep+s
            adev=adev+abs(s)
            p=s*s
            rvar=rvar+p
            p=p*s
            rskew=rskew+p
            p=p*s
            rcurt=rcurt+p
   12   continue
        adev=adev/n
        var=(rvar-ep**2/n)/(n-1)
        if ( var.le.0 ) then
            sdev = 0
        else
            sdev=sqrt(var)
        endif
        if ( lwrite ) print *,'adev,var,sdev = ',adev,var,skew
        if(sdev**3.ne.0.)then
            skew=rskew/(n*sdev**3)
        else
***     print *,'moment: warning: found sdev~0 ',n,data
            skew = 3e33
        endif
        if(dble(var)**2.ne.0.)then
            curt=dble(rcurt)/(n*dble(var)**2)-3.
        else
***     print *,'moment: warning: found var~0 ',n,data
            curt = 3e33
        endif
        if ( lwrite ) print *,'skew,curt     = ',skew,curt
        return
        END
