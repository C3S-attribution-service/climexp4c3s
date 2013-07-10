C  (C) Copr. 1986-92 Numerical Recipes Software +.-).
      SUBROUTINE chsone(bins,ebins,nbins,knstrn,df,chsq,prob)
      INTEGER knstrn,nbins
      REAL chsq,df,prob,bins(nbins),ebins(nbins)
CU    USES gammq
      INTEGER j
      REAL gammq
      df=nbins-knstrn
      chsq=0.
      do 11 j=1,nbins
          if(ebins(j).le.0.)then
              if ( bins(j).ne.0 ) then
                  write(0,*) 'chsone: ebins(',j,') <= 0 ',ebins(j)
                  write(*,*) 'chsone: ebins(',j,') <= 0 ',ebins(j)
                  call abort
              endif
          else
              chsq=chsq+(bins(j)-ebins(j))**2/ebins(j)
          endif
11    continue
      prob=gammq(0.5*df,0.5*chsq)
      return
      END
