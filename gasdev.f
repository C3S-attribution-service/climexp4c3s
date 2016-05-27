C  (C) Copr. 1986-92 Numerical Recipes Software +.-).
      FUNCTION gasdev(iseed)
      REAL gasdev
CU    USES ran1
      INTEGER iset,nerr,ntot
      REAL fac,gset,rsq,v1,v2
      SAVE iset,gset,ntot,nerr
      DATA iset,nerr,ntot/0,0,0/
      if (iset.eq.0) then
1       continue
        call random_number(v1)
        v1=2*v1-1
        call random_number(v2)
        v2=2*v2-1
        ntot = ntot + 1
        if ( v1.eq.v2 ) then
            nerr = nerr + 1
            if ( nerr.gt.ntot/10 ) then
                write(0,*) 'gasdev: error: compiler optimisations'
            endif
        endif
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END
