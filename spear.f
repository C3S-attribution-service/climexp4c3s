C  (C) Copr. 1986-92 Numerical Recipes Software +.-).
      SUBROUTINE spearx(data1,data2,n,wksp1,wksp2,d,zd,probd,rs,probrs,
     +        sum,a1,s1,a2,s2)
      INTEGER n,nover
      REAL d,probd,probrs,sum,rs,zd,data1(n),data2(n),wksp1(n),wksp2(n)
      real a1,s1,a2,s2
CU    USES betai,crank,erfcc,sort2
      INTEGER j
      REAL aved,df,en,en3n,fac,sf,sg,t,vard,betai,erfcc
	logical lcheck
	parameter(lcheck=.false.)
      do 11 j=1,n
        wksp1(j)=data1(j)
        wksp2(j)=data2(j)
11    continue
      call sort2(n,wksp1,wksp2)
	if ( lcheck) then
	    do j=1,n-1
	    	if ( wksp1(j).gt.wksp1(j+1) ) then
		    print *,'spearx: error: unsorted array: ',
     +			wksp1(j),wksp1(j+1)
		endif
	    enddo
	endif
      if ( mod(n,2).eq.0 ) then
          a1 = (wksp1(n/2) + wksp1(n/2+1))/2
          if ( mod(n,4).eq.0 ) then
              s1 = (wksp1(3*n/4) + 3*wksp1(3*n/4+1) 
     +              - 3*wksp1(n/4) - wksp1(n/4+1))/4
          else
              s1 = (3*wksp1((3*n+2)/4) + wksp1((3*n+2)/4+1) 
     +              - wksp1((n-2)/4) - 3*wksp1((n-2)/4+1))/4
          endif
      else
          a1 = wksp1((n+1)/2)
          if ( mod(n+1,4).eq.0 ) then
              s1 = wksp1(3*(n+1)/4) - wksp1((n+1)/4)
          else
              s1 = (wksp1((3*n+1)/4) + wksp1((3*n+1)/4+1) 
     +              - wksp1((n-1)/4) - wksp1((n-1)/4+1))/2
          endif
      endif
      s1 = n*s1**2
      call crank(n,wksp1,sf)
      call sort2(n,wksp2,wksp1)
	if ( lcheck ) then
	    do j=1,n-1
	    	if ( wksp2(j).gt.wksp2(j+1) ) then
		    print *,'spearx: error: unsorted array2: ',
     +			wksp2(j),wksp2(j+1)
	    	endif
	    enddo
	endif
      if ( mod(n,2).eq.0 ) then
          if ( lcheck ) print *,'mod(n,2).eq.0'
          a2 = (wksp2(n/2) + wksp2(n/2+1))/2
          if ( mod(n,4).eq.0 ) then
              s2 = (wksp2(3*n/4) + 3*wksp2(3*n/4+1) 
     +              - 3*wksp2(n/4) - wksp2(n/4+1))/4
          else
              s2 = (3*wksp2((3*n+2)/4) + wksp2((3*n+2)/4+1) 
     +              - wksp2((n-2)/4) - 3*wksp2((n-2)/4+1))/4
          endif
      else
          a2 = wksp2((n+1)/2)
          if ( mod(n+1,4).eq.0 ) then
              s2 = wksp2(3*(n+1)/4) - wksp2((n+1)/4)
          else
              s2 = (wksp2((3*n+1)/4) + wksp2((3*n+1)/4+1) 
     +              - wksp2((n-1)/4) - wksp2((n-1)/4+1))/2
          endif
      endif
      s2 = n*s2**2
      call crank(n,wksp2,sg)
      d=0.
      do 12 j=1,n
        d=d+(wksp1(j)-wksp2(j))**2
12    continue
      en=n
      en3n=en**3-en
      aved=en3n/6.-(sf+sg)/12.
      fac=(1.-sf/en3n)*(1.-sg/en3n)
      vard=((en-1.)*en**2*(en+1.)**2/36.)*fac
      zd=(d-aved)/sqrt(vard)
	if ( sum.ge.1e33 ) then
	    probd = 3e33
	else
	    probd=erfcc(abs(zd)/1.4142136)
	endif
	if ( fac.eq.0 ) then
	    rs = 3e33
	    probrs = 3e33
	    return
	else
      	    rs=(1.-(6./en3n)*(d+(sf+sg)/12.))/sqrt(fac)
	endif
	if ( lcheck ) then
	    print *,'spear: d =    ',d
	    print *,'spear: en =   ',en
	    print *,'spear: en3n = ',en3n
	    print *,'spear: sf =   ',sf
	    print *,'spear: sg =   ',sg
	    print *,'spear: fac =  ',fac
	    print *,'spear: rs =   ',rs
	endif
      fac=(1.+rs)*(1.-rs)
      if(fac.gt.0.)then
        t=rs*sqrt((en-2.)/fac)
        df=en/sum-2.
        probrs=betai(0.5*df,0.5,df/(df+t**2))
      else
        probrs=0.
      endif
      return
      END
