
      program gauss
c-----------------------------------------------------------------------
c low-pass filters data; input and output as a grads file
c-----------------------------------------------------------------------
      parameter(nlat=32,nlon=64,itot=999)
      dimension te1(nlat,nlon,itot),te2(nlat,nlon,itot),
     *  te3(nlat,nlon,itot),te4(nlat,nlon,itot),te5(nlat,nlon,itot)
      dimension te1g(nlat,nlon,itot),te2g(nlat,nlon,itot),  
     *  te3g(nlat,nlon,itot),te4g(nlat,nlon,itot),te5g(nlat,nlon,itot)     
c
      open(52,file='atmos/atmc_seas.dat',form='unformatted',
     *                        access='direct',recl=32*64)
      open(54,file='atmos/atmc_g50.dat',form='unformatted',
     *                        access='direct',recl=32*64)
c
c read seasonal data
c
      do it=1,itot
        read(52) ((te1(i,j,it),j=1,nlon), i=1,nlat)
        read(52) ((te2(i,j,it),j=1,nlon), i=1,nlat)
        read(52) ((te3(i,j,it),j=1,nlon), i=1,nlat)
        read(52) ((te4(i,j,it),j=1,nlon), i=1,nlat)
          do i=1,nlat
            do j=1,nlon
              te5(i,j,it)=(te1(i,j,it)+te2(i,j,it)+
     *             te3(i,j,it)+te4(i,j,it))/4.
            enddo             
          enddo
      enddo     
c
c apply iw-year Gaussian filter to sT, wT and aT
c iw can be 10, 25, 50, 100 years
c
      iw=50
      iy1=1
      iy2=999
      write(6,*) iy1,iy2, iy2-iy1
      call filter(iw,iy1,iy2,te1,iw2,te1g)
      call filter(iw,iy1,iy2,te3,iw2,te3g)
      call filter(iw,iy1,iy2,te5,iw2,te5g)
      write(6,*) 'G filter', iw,iw2
      write(6,*) iy1+iw2,iy2-iw2, iy2-iy1-2*iw2
c
      do it=iy1+iw2,iy2-iw2
        write(54) ((te1g(i,j,it),j=1,nlon), i=1,nlat)
        write(54) ((te3g(i,j,it),j=1,nlon), i=1,nlat)
        write(54) ((te5g(i,j,it),j=1,nlon), i=1,nlat)
      enddo     

      stop
      end
c
c====================================================================
      subroutine filter(iw,iy1,iy2,te,iw2,teg)
c Applies a iw-year gaussian filter to data (help) covering the 
c years iy1 to iy2. The filtered data (helpsum) cover iy2+iw2 to iy2-iw2
c Gaussian weights are assigned to the array w; iw can be 10,25,50,100 years
c
      parameter(nlat=32,nlon=64,nh=3,itot=999)
      dimension te(nlat,nlon,itot),teg(nlat,nlon,itot),
     *        sum(nlat,nlon),w(-36:36)
c
      if (iw.eq.10) then
        iv=9
        w(0)=0.2408    
        w(1)=0.2011   
        w(2)=0.1172   
        w(3)=0.0477   
        w(4)=0.0136
        do j=1,4
          w(-j)=w(j)
        enddo
      endif
c
      if (iw.eq.25) then
        iv=19
        w(0)=0.0979    
        w(1)=0.0951    
        w(2)=0.0873   
        w(3)=0.0756    
        w(4)=0.0618    
        w(5)=0.0477    
        w(6)=0.0347    
        w(7)=0.0239    
        w(8)=0.0155    
        w(9)=0.0095    
        do j=1,9
          w(-j)=w(j)
        enddo
      endif
c
      if (iw.eq.50) then
        iv=37
        w(0)=0.0492
        w(1)=0.0488
        w(2)=0.0478
        w(3)=0.0461
        w(4)=0.0438
        w(5)=0.0411
        w(6)=0.0379
        w(7)=0.0345
        w(8)=0.0310
        w(9)=0.0274
        w(10)=0.0239
        w(11)=0.0206
        w(12)=0.0174
        w(13)=0.0146
        w(14)=0.0120
        w(15)=0.0097
        w(16)=0.0078
        w(17)=0.0061
        w(18)=0.0049
        do j=1,18
          w(-j)=w(j)
        enddo
      endif
c
      if (iw.eq.100) then
        iv=73
        w(0)=0.02464
        w(1)=0.02459
        w(2)=0.02446
        w(3)=0.02424
        w(4)=0.02394
        w(5)=0.02356
        w(6)=0.02309
        w(7)=0.02256
        w(8)=0.02196
        w(9)=0.02129
        w(10)=0.02058
        w(11)=0.01982
        w(12)=0.01902
        w(13)=0.01818
        w(14)=0.01731
        w(15)=0.01643
        w(16)=0.01554
        w(17)=0.01464
        w(18)=0.01375
        w(19)=0.01286
        w(20)=0.01199
        w(21)=0.01114
        w(22)=0.01031
        w(23)=0.00950
        w(24)=0.00874
        w(25)=0.00800
        w(26)=0.00730
        w(27)=0.00663
        w(28)=0.00601
        w(29)=0.00542
        w(30)=0.00488
        w(31)=0.00437
        w(32)=0.00390
        w(33)=0.00347
        w(34)=0.00308
        w(35)=0.00272
        w(36)=0.00239
        do j=1,36
          w(-j)=w(j)
        enddo
      endif

      iv2=1+(iv-1)/2
      do k=iy1,iy2-(iv-1)
        do i=1,nlat
        do j=1,nlon
          sum(i,j)=0.
          do m=1,iv
            sum(i,j)=sum(i,j)+w(m-iv2)*te(i,j,k+m-1)
          enddo
          teg(i,j,k+iv2-1)=sum(i,j)
        enddo
        enddo
      enddo
      iw2=iv2-1
c
      return
      end
