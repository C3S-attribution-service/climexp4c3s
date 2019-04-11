subroutine spearx(data1,data2,n,wksp1,wksp2,d,zd,probd,rs,probrs,sum,a1,s1,a2,s2)
!
!   wrapper for extended spear call used in Climate Explorer
!
    use fgsl
    implicit none
    integer,intent(in) :: n
    real,intent(in) :: data1(n),data2(n),wksp1(1),wksp2(1),sum
    real,intent(out) :: d,zd,probd,rs,probrs,a1,s1,a2,s2
    integer(fgsl_size_t) :: dim
    real(fgsl_double),allocatable :: ddata1(:),ddata2(:),work(:)
    real(fgsl_double) :: drs,fac,t,arg1,arg2,arg3,dt,df,da1,ds1,da2,ds2
!
!   most output is not used, set to undef
    d = 3e33
    zd = 3e33
    probd = 3e33
    rs = 3e33
    probrs = 3e33
    a1 = 3e33
    s1 = 3e33
    a2 = 3e33
    s2 = 3e33
    if ( n <= 2 ) return
!
!   rank correlation coefficient
!
    dim = n
    allocate(ddata1(n),ddata2(n),work(2*n))
    ddata1 = data1
    ddata2 = data2
    drs = fgsl_stats_spearman(ddata1,1_fgsl_size_t,ddata2,1_fgsl_size_t,dim,work)
    da1 = fgsl_stats_mean(ddata1,1_fgsl_size_t,dim)
    ds1 = fgsl_stats_variance_m(ddata1,1_fgsl_size_t,dim,da1)
    da2 = fgsl_stats_mean(ddata2,1_fgsl_size_t,dim)
    ds2 = fgsl_stats_variance_m(ddata2,1_fgsl_size_t,dim,da2)
    deallocate(ddata1,ddata2,work)
    rs = drs
    a1 = da1
    s1 = sqrt(ds1)
    a2 = da2
    s2 = sqrt(ds2)
    fac = (1-drs)*(1+drs)
    if ( sum < 1 ) then
        write(0,*) 'spearx: error: sum < 1 ',sum
        call exit(-1)
    end if
    if ( fac == 0 ) then ! r=±1
        probrs = 0
    else if ( n > 100 ) then
        dt = drs*sqrt((n/sum-2)/fac)
        probrs = fgsl_sf_erfc(dt/sqrt(2d0))    
    else ! Eq 14.6.2 in the blue book
        dt = drs*sqrt((n-2)/fac)
        df = n/sum - 2. ! sum represents how much samples are dependent
        arg1 = df/2
        arg2 = 0.5d0
        arg3 = df/(df+dt**2)
        probrs = fgsl_sf_beta_inc(arg1,arg2,arg3)
    end if
end subroutine spearx