subroutine printcorr(dindx,ddata,lfirst,dddata,yrmo,n,n0,j1,j2 &
    ,month,nperyear,lag,string,lboot,lprint,result,dresult,prob)

!   actually perform the correlations based on the linear ararys
!   and print the result (if lprint)

    implicit none
    include 'param.inc'
    include 'getopts.inc'
    integer :: n,n0,month,nperyear,lag,j1,j2,yrmo(3,n)
    real :: dindx(n),ddata(n),dddata(n),result,dresult(-2:2),prob
    logical :: lfirst(n)
    character*(*) :: string
    logical :: lboot,lprint
    integer :: i,j,nunit,ilen
    integer :: iindx(ndata),ks,ks1,ks2,l,n1,ndecor
    integer :: ncont(0:3,0:3)
    real :: zd,pt,adat1,avar1,adat2,avar2,d,rmse,mae,a1,a2
    real :: chi2,xn,x(ndata),sig(ndata),u(ndata,2),v(2,2),w(2),da(2,2),aa(ndata),bb(ndata)
    real :: df,sum,a(2),q,r,z,adata,sxx,aindx,syy,sxy,probd,bias,xneff,chieff
    real,allocatable :: bb1(:,:),wdata(:),windx(:),w1(:),w2(:)
    real :: gammq
    external gammq,findx

    if ( n > ndata ) then
        write(0,*) 'getcorr: increase ndata to ',n
        write(*,*) 'getcorr: increase ndata to ',n
        call exit(-1)
    end if

    result = 3e33
    if ( lks ) then

!       Kolmogorov-Smirnov test, NOT PORTED YET, NOT USED IN CLIMEXP

        if ( kscut < -1e33 ) then
            ks1 = 5
            ks2 = n-4
        else
            ks1 = 0
            ks2 = 0
        end if
        call ffsort(dindx,iindx,n)
        do ks=ks1,ks2
!           NumRec vernaggelt (sorteert) dddata...
            do l=1,n
                dddata(l) = ddata(iindx(l))
            end do
            if ( ks1 /= 0 ) then
                kscut = (dindx(iindx(ks-1)) + dindx(iindx(ks)))/2
            end if
            do l=1,n
                if ( dindx(iindx(l)) > kscut ) goto 730
            end do
        730 n1 = l-1
            if ( lwrite ) then
                print *,'Found ',n1,' points below ',kscut
                if ( n1 < 20 ) print '(f10.2)',(dddata(l),l=1,n1)
                print *,'Found ',n-n1,' points above ',kscut
                if ( n-n1 < 20 ) print '(f10.2)',(dddata(l),l=n1+1,n)
            end if
            if ( n >= 4 .and. n-n1 >= 4 ) then
!               NumRec routine
!               NumRec routines
                write(0,*) 'printcorr: error: KS-test not yet implemented'
                call exit(-1)
                !!!call ttestx(dddata,n1,dddata(n1+1),n-n1,zd,pt,adat1,avar1,adat2,avar2)
                !!!call kstwo(dddata,n1,dddata(n1+1),n-n1,d,prob)
            else
                print *,'not enough data below/above cutoff:',n1,n-n1
                goto 800
            end if
            if ( lprint ) then
                1001 format(a8,f8.2,2f8.3,i5,2f8.2,i5,2f8.2)
                print 1001,string,kscut,100*(1-pt),100*(1-prob),n1,adat1,sqrt(avar1),n-n1 &
                    ,adat2,sqrt(avar2)
                if ( dump ) write(10,1001) ' ',kscut,100*(1-pt),100*(1-prob),n1,adat1 &
                    ,sqrt(avar1),n-n1,adat2,sqrt(avar2)
            end if
            result = pt
            do i=-2,2
                dresult(i) = 0
            end do
        end do
    else if ( lconting ) then

!       contingency tables  - only 3x3 at the moment

        if ( n <= 0 ) then
            print *,'no data for contingency table'
        else
            do j=0,3
                do i=0,3
                    ncont(i,j) = 0
                end do
            end do
            do l=1,n
                if ( dindx(l) < minindx ) then
                    i = 1
                else if ( dindx(l) < maxindx ) then
                    i = 2
                else
                    i = 3
                end if
                if ( ddata(l) < mindata ) then
                    j = 1
                else if ( ddata(l) < maxdata ) then
                    j = 2
                else
                    j = 3
                end if
                ncont(i,j) = ncont(i,j) + 1
            end do
            do j=1,3
                do i=1,3
                    ncont(0,j) = ncont(0,j) + ncont(i,j)
                    ncont(i,0) = ncont(i,0) + ncont(i,j)
                end do
            end do
            ncont(0,0) = n
!           compute significance (NumRec p.624)
            chi2 = 0
            do i=1,3
                do j=1,3
                    xn = ncont(i,0)*ncont(0,j)/real(n)
                    if ( xn > 0 ) then
                        chi2 = chi2 + (ncont(i,j)-xn)**2/xn
                    end if
                end do
            end do
            if ( month == 0 ) then
                df = 4/(max(lsum,lsum2) + decor + max(1,-ndiff)-1)
            else
                df = 4/(1 + (max(lsum,lsum2)-1)/nperyear &
                    + decor/nperyear)/real(max(1,1-ndiff))
            end if
            prob = gammq(0.5*df,0.5*chi2)
!           print results
            if ( lprint ) then
                print '(a,2f12.2)','cutoff data            ',mindata,maxdata
           1010 format(a8,f10.2,4(i5,1x,'(',i3,'%)'),f7.2,i4)
           1011 format(a8,a10,4(i5,1x,'(',i3,'%)'),f7.2,i4)
                if ( ncont(3,0) > 0 ) &
                print 1011,string,'          ',(ncont(3,j) &
                    ,nint(100*ncont(3,j)/real(n)),j=1,3),ncont(3,0 &
                    ),nint(100*ncont(3,0)/real(n))
                if ( ncont(2,0) > 0 ) &
                print 1010,string,maxindx,(ncont(2,j),nint(100 &
                    *ncont(2,j)/real(n)),j=1,3),ncont(2,0) &
                    ,nint(100*ncont(2,0)/real(n))
                if ( ncont(1,0) > 0 ) &
                print 1010,string,minindx,(ncont(1,j),nint(100 &
                    *ncont(1,j)/real(n)),j=1,3),ncont(1,0) &
                    ,nint(100*ncont(1,0)/real(n))
                print 1011,string,'sums      ', &
                    (ncont(0,j),nint(100*ncont(0,j)/real(n)),j=1,3 &
                    ),ncont(0,0),nint(100*ncont(0,0)/real(n)),100 &
                    *(1-prob),lag
                print *
           1015 format(2a,f7.2,i4)
                if ( dump ) write(10,1015)'# ',string,100*(1-prob),lag
                if ( plot ) then
               1012 format(i3,i5,g10.2,16i6,4g12.2)
                    write(11,1012) month,lag,prob,ncont,mindata &
                        ,maxdata,minindx,maxindx
                    if ( month == 0 .and. m1 /= m2) write(11,'(a)')
                    if ( lag1 /= lag2 .and. lag == lag2 ) write(11,'(a)')
                end if
            end if
        end if
    else if ( nfittime > 0 ) then
        if ( n < max(minnum,5) ) then
            print *,'not enough data for time fit'
            goto 800
        end if

!       fit to data = a*d(data)/dt + b*indx

        do j=1,n
            x(j) = j
        end do
        do j=1,n
            sig(j) = 1
        end do
        call svdfit(x,ddata,sig,n,a,2,u,v,w,ndata,2,chi2,findx)
        call svdvar(v,2,2,w,da,2)
        do j=1,n
            x(j) = a(1)*dddata(j) + a(2)*dindx(j)
        end do
        if ( month == 0 ) then
            df = (n-n0)/(max(lsum,lsum2) + decor) &
            /real(max(1,1-ndiff)) - 3
        else
            df = (n-n0)/(1 + (max(lsum,lsum2)-1)/nperyear + &
            decor/nperyear)/real(max(1,1-ndiff)) - 3
        end if
        call pearsncross(ddata,x,n,r,prob,z,adata,sxx,aindx,syy,sxy,df,ncrossvalidate)
        if ( lprint ) then
            1003 format(a8,f7.3,f7.2,i5,f10.4,f8.4,f10.4,f8.4,f6.2,i4)
            if ( da(1,1) == 0 ) da(1,1) = 1e-33
            if ( da(2,2) == 0 ) da(2,2) = 1e-33
            print 1003,string,r,100*(1-prob),n,a(1),sqrt(da(1,1) &
                /chi2),a(2),sqrt(da(2,2)/chi2),da(1,2)/sqrt(da(1,1)*da(2,2)),lag
            if ( plot ) then
                write(11,1020) month,lag,r,prob,n-n0,adata &
                    ,sqrt(sxx/(n-1)),aindx,sqrt(syy/(n-1)), &
                    a(1),da(1,1),a(2),da(2,2)
                if ( month == 0 .and. m1 /= m2) write(11,'(a)')
                if ( lag1 /= lag2 .and. lag == lag2 ) write(11,'(a)')
            end if
        end if
    else                    ! not KS or conting or fittime

!       correlations

        if ( n >= max(minnum,3) ) then
            if ( lnooverlap ) then
                df = (n-n0) - 2
            else if ( j1 /= j2 ) then
                df = (n-n0)/(max(lsum,lsum2) + decor)/ &
                real(max(1,1-ndiff,1-ndiff2)) - 2
            else
                df = (n-n0)/(1+(max(lsum,lsum2)-1)/nperyear &
                + decor/nperyear)/real(max(1,1-ndiff,1-ndiff2)) - 2
                if ( lwrite ) then
                    print *,'n,n0             = ',n,n0
                    print *,'lsum,lsum2       = ',lsum,lsum2
                    print *,'decor            = ',decor
                    print *,'1-ndiff,1-ndiff2 = ',1-ndiff,1-ndiff2
                    print *,'              df = ',df
                end if
            end if
            call getautocor1(dindx,ddata,n,a1,a2,lwrite)
!           if autocorrelation is significantly different from zero and does not increase df
            if ( a1 < 1e33 .and. a1 > exp(-1.) .and. &
            a1 > 1/sqrt(real(n)) ) then
                df = (n-n0-2)*(-log(a1))
                if ( lwrite ) then
                    print *,'from lag-1    df = ',df
                end if
            else
                df = (n-n0-2)
            end if
            allocate(w1(ndata),w2(ndata))
            if ( lrank ) then
!               NumRec
                if ( j1 /= j2 ) then
                    sum = max(lsum,lsum2) + decor
                else
                    sum = 1 + (max(lsum,lsum2)-1)/nperyear + decor/nperyear
                end if
                call spearx(ddata,dindx,n,w1,w2,d,zd,probd,r,prob &
                    ,sum,adata,sxx,aindx,syy)
                result = r
                if ( lboot ) then
                    ndecor = int(n/df+0.001)
                    call bootstrap(ddata,dindx,lfirst,u(1,1),u(1,2) &
                        ,n,dresult,2,ndecor,w1,w2,ncrossvalidate)
                end if
            else
                if ( lrmse .or. lmae ) then
                    allocate(wdata(n))
                    allocate(windx(n))
                    wdata = 1
                    windx = 1
                    ndecor = int(n/df+0.001)
                    if ( lrmse ) then
                        call getrms(ddata,wdata,dindx,windx,n,bias, &
                            adata,sxx,aindx,syy,rmse)
                        result = rmse
                        if ( lboot ) then
                            call bootstrap(ddata,dindx,lfirst, &
                                u(1,1),u(1,2),n,dresult,3,ndecor, &
                                w1,w2,ncrossvalidate)
                        end if
                    else if ( lmae ) then
                        call getmae(ddata,wdata,dindx,windx,n,bias, &
                            adata,sxx,aindx,syy,mae)
                        result = mae
                        if ( lboot ) then
                            call bootstrap(ddata,dindx,lfirst, &
                                u(1,1),u(1,2),n,dresult,4,ndecor, &
                                w1,w2,ncrossvalidate)
                        end if
                    else
                        write(0,*) 'printcorr: error hbgwkfwghv'
                        call exit(-1)
                    end if
                    sxx = sxx*n ! agree with stupid definition is pearsnx...
                    syy = syy*n
                    deallocate(wdata)
                    deallocate(windx)
                else
!                   NumRecp
                    call pearsncross(ddata,dindx,n,r,prob,z,adata, &
                        sxx,aindx,syy,sxy,df,ncrossvalidate)
                    result = r
                    call fitcross(dindx,ddata,n,sig,0,a(2),a(1), &
                        da(2,2),da(1,1),chi2,q,ncrossvalidate, &
                        aa,bb, .true. )
!                   write out bb
                    if ( lprint .and. ncrossvalidate > 0 .and. &
                            bbfile /= ' ' ) then
                        call rsunit(nunit)
                        ilen = index(bbfile,'.dat') - 1
                        if ( ilen <= 0 ) ilen = len_trim(bbfile)
                        open(nunit,file=bbfile(:ilen)//'_b2.dat')
!                       overwrite if there happens to one already
                        write(nunit,'(a)') '# cross-validated '// &
                            'regreesion coefficients B (regr)'
                        allocate(bb1(nperyear,yrmo(1,1):yrmo(1,n)))
                        bb1 = 3e33
                        do i=1,n
                            bb1(yrmo(2,i),yrmo(1,i)) = bb(i)
                        end do
                        call printdatfile(nunit,bb1,nperyear, &
                            nperyear,yrmo(1,1),yrmo(1,n))
                        close(nunit)
                        open(nunit,file=bbfile(:ilen)//'_a.dat')
!                       overwrite if there happens to one already
                        write(nunit,'(a)') '# cross-validated '// &
                            'regreesion coefficients A (const)'
                        bb1 = 3e33
                        do i=1,n
                            bb1(yrmo(2,i),yrmo(1,i)) = aa(i)
                        end do
                        call printdatfile(nunit,bb1,nperyear, &
                            nperyear,yrmo(1,1),yrmo(1,n))
                        close(nunit)
                        deallocate(bb1)
                    end if
!!!                 print *,'df was ',df
                    if ( lboot ) then
                        ndecor = int(n/df+0.001)
                        if ( lwrite ) print *,'calling bootstrap with ndecor = ',ndecor
                        call bootstrap(ddata,dindx,lfirst, &
                            u(1,1),u(1,2),n,dresult,1,ndecor, &
                            w1,w2,ncrossvalidate)
                    end if
                end if
            end if
            deallocate(w1,w2)
        else
            if ( lprint ) then
                print *,'printcorr: not enough data for correlation ',n,max(minnum,2)
                if ( plot ) then
                    write(11,'(a)')
                end if
            end if
            goto 800
        end if
        if ( lprint ) then
       1000 format(a8,f7.3,f7.2,i6,g10.2e1,g8.2e1,g10.2e1,g8.2e1,i4,5f7.3)
       2000 format(5a,i4,a,f7.3,a,f7.4,a,i6,a,f6.2,a,f6.2)
            if ( lboot ) then
                if ( lweb ) then
                    print 2000,'<tr><td>',string, &
                        '</td><td align="right">',corrmonths, &
                        '</td><td align="right">',lag, &
                        '</td><td align="right">',result, &
                        '</td><td align="right">',prob, &
                        '</td><td align="right">',n-n0, &
                        '</td><td align="right">',dresult(-2), &
                        '...',dresult(2), &
                        '</td></tr>'
                else
                    print 1000,string,result,100*(1-prob),n-n0,adata &
                        ,sqrt(sxx/(n-1)),aindx,sqrt(syy/(n-1)),lag &
                        ,dresult
                end if
            else
                if ( lweb ) then
                    print 2000,'<tr><td>',string, &
                        '</td><td align="right">',corrmonths, &
                        '</td><td align="right">',lag, &
                        '</td><td align="right">',result, &
                        '</td><td align="right">',prob, &
                        '</td><td align="right">',n-n0, &
                        '</td></tr>'
                else
                    print 1000,string,result,100*(1-prob),n-n0,adata &
                        ,sqrt(sxx/(n-1)),aindx,sqrt(syy/(n-1)),lag
                end if
            end if
       1100 format(a,a8,f7.3,f7.2,i6,g10.2e1,g8.2e1,g10.2e1,g8.2e1,i4,2g10.2e1)
            if ( dump ) then
                if ( df > 0 ) then
                    xneff = n-2*(n-2)/df
                else
                    xneff = -3e33
                end if
                if ( xneff > 0 ) then
                    chieff = sqrt(chi2/xneff)
                else
                    chieff = -999.9
                end if
                write(10,1100)'# ',string,r,100*(1-prob),n &
                    -n0,adata,sqrt(sxx/(n-1)),aindx,sqrt(syy/(n-1)) &
                    ,lag,chieff,rmse
            end if
            if ( plot ) then
           1002 format(i3,i5,f8.3,g10.2,i6,g12.4e2,g9.3e1,g12.4e1,g9.3e1,5f7.3,4g12.4e1)
           1020 format(i3,i5,f8.3,g10.2,i6,g12.4e2,g9.3e1,g12.4e1,g9.3e1,4g12.4e1)
                if ( lboot ) then
                    write(11,1002) month,lag,result,prob,n-n0,adata &
                        ,sqrt(sxx/(n-1)),aindx,sqrt(syy/(n-1)) &
                        ,dresult,a(1),da(1,1),a(2),da(2,2)
                else
                    write(11,1020) month,lag,result,prob,n-n0,adata &
                        ,sqrt(sxx/(n-1)),aindx,sqrt(syy/(n-1)), &
                        a(1),da(1,1),a(2),da(2,2)
                end if
                if ( month == 0 .and. m1 /= m2) write(11,'(a)')
                if ( lag1 /= lag2 .and. lag == lag2 ) &
                write(11,'(a)')
            end if
        end if
    end if
800 continue
end subroutine printcorr

subroutine findx(xi,f,n)

!       used by the multiple-parameter fitting routine (lfittime)

    implicit none
    include 'param.inc'
    integer :: n
    real :: xi,f(n)
    real :: dddata(ndata),dindx(ndata)
    common /c_findx/ dddata,dindx
    integer :: i,j

    if ( n /= 2 ) goto 901
    i = nint(xi)
    if ( abs(xi-i) > 0.01 ) goto 902
    f(1) = dddata(i)
    f(2) = dindx(i)
    return
901 print *,'findx: should be called with n=2, not ',n
    call exit(-1)
902 print *,'findx: wrong input! ',xi
    call exit(-1)
end subroutine findx
