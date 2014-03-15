        subroutine printcorr(dindx,ddata,lfirst,dddata,yrmo,n,n0,j1,j2
     +       ,month,nperyear,lag,string,lboot,lprint,result,dresult,prob
     +       )
*
*       actually perform the correlations based on the linear ararys
*       and print the result (if lprint)
*
        implicit none
#include "param.inc"
#include "getopts.inc"
        integer n,n0,month,nperyear,lag,j1,j2,yrmo(2,n)
        real dindx(n),ddata(n),dddata(n),result,dresult(-2:2),prob
        logical lfirst(n)
        character*(*) string
        logical lboot,lprint
        integer i,j,nunit,ilen
        integer iindx(ndata),ks,ks1,ks2,l,n1,ndecor
        integer ncont(0:3,0:3)
        real zd,pt,adat1,avar1,adat2,avar2,d,rmse,mae,a1,a2
        real chi2,xn,x(ndata),sig(ndata),u(ndata,2),v(2,2),w(2),da(2,2)
     +       ,aa(ndata),bb(ndata)
        real df,sum,a(2),w1(ndata),w2(ndata),q,r,z,adata,sxx
     +       ,aindx,syy,sxy,probd,bias,xneff,chieff
        real,allocatable :: bb1(:,:),wdata(:),windx(:)
        real gammq
        external gammq,findx
*
        if ( n.gt.ndata ) then
            write(0,*) 'getcorr: increase ndata to ',n
            write(*,*) 'getcorr: increase ndata to ',n
            call abort
        endif
*       
        result = 3e33
        if ( lks ) then
*       
*       Kolmogorov-Smirnov test
*       
            if ( kscut.lt.-1e33 ) then
                ks1 = 5
                ks2 = n-4
            else
                ks1 = 0
                ks2 = 0
            endif
            call ffsort(dindx,iindx,n)
            do ks=ks1,ks2
*       NumRec vernaggelt (sorteert) dddata...
                do l=1,n
                    dddata(l) = ddata(iindx(l))
                enddo
                if ( ks1.ne.0 ) then
                    kscut = (dindx(iindx(ks-1))
     +                    +dindx(iindx(ks)))/2
                endif
                do l=1,n
                    if ( dindx(iindx(l)).gt.kscut ) goto 730
                enddo
  730           n1 = l-1
                if ( lwrite ) then
                    print *,'Found ',n1,' points below ',kscut
                    if ( n1.lt.20 ) print '(f10.2)',(dddata(l),l=1,n1)
                    print *,'Found ',n-n1,' points above ',kscut
                    if ( n-n1.lt.20 ) print '(f10.2)',(dddata(l),l=n1+1
     +                    ,n)
                endif
                if ( n.ge.4 .and. n-n1.ge.4 ) then
*                   NumRec routine
                    call ttestx(dddata,n1,dddata(n1+1),n-n1,zd,pt,adat1
     +                    ,avar1,adat2,avar2)
*                   NumRec routine
                    call kstwo(dddata,n1,dddata(n1+1),n-n1,d,prob)
                else
                    print *,'not enough data below/'//
     +                    'above cutoff:',n1,n-n1
                    goto 800
                endif
                if ( lprint ) then
 1001               format(a8,f8.2,2f8.3,i5,2f8.2,i5,2f8.2)
                    print 1001,string,kscut,100*(1-pt),
     +                    100*(1-prob),n1,adat1,sqrt(avar1),n-n1
     +                    ,adat2,sqrt(avar2)
                    if ( dump ) write(10,1001) ' ',kscut,
     +                    100*(1-pt),100*(1-prob),n1,adat1
     +                    ,sqrt(avar1),n-n1,adat2,sqrt(avar2)
                endif
                result = pt
                do i=-2,2
                    dresult(i) = 0
                enddo
            enddo
        elseif ( lconting ) then
*
*           contingency tables  - only 3x3 at the moment
*
            if ( n.le.0 ) then
                print *,'no data for contingency table'
            else
                do j=0,3
                    do i=0,3
                        ncont(i,j) = 0
                    enddo
                enddo
                do l=1,n
                    if ( dindx(l).lt.minindx ) then
                        i = 1
                    elseif ( dindx(l).lt.maxindx ) then
                        i = 2
                    else
                        i = 3
                    endif
                    if ( ddata(l).lt.mindata ) then
                        j = 1
                    elseif ( ddata(l).lt.maxdata ) then
                        j = 2
                    else
                        j = 3
                    endif
                    ncont(i,j) = ncont(i,j) + 1
                enddo
                do j=1,3
                    do i=1,3
                        ncont(0,j) = ncont(0,j) + ncont(i,j)
                        ncont(i,0) = ncont(i,0) + ncont(i,j)
                    enddo
                enddo
                ncont(0,0) = n
*               compute significance (NumRec p.624)
                chi2 = 0
                do i=1,3
                    do j=1,3
                        xn = ncont(i,0)*ncont(0,j)/real(n)
                        if ( xn.gt.0 ) then
                            chi2 = chi2 + (ncont(i,j)-xn)**2/xn
                        endif
                    enddo
                enddo
                if ( month.eq.0 ) then
                    df = 4/(max(lsum,lsum2) + decor + max(1,-ndiff)-1)
                else
                    df = 4/(1 + (max(lsum,lsum2)-1)/nperyear
     +                    + decor/nperyear)/real(max(1,1-ndiff))
                endif
                prob = gammq(0.5*df,0.5*chi2)
*               print results
                if ( lprint ) then
                    print '(a,2f12.2)','cutoff data            ',mindata
     +                    ,maxdata
 1010               format(a8,f10.2,4(i5,1x,'(',i3,'%)'),f7.2,i4)
 1011               format(a8,a10,4(i5,1x,'(',i3,'%)'),f7.2,i4)
                    if ( ncont(3,0).gt.0 )
     +                    print 1011,string,'          ',(ncont(3,j)
     +                    ,nint(100*ncont(3,j)/real(n)),j=1,3),ncont(3,0
     +                    ),nint(100*ncont(3,0)/real(n))
                    if ( ncont(2,0).gt.0 )
     +                    print 1010,string,maxindx,(ncont(2,j),nint(100
     +                    *ncont(2,j)/real(n)),j=1,3),ncont(2,0)
     +                    ,nint(100*ncont(2,0)/real(n))
                    if ( ncont(1,0).gt.0 )
     +                    print 1010,string,minindx,(ncont(1,j),nint(100
     +                    *ncont(1,j)/real(n)),j=1,3),ncont(1,0)
     +                    ,nint(100*ncont(1,0)/real(n))
                    print 1011,string,'sums      ',
     +                    (ncont(0,j),nint(100*ncont(0,j)/real(n)),j=1,3
     +                    ),ncont(0,0),nint(100*ncont(0,0)/real(n)),100
     +                    *(1-prob),lag
                    print *
 1015               format(2a,f7.2,i4)
                    if ( dump ) write(10,1015)'# ',string,100*(1-prob)
     +                    ,lag
                    if ( plot ) then
 1012                   format(i3,i5,g10.2,16i6,4g12.2)
                        write(11,1012) month,lag,prob,ncont,mindata
     +                        ,maxdata,minindx,maxindx
                        if ( month.eq.0 .and. m1.ne.m2)
     +                        write(11,'(a)')
                        if ( lag1.ne.lag2 .and. lag.eq.lag2 )
     +                        write(11,'(a)')
                    endif
                endif
            endif
        elseif ( nfittime.gt.0 ) then
            if ( n.lt.max(minnum,5) ) then
                print *,'not enough data for time fit'
                goto 800
            endif
*
*           fit to data = a*d(data)/dt + b*indx
*
            do j=1,n
                x(j) = j
            enddo
            do j=1,n
                sig(j) = 1
            enddo
            call svdfit(x,ddata,sig,n,a,2,u,v,w,ndata,2,chi2,findx)
            call svdvar(v,2,2,w,da,2)
            do j=1,n
                x(j) = a(1)*dddata(j) + a(2)*dindx(j)
            enddo
            if ( month.eq.0 ) then
                df = (n-n0)/(max(lsum,lsum2) + decor)
     +                /real(max(1,1-ndiff)) - 3
            else
                df = (n-n0)/(1 + (max(lsum,lsum2)-1)/nperyear +
     +                decor/nperyear)/real(max(1,1-ndiff)) - 3
            endif
            call pearsncross(ddata,x,n,r,prob,z,adata,sxx,aindx,syy,sxy
     +           ,df,ncrossvalidate)
            if ( lprint ) then
 1003           format(a8,f7.3,f7.2,i5,f10.4,f8.4,f10.4,f8.4,f6.2,i4)
                if ( da(1,1).eq.0 ) da(1,1) = 1e-33
                if ( da(2,2).eq.0 ) da(2,2) = 1e-33
                print 1003,string,r,100*(1-prob),n,a(1),sqrt(da(1,1)
     +                /chi2),a(2),sqrt(da(2,2)/chi2),da(1,2)/sqrt(da(1,1
     +                )*da(2,2)),lag
                if ( plot ) then
                    write(11,1020) month,lag,r,prob,n-n0,adata
     +                    ,sqrt(sxx/(n-1)),aindx,sqrt(syy/(n-1)),
     +                    a(1),da(1,1),a(2),da(2,2)
                    if ( month.eq.0 .and. m1.ne.m2) write(11,'(a)')
                    if ( lag1.ne.lag2 .and. lag.eq.lag2 )
     +                    write(11,'(a)')
                endif
            endif
        else                    ! not KS or conting or fittime
*
*       correlations
*
            if ( n.ge.max(minnum,3) ) then
                if ( lnooverlap ) then
                    df = (n-n0) - 2
                else if ( j1.ne.j2 ) then
                    df = (n-n0)/(max(lsum,lsum2) + decor)/
     +                   real(max(1,1-ndiff,1-ndiff2)) - 2
                else
                    df = (n-n0)/(1+(max(lsum,lsum2)-1)/nperyear
     +                   + decor/nperyear)/real(max(1,1-ndiff,1
     +                   -ndiff2)) - 2
                    if ( lwrite ) then
                        print *,'n,n0             = ',n,n0
                        print *,'lsum,lsum2       = ',lsum,lsum2
                        print *,'decor            = ',decor
                        print *,'1-ndiff,1-ndiff2 = ',1-ndiff,1
     +                       -ndiff2
                        print *,'              df = ',df
                    end if
                endif
                call getautocor1(dindx,ddata,n,a1,a2,lwrite)
!               if autocorrelation is significantly different from zero and does not increase df
                if ( a1.lt.1e33 .and. a1.gt.exp(-1.) .and. 
     +               a1.gt.1/sqrt(real(n)) ) then
                    df = (n-n0-2)*(-log(a1))
                    if ( lwrite ) then
                        print *,'from lag-1    df = ',df
                    end if
                else
                    df = (n-n0-2)
                end if
                if ( lrank ) then
*                   NumRec
                    if ( j1.ne.j2 ) then
                        sum = max(lsum,lsum2) + decor
                    else
                        sum = 1 + (max(lsum,lsum2)-1)
     +                        /nperyear + decor/nperyear
                    endif
                    call spearx(ddata,dindx,n,w1,w2,d,zd,probd,r,prob
     +                    ,sum,adata,sxx,aindx,syy)
                    result = r
                    if ( lboot ) then
                        ndecor = int(n/df+0.001)
                        call bootstrap(ddata,dindx,lfirst,u(1,1),u(1,2)
     +                       ,n,dresult,2,ndecor,w1,w2,ncrossvalidate)
                    endif
                else
                    if ( lrmse .or. lmae ) then
                        allocate(wdata(n))
                        allocate(windx(n))
                        wdata = 1
                        windx = 1
                        ndecor = int(n/df+0.001)
                        if ( lrmse ) then
                            call getrms(ddata,wdata,dindx,windx,n,bias,
     +                           adata,sxx,aindx,syy,rmse)
                            result = rmse
                            if ( lboot ) then
                                call bootstrap(ddata,dindx,lfirst,
     +                               u(1,1),u(1,2),n,dresult,3,ndecor,
     +                               w1,w2,ncrossvalidate)
                            endif
                        else if ( lmae ) then
                            call getmae(ddata,wdata,dindx,windx,n,bias,
     +                           adata,sxx,aindx,syy,mae)
                            result = mae
                            if ( lboot ) then
                                call bootstrap(ddata,dindx,lfirst,
     +                               u(1,1),u(1,2),n,dresult,4,ndecor,
     +                               w1,w2,ncrossvalidate)
                            endif
                        else
                            write(0,*) 'printcorr: error hbgwkfwghv'
                            call abort
                        end if
                        sxx = sxx*n ! agree with stupid definition is pearsnx...
                        syy = syy*n
                        deallocate(wdata)
                        deallocate(windx)
                    else
*                       NumRecp
                        call pearsncross(ddata,dindx,n,r,prob,z,adata,
     +                       sxx,aindx,syy,sxy,df,ncrossvalidate)
                        result = r
                        call fitcross(dindx,ddata,n,sig,0,a(2),a(1),
     +                       da(2,2),da(1,1),chi2,q,ncrossvalidate,
     +                       aa,bb,.true.)
!                       write out bb
                        if ( lprint .and. ncrossvalidate.gt.0 .and. 
     +                       bbfile.ne.' ' ) then
                            call rsunit(nunit)
                            ilen = index(bbfile,'.dat') - 1
                            if ( ilen.le.0 ) ilen = len_trim(bbfile)
                            open(nunit,file=bbfile(:ilen)//'_b2.dat')
!                           overwrite if there happens to one already
                            write(nunit,'(a)') '# cross-validated '//
     +                           'regreesion coefficients B (regr)'
                            allocate(bb1(nperyear,yrmo(1,1):yrmo(1,n)))
                            bb1 = 3e33
                            do i=1,n
                                bb1(yrmo(2,i),yrmo(1,i)) = bb(i)
                            end do
                            call printdatfile(nunit,bb1,nperyear,
     +                           nperyear,yrmo(1,1),yrmo(1,n))
                            close(nunit)
                            open(nunit,file=bbfile(:ilen)//'_a.dat')
!                           overwrite if there happens to one already
                            write(nunit,'(a)') '# cross-validated '//
     +                           'regreesion coefficients A (const)'
                            bb1 = 3e33
                            do i=1,n
                                bb1(yrmo(2,i),yrmo(1,i)) = aa(i)
                            end do
                            call printdatfile(nunit,bb1,nperyear,
     +                           nperyear,yrmo(1,1),yrmo(1,n))
                            close(nunit)
                            deallocate(bb1)
                        end if
!!!                     print *,'df was ',df
                        if ( lboot ) then
                            ndecor = int(n/df+0.001)
                            if ( lwrite ) print *,
     +                           'calling bootstrap with ndecor = ',
     +                           ndecor
                            call bootstrap(ddata,dindx,lfirst,
     +                           u(1,1),u(1,2),n,dresult,1,ndecor,
     +                           w1,w2,ncrossvalidate)
                        endif
                    end if
                endif
            else
                if ( lprint ) then
                    print *
     +                    ,'printcorr: not enough data for correlation '
     +                    ,n,max(minnum,2)
                    if ( plot ) then
                        write(11,'(a)')
                    endif
                endif
                goto 800
            endif
            if ( lprint ) then
 1000           format(a8,f7.3,f7.2,i6,g10.2e1,g8.2e1,g10.2e1,g8.2e1,i4
     +               ,5f7.3)
 2000           format(5a,i4,a,f7.3,a,f7.4,a,i6,a,f6.2,a,f6.2)
                if ( lboot ) then
                    if ( lweb ) then
                        print 2000,'<tr><td>',string,
     +                       '</td><td align="right">',corrmonths,
     +                       '</td><td align="right">',lag,
     +                       '</td><td align="right">',result,
     +                       '</td><td align="right">',prob,
     +                       '</td><td align="right">',n-n0,
     +                       '</td><td align="right">',dresult(-2),
     +                       '...',dresult(2),
     +                       '</td></tr>'
                    else
                        print 1000,string,result,100*(1-prob),n-n0,adata
     +                    ,sqrt(sxx/(n-1)),aindx,sqrt(syy/(n-1)),lag
     +                    ,dresult
                    endif
                else
                    if ( lweb ) then
                        print 2000,'<tr><td>',string,
     +                       '</td><td align="right">',corrmonths,
     +                       '</td><td align="right">',lag,
     +                       '</td><td align="right">',result,
     +                       '</td><td align="right">',prob,
     +                       '</td><td align="right">',n-n0,
     +                       '</td></tr>'
                    else
                        print 1000,string,result,100*(1-prob),n-n0,adata
     +                    ,sqrt(sxx/(n-1)),aindx,sqrt(syy/(n-1)),lag
                    endif
                endif
 1100           format(a,a8,f7.3,f7.2,i6,g10.2e1,g8.2e1,g10.2e1,g8.2e1
     +               ,i4,2g10.2e1)
                if ( dump ) then
                    xneff = n-2*(n-2)/df
                    if ( xneff.gt.0 ) then
                        chieff = sqrt(chi2/xneff)
                    else
                        chieff = -999.9
                    end if
                    write(10,1100)'# ',string,r,100*(1-prob),n
     +                -n0,adata,sqrt(sxx/(n-1)),aindx,sqrt(syy/(n-1))
     +                ,lag,chieff,rmse
                end if
                if ( plot ) then
 1002               format(i3,i5,f8.3,g10.2,i6,g12.4e2,g9.3e1,g12.4e1,g9
     +                    .3e1,5f7.3,4g12.4e1)
 1020               format(i3,i5,f8.3,g10.2,i6,g12.4e2,g9.3e1,g12.4e1,g9
     +                    .3e1,4g12.4e1)
                    if ( lboot ) then
                        write(11,1002) month,lag,result,prob,n-n0,adata
     +                        ,sqrt(sxx/(n-1)),aindx,sqrt(syy/(n-1))
     +                        ,dresult,a(1),da(1,1),a(2),da(2,2)
                    else
                        write(11,1020) month,lag,result,prob,n-n0,adata
     +                        ,sqrt(sxx/(n-1)),aindx,sqrt(syy/(n-1)),
     +                        a(1),da(1,1),a(2),da(2,2)
                    endif
                    if ( month.eq.0 .and. m1.ne.m2) write(11,'(a)')
                    if ( lag1.ne.lag2 .and. lag.eq.lag2 )
     +                    write(11,'(a)')
                endif
            endif
        endif
  800   continue
        end

        subroutine findx(xi,f,n)
*       
*       used by the multiple-parameter fitting routine (lfittime)
*       
        implicit none
#include "param.inc"
        integer n
        real xi,f(n)
        real dddata(ndata),dindx(ndata)
        common /c_findx/ dddata,dindx
        integer i,j
*       
        if ( n.ne.2 ) goto 901
        i = nint(xi)
        if ( abs(xi-i).gt.0.01 ) goto 902
        f(1) = dddata(i)
        f(2) = dindx(i)
        return
  901   print *,'findx: should be called with n=2, not ',n
        call abort
  902   print *,'findx: wrong input! ',xi
        call abort
        end

