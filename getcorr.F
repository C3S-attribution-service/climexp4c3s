        subroutine getcorr(dindx,ddata,dddata,n,n0,month,nperyear,lag
     +        ,string,lboot)
*
*       actually perform the correlations based on the linear ararys
*
        implicit none
#include "param.inc"
#include "getopts.inc"
        integer n,n0,month,nperyear,lag
        real dindx(n),ddata(n),dddata(n)
        character*(*) string
        logical lboot
        integer i,j
        integer iindx(ndata),ks,ks1,ks2,l,n1
        integer ncont(0:3,0:3)
        real zd,pt,adat1,avar1,adat2,avar2,prob,d
        real chi2,xn,x(ndata),sig(ndata),u(ndata,2),v(2,2),w(2),da(2,2)
        real df,sum,a(2),w1(ndata),w2(ndata),q,r,rr,z,adata,sxx,aindx
     +        ,syy,sxy,probd,aa(ndata),bb(ndata)
        real gammq
        external gammq,findx
*
        if ( n.gt.ndata ) then
            write(0,*) 'getcorr: increase ndata to ',n
            write(*,*) 'getcorr: increase ndata to ',n
            call abort
        endif
*
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
*       NumRec routine
                    call ttestx(dddata,n1,dddata(n1+1),n-n1,zd,pt,adat1
     +                    ,avar1,adat2,avar2)
*       NumRec routine
                    call kstwo(dddata,n1,dddata(n1+1),n-n1,d,prob)
                else
                    print *,'not enough data below/'//
     +                    'above cutoff:',n1,n-n1
                    goto 800
                endif
 1001           format(a,f8.2,2f8.3,i5,2f8.2,i5,2f8.2)
                print 1001,string,kscut,100*(1-pt),
     +                100*(1-prob),n1,adat1,sqrt(avar1),n-n1
     +                ,adat2,sqrt(avar2)
                if ( dump ) write(10,1001) ' ',kscut,
     +                100*(1-pt),100*(1-prob),n1,adat1
     +                ,sqrt(avar1),n-n1,adat2,sqrt(avar2)
            enddo
        elseif ( lconting ) then
*
*       contingency tables  - only 3x3 at the moment
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
*       compute significance (NumRec p.624)
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
*       print results
                print '(a,2f12.2)','cutoff data            ',mindata
     +                ,maxdata
 1010           format(a,f10.2,4(i5,1x,'(',i3,'%)'),f7.2,i4)
 1011           format(a,a10,4(i5,1x,'(',i3,'%)'),f7.2,i4)
                if ( ncont(3,0).gt.0 )
     +                print 1011,string,'          ',(ncont(3,j)
     +                ,nint(100*ncont(3,j)/real(n)),j=1,3),ncont(3,0)
     +                ,nint(100*ncont(3,0)/real(n))
                if ( ncont(2,0).gt.0 )
     +                print 1010,string,maxindx,(ncont(2,j),nint(100
     +                *ncont(2,j)/real(n)),j=1,3),ncont(2,0),nint(100
     +                *ncont(2,0)/real(n))
                if ( ncont(1,0).gt.0 )
     +                print 1010,string,minindx,(ncont(1,j),nint(100
     +                *ncont(1,j)/real(n)),j=1,3),ncont(1,0),nint(100
     +                *ncont(1,0)/real(n))
                print 1011,string,'sums      ',
     +                (ncont(0,j),nint(100*ncont(0,j)/real(n)),j=1,3)
     +                ,ncont(0,0),nint(100*ncont(0,0)/real(n)),100*(1
     +                -prob),lag
                print *
 1015           format(a,f7.2,i4)
                if ( dump ) write(10,1015)'# '//string,100*(1-prob)
     +                ,lag
                if ( plot ) then
 1012               format(i3,i5,g10.2,16i6,4g12.2)
                    write(11,1012) month,lag,prob,ncont,mindata,maxdata
     +                    ,minindx,maxindx
                    if ( month.eq.0 .and. m1.ne.m2)
     +                    write(11,'(a)')
                    if ( lag1.ne.lag2 .and. lag.eq.lag2 )write(11,'(a)')
                endif
            endif
        elseif ( nfittime.gt.0 ) then
            if ( n.le.4 ) then
                print *,'not enough data for time fit'
                goto 800
            endif
*
*       fit to data = a*d(data)/dt + b*indx
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
 1003       format(a,f6.2,f7.2,i5,f10.4,f8.4,f10.4,f8.4,f6.2,i4)
            print 1003,string,r,100*(1-prob),n,a(1),sqrt(da(1,1)
     +            /chi2),a(2),sqrt(da(2,2)/chi2),da(1,2)/sqrt(da(1,1)
     +            *da(2,2)),lag
        else                    ! not KS or conting or fittime
*
*       correlations
*
            if ( n.gt.2 ) then
                if ( lrank ) then
*       NumRec
                    if ( month.eq.0 ) then
                        sum = max(lsum,lsum2) + decor
                    else
                        sum = 1 + (max(lsum,lsum2)-1)
     +                        /nperyear + decor/nperyear
                    endif
                    call spearx(ddata,dindx,n,w1,w2,d,zd,probd,r,prob
     +                    ,sum,adata,sxx,aindx,syy)
                    if ( lboot ) then
                        call bootstrap(ddata,dindx,u(1,1),u(1,2),n,rr,2
     +                        ,w1,w2,ncrossvalidate)
                    endif
                else
*       NumRecp
                    if ( month.eq.0 ) then
                        df = (n-n0)/(max(lsum,lsum2) + decor)/
     +                        real(max(1,1-ndiff)) - 2
                    else
                        df = (n-n0)/(1+(max(lsum,lsum2)-1)/nperyear
     +                        + decor/nperyear)/real(max(1,1-ndiff)) - 2
                    endif
                    call pearsncross(ddata,dindx,n,r,prob,z,adata,sxx
     +                   ,aindx,syy,sxy,df,ncrossvalidate)
                    if ( lboot ) then
                        call bootstrap(ddata,dindx,u(1,1),u(1,2),n,rr,1
     +                        ,w1,w2,ncrossvalidate)
                    endif
                    if ( sxx.ne.0 .and. syy.ne.0 ) then
                        call fitcross(dindx,ddata,n,sig,0,a(2),a(1),
     +                       da(2,2),da(1,1),chi2,q,ncrossvalidate,aa,bb
     +                       ,.true.)
                    endif
                endif
            else
                print *,'not enough data for correlation'
                goto 800
            endif
 1000       format(a,f6.2,f7.2,i6,g10.2e1,g8.2e1,g10.2e1,g8.2e1,i4,5f6.2
     +            )
            if ( lboot ) then
                print 1000,string,r,100*(1-prob),n-n0,adata
     +                ,sqrt(sxx/n),aindx,sqrt(syy/n),lag,rr
            else
                print 1000,string,r,100*(1-prob),n-n0,adata
     +                ,sqrt(sxx/n),aindx,sqrt(syy/n),lag
            endif
            if ( dump ) write(10,1000)'# '//string,r,100*(1-prob),n
     +            -n0,adata,sqrt(sxx/n),aindx,sqrt(syy/n),lag
            if ( plot ) then
 1002           format(i3,i5,f6.2,g10.2,i6,g12.4e1,g9.3e1,g12.4e1,g9.3e1
     +                ,5f6.2)
 1020           format(i3,i5,f6.2,g10.2,i6,g12.4e1,g9.3e1,g12.4e1,g9.3e1
     +                ,4g12.4e1)
                if ( lboot ) then
                    write(11,1002) month,lag,r,prob,n-n0,adata
     +                    ,sqrt(sxx/n),aindx,sqrt(syy/n),rr
                else
                    write(11,1020) month,lag,r,prob,n-n0,adata
     +                    ,sqrt(sxx/n),aindx,sqrt(syy/n),a(1),da(1,1)
     +                    ,a(2),da(2,2)
                endif
                if ( month.eq.0 .and. m1.ne.m2) write(11,'(a)')
                if ( lag1.ne.lag2 .and. lag.eq.lag2 )
     +                write(11,'(a)')
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

