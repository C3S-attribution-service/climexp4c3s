subroutine printpoint(i,f,ntype,x,s,yr)
    implicit none
    integer :: i,ntype,yr
    real :: f,x,s
    if ( ntype == 2 ) then
        if ( f > 0 ) then
            print '(i8,4g22.6,i12)',i,-log(-log(f)),x,s,1/(1-f),yr
        endif
    elseif ( ntype == 3 ) then
        print '(i8,4g22.6,i12)',i,-log(1-f),x,s,1/(1-f),yr
    elseif ( ntype == 4 ) then
        print '(i8,4g22.6,i12)',i,sqrt(-log(1-f)),x,s,1/(1-f),yr
    else
        write(0,*) 'histogram: error: unknown ntype ',ntype
        call exit(-1)
    endif
end subroutine printpoint
    

subroutine printreturnvalue(ntype,t,t25,t975,lweb)

!   print return times

    implicit none
    integer :: ntype
    real :: t(10),t25(10),t975(10)
    logical :: lweb
    integer :: i

    if ( ntype == 2 .or. ntype == 3 .or. ntype == 4 ) then ! extreme value  plot
        if ( lweb ) then
            do i=1,10,3
                print '(a,i5,a,f16.3,a,f16.3,a,f16.3,a)' &
                    ,'# <tr><td>return value ',10**(1+i/3) &
                    ,' yr</td><td>',t(i),'</td><td>',t25(i) &
                    ,' ... ',t975(i),'</td></tr>'
                call print3untransf(t(i),t25(i),t975(i),0)
            enddo
        else
            do i=1,4
                print '(a,i5,a,f16.3,a,2f16.3)' &
                    ,'# value for return period ',10**i,' year: ' &
                    ,t(i),' 95% CI ',t25(i),t975(i)
                call print3untransf(t(i),t25(i),t975(i),0)
            enddo
        endif
    endif
end subroutine printreturnvalue


subroutine printcovreturnvalue(ntype,t,t25,t975,yr1a,yr2a,lweb,plot,assume,lnone)

!   print return values for a set of fixed return times

    implicit none
    integer :: ntype,yr1a,yr2a
    real :: t(10,3),t25(10,3),t975(10,3)
    logical :: lweb,plot,lnone
    character assume*(*)
    integer :: i
    logical lprintreturnvalue
    character units*30

    lprintreturnvalue = .false.
    if ( ntype == 2 .or. ntype == 3 .or. ntype == 4 ) then ! extreme value  plot
        if ( assume == 'scale ' ) then
            write(units,'(a,i4)') ' % relative to ',yr1a
        else
            units = ' '
        end if
        if ( lweb ) then
            do i=1,10,3
                if ( lnone ) then
                    if ( lprintreturnvalue ) then
                        print '(a,i5,a,f16.3,a,f16.3,a,f16.3,a)' &
                            ,'# <tr><td colspan=2>return value ',10**(1+i/3) &
                            ,' yr</td><td>',t(i,1),'</td><td>',t25(i,1) &
                            ,' ... ',t975(i,1),'</td></tr>'
                        call print3untransf(t(i,1),t25(i,1),t975(i,1),-1)
                    end if
                else
                    if ( lprintreturnvalue ) then
                        print '(a,i5,a,i4,a,f16.3,a,f16.3,a,f16.3,a)' &
                            ,'# <tr><td>return value ',10**(1+i/3) &
                            ,' yr</td><td>',yr1a,'</td><td>' &
                            ,t(i,1),'</td><td>',t25(i,1) &
                            ,' ... ',t975(i,1),'</td></tr>'
                        print '(a,i4,a,f16.3,a,f16.3,a,f16.3,a)' &
                            ,'# <tr><td>&nbsp;</td><td>',yr2a &
                            ,'</td><td>',t(i,2),'</td><td>',t25(i,2) &
                            ,' ... ',t975(i,2),'</td></tr>'
                    end if
                    if ( i == 10 .or. assume == 'both' ) then
                        print '(3a,f16.3,a,f16.3,a,f16.3,a)' &
                            ,'# <tr><td>return values</td><td>diff',trim(units), &
                            '</td><td>',t(i,3),'</td><td>',t25(i,3) &
                            ,' ... ',t975(i,3),'</td></tr>'
                    end if
                    if ( lprintreturnvalue ) then
                        call print3untransf(t(i,1),t25(i,1),t975(i,1),yr1a)
                        call print3untransf(t(i,2),t25(i,2),t975(i,2),yr2a)
                    end if
                end if
            enddo
            if ( plot ) then ! output for stationlist
                do i=1,10
                    if ( lnone ) then
                        write(11,'(4g20.4,i6,a)') t(i,1),t25(i,1), &
                            t975(i,1),10**(1+i/3.),0,' return value'
                    else
                        write(11,'(4g20.4,i6,a)') t(i,1),t25(i,1), &
                            t975(i,1),10**(1+i/3.),yr1a,' return value'
                        write(11,'(4g20.4,i6,a)') t(i,2),t25(i,2), &
                            t975(i,2),10**(1+i/3.),yr2a,' return value'
                    end if
                end do
            end if
        else
            do i=1,4
                if ( lnone ) then
                    print '(a,i5,a,f16.3,a,2f16.3)' &
                        ,'# value for return period ',10**i &
                        ,': ' &
                        ,t(i,1),' 95% CI ',t25(i,1),t975(i,1)
                    call print3untransf(t(i,1),t25(i,1),t975(i,1),0)
                else
                    print '(a,i5,a,i4,a,f16.3,a,2f16.3)' &
                        ,'# value for return period ',10**i &
                        ,' year at yr=',yr1a,': ' &
                        ,t(i,1),' 95% CI ',t25(i,1),t975(i,1)
                    print '(a,i5,a,i4,a,f16.3,a,2f16.3)' &
                        ,'# value for return period ',10**i &
                        ,' year at yr=',yr2a,': ' &
                        ,t(i,2),' 95% CI ',t25(i,2),t975(i,2)
                    print '(a,i5,3a,f16.3,a,2f16.3)' &
                        ,'# value for return period ',10**i &
                        ,' year, difference',trim(units),': ' &
                        ,t(i,3),' 95% CI ',t25(i,3),t975(i,3)
                    call print3untransf(t(i,1),t25(i,1),t975(i,1),yr1a)
                    call print3untransf(t(i,2),t25(i,2),t975(i,2),yr2a)
                end if
            enddo
        endif
    endif
end subroutine printcovreturnvalue


subroutine printreturntime(year,xyear,tx,tx25,tx975,lweb)

!   print return time of year

    implicit none
    integer :: year
    real :: xyear,tx,tx25,tx975
    logical :: lweb
            
    if ( xyear < 1e33 ) then
        if ( lweb ) then
            if ( year == 9999 ) then
                print '(a,g16.5,a,g16.5,a,g16.5,a,g16.5,a)' &
                    ,'# <tr><td>return period ',xyear,'</td><td>' &
                    ,tx,'</td><td>',tx25,' ... ',tx975,'</td></tr>'
                print '(a,g16.5,a,g16.5,a,g16.5,a,g16.5,a)' &
                    ,'# <tr><td>probability ',xyear,'</td><td>' &
                    ,1/tx,'</td><td>',1/tx975,' ... ',1/tx25,'</td></tr>'
            else
                print '(a,g16.5,a,i5,a,g16.5,a,g16.5,a,g16.5,a)' &
                    ,'# <tr><td>return period ',xyear,'(',year &
                    ,')</td><td>' &
                    ,tx,'</td><td>',tx25,' ... ',tx975,'</td></tr>'
                print '(a,g16.5,a,i5,a,g16.5,a,g16.5,a,g16.5,a)' &
                    ,'# <tr><td>probability ',xyear,'(',year &
                    ,')</td><td>' &
                    ,1/tx,'</td><td>',1/tx975,' ... ',1/tx25,'</td></tr>'
            end if
        else
            print '(a,f16.5,a,i4,a,f16.5,a,2f18.5)','# return time ' &
                ,xyear,' (',year,') = ',tx,' 95% CI ',tx25,tx975
        endif
    endif
end subroutine printreturntime


subroutine printcovreturntime(year,xyear,idmax,tx,tx25,tx975,yr1a,yr2a,yr2b,lweb,plot,assume,lnone,i12)

!   print return time of year at cov1 and cov2 (and optionally cov3)

    implicit none
    integer :: year,yr1a,yr2a,yr2b,i12
    real :: xyear,tx(4),tx25(4),tx975(4)
    logical :: lweb,plot,lnone
    character idmax*(*),assume*(*)
    integer :: i,nj
    character atx(4)*16,atx25(4)*16,atx975(4)*16, &
        ainvtx(4)*16,ainvtx25(4)*16,ainvtx975(4)*16

    if ( xyear < 1e33 ) then
        if ( yr2b > 0 ) then
            nj = 4
        else
            nj = 3
        end if
        do i=1,nj
            call val_or_inf(atx(i),tx(i),lweb)
            call val_or_inf(atx25(i),tx25(i),lweb)
            call val_or_inf(atx975(i),tx975(i),lweb)
            call val_or_inf(ainvtx(i),1/tx(i),lweb)
            if ( tx25(i) == 0 ) tx25(i) = 1e-20
            call val_or_inf(ainvtx25(i),1/tx25(i),lweb)
            call val_or_inf(ainvtx975(i),1/tx975(i),lweb)
        end do
        if ( lweb ) then
            if ( .not. lnone ) then
                if ( i12 == 1 ) then
                    if ( idmax == ' ' ) then
                        print '(a,i5,a,g16.5,a,i5,7a)' &
                            ,'# <tr><td><!--atr2-->return period event ',year, &
                            ' (value ',xyear,')</td><td>',yr1a,'</td><td>',atx(1), &
                            '</td><td>',atx25(1),' ... ',atx975(1),'</td></tr>'
                    else
                        print '(a,i5,a,g16.5,3a,i5,7a)' &
                            ,'# <tr><td><!--atr2-->return period event ',year, &
                            ' (value ',xyear,' at ',trim(idmax),')</td><td>',yr1a,'</td><td>',atx(1), &
                            '</td><td>',atx25(1),' ... ',atx975(1),'</td></tr>'
                    end if
                    print '(a,i5,7a)' &
                        ,'# <tr><td>probability</td><td>',yr1a,'</td><td>',ainvtx(1), &
                        '</td><td>',ainvtx975(1),' ... ',ainvtx25(1),'</td></tr>'
                    if ( idmax == ' ' ) then
                        print '(a,i5,a,g16.5,a,i5,7a)' &
                            ,'# <tr><td><!--atr1-->return period event ',year, &
                            ' (value ',xyear,')</td><td>',yr2a,'</td><td>',atx(2), &
                            '</td><td>',atx25(2),' ... ',atx975(2),'</td></tr>'
                    else
                        print '(a,i5,a,g16.5,3a,i5,7a)' &
                            ,'# <tr><td><!--atr1-->return period ',year, &
                            ' (value ',xyear,' at ',trim(idmax),')</td><td>',yr2a,'</td><td>',atx(2), &
                            '</td><td>',atx25(2),' ... ',atx975(2),'</td></tr>'
                    end if
                    print '(a,i5,7a)' &
                        ,'# <tr><td>probability</td><td>' &
                        ,yr2a,'</td><td>',ainvtx(2), &
                        '</td><td>',ainvtx975(2),' ... ',ainvtx25(2),'</td></tr>'
                    print '(8a)' &
                            ,'# <tr><td><!--atra-->probability ratio</td><td>&nbsp;', &
                            '</td><td>',atx(3),'</td><td>',atx25(3), &
                            ' ... ',atx975(3),'</td></tr>'
                    if ( tx(3) < 1 ) then
                        print '(8a)' &
                            ,'# <tr><td>inverse probability ratio</td><td>&nbsp;', &
                            '</td><td>',ainvtx(3),'</td><td>',ainvtx975(3), &
                            ' ... ',ainvtx25(3),'</td></tr>'
                    end if
                else
                    if ( yr2b > 0 ) then
                        if ( idmax == ' ' ) then
                            print '(a,i5,a,g16.5,a,i5,7a)' &
                                ,'# <tr><td>return period event ',year, &
                                ' (value ',xyear,')</td><td>',yr2b,'</td><td>',atx(4), &
                                '</td><td>',atx25(4),' ... ',atx975(4),'</td></tr>'
                        else
                            print '(a,i5,a,g16.5,3a,i5,7a)' &
                                ,'# <tr><td>return period event ',year, &
                                ' (value ',xyear,' at ',trim(idmax),')</td><td>',yr2b,'</td><td>',atx(4), &
                                '</td><td>',atx25(4),' ... ',atx975(4),'</td></tr>'
                        end if
                        print '(a,i5,7a)' &
                            ,'# <tr><td>probability</td><td>',yr2b,'</td><td>',ainvtx(4), &
                            '</td><td>',ainvtx975(4),' ... ',ainvtx25(4),'</td></tr>'
                    end if
                end if
            else
                print '(a,g16.5,7a)' &
                    ,'# <tr><td colspan=2><!--atr2-->return period ',xyear, &
                    '</td><td>',atx(1), &
                    '</td><td>',atx25(1),' ... ',atx975(1),'</td></tr>'
                print '(7a)' &
                    ,'# <tr><td colspan=2><!--atr2-->probability </td><td>',ainvtx(1), &
                    '</td><td>',ainvtx975(1),' ... ',ainvtx25(1),'</td></tr>'
            end if
        else
            if ( .not. lnone ) then
                print '(a,i4,a,i4,5a)' &
                    ,'# return time ',year,' at yr=',yr1a &
                    ,' = ',atx(1),' 95% CI ',atx25(1),atx975(1)
                print '(a,i4,a,i4,5a)' &
                    ,'# return time ',year,' at yr=',yr2a &
                    ,' = ',atx(2),' 95% CI ',atx25(2),atx975(2)
                print '(a,i4,5a)','# return time ',year,' ratio = ' &
                    ,atx(3),' 95% CI ',atx25(3),atx975(3)
                if ( yr2b > 0 ) then
                    print '(a,i4,a,i4,5a)' &
                        ,'# return time ',year,' at yr=',yr2b &
                        ,' = ',atx(4),' 95% CI ',atx25(4),atx975(4)
                end if
            else
                print '(a,g16.5,7a)','# return time ',xyear,' = ',atx(1),' 95% CI ',atx25(1),atx975(1)                
            end if
        end if
        if ( plot ) then    ! output for stationlist
            write(11,'(3g20.4,i6,a)') tx(1),tx25(1),tx975(1),yr1a,' return time'
            write(11,'(3g20.4,i6,a)') tx(2),tx25(2),tx975(2),yr2a,' return time'
            write(11,'(3g20.4,a)') tx(3),tx25(3),tx975(3),' ratio'
        end if
    endif
end subroutine printcovreturntime


subroutine printcovpvalue(txtx,nmc,nens,lweb,plot)

!   print out the p-value at which the ratio of return time is unequal to 1 (FAR unequal to 0)

    implicit none
    integer :: nmc,nens
    real :: txtx(nmc,3)
    logical :: lweb,plot
    integer :: i
    real :: p,one
    one = 1
    call invgetcut(p,one,nens,txtx(1,3))
    if ( p > 0.5 ) p = 1-p
    if ( lweb ) then
        print '(2a,f7.4,2a)','# <tr><td><i>p</i>-value probability ratio (one-sided)', &
            '</td><td>&#8800; 1</td><td>',p,'</td><td>&nbsp;</td>', &
            '</tr>'
    else
        print '(a,f7.4,a)','# p-value for ratio/=1 (one-sided) ',p
    end if
    if ( plot ) then    ! output for stationlist
        write(11,'(g20.4,a)') p,' p-value (one-sided)'
    end if
end subroutine printcovpvalue


subroutine plotreturnvalue(ntype,t25,t975,n)

!   print return value at a few times

    implicit none
    integer :: n
    real :: t25(10),t975(10)
    logical :: lweb
    integer :: i,ntype,init
    save init
    real :: x,f
    character fx*10,ff*10
    data init /0/

    if ( ntype == 2 .or. ntype == 3 .or. ntype == 4 ) then ! extreme value  plot
        if ( init == 0 ) then
            init = 1
            if ( ntype == 2 ) then
                fx = 'Gumbel(T)'
            else if ( ntype == 3 ) then
                fx = 'log(T)'
            else if ( ntype == 4 ) then
                fx = 'sqrtlog(T)'
            else
                fx = '???'
            end if
            ff = 'fit'      ! maybe later propagate nfit and make beautiful...
            print '(5a)','#     n            ',fx, &
                '            Y                     ',ff, &
                '            T              date'
        end if
        do i=1,10
            if ( mod(i,3) == 1 ) then
                x = 1 + i/3
            elseif ( mod(i,3) == 2 ) then
                x = 1 + log10(2.) + i/3
            else
                x = log10(5.) + i/3
            endif
            x = x + log10(real(n))
            f = 1 - 10.**(-x)
            call printpoint(0,f,ntype,-999.9,t25(i),0)
        enddo
        print '(a)'
        do i=1,10
            if ( mod(i,3) == 1 ) then
                x = 1 + i/3
            elseif ( mod(i,3) == 2 ) then
                x = 1 + log10(2.) + i/3
            else
                x = log10(5.) + i/3
            endif
            x = x + log10(real(n))
            f = 1 - 10.**(-x)
            call printpoint(0,f,ntype,-999.9,t975(i),0)
        enddo
        print '(a)'
    endif
end subroutine plotreturnvalue


subroutine val_or_inf(atx,tx,lweb)
    implicit none
    real :: tx
    character atx*(*)
    logical :: lweb
    if ( abs(tx) < 1e19 .and. abs(tx) > 1e-19 ) then
        write(atx,'(g16.5)') tx
    else
        if ( abs(tx) >= 1e33 .or. abs(tx) < 1e-33 ) then
            atx = 'undefined'
        else if ( abs(tx) >= 1e19 ) then
            if ( lweb ) then
                atx = '&infin;'
            else
                atx = 'infinity'
            end if
        else
            atx = '0'
        end if
        if ( tx < 0 ) then
            atx = '-'//trim(atx)
        end if
    end if
end subroutine val_or_inf


subroutine plot_ordered_points(xx,xs,yrs,ntot,ntype,nfit, &
    frac,a,b,xi,j1,j2,minindx,mindata,pmindata, &
    yr2a,xyear,snorm,lchangesign,lwrite,last)

!   Gumbel or (sqrt) logarithmic plot

    implicit none
    integer :: ntot,ntype,nfit,j1,j2,yr2a
    integer :: yrs(0:ntot)
    real :: xx(ntot),xs(ntot),frac,a,b,xi
    real :: minindx,mindata,pmindata,xyear,snorm
    logical :: lchangesign,lwrite,last
    integer :: i,j,it,nzero
    real :: f,ff,s,x,z,sqrt2,tmax
    real,save :: smin=3e33,smax=-3e33
    character string*100
    logical :: lprinted
    real :: erf
    real,external :: gammp,gammq,invcumgamm,invcumpois,serfi
    real :: scalingpower
    common /c_scalingpower/ scalingpower

    if ( lwrite ) then
        print *,'plot_ordered_points: xyear = ',xyear
        print *,'                    a,b,xi = ',a,b,xi
        print *,'  minindx,mindata,pmindata = ',minindx,mindata,pmindata
    end if
    call keepalive1('Output',0,ntot+100)
    call nrsort(ntot,xx)
    do nzero = 1,ntot
        if ( xx(nzero) /= 0 ) exit
    end do
    nzero = nzero - 1
    do i=1,ntot+100
        if ( mod(i,1000) == 0 ) call keepalive1('Output',i,ntot+100)
    
    !   search unsorted year that belongs to this point
    !   attention: quadratic algorithm!
    
        if ( i <= ntot ) then
            if ( xx(i) < 1e33 ) then
                do j=1,ntot
                    if ( xs(j) == xx(i) ) goto 790
                end do
                print *,'warning: cannot find year for x = ',x
            end if
        end if
        j = 0 ! yrs(0) = 0 is orinted
790     continue
        if ( j > 0 ) xs(j) = 3e33

        if ( i <= ntot ) then
            f = real(i)/real(ntot+1)
            if ( xx(i) < 1e33 ) then
                x = xx(i)
            else
                x = -999.9
            end if
        else
            f = 1 - 1/real(ntot+1)*0.9**(i-ntot)
            if ( abs(f-1) < 1e-6 ) goto 800
            x = -999.9
        endif
        if ( frac /= 1 ) f = 1-(1-f)*frac
        if ( nfit == 0 ) then
            s = -999.9
        elseif ( nfit == 1 ) then
        !   Poisson distribution - only the last point of a bin makes sense
            if ( i > ntot .or. xx(min(ntot,i)) /= xx(min(ntot,i+1)) ) then
                f = snorm*f
                if ( minindx > 0 ) then
                    f = f + gammq(minindx+0.5,a)
                endif
                s = invcumpois(f,1-f,a)
            else
                s = -999.9
            endif
        elseif ( nfit == 2 ) then
        !   Gaussian distribution
            sqrt2 = sqrt(2.)
            ff = 2*snorm*f
            if ( minindx > -1e33 ) then
                ff = ff + erf((minindx-a)/(sqrt2*b))
            else
                ff = ff - 1
            endif
        !   netlib routine
            z = serfi(ff)
            s = a + sqrt2*b*z
            if ( lchangesign ) s = -s
        elseif ( nfit == 3 ) then
        !   Gamma distribution plus possible delta at zero
            if ( i <= nzero ) then
                s = 0
            else
                ff = (f*ntot - nzero)/(ntot - nzero)
                ff = snorm*ff
                if ( minindx > 0 ) then
                    ff = ff + gammp(a,minindx/b)
                endif
                if ( lchangesign ) ff = 1 - ff
                s = invcumgamm(ff,1-ff,a,abs(b))
            end if
        elseif ( nfit == 4 ) then
        !   Gumbel distribution
            s = snorm*f
            if ( minindx > -1e33 ) then
                s = s + exp(-exp(-(minindx-a)/b))
            endif
            s = a - b*log(-log(s))
            if ( lchangesign ) s = -s
        elseif ( nfit == 5 ) then
        !   GEV distribution
            s = snorm*f
            if ( minindx > -1e33 ) then
                s = s + exp(-(1+(minindx-a)/b)**(-1/xi))
            endif
            if ( xi == 0 ) then
                s = a - b*log(-log(s))
            else
                s = a + b/xi*((-log(s))**(-xi)-1)
            endif
            if ( lchangesign ) s = -s
        elseif ( nfit == 6 ) then
        !   GPD distribution - only in the tail
            if ( f < pmindata/100 ) then
                s = -999.9
            else
                s = (f - pmindata/100)/(1 - pmindata/100)
                if ( abs(xi) < 1e-4 ) then
                    s = b*(-log(1-s) + 0.5*xi*(log(1-s))**2)
                else
                    s = b/xi*((1-s)**(-xi) - 1)
                endif
                s = mindata + s
                if ( lchangesign ) s = -s
            endif
        else
            write(0,*) 'histogram: error: unknown distribution ',nfit
            call exit(-1)
        endif
    
        if ( lchangesign ) then
            if ( x /= -999.9 .and. x < 1e33 ) x = -x
        endif
        call printpoint(i,f,ntype,x,s,yrs(j))
        smin = min(s,smin)
        smax = max(s,smax)
        if ( i > ntot .and. (1-f)*(j2-j1+1) < 0.0001 ) goto 800
    enddo
800 continue
    if ( last ) then
        print '(a)'
        print '(a)'
        if ( xyear /= 3e33 ) then
            call printpoint(0,1/real(ntot+1),ntype,-999.9,xyear,10000*yr2a)
            call printpoint(0,f,ntype,-999.9,xyear,10000*yr2a)
        else
            call printpoint(0,1/real(ntot+1),ntype,-999.9,-999.9,10000*yr2a)
            call printpoint(0,f,ntype,-999.9,-999.9,10000*yr2a)
        endif
        call print_xtics(6,ntype,ntot,j1,j2)
        if ( scalingpower /= 1 ) then
            call printy2tics(6,smin,smax,scalingpower)
        end if
    end if
end subroutine plot_ordered_points


subroutine print_bootstrap_message(ndecor,j1,j2)
    implicit none
    integer :: ndecor,j1,j2
    if ( ndecor == 1 ) then
        print '(3a)' &
            ,'# The error margins were computed with a ' &
            ,'bootstrap method that assumes all points are ' &
            ,'temporally independent.'
    elseif ( j1 /= j2 ) then
        print '(2a,i4,a)' &
            ,'# The error margins were computed with a ' &
            ,'moving block bootstrap with block size ' &
            ,ndecor,' months.'
    else
        print '(2a,i4,a)' &
            ,'# The error margins were computed with a ' &
            ,'moving block bootstrap with block size ' &
            ,ndecor,' years.'
    endif
end subroutine print_bootstrap_message


subroutine print_xtics(unit,ntype,ntot,j1,j2)

!   convert Gumbel/log variates to return periods

    implicit none
    integer :: unit,ntype,ntot,j1,j2
    integer :: tmax,it,i
    real :: f,t,xmax
    character string*80
    logical :: lprinted

    tmax = max(10000,int((ntot+1)/real(j2-j1+1))) ! in years
    if ( ntype == 2 ) then
        write(unit,'(a,f8.4,a)') '#@ set xrange [:', &
            -log(-log(1-1/(tmax*real(j2-j1+1)+1))),']'
    elseif ( ntype == -2 ) then
        t = (tmax*real(j2-j1+1)+1)/10 ! 10 times smaller range as we go two sides
        xmax = log(t)
        write(unit,'(a,f8.4,a,f8.4,a)') '#@ set xrange [',-xmax,':',xmax,']'
    elseif ( ntype == 3 ) then
        write(unit,'(a,f8.4,a)') '#@ set xrange [:', &
            -log(1/(tmax*real(j2-j1+1)+1)),']'
    elseif ( ntype == 4 ) then
        write(unit,'(a,f8.4,a)') '#@ set xrange [:', &
            sqrt(-log(1/(tmax*real(j2-j1+1)+1))),']'
    else
        write(0,*) 'histogram: error: unknown ntype ',ntype
        call exit(-1)
    endif
    write(unit,'(a)') '#@ set xtics (\\'
    it = 1
    do i=1,100
        f = 1-1/(it*real(j2-j1+1))
        t = it*real(j2-j1+1)
        write(string,'(a,i1,a)') '(a,i',1+int(log10(real(it))),',a,f8.4,a)'
        lprinted = .TRUE. 
        if ( ntype == 2 ) then
            if ( f > 0 ) then
                if ( it > 10 .and. mod(i,3) /= 1 ) then
                    write(unit,'(a,f8.4,a)') '#@ "" ',-log(-log(f)),'\\'
                else
                    write(unit,string) '#@ "',it,'" ',-log(-log(f)),'\\'
                endif
            else
                lprinted = .FALSE. 
            endif
        elseif ( ntype == -2 ) then
            if ( mod(i,3) /= 1 ) then
                write(unit,'(a,f8.4,a)') '#@ "" ',log(t),',\\'
                write(unit,'(a,f8.4,a)') '#@ "" ',log(t),'\\'
            else
                write(unit,string) '#@ "',it,'" ',log(t),',\\'
                write(unit,string) '#@ "1/',-it,'" ',log(t),'\\'
            endif
        elseif ( ntype == 3 ) then
            if ( it > 10 .and. mod(i,3) /= 1 ) then
                write(unit,'(a,f8.4,a)') '#@ "" ',-log(1-f),'\\'
            else
                write(unit,string) '#@ "',it,'" ',-log(1-f),'\\'
            endif
        elseif ( ntype == 4 ) then
            if ( it > 10 .and. mod(i,3) /= 1 ) then
                write(unit,'(a,f8.4,a)') '#@ "" ',sqrt(-log(1-f)),'\\'
            else
                write(unit,string) '#@ "',it,'" ',sqrt(-log(1-f)),'\\'
            endif
        else
            write(0,*) 'histogram: error: unknown ntype ',ntype
            call exit(-1)
        endif
        if ( mod(i,3) == 2 ) then
            it = nint(2.5*it)
        else
            it = 2*it
        endif
        if ( it > tmax ) goto 801
        if ( lprinted ) write(unit,'(a)') '#@ ,\\'
    enddo
801 continue
    write(unit,'(a)') '#@ )'
end subroutine print_xtics


subroutine printy2tics(unit,smin,smax,scalingpower)

!   print the gnuplot command to replace the RH Y-axis with the original units

    implicit none
    integer :: unit
    real :: smin,smax,scalingpower
    integer :: i,j,n,pow,i1,i2
    real :: delta,s
    character string*100

    if ( scalingpower == 1 ) return
    if ( smin < 0 ) smin = 0
    if ( scalingpower == 0 ) then
        smin = 10.**smin
        smax = 10.**smax
    else
        smin = smin**(1/scalingpower)
        smax = smax**(1/scalingpower)
    end if
    delta = (smax-smin)/5
!   round to 1,2,5*10^pow
    pow = int(log10(delta))
    delta = delta/10**pow
    if ( delta < 2 ) then
        delta = 1
    else if ( delta < 5 ) then
        delta = 2
    else
        delta = 5
    end if
    delta = delta*10**pow
    i1 = int(smin/delta)
    i2 = int(smax/delta) + 1

    write(unit,'(a)') '#@ set ytics nomirror'
    write(unit,'(a)') '#@ set y2tics (\\'
    do i = i1,i2
        if ( nint(i*delta) >= 1 .or. i == 0 ) then
            if ( i == 0 ) then
                string = '(a,i1,a,f15.4,a)'
            else
                n = 1+int(0.01+log10(i*delta))
                n = max(1,min(9,n))
                write(string,'(a,i1,a)') '(a,i',n,',a,f15.4,a)'
            end if
            if ( scalingpower == 0 ) then
                if ( i /= 0 ) then
                    write(unit,string) '#@ "',nint(i*delta), &
                        '" ',log10(i*delta),'\\'
                end if
            else
                write(unit,string) '#@ "',nint(i*delta), &
                    '" ',(i*delta)**(scalingpower),'\\'
            end if
        else
            n = int(-0.01-log10(i*delta))
            n = max(1,min(n,7))
            write(string,'(a,i1,a,i1,a)') '(a,f',n+2,'.',n, &
                ',a,f15.4,a)'
            if ( scalingpower == 0 ) then
                write(unit,string) '#@ "',i*delta, &
                    '" ',log10(i*delta),'\\'
            else
                write(unit,string) '#@ "',i*delta, &
                    '" ',(i*delta)**(scalingpower),'\\'
            end if
        end if
        if ( i < i2 ) then
            write(unit,'(a)') '#@ ,\\'
        end if
    end do
    write(unit,'(a)') '#@ )'
end subroutine printy2tics


subroutine plot_tx_cdfs(txtx,nmc,nens,ntype,j1,j2)

!   make a file to plot CDFs

    integer :: nmc,nens,ntype,j1,j2
    real :: txtx(nmc,4)
    integer :: iens,j
    real :: logtxtx(4)
    write(10,'(a)') '# fraction, return values in the past '// &
        'climate, current climate and difference (sorted), '
    if ( ntype == 2 ) then
        write(10,'(a)') '# repeated as gumbel transforms'
    else if ( ntype == 3 ) then
        write(10,'(a)') '# repeated as logarithmic transforms'
    else if ( ntype == 4 ) then
        write(10,'(a)') '# repeated as sqrt-logarithmic transforms'
    else
        write(0,*) 'plot_tx_cdfs: error: unknown value for ntype ',ntype
        call exit(-1)
    end if
    if ( nens > nmc ) then
        write(0,*) 'plot_tx_cdfs: internal error: nens>nmc ',nens,nmc
        nens = nmc
    end if
    do iens=1,nens
        do j=1,4
            if ( j == 3 ) cycle
            if ( ntype == 2 ) then
                if ( txtx(iens,j) < 1e20 .and. txtx(iens,j) > 1 ) then
                    logtxtx(j) = -log(1-1/((j2-j1+1)*txtx(iens,j)))
                    if ( logtxtx(j) > 0 ) then
                        logtxtx(j) = -log(logtxtx(j))
                    else
                        logtxtx(j) = 3e33
                    end if
                else
                    logtxtx(j) = 1e20
                end if
            else if ( ntype == 3 ) then
                if ( txtx(iens,j) < 1e20 .and. txtx(iens,j) > 0 ) then
                    logtxtx(j) = log(txtx(iens,j)*(j2-j1+1))
                else
                    logtxtx(j) = 1e20
                end if
            else if ( ntype == 4 ) then
                if ( txtx(iens,j) < 1e20 .and. txtx(iens,j) > 1 ) then
                    logtxtx(j) = sqrt(log(txtx(iens,j)*(j2-j1+1)))
                else
                    logtxtx(j) = 1e20
                end if
            end if
        end do
        !   simple logscale for the time being
        if ( txtx(iens,3) < 1e20 .and. txtx(iens,3) > 0 ) then
            logtxtx(3) = log(txtx(iens,3)*(j2-j1+1))
        else
            logtxtx(3) = 1e20
        end if
        write(10,'(f6.4,8g16.5)') iens/real(nens+1), &
            (txtx(iens,j),j=1,3),(logtxtx(j),j=1,3), &
            (txtx(iens,j),j=4,4),(logtxtx(j),j=4,4)
    end do
end subroutine plot_tx_cdfs


subroutine getreturnlevels(a,b,xi,alpha,beta,cov1,cov2,cov3,covreturnlevel,j1,j2,assume,t)
    implicit none
    integer :: j1,j2
    real :: a,b,xi,alpha,beta,cov1,cov2,cov3,t(10,4)
    character assume*(*)
    real,external :: covreturnlevel
    integer :: i
    real :: x,xx
    do i=1,10
        if ( mod(i,3) == 1 ) then
            x = 1 + i/3
        elseif ( mod(i,3) == 2 ) then
            x = 1 + log10(2.) + i/3
        else
            x = log10(5.) + i/3
        endif
        x = x + log10(real(j2-j1+1))
        xx = x ! some routines modify x :-(
        t(i,1) = covreturnlevel(a,b,xi,alpha,beta,xx,cov1)
        xx = x
        t(i,2) = covreturnlevel(a,b,xi,alpha,beta,xx,cov2)
        if ( cov3 < 1e33 ) then
            xx = x
            t(i,4) = covreturnlevel(a,b,xi,alpha,beta,xx,cov3)
        else
            t(i,4) = 3e33
        end if
        t(i,3) = t(i,2) - t(i,1)
        if ( assume == 'scale' ) then ! relative change in % relative to *past*
            t(i,3) = 100*t(i,3)/t(i,1)
        end if
    enddo
end subroutine getreturnlevels


subroutine getreturnyears(a,b,xi,alpha,beta,xyear,cov1,cov2,cov3,covreturnyear,j1,j2,tx, &
        lchangesign,lwrite)
    implicit none
    integer :: j1,j2
    real :: a,b,xi,alpha,beta,xyear,cov1,cov2,cov3,tx(4)
    logical :: lchangesign,lwrite
    real,external :: covreturnyear
    tx(1) = covreturnyear(a,b,xi,alpha,beta,xyear,cov1,lchangesign)
    tx(2) = covreturnyear(a,b,xi,alpha,beta,xyear,cov2,lchangesign)
    if ( cov3 < 1e33 ) then
        tx(4) = covreturnyear(a,b,xi,alpha,beta,xyear,cov3,lchangesign)
    end if
    if ( tx(2) > 1e19 ) then
        if ( tx(1) > 1e19 ) then
            tx(3) = 3e33
        else
            tx(3) = 1e20
        end if
    else
        if ( tx(1) > 1e19 ) then
            tx(3) = 1e20
        else
            tx(3) = tx(1) / tx(2)
        end if
    end if
    if ( lwrite ) then
        print *,'return time = ',tx
    end if
end subroutine getreturnyears


subroutine getabfromcov(a,b,alpha,beta,cov,aa,bb)
    implicit none
    real :: a,b,alpha,beta,cov,aa,bb
    real :: arg
    character cassume*5
    common /fitdata4/ cassume
    integer :: ncur
    real :: restrain
    logical :: llwrite,llchangesign
    common /fitdata2/ restrain,ncur,llwrite,llchangesign
            
    if ( alpha > 1e33 ) then
        aa = a
        bb = b
    else if ( cassume == 'shift' ) then
        aa = a + alpha*cov
        bb = b
    else if ( cassume == 'scale' ) then
        arg = alpha*cov/a
        if ( arg < 70 .and. a < 1e33 .and. b < 1e33 ) then
            aa = a*exp(alpha*cov/a)
            bb = b*aa/a
        else
            aa = 3e33
            bb = 3e33
        end if
    else if ( cassume == 'both' ) then
        aa = a + alpha*cov
        bb = b + beta*cov
    else
        write(0,*) 'getabfromcov: error: unknown value for assume: ',cassume
        call exit(-1)
    end if
end subroutine getabfromcov


subroutine printab(restrain,lnone,lweb)
    implicit none
    real,intent(in) :: restrain
    logical,intent(in) :: lnone,lweb
    character :: string*80
    character cassume*5
    common /fitdata4/ cassume

    string = ' '
    if ( restrain /= 0 ) then 
        if ( lweb ) then
            write(string,'(a,f4.1)') ' and a Gaussian penalty on &xi; of width ',restrain/2
        else
            write(string,'(a,f4.1)') ' and a Gaussian penalty on xi of width ',restrain/2
        end if
    end if
    if ( lweb ) then
        if ( lnone ) then
            if ( string /= ' ' ) then
                print '(a)','# <tr><td colspan="4">'//trim(string)//'</td></tr>'
            endif
        else if ( cassume == 'shift' ) then
            print '(a)','# <tr><td colspan="4">'// &
                'with &mu;''= &mu;+&alpha;T '// &
            '   and &sigma;'' = &sigma;'//trim(string)//'</td></tr>'
        else if ( cassume == 'scale' ) then
            print '(a)','# <tr><td colspan="4">'// &
                'with &mu;'' = &mu; exp(&alpha;T/&mu;) and '// &
                '&sigma;'' = &sigma; exp(&alpha;T/&mu;)'//trim(string)//'</td></tr>'
        else if ( cassume == 'both' ) then
            print '(a)','# <tr><td colspan="4">'// &
                'with &mu;''= &mu;+&alpha;T '// &
                'and &sigma;'' = &sigma;+&beta;T'//trim(string)//'</td></tr>'
        else
            write(0,*) 'printab: error: unknow value for assume ',cassume
        end if
    else
        if ( lnone ) then
            if ( string /= ' ' ) then
                print '(a)',trim(string)
            end if
        else if ( cassume == 'shift' ) then
            print '(a)','a'' = a+alpha*T, b''= b'//trim(string)
        else if ( cassume == 'scale' ) then
            print '(a)','a'' = a*exp(alpha*T/a), b''= b*a''/a'//trim(string)
        else if ( cassume == 'both' ) then
            print '(a)','a'' = a+alpha*T, b''= b+beta*T'//trim(string)
        else
            write(0,*) 'printab: error: unknow value for assume ',cassume
        end if
    end if
end subroutine printab


subroutine adjustyy(ntot,xx,assume,a,b,alpha,beta,cov,yy,zz,aaa,bbb,lchangesign,lwrite)

!   input: xx,assume,a,b,alpha,beta,cov
!   output: yy,zz,aaa,bbb
!   flag: lwrite

    implicit none
    integer :: ntot
    real :: xx(2,ntot),yy(ntot),zz(ntot)
    real :: a,b,alpha,beta,cov,aaa,bbb
    character assume*(*)
    logical :: lchangesign,lwrite
    integer :: i
    real :: arg

    if ( lwrite ) then
        print *,'adjustyy: input ',assume
        print *,'a,b,alpha,beta,cov = ',a,b,alpha,beta,cov
        do i=1,min(5,ntot)
            print *,i,xx(1,i),xx(2,i)
        end do
    end if
    do i=1,ntot
        yy(i) = xx(1,i)
        zz(i) = xx(2,i)
    end do
    if ( alpha > 1e33 ) then
        aaa = a
        bbb = b
    else if ( assume == 'shift' ) then
        do i=1,ntot
            yy(i) = yy(i) - alpha*(zz(i)-cov)
        end do
        aaa = a+cov*alpha
        bbb = b
    else if ( assume == 'scale' ) then
        do i=1,ntot
            arg = -alpha*(zz(i)-cov)/a
            if ( arg < 70 ) then
                yy(i) = yy(i)*exp(arg)
            else
                yy(i) = 3e33
            end if
        end do
        arg = alpha*cov/a
        if ( arg < 70 ) then
            aaa = a*exp(arg)
            bbb = b*exp(arg)
        else
            aaa = 3e33
            bbb = 3e33
        end if
    else if ( assume == 'both' ) then
        write(0,*) 'adjustyy: error: not yet ready'
        write(*,*) 'adjustyy: error: not yet ready'
        call exit(-1)
    else
        write(0,*) 'adjustyy: error: unknown assumption ',assume
        write(*,*) 'adjustyy: error: unknown assumption ',assume
        call exit(-1)
    end if
end subroutine adjustyy


subroutine write_obscov(xx,yrs,ntot,xmin,cov2,xyear,year,offset,lchangesign)
    implicit none
    integer :: ntot,yrs(0:ntot),year
    real :: xx(2,ntot),xmin,cov2,xyear,offset
    logical :: lchangesign
    integer :: i,is
    logical :: lopen
    character string*1000,arg*250
    inquire(unit=15,opened=lopen)
    if ( lopen ) then
        if ( lchangesign ) then
            is = -1
        else
            is = +1
        end if
        string = ' '
        do i=0,command_argument_count()
            call get_command_argument(i,arg)
            string = trim(string)//' '//arg(index(arg,'/',.TRUE.)+1:)
        end do
        write(15,'(a)') '# ' // trim(string)
        write(15,'(a)') '# covariate  value'
        do i=1,ntot
            if ( xx(1,i) > xmin ) then
                write(15,'(2g20.6,i11)') xx(2,i)+offset,is*xx(1,i),yrs(i)
            end if
        end do
        write(15,'(a)')
        write(15,'(a)')
        if ( xyear < 1e33 ) then
            write(15,'(2g20.6,i11)') cov2+offset,is*xyear,year
        else
            write(15,'(g20.6,a,i11)') cov2+offset,'-999.900',0
        end if
    end if
end subroutine write_obscov


subroutine write_threshold(cmin,cmax,a,b,xi,alpha,beta,offset,lchangesign,covreturnlevel)
    implicit none
    real :: cmin,cmax,a,b,xi,alpha,beta,offset
    real,external :: covreturnlevel
    logical :: lchangesign
    integer :: i,is
    real :: c,aa,bb,x6,x40,v6,v40
    logical :: lopen
    inquire(unit=15,opened=lopen) ! no flag in getopts.inc
    if ( lopen ) then
        if ( lchangesign ) then
            is = -1
        else
            is = +1
        end if
        write(15,'(a)')
        write(15,'(a)')
        write(15,'(a)') '# covariate threshold/position value(6yr) val(40yr)'
        do i=0,100
            c = cmin + (cmax-cmin)*i/100.
            call getabfromcov(a,b,alpha,beta,c,aa,bb)
            x6 = log10(6.)
            v6 = covreturnlevel(a,b,xi,alpha,beta,x6,c)
            x40 = log10(40.)
            v40 = covreturnlevel(a,b,xi,alpha,beta,x40,c)
            write(15,'(4g20.6)') c+offset,is*aa,is*v6,is*v40
        end do
    endif
end subroutine write_threshold


subroutine write_dthreshold(cov1,cov2,cov3,acov,offset,lchangesign)
    implicit none
    real :: cov1,cov2,cov3,acov(3,3),offset
    logical :: lchangesign
    integer :: i,j,is
    real :: a(3)
    logical :: lopen
    
    inquire(unit=15,opened=lopen) ! no flag in getopts.inc
    if ( lopen ) then
        is = +1 ! pre-flipped
        write(15,'(a)')
        write(15,'(a)')
        write(15,'(a)') '# covariate locationparameter/threshold lowerbound upperbound'
        do j=1,3
            if ( acov(j,1) < 1e33 ) then
                a(j) = is*acov(j,1)
            else
                a(j) = -999.9
            end if
        end do
        write(15,'(4g20.6)') cov1+offset,a
        do j=1,3
            if ( acov(j,2) < 1e33 ) then
                a(j) = is*acov(j,2)
            else
                a(j) = -999.9
            end if
        end do
        write(15,'(4g20.6)') cov2+offset,a
        if ( cov3 < 1e33 ) then
            do j=1,3
                if ( acov(j,3) < 1e33 ) then
                    a(j) = is*acov(j,3)
                else
                    a(j) = -999.9
                end if
            end do
            write(15,'(4g20.6)') cov3+offset,a
        end if
    endif
end subroutine write_dthreshold