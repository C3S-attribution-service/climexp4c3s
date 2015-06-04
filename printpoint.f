*  #[ printpoint:
        subroutine printpoint(i,f,ntype,x,s,yr)
        implicit none
        integer i,ntype,yr
        real f,x,s
        if ( ntype.eq.2 ) then
            if ( f.gt.0 ) then
                print '(i8,4g22.6,i12)',i,
     +               -log(-log(f)),x,s,1/(1-f),yr
            endif
        elseif ( ntype.eq.3 ) then
            print '(i8,4g22.6,i12)',i,-log(1-f),x,s,1/(1-f),yr
        elseif ( ntype.eq.4 ) then
            print '(i8,4g22.6,i12)',i,sqrt(-log(1-f)),x,s,1/(1-f),yr
        else
            write(0,*) 'histogram: error: unknown ntype '
     +           ,ntype
            call abort
        endif
        end
*  #] printpoint:
*  #[ printreturnvalue:
        subroutine printreturnvalue(ntype,t,t25,t975,lweb)
*
*       print return times
*
        implicit none
        integer ntype
        real t(10),t25(10),t975(10)
        logical lweb
        integer i
*
        if ( ntype.eq.2 .or. ntype.eq.3 .or. ntype.eq.4 ) then ! extreme value  plot
            if ( lweb ) then
                do i=1,10,3
                    print '(a,i5,a,f16.3,a,f16.3,a,f16.3,a)'
     +                   ,'# <tr><td>return value ',10**(1+i/3)
     +                   ,' yr</td><td>',t(i),'</td><td>',t25(i)
     +                   ,' ... ',t975(i),'</td></tr>'
                    call print3untransf(t(i),t25(i),t975(i),0)
                enddo
            else
                do i=1,4
                    print '(a,i5,a,f16.3,a,2f16.3)'
     +                   ,'# value for return period ',10**i,' year: '
     +                   ,t(i),' 95% CI ',t25(i),t975(i)
                    call print3untransf(t(i),t25(i),t975(i),0)
                enddo
            endif
        endif
        end
*  #] printreturnvalue:
*  #[ printcovreturnvalue:
        subroutine printcovreturnvalue(ntype,t,t25,t975,yr1a,yr2a,lweb,
     +       plot)
*
*       print return times
*
        implicit none
        integer ntype,yr1a,yr2a
        real t(10,3),t25(10,3),t975(10,3)
        logical lweb,plot
        integer i
*
        if ( ntype.eq.2 .or. ntype.eq.3 .or. ntype.eq.4 ) then ! extreme value  plot
            if ( lweb ) then
                do i=1,10,3
                    print '(a,i5,a,i4,a,f16.3,a,f16.3,a,f16.3,a)'
     +                   ,'# <tr><td>return value ',10**(1+i/3)
     +                   ,' yr</td><td>',yr1a,'</td><td>'
     +                   ,t(i,1),'</td><td>',t25(i,1)
     +                   ,' ... ',t975(i,1),'</td></tr>'
                    print '(a,i4,a,f16.3,a,f16.3,a,f16.3,a)'
     +                   ,'# <tr><td>&nbsp;</td><td>',yr2a
     +                   ,'</td><td>',t(i,2),'</td><td>',t25(i,2)
     +                   ,' ... ',t975(i,2),'</td></tr>'
                    if ( i.eq.10 ) then ! they are all the same by construction...
                        print '(a,f16.3,a,f16.3,a,f16.3,a)'
     +                       ,'# <tr><td>&nbsp;</td><td>diff</td><td>'
     +                       ,t(i,3),'</td><td>',t25(i,3)
     +                       ,' ... ',t975(i,3),'</td></tr>'
                    end if
                    call print3untransf(t(i,1),t25(i,1),t975(i,1),yr1a)
                    call print3untransf(t(i,2),t25(i,2),t975(i,2),yr2a)
                enddo
                if ( plot ) then ! output for stationlist
                    do i=1,10
                        write(11,'(4g20.4,i6,a)') t(i,1),t25(i,1),
     +                       t975(i,1),10**(1+i/3),yr1a,' return value'
                        write(11,'(4g20.4,i6,a)') t(i,2),t25(i,2),
     +                       t975(i,2),10**(1+i/3),yr2a,' return value'
                    end do
                end if
            else
                do i=1,4
                    print '(a,i5,a,i4,a,f16.3,a,2f16.3)'
     +                   ,'# value for return period ',10**i
     +                   ,' year at yr=',yr1a,': '
     +                   ,t(i,1),' 95% CI ',t25(i,1),t975(i,1)
                    print '(a,i5,a,i4,a,f16.3,a,2f16.3)'
     +                   ,'# value for return period ',10**i
     +                   ,' year at yr=',yr2a,': '
     +                   ,t(i,2),' 95% CI ',t25(i,2),t975(i,2)
                    print '(a,i5,a,f16.3,a,2f16.3)'
     +                   ,'# value for return period ',10**i
     +                   ,' year, difference: '
     +                   ,t(i,3),' 95% CI ',t25(i,3),t975(i,3)
                    call print3untransf(t(i,1),t25(i,1),t975(i,1),yr1a)
                    call print3untransf(t(i,2),t25(i,2),t975(i,2),yr2a)
                enddo
            endif
        endif
        end
*  #] printcovreturnvalue:
*  #[ printreturntime:
        subroutine printreturntime(year,xyear,tx,tx25,tx975,lweb)
*
*       print return time of year
*
        implicit none
        integer year
        real xyear,tx,tx25,tx975
        logical lweb
        
        if ( xyear.lt.1e33 ) then
            if ( lweb ) then
                if ( year.eq.9999 ) then
                    print '(a,g16.5,a,g16.5,a,g16.5,a,g16.5,a)'
     +                   ,'# <tr><td>return period ',xyear,'</td><td>'
     +                   ,tx,'</td><td>',tx25,' ... ',tx975,'</td></tr>'
                else
                    print '(a,g16.5,a,i5,a,g16.5,a,g16.5,a,g16.5,a)'
     +                   ,'# <tr><td>return period ',xyear,'(',year
     +                   ,')</td><td>'
     +                   ,tx,'</td><td>',tx25,' ... ',tx975,'</td></tr>'
                end if
            else
                print '(a,f16.5,a,i4,a,f16.5,a,2f18.5)','# return time '
     +               ,xyear,' (',year,') = ',tx,' 95% CI ',tx25,tx975
            endif
        endif
        end
*  #] printreturntime:
*  #[ printcovreturntime:
        subroutine printcovreturntime(year,xyear,idmax,tx,tx25,tx975,
     +       yr1a,yr2a,lweb,plot)
*
*       print return time of year at cov1 and cov2
*
        implicit none
        integer year,yr1a,yr2a
        real xyear,tx(3),tx25(3),tx975(3)
        logical lweb,plot
        character idmax*(*)
        integer i
        character atx(3)*16,atx25(3)*16,atx975(3)*16

        if ( xyear.lt.1e33 ) then
            do i=1,3
                call val_or_inf(atx(i),tx(i),lweb)
                call val_or_inf(atx25(i),tx25(i),lweb)
                call val_or_inf(atx975(i),tx975(i),lweb)
            end do
            if ( lweb ) then
                print '(a,i5,a,i5,7a)'
     +               ,'# <tr><td><!--atr2-->return period ',year,
     +               '</td><td>',yr1a,'</td><td>',atx(1),
     +               '</td><td>',atx25(1),' ... ',atx975(1),'</td></tr>'
                print '(a,g16.5,a,i5,7a)'
     +               ,'# <tr><td><!--atr1-->(value ',xyear,')</td><td>'
     +               ,yr2a,'</td><td>',atx(2),
     +               '</td><td>',atx25(2),' ... ',atx975(2),'</td></tr>'
                if ( idmax.eq.' ' ) then
                    print '(8a)'
     +               ,'# <tr><td><!--atra-->&nbsp;</td><td>ratio',
     +               '</td><td>',atx(3),'</td><td>',atx25(3),
     +               ' ... ',atx975(3),'</td></tr>'
                else
                    print '(10a)'
     +               ,'# <tr><td><!--atra-->(at ',trim(idmax),
     +               ')</td><td>ratio',
     +               '</td><td>',atx(3),'</td><td>',atx25(3),
     +               ' ... ',atx975(3),'</td></tr>'
                end if
            else
                print '(a,i4,a,i4,5a)'
     +               ,'# return time ',year,' at yr=',yr1a
     +               ,' = ',atx(1),' 95% CI ',atx25(1),atx975(1)
                print '(a,i4,a,i4,5a)'
     +               ,'# return time ',year,' at yr=',yr2a
     +               ,' = ',atx(2),' 95% CI ',atx25(2),atx975(2)
                print '(a,i4,5a)'
     +               ,'# return time ',year,' ratio = '
     +               ,atx(3),' 95% CI ',atx25(3),atx975(3)
            endif
            if ( plot ) then    ! output for stationlist
                write(11,'(3g20.4,i6,a)') tx(1),tx25(1),tx975(1),yr1a,
     +               ' return time'
                write(11,'(3g20.4,i6,a)') tx(2),tx25(2),tx975(2),yr2a,
     +               ' return time'
                write(11,'(3g20.4,a)') tx(3),tx25(3),tx975(3),' ratio'
            end if
        endif
        end
*  #] printcovreturntime:
*  #[ printcovpvalue:
        subroutine printcovpvalue(txtx,nmc,nens,lweb)
!
!       print out the p-value at which the ratio of return time is unequal to 1 (FAR unequal to 0)
!
        implicit none
        integer nmc,nens
        real txtx(nmc,3)
        logical lweb
        integer i
        real p,one
        one = 1
        call invgetcut(p,one,nens,txtx(1,3))
        if ( p.gt.0.5 ) p = 1-p
        if ( lweb ) then
            print '(2a,f7.4,2a)','<tr><td><i>p</i>-value (one-sided)',
     +           '</td><td>&#8800; 1</td><td>',p,'</td><td>&nbsp;</td>',
     +           '</tr>'
        else
            print '(a,f7.4,a)','# p-value for ratio/=1 (one-sided) ',p
        end if
        end subroutine
*  #] printcovpvalue:
*  #[ plotreturnvalue:
        subroutine plotreturnvalue(ntype,t25,t975,n)
*
*       print return value at a few times
*
        implicit none
        integer n
        real t25(10),t975(10)
        logical lweb
        integer i,ntype,init
        save init
        real x,f
        character fx*10,ff*10
        data init /0/
*
        if ( ntype.eq.2 .or. ntype.eq.3 .or. ntype.eq.4 ) then ! extreme value  plot
            if ( init.eq.0 ) then
                init = 1
                if ( ntype.eq.2 ) then
                    fx = 'Gumbel(T)'
                else if ( ntype.eq.3 ) then
                    fx = 'log(T)'
                else if ( ntype.eq.4 ) then
                    fx = 'sqrtlog(T)'
                else
                    fx = '???'
                end if
                ff = 'fit'      ! maybe later propagate nfit and make beautiful...
                print '(5a)','#     n            ',fx,
     +               '            Y                     ',ff,
     +               '            T              date'
            end if
            do i=1,10
                if ( mod(i,3).eq.1 ) then
                    x = 1 + i/3
                elseif ( mod(i,3).eq.2 ) then
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
                if ( mod(i,3).eq.1 ) then
                    x = 1 + i/3
                elseif ( mod(i,3).eq.2 ) then
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
        end
*  #] plotreturnvalue:
*  #[ val_or_inf:
        subroutine val_or_inf(atx,tx,lweb)
        implicit none
        real tx
        character atx*(*)
        logical lweb
        if ( abs(tx).lt.1e19 ) then
            write(atx,'(g16.5)') tx
        else
            if ( lweb ) then
                atx = '&infin;'
            else
                atx = 'infinity'
            end if
            if ( tx.lt.0 ) then
                atx = '-'//atx
            end if
        end if
        end
*  #] val_or_inf:
*  #[ plot_ordered_points:
        subroutine plot_ordered_points(xx,xs,yrs,ntot,ntype,nfit,
     +       frac,a,b,xi,j1,j2,minindx,mindata,pmindata,
     +       yr2a,xyear,snorm,lchangesign,lwrite,last)
*
*       Gumbel or (sqrt) logarithmic plot
*
        implicit none
        integer ntot,ntype,nfit,j1,j2,yr2a
        integer yrs(0:ntot)
        real xx(ntot),xs(ntot),frac,a,b,xi
        real minindx,mindata,pmindata,xyear,snorm
        logical lchangesign,lwrite,last
        integer i,j,ier,it
        real f,ff,s,x,z,sqrt2,tmax
        real,save :: smin=3e33,smax=-3e33
        character string*100
        logical lprinted
        real erf
        real,external :: gammp,gammq,invcumgamm,invcumpois
        real scalingpower
        common /c_scalingpower/ scalingpower

        if ( lwrite ) then
            print *,'plot_ordered_points: xyear = ',xyear
        end if
        call keepalive1('Output',0,ntot+100)
        call nrsort(ntot,xx)
        do i=1,ntot+100
            if ( mod(i,1000).eq.0 ) call keepalive1('Output',i,ntot+100)
*
*           search unsorted year that belongs to this point
*           attention: quadratic algorithm!
*
            if ( i.le.ntot ) then
                do j=1,ntot
                    if ( xs(j).eq.xx(i) ) goto 790
                enddo
                print *,'warning: cannot find year for x = ',x
            endif
            j = 0
 790        continue
            if ( j.gt.0 ) xs(j) = 3e33

            if ( i.le.ntot .and. j.gt.0 ) then
                f = real(i)/real(ntot+1)
                x = xx(i)
            else
                f = 1 - 1/real(ntot+1)*0.9**(i-ntot)
                if ( abs(f-1).lt.1e-6 ) goto 800
                x = -999.9
            endif
            if ( frac.ne.1 ) f = 1-(1-f)*frac
            if ( nfit.eq.0 ) then
                s = -999.9
            elseif ( nfit.eq.1 ) then
*               Poisson distribution - only the last point
*               of a bin makes sense
                if ( i.gt.ntot .or. 
     +               xx(min(ntot,i)).ne.xx(min(ntot,i+1)) ) then
                    f = snorm*f
                    if ( minindx.gt.0 ) then
                        f = f + gammq(minindx+0.5,a)
                    endif
                    s = invcumpois(f,a)
                else
                    s = -999.9
                endif
            elseif ( nfit.eq.2 ) then
*               Gaussian distribution
                sqrt2 = sqrt(2.)
                ff = 2*snorm*f
                if ( minindx.gt.-1e33 ) then
                    ff = ff + erf((minindx-a)/(sqrt2*b))
                else
                    ff = ff - 1
                endif
*               IMSL routine
                call merfi(ff,z,ier)
                s = a + sqrt2*b*z
            elseif ( nfit.eq.3 ) then
*               Gamma distribution
                f = snorm*f
                if ( minindx.gt.0 ) then
                    f = f + gammp(a,minindx/b)
                endif
                s = invcumgamm(f,a,b)
            elseif ( nfit.eq.4 ) then
*               Gumbel distribution
                s = snorm*f
                if ( minindx.gt.-1e33 ) then
                    s = s + exp(-exp(-(minindx-a)/b))
                endif
                s = a - b*log(-log(s))
            elseif ( nfit.eq.5 ) then
*               GEV distribution
                s = snorm*f
                if ( minindx.gt.-1e33 ) then
                    s = s + exp(-(1+(minindx-a)/b)**(-1/xi))
                endif
                if ( xi.eq.0 ) then
                    s = a - b*log(-log(s))
                else
                    s = a + b/xi*((-log(s))**(-xi)-1)
                endif
            elseif ( nfit.eq.6 ) then
*               GPD distribution - only in the tail
                if ( f.lt.pmindata/100 )
     +               then
                    s = -999.9
                else
                    s = (f - pmindata/100)/(1 - pmindata/100)
                    if ( abs(xi).lt.1e-4 ) then
                        s = b*(-log(1-s) + 0.5*xi*(log(1-s))**2)
                    else
                        s = b/xi*((1-s)**(-xi) - 1)
                    endif
                    s = mindata + s
                endif
            else
                write(0,*) 'histogram: error: '//
     +               'unknown distribution ',nfit
                call abort
            endif
*           
            if ( lchangesign ) then
                if ( x.ne.-999.9 ) x = -x
            endif
            call printpoint(i,f,ntype,x,s,yrs(j))
            smin = min(s,smin)
            smax = max(s,smax)
            if ( i.gt.ntot .and. (1-f)*(j2-j1+1).lt.0.0001 ) 
     +           goto 800
        enddo
 800    continue
        if ( last ) then
            print '(a)'
            print '(a)'
            if ( xyear.ne.3e33 ) then
                call printpoint(0,1/real(ntot+1),ntype,-999.9,xyear
     +               ,10000*yr2a)
                call printpoint(0,f,ntype,-999.9,xyear,10000*yr2a)
            else
                call printpoint(0,1/real(ntot+1),ntype,-999.9,-999.9
     +               ,10000*yr2a)
                call printpoint(0,f,ntype,-999.9,-999.9,10000*yr2a)
            endif
            call print_xtics(6,ntype,ntot,j1,j2)
            if ( scalingpower.ne.1 ) then
                call printy2tics(6,smin,smax,scalingpower)
            end if
        end if
        end
*  #] plot_ordered_points:
*  #[ print_bootstrap_message:
        subroutine print_bootstrap_message(ndecor,j1,j2)
        implicit none
        integer ndecor,j1,j2
        if ( ndecor.eq.1 ) then
            print '(3a)'
     +           ,'# The error margins were computed with a '
     +           ,'bootstrap method that assumes all points are '
     +           ,'independent'
        elseif ( j1.ne.j2 ) then
            print '(2a,i4,a)'
     +           ,'# The error margins were computed with a '
     +           ,'moving block bootstrap with block size '
     +           ,ndecor,' months'
        else
            print '(2a,i4,a)'
     +           ,'# The error margins were computed with a '
     +               ,'moving block bootstrap with block size '
     +           ,ndecor,' years'
        endif
        end
*  #] print_bootstrap_message:
*  #[ print_xtics:
        subroutine print_xtics(unit,ntype,ntot,j1,j2)
* 
*       convert Gumbel/log variates to return periods
*
        implicit none
        integer unit,ntype,ntot,j1,j2
        integer tmax,it,i
        real f,t,xmax
        character string*80
        logical lprinted

        tmax = max(10000,int((ntot+1)/real(j2-j1+1))) ! in years
        if ( ntype.eq.2 ) then
            write(unit,'(a,f8.4,a)') '#@ set xrange [:',
     +           -log(-log(1-1/(tmax*real(j2-j1+1)+1))),']'
        elseif ( ntype.eq.-2 ) then
            t = (tmax*real(j2-j1+1)+1)/10 ! 10 times smaller range as we go two sides
            xmax = log(t)
            write(unit,'(a,f8.4,a,f8.4,a)') '#@ set xrange [',
     +           -xmax,':',xmax,']'
        elseif ( ntype.eq.3 ) then
            write(unit,'(a,f8.4,a)') '#@ set xrange [:',
     +           -log(1/(tmax*real(j2-j1+1)+1)),']'
        elseif ( ntype.eq.4 ) then
            write(unit,'(a,f8.4,a)') '#@ set xrange [:',
     +           sqrt(-log(1/(tmax*real(j2-j1+1)+1))),']'
        else
            write(0,*) 'histogram: error: unknown ntype ',ntype
            call abort
        endif
        write(unit,'(a)') '#@ set xtics (\\'
        it = 1
        do i=1,100
            f = 1-1/(it*real(j2-j1+1))
            t = it*real(j2-j1+1)
            write(string,'(a,i1,a)') '(a,i',
     +           1+int(log10(real(it))),',a,f8.4,a)'
            lprinted = .true.
            if ( ntype.eq.2 ) then
                if ( f.gt.0 ) then
                    if ( it.gt.10 .and. mod(i,3).ne.1 ) then
                        write(unit,'(a,f8.4,a)') '#@ "" ',-log(-log(f)
     +                       ),'\\'
                    else
                        write(unit,string) '#@ "',it,'" ',-log(-log(f)
     +                       ),'\\'
                    endif
                else
                    lprinted = .false.
                endif
            elseif ( ntype.eq.-2 ) then
                if ( mod(i,3).ne.1 ) then
                    write(unit,'(a,f8.4,a)') '#@ "" ',
     +                   log(t),',\\'
                    write(unit,'(a,f8.4,a)') '#@ "" ',
     +                   log(t),'\\'
                else
                    write(unit,string) '#@ "',it,'" ',
     +                   log(t),',\\'
                    write(unit,string) '#@ "1/',-it,'" ',
     +                   log(t),'\\'
                endif
            elseif ( ntype.eq.3 ) then
                if ( it.gt.10 .and. mod(i,3).ne.1 ) then
                    write(unit,'(a,f8.4,a)') '#@ "" ',-log(1-f),
     +                   '\\'
                else
                    write(unit,string) '#@ "',it,'" ',-log(1-f),
     +                   '\\'
                endif
            elseif ( ntype.eq.4 ) then
                if ( it.gt.10 .and. mod(i,3).ne.1 ) then
                    write(unit,'(a,f8.4,a)') '#@ "" ',sqrt(-log(1-f)),
     +                   '\\'
                else
                    write(unit,string) '#@ "',it,'" ',sqrt(-log(1-f)),
     +                   '\\'
                endif
            else
                write(0,*) 'histogram: error: unknown ntype '
     +               ,ntype
                call abort
            endif
            if ( mod(i,3).eq.2 ) then
                it = nint(2.5*it)
            else
                it = 2*it
            endif
            if ( it.gt.tmax ) goto 801
            if ( lprinted ) write(unit,'(a)') '#@ ,\\'
        enddo
 801    continue
        write(unit,'(a)') '#@ )'
        end
*  #] print_xtics:
*  #[ print_y2tics:
        subroutine printy2tics(unit,smin,smax,scalingpower)
!
!       print the gnuplot command to replace the RH Y-axis with the original units
!
        implicit none
        integer unit
        real smin,smax,scalingpower
        integer i,j,n,pow,i1,i2
        real delta,s
        character string*100

        if ( scalingpower.eq.1 ) return
        if ( smin.lt.0 ) smin = 0
        if ( scalingpower.eq.0 ) then
            smin = exp(smin)
            smax = exp(smax)
        else
            smin = smin**(1/scalingpower)
            smax = smax**(1/scalingpower)
        end if
        delta = (smax-smin)/5
        ! round to 1,2,5*10^pow
        pow = int(log10(delta))
        delta = delta/10**pow
        if ( delta.lt. 2 ) then
            delta = 1
        else if ( delta.lt.5 ) then
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
            if ( nint(i*delta).ge.1 .or. i.eq.0 ) then
                if ( i.eq.0 ) then
                    string = '(a,i1,a,f15.4,a)'
                else
                    n = 1+int(0.01+log10(i*delta))
                    n = max(1,min(9,n))
                    write(string,'(a,i1,a)') '(a,i',n,',a,f15.4,a)'
                end if
                if ( scalingpower.eq.0 ) then
                    write(unit,string) '#@ "',nint(i*delta),
     +                   '" ',log(i*delta),'\\'
                else
                    write(unit,string) '#@ "',nint(i*delta),
     +               '" ',(i*delta)**(scalingpower),'\\'
                end if
            else
                n = int(-0.01-log10(i*delta))
                n = max(1,min(n,7))
                write(string,'(a,i1,a,i1,a)') '(a,f',n+2,'.',n,
     +               ',a,f15.4,a)'
                if ( scalingpower.eq.0 ) then
                    write(unit,string) '#@ "',i*delta,
     +                   '" ',log(i*delta),'\\'
                else
                    write(unit,string) '#@ "',i*delta,
     +               '" ',(i*delta)**(scalingpower),'\\'
                end if
            end if
            if ( i.lt.i2 ) then
                write(unit,'(a)') '#@ ,\\'
            end if
        end do
        write(unit,'(a)') '#@ )'
        end subroutine
*  #] print_y2tics:
*  #[ plot_tx_cdfs:
        subroutine plot_tx_cdfs(txtx,nmc,nens,ntype,j1,j2)
!
!       make a file to plot CDFs
!
        integer nmc,nens,ntype,j1,j2
        real txtx(nmc,3)
        integer iens,j
        real logtxtx(3)
        integer iargc
        write(10,'(a)') '# fraction, return values in the past '//
     +       'climate, current climate and difference (sorted), '
        if ( ntype.eq.2 ) then
            write(10,'(a)') '# repeated as gumbel transforms'
        else if ( ntype.eq.3 ) then
            write(10,'(a)') '# repeated as logarithmic transforms'
        else if ( ntype.eq.4 ) then
            write(10,'(a)') '# repeated as sqrt-logarithmic transforms'
        else
            write(0,*) 'plot_tx_cdfs: error: unknown value for ntype ',
     +           ntype
            call abort
        end if
        do iens=1,nens
            do j=1,2
                if ( ntype.eq.2 ) then
                    if ( txtx(iens,j).lt.1e20 .and. txtx(iens,j).gt.1 )
     +                   then
                        logtxtx(j) = -log(1-1/((j2-j1+1)*txtx(iens,j)))
                        if ( logtxtx(j).gt.0 ) then
                            logtxtx(j) = -log(logtxtx(j))
                        else
                            logtxtx(j) = 3e33
                        end if
                    else
                        logtxtx(j) = 1e20
                    end if
                else if ( ntype.eq.3 ) then
                    if ( txtx(iens,j).lt.1e20 .and. txtx(iens,j).gt.0 )
     +                   then
                        logtxtx(j) = log(txtx(iens,j)*(j2-j1+1))
                    else
                        logtxtx(j) = 1e20
                    end if                        
                else if ( ntype.eq.4 ) then
                    if ( txtx(iens,j).lt.1e20 .and. txtx(iens,j).gt.1 )
     +                   then
                        logtxtx(j) = sqrt(log(txtx(iens,j)*(j2-j1+1)))
                    else
                        logtxtx(j) = 1e20
                    end if                        
                end if
            end do
            ! simple logscale for the time being
            logtxtx(3) = log(txtx(iens,3))
            write(10,'(f6.4,6g16.5)') iens/real(nens+1),
     +           (txtx(iens,j),j=1,3),(logtxtx(j),j=1,3)
        end do
        ! DO NOT write xtics in return times for easier reading of the plot
        !!!call print_xtics(10,-ntype,ntot,j1,j2)
        end
*  #] plot_tx_cdfs:
*  #[ getreturnlevels:
        subroutine getreturnlevels(a,b,xi,alpha,beta,cov1,cov2,
     +       covreturnlevel,j1,j2,t)
        implicit none
        integer j1,j2
        real a,b,xi,alpha,beta,cov1,cov2,t(10,3)
        real,external :: covreturnlevel
        integer i
        real x,xx
        do i=1,10
            if ( mod(i,3).eq.1 ) then
                x = 1 + i/3
            elseif ( mod(i,3).eq.2 ) then
                x = 1 + log10(2.) + i/3
            else
                x = log10(5.) + i/3
            endif
            x = x + log10(real(j2-j1+1))
            xx = x ! some routines modify x :-(
            t(i,1) = covreturnlevel(a,b,xi,alpha,beta,xx,cov1)
            xx = x
            t(i,2) = covreturnlevel(a,b,xi,alpha,beta,xx,cov2)
            t(i,3) = t(i,2) - t(i,1)
        enddo
        end
*  #] getreturnlevels:
*  #[ getreturnyears:
        subroutine getreturnyears(a,b,xi,alpha,beta,xyear,cov1,cov2,
     +       covreturnyear,j1,j2,tx,lwrite)
        implicit none
        integer j1,j2
        real a,b,xi,alpha,beta,xyear,cov1,cov2,tx(3)
        logical lwrite
        real,external :: covreturnyear
        tx(1) = covreturnyear(a,b,xi,alpha,beta,xyear,cov1)
        tx(2) = covreturnyear(a,b,xi,alpha,beta,xyear,cov2)
        if ( tx(2).gt.1e19 ) then
            if ( tx(1).gt.1e19 ) then
                tx(3) = 3e33
            else
                tx(3) = 1e20
            end if
        else
            if ( tx(1).gt.1e19 ) then
                tx(3) = 1e20
            else
                tx(3) = tx(1) / tx(2)
            end if
        end if
        if ( lwrite ) then
            print *,'return time = ',tx
        end if
        end
*  #] getreturnyears:
*  #[ getabfromcov:
        subroutine getabfromcov(a,b,alpha,beta,cov,aa,bb)
        implicit none
        real a,b,alpha,beta,cov,aa,bb
        character cassume*5
        common /fitdata4/ cassume
        integer ncur
        real restrain
        logical llwrite
        common /fitdata2/ restrain,ncur,llwrite
        
        if ( cassume.eq.'shift' ) then
            aa = a + alpha*cov
            bb = b
        else if ( cassume.eq.'scale' ) then
            aa = a*exp(alpha*cov/a)
            bb = b*aa/a
        else if ( cassume.eq.'both' ) then
            aa = a + alpha*cov
            bb = b + beta*cov
        else
            write(0,*) 'getabfromcov: error: unknown value for '//
     +           'assume: ',cassume
            call abort
        end if
        end
*  #] getabfromcov:
*  #[ printab:
        subroutine printab(lweb)
        implicit none
        logical lweb
        character cassume*5
        common /fitdata4/ cassume

        if ( lweb ) then
            if ( cassume.eq.'shift' ) then
                print '(a)','# <tr><td colspan="4">'//
     +               'with a''= a+&alpha;T '//
     +               'and b'' = b</td></tr>'
            else if ( cassume.eq.'scale' ) then
                print '(a)','# <tr><td colspan="4">'//
     +               'with a'' = a exp(&alpha;T/a) '//
     +               'and b'' = b exp(&alpha;T/a)</td></tr>'
            else if ( cassume.eq.'both' ) then
                print '(a)','# <tr><td colspan="4">'//
     +               'with a''= a+&alpha;T '//
     +               'and b'' = b+&beta;T</td></tr>'
            else
                write(0,*) 'printpoint: error: unknow value for assume '
     +               ,cassume
            end if
        else
            if ( cassume.eq.'shift' ) then
                print '(a)','a'' = a+alpha*T, b''= b'
            else if ( cassume.eq.'scale' ) then
                print '(a)','a'' = a*exp(alpha*T/a), b''= b*a''/a'
            else if ( cassume.eq.'both' ) then
                print '(a)','a'' = a+alpha*T, b''= b+beta*T'
            else
                write(0,*) 'fitgaucov: error: unknow value for assume ',
     +               cassume
            end if
        end if
        end
*  #] printab:
*  #[ adjustyy:
        subroutine adjustyy(ntot,xx,assume,a,b,alpha,beta,cov,
     +       yy,zz,aaa,bbb,lchangesign,lwrite)
!
!       input: xx,assume,a,b,alpha,beta,cov
!       output: yy,zz,aaa,bbb
!       flag: lwrite
!
        implicit none
        integer ntot
        real xx(2,ntot),yy(ntot),zz(ntot)
        real a,b,alpha,beta,cov,aaa,bbb
        character assume*(*)
        logical lchangesign,lwrite
        integer i
!
        do i=1,ntot
            yy(i) = xx(1,i)
            zz(i) = xx(2,i)
        end do
        if ( assume.eq.'shift' ) then
            if ( lchangesign ) then
                do i=1,ntot
                    yy(i) = yy(i) + alpha*(zz(i)-cov)
                end do
            else
                do i=1,ntot
                    yy(i) = yy(i) - alpha*(zz(i)-cov)
                end do
            end if
            aaa = a+cov*alpha
            bbb = b
        else if ( assume.eq.'scale' ) then
            if ( lchangesign ) then
                do i=1,ntot
                    yy(i) = yy(i)*exp(alpha*(zz(i)-cov)/a)
                end do
            else
                do i=1,ntot
                    yy(i) = yy(i)*exp(-alpha*(zz(i)-cov)/a)
                end do
            end if
            aaa = a*exp(alpha*cov/a)
            bbb = b*exp(alpha*cov/a)
        else if ( assume.eq.'both' ) then
            write(0,*) 'adjustyy: error: not yet ready'
            write(*,*) 'adjustyy: error: not yet ready'
            call abort
        else
            write(0,*) 'adjustyy: error: unknown assumption ',assume
            write(*,*) 'adjustyy: error: unknown assumption ',assume
            call abort            
        end if

        end
*  #] adjustyy:
*  #[ write_obscov:
        subroutine write_obscov(xx,yrs,ntot,xmin,cov2,xyear,year,offset,
     +       lchangesign)
        implicit none
        integer ntot,yrs(0:ntot),year
        real xx(2,ntot),xmin,cov2,xyear,offset
        logical lchangesign
        integer i,is
        logical lopen
        character string*1000,arg*250
        integer iargc
        inquire(unit=15,opened=lopen)
        if ( lopen ) then
            if ( lchangesign ) then
                is = -1
            else
                is = +1
            end if
            string = ' '
            do i=0,iargc()
                call getarg(i,arg)
                string = trim(string) // ' '
     +               // arg(index(arg,'/',.true.)+1:)
            end do
            write(15,'(a)') '# ' // trim(string)
            write(15,'(a)') '# covariate  value'
            do i=1,ntot
                if ( xx(1,i).gt.xmin ) then
                    write(15,'(2g20.6,i11)') xx(2,i)+offset,is*xx(1,i),
     +                   yrs(i)
                end if
            end do
            write(15,'(a)')
            write(15,'(a)')
            if ( xyear.lt.1e33 ) then
                write(15,'(2g20.6,i11)') cov2+offset,is*xyear,year
            else
                write(15,'(g20.6,a,i11)') cov2+offset,'-999.900',0
            end if
        end if
        end subroutine
*  #] write_obscov:
*  #[ write_threshold
        subroutine write_threshold(cmin,cmax,a,b,alpha,beta,offset,
     +       lchangesign)
        implicit none
        real cmin,cmax,a,b,alpha,beta,offset
        logical lchangesign
        integer i,is
        real c,aa,bb
        logical lopen
        inquire(unit=15,opened=lopen) ! no flag in getopts.inc
        if ( lopen ) then
            if ( lchangesign ) then
                is = -1
            else
                is = +1
            end if
            write(15,'(a)')
            write(15,'(a)')
            write(15,'(a)') '# covariate threshold/position scale'
            do i=0,100
                c = cmin + (cmax-cmin)*i/100.
                call getabfromcov(a,b,alpha,beta,c,aa,bb)
                write(15,'(3g20.6)') c+offset,is*aa,is*bb
            end do
        endif
        end subroutine
*  #] write_threshold:
*  #[ write_dthreshold
        subroutine write_dthreshold(cov1,cov2,acov,offset,lchangesign)
        implicit none
        real cov1,cov2,acov(3,2),offset
        logical lchangesign
        integer i,j,is
        logical lopen
        inquire(unit=15,opened=lopen) ! no flag in getopts.inc
        if ( lopen ) then
            is = +1 ! pre-flipped
            write(15,'(a)')
            write(15,'(a)')
            write(15,'(a)')
     +           '# covariate threshold lowerbound upperbound'
            write(15,'(4g20.6)') cov1+offset,(is*acov(j,1),j=1,3)
            write(15,'(4g20.6)') cov2+offset,(is*acov(j,2),j=1,3)
        endif
        end subroutine
*  #] write_threshold:
