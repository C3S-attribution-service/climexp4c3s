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
     +                   ,'...',t975(i),'</td></tr>'
                    call print3untransf(t(i),t25(i),t975(i))
                enddo
            else
                do i=1,4
                    print '(a,i5,a,f16.3,a,2f16.3)'
     +                   ,'# value for return period ',10**i,' year: '
     +                   ,t(i),' 95% CI ',t25(i),t975(i)
                    call print3untransf(t(i),t25(i),t975(i))
                enddo
            endif
        endif
        end
*  #] printreturnvalue:
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
                print '(a,i5,a,g16.5,a,g16.5,a,g16.5,a)'
     +               ,'# <tr><td>return period ',year,'</td><td>',tx
     +               ,'</td><td>',tx25,'...',tx975,'</td></tr>'
            else
                print '(a,i4,a,f16.5,a,2f16.5)','# return time ',year
     +               ,' = ',tx,' 95% CI ',tx25,tx975
            endif
        endif
        end
*  #] printreturntime:
*  #[ plotreturnvalue:
        subroutine plotreturnvalue(ntype,t25,t975,n)
*
*       print return value at a few times
*
        implicit none
        integer ntype,n
        real t25(10),t975(10)
        logical lweb
        integer i
        real x,f
*
        if ( ntype.eq.2 .or. ntype.eq.3 .or. ntype.eq.4 ) then ! extreme value  plot
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
