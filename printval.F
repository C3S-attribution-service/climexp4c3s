*  #[ printval:
      subroutine printval(string,val,plot,iuplot,lweb,lchangesign)
      implicit none
      integer iuplot
      real val
      logical plot,lweb,lchangesign
      character string*(*)
      real x
      x = val
      if ( lchangesign ) x = -x
      if ( lweb ) then
          print '(3a,g16.4,10a)','# <tr><td>',string(3:),'</td><td>',x
     +         ,'</td><td>&nbsp;</td></tr>'
      else
          print '(a,6g16.4)',string,x
      endif
      if ( plot ) write(iuplot,'(6g16.4)') x
      end
*  #] printval:
*  #[ printvalerr:
      subroutine printvalerr(string,array,mboot,plot,iuplot,lweb
     +       ,lchangesign)
      implicit none
      integer iuplot,mboot
      real array(0:mboot)
      logical plot,lweb,lchangesign
      character string*(*)
      real x,x025,x160,x500,x840,x975
      if ( mboot.eq.0 ) then
          call printval(string,array(0),plot,iuplot,lweb,lchangesign)
      else
          x = array(0)
          call getcut1(x025, 2.5,mboot,array(1),.false.)
          call getcut1(x160,16.0,mboot,array(1),.false.)
          call getcut1(x500,50.0,mboot,array(1),.false.)
          call getcut1(x840,84.0,mboot,array(1),.false.)
          call getcut1(x975,97.5,mboot,array(1),.false.)
          if ( lchangesign ) then
              x = -x
              call swapminus(x025,x975)
              call swapminus(x160,x840)
              x500 = -x500
          endif
          if ( lweb ) then
              print '(3a,7(g16.6,a))','# <tr><td>',string(3:),
     +             '</td><td>',x,'&plusmn;',(x975-x025)/2
     +             ,'</td><td>',x025,'...',x975,'</td></tr>'
          else
              print '(a,7g16.6)',string,x,x025,x160,x500,x840,x975,
     +             (x975-x025)/2
          endif
          if ( plot ) write(iuplot,'(7g16.6)')
     +         x,x025,x160,x500,x840,x975
      endif
      end
*  #] printvalerr:
*  #[ printuntransf:
        subroutine printuntransf(xin)
        implicit none
#include "getopts.inc"
        real xin
        real x
        x = xin
        if ( logscale .or. sqrtscale .or. squarescale .or. twothirdscale
     +       ) then
            if ( logscale ) x = 10.**(x)
            if ( sqrtscale ) x = x**2
            if ( squarescale ) x = sqrt(x)
            if ( twothirdscale ) x = x**(3./2.)
            if ( lweb ) then
                print '(a,f16.2)','# w/o transf.',x
            else
                print '(a,f16.2,a)','# <tr><td>w/o transf.</td><td>',x
     +               ,'</td><td>&nbsp;</td></tr>'
            endif
        endif
        end
*  #] printuntransf:
*  #[ print3untransf:
        subroutine print3untransf(x1in,x2in,x3in,year)
        implicit none
#include "getopts.inc"
        integer year
        real x1in,x2in,x3in
        real x(3)
        integer i
        if ( logscale .or. sqrtscale .or. squarescale .or.
     +       cubescale .or. twothirdscale ) then
            do i=1,3
                if ( i.eq.1 ) then
                    x(i) = x1in
                elseif ( i.eq.2 ) then
                    x(i) = x2in
                else
                    x(i) = x3in
                endif
                if ( logscale ) x(i) = 10.**(x(i))
                if ( sqrtscale ) x(i) = x(i)**2
                if ( cubescale ) x(i) = x(i)**3
                if ( squarescale ) x(i) = sqrt(x(i))
                if ( twothirdscale ) x(i) = x(i)**(3./2.)
            enddo
            if ( .not.lweb ) then
                print '(a,f16.3,a,f16.3,a,f16.3,a)','# w/o transf.',x(1)
     +               ,'(',x(2),'...',x(3),')'
            else if ( year.eq.0 ) then
                print '(a,g16.5,a,g16.5,a,g16.5,a)',
     +               '# <tr><td>w/o transformation</td><td>',x(1),
     +               '</td><td>',x(2),'...',x(3),'</td></tr>'
            else
                print '(a,i4,a,g16.5,a,g16.5,a,g16.5,a)',
     +               '# <tr><td>w/o transformation</td><td>',year,
     +               '</td><td>',x(1),
     +               '</td><td>',x(2),'...',x(3),'</td></tr>'
            end if
        endif
        end
*  #] print3untransf:
*  #[ swapminus:
        subroutine swapminus(a,b)
        real a,b
        real c
        c = a
        a = -b
        b = -c
        end
*  #] swapminus:
