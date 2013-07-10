        program plotdat
*
*       get nperyear
*
        implicit none
	include 'param.inc'
        integer nperyear
        real data(npermax,yrbeg:yrend)
        character line*255,var*40,units*20
        logical lwrite
	integer iargc
*
*       init
*
        lwrite = .false.
        if ( iargc().eq.0 ) then
            print *,'usage: getnperyear datafile'
            stop
        endif
*
        call getarg(1,line)
        call readseries(line,data(1,yrbeg),npermax,yrbeg,yrend,nperyear
     +       ,var,units,.false.,lwrite)
        if ( nperyear.lt.10 ) then
            print '(i1)',nperyear
        elseif ( nperyear.lt.100 ) then
            print '(i2)',nperyear
        elseif ( nperyear.lt.1000 ) then
            print '(i3)',nperyear
        elseif ( nperyear.lt.10000 ) then
            print '(i4)',nperyear
        endif
        end

