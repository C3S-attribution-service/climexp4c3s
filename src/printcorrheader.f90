subroutine printcorrheader

!   prints the header for the correlation tables

    implicit none
    include 'getopts.inc'
    integer :: ndash,i
    character :: stationstring*255,width*20,var*4

    if ( lrmse ) then
        var = 'RMSE'
    else if ( lmae ) then
        var = 'MAE '
    else
        var = 'corr'
    end if
    if ( lweb ) then
        print '(a)','<table class="realtable" width=452 border=0'// &
            ' cellpadding=0 cellspacing=0>'
    else
        if ( lbootstrap ) then
            ndash=103
        else
            ndash=68
        endif
        print '(1000a)',('=',i=1,ndash)
    endif
    if ( lks ) then
        print '(a)', &
            ' index    cutoff  t-test KS             below'// &
            '                above'
        print '(a)', &
            '  name            sign.  sign.    #    mean    s.d.'// &
            '    #    mean    s.d.'
    elseif ( lconting ) then
        print '(a)', &
            ' index      cutoff       below      normal'// &
            '       above         sum  sign. lag'
    elseif ( nfittime > 0 ) then
        print '(3a)', &
            ' index    ',var,'   sign.   no       relaxation'// &
            '        relation     corr lag'
        print '(a)','  name              %  pnts '// &
            '    regr.   error     regr.   error'
    elseif ( lbootstrap ) then
        if ( lweb ) then
            if ( namestring == ' ' ) then
                stationstring = '&nbsp;'
                width = ' '
            else
                stationstring = namestring
                width = 'width=220'
            endif
            print '(10a)','<tr><th ',trim(width),'>', &
                trim(stationstring),'</th>'// &
                '<th align="right">months</th>'// &
                '<th align="right">lag</th>'// &
                '<th align="right">'//trim(var)//'</th>'// &
                '<th align="right">p</th>'// &
                '<th align="right">no</th>'// &
                '<th align="right">95% CI</th></tr>'
        else
            print '(3a)', &
                ' index     ',var,'   sign.  no          data'// &
                '             index      lag          bootstrap'
            print '(a)', &
                '  name              % pnts      mean    s.d.      mean    s.d.' &
                //'         2.5%    16%    50%    84%  97.5%'
        endif
    else
        if ( lweb ) then
            if ( namestring == ' ' ) then
                stationstring = '&nbsp;'
            else
                stationstring = namestring
            endif
            print '(4a)','<tr><th width=220>',trim(stationstring), &
                '</th>', &
                '<th align="right">months</th>'// &
                '<th align="right">lag</th>'// &
                '<th align="right">'//trim(var)//'</th>'// &
                '<th align="right">p</th>'// &
                '<th align="right">no</th>'
        else
            print '(3a)', &
                ' index     ',var,'   sign.  no          data'// &
                '             index      lag'
            print '(a)', &
                '  name              % pnts      mean    s.d.      mean    s.d.'
        endif
    endif
    if ( .not. lweb ) then
        print '(1000a)',('=',i=1,ndash)
    endif
end subroutine printcorrheader

subroutine printcorrfooter

!   prints the footer for the correlation tables

    implicit none
    include 'getopts.inc'
    integer :: ndash,i
    if ( lweb ) then
        print '(a)','</table>'
    else
        if ( lbootstrap ) then
            ndash=103
        else
            ndash=68
        endif
        print '(1000a)',('=',i=1,ndash)
    endif
end subroutine printcorrfooter
