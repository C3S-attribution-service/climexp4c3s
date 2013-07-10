        program untransform
!
!       remove the transformation in a gumbel or log plot
!
        implicit none
        integer i,n,itrans,yyyymm,nwords
        real xx(4)
        character line*256,flag*10
        integer iargc,getnumwords
!
        if ( iargc().ne.3 ) then
            write(0,*)
     +           'usage: untransform true|false true|false true|false'
            write(0,*) '       transforms a log|sqtr|square back'
            stop
        endif
!
        itrans = 0
        do i=1,3
            call getarg(i,flag)
            if ( flag.ne.'false' .and. flag.ne.'off' ) then
                if ( itrans.ne.0 ) then
                    write(0,*) 'untransform: error: can only handle'//
     +                   ' one transformation at a a time'
                    call abort
                endif
                itrans = i
            endif
        enddo
 100    continue
        read(*,'(a)',end=800) line
        if ( line(1:1).eq.'<' ) then
            print '(2a)','# ',trim(line)
            goto 100
        endif
        if ( line(1:1).eq.'#' .or. line.eq.' ' ) then
            print '(a)',trim(line)
            goto 100
        endif
        if ( index(line,'amoeba').ne.0 ) goto 100
        nwords = getnumwords(line)
        if ( nwords.eq.5 ) then
            read(line,*) n,xx
        elseif ( nwords.eq.6 ) then
            read(line,*) n,xx,yyyymm
        else
            write(0,'(2a)') 'untransform: funny line: ',trim(line)
            write(*,'(2a)') '# untransform: funny line: ',trim(line)
        endif
        if ( itrans.eq.1 ) then
            if ( xx(2).ne.-999.9 ) xx(2) = 10.**xx(2)
            if ( xx(3).ne.-999.9 ) xx(3) = 10.**xx(3)
        elseif ( itrans.eq.2 ) then
            if ( xx(2).ne.-999.9 ) xx(2) = xx(2)**2
            if ( xx(3).ne.-999.9 ) xx(3) = xx(3)**2
        elseif ( itrans.eq.3 ) then
            if ( xx(2).ne.-999.9 ) then
                if ( xx(2).ge.0 ) then
                    xx(2) = sqrt(xx(2))
                else
                    write(0,*) 'untransform: error: neagtive number: '
     +                   ,xx(2)
                endif
            endif
            if ( xx(3).ne.-999.9 ) then
                if ( xx(3).ge.0 ) then
                    xx(3) = sqrt(xx(3))
                else
                    write(0,*) 'untransform: error: neagtive number: '
     +                   ,xx(3)
                endif
            endif
        endif
        if ( nwords.eq.5 ) then
            write(*,'(i8,4g22.6,i9)') i,xx
        elseif ( nwords.eq.6 ) then
            write(*,'(i8,4g22.6,i9)') i,xx,yyyymm
        endif
        goto 100
 800    continue
        end
