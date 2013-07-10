*  #[ filloutens:
        subroutine filloutens(file,iens)
        implicit none
        character*(*) file
        integer iens
        integer i,j
        character formaat*4,teken*1
        integer llen
        external llen
*
        i = index(file,'%')
        if ( i.eq.0 ) then
            i = index(file,'++')
            if ( i.eq.0 ) then
                write(0,*)
     +                'filloutens: error: cannot find + or % in file '
     +                ,trim(file)
                call abort
            endif
            teken = '+'
        else
            teken = '%'
        endif
        j = i
   10   continue
        j = j + 1
        if ( j.le.len(file) ) then
                if ( file(j:j).eq.teken ) goto 10
        endif
        j = j - 1
***        write(formaat,'(a,i1,a,i1)') 'i',j-i+1,'.',j-i+1
***        write(file(i:j),formaat) iens
	if ( j.eq.i ) then
            write(file(i:j),'(i1.1)') iens
        elseif ( j.eq.i+1 ) then
            write(file(i:j),'(i2.2)') iens
        elseif ( j.eq.i+2 ) then
            write(file(i:j),'(i3.3)') iens
        elseif ( j.eq.i+3 ) then
            write(file(i:j),'(i4.4)') iens
        else
	    write(0,*) 'filloutens: iens too big ',iens,10**(j-i+1)-1
	    call abort
        endif
        end
*  #] filloutens:
*  #[ filloutleadens:
        subroutine filloutleadens(file,lead,iens)
        implicit none
        character*(*) file
        integer iens,lead
        integer i,j
        character formaat*4,teken*1
        integer llen
        external llen
*
        i = index(file,'%%')
        if ( i.eq.0 ) then
            i = index(file,'++')
            if ( i.eq.0 ) then
                write(0,*)
     +                'filloutleadens: error: '//
     +                'cannot find ++ or %% in file ',file(1:llen(file))
                call abort
            endif
            teken = '+'
        else
            teken = '%'
        endif
        j = i
   10   continue
        j = j + 1
        if ( j.le.len(file) ) then
                if ( file(j:j).eq.teken ) goto 10
        endif
        j = j - 1
***        write(formaat,'(a,i1,a,i1)') 'i',j-i+1,'.',j-i+1
***        write(file(i:j),formaat) iens
	if ( iens.lt.10 .and. j.eq.i ) then
            write(file(i:j),'(i1.1)') iens
        elseif ( iens.lt.100 .and. j.eq.i+1 ) then
            write(file(i:j),'(i2.2)') iens
        elseif ( iens.lt.1000 .and. j.eq.i+2 ) then
            write(file(i:j),'(i3.3)') iens
        elseif ( iens.lt.10000 .and. j.eq.i+3 ) then
            write(file(i:j),'(i4.4)') iens
        else
	    write(0,*) 'filloutleadens: iens too large ',iens,10**(j-i)
	    call abort
        endif
*       
        j = i-1
        i = index(file(1:j),'%')
        if ( i.eq.0 ) then
            i = index(file(1:j),'+_')
            if ( i.eq.0 ) then
                goto 900
            endif
            teken = '+'
        else
            teken = '%'
        endif
        j = i
   20   continue
        j = j + 1
        if ( j.le.len(file) ) then
                if ( file(j:j).eq.teken ) goto 20
        endif
        j = j - 1
***        write(formaat,'(a,i1,a,i1)') 'i',j-i+1,'.',j-i+1
***        write(file(i:j),formaat) lead
	if ( lead.lt.10 .and. j.eq.i ) then
            write(file(i:j),'(i1)') lead
        elseif ( lead.lt.10 .and. j.eq.i+1 ) then
            write(file(i:j),'(2i1)') 0,lead
        elseif ( lead.lt.100 .and. j.eq.i+1 ) then
            write(file(i:j),'(i2)') lead
        else
	    write(0,*) 'filloutleadens: cannot handle lead>99 yet'
	    call abort
        endif
  900   continue
        end
*  #] filloutleadens:
