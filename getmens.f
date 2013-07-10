        subroutine getmens(file,mens,nensmax)
        implicit none
        integer mens,nensmax
        character*(*) file
        character*255 ensfile
        logical lexist
*
        mens = -1
        if ( index(file,'%').ne.0 .or. index(file,'++').ne.0 ) then
            do mens=0,nensmax
                ensfile = file
                call filloutens(ensfile,mens)
                inquire(file=endfile,exist=lexist)
                if ( .not.lexist ) goto 100
            enddo
  100       continue
            mens = mens - 1
            call filloutens(file,0)
        endif
        end
