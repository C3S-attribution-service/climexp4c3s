        subroutine checklsmask(lsmask,nx,ny,lwrite)
*
*       normalize the land/sea mask to out convention 0=sea, 1=land
*       based on too much experience
*
        implicit none
        integer nx,ny
        real lsmask(nx,ny)
        logical lwrite
        integer i,j
        real w
*
        w = 0
        do j=1,ny
            do i=1,nx
                if ( abs(lsmask(i,j)).lt.0.0001 ) lsmask(i,j) = 0
                if ( abs(lsmask(i,j)).gt.1e10 ) lsmask(i,j) = 0
                w = max(w,lsmask(i,j))
            enddo
        enddo
        if ( abs(w-100).lt.1e-2 ) then
            do j=1,ny
                do i=1,nx
                    lsmask(i,j) = lsmask(i,j)/100
                enddo
            enddo
            w = w/100
        endif
        if ( abs(w-1).gt.1e-4 ) then
            write(0,*) 'checklsmask: error: max LS mask = ',w
            write(*,*) 'checklsmask: error: max LS mask = ',w
            call abort
        endif
        if ( lwrite ) then
            print *,'LS mask'
            do j=ny,1,-1
                print '(i4,1x,1000i1)',j,
     +               (nint(lsmask(i,j)),i=1,nx)
            enddo
        endif
        end
