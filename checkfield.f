        subroutine checkfield(field,nxf,nyf,npermax,firstyr,lastyr,
     +       nx,ny,nperyear,yr1,yr2,nens1,nens2,lnomissing)
        implicit none
        integer nxf,nyf,npermax,firstyr,lastyr,nx,ny,nperyear,yr1,yr2
     +       ,nens1,nens2
        real field(nxf,nyf,npermax,firstyr:lastyr,0:nens2)
        logical lnomissing
        integer i,j,mo,yr,iens
*
        if ( .not.lnomissing ) return
        do iens=nens1,nens2
            do yr=yr1,yr2-1
                do mo=1,nperyear
                    do j=1,ny
                        do i=1,nx
                            if ( field(i,j,mo,yr,iens).gt.1e20 ) then
                                write(0,*)
     +                               'checkfield: error: lnomissing'//
     +                               ' is true but field(',i,j,mo,yr
     +                               ,iens,') = ',field(i,j,mo,yr,iens)
                                call abort
                            endif
                        enddo
                    enddo
                enddo
            enddo
        enddo
        end
