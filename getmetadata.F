        subroutine getmetadata(infile,mens1,mens,ncid,datfile,nxmax,nx
     +       ,xx,nymax,ny,yy,nzmax,nz,zz,lz,nt,nperyear,firstyr,firstmo
     +       ,ltime,undef,endian,title,history,nvarmax,nvars,vars,ivars
     +       ,lvars,svars,units,cell_methods,lwrite)
!
!       get metadata of the field in file infile, either grads or netcdf
!
        implicit none
#include "netcdf.inc"
        integer mens1,mens,ncid,nxmax,nx,nymax,ny,nzmax,nz,nt,nperyear
     +       ,firstyr,firstmo,endian,nvarmax,nvars,ivars(6,nvarmax)
        real xx(nxmax),yy(nymax),zz(nzmax),undef
        character infile*(*),datfile*(*),lz(3)*(*),ltime*(*),title*(*)
     +       ,history*(*),vars(nvarmax)*(*),lvars(nvarmax)*(*)
     +       ,svars(nvarmax)*(*),units(nvarmax)*(*),
     +       cell_methods(nvarmax)*(*)
        logical lwrite
        integer nensmax,ntmax
        parameter(nensmax=120,ntmax=60000)
        integer status,i,nens1,nens2,it
        character file*255
        logical ensemble,lexist,tdefined(60000)

        if ( lwrite ) then
            print *,'getmetadata: infile = ',trim(infile)
        endif
        file = infile
        if ( index(file,'%').gt.0 .or. index(file,'++').gt.0 ) then
            ensemble = .true.
            call filloutens(file,0)
            inquire(file=file,exist=lexist)
            if (.not.lexist ) then
                mens1 = 1
                file = infile
                call filloutens(file,1)
            else
                mens1 = 0
            endif
        else
            ensemble = .false.
            mens1 = 0
            mens = 0
        endif
        if ( lwrite ) then
            print *,'mens1 = ',mens1
            print *,'getmetadata: nf_opening file ',trim(file)
        endif
        status = nf_open(file,nf_nowrite,ncid)
        if ( status.ne.nf_noerr ) then
            call parsectl(file,datfile,nxmax,nx,xx,nymax,ny,yy,
     +           nzmax,nz,zz,nt,nperyear,firstyr,firstmo,undef,endian,
     +           title,nvarmax,nvars,vars,ivars,lvars,units)
            nz = max(1,ivars(1,1))
            ncid = -1
            if ( ensemble ) then
                do mens=1,nensmax
                    file = infile
                    call filloutens(file,mens)
                    inquire(file=file,exist=lexist)
                    if ( .not.lexist ) goto 100
                enddo
  100           continue
                mens = mens - 1
            endif
            lz = ' '
            ltime  = ' '
            history = ' '
            svars(1:nvars) = ' '
            cell_methods = ' '
        else
            datfile = file
            call ensparsenc(file,ncid,nxmax,nx,xx,nymax,ny,yy,nzmax,
     +           nz,zz,lz,nt,nperyear,firstyr,firstmo,ltime,tdefined
     +           ,ntmax,nens1,nens2,undef,title,history,nvarmax,nvars
     +           ,vars,ivars,lvars,svars,units,cell_methods)
            do it=1,nt
                if ( .not.tdefined(it) ) then
                    write(0,*) 'getmetadat: error: cannot handle ',
     +                   'holes in time axis yet'
                    call abort
                end if
            end do
            if ( .not.ensemble ) then
                mens1 = nens1
                mens  = nens2
            endif
            if ( ensemble ) then
                do mens=1,nensmax
                    file = infile
                    call filloutens(file,mens)
                    status = nf_open(file,nf_nowrite,i)
                    if ( status.ne.nf_noerr ) goto 200
                    status = nf_close(i)
                enddo
 200            continue
                mens = mens - 1
            endif
        endif
        if ( ensemble ) then
            i = 1 + index(infile,'/',.true.)
            write(0,*) 'located ',mens-mens1+1
     +               ,' ensemble members of ',trim(infile(i:)),'<br>'
        endif
        end
