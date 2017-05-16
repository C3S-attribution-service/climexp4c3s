        subroutine readfield(ncid,infile,datfile,field,nxf,nyf
     +       ,nzf,npermax,fyr,lyr,nens1,nens2,nx,ny,nz,nperyear,yr1,yr2
     +       ,firstyr,firstmo,nt,undef,endian,vars,units,lstandardunits
     +       ,lwrite)
*
*       read field, change absent values to our convention
*
        implicit none
#include "netcdf.inc"
        integer ncid,nxf,nyf,nzf,npermax,fyr,lyr,nx,ny,nz,nperyear
     +       ,yr1,yr2,nens1,nens2,firstyr,firstmo,nt,endian
        real field(nxf,nyf,nzf,npermax,fyr:lyr,0:nens2)
        real undef
        logical lstandardunits,lwrite
        character infile*(*),datfile*(*),vars(1)*(*),units(1)*(*)
        integer i,j,iens,fiyr,fimo,ix,iy,iz,it,nvars,jvars(6,1),status
        real xx(nxf),yy(nyf),zz(nzf)
        character file*255,line*255,tmpvars(1)*60,lvars(1)*100,
     +       tmpunits(1)*80,title*512

        if ( lwrite ) then
            print *,'readfield: field has dims ',nxf,nyf,nzf,npermax,fyr
     +           ,lyr,nens2
            print *,'           ncid    = ',ncid
            print *,'           nens1,2 = ',nens1,nens2
        end if
        if ( firstyr.gt.10000 ) then
            write(0,*) 'readfield: error while reading ',trim(infile)
            write(0,*) 'readfield: error: wrong value for firstyr '
     +           ,firstyr
            call abort
        end if
        call keepalive(0,0)
        i = index(infile,'%') + index(infile,'++')
        if ( i.gt.0 ) then
            write(0,*)'Using ensemble members ',nens1,' to ',nens2,'<p>'
        end if
        if ( ncid.eq.-1 ) then
            do iens=nens1,nens2
                if ( iens.gt.nens1 ) then
                    file = infile
                    call filloutens(file,iens)
                    open(1,file=file,status='old')
                    read(1,'(a)') line
                    j = index(line,'^')
                    if ( j.ne.0 ) then
                        do i=len(infile),1,-1
                            if ( infile(i:i).eq.'/' ) goto 110
                        enddo
 110                    continue
                        datfile = file(1:i)//line(j+1:)
                    else
                        datfile = line(6:)
                    endif
                endif
                call zreaddatfile(datfile,field(1,1,1,1,fyr,iens),
     +               nxf,nyf,nzf,nx,ny,nz,nperyear,fyr,lyr,
     +               firstyr,firstmo,nt,undef,endian,lwrite,
     +               yr1,yr2,1,1)
                call keepalive(iens,nens2-nens1+1)
            end do
        else
            do iens=nens1,nens2
                file = infile
                if ( i.gt.0 ) call filloutens(file,iens)
                call parsenc(file,ncid,nxf,ix,xx,nyf,iy,yy
     +               ,nzf,iz,zz,it,nperyear,fiyr,fimo,undef
     +               ,title,1,nvars,tmpvars,jvars,lvars,tmpunits)
                if ( ix.ne.nx .or. iy.ne.ny .or. iz.ne.nz ) then
                    write(0,*) 'readfield: error: unequal grids'
                    write(0,*) ix,iy,iz,nx,ny,nz
                    call abort
                end if
                if ( fiyr.gt.10000 ) then
                    write(0,*) 'readfield: error while reading '
     +                   ,trim(infile)
                    write(0,*) 'readfield: error: wrong value for fiyr '
     +                   ,fiyr
                    call abort
                end if
                call zreadncfile(ncid,field(1,1,1,1,fyr,iens),
     +               nxf,nyf,nzf,ix,iy,iz,nperyear,fyr,lyr,
     +               fiyr,fimo,it,undef,lwrite,yr1,yr2,jvars)
                call keepalive(iens,nens2-nens1+1)
            enddo
        endif
        if ( lstandardunits ) then
*           convert to standard units
            do iens=nens1,nens2
                tmpunits(1) = units(1) ! they are otherwise adjusted
                call makestandardfield(field(1,1,1,1,fyr,iens),nx,ny
     +               ,nz,nperyear,fyr,lyr,nx,ny,1,nperyear
     +               ,max(fyr,yr1),min(lyr,yr2),vars(1)
     +               ,tmpunits,lwrite)
            enddo
            units(1) = tmpunits(1)
        end if
        end
