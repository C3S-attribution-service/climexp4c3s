        program seriesanomal
!
!       compute the anomaly series
!
        implicit none
#include "param.inc"
        integer yr1,yr2,nperyear,mens1,mens,iens
        real,allocatable :: data(:,:,:)
        logical lwrite,lstandardunits
        character file*255,outfile*255,ensfile*255,var*20,units*40,
     +       line*255
        integer iargc
        lwrite = .false.
        lstandardunits = .false.

        allocate(data(npermax,yrbeg:yrend,0:nensmax))
        if ( iargc().lt.2 ) then
            print *,'usage: seriesanomal [yr1 yr2] infile outfile'
            call abort
        endif

        if ( iargc().eq.4 ) then
            call getarg(1,file)
            read(file,*) yr1
            call getarg(2,file)
            read(file,*) yr2
            call getarg(3,file)
            call getarg(4,outfile)
        else
            yr1 = yrbeg
            yr2 = yrend
            call getarg(1,file)
            call getarg(2,outfile)
        endif
        call readensseries(file,data,npermax,yrbeg,yrend,nensmax
     +        ,nperyear,mens1,mens,var,units,lstandardunits,lwrite)
        call ensanomal(data,npermax,nperyear,yrbeg,yrend,mens1
     +       ,mens,yr1,yr2)
        do iens=mens1,mens
            if ( mens1.eq.mens ) then
                open(1,file=file)
                open(2,file=outfile)
            else
                ensfile=file
                call filloutens(ensfile,iens)
                open(1,file=ensfile)
                ensfile=outfile
                call filloutens(ensfile,iens)
                open(2,file=ensfile)
            endif
            do
                read(1,'(a)') line
                if ( line(1:1).ne.'#' ) exit
                if ( line(3:).ne.' ' ) write(2,'(a)') trim(line)
            enddo
            close(1)
            write(2,'(a,i4,a,i4,a,i2,a,i2)')
     +           '# anomalies with respect to ',yr1,'-',yr2,
     +           ' ensemble members ',mens1,'-',mens
            call printdatfile(2,data(1,yrbeg,iens),npermax,nperyear
     +           ,yrbeg,yrend)
            close(2)
        enddo
        end
