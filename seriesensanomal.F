        program seriesensanomal
!
!       compute the anomaly series
!
        implicit none
#include "param.inc"
        integer nperyear,mens1,mens,iens,nens1
        real data(npermax,yrbeg:yrend,0:nensmax)
        logical lwrite,lstandardunits,lexist
        character file*255,ensfile*255,var*40,units*40,line*255
        integer iargc
        lwrite = .false.
        lstandardunits = .false.

        if ( iargc().ne.3 ) then
            print *,'usage: seriesensanomal ensmember ensfile dummy'
            call abort
        endif

        call getarg(1,ensfile)
        read(ensfile,*) nens1
        call getarg(2,ensfile)        

        call readensseries(ensfile,data,npermax,yrbeg,yrend,nensmax
     +        ,nperyear,mens1,mens,var,units,lstandardunits,lwrite)
        call anomalensemble(data,npermax,nperyear,yrbeg,yrend,
     +           yrbeg,yrend,mens1,mens)
        call filloutens(ensfile,nens1)
        open(1,file=ensfile,status='old')
        do
            read(1,'(a)') line
            if ( line(1:1).ne.'#' ) exit
            if ( line(3:).ne.' ' ) write(*,'(a)') trim(line)
        enddo
        close(1)
        write(*,'(a)') '# anomalies with respect to the ensemble mean'
        call printdatfile(6,data(1,yrbeg,nens1),npermax,nperyear
     +       ,yrbeg,yrend)
        end
