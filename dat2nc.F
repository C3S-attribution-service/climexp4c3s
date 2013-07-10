        program dat2nc
!
!       program to convert a .dat-style time series into a netCDF file 
!       using Camiel's classes
!
        implicit none
#include "param.inc"
        integer             :: yr,mo,nperyear,i,j,n
        real*4              :: data(npermax,yrbeg:yrend),x,y,x0,y0
        character           :: file*256,since*32,var*40,lvar*80,
     +       description*128
        character           :: type*1,comment*10000,line*80,units*20
        logical             :: lwrite,lstandardunits
        integer             :: iargc
        integer,external    :: llen,leap
        lwrite = .FALSE.
*
*        read data - easy f77
*
        if ( iargc().ne.4 ) then
            print *,'usage: dat2nc infile type name outfile ',iargc()
            stop
        endif
        call getarg(1,file)
        call makeabsent(data,npermax,yrbeg,yrend)
        if ( lwrite ) print *,'reading from ',trim(file)
        call readseries(file,data,npermax,yrbeg,yrend,nperyear,var,
     +           units,lstandardunits,lwrite)
        open(1,file=file)
        comment = ' '
        lvar = ' '
        do i=1,100
            read(1,'(a)') description
            if ( description(1:1).ne.'#' ) exit
            j = index(description,']')
            if ( j.eq.0 ) then
                comment(2+len_trim(comment):) = trim(description(3:))
            else
                lvar = description(j+2:)
            end if
        end do
        close(1)
        comment(2+llen(comment):) = 'via the KNMI Climate Explorer'
        comment(2+llen(comment):) = '(http://climexp.knmi.nl)'
        if ( var.eq.' ' ) then
            call getarg(2,file)
            if ( file.eq.'i') then
                var = 'index'
            elseif ( file.eq.'t' ) then
                var = 'temperature'
            elseif ( file.eq.'p' ) then
                var = 'precipitation'
            elseif ( file.eq.'s' ) then
                var = 'pressure'
            elseif ( file.eq.'l' ) then
                var = 'sealevel'
            else
                var = 'unknown'
            endif
        endif
        call getarg(3,description)
        do i=1,len(description)
            if ( description(i:i).eq.'_') description(i:i) = ' '
        enddo
        call getarg(4,file)
!
!       and write out
!
        call writencseries(file,data,npermax,yrbeg,yrend,nperyear
     +       ,description,comment,var,lvar,units)
        end
