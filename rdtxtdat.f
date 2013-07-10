        subroutine rdtxtdat(file,data,yrbeg,yrend)
*
*       reads the data in file 'file' into the array 'data'
*       the data is assumed to be in either one of two formats:
*       - if the name contains the string 'txt' then it is a 
*         NCDC formatted file with one month per line
*       - otherwise it has 12 months per line
*
        implicit none
        integer yrbeg,yrend
        real data(12,yrbeg:yrend)
        character file*(*)
        integer i,j,year,month,unit
        real adata,absent
        parameter (absent=3e33)
        logical txtdat
        character line*128
*
        if ( index(file,'txt').ne.0 ) then
            txtdat = .TRUE.
        else
            txtdat = .FALSE.
        endif
        if ( txtdat ) then
            if ( file.ne.'-' ) then
                call rsunit(unit)
                open(unit,file=file,status='old')
            else
                unit = 5
            endif
*           skip header
            do i=1,5
                read(unit,'(a)') file
            enddo
*           data is in the format 
*             (blank line)
*             yyyy
*             text
*             mm value      12 times, value can be 'NO DATA'
*
            do i=1,10000
   10           read(unit,'(a)',err=900,end=100) line
                if ( line.eq.' ' ) goto 10
                read(line,*,err=900,end=100) year
   20           read(unit,'(a)',err=900,end=100) line
                if ( line.eq.' ' ) goto 20
                do j=1,12
   30               read(unit,'(a)',err=900,end=100) line
                    if ( line.eq.' ' ) goto 30
                    if ( index(line,'NO DATA').ne.0 ) then
                        read(line,*,err=900) month
                        adata = absent
                    else
                        read(line,*,err=900) month,adata
                        if ( adata.eq.-99.9 ) adata = absent
                    endif
                    if ( month.ne.j ) then
                        print *,'error reading station data ',line
                    endif
                    data(month,year) = adata
                enddo
            enddo
  100       continue
            if ( unit.ne.5 ) close(unit)
        else
            call readdat(data,yrbeg,yrend,file)
        endif
        return
*
*       error messages
*
  900   print *,'error reading station data ',line
        if ( txtdat ) then
            print*,'expecting: <lf>yyyy<lf>text<lf>12*(month value<lf>)'
        else
            print *,'expecting: yyyy val1 ...val12'
        endif
        call abort
        end
