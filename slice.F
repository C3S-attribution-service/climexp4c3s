        program slice
*
*       slice a 2D field from a 3D GrADS datafile
*
        implicit none
        integer nxmax,nymax,nzmax
        parameter(nxmax=360,nymax=180,nzmax=100)
        integer i,j,k,ninfile,noutfile,nx,ny,nz,nt
        logical xrev,yrev,zrev
        real value,xx(nxmax),yy(nymax),zz(nzmax),field(nxmax,nymax,nzmax
     +        )
        character string*256,infile*256,outfile*256,direction*1
        integer iargc,llen
        common /c_slice/ field
*
*       usage
        if ( iargc().ne.4 ) then
            print *,'usage: slice infile[.ctl|.dat] [lon|lat|lev] value'
     +            //' outfile'
            stop
        endif
*       
*       arguments
        call getarg(2,string)
        if ( string(1:3).eq.'lon' ) then
            direction = 'x'
        elseif ( string(1:3).eq.'lat' ) then
            direction = 'y'
        elseif ( string(1:3).eq.'lev' ) then
            direction = 'z'
        else
            write(0,*) 'error: unrecognized argument, '//
     +            'expected lon,lat or lev but got ',string
            call abort
        endif
        call getarg(3,string)
        read(string,*,err=901) value
*
*       open files
        call getarg(1,infile)
        call getarg(2,outfile)
        call getnfile(outfile,noutfile)
        open(3,file=outfile(1:noutfile)//'.ctl',status='unknown')
        call parsectl(########)

        stop
  901   write(0,*) 'error reading value from ',string
        end
