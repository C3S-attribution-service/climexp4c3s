      program daily2longerfield
*       
*     compute aggregate quantities from daily data
*     input: daily field
*     output: yearly/monthly/10-dy/5-dy field
*     G.J. van Oldenborgh, KNMI, 2004
*     
#include "recfac.h"
      integer nxmax,nymax,nzmax,mpermax,npermax
      parameter (nxmax=720,nymax=360,nzmax=1,mpermax=366,npermax=73)
      integer nperyear,nperyearnew,yr,mo,dy,i,j,n,dpm(12),itype,
     +     jold,jj,nx,ny,nz,nt,firstyr,firstmo,endian,nvars,ivars(2,1),
     +     iout,irec,ix,iy
      integer*2 nn(nxmax,nymax,mpermax)
      real oldfield(nxmax,nymax,mpermax),newfield(nxmax,nymax,npermax),
     +     oldclim(nxmax,nymax,mpermax),newclim(nxmax,nymax,npermax),
     +     oldseries(mpermax),newseries(npermax),
     +     xx(nxmax),yy(nymax),zz(nzmax),undef,
     +     cut(mpermax),pcut,normdata(npermax),s
      character infile*256,datfile*256,outfile*256,string*40,
     +     lgt*1,oper*3,title*80,vars(1)*20,lvars(1)*80
      logical normal,lwrite,llwrite
      integer iargc,llen
      data dpm/ 31,29,31,30,31,30,31,31,30,31,30,31/
      lwrite = .false.
      if ( iargc().lt.2 ) then
          print *,'usage: daily2longerfield infile.ctl nperyearnew '//
     +         'mean|sd|sum|abo|bel|min|max|num [<> val[%|p]] '//
     +         'outfile.ctl'
          print *,'(more options will come as requested)'
          stop
      endif
*
*     read operation
*
      call getarg(2,string)
      read(string,*,err=901) nperyearnew
      if ( nperyearnew.gt.npermax ) then
          write(0,*) 'daily2longerfield: error: nperyearnew = '
     +         ,nperyearnew,' not yet supported'
          write(*,*) 'daily2longerfield: error: nperyearnew = '
     +         ,nperyearnew,' not yet supported'
          call abort
      endif
      call getarg(3,oper)
      if ( oper.ne.'mea' .and. oper.ne.'sd ' .and.
     +     oper.ne.'min' .and. oper.ne.'max' .and.
     +     oper.ne.'num' .and. oper.ne.'sum' .and.
     +     oper.ne.'bel' .and. oper.ne.'abo' ) then
          write(0,*) 'daily2longerfield: error: unknown operation '
     +         ,oper
          call abort
      endif
      normal = .false.
      if ( iargc().ge.6 ) then
          call getarg(4,lgt)
          if ( lgt.eq.'l' ) lgt = '<'
          if ( lgt.eq.'g' ) lgt = '>'
          if ( lgt.ne.'<' .and. lgt.ne.'>' ) then
              write(0,*)
     +             'daily2longerfield: error: unknown comparison ',lgt
              call abort
          endif
          call getarg(5,string)
          read(string,*,err=902) pcut
          if ( index(string,'%').ne.0 .or. index(string,'p').ne.0 ) 
     +         then
              write(0,*) 'percentage cuts not yet implemented'
              write(*,*) 'percentage cuts not yet implemented'
              call abort
          else
              do j=1,nperyear
                  cut(j) = pcut
              enddo
          endif
      else
          lgt = ' '
      endif
*
*     get field parameters
*
      call getarg(1,infile)
      call parsectl(infile,datfile,nxmax,nx,xx,nymax,ny,yy,nzmax,nz
     +            ,zz,nt,nperyear,firstyr,firstmo,undef,endian,title,1
     +            ,nvars,vars,ivars,lvars,units)
      itype = 0
      call getarg(iargc(),outfile)
      i = index(outfile,'.ctl')
      outfile(i:) = '.grd'
      call rsunit(iout)
      open(iout,file=outfile,form='unformatted',access='direct',
     +     recl=recfa4*nx*ny)
      irec = 0
*
*     treat each year individually - we do not have enough memory to
*     read the whole field in one go
*
*     first read the climatology if it exists
*
      open(1,file='oldclim.grd',status='old',err=100,access='direct',
     +     recl=nxmax*nymax*mpermax*recfa4)          
      if ( lwrite ) print *,'reading climatologies from file'
      read(1,rec=1) oldclim
      close(1)
      goto 200
 100  continue
      stop
      if ( lwrite ) print *,'setting oldclim to zero'
      do i=1,nperyear
          do iy=1,ny
              do ix=1,nx
                  oldclim(ix,iy,i) = 0
              enddo
          enddo
      enddo
      if ( lwrite ) print *,'setting nn to zero'
      do i=1,nperyear
          do iy=1,ny
              do ix=1,nx
                  nn(ix,iy,i) = 0
              enddo
          enddo
      enddo
      do yr=firstyr,firstyr+(nt+firstmo-2)/nperyear
          call keepalive(yr-firstyr+1,(nt+firstmo-2)/nperyear+1)
          if ( lwrite ) print *,'reading ',yr,' for climatology'
          call readdatfile(datfile,oldfield,nxmax,nymax,nx,ny,nperyear
     +         ,yr,yr,firstyr,firstmo,nt,undef,endian,lwrite,yr,yr,1,1)
          do i=1,nperyear
              if ( lwrite ) print *,'computing climatology field ',i
              do iy=1,ny
                  do ix=1,nx
                      if ( oldfield(ix,iy,i).lt.1e33 ) then
                          nn(ix,iy,i) = nn(ix,iy,i) + 1
                          oldclim(ix,iy,i) = oldclim(ix,iy,i) + 
     +                         oldfield(ix,iy,i)
                      endif
                  enddo
              enddo
          enddo
      enddo
      if ( lwrite ) print *,'normalizing climatology field'
      do i=1,nperyear
          do iy=1,ny
              do ix=1,nx
                  if ( nn(ix,iy,i).gt.0 ) then
                      oldclim(ix,iy,i) = oldclim(ix,iy,i)/nn(ix,iy,i)
                  else
                      oldclim(ix,iy,i) = 3e33
                  endif
              enddo
          enddo
      enddo
      if ( .true. ) then
          open(1,file='oldclim.grd',access='direct',
     +         recl=recfa4*nxmax*nymax*mpermax)
          write(1,rec=1) oldclim
          close(1)
      endif
 200  continue
*
*     interpolate the climatology
*
      do iy=1,ny
          do ix=1,nx
              do i=1,nperyear
                  oldseries(i) = oldclim(ix,iy,i)
              enddo
              call allday2period(
     +             oldseries,mpermax,nperyear,
     +             newseries,npermax,nperyearnew,
     +             1,1,oper,lgt,cut,itype,.false.)
              do i=1,nperyearnew
                  newclim(ix,iy,i) = newseries(i)
              enddo
          enddo
      enddo
      if ( .false. ) then
          open(1,file='newclim.grd',access='direct',
     +         recl=recfa4*nxmax*nymax*npermax)
          write(1,rec=1) newclim
          close(1)
      endif
*
*     next do the summing etc.
*
      do yr=firstyr,firstyr+(nt+firstmo-2)/nperyear
          call keepalive(yr-firstyr+1,(nt+firstmo-2)/nperyear+1)
          if ( lwrite ) print *,yr
          call readdatfile(datfile,oldfield,nxmax,nymax,nx,ny,nperyear
     +         ,yr,yr,firstyr,firstmo,nt,undef,endian,lwrite,yr,yr,1,1)
          do iy=1,ny
              do ix=1,nx
                  do i=1,nperyear
                      if ( oldfield(ix,iy,i).lt.1e33 .and. 
     +                     oldclim(ix,iy,i).lt.1e33 ) then
                          oldseries(i) = oldfield(ix,iy,i) -
     +                         oldclim(ix,iy,i)
                      else
                          oldseries(i) = 3e33
                      endif
                  enddo
                  if ( .false. ) then
                  if ( iy.gt.100 .and. oldseries(73).ne.0 .and.
     +                 oldseries(73).lt.1e33 ) then
                      print *,'oldseries,oldfield,oldclim = ',ix,iy,i
                      do i=1,nperyear
                          print *,i,oldseries(i),oldfield(ix,iy,i),
     +                         oldclim(ix,iy,i)
                      enddo
                      llwrite = lwrite
                  endif
                  else
                      llwrite = .false.
                  endif
                  call allday2period(
     +                 oldseries,mpermax,nperyear,
     +                 newseries,npermax,nperyearnew,
     +                 yr,yr,oper,lgt,cut,itype,llwrite)
                  if ( llwrite ) then
                      stop
                  endif
*                 
*                 copy back
*
                  do i=1,nperyearnew
                      if ( newseries(i).lt.1e33 .and. 
     +                     newclim(ix,iy,i).lt.1e33 ) then
                          newfield(ix,iy,i) = newseries(i) + 
     +                         newclim(ix,iy,i)
                      else
                          newfield(ix,iy,i) = 3e33
                      endif
                  enddo
              enddo             ! ix
          enddo                 ! iy
*
*         save output
*         
          do n = 1,nperyearnew
              irec = irec + 1
              write(iout,rec=irec) ((newfield(ix,iy,n),ix=1,nx),iy=1,ny)
          enddo
      enddo                     ! yr
      datfile = outfile
      i = index(outfile,'.grd')
      outfile(i:) = '.ctl'
      call writectl(outfile,datfile,nx,xx,ny,yy,nz,zz,
     +     irec,nperyearnew,firstyr,1,3e33,title,nvars,vars,ivars
     +     ,lvars,units)
      goto 999
*
*     error messages
*     
      goto 999
 901  write(0,*) 'daily2longerfield: expecting nperyearnew, not ',string
      call abort
 902  write(0,*) 'daily2longerfield: expecting value[%|p], not ',string
      call abort
 999  continue
      end
