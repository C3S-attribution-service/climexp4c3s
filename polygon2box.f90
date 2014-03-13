program polygon2box
  !
  ! get the boundingbox for an (SREX) polygon
  !
  implicit none
  integer i
  character line*80,file*255
  real x,y,xmin,xmax,ymin,ymax,xold

  call getarg(1,file)
  if ( file.eq.' ' ) then
     write(0,*) 'usage: polygon2box file.txt'
     stop
  end if
  xmin = 3e33
  xmax = -3e33
  ymin = 3e33
  ymax = -3e33
  xold = 3e33
  open(1,file=trim(file),status='old')
1 continue
  read(1,'(a)',end=800) line
  if ( line(1:1).eq.'#' .or. line.eq.' ' ) goto 1
  if ( index(line,'NaN').ne.0 .or. index(line,'nan').ne.0 ) goto 1
  read(line,*) x,y
  if ( xold.lt.1e33 ) then
    if ( abs(x-xold-360).lt.abs(x-xold) ) then
        x = x - 360
    else if ( abs(x-xold+360).lt.abs(x-xold) ) then
        x = x + 360
    end if
  end if
  xold = x
  xmin = min(x,xmin)
  xmax = max(x,xmax)
  ymin = min(y,ymin)
  ymax = max(y,ymax)
  goto 1
800 continue
  print '(a,f10.5)','xmin=',xmin
  print '(a,f10.5)','xmax=',xmax
  print '(a,f10.5)','ymin=',ymin
  print '(a,f10.5)','ymax=',ymax
end program polygon2box
