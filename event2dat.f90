program event2dat
  ! generates a standard dat files with zeros everywhere excpt in a list of events, which are put to one
  implicit none
  integer yr,mo,yrbeg,yrend,nperyear
  real val
  real,allocatable :: data(:,:)
  integer iargc
  character file*255,string*80

  if ( iargc().ne.4 ) then
     print *,'usage: event2dat infile yrbeg yrend nperyear'
     print *,'       sets dates in infile to 1, others to 0'
     stop
  end if
  
  call getarg(1,file)
  call getarg(2,string)
  read(string,*) yrbeg
  call getarg(3,string)
  read(string,*) yrend
  call getarg(4,string)
  read(string,*) nperyear
  if ( nperyear.le.1 .or. nperyear.gt.12 ) then
     write(0,*) 'dat2event: can only handle monthly o seasonal data at the moment,'
     write(0,*) '           not nperyear = ', nperyear
     call abort
  end if
  allocate(data(nperyear,yrbeg:yrend))
  data = 0

  open(1,file=trim(file),status='old')
  do
     read(1,'(a)',end=100) string
     if ( string(1:1).eq.'#' ) then
        print '(a)',trim(string)
     else
        read(string,*,end=100) yr,mo,val
        if ( yr.lt.yrbeg .or. yr.gt.yrend ) then
           write(0,*) 'event2dat: warning: skipping year ',yr
        else if ( mo.lt.1 .or. mo.gt.nperyear ) then
           write(0,*) 'event2dat: error: wrong month ',mo
           call abort
        else
           data(mo,yr) = val
        end if
     end if
  end do
100 continue
  do yr=yrbeg,yrend
     print '(i4,366f6.2)',yr,(data(mo,yr),mo=1,nperyear)
  end do

end program
