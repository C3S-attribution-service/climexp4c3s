        program plotpdf
*
*       program to read two datasets and plot the PDF in the following
*       way: sort according to indices, take the first n, print 
*       average index vs percentiles @ 1,...,n-1; shift one in n.
*
        implicit none
        integer yrbeg,yrend,nmax
        parameter (yrbeg=1749,yrend=2020,nmax=100)
        integer i,j,k,ii(12*(yrend-yrbeg+1)),jj(nmax),m,n,lsum,lag,lskip
     +        ,m1,m2
        real data(12,yrbeg:yrend),indx(12,yrbeg:yrend),iindx(12*(yrend
     +        -yrbeg+1)),ddata(12*(yrend-yrbeg+1)),s,y(nmax)
        character line*255,oper*1
        integer iargc
*
        if ( iargc().lt.3 ) then
            print '(2a)','usage: plotpdf dataset index n [month n] '//
     +            '[sum n] [lag n]'
            call abort
        endif
        call makeabsent(data,12,yrbeg,yrend)
        call makeabsent(indx,12,yrbeg,yrend)
        lsum = 1
        lag = 0
        oper = '+'
        m1 = 1
        m2 = 12
        call getarg(1,line)
        call rdtxtdat(line,data,yrbeg,yrend)
        call getarg(2,line)
        call rdtxtdat(line,indx,yrbeg,yrend)
        call getarg(3,line)
        read(line,*) n
        print '(a,i4)','# using n = ',n
        if ( n.gt.nmax ) then
            print *,'Recompile with nmax at least as large as ',n
            call abort
        endif
        lskip = 0
        do i=4,iargc()
            if ( lskip.gt.0 ) then
                lskip = lskip - 1
            else
                call getarg(i,line)
                if ( line(1:3).eq.'sum' ) then
                    call getarg(i+1,line)
                    read(line,*,err=906) lsum
                    print '(a,i3,a)','# summing ',lsum,' months'
                    oper = '+'
                    lskip = 1
                endif
                if ( line(1:3).eq.'lag' ) then
                    call getarg(i+1,line)
                    read(line,*,err=905) lag
                    print '(a,i3)','# using lag ',lag
                    lskip = 1
                endif
                if ( line(1:3).eq.'mon' ) then
                    call getarg(i+1,line)
                    read(line,*,err=905) m1
                    m2 = m1
                    print '(a,i3)','# using month ',m1
                    lskip = 1
                endif
            endif
        enddo
*
        call sumit(data,12,yrbeg,yrend,lsum,oper)
        call sumit(indx,12,yrbeg,yrend,lsum,'+')
*       
*       make linear arrays without absent values
        m = 0
        do i=yrbeg,yrend
            do j=m1,m2
                if ( indx(j,i).lt.1e33 .and. data(j,i).lt.1e33 ) then
                    m = m + 1
                    iindx(m) = indx(j,i)
                    ddata(m) = data(j,i)
                endif
            enddo
        enddo
*
        call ffsort(iindx,ii,m)
        do k=0,m-n
*
*           compute average index
            s = 0
            do i=1,n
                s = s + iindx(ii(k+i))
            enddo
            s = s/n
*
*           get percentiles
            do i=1,n
                y(i) = ddata(ii(k+i))
            enddo
            call ffsort(y,jj,n)
            print '(999f10.4)',s/lsum,((y(jj(i))+y(jj(i+1)))/2,i=1,n-1)
        enddo
  800   continue
*
*       error messages
*
        stop
  905   print *,'error reading lag value ',line
        stop
  906   print *,'error reading sum value ',line
        stop
  907   print *,'error reading begin value ',line
        stop
  908   print *,'error reading end value ',line
        end
