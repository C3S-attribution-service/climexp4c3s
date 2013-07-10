        program multifit
*
*       read a correlate-dumpfile and 
*       do a simultaneous fit to data with n predictors
*       and write the resulting optimal compound index on stdout
*       in dump format
*
*       Uses Numerical Recipes (I am lazy)
*
        implicit none
        integer indxmx,ndatmx
        parameter (indxmx=12,ndatmx=12000)
        integer i,j,k,n,year(ndatmx),month(ndatmx)
        real data(ndatmx),indx(ndatmx,indxmx),chisq,sig(ndatmx),
     +        x(ndatmx),u(ndatmx,indxmx),v(indxmx,indxmx),w(indxmx),
     +        a(indxmx),da(indxmx,indxmx),nindx(ndatmx)
        common /cindx/ indx
        logical lwrite
        character string*255,indname(indxmx)*20
        integer iargc,llen
        external findx
*
*       process arguments
*
        lwrite = .FALSE.
        n = iargc()
        if ( n.lt.2 ) then
            print *,'usage: multifit dumpfile index1 index2 ...'
            stop
        endif
        call getarg(1,string)
        i = index(string,' ')
        if ( lwrite ) print *,'opening old dumpfile ',string(1:i-1)
        open(1,file=string,status='old')
        string(i:i) = '1'
        if ( lwrite ) print *,'opening new dumpfile ',string(1:i)
        open(2,file=string)
        indname(1) = 'constant'
        do i=2,n
            call getarg(i,indname(i))
        enddo
*
*       read data
*       
        i = 0
  100   continue
        read (1,'(a)',end=200,err=900) string
        if ( string(1:1).eq.'#' ) then
            if ( lwrite ) print '(a)',string(1:llen(string))
            write(2,'(a)') string(1:llen(string))
            goto 100
        endif
        if ( string.eq.' ' ) goto 100
        if ( index(string,'999.9').ne.0 ) goto 100
        i = i + 1
        read(string,*,err=901)(indx(i,j),j=1,n-1),data(i),year(i)
     +        ,month(i)
        goto 100
  200   continue
        close(1)
*
*       call fit routine
*
        if ( i.lt.n ) goto 902
        do j=1,i
            x(j) = j
        enddo
        do j=1,i
            sig(j) = 1
        enddo
        if ( lwrite ) print *,'calling svdfit'
        call svdfit(x,data,sig,i,a,n,u,v,w,ndatmx,indxmx,chisq,findx)
        print '(a,f12.2,a,i5,a,f12.4)','chisq/ndf = ',chisq,'/',i-n
     +        ,' = ',chisq/(i-n)
*
*       get covariances of fit parameters
*
        if ( lwrite ) print *,'calling svdvar'
        call svdvar(v,n,indxmx,w,da,indxmx)
*
*       print output
*
        print '(a)','Best fit:'
        do j=1,n
            print '(a,g16.6,a,g16.6)',indname(j),a(j),' +/- ',
     +            sqrt(da(j,j))
        enddo
        print '(a)','Correlation matrix:'
        do j=1,n
            print '(a,20f7.3)',indname(j),(da(k,j)/sqrt(da(j,j)*da(k,k)
     +            ),k=1,j)
        enddo
*
*       construct new index (without constant term)
*
        do j=1,i
            nindx(j) = 0
            do k=2,n
                nindx(j) = nindx(j) + a(k)*indx(j,k-1)
            enddo
        enddo
*
*       write new dumpfile
*       
        do j=1,i
            write(2,'(2f12.4,2i5)') nindx(j),data(j),year(j),month(j)
        enddo
        write(2,'(a,i2.2,a,f10.6)') ('# a',k,'=',a(k),k=1,n)
        close(2)
*
*       error messages
*
        goto 999
  900   print *,'error reading dumpfile at ',string
        call abort
  901   print *,'error reading data from string ',string
        call abort
  902   print *,'not enough data points to fit ',i,n
        call abort
  999   continue
        end
*        
        integer function llen(a)
        character*(*) a
        do 10 i=len(a),1,-1
            if(a(i:i).ne.'?' .and. a(i:i).ne.' ')goto 20
   10   continue
        llen=len(a)
   20   continue
        llen = i
        end
*       
        subroutine findx(xi,f,n)
        implicit none
        integer n
        real xi,f(n)
        integer indxmx,ndatmx
        parameter (indxmx=12,ndatmx=12000)
        real indx(ndatmx,indxmx)
        common /cindx/ indx
        integer i,j
*
        i = nint(xi)
        if ( abs(xi-i).gt.0.01 ) print *,'findx: wrong input! ',xi
        if ( n.gt.indxmx ) print *,'findx: wrong input! ',n,indxmx
        f(1) = 1
        do j=2,n
            f(j) = indx(i,j-1)
        enddo
        end
