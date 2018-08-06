program svd

!       program make an SVD of two fields, modelled on
!       correlatefieldfield and eof.
!       this version only accepts a single region
!       and a single starting month

    implicit none
    include 'params.h'
    include 'getopts.inc'
    include 'netcdf.inc'
    integer,parameter :: nsvdmax=20,recfa4=4
    integer :: i,ii,j,jj,k,m,n,iens,mo,yr,lwork,info,ivars(2,nsvdmax), &
        nsvd,nxy, &
        ncid1,nx1,ny1,nz1,nt1,nperyear1,firstyr1,firstmo1,nvars1, &
        ivars1(6,1),endian1,mens11,mens1,nxf1,nyf1,nzf1, &
        nens11,nens12,nxymax1, &
        ncid2,nx2,ny2,nz2,nt2,nperyear2,firstyr2,firstmo2,nvars2, &
        ivars2(6,1),endian2,mens12,mens2,nxf2,nyf2,nzf2, &
        nens21,nens22,nxymax2, &
        status,mens,firstyr,firstmo,nperyear,nt,lastyr, &
        i1,j1,n1,nn1,jx1,jy1,jz1,x11,y11,z11,x12,y12,z12, &
        i2,j2,n2,nn2,jx2,jy2,jz2,x21,y21,z21,x22,y22,z22
    integer,allocatable :: ijk1(:,:),ijk2(:,:)
    real :: xx1(nxmax),yy1(nymax),zz1(nzmax), &
        xx2(nxmax),yy2(nymax),zz2(nzmax), &
        wx1(nxmax),wy1(nymax),wz1(nzmax), &
        wx2(nxmax),wy2(nymax),wz2(nzmax), &
        undef1,undef2,tarray(3),s,f1,f2,t,a,b,scratch(1),w
    real,allocatable :: &
        field1(:,:,:,:,:,:), &
        field2(:,:,:,:,:,:), &
        cov(:,:),ss(:),u(:,:),vt(:,:),work(:), &
        svd1(:,:,:,:),svd2(:,:,:,:), &
        pattern1(:,:,:,:,:),pattern2(:,:,:,:,:), &
        pc1(:,:,:,:),pc2(:,:,:,:)
    real,allocatable :: fxy1(:,:,:),fxy2(:,:,:)
    character :: infile1*256,infile2*256,datfile1*256,datfile2*256, &
        title1*256,vars1(1)*40,svars1(1)*20,lvars1(1)*80, &
        units1(1)*60,cell_methods1*100,metadata1(2,100)*2000, &
        title2*256,vars2(1)*40,svars2(1)*20,lvars2(1)*80, &
        units2(1)*60,cell_methods2*100,metadata2(2,100)*2000, &
        title*256,dir*255,outfile1*256,outfile2*256,yesno,file*256
    character :: lz1(3)*20,lz2(3)*20,ltime1*20,ltime2*20,history1*1024 &
        ,history2*1024,string*80
    character vars(nsvdmax)*10,lvars(nsvdmax)*80,units(nsvdmax)
    logical :: lexist1,lexist2,xrev1,yrev1,xwrap1,xrev2,yrev2,xwrap2
    real :: etime

!       check arguments

    lwrite = .false. 
    n = command_argument_count()
    if ( n < 4 ) then
        print *,'usage: svd '// &
            'field1.[ctl|nc] field2.[ctl|nc] nsvd '// &
            '[lag n[:m]] [sum|ave|max|min|sel n] '// &
            '[log|sqrt|rank] '// &
            '[minfac r] [minnum n] [begin yr] [end yr] '// &
            '[diff [nyr]] [detrend] '// &
            'outfield1.ctl outfield2.ctl'
        call exit(-1)
    end if
    call get_command_argument(1,infile1)
    call getmetadata(infile1,mens11,mens1,ncid1,datfile1,nxmax,nx1 &
    ,xx1,nymax,ny1,yy1,nzmax,nz1,zz1,lz1,nt1,nperyear1,firstyr1 &
    ,firstmo1,ltime1,undef1,endian1,title1,history1,1,nvars1 &
    ,vars1,ivars1,lvars1,svars1,units1,cell_methods1,metadata1 &
    ,lwrite)

    call get_command_argument(2,infile2)
    call getmetadata(infile2,mens12,mens2,ncid2,datfile2,nxmax,nx2 &
    ,xx2,nymax,ny2,yy2,nzmax,nz2,zz2,lz2,nt2,nperyear2,firstyr2 &
    ,firstmo2,ltime2,undef2,endian2,title2,history2,1,nvars2 &
    ,vars2,ivars2,lvars2,svars2,units2,cell_methods2,metadata2 &
    ,lwrite)
    call get_command_argument(3,string)
    read(string,*,err=10) nsvd
    i1 = 4
    if ( nsvd > nsvdmax ) then
        write(0,*) 'svd: error: max number of SVDs is ',nsvdmax
        call exit(-1)
    end if
    goto 11
 10 continue
    i1 = 3
    nsvd = 4
 11 continue
    print '(a,i3,a)','Computing ',nsvd,' SVDs'

    call get_command_argument(command_argument_count()-1,outfile1)
    i = index(outfile1,'.ctl')
    if ( i == 0 ) then
        write(0,*) 'svd: sorry. netcdf output not yet supported'
        call exit(-1)
    end if
    call get_command_argument(command_argument_count(),outfile2)
    i = index(outfile2,'.ctl')
    if ( i == 0 ) then
        write(0,*) 'svd: sorry. netcdf output not yet supported'
        call exit(-1)
    end if
    inquire(file=trim(outfile1),exist=lexist1)
    inquire(file=trim(outfile2),exist=lexist2)
    if ( lexist1 .or. lexist2 ) then
        print *,'output file ',trim(outfile1),' or ',trim(outfile2), &
        ' already exists, overwrite? [y/n]'
        read(*,'(a)') yesno
        if (  yesno /= 'y' .and. yesno /= 'Y' .and. &
        yesno /= 'j' .and. yesno /= 'J' ) then
            call exit(-1)
        end if
        if ( lexist1 ) then
            open(1,file=trim(outfile1))
            close(1,status='delete')
        end if
        if ( lexist2 ) then
            open(1,file=trim(outfile2))
            close(1,status='delete')
        end if
    end if

    firstyr = max(firstyr1,firstyr2)
    if ( nperyear1 /= nperyear2 ) then
        write(0,*) 'svd: error: cannot interpolate in time yet',nperyear1,nperyear2
        write(*,*) 'svd: error: cannot interpolate in time yet',nperyear1,nperyear2
        call exit(-1)
    end if
    nperyear = nperyear1
    lastyr = min(firstyr1 + (nt1-1)/nperyear,firstyr2 + (nt2-1) &
    /nperyear)

!       save time on the initialization - but not too much.
    nt = nperyear*(lastyr-firstyr+1)
    n = command_argument_count()
    mens1 = min(mens11,mens12)
    mens = max(mens1,mens2)
    call getopts(i1,n-2,nperyear,yrbeg,yrend, .true. ,mens1,mens2)
    if ( lag1 < 0 ) print *,'(field1 leading field2)'
    if ( lag2 > 0 ) print *,'(field2 leading field1)'
    if ( dump ) write(0,*)'svd: dump not supported'
    if ( plot ) write(0,*)'svd: plot not supported'
    if ( lks ) write(0,*)'svd: K-S not supported'
    if ( lconting ) write(0,*)'svd: contingency '// &
    'tables not supported'
    if ( composite ) write(0,*)'composites not yet supported'
    do i=1,indxuse
        if ( lincl(i) ) write(0,*)'svd: what do ', &
        'you mean with ',strindx(i),'?'
    end do
    if ( fix2 ) then
        write(0,*) 'svd: fix2 not yet supported'
        call exit(-1)
    end if
    yr1 = max(yr1,firstyr,firstyr - (min(lag1,lag2)+nperyear-1) &
    /nperyear)
    yr2 = min(yr2,lastyr,lastyr - (max(lag1,lag2)-nperyear+1) &
    /nperyear)
    if ( mens1 == 0 ) then
        nens11 = 0
        nens21 = 0
    else
        nens11 = max(nens1,mens11)
        nens21 = min(nens2,mens1)
    end if
    if ( mens2 == 0 ) then
        nens12 = 0
        nens22 = 0
    else
        nens12 = max(nens1,mens12)
        nens22 = min(nens2,mens2)
    end if
    if ( lwrite ) then
        print *,'svd: correlating ',trim(infile1)
        print *,'     with ',trim(infile2)
        print *,'     years: ',yr1,yr2
    end if

!       allocate fields

    firstyr = yr1
    lastyr = yr2
    mens1 = min(mens1,nens2)
    mens2 = min(mens2,nens2)
    mens = max(mens1,mens2)
    if ( lwrite ) print *,'allocating field1 ',nx1,ny1,nz1,nperyear &
    ,firstyr,lastyr,mens1,mens2,mens
    nxf1 = nx1
    nyf1 = ny1
    nzf1 = nz1
    allocate(field1(nx1,ny1,nz1,nperyear,firstyr:lastyr,0:mens1))
    if ( lwrite ) print *,'allocating field2 ',nx2,ny2,nz2,nperyear &
    ,firstyr,lastyr,mens1,mens2,mens
    nxf2 = nx2
    nyf2 = ny2
    nzf2 = nz2
    allocate(field2(nx2,ny2,nz2,nperyear,firstyr:lastyr,0:mens2))
    allocate(fxy1(nperyear,firstyr:lastyr,0:mens))
    allocate(fxy2(nperyear,firstyr:lastyr,0:mens))

!       read kill info file (in the climexp) and add own PID

    call killfile(dir,title,file,0)

!       init

    if ( m1 /= m2 ) then
        print *,'Sorry, can only handle a single starting month'
        print *,'months:',m1,m2
        call exit(-1)
    end if
    if ( lag1 /= lag2 ) then
        print *,'Sorry, can only handle a single lag'
        print *,'lags:',lag1,lag2
        call exit(-1)
    end if
    if ( lag2 /= lag1 ) then
        print *,'Sorry, only makes sense for lag1=lag2'
        call exit(-1)
    end if
    print *,'init'
    do iens=0,mens
        call makeabsent(fxy1(1,firstyr,iens),nperyear,firstyr,lastyr &
        )
        call makeabsent(fxy2(1,firstyr,iens),nperyear,firstyr,lastyr &
        )
    end do

!       compute minfac if it has not been set explicitly

    if ( minfac < 0 .and. minnum < 0 ) then
    !           heuristic, gives 0.25 for 150 yrs, 0.5 for 50 yrs, 0.75 for
    !           20yrs
        minfac = max(0.1, &
        min(0.6, &
        &            1.5-log(1+real(min(nt,nperyear*(yr2-yr1+1))-1) &
        /nperyear)/4))
    end if
    write(0,'(a,i2,a)') 'Requiring at least ', &
    nint(100*minfac),'% valid points<p>'

!       read fields

    do iens=nens11,nens21
        call keepalive(1+iens-nens11,2+nens21-nens11+nens22-nens12)
        if ( ncid1 == -1 ) then
            if ( iens > nens11 ) then
                file = infile1
                call filloutens(file,iens)
                call parsectl(file,datfile1,nxmax,nx1,xx1,nymax,ny1 &
                ,yy1,nzmax,nz1,zz1,nt1,nperyear1,firstyr1 &
                ,firstmo1,undef1,endian1,title,1,nvars1,vars1 &
                ,ivars1,lvars1,units1)
            end if
            call zreaddatfile(datfile1,field1(1,1,1,1,firstyr,iens) &
            ,nx1,ny1,nz1,nx1,ny1,nz1,nperyear,firstyr,lastyr &
            ,firstyr1,firstmo1,nt1,undef1,endian1,lwrite,yr1 &
            ,yr2,1,1)
        else
            if ( iens > nens11 ) then
                file = infile1
                call filloutens(file,iens)
                call parsenc(file,ncid1,nxmax,nx1,xx1,nymax,ny1 &
                ,yy1,nzmax,nz1,zz1,nt1,nperyear1,firstyr1 &
                ,firstmo1,undef1,title1,1,nvars1,vars1,ivars1 &
                ,lvars1,units1)
            end if
            call zreadncfile(ncid1,field1(1,1,1,1,firstyr,iens) &
            ,nx1,ny1,nz1,nx1,ny1,nz1,nperyear,firstyr,lastyr &
            ,firstyr1,firstmo1,nt1,undef1,lwrite,yr1,yr2,ivars1 &
            )
        end if
    end do
    if ( lwrite ) then
        print *,'field1 @ 0,60N'
        call dump060(xx1,yy1,zz1,field1,nxf1,nyf1,nzf1,nx1,ny1,nz1 &
        ,nperyear,firstyr,lastyr)
    end if
    do iens=nens12,nens22
        call keepalive(1+iens-nens12+1+nens21-nens11, &
        &            2+nens21-nens11+nens22-nens12)
        if ( ncid2 == -1 ) then
            if ( iens > nens12 ) then
                file=infile2
                call filloutens(file,iens)
                call parsenc(file,ncid2,nxmax,nx2,xx2,nymax,ny2 &
                ,yy2,nzmax,nz2,zz2,nt2,nperyear2,firstyr2 &
                ,firstmo2,undef2,title2,1,nvars2,vars2,ivars2 &
                ,lvars2,units2)
            end if
            call zreaddatfile(datfile2,field2(1,1,1,1,firstyr,iens) &
            ,nx2,ny2,nz2,nx2,ny2,nz2,nperyear,firstyr,lastyr &
            ,firstyr2,firstmo2,nt2,undef2,endian2,lwrite,yr1 &
            ,yr2,1,1)
        else
            if ( iens > nens12 ) then
                file=infile2
                call filloutens(file,iens)
                call parsenc(file,ncid2,nxmax,nx2,xx2,nymax,ny2 &
                ,yy2,nzmax,nz2,zz2,nt1,nperyear,firstyr2 &
                ,firstmo2,undef2,title,1,nvars2,vars2,ivars2 &
                ,lvars2,units2)
            end if
            call zreadncfile(ncid2,field2(1,1,1,1,firstyr,iens) &
            ,nx2,ny2,nz2,nx2,ny2,nz2,nperyear,firstyr,lastyr &
            ,firstyr2,firstmo2,nt2,undef2,lwrite,yr1,yr2,ivars2 &
            )
        end if
    end do
    if ( lwrite ) then
        print *,'field2 @ 0,60N'
        call dump060(xx2,yy2,zz2,field2,nxf2,nyf2,nzf2,nx2,ny2,nz2 &
        ,nperyear,firstyr,lastyr)
    end if
    call keepalive(0,0)
    write(0,'(a,f8.2,a)') 'Averaging, shifting, cutting, time: ' &
    ,etime(tarray),'<br>'

!	get boundaries in grid points of field1

    call getxyprop(xx1,nx1,yy1,ny1,xrev1,yrev1,xwrap1)
    call getlatlonwindow(lat1,lat2,lon1,lon2,xx1,nx1,xwrap1,avex,yy1 &
    ,ny1,avey,x11,x12,y11,y12,lwrite)
    call getlevwindow(lev1,lev2,zz1,nz1,z11,z12,lwrite)

!       average, cut out window - everything to make the arrays smaller

    if ( nz1 > 1 ) then
        write(0,*) 'not yet ready for Z'
        call exit(-1)
    else
        call enscutoutwindow(x11,x12,y11,y12,xx1,nx1,xwrap1,xrev1 &
        ,avex,yy1,ny1,avey,wx1,wy1,field1,nxf1,nyf1,nens1,nens2 &
        ,nperyear,firstyr,lastyr,yr1,yr2,lwrite)
    end if
    call getweights('z',zz1,wz1,nz1, .false. ,lwrite)

!	get boundaries in grid points of field2

    call getxyprop(xx2,nx2,yy2,ny2,xrev2,yrev2,xwrap2)
    call getlatlonwindow(altlat1,altlat2,altlon1,altlon2,xx2,nx2 &
    ,xwrap2,altavex,yy2,ny2,altavey,x21,x22,y21,y22,lwrite)
    call getlevwindow(altlev1,altlev2,zz2,nz2,z21,z22,lwrite)

!       average, cut out window - everything to make the arrays smaller

    write(0,'(a,f8.2,a)') 'Averaging, shifting, cutting, time: ' &
    ,etime(tarray),'<br>'
    if ( nz2 > 1 ) then
        write(0,*) 'not yet ready for Z'
        call exit(-1)
    else
        call enscutoutwindow(x21,x22,y21,y22,xx2,nx2,xwrap2,xrev2 &
        ,altavex,yy2,ny2,altavey,wx2,wy2,field2,nxf2,nyf2,nens1 &
        ,nens2,nperyear,firstyr,lastyr,yr1,yr2,lwrite)
    end if
    call getweights('z',zz2,wz2,nz2, .false. ,lwrite)

!       compute covariance matrix

    call processoptions(field1,nxf1,nyf1,nzf1,nperyear,firstyr &
    ,lastyr,mens,nx1,ny1,nz1,nperyear,1)
    call processoptions(field2,nxf2,nyf2,nzf2,nperyear,firstyr &
    ,lastyr,mens,nx2,ny2,nz2,nperyear,2)

!       estimate time to completion

    nxy=min(nx1*ny1*nz1,nx2*ny2*nz2)
    if ( m1 == 0 ) then
        n = nt
    else
        n = (yr2-yr1+1)*(m2+lsel-m1)
    end if
!       these are ancient values, I'll have to re-time them
    a = 4*40e-9
    b = 25e-9
    t = a*n*(nens2-nens1+1)*real(nxy)**2 + b*real(nxy)**3
    t = t/3 ! for the new server
    write(0,'(a,f8.0,a)') 'Estimated time to completion: ',t &
    ,'s<p>'

!       covariance

    write(0,'(a,f8.2,a)') 'Computing covariance matrix, time: ' &
    ,etime(tarray),'<br>'
    nxymax1=nx1*ny1*nz1
    nxymax2=nx2*ny2*nz2
    if ( lwrite ) print *,'allocating cov(',nxymax1,nxymax2,')'
    allocate(cov(nxymax1,nxymax2))
    allocate(ijk1(3,nxymax1))
    allocate(ijk2(3,nxymax2))
    cov = 3e33
    n1 = 0
    n2 = 0
    do jx1=1,nx1
        do jy1=1,ny1
            do jz1=z11,z12
                call keepalive(jz1+(jy1-1+(jx1-1)*ny1)*nz1, &
                nx1*ny1*nz1)
                do jx2=1,nx2
                    do jy2=1,ny2
                        do jz2=z21,z22
                            if ( lwrite ) print '(a,3i4,a,3i4)', &
                            'computing covariance between ', &
                            jx1,jy1,jz1,' and ',jx2,jy2,jz2
                            s = 0
                            n = 0
                            mo=m1
                            do iens=nens1,nens2
                                do yr=yr1,yr2
                                    if ( mo == 0 ) then
                                        j1 = 1
                                        j2 = nperyear
                                    else
                                        j1 = mo
                                        j2 = mo + lsel - 1
                                    end if
                                    do jj=j1,j2
                                        j = jj
                                        call normon(j,yr,i,nperyear)
                                        if ( i < firstyr .or. &
                                        i > lastyr ) &
                                        goto 710
                                        m = j-lag1
                                        call normon(m,i,ii,nperyear)
                                        if ( ii < firstyr .or. &
                                        ii > lastyr ) &
                                        goto 710
                                        f1 = field1(jx1,jy1,jz1,j,i &
                                        ,iens)
                                        f2 = field2(jx2,jy2,jz2,m,ii &
                                        ,iens)
                                        if ( f1 < 1e33 .and. ( &
                                        ((f1 <= maxdata).eqv. &
                                        (f1 >= mindata)) .eqv. &
                                        (maxdata >= mindata) ) &
                                         .and. &
                                        f2 < 1e33 .and. ( &
                                        ((f2 <= maxindx).eqv. &
                                        (f2 >= minindx) .eqv. &
                                        (maxindx >= minindx) ) &
                                        ) ) then
                                            n = n + 1
                                            s = s + f1*f2
                                        end if
                                        710 continue
                                    end do ! jj
                                end do ! yr
                            end do ! iens
                            if ( lwrite ) then
                                if ( mo == 0 ) then
                                    print *,'Comparing n=',n &
                                    ,' with minfac*N = ',minfac &
                                    ,min(nt,nperyear*(yr2-yr1+1 &
                                    ))
                                else
                                    print *,'Comparing n=',n &
                                    ,' with minfac*N = ',minfac &
                                    ,min(nt/nperyear,yr2-yr1+1) &
                                    *lsel
                                end if
                            end if
                            if ( mo == 0 .and. n < minfac*min(nt,nperyear*(yr2-yr1+1)) .or. &
                                mo /= 0 .and. n < minfac* &
                                min(nt/nperyear,yr2-yr1+1)*lsel .or. n < minnum ) then
                                if ( lwrite ) print '(a,6i5,a,2i6)' &
                                ,'not enough valid points at ' &
                                ,jx1,jy1,jz1,jx2,jy2,jz2,': ',n &
                                ,nt
                                goto 790
                            end if
                        !                               first coordinate
                            if ( n1 == 0 ) then
                                n1 = 1
                                ijk1(1,n1) = jx1
                                ijk1(2,n1) = jy1
                                ijk1(3,n1) = jz1
                            elseif ( ijk1(1,n1) /= jx1 .or. &
                                ijk1(2,n1) /= jy1 .or. &
                                ijk1(3,n1) /= jz1 ) then
                                n1 = n1 + 1
                                if ( n1 > nxymax1 ) then
                                    write(0,*) 'svd: error: too '// &
                                    'many points in field1'
                                    write(*,*) 'svd: error: too '// &
                                    'many points in field1'
                                    call exit(-1)
                                end if
                                ijk1(1,n1) = jx1
                                ijk1(2,n1) = jy1
                                ijk1(3,n1) = jz1
                            end if
                        !                               second coordinate
                            do i2=n2,1,-1
                                if ( ijk2(1,i2) == jx2 .and. &
                                ijk2(2,i2) == jy2 .and. &
                                ijk2(3,i2) == jz2 ) &
                                goto 780
                            end do
                            n2 = n2 + 1
                            if ( n2 > nxymax2 ) then
                                write(0,*) 'svd: error: too '// &
                                'many points in field2', &
                                jx2,jy2,jz2
                                write(*,*) 'svd: error: too '// &
                                'many points in field2', &
                                jx2,jy2,jz2
                                call exit(-1)
                            end if
                            ijk2(1,n2) = jx2
                            ijk2(2,n2) = jy2
                            ijk2(3,n2) = jz2
                            i2 = n2
                            780 continue
                            cov(n1,i2) = sqrt(wx1(jx1)*wx2(jx2)* &
                            wy1(jy1)*wy2(jy2)*wz1(jz1)*wz2(jz2) &
                            )*s/(n-1)
                            if (lwrite) print &
                            '(a,6i4,a,i6,a,2i6,a,6f9.4)' &
                            ,'point',jx1,jy1,jz1,jx2,jy2,jz2 &
                            ,' OK (',n,'): cov(',n1,i2,') = ' &
                            ,cov(n1,i2)
                            790 continue ! valid point
                        end do ! nz2
                    end do   ! ny2
                end do       ! nx2
            end do           ! nz1
        end do               ! ny1
    end do                   ! nx1
    print *,'I now have a covariance matrix of ',n1,' by ',n2

!       I now have a covariance matrix with invalid data
!       get rid of the rows and columns with too many undefines
!       and put the rest of the undefineds to zero (I do not know a
!       better value)

    if ( n1 == 0 .or. n2 == 0 ) then
        write(0,*) 'svd: error: no valid data found. '// &
        'Please check your input.'
        write(0,*) 'Common errors are selecting a land area when '// &
        'there is only data over sea or vice vera, a time '// &
        'period with data or a season without data.'
        call exit(-1)
    end if
    801 continue
    do i1=1,n1
        do i2=1,n2
            if ( cov(i1,i2) > 1e30 ) then
                nn1 = 0
                do j1=1,n1
                    if ( cov(j1,i2) > 1e33 ) nn1 = nn1 + 1
                end do
                nn2 = 0
                do j2=1,n2
                    if ( cov(j1,i2) > 1e33 ) nn2 = nn2 + 1
                end do
                if ( nn1 > n1/2 .or. nn2 > n2/2 ) then
                    if ( nn1/real(n1) > nn2/real(n2) ) then
                        if ( lwrite ) print *,'deleting row ',i1 &
                        ,' with ',nn1,'/',n1,' undefs'
                        n1 = n1 - 1
                        do j1=i1,n1
                            do k=1,3
                                ijk1(k,j1) = ijk1(k,j1+1)
                            end do
                            do j2=1,n2
                                cov(j1,j2) = cov(j1+1,j2)
                            end do
                        end do
                    else
                        if ( lwrite ) print *,'deleting col ',i2 &
                        ,' with ',nn2,'/',n2,' undefs'
                        n2 = n2 - 1
                        do j2=i2,n2
                            do k=1,3
                                ijk2(k,j2) = ijk2(k,j2+1)
                            end do
                            do j1=1,n1
                                cov(j1,j2) = cov(j1,j2+1)
                            end do
                        end do
                    end if
                    goto 801 ! for safety
                end if
            end if
        end do
    end do
    n = 0
    do i1=1,n1
        do i2=1,n2
            if ( cov(i1,i2) > 1e30 ) then
                n = n + 1
                cov(i1,i2) = 0
            end if
        end do
    end do
    print '(a,i6,a)','svd: set ',n,' undefined covariances to zero'
    print *,'I now have a covariance matrix of ',n1,' by ',n2

!       SVD

    write(0,'(a,f8.2,a)') 'Computing SVD decomposition, time: ' &
    ,etime(tarray),'<br>'

    n = min(n1,n2)
    if ( 4*nxf1*nxf2*nzf1*n > 3e9 ) then
        write(0,*) 'svd: error: too much memory needed (', 4*nxf1 &
        *nxf2*nzf1*n/1e9,' GB). Please reduce the size of ', &
        'the problem by averaging over grid boxes.'
        write(*,*) 'svd: error: too much memory needed (', 4*nxf1 &
        *nxf2*nzf1*n/1e9,' GB). Please reduce the size of ', &
        'the problem by averaging over grid boxes.'
        call exit(-1)
    end if
    nsvd = min(n,nsvd)
    allocate(ss(n))
    allocate(u(n1,n))
    allocate(vt(n,n2))
    if ( .true. ) then
        lwork = -1
        if ( lwrite ) print *,'calling sgesvd to compute lwork'
        call sgesvd('S','S',n1,n2,cov,nxymax1,ss,u,n1,vt,n,scratch &
        ,lwork,info)
        lwork = nint(scratch(1))
        if ( lwrite ) print *,'sgesvd tells me I need lwork = ' &
        ,lwork
    else
        lwork = 50*max(n1,n2)
        if ( lwrite ) print *,'I guess that lwork = ',lwork, &
        ' is enough'
    end if
    allocate(work(lwork))
    if ( lwrite ) print *,'calling sgesvd'
    call sgesvd('S','S',n1,n2,cov,nxymax1,ss,u,n1,vt,n,work,lwork &
    ,info)
    if ( info /= 0 ) then
        write(0,*) 'svd: error in sgesvd, info = ',info
        call exit(-1)
    end if
    deallocate(work)
    write(0,'(a,f8.2)')'Eigenvalues, time: ',etime(tarray)
    write(0,'(a)')'<table class="realtable" '// &
    'border=0 cellpadding=0 cellspacing=0>'// &
    '<tr><th>#</th><th>eigenvalue</th>'
    write(0,'(a)')'<th>explained variance</th>'// &
    '<th>cumulative</th></tr>'
    t = 0
    do i=1,n
        t = t + ss(i)
    end do
    s = 0
    do i=1,nsvd
        s = s + ss(i)/t
        write(0,'(a,i3,a,g16.5,a,f8.2,a,f8.2,a)') &
        '<tr><td>',i,'</td><td>',ss(i) &
        ,'</td><td align=right>',ss(i)/t*100 &
        ,'%</td><td align=right>',s*100 &
        ,'%</td></tr>'
    end do
    write(0,'(a)')'</table>'

!       get rid of too small modes

    do i=1,nsvd
        if ( ss(i)/t < 1e-5 ) then
            nsvd = i-1
            write(0,'(a,i3)') 'number of SVDs reduced to ',nsvd
            exit
        end if
    end do

!       unpack the left singular vectors

    if ( lwrite ) then
        print *,'allocating svd1(',nxf1,nyf1,nzf1,n,')',nxf1*nyf1 &
        *nzf1*n
        print *,'allocating pattern1(',nxf1,nyf1,nzf1,nperyear, &
        '0:1)',nxf1*nyf1*nzf1*nperyear*2
    end if
    allocate(svd1(nxf1,nyf1,nzf1,n))
    allocate(pattern1(nxf1,nyf1,nzf1,nperyear,0:1))
    svd1 = 3e33
    do i=1,nsvd             ! number of vectors
        do i1=1,n1           ! point of vector
            jx1 = ijk1(1,i1)
            jy1 = ijk1(2,i1)
            jz1 = ijk1(3,i1)
            w = wx1(jx1)*wy1(jy1)*wz1(jz1)
            if ( w > 0 ) then
                svd1(jx1,jy1,jz1,i) = u(i1,i)/sqrt(w)
            else
                svd1(jx1,jy1,jz1,i) = 3e33
            end if
        end do
    end do
    if ( normalization > 0 ) then
        call normsvd1(svd1,nxf1,nyf1,nzf1,nx1,ny1,nz1,n &
        ,normalization)
    end if

!       get corresponding time series

    if ( lwrite ) print *,'allocating pc1(',nperyear,yr1,yr2,nens2 &
    ,nsvd,')',nperyear*(yr2-yr1+1)*(1+nens2)*nsvd
    allocate(pc1(nperyear,yr1:yr2,0:nens2,nsvd))
    pc1 = 3e33
    do i=1,nsvd
    !           construct a pattern file suitable for use in project
        pattern1 = 3e33
        if ( m1 /= 0 ) then
            pattern1(1:nx1,1:ny1,1:nz1,m1,1) = &
            svd1(1:nx1,1:ny1,1:nz1,i)
        else
            pattern1(1:nx1,1:ny1,1:nz1,nperyear,0) = &
            svd1(1:nx1,1:ny1,1:nz1,i)
        end if
        if ( lwrite ) print *,'calling project3'
        call project3(pc1(1,yr1,0,i),nperyear, &
        yr1,yr2,nens1,nens2,xx1,nx1,yy1,ny1,zz1,nz1, &
        field1,pattern1,nxf1,nyf1,nzf1,nperyear,firstyr,lastyr, &
        m1,minfac,lnomissing,(m1 /= 0),anom,lwrite)
    end do
    if ( normalization < 0 ) then
        if ( lwrite ) print *,'calling normsvd2'
        call normsvd2(pc1,svd1,nperyear,yr1,yr2,m1,m2, &
        nens1,nens2,nsvd,nx1,ny1,nz1,nx1,ny1,nz1, &
        normalization,lnomissing,lwrite)
    end if

!       write to file

    if ( lwrite ) print *,'writing output'
    do i=1,nsvd
        write(vars(i),'(a,i2.2)') 'svd',i
        write(lvars(i),'(a,i3)') 'left singular vector number ',i
        ivars(1,i) = nz1
        ivars(2,i) = 99
        if ( normalization >= 0 ) then
            units(i) = '1'
        else
            units(i) = units1(1)
        end if
    end do
    title1 = 'Left singular vactors of SVD of '//trim(infile1)// &
    ' and '//trim(infile2)
    if ( abs(xx1(1)-360-lon1) < abs(xx1(1)-lon1) ) then
        xx1(1:nx1) = xx1(1:nx1) - 360
    elseif (  abs(xx1(1)+360-lon1) < abs(xx1(1)-lon1) ) then
        xx1(1:nx1) = xx1(1:nx1) + 360
    end if
    ii = index(outfile1,'.ctl')
    if ( ii /= 0 ) then
        file = outfile1(1:ii-1)//'.grd'
        call writectl(outfile1,file,nx1,xx1,ny1,yy1,nz1,zz1, &
        &            1,nperyear,1,max(1,m1),3e33,title1,nsvd,vars &
        ,ivars,lvars,units)
        open(1,file=trim(file),access='direct',form='unformatted' &
        ,recl=recfa4*nx1*ny1*nz1)
        do i=1,nsvd
            write(1,rec=i) (((svd1(jx1,jy1,jz1,i), &
            jx1=1,nx1),jy1=1,ny1),jz1=1,nz1)
        end do
        close(1)
    else
        write(0,*) 'svd: sorry, netcdf output not yet ready'
        call exit(-1)
    end if
    do i=1,nsvd
        do iens=nens1,nens2
            if ( nens1 == nens2 ) then
                write(outfile1(ii:),'(a,i2.2,a)') '_',i,'.dat'
            else
                write(outfile1(ii:),'(a,i2.2,a,i2.2,a)') &
                '_',i,'_',iens,'.dat'
            end if
            open(1,file=trim(outfile1))
            write(1,'(a,i2,a)') '# Time series of left singular vactor ',i &
                ,' of SVD of '//trim(infile1)//' and '//trim(infile2)
            if ( nens1 /= nens2 ) then
                write(1,'(a,i2)') '# ensemble member ',iens
            end if
            if ( normalization <= 0 ) then
                units(i) = '1'
            else
                units(i) = units1(1)
            end if
            write(1,'(a,i2.2,3a)') '# PCl',i,' [',trim(units(i)) &
            ,'] left PC1 of SVD'
            call printdatfile(1,pc1(1,yr1,iens,i),nperyear &
            ,nperyear,yr1,yr2)
            close(1)
        end do
    end do

!       unpack the right singular vectors and write to file

    allocate(svd2(nxf2,nyf2,nzf2,n))
    allocate(pattern2(nxf2,nyf2,nzf2,nperyear,0:1))
    svd2 = 3e33
    do i=1,nsvd             ! number of vectors
        do i2=1,n2          ! point of vector
            jx2 = ijk2(1,i2)
            jy2 = ijk2(2,i2)
            jz2 = ijk2(3,i2)
            w = wx2(jx2)*wy2(jy2)*wz2(jz2)
            if ( w > 0 ) then
                svd2(jx2,jy2,jz2,i) = vt(i,i2)/sqrt(w)
            else
                svd2(jx2,jy2,jz2,i) = 3e33
            end if
        end do
    end do
    if ( normalization > 0 ) then
        call normsvd1(svd2,nxf2,nyf2,nzf2,nx2,ny2,nz2,n &
        ,normalization)
    end if

!       get corresponding time series

    allocate(pc2(nperyear,yr1:yr2,0:nens2,nsvd))
    pc2 = 3e33
    if ( m /= 0 ) then
        m = m1 - lag1
        if ( m > nperyear ) m = 1 + mod(m-1,nperyear)
        if ( m <= 0 ) m = nperyear + mod(m,nperyear)
    end if
    do i=1,nsvd
    !           construct a pattern file suitable for use in project
        pattern2 = 3e33
        if ( m1 /= 0 ) then
            pattern2(1:nx2,1:ny2,1:nz2,m,1) = &
            svd2(1:nx2,1:ny2,1:nz2,i)
        else
            pattern2(1:nx2,1:ny2,1:nz2,nperyear,0) = &
            svd2(1:nx2,1:ny2,1:nz2,i)
        end if
        call project3(pc2(1,yr1,0,i),nperyear, &
        yr1,yr2,nens1,nens2,xx2,nx2,yy2,ny2,zz2,nz2, &
        field2,pattern2,nxf2,nyf2,nzf2,nperyear,firstyr,lastyr, &
        m,minfac,lnomissing,(m1 /= 0),anom,lwrite)
    end do
    if ( normalization < 0 ) then
        call normsvd2(pc2,svd2,nperyear,yr1,yr2,m1,m2, &
        nens1,nens2,nsvd,nxf2,nyf2,nzf2,nx2,ny2,nz2, &
        normalization,lnomissing,lwrite)
    end if
    do i=1,nsvd
        write(vars(i),'(a,i2.2)') 'svd',i
        write(lvars(i),'(a,i3)') 'right singular vector number ',i
        ivars(1,i) = nz2
        ivars(2,i) = 99
        if ( normalization >= 0 ) then
            units(i) = '1'
        else
            units(i) = units2(1)
        end if
    end do
    title2 = 'Right singular vactors of SVD of '//trim(infile1)// &
    ' and '//trim(infile2)
    if ( abs(xx2(1)-360-altlon1) < abs(xx2(1)-lon1) ) then
        xx2(1:nx2) = xx2(1:nx2) - 360
    elseif (  abs(xx2(1)+360-lon1) < abs(xx2(1)-lon1) ) then
        xx2(1:nx2) = xx2(1:nx2) + 360
    end if
    if ( m1 /= 0 ) then
        m = m1 - lag1
        i = 1
        call normon(m,i,yr,nperyear)
    else
        m = 0
        yr = 1
    end if
    ii = index(outfile2,'.ctl')
    if ( i /= 0 ) then
        file = outfile2(1:ii-1)//'.grd'
        call writectl(outfile2,file,nx2,xx2,ny2,yy2,nz2,zz2, &
        &            1,nperyear,yr,max(1,m),3e33,title2,nsvd,vars &
        ,ivars,lvars,units)
        open(1,file=trim(file),access='direct',form='unformatted' &
        ,recl=recfa4*nx2*ny2*nz2)
        do i=1,nsvd
            write(1,rec=i) (((svd2(jx2,jy2,jz2,i), &
            jx2=1,nx2),jy2=1,ny2),jz2=1,nz2)
        end do
        close(1)
    else
        write(0,*) 'svd: sorry, netcdf output not yet ready'
        call exit(-1)
    end if
    do i=1,nsvd
        do iens=nens1,nens2
            if ( nens1 == nens2 ) then
                write(outfile2(ii:),'(a,i2.2,a)') '_',i,'.dat'
            else
                write(outfile2(ii:),'(a,i2.2,a,i2.2,a)') &
                '_',i,'_',iens,'.dat'
            end if
            open(1,file=trim(outfile2))
            write(1,'(a,i2,a)') '# Time series of reight singular vactor ',i &
            ,' of SVD of '//trim(infile1)//' and '//trim(infile2)
            if ( nens1 /= nens2 ) then
                write(1,'(a,i2)') '# ensemble member ',iens
            end if
            if ( normalization <= 0 ) then
                units(i) = '1'
            else
                units(i) = units2(1)
            end if
            write(1,'(a,i2.2,3a)') '# PCl',i,' [',trim(units(i)) &
            ,'] right PC1 of SVD'
            call printdatfile(1,pc2(1,yr1,iens,i),nperyear &
            ,nperyear,yr1,yr2)
            close(1)
        end do
    end do
end program svd

subroutine dump060(xx,yy,zz,field,nxf,nyf,nzf,nx,ny,nz &
    ,nperyear,firstyr,lastyr)

!       dumps the field at 0,60N

    implicit none
    include 'params.h'
    integer :: nxf,nyf,nzf,nx,ny,nz,nperyear,firstyr,lastyr
    real :: xx(nx),yy(ny),zz(nz), &
    field(nxf,nyf,nzf,nperyear,firstyr:lastyr)
    integer :: x1,x2,y1,y2,i,j,yr,mo
    real :: lon1,lat1,lon2,lat2,lon1c,lat1c,lon2c,lat2c, &
    data(npermax,yrbeg:yrend)
            
    lon1 = 0
    lat1 = 60
    lon2 = 0
    lat2 = 60
    call getlonwindow(lon1,lon2,x1,x2,xx,nx,lon1c,lon2c, .false. )
    call getlatwindow(lat1,lat2,y1,y2,yy,ny,lat1c,lat2c, .false. )
    if ( lon1c > 1e33 .or. lat1c >= 1e33 ) then
        x1 = 1
        lon1c = xx(1)
        y1 = 1
        lat1c = yy(1)
    end if
    print *,'cutting out longitude ',x1,x2,lon1c,lon2c
    print *,'cutting out latitude  ',y1,y2,lat1c,lat2c
    call makeabsent (data,npermax,yrbeg,yrend)
    do yr=max(yrbeg,firstyr),min(lastyr,yrend)
        do mo=1,nperyear
            data(mo,yr) = field(x1,y1,1,mo,yr)
        end do
    end do
    call printdatfile(6,data,npermax,nperyear,yrbeg,yrend)
end subroutine dump060

subroutine processoptions(field,nxf,nyf,nzf,mpermax,firstyr &
    ,lastyr,mens,nx,ny,nz,nperyear,i12)
    implicit none
    include 'params.h'
    include 'getopts.inc'
    integer :: nxf,nyf,nzf,mpermax,nx,ny,nz,nperyear,firstyr, &
        lastyr,mens,i12
    real :: field(nxf,nyf,nzf,mpermax,firstyr:lastyr,0:mens)
    integer :: jx,jy,jz,iens,i,j
    real :: fxy(npermax,yrbeg:yrend)

    do iens=nens1,nens2
        do jz=1,nz
            do jy=1,ny
                call keepalive(jy+(jz-1)*ny,ny*nz)
                do jx=1,nx
                
                !                       create 1-D series from field
                
                    call makeabsent(fxy,npermax,yrbeg,yrend)
                    do i=yr1,yr2
                        do j=1,nperyear
                            if ( field(jx,jy,jz,j,i,iens) < 1e30 ) &
                            then
                                fxy(j,i) = field(jx,jy,jz,j,i,iens)
                            end if
                        end do
                    end do
                
                !                       take monthly anomalies
                
                    if ( i12 == 1 .and. mdiff > 0 ) then
                        call mdiffit(fxy,npermax,nperyear,yrbeg &
                        ,yrend,mdiff)
                    end if
                    if ( i12 == 2 .and. mdiff2 > 0 ) then
                        call mdiffit(fxy,npermax,nperyear,yrbeg &
                        ,yrend,mdiff2)
                    end if
                
                !                       sum
                
                    if ( i12 == 1 .and. lsum > 1 ) then
                        call sumit(fxy,npermax,nperyear,yrbeg,yrend &
                        ,lsum,'v')
                    end if
                    if ( i12 == 2 .and. lsum2 > 1 ) then
                        call sumit(fxy,npermax,nperyear,yrbeg,yrend &
                        ,lsum2,'v')
                    end if
                
                !                       log, sqrt
                
                    if ( i12 == 1 .and. logscale .or. &
                    i12 == 2 .and. logfield ) then
                        call takelog(fxy,npermax,nperyear,yrbeg &
                        ,yrend)
                    end if
                    if ( i12 == 1 .and. sqrtscale .or. &
                    i12 == 2 .and. logfield ) then
                        call takesqrt(fxy,npermax,nperyear,yrbeg &
                        ,yrend)
                    end if
                
                !                       detrend
                
                    if ( ldetrend ) then
                        if ( lwrite ) print *,'Detrending field'
                        if ( lag1 == 0 .and. lag2 == 0 .or. m1 == 0 &
                         .or. lsel == nperyear ) then
                            call detrend(fxy,npermax,nperyear, &
                            yrbeg,yrend,yr1,yr2,m1,m2,lsel)
                        else
                            call detrend(fxy,npermax,nperyear, &
                            yrbeg,yrend,yr1,yr2,1,12,lsel)
                        end if
                    end if
                
                !                       differentiate
                
                    if ( ndiff > 0 ) then
                        if ( lwrite ) print *,'Taking differences'
                        call diffit(fxy,npermax,nperyear,yrbeg &
                        ,yrend,ndiff)
                    end if
                
                !                       normalize to s.d.
                
                    if ( lnormsd ) then
                        call normsd(fxy,npermax,nperyear,yrbeg,yrend &
                        ,yr1,yr2)
                    else
                        call anomal(fxy,npermax,nperyear,yrbeg &
                        ,yrend,yr1,yr2)
                    end if
                
                !                       copy back to field
                
                    do i=yr1,yr2
                        do j=1,nperyear
                            field(jx,jy,jz,j,i,iens) = fxy(j,i)
                        end do
                    end do
                end do       ! iens
            end do           ! jz
        end do               ! jy
    end do                   ! jx
end subroutine processoptions
