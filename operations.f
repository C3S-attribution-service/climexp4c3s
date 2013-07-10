*  #[ takelog:
        subroutine takelog(data,npermax,nperyear,yrbeg,yrend)
*
*       take logarithm of time series
*       replaces zero's by 1/5 of minimum value
*
        implicit none
        integer yrbeg,yrend,npermax,nperyear
        real data(npermax,yrbeg:yrend)
        integer i,j,ii,jj
        integer,save :: nfirstzero
        real minval
        data nfirstzero /0/
*
        minval = 3e33
        do i=yrbeg,yrend
            do j=1,nperyear
                if ( data(j,i).gt.0 ) minval = min(minval,data(j,i))
            enddo
        enddo
        do i=yrbeg,yrend
            do j=1,nperyear
                if ( data(j,i).lt.1e33 ) then
                    if ( data(j,i).eq.0 ) then
                        if ( minval.gt.1e33 ) then
*                           this happens when I only have 0
                            data(j,i) = 3e33
                        else
                            if ( nfirstzero.lt.3 ) then
                                print '(a,g14.6)',
     +                               '# takelog: replaced zero by '
     +                               ,minval/2
                                nfirstzero = nfirstzero + 1
                            endif
                            data(j,i) = minval/2
                        endif
                    endif
                    if ( data(j,i).gt.0 ) then
                        data(j,i) = log10(data(j,i))
                    else
                        data(j,i) = 3e33
                    endif
                endif
            enddo
        enddo
*
        end
*  #] takelog:
*  #[ takesqrt:
        subroutine takesqrt(data,npermax,nperyear,yrbeg,yrend)
*
*       take sqrt of time series
*
        implicit none
        integer yrbeg,yrend,npermax,nperyear
        real data(npermax,yrbeg:yrend)
        integer i,j
*
        do i=yrbeg,yrend
            do j=1,nperyear
                if ( data(j,i).lt.1e33 ) then
                    if ( data(j,i).gt.0 ) then
                        data(j,i) = sqrt(data(j,i))
                    else
                        data(j,i) = 3e33
                    endif
                endif
            enddo
        enddo
*
        end
*  #] takesqrt:
*  #[ takeexp:
        subroutine takeexp(data,npermax,nperyear,yrbeg,yrend)
*
*       take exp of time series
*
        implicit none
        integer yrbeg,yrend,npermax,nperyear
        real data(npermax,yrbeg:yrend)
        integer i,j
        integer,save :: n
        data n /0/
*
        do i=yrbeg,yrend
            do j=1,nperyear
                if ( data(j,i).lt.1e33 ) then
                    if ( data(j,i).lt.30 ) then
                        data(j,i) = 10.**data(j,i)
                    else
                        n = n + 1
                        if ( n.lt.5 ) then
                            write(0,*) 'takeexp: unlikely value ',
     +                           data(j,i),' set to undefined'
                        endif
                        data(j,i) = 3e33
                    endif
                endif
            enddo
        enddo
*
        end
*  #] takeexp:
*  #[ takesquare:
        subroutine takesquare(data,npermax,nperyear,yrbeg,yrend)
*
*       take square of time series
*
        implicit none
        integer yrbeg,yrend,npermax,nperyear
        real data(npermax,yrbeg:yrend)
        integer i,j
*
        do i=yrbeg,yrend
            do j=1,nperyear
                if ( data(j,i).lt.1e33 ) then
                    data(j,i) = data(j,i)**2
                endif
            enddo
        enddo
*
        end
*  #] takesquare:
*  #[ takeinv:
        subroutine takeinv(data,npermax,nperyear,yrbeg,yrend)
*
*       take inverse of time series
*
        implicit none
        integer yrbeg,yrend,npermax,nperyear
        real data(npermax,yrbeg:yrend)
        integer i,j
*
        do i=yrbeg,yrend
            do j=1,nperyear
                if ( data(j,i).lt.1e33 ) then
                    if ( abs(data(j,i)).gt.1e-33 ) then
                        data(j,i) = 1/data(j,i)
                    else
                        data(j,i) = 3e33
                    end if
                endif
            enddo
        enddo
*
        end
*  #] takeinv:
*  #[ changesign:
        subroutine changesign(data,npermax,nperyear,yrbeg,yrend)
*
*       change sign of time series
*
        implicit none
        integer yrbeg,yrend,npermax,nperyear
        real data(npermax,yrbeg:yrend)
        integer i,j
*
        do i=yrbeg,yrend
            do j=1,nperyear
                if ( data(j,i).lt.1e33 ) then
                    data(j,i) = -data(j,i)
                endif
            enddo
        enddo
*
        end
*  #] changesign:
*  #[ takefieldlog:
        subroutine takefieldlog(data,nxmax,nymax,npermax,nperyear,nx,ny
     +        ,yrbeg,yrend)
*
*       take logarithm of field
*       replaces zero's by 1/5 of minimum value
*
        implicit none
        integer nxmax,nymax,yrbeg,yrend,npermax,nperyear,nx,ny
        real data(nxmax,nymax,npermax,yrbeg:yrend)
        integer i,j,ix,iy
        logical lfirstzero
        real minval
*
        minval = 3e33
        do i=yrbeg,yrend
            do j=1,nperyear
                do iy=1,ny
                    do ix=1,nx
                        if ( data(ix,iy,j,i).gt.0 ) minval = min(minval
     +                        ,data(ix,iy,j,i))
                    enddo
                enddo
            enddo
        enddo
        lfirstzero = .TRUE.
        do i=yrbeg,yrend
            do j=1,nperyear
                do iy=1,ny
                    do ix=1,nx
                        if ( data(ix,iy,j,i).lt.1e33 ) then
                            if ( data(ix,iy,j,i).eq.0 ) then
                                if ( lfirstzero ) then
                                    write(0,'(a,f12.2,a)')
     +                                   'takefieldlog: replaced '//
     +                                   'zero by ',minval/2,'<br>'
                                    lfirstzero = .FALSE.
                                endif
                                data(ix,iy,j,i) = minval/2
                            endif
                            if ( data(ix,iy,j,i).gt.0 ) then
                                data(ix,iy,j,i) = log10(data(ix,iy,j,i))
                            else
                                data(ix,iy,j,i) = 3e33
                            endif
                        endif
                    enddo
                enddo
            enddo
        enddo
*
        end
*  #] takefieldlog:

