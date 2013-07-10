        program get_area
*
*       get the average from a GrADS data file of a certain area
*
        implicit none
#include "recfac.h"
        integer nxmax,nymax
        parameter (nxmax=500,nymax=250)
        integer i,j,k,n,nx,ny,nt,xl,x1,x2,y1,y2,yrbegin,year,month,iskip
     +        ,nline,sgn,ii,jj
        real field(nxmax,nymax),minfac,ave(12)
        logical lwrite
        character string*255,dipole*2,file*255
        integer iargc

        minfac = 0.1
        lwrite = .FALSE.
        if ( iargc().lt.9 ) then
            print *,'usage: get_area file nx ny nt yrbegin x1 x2 y1 y2 '
     +            //'[minfac r]'
        endif
        
        call getarg(2,string)
        read(string,*) nx
        call getarg(3,string)
        read(string,*) ny
        call getarg(4,string)
        read(string,*) nt
        call getarg(5,string)
        read(string,*) yrbegin

        call getarg(1,file)
        print *,'using file ',file(1:index(file,' ')-1)
        open(1,file=file,form='unformatted',access='direct',recl
     +        =recfa4*nx*ny,status='old')
        if ( lwrite ) print *,'recl = ',recfa4,'*',nx,'*',ny,'=',recfa4
     +        *nx*ny
*
        call getarg(6,string)
        read(string,*) x1
        call getarg(7,string)
        read(string,*) x2
        call getarg(8,string)
        read(string,*) y1
        call getarg(9,string)
        read(string,*) y2
        print *,'cutting out region x=',x1,x2,', y=',y1,y2

        minfac = 0.1
        iskip = 0
        dipole = '  '
        nline = 3
        do i=10,iargc()
            if ( iskip.gt.0 ) then
                iskip = iskip - 1
            else
                call getarg(i,string)
                if ( string(1:6).eq.'minfac' ) then
                    call getarg(i+1,string)
                    read(string,*) minfac
                    iskip = 1
                    print *,'using minimal fraction of valid points '
     +                    ,minfac
                    nline = nline + 1
                endif
                if ( string(1:6).eq.'minnum' ) then
                    call getarg(i+1,string)
                    read(string,*) minnum
                    iskip = 1
                    print *,'using minimal number of valid points '
     +                    ,minnum
                    nline = nline + 1
                endif
                if ( string(1:6).eq.'dipole' ) then
                    call getarg(i+1,dipole)
                    iskip = 1
                    print *,'using dipole ',dipole
                    nline = nline + 1
                endif
            endif
        enddo
        do i=nline,5
            print *,' '
        enddo

        month = 0
        year = yrbegin
        do k=1,nt
            if ( lwrite ) print *,'Reading record ',k
            read(1,rec=k) ((field(i,j),i=1,nx),j=1,ny)
            month = month + 1
            if ( month.ge.13 ) then
                print '(i5,12f8.2)',year,ave
                month = month - 12
                year = year + 1
                do i=1,12
                    ave(i) = 999.9
                enddo
            endif
            ave(month) = 0
            n = 0
            do j=y1,y2
                do i=x1,x2
                    if ( i.le.0 ) then
                        ii = i + (-i/nx+1)*nx
                    elseif ( i.gt.nx ) then
                        ii = 1 + mod(i-1,nx)
                    else
                        ii = i
                    endif
                    if ( index(file,'jones').ne.0 .and. field(ii,j).lt.
     +                    32 .or. field(ii,j).lt.1e33 ) then
                        if (  dipole.eq.'  ' .or.
     +                        dipole.eq.'ll' .and. i.lt.(x1+x2)/2. .or.
     +                        dipole.eq.'rr' .and. i.gt.(x1+x2)/2. .or.
     +                        dipole.eq.'tt' .and. j.lt.(y1+y2)/2. .or.
     +                        dipole.eq.'bb' .and. j.gt.(y1+y2)/2. .or.
     +                        (dipole.eq.'bl' .or. dipole.eq.'lb').and.
     +                        (i-x1)*(y1-y2).lt.(x2-x1)*(j-y2) .or.
     +                        (dipole.eq.'tr' .or. dipole.eq.'rt').and.
     +                        (i-x1)*(y1-y2).gt.(x2-x1)*(j-y2) .or.
     +                        (dipole.eq.'br' .or. dipole.eq.'rb').and.
     +                        (i-x1)*(y2-y1).gt.(x2-x1)*(j-y1) .or.
     +                        (dipole.eq.'tl' .or. dipole.eq.'lt').and.
     +                        (i-x1)*(y2-y1).lt.(x2-x1)*(j-y1) ) then
                            sgn = +1
                        else
                            sgn = -1
                        endif
                        if ( lwrite ) print *,'adding field(',ii,j
     +                        ,') = ',field(ii,j)
                        ave(month) = ave(month) + sgn*field(ii,j)
                        n = n + 1
                    elseif ( abs(field(ii,j)-32.767).gt.1e-4 .and. 
     +                        field(ii,j).lt.1e33 ) then
                        print *,'found large value ',field(ii,j)
     +                        ,' assumed invalid'
                    endif
                enddo
            enddo
            if ( n.lt.minfac*(x2-x1+1)*(y2-y1+1) .or. n.lt.minnum 
     +          ) then
                ave(month) = -999.9
            else
                ave(month) = ave(month)/n
            endif
        enddo
        print '(i5,12f8.2)',year,ave

        end
