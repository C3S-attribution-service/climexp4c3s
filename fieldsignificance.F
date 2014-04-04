        program fieldsignficance
!
!       given the output of correlatefield, compute the field
!       signficance
!
        implicit none
#include "netcdf.inc"
#include "params.h"
#include "recfac.h"
#include "getopts.inc"
        integer nvarmax,ncuts
        parameter (nvarmax=30,ncuts=9)
        integer i,j,itime,ivar,iivar,irec,ix,ix1,ixm1,iy,iz,nd,np,x1,x2
     +       ,y1,y2,ncut,icut
        integer mens1,mens,ncid,nx,ny,nz,nt,nperyear,firstmo,firstyr
     +       ,nvars,endian,ivars(6,nvarmax),nxx,nyy
        integer,allocatable :: npairs(:)
        real s,s0,p,w,d,pi,deg2rad,q,qq,qcut,qcut1,dx,dy,area
        real,allocatable :: prob(:,:,:),corr(:,:,:),pairs(:,:,:)
     +       ,autocor(:,:),invprob(:,:,:)
        real xx(nxmax),yy(nymax),zz(nzmax),undef,wx(nxmax),wy(nymax)
     +       ,cuts(ncuts),pcuts(6)
        character infile*255,datfile*255,string*20,lz(3)*20,ltime*100,
     +       title*1000,history*10000,vars(nvarmax)*40,
     +       lvars(nvarmax)*200,svars(nvarmax)*200,units(nvarmax)*100,
     +       cell_methods(nvarmax)*128,all*6
        logical xrev,yrev,xwrap
        integer iargc
        real df,pcut
        common /qcom/ df,pcut
        lwrite = .false.
        data cuts /0.1,0.2,0.5,1.,2.,5.,10.,20.,50./
        data pcuts /0.05,0.02,0.01,0.005,0.002,0.001/ ! the return times 20,50,100,200,500,1000 years
!
!       process options
!
        call getarg(3,string)
        if ( string.eq.' ' ) then
            write(0,*) 'usage: fieldsignificance file timestep pcut% ',
     +           'lon1 x1 lon2 x2 lat1 y1 lat2 y2'
            call abort
        end if
        all = string
        if ( string(1:3).eq.'all' ) then
            pcut = -1
        else
            read(string,*,err=901) pcut
            pcut = pcut/100
            if ( pcut.eq.1 ) pcut = 0.1 ! choose an arbitrary cut
        end if
        if ( lwrite ) print *,'read pcut = ',pcut
        call getarg(2,string)
        read(string,*,err=902) itime
        if ( lwrite ) print *,'read itime = ',itime
        call getarg(1,infile)
!
!       read data
!
        call getmetadata(infile,mens1,mens,ncid,datfile,nxmax,nx
     +       ,xx,nymax,ny,yy,nzmax,nz,zz,lz,nt,nperyear,firstyr,firstmo
     +       ,ltime,undef,endian,title,history,nvarmax,nvars,vars,ivars
     +       ,lvars,svars,units,cell_methods,lwrite)
        if ( mens1.gt.0 ) then
            write(0,*) 'fieldsignficance: cannot handle ensembles',mens1
            call abort
        end if
        do ivar=1,nvars
            if ( vars(ivar)(1:4).eq.'prob' .or.
     +           index(vars(ivar),'_rt').ne.0 ) then
                if ( all.eq.'all_lo' ) then
                    if ( vars(ivar)(1:3).eq.'lo_' ) goto 101
                else if ( all.eq.'all_hi' ) then
                    if ( vars(ivar)(1:3).eq.'hi_' ) goto 101
                else
                    goto 101
                end if
            end if
        end do
        write(*,*) 'fieldsignficance: error: cannot find variable '//
     +       'prob* or *_rt*'
        goto 999
 101    continue
        if ( vars(ivar)(1:4).eq.'prob' ) then
            do iivar=1,nvars
                if ( vars(iivar)(1:3).eq.'cor' .or.
     +               vars(iivar)(1:3).eq.'dif' .or.
     +               vars(iivar)(1:6).eq.'reldif' ) goto 102
            end do
            write(*,*) 'fieldsignficance: error: cannot find variable'//
     +           ' cor*'
            goto 999
 102        continue
        else
            iivar = -999
        end if
        if ( itime.le.0 .or. itime.gt.nt ) then
            write(0,*) 'fieldignificance: error: itime < 0 or itime > ',
     +           'nt: ',itime,nt
            call abort
        end if
        call getopts(4,iargc(),nperyear,yrbeg,yrend,.true.,mens1,mens)
        allocate(prob(nx,ny,nz))
        allocate(invprob(nx,ny,nz))
        allocate(corr(nx,ny,nz))
        if ( ncid.eq.-1 ) then
            open(1,file=trim(datfile),form='unformatted',
     +           access='direct',recl=recfa4*nx*ny*nz,status='old')
            irec = ivar + (itime-1)*nvars
            read(1,rec=irec)
     +           (((prob(ix,iy,iz),ix=1,nx),iy=1,ny),iz=1,nz)
            if ( iivar.ne.-999 ) then
                irec = iivar + (itime-1)*nvars
                read(1,rec=irec)
     +               (((corr(ix,iy,iz),ix=1,nx),iy=1,ny),iz=1,nz)
                close(1)
            else
                corr = 1
            end if
        else
            call readonencfield(ncid,prob,nx,ny,nz,itime,ivars(1,ivar))
            if ( iivar.ne.-999 ) then
                call readonencfield(ncid,corr,nx,ny,nz,itime,
     +               ivars(1,iivar))
            else
                corr = 1
            end if
        end if
        if (  vars(ivar)(1:4).ne.'prob' ) then
            forall(ix=1:nx,iy=1:ny,iz=1:nz, prob(ix,iy,iz).lt.1e33) 
     +           prob(ix,iy,iz) = 1/prob(ix,iy,iz)
        end if
!
!       compute area with prob < pcut
!
        call getxyprop(xx,nx,yy,ny,xrev,yrev,xwrap)
        call getlatlonwindow(lat1,lat2,lon1,lon2,xx,nx,xwrap,avex,yy,ny
     +        ,avey,x1,x2,y1,y2,lwrite)
        call getweights('x',xx,wx,nx,xwrap,lwrite)
        call getweights('y',yy,wy,ny,.FALSE.,lwrite)
        if ( pcut.eq.-1 ) then
            ncut = 6
        else
            ncut = 1
        end if
        do icut=1,ncut
            if ( ncut.gt.1 ) then
                pcut = pcuts(icut)
            end if
            s = 0
            s0 = 0
            do iz=1,nz
                do iy=y1,y2
                    do ix1=x1,x2
                        ix = ix1
                        if ( xwrap .and. ix.lt.1 ) ix = ix + nx
                        if ( xwrap .and. ix.gt.nx ) ix = ix - nx
                        if ( abs(prob(ix,iy,iz)/undef-1).gt.1e-6 ) then
                            w = wx(ix)*wy(iy)
                            s0 = s0 + w
                            if ( prob(ix,iy,iz).lt.pcut ) then
                                s = s + w
                            end if
                        end if
                    end do
                end do
            end do
            area = s0
            if ( s.eq.0 ) then
                write(0,'(2a)') 'Cannot find enough data to compute ',
     +               'field significance '
                goto 999
            end if
            p = s/s0
            if ( nz.gt.1 ) then
                write(0,'(a)') '(Giving equal weight to each layer.)'
                goto 904        ! cannot handle vertical correlations yet
            end if
            if ( .false. .and. p.le.pcut ) then
!               save time
                write(0,'(a)')
     +               'Statistically, this is indistinguishable from '
     +               ,'pure chance (p&gt;0.5).'
                goto 999
            end if
!           not enough data...
            if ( x1.eq.x2 .or. y1.eq.y2 ) goto 999
            if ( icut.eq.1 ) then
!
!               compute number of degrees of freedom
!               first map the probabilities to th real axis using the sign of the
!               correlation. I think that we should emphasise the spatial scale of 
!               the small p-values more than the scale of the big ones so I changed
!               it to log(1/p) on 19-06-2013
!               
                do iz=1,nz
                    do iy=1,ny
                        do ix=1,nx
                            if ( prob(ix,iy,iz).lt.1e33 ) then
                                if ( corr(ix,iy,iz).lt.0 ) then
                                    invprob(ix,iy,iz) = 
     +                                   -log(1/prob(ix,iy,iz))
                                else
                                    invprob(ix,iy,iz) =
     +                                   log(1/prob(ix,iy,iz))
                                end if
                            end if
                        end do
                    end do
                end do
!
!               next compute the parameters for getdecor
!
                s = 0
                do ix1=x1+1,x2
                    ix = ix1
                    ixm1 = ix-1
                    if ( xwrap ) then
                        if ( ix.le.0 ) ix = ix + nx
                        if ( ix.gt.nx ) ix = ix - nx
                        if ( ixm1.le.0 ) ixm1 = ixm1 + nx
                        if ( ixm1.gt.nx ) ixm1 = ixm1 - nx
                    end if
                    
                    d = min(abs(xx(ix) - xx(ixm1)),
     +                   abs(xx(ix) - xx(ixm1) + 360),
     +                   abs(xx(ix) - xx(ixm1) - 360))
                    s = s + d
                end do
                do iy=y1+1,y2
                    s = s + abs(yy(iy) - yy(iy-1))
                end do
                d = s/(x2-x1+y2-y1) ! average bin size
                nxx = min(20,x2-x1+1)
                nyy = min(20,y2-y1+1)
                nd = min(50,x2-x1+1,y2-y1+1) ! number of bins for autocor function
                np = (nxx*nyy)*(min(nxx,2*nd+1)*min(nyy,2*nd+1))/2+1
                if ( lwrite ) then
                    print *,'s,n = ',s,(x2-x1)+(y2-y1)
                    print *,'d = ',d
                    print *,'nxx = ',nxx
                    print *,'nyy = ',nyy
                    print *,'nd  = ',nd
                    print *,'np  = ',np
                    print *,'size = ',np*(nd+1)*2,2.*np*(nd+1)
                end if
!           
!               and call getcdecor
!
                allocate(pairs(np,0:nd,2))
                allocate(npairs(0:nd))
                allocate(autocor(0:nd,2))
                call getdecor(decor,invprob,xx,nx,nxx,xwrap,yy,ny,nyy,
     +               x1,x2,y1,y2,d,pairs,npairs,autocor,nd,np,lwrite)
                pi = 4*atan(1.)
                deg2rad = atan(1.)/45.
!               decorrelation disk (assuming equal lat/lon decorrelation lengths
!               etc)
                d = pi*decor**2
!               degrees of freedom.  This is not correct if the grid points with
!               valid data are scattered rather than contiguous.
                df = max(1.,area/d)
                if ( .false. .and. df.lt.2 ) then
                    write(0,'(a)') 'This is too small to compute a ',
     +                   'field significance.'
                    goto 999
                end if
            end if              ! icut.eq.1
!           
!           estimate field significance
!           
            qcut1 = 0
            do i=1,ncuts
                qcut = cuts(i)/100.
                qq = q(qcut)
                if ( lwrite ) print *,'pcut,qcut,df,q = ',pcut,qcut,df,
     +               '=>',qq
                if ( p.gt.qq ) then
                    exit
                end if
                qcut1 = qcut
            end do
            if ( i.gt.ncuts ) then
                qcut = 1
            end if
            if ( ncut.eq.1 ) then ! normal map
                if ( qcut.le.0.25 ) then
                    write(0,'(a)') 'Statistically, there is'
                    if ( qcut.le.0.001 ) then
                        write(0,'(a)') 'almost certainly'
                    else if ( qcut.le.0.01 ) then
                        write(0,'(a)') 'very likely'
                    else if ( qcut.le.0.1 ) then
                        write(0,'(a)') 'likely'
                    else
                        write(0,'(a)') 'maybe'
                    end if
                    write(0,'(2a)')'a significant connection in the map'
                    if ( lweb ) then
                        write(0,'(a,f5.1,a)') '(p<sub>field</sub> &lt; '
     +                       ,qcut*100,'%).'
                    else
                        write(0,'(a,f5.1,a)')'(pfield < ',qcut*100,'%).'
                    end if
                else
                    write(0,'(2a)') 'Statistically, the map is ',
     +                   'indistinguishable from random noise'
                    write(0,'(a,f5.1,a)') '(pfield > ',100*qcut1,'%).'
                end if
                i = index(infile,'/',.true.)
                j = index(infile,'.',.true.)
                if ( j.eq.0 ) j = 1 + len_trim(infile)
                if ( lweb ) then
                    write(0,'(4a,i4.4,a)') '<a href="javascript:',
     +                   'hidden_info_switch(''hidden_',infile(i+1:j-1),
     +                   '_',itime,''');">Details...</a></div>'
                    write(0,'(3a,i4.4,a)') '<div id="hidden_',
     +                   infile(i+1:j-1),'_',itime,
     +                   '" style="display: none;">'
                end if
                write(0,'(a,f6.2,a,f6.2,a)')
     +               'The fraction of the map with p&lt;',100*pcut,
     +               '% is ',100*p,'%.'
            else
                write(0,'(a,i6,a,f6.2,a)')
     +               'The fraction of the map with return time &gt;',
     +               nint(1/pcut),' yr is ',100*p,'%.'
            end if
            if ( icut.eq.1 ) then
                write(0,'(a,f5.1,a,f4.0,a,i5,a)')
     +               'With an estimated decorrelation scale of ',decor
     +               ,'&deg; and (',sqrt(area)
     +               ,'&deg;)&sup2; with data there are about',
     +               2*nint(df),' degrees of freedom in the map.'
            end if
            write(0,'(a,f5.1,a,f5.1,a)')
     +           'This gives a field significance of ',100*qcut1,
     +           '% &lt; pfield &lt; ',100*qcut,'%.'
!!!         write(0,'(a,f4.1,a)') '(with an error margin of about ',
!!!         +       20*qcut,'%).'
            if ( lweb ) then
                if ( ncut.eq.1 ) then
                    write(0,'(a)') '</div>'
                else
                    write(0,'(a)') '<br>'
                end if
            end if
        end do ! icut
!
!       error messages
!
        goto 999
 901    write(0,*) 'fieldsignficance: error reading pcut from '
     +       ,trim(string)
        call abort
 902    write(0,*) 'fieldsignficance: error reading itime from '
     +       ,trim(string)
        call abort
 904    write(0,*) 'Cannot compute field significance for this grid.'
 999    continue
        end

        subroutine getdecor(decor,field,xx,nx,nxx,xwrap,yy,ny,nyy,
     +       x1,x2,y1,y2,d,pairs,npairs,autocor,nd,np,lwrite)
!
!       compute one spatial decorrelation of (part of) a 2D field
!
        implicit none
        integer nx,nxx,ny,nyy,nd,np,x1,x2,y1,y2,npairs(0:nd)
        real decor,field(nx,ny),xx(nx),yy(ny),d,pairs(np,0:nd,2)
     +       ,autocor(0:nd,2)
        logical xwrap,lwrite
        integer i,j,ii,jj,ip,id,idd,ni,nii,i1,j1,ii1,jj1
        integer itoobig,ntoobig,jtoobig,mtoobig
        real dx,deg2rad,lon1,lon2,lat1,lat2,sxx,syy,sxy,e
        external findx
        itoobig = 0
        jtoobig = 0
        ntoobig = 1
        mtoobig = 1
!
!       get pairs for autocorrelation
!
        deg2rad = atan(1.)/45.
        npairs = 0
        do j1=1,nyy
            j = y1 + nint(real(j1-1)*(y2-y1)/(nyy-1))
            if ( j.le.0  ) j = j + ny
            if ( j.gt.ny ) j = j - ny
            if ( lwrite ) print *,'j = ',j
            ni = int(1/cos(yy(j)*deg2rad))
            do i1=1,nxx
                i = x1 + nint(real(i1-1)*(x2-x1)/(nxx-1))
                if ( i.le.0  ) i = i + nx
                if ( i.gt.nx ) i = i - nx
                if ( field(i,j).gt.1e33 ) cycle
                if ( abs(field(i,j)).gt.1e17 ) then
                    itoobig = itoobig + 1
                    if ( itoobig.ge.ntoobig ) then
                        ntoobig = 2*ntoobig
                        print '(2a,2i5,a,g20.4)',
     +                       '# Found too large number ',
     +                       'in field(',i,j,' = ',field(i,j)
                    end if
                    cycle
                end if
                if ( lwrite ) print *,' i = ',i
                do jj1=max(j1,y1),y2
                    jj = jj1
                    if ( jj.le.0  ) jj = jj + ny
                    if ( jj.gt.ny ) jj = jj - ny
                    if ( abs(yy(j)-yy(jj)).gt.(nd+1)*d ) cycle
                    nii = int(1/cos(yy(jj)*deg2rad))
                    do ii1=max(i1,x1),x2,nii
                        ii = ii1
                        if ( ii.le.0  ) ii = ii + nx
                        if ( ii.gt.nx ) ii = ii - nx
                        if ( field(ii,jj) .gt.1e33 ) cycle
                        if ( abs(field(ii,jj)).gt.1e17 ) then
                            jtoobig = jtoobig + 1
                            if ( jtoobig.ge.mtoobig ) then
                                mtoobig = 2*mtoobig
                                print '(2a,2i5,a,g20.4)',
     +                               '# Found too large number 2 ',
     +                               'in field(',ii,jj,' = ',
     +                               field(ii,jj)
                            end if
                            cycle
                        end if
                        lon1 = xx(i)*deg2rad
                        lat1 = yy(j)*deg2rad
                        lon2 = xx(ii)*deg2rad
                        lat2 = yy(jj)*deg2rad
!                       compute great circle distance
!                       http://en.wikipedia.org/wiki/Great-circle_distance
                        dx = atan2( sqrt( (cos(lat2)*sin((lon2-lon1)))
     +                       **2+ (cos(lat1)*sin(lat2)- sin(lat1)
     +                       *cos(lat2)*cos(lon2-lon1))**2 )
     +                       , sin(lat1)*sin(lat2)
     +                       + cos(lat1)*cos(lat2)*cos(lon2-lon1) )
!                       bin
                        id = nint(dx/deg2rad/d)
                        if ( id.lt.0 .or. id.gt.nd ) cycle
!                       add to list
                        if ( npairs(id).ge.np ) then
!!!                            write(0,*) 'too many pairs'
                            cycle
                        endif
                        npairs(id) = npairs(id) + 1
                        pairs(npairs(id),id,1) = field(i,j)
                        pairs(npairs(id),id,2) = field(ii,jj)
                        if ( .false. .and. lwrite ) print
     +                       '(2f7.1,f7.3,2f7.1,f7.3,i3,i4)',xx(i),yy(j)
     +                       ,field(i,j),xx(ii),yy(jj),field(ii,jj),id
     +                       ,npairs(id)
                    enddo
                enddo
            enddo
        enddo
!
!       compute autocorrelation relative to (1,1) and error estimate
!
        autocor = 0
        do id=0,nd
            if ( npairs(id).lt.4 ) then
                autocor(id,1:2) = 3e33
                cycle
            endif
            sxx = 0
            syy = 0
            sxy = 0
            do i=1,npairs(id)
                sxx = sxx + (1-pairs(i,id,1))**2
                syy = syy + (1-pairs(i,id,2))**2
                sxy = sxy + (1-pairs(i,id,1))*(1-pairs(i,id,2))
            end do
            if ( sxx.eq.0 .or. syy.eq.0 ) then
                autocor(id,1) = 3e33
                autocor(id,2) = 3e33
            else
                autocor(id,1) = sxy/(sqrt(sxx*syy))
                autocor(id,2) = autocor(id,1)/sqrt(real(npairs(id))) ! roughly
            end if
        enddo
        if ( lwrite ) then
            print *,'autocor = '
            do id=0,nd
                print *,id,autocor(id,1),autocor(id,2),npairs(id)
            enddo
        endif
!
!       determine distance where it goes through 1/e
!
        e = exp(-1.)
        decor = nd*d            ! default
        do id=0,nd-1
            if ( autocor(id,1).lt.1e33 ) then
                do idd=id,nd
                    if ( autocor(idd,1).lt.1e33 ) then
                        if ( autocor(id,1).gt.e .and. 
     +                       autocor(idd,1).le.e ) then
                            decor = ((autocor(id,1)-e)*idd 
     +                           + (e-autocor(idd,1))*id)
     +                           /(autocor(id,1) - autocor(idd,1))
                            decor = decor*d
                            goto 999
                        end if
                    end if
                end do
            endif
        enddo
 999    continue
!
        end

        real function q(qcut)
!
!       computes q(qcut;df,pcut) in fraction
!       see Telecon/log for details
!
        implicit none
        real qcut
        logical lwrite
        real df,pcut
        common /qcom/ df,pcut
        real aaa,bbb,ccc,aa,bb,a,qq
        lwrite = .false.

        if ( lwrite ) print *,'q: qcut,pcut,df = ',qcut,pcut,df
        aaa = 0.158147
        bbb = -0.0565358
        ccc = 0.00218777
!       note that the fits were done for qcut in percentage points.
        aa = aaa + bbb*log(100*qcut) + ccc*log(100*qcut)**2

        aaa = 0.584773
        bbb = 0.052409
        ccc = 0.00581538
        bb = aaa + bbb*log(100*qcut) + ccc*log(100*qcut)**2
        if ( lwrite ) print *,'aa,bb = ',aa,bb
        
!       again, I used percentages in the fit
        a = aa*(100*pcut)**bb

        qq = a/sqrt(df)
        if ( lwrite ) print *,'a,qq = ',a,qq

        q = pcut + (1-pcut)*(1d0-exp(-dble(qq)))
        if ( lwrite ) print *,'q: q = ',q

        end


        subroutine readonencfield(ncid,data,nx,ny,nz,itime,jvars)
!
!       read one XYZ slice from var described by jvars at time itime
!       from netcdf file open at ncid, assuming all dimensions are filled
!
        implicit none
#include "netcdf.inc"
#include "getopts.inc"
        integer ncid,nx,ny,nz,itime,jvars(6)
        real data(nx,ny,nz)
        integer i,k,start(4),count(4),status

        k = 0
        if ( jvars(2).gt.0 ) then
            k = k + 1
            start(k) = 1
            count(k) = max(1,nx)
        endif
        if ( jvars(3).gt.0 ) then
            k = k + 1
            start(k) = 1
            count(k) = max(1,ny)
        endif
        if ( jvars(4).gt.0 ) then
            k = k + 1
            start(k) = 1
            count(k) = max(1,nz)
        endif
        if ( jvars(5).gt.0 ) then
            k = k + 1
            start(k) = itime
            count(k) = 1
        else
            write(0,*) 'zreadncfile: error: time undefined in ncid '
     +           ,ncid
            call abort
        endif
        if ( lwrite ) then
            print *,'readonencfield: startvec = ',(start(i),i=1,k)
            print *,'                countvec = ',(count(i),i=1,k)
        endif
        status = nf_get_vara_real(ncid,jvars(1),start,count,data)
        if ( status.ne.nf_noerr ) call handle_err(status
     +       ,'readonencfield: nf_get_vara_real: ')

        end
