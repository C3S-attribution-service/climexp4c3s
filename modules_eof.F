      program eoff
**    parameter (mxy=11000,mxy2=6700,msteps=6000)
      parameter (mxy=11000,mxy2=2000,msteps=6000)
#include "absent.h"
      dimension f(mxy,2),cov(mxy2,mxy2),eval(mxy),mask(mxy),iadres(mxy2)
     ,         ,tserie(msteps)
      logical normpc,latico
      character*128 dump
      call readat(f,mask,ndim,ndim2,nsteps,iadres,npc,normpc,latico,
     ,                                         long,lati,pl0,pdl,dump)
      if (ndim2.gt.mxy2) stop 'increase mxy2'
      if (nsteps.gt.msteps) stop 'increase msteps'
      call eof(f,iadres,ndim,ndim2,nsteps,cov,eval,long,lati,pl0,pdl,
     ,                                               mask,latico,absent)
      call pricom(f,ndim,ndim2,nsteps,cov,eval,iadres,npc,normpc,latico,
     ,                           long,lati,pl0,pdl,dump,absent,tserie)
      end
      subroutine pricom(f,ndim,ndim2,nsteps,cov,eval,iadres,npc,normpc,
     ,                    latico,long,lati,pl0,pdl,dump,absent,tserie)
      dimension f(ndim,2),cov(ndim2,ndim2),eval(ndim),iadres(ndim2),
     ,          tserie(*)
      logical normpc,latico
      character*128 dump
      if (npc.le.0) return
      if (npc.gt.ndim2) then
         npc = ndim2
         write(6,'(a24,i4)') 'The # PC''s is reduced to',npc
      end if
      call zero(eval(ndim2+1),ndim-ndim2)
      write(2) -1,503,9999,ndim
      write(2) eval
      npc2 = npc
      nfield = ndim / (long*lati)
      do 10 i=1,npc
         call scatter(cov(1,i),f(1,2),iadres,ndim,ndim2,absent)
         if (latico) then
            call correc(f(1,2),long,lati,nfield,pl0,pdl,'fro',absent)
            call norm(f(1,2),ndim,absent)
         end if
         if (normpc) then
            if (eval(i).gt.0.0) then
               eval(i) = 1./sqrt(eval(i))
               call scale(f(1,2),eval(i),ndim,absent)
            else
               npc2 = min(npc2,i)
               goto 20
            end if
         else
            eval(i) = 1.0
         end if
         write(2) i,503,9999,ndim
         write(2) (f(k,2),k=1,ndim)
10    continue
20    if (dump(1:7).ne.'dumppc ') return
      if (npc.ne.npc2) print*,'!!! WARNING !!!, npc reduced to:',npc 
      npc = npc2
      absen = absent * .9
      if (dump(1:7).eq.'dumppc ') then
         do 50 j=1,npc
            rewind (1)
            write(2) j,504,9999,nsteps
            do 40 i=1,nsteps
               read(1) i1,i2,i3,i4
               read(1) (f(k,1),k=1,i4)
               if (latico)call correc(f,long,lati,nfield,pl0,pdl,'to ',
     ,                                                          absent)
               do 30 k=1,ndim2
30             f(k,1) = f(iadres(k),1)
               vv = vpnorm(f(1,1),cov(1,j),ndim2,absent)     
               if (vv.lt.absen) then
                  tserie(i) = vv * eval(j)
               else
                  tserie(i) = absent
               end if
40          continue
            write(2) (tserie(i),i=1,nsteps)
50       continue
      else
         rewind (1)
         do 80 i=1,nsteps
            read(1) i1,i2,i3,i4
            read(1) (f(k,1),k=1,i4)
            if (latico)call correc(f,long,lati,nfield,pl0,pdl,'to ',
     ,                                                       absent)
            do 60 k=1,ndim2
60          f(k,1) = f(iadres(k),1)
            write(2) i1,504,9999,npc
            do 70 j=1,npc
               vv = vpnorm(f(1,1),cov(1,j),ndim2,absent)     
               if (vv.lt.absen) then
                  tserie(j) = vv * eval(j)
               else
                  tserie(j) = absent
               end if
70          continue
            write(2) (tserie(j),j=1,npc)
80       continue
      endif
      close (2)
      return
      end
      function vpnorm(a,b,ndim,absent)
      dimension a(ndim),b(ndim)
      bnorm = 0.0
      vpnorm = 0.0
      absen = absent * .9
      do 10 i=1,ndim
         if (a(i).lt.absen) then
            bnorm = bnorm + b(i) * b(i)
            vpnorm = vpnorm + a(i) * b(i)
         end if
10    continue
      if (bnorm.eq.0.0) then
         vpnorm = absent
      else
         vpnorm = vpnorm / sqrt(bnorm)
      end if
      return
      end 
      subroutine scale(a,scal,ndim,absent)
      dimension a(ndim)
      absen = absent * .9
      do 10 i=1,ndim
10    if (a(i).lt.absen) a(i) = a(i) * scal
      return
      end 
      subroutine scatter(a,b,iadres,ndim,ndim2,absent)
      dimension a(ndim2),b(ndim),iadres(ndim2)
      do 10 i=1,ndim
10    b(i) = absent
      do 20 i=1,ndim2
20    b(iadres(i)) = a(i)
      return
      end
      subroutine compres(a,b,n,new,absent)
      dimension a(n),b(*)
      new = 0
      absen = absent * .9
      do 10 i=1,n
         if (a(i).lt.absen) then
            new = new + 1
            b(new) = a(i)
         end if
10    continue
      return
      end
      subroutine readat(f,mask,ndim,ndim2,nsteps,iadres,npc,normpc,
     ,                             latico,long,lati,pl0,pdl,dump)
#include "absent.h"
      dimension f(*),mask(*),iadres(*)
      character*128 infile,outfile,mode,mode2,dump,ainput
      logical normpc,latico
      mode2 = 'nocorrection'
      call inputt(5,numdat)
      infile = ainput()
      outfile = ainput()
      minper = iinput()
      npc = iinput()
      mode = ainput()
      mode2 = ainput()
      dump = ainput()
      write(6,11) infile,outfile,minper,npc,mode,mode2,absent
11    format('--> EOF program <--',/,
     ,      'Input  file:',a,/,'Output file:',a,/,
     ,       'Minimum % data required:',i4,/,
     ,       '# principle components:',i4,/,
     ,       'Normalisation mode:',a,/,
     ,       'Correction mode:',a,/,'Absent value:',e8.2)
      if (mode(1:6).eq.'normpc') then
         print*,'The principle components are normalised'
         normpc = .true.
      else
         normpc = .false.
      end if
      if (dump(1:7).eq.'dumppc ') then
         print*,'Principle components will be dumped per EOF'
      else if (dump(1:7).eq.'dumppcd') then
         print*,'Principle components will be dumped per timestep'
      else
         print*,'Principle components will NOT be dumped'
      end if
      open(1,file=infile,form='unformatted')
      open(2,file=outfile,form='unformatted')
      nsteps = 0
      i0 = 1
10    read(1,end=999) i1,i2,i3,i4
      if (nsteps.eq.0) call izero(mask,i4)
      read(1) (f(i),i=1,i4)
      do 20 i=1,i4
         if (f(i).ge.absent*.9) then
            mask(i) = mask(i)+1
         end if
20    continue
      nsteps = nsteps + 1
      goto 10
999   write(6,44) nsteps
44    format('# steps read:',i8)
      rewind (1)
      ndim = i4
      ndim2 = 0
      do 30 i=1,ndim
         if (nsteps-mask(i).ge.minper*nsteps/100) then
            ndim2 = ndim2 + 1
            iadres(ndim2) = i
            mask(ndim2) = mask(i)
         end if
30    continue
      write(6,55) ndim-ndim2,100-minper
55    format(i4,' points have more than',i4,'% absent values')
      if (mode2(1:7).eq.'laticor') then
         print*,'Signal will undergo latitudinal correction'
         latico = .true.
         call inputt(5,numdat)
         long = iinput()
         lati = iinput()
         rla0 = finput()
         dla = finput()
         pi = 3.1415926535898
         pl0 = rla0*pi/180.
         pdl =  dla*pi/180.
         write(6,22) long,lati,rla0,dla,
     ,               rla0+(lati-1)*dla,
     ,               (sqrt(cos(pl0+(j-1)*pdl)),j=1,lati)
22       format('# longitude:',i4,' # latitude:',i4,/,
     ,          'latitude of first row:',f10.4,' increment:',f10.4,/,
     ,          'latitude of last row:',f10.4,/,
     ,          'Latitudinal weights:',(t21,10f6.2))
         if ((ndim.ne.long*lati).and.(mod(ndim,long*lati).eq.0)) then
            write(6,33) ndim/(long*lati)
33          format('---------------------------------------------------'
     ,          ,/,'-->> The field is assumed to be a vector field <<--'
     ,      ,/,'-->> in',i2,' dimensions                           <<--'
     ,         ,/,'---------------------------------------------------')
         else if (long*lati.ne.ndim) then
            print*,'long,lati,ndim:',long,lati,ndim
            stop 'inconsistency long/lati <-> dimension'
         end if
      else
         latico = .false.
         long = 1
         lati = 1
      end if
      return
      end
      subroutine zero(a,n)
      dimension a(*)
      do 10 i=1,n
10    a(i) = 0.0
      return
      end
      subroutine izero(ia,n)
      dimension ia(*)
      do 10 i=1,n
10    ia(i) = 0
      return
      end
      subroutine eigen(a,eval,n)
      parameter (mn=10000)
      dimension a(n,n),eval(n),fv1(mn)
      if (n.gt.mn) stop 'increase mn in eigen'
*org  call tred2(n,n,a,eval,fv1,a)
*org  call tql2(n,n,eval,fv1,a,ierr)
      call xtred2(n,n,a,eval,fv1,a)
      call xtql2(n,n,eval,fv1,a,ierr)
      if (ierr.ne.0) stop 'error in eigen'
      do 10 i=1,n/2
         ee = eval(i)
         eval(i) = eval(n-i+1)
         eval(n-i+1) = ee
         call fmove(a(1,i),fv1,n)
         call fmove(a(1,n-i+1),a(1,i),n)
         call fmove(fv1,a(1,n-i+1),n)
10    continue
      return
      end
      subroutine xtql2(nm,n,d,e,z,ierr)
c
cc    Eispack routine for determining the eigenvalues and eigenvectors of
cc    a symmetric tridiagonal matrix
c
      dimension d(n),e(n),z(nm,n)
      eps = 1.0e-20
      ierr = 0
      if (n.eq.1) goto 1001
      do 100 i=2,n
100   e(i-1) = e(i)
      f = 0.0e0
      b = 0.0e0
      e(n) = 0.0e0
      do 240 l=1,n
         j = 0
         h = eps * (abs(d(l))+abs(e(l)))
         if (b.lt.h) b = h
         do 110 m=l,n
            if (abs(e(m)).le.b) goto 120
110      continue
120      if (m.eq.l) goto 220
130      if (j.eq.30) goto 1000
         j = j + 1
         l1 = l + 1
         g = d(l)
         p = (d(l1) - g) / (2.0 * e(l))
         r = sqrt(p*p+1.0)
         d(l) = e(l) / (p+ sign(r,p))
         h = g - d(l)
         do 140 i=l1,n
140      d(i) = d(i) - h
         f = f + h
         p = d(m)
         c = 1.0
         s = 0.0
         mml = m - l
         do 200 ii=1,mml
            i = m - ii
            g = c * e(i)
            h = c * p
            if (abs(p).lt.abs(e(i))) goto 150
            c = e(i) / p
            r = sqrt(c*c+1.0)
            e(i+1) = s * p * r
            s = c / r
            c = 1.0 / r
            goto 160
150         c = p / e(i)
            r = sqrt(c*c+1.0)
            e(i+1) = s * e(i) * r
            s = 1.0 / r
            c = c * s
160         p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
            do 180 k=1,n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
180         continue
200      continue
         e(l) = s * p
         d(l) = c * p
         if (abs(e(l)).gt.b) goto 130
220      d(l) = d(l) + f
240   continue
      do 300 ii=2,n
         i = ii - 1
         k = i
         p = d(i)
         do 260 j=ii,n
            if (d(j).ge.p) goto 260
            k = j
            p = d(j)
260      continue
         if (k.eq.i) goto 300
         d(k) = d(i)
         d(i) = p
         do 280 j=1,n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
280      continue
300   continue
      goto 1001
1000  ierr = l
1001  return
      end
      subroutine xtred2(nm,n,a,d,e,z)
      dimension a(nm,n),d(n),e(n),z(nm,n)
      do 100 i=1,n
      do 100 j=1,i
100   z(i,j) = a(i,j)
      if (n.eq.1) goto 320
      do 300 ii=2,n
         i = n + 2 - ii
         l = i - 1
         h = 0.0
         scale = 0.0
         if (l.lt.2) goto 130
         do 120 k=1,l
120      scale = scale + abs(z(i,k))
         if (scale.ne.0.0) goto 140
130      e(i) = z(i,l)
         goto 290
140      do 150 k=1,l
            z(i,k) = z(i,k) / scale
            h = h + z(i,k) * z(i,k)
150      continue
         f = z(i,l)
         g = -sign(sqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         z(i,l) = f - g
         f = 0.0
         do 240 j=1,l
            z(j,i) = z(i,j) / h
            g = 0.0
            do 180 k=1,j
180         g = g + z(j,k) * z(i,k)
            jp1 = j + 1
            if (l.lt.jp1) goto 220
            do 200 k=jp1,l
200         g = g + z(k,j) * z(i,k)
220         e(j) = g / h
            f = f + e(j) * z(i,j)
240      continue
         hh = f / (h+h)
         do 260 j=1,l
            f = z(i,j)
            g = e(j) - hh * f
            e(j) = g
            do 260 k=1,j
               z(j,k) = z(j,k) - f * e(k) - g * z(i,k)
260      continue
290      d(i) = h
300   continue 
320   d(1) = 0.0
      e(1) = 0.0
      do 500 i=1,n
         l = i - 1
         if (d(i).eq.0.0) goto 380
         do 360 j=1,l
            g = 0.0
            do 340 k=1,l
340         g = g + z(i,k) * z(k,j)
            do 360 k=1,l
               z(k,j) = z(k,j) - g * z(k,i)
360      continue
380      d(i) = z(i,i)
         z(i,i) = 1.0
         if (l.lt.1) goto 500
         do 400 j=1,l
            z(i,j) = 0.0
            z(j,i) = 0.0
400      continue
500   continue
      return
      end          
      subroutine fmove(a,b,n)
      dimension a(*),b(*)
      do 10 i=1,n
10    b(i) = a(i)
      return
      end     
      subroutine covar(f,iadres,b,nra,nra2,nca,long,lati,pl0,pdl,mask,
     ,                                                  latico,absent)
      dimension f(nra),iadres(nra2),b(nra2,nra2),mask(nra)
      logical latico
      absen = absent*.9
      call zero(b,nra2*nra2)
      call izero(mask,nra2)
      nfield = nra/(long*lati)
      do 40 k=1,nca
         read(1) i1,i2,i3,i4
         read(1) f
         if (latico)call correc(f,long,lati,nfield,pl0,pdl,'to ',absent)
         do 10 i=1,nra2
10       f(i) = f(iadres(i))
         do 30 j=1,nra2
            if (f(j).lt.absen) then
               do 20 i=1,j-1
                  if (f(i).le.absen) then
                     b(i,j) = b(i,j) + f(i) * f(j)
                     b(j,i) = b(j,i) + 1
                  endif
20             continue
               b(j,j) = b(j,j) + f(j) * f(j)
               mask(j) = mask(j) + 1
            endif
30       continue
40    continue
      do 60 j=1,nra2
         do 50 i=1,j-1
            b(i,j) = b(i,j) / b(j,i)
50       b(j,i) = b(i,j) 
         b(j,j) = b(j,j) / mask(j)
60    continue
      return
      end
      subroutine eof(f,iadres,ndim,ndim2,nsteps,cov,eval,long,lati,pl0,
     ,                                           pdl,mask,latico,absent)
      dimension f(ndim),iadres(ndim2),cov(ndim2,ndim2),eval(ndim2),
     ,          mask(ndim)
      logical latico
      call covar(f,iadres,cov,ndim,ndim2,nsteps,long,lati,pl0,pdl,
     ,                                               mask,latico,absent)
      call eigen(cov,eval,ndim2)
      sum = 0.0
      do 30 i=1,ndim2
30    sum = sum + eval(i)
      print*,'sum eval:',sum
      print*,'eigenvalues :'
      write(6,'(8e10.3)') (eval(i),i=1,ndim2)
      print*,'relative importance :'
      write(6,'(10f8.3)') (eval(i)*100./sum,i=1,ndim2)
      return
      end
      subroutine correc(f,long,lati,nfield,rla0,dla,mode,absent)
      dimension f(long,lati,nfield)
      character*3 mode
      absen = absent * .9
      if (mode.eq.'to ') then
         do 20 j=1,lati
            factor = sqrt(cos(rla0+(j-1)*dla))
            do 10 k=1,nfield
            do 10 i=1,long
10          if (f(i,j,k).lt.absen) f(i,j,k) = f(i,j,k) * factor 
20       continue
      else if (mode.eq.'fro') then
         do 40 j=1,lati
            factor = 1./sqrt(cos(rla0+(j-1)*dla))
            do 30 k=1,nfield
            do 30 i=1,long
30          if (f(i,j,k).lt.absen) f(i,j,k) = f(i,j,k) * factor 
40       continue
      else
         stop 'unkown mode'
      end if
      return
      end
      subroutine norm(a,n,absent)
      dimension a(n)
      absen = absent*.9
      rnorm = 0.0
      do 10 i=1,n
10    if (a(i).lt.absen) rnorm = rnorm + a(i) * a(i)
      if (rnorm.le.0.0) stop 'zero norm in norm'
      rnorm = 1./sqrt(rnorm)
      do 20 i=1,n
20    if (a(i).lt.absen) a(i) = a(i) * rnorm
      return
      end
      function vprod(a,b,n,absent)
      dimension a(*),b(*)
      vprod = 0.0
      absen = absent * .9
      do 10 i=1,n
10    if ((a(i).lt.absen).and.(b(i).lt.absen)) vprod = vprod + a(i)*b(i)
      return
      end
