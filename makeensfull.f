        subroutine makeensfull(ndup,nperyear,data,npermax,yrbeg,yrend
     +        ,nens1,nens2,validens,lwrite)
*
*       copy ensemble members so that there is the same
*       number of valid ones at every time step
*
        implicit none
        integer npermax,nperyear,ndup(nperyear),yrbeg,yrend,nens1,nens2
     +        ,validens(nens2+1),n
        real data(npermax,yrbeg:yrend,0:nens2)
        logical lwrite
        integer i,j,k,l,mens,iens,jens,seed(100)
        real xran
*
*	reset random number generator to get identical numbers each time
*	the routine is called
        call random_seed(size=i)
        n=1234567
        do j=1,i
            n = n*j
            seed(j) = n
        enddo
	call random_seed(put=seed(1:i))
*
        do j=1,nperyear
            do i=yrbeg,yrend
                mens = 0
                do iens=nens1,nens2
                    if ( data(j,i,iens).lt.1e33 ) then
                        mens = mens + 1
                        validens(mens) = iens
                    endif
                enddo
                if ( mens.gt.1 .and. mens.le.(nens2-nens1) ) then
*                   some - but not all - are undefined
*                   copy as many whole sets in as fit
                    if ( lwrite ) print '(a,i4,i2,a,i4,a,2i4)','At ',i,j
     +                    ,' only found ',mens
     +                    ,' ensemble members out of ',nens1,nens2
                    k = (nens2-nens1+1)/mens
                    iens=nens1
                    do l=1,k-1
                        do jens=1,mens
  300                       continue
                            if ( data(j,i,iens).lt.1e33 ) then
                                iens = iens + 1
                                goto 300
                            endif
                            if ( lwrite ) print *,'copying member '
     +                            ,validens(jens),' to ',iens,data(j,i
     +                            ,validens(jens)),data(j,i,iens)
                            ndup(j) = ndup(j) + 1
                            data(j,i,iens) = data(j,i,validens(jens))
                        enddo
                    enddo
*       and fill out the rest with random members
                    do l=nens1+k*mens,nens2
*       check...
                        do k=1,mens
                            if ( validens(k).ge.0 ) goto 305
                        enddo
                        write(0,*) 'makeensfull: error: no more members'
                        call abort
  305                   continue
  310                   continue
                        call random_number(xran)
                        k = 1 + int(mens*xran)
                        if ( k.le.0 .or. k.gt.mens .or. validens(k).lt.0
     +                        ) goto 310
  320                   continue
                        if ( data(j,i,iens).lt.1e33 ) then
                            iens = iens + 1
                            goto 320
                        endif
                        if ( lwrite ) print *,'copying member '
     +                        ,validens(k),' to ',iens,data(j,i
     +                        ,validens(k)),data(j,i,iens)
                        ndup(j) = ndup(j) + 1
                        data(j,i,iens) = data(j,i,validens(k))
                        validens(k) = -1
                    enddo
                endif
            enddo
        enddo
        end
