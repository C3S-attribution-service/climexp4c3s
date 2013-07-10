        program dat2indices
THIS DOES NOT WORK
*
*       convert standar .dat files (obtained with get_area)
*       to sst*.indices format of NINO3 indices
*
        implicit none
        integer yrbeg,yrend
        parameter(yrbeg=1700,yrend=2020)
        integer ii(2)
        real nino(12,yrbeg:yrend,2:5),xx(8)
        character string*128
*
        call readdat(nino(1,yrbeg,3),yrbeg,yrend,'kaplan_nino3.dat')
        call readdat(nino(1,yrbeg,4),yrbeg,yrend,'kaplan_nino4.dat')
        call readdat(nino(1,yrbeg,5),yrbeg,yrend,'kaplan_nino34.dat')
*
        open(1,file='sstkap.indices',status='old')
        read(1,'(a)') string
        print '(a)',string
  100   continue
        read(1,*,err=900,end=200) ii,xx
        if ( abs(xx(4)-nino(ii(2),ii(1),3)).gt.0.015 ) then
            print *,'inconsistent NINO3 index at ',ii,': ',
     +            xx(4),nino(ii(2),ii(1),3)
        endif
        xx(6) = nino(ii(2),ii(1),4)
        xx(8) = nino(ii(2),ii(1),5)
        print '(i5,i3,8f8.2)',ii,xx
        goto 100
  200   continue
        stop
*
  900   print *,'error reading file sstkap.indices'
        end
