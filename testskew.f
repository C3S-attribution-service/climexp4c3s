        program testskew
!
!       test makeskew
!
        implicit none
        integer nmax
        parameter (nmax=1000000)
        integer i
        real x(nmax),s,ave,adev,sdev,var,skew,curt

        print *,'skew?'
        read *,s
        do i=1,nmax
            call makeskew(x(i),s)
        end do

        call moment(x,nmax,ave,adev,sdev,var,skew,curt)
        print *,'ave = ',ave,', should be 0'
        print *,'var = ',var,', should be 1'
        print *,'skew= ',skew,', should be ',s
        end
