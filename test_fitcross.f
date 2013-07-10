        program test_fitcross
        implicit none
        integer i,j,n,nc,n1,nc1
        real xx(100),yy(100),aa(100),bb(100),sig(1),a,b,siga,sigb,chi2,q
        real xx1(100),yy1(100)

        print *,'n?'
        read *,n
        if ( n.gt.100 ) then
            write(0,*) 'cannot handle n>100 '
            call abort
        end if
        do i=1,n
            call random_number(xx(i))
            call random_number(yy(i))
        end do
        print *,'nc?'
        read *,nc
        call fitcross(xx,yy,n,sig,0,a,b,siga,sigb,chi2,q,nc,
     +       aa,bb,.true.)
        
        do i=1,n
            n1 = 0
            do j=1,n
                if ( j.lt.i-(nc-1)/2 ) then
                    xx1(j) = xx(j)
                    yy1(j) = yy(j)
                    n1 = n1 + 1
                elseif ( j.gt.i+nc/2 ) then
                    nc1 = min(i+nc/2,nc)
                    xx1(j-nc1) = xx(j)
                    yy1(j-nc1) = yy(j)
                    n1 = n1 + 1
                end if
            end do
            call fitcross(xx1,yy1,n1,sig,0,a,b,siga,sigb,chi2,q,0,aa
     +           ,bb,.false.)
            print *,i,a,aa(i),a-aa(i),b,bb(i),b-bb(i)
        end do

        end
