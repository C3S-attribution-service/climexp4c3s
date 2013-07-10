        subroutine getrms(f1,w1,f2,w2,n,bias,mean1,var1,mean2,var2,rmse)
        implicit none
        integer n
        real f1(n),w1(n),f2(n),w2(n),bias,mean1,var1,mean2,var2,rmse
        integer i
        real ww1,ww2,ww12
        mean1 = 0
        mean2 = 0
        ww1 = 0
        ww2 = 0
        do i=1,n
            mean1 = mean1 + w1(i)*f1(i)
            ww1 = ww1 + w1(i)
            mean2 = mean2 + w2(i)*f2(i)
            ww2 = ww2 + w2(i)
        enddo
        mean1 = mean1/ww1
        mean2 = mean2/ww2
        bias = mean1 - mean2
        var1 = 0
        var2 = 0
        rmse = 0
        ww12 = 0
        do i=1,n
            var1 = var1 + w1(i)*(f1(i)-mean1)**2
            var2 = var2 + w2(i)*(f2(i)-mean2)**2
            rmse = rmse + sqrt(w1(i)*w2(i))*(f1(i)-mean1-f2(i)+mean2)**2
            ww12 = ww12 + sqrt(w1(i)*w2(i))
        enddo
        var1 = var1/ww1*n/(n-1)
        var2 = var2/ww2*n/(n-1)
        rmse = sqrt(rmse/ww12*n/(n-1))
        end
