        subroutine getdf(df,month,n,ndup,decor,lsum,ndiff,nperyear)
        implicit none
        real df,decor
        integer month,n,ndup,lsum,ndiff,nperyear
        if ( month.eq.0 ) then
            df = (n-ndup)/(lsum + decor)/real(ndiff) - 2
        else
            df = (n-ndup)/(1 + (lsum-1)/nperyear + decor/nperyear)/ndiff
     +           - 2
        endif
        end
