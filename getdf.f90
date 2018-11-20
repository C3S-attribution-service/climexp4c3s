subroutine getdf(df,month,n,ndup,decor,lsum,ndiff,nperyear)
    implicit none
    integer,intent(in) :: month,n,ndup,lsum,ndiff,nperyear
    real,intent(in) :: decor
    real,intent(out) :: df
    if ( month == 0 ) then
        df = (n-ndup)/(lsum + decor)/real(ndiff) - 2
    else
        df = (n-ndup)/(1 + (lsum-1)/nperyear + decor/nperyear)/ndiff - 2
    endif
end subroutine getdf
