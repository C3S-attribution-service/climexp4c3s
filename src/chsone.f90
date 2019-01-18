subroutine chsone(bins,ebins,nbins,nconstraints,df,chi2,prob)
!
!   compute \chi^2 statistic of a histohgram compared with the expected distribution
!   and p-value that the histogram is drawn from the expected distribution.
!   uses (F)GSL rather than Numerical Recipes
!
    use fgsl
    implicit none
    integer :: nconstraints,nbins
    real :: chi2,df,prob,bins(nbins),ebins(nbins)
    integer :: i
    real(fgsl_double) :: arg1,arg2

    df = nbins - nconstraints
    chi2 = 0
    do i=1,nbins
        if ( ebins(i) > 0 ) then
            chi2 = chi2 + (bins(i)-ebins(i))**2/ebins(i)
        else if ( bins(i) /= 0 ) then
            write(0,*) 'chsone: error: ebins(',i,') <= 0 ',ebins(i)
            call exit(-1)
        end if
    end do
    arg1 = dble(df)/2
    arg2 = dble(chi2)/2
    prob = fgsl_sf_gamma_inc_Q(arg1,arg2)
end subroutine chsone