subroutine svdfit(x,y,sig,ndata,a,ma,u,v,w,mp,np,chisq,func)
!
!   wrapper that calls GSL with interface of Numerical recipes routine
!   u,v,w are working spaces that are not used in this version, v is used to return the covariance matrix
!
    use fgsl
    implicit none
    integer :: ndata,ma,mp,np
    real :: x(ndata),y(ndata),sig(ndata),a(ma),u(mp,np),v(np,np),w(np),chisq
    external :: func
    integer :: i,j
    real,allocatable :: vector(:)
    real*8,allocatable :: dmat(:,:),dy(:),fitparams(:),dcov(:,:)
    integer(fgsl_int) :: status
    integer(fgsl_size_t) :: n,p
    type(fgsl_multifit_linear_workspace) :: workspace
    type(fgsl_vector) :: yvec,cvec
    type(fgsl_matrix) :: xmat,cov
    real(fgsl_double) :: chisquare
    
    n = ndata
    p = ma
    workspace = fgsl_multifit_linear_alloc(n, p)
    allocate(vector(ma),dmat(ma,ndata),fitparams(ma),dcov(ma,ma),dy(ndata))
    do i=1,n
        call func(real(i),vector,ma) ! single precision arguments
        do j=1,ma
            dmat(j,i) = vector(j)
        end do
    end do
    dy = y ! convert to double preicision
    xmat = fgsl_matrix_init(1.0_fgsl_double)
    status = fgsl_matrix_align(dmat,p,p,n,xmat)
    yvec = fgsl_vector_init(1.0_fgsl_double)
    status = fgsl_vector_align(dy,n,yvec,n,0_fgsl_size_t,1_fgsl_size_t)
    cvec = fgsl_vector_init(1.0_fgsl_double)
    status = fgsl_vector_align(fitparams,p,cvec,p,0_fgsl_size_t,1_fgsl_size_t)
    cov = fgsl_matrix_init(1.0_fgsl_double)
    status = fgsl_matrix_align(dcov,p,p,p,cov)
    do i=1,ndata
        if ( sig(i) /= 1 ) then
            write(0,*) 'svdfit: error: wrapper cannot handle weights yet'
            call exit(-1)
        end if
    end do
    status = fgsl_multifit_linear(xmat, yvec, cvec, cov, chisquare, workspace)
    a = fitparams ! should be buried inside cvec
    do j=1,ma
        do i=1,ma
            v(i,j) = dcov(i,j)
        end do
    end do
    chisq = chisquare
    ! Numerical Recipes gives varianbces that are a factor 1/(chi2/dof) larger
    v = v/(chisq/(n-ma))
    call fgsl_matrix_free(xmat)
    call fgsl_vector_free(yvec)
    call fgsl_vector_free(cvec)
    call fgsl_matrix_free(cov)
    deallocate(vector,dmat,fitparams,dcov)
    call fgsl_multifit_linear_free(workspace)
end subroutine svdfit

subroutine svdvar(v,ma,np,w,cvm,ncvm)
!
!   just copy the covariance matrix that was computed by fgsl_multifit_linear
!
    integer :: ma,np,ncvm
    real :: v(np,np),w(np),cvm(ncvm,ncvm)
    integer :: i,j
    
    do j=1,np
        do i=1,np
            cvm(i,j) = v(i,j)
        end do
    end do
end subroutine svdvar