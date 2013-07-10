        module loess_d_mod
        contains

        subroutine condition(surface, new_stat, trace_hat_in, surf_stat)
        implicit none
        character surface*(*), new_stat*(*), trace_hat_in*(*),
     +       surf_stat*(*)

	if ( surface.eq.'interpolate' ) then
            if ( new_stat.eq.'none' ) then
                surf_stat = 'interpolate/none'
            else if ( new_stat.eq."exact" ) then
                surf_stat = 'interpolate/exact'
            else if ( new_stat.eq.'approximate' ) then
                if ( trace_hat_in.eq.'approximate' ) then
                    surf_stat = 'interpolate/2.approx'
                else if( trace_hat_in.eq.'exact' ) then
                    surf_stat = 'interpolate/1.approx'
		end if
            end if
	else if ( surface.eq.'direct' ) then
            if ( new_stat.eq.'none' ) then
                surf_stat = 'direct/none'
            else if ( new_stat.eq.'exact' ) then
                surf_stat = 'direct/exact'
            else if ( new_stat.eq.'approximate' ) then
                surf_stat = 'direct/approximate'
            end if
	end if
        if ( .false.) print *,'surf_stat = ',trim(surf_stat)
        end subroutine

        subroutine loess_raw(y, x, weights, robust, d, n, max_kd, span,
     +       degree, nonparametric, order_drop_sqr, sum_drop_sqr, cell,
     +       surf_stat, surface, parameter, a, xi, vert, vval, diagonal,
     +       trL, one_delta,two_delta, setLf)
        implicit none
        integer d, n, max_kd, parameter(7), a(max_kd), degree,
     +       nonparametric, order_drop_sqr(8), sum_drop_sqr
        double precision y(n), x(n,d), weights(n), robust(n), span, cell
     +       , surface(n),xi(max_kd), vert(max_kd), vval(max_kd),
     +       diagonal(n), trL,one_delta, two_delta
        character surf_stat*(*)
        logical setLf
        integer nsing,i,k
        integer liv,lv,tau
        integer,allocatable :: iv(:)
        double precision,allocatable :: v(:),hat_matrix(:,:),LL(:,:)

        trL = 0
        call loess_workspace(d, n, span, degree, nonparametric,
     +       order_drop_sqr, sum_drop_sqr, setLf, iv, v, lv, liv, tau)
        v(2) = cell
	if ( surf_stat.eq.'interpolate/none' ) then
            call lowesb(x, y, robust, 0, 0, iv, liv, lv, v);
            call lowese(iv, liv, lv, v, n, x, surface)
            call loess_prune(parameter, a, max_kd, xi, vert, vval, iv, v
     +           , lv,liv)
	else if ( surf_stat.eq.'direct/none' ) then
            call lowesf(x, y, robust, iv, liv, lv, v, n, x, 0, 0,
     +           surface)
	else if ( surf_stat.eq.'interpolate/1.approx' ) then
            call lowesb(x, y, weights, diagonal, 1, iv, liv, lv, v)
            call lowese(iv, liv, lv, v, n, x, surface);
            nsing = iv(30)
            do i=1,n
                trL = trL + diagonal(i)
            end do
            call lowesa(trL, n, d, tau, nsing, one_delta, two_delta)
            call loess_prune(parameter, a, max_kd, xi, vert, vval, iv, v
     +           , lv,liv)
        else if ( surf_stat.eq.'interpolate/2.approx' ) then
            call lowesb(x, y, robust, 0, 0, iv, liv, lv, v)
            call lowese(iv, liv, lv, v, n, x, surface)
            nsing = iv(30)
            call ehg196(tau, d, span, trL)
            call lowesa(trL, n, d, tau, nsing, one_delta, two_delta)
            call loess_prune(parameter, a, max_kd, xi, vert, vval, iv, v
     +           , lv,liv)
	else if ( surf_stat.eq.'direct/approximate' ) then
            call lowesf(x, y, weights, iv, liv, lv, v, n, x, diagonal, 1
     +           , surface)
            nsing = iv(30)
            do i=1,n
                trL = trL + diagonal(i)
            end do
            call lowesa(trL, n, d, tau, nsing, one_delta, two_delta)
	else if ( surf_stat.eq.'interpolate/exact' ) then
            allocate(hat_matrix(n,n))
            hat_matrix = 0
            allocate(LL(n,n))
            LL = 0
            call lowesb(x, y, weights, diagonal, 1, iv, liv, lv, v)
            call lowesl(iv, liv, lv, v, n, x, hat_matrix)
            call lowesc(n, hat_matrix, LL, trL, one_delta, two_delta)
            call lowese(iv, liv, lv, v, n, x, surface)
            call loess_prune(parameter, a, max_kd, xi, vert, vval, iv, v
     +           , lv,liv)
            deallocate(hat_matrix)
            deallocate(LL)
	else if ( surf_stat.eq.'direct/exact' ) then
            allocate(hat_matrix(n,n))
            hat_matrix = 0
            allocate(LL(n,n))
            LL = 0
            call lowesf(x, y, weights, iv, liv, lv, v, n, x, hat_matrix,
     +           2, surface)
            call lowesc(n, hat_matrix, LL, trL, one_delta, two_delta)
            do i=1,n
                diagonal(i) = hat_matrix(i,i)
            end do
            deallocate(hat_matrix)
            deallocate(LL)
	end if
        deallocate(iv)
        deallocate(v)
        end subroutine

        subroutine loess_workspace(d, n, span, degree, nonparametric,
     +       order_drop_sqr, sum_drop_sqr, setLf, iv, v, lv, liv, tau)
        implicit none
        integer d, n, degree, nonparametric, order_drop_sqr(8),
     +       sum_drop_sqr, lv, liv, tau
        double precision span
        logical setLf
        integer,allocatable :: iv(:)
        double precision,allocatable :: v(:)
	integer tau0, nvmax, nf, version, i
        version = 106

        if ( .false. ) then
            print *,'d = ',d
            print *,'n = ',n
            print *,'span = ',span
            print *,'degree = ',degree
            print *,'nonparametric = ',nonparametric
            print *,'sum_drop_sqr = ',sum_drop_sqr
            do i=1,d
                print *,'order_drop_sqr(',i,') = ',order_drop_sqr(i)
            end do
        end if

	nvmax = max(200, n)
        nf = min(n, int(n*span))
        if ( degree.gt.1 ) then
            tau0 = (d+2)*(d+1)/2
        else
            tau0 = d+1
        end if
        tau = tau0 - sum_drop_sqr
        lv = 50 + (3*d + 3)*nvmax + n + (tau0 + 2)*nf
	liv = 50 + (2**d + 4)*nvmax + 2*n
	if ( setLf ) then
            lv = lv + (D + 1)*nf*nvmax
            liv = liv + nf*nvmax
	end if
        if ( .false. ) print *,'liv,lv,tau,nf,nvmax,setLf = ',liv,lv,tau
     +       ,nf,nvmax,setLf
        allocate(iv(liv))
        allocate(v(lv))
        call lowesd(version, iv, liv, lv, v, d, n, span, degree, nvmax,
     +       setLf)
        iv(33) = nonparametric
        do i=1,d
            iv(i+40) = order_drop_sqr(i)
        end do
        end subroutine

        subroutine loess_prune(parameter, a, max_kd, xi, vert, vval, iv,
     +       v, lv, liv)
        implicit none
        integer max_kd, parameter(7), a(max_kd), liv, lv, iv(liv)
        double precision xi(max_kd), vert(max_kd), vval(max_kd), v(lv)
        integer d, vc, a1, v1, xi1, vv1, nc, nv, nvmax, i, j, k
	
	d = iv(2)
	vc = iv(4) - 1
	nc = iv(5)
	nv = iv(6)
	a1 = iv(7) - 1
	v1 = iv(11) - 1
	xi1 = iv(12) - 1
	vv1 = iv(13) - 1
	nvmax = iv(14)

        do i=1,5
            parameter(i) = iv(i+1)
        end do
	parameter(6) = iv(22) - 1
	parameter(7) = iv(15) - 1

        do i=1,d
            k = nvmax*(i-1)
            vert(i) = v(1 + v1 + k)
            vert(i + d) = v(1 + v1 + vc + k)
	end do
        do i=1,nc
            xi(i) = v(xi1 + i)
            a(i) = iv(a1 + i)
	end do
	k = (d + 1)*nv
        do i=1,k
            vval(i) = v(vv1 + i)
        end do
        end subroutine

        end module

        subroutine loess(fitted_values,fitted_residuals,pseudovalues
     +       ,diagonalout,robustout,divisorout,enp,s,one_delta,two_delta
     +       ,yin,xin,weightsin,n,d,degree,infamily,spanin,cellin)
!
!       translation into Fortran of the C routine loess.c from netlib
!       no structs for the time being
!       degree = 1,2: linear or quadrtic local fit
!       infamily = 'gaussian','symmetric'
!
        use loess_d_mod
        implicit none
        integer n,d,degree
        real fitted_values(n),fitted_residuals(n),pseudovalues(n)
     +       ,diagonalout(n),robustout(n),divisorout(d),yin(n),xin(n,d),
     +       weightsin(n),spanin,cellin,enp,s,one_delta,two_delta
        character infamily*(*)
        integer i,j,max_kd,iterations,nonparametric,sum_drop_sqr
     +       ,sum_parametric
        double precision trace_hat_out,new_cell,trL_tmp,d1_tmp,d2_tmp
     +       ,sum,mean,delta1,delta2,sum_squares,trL,span,cell
        integer,allocatable :: order_parametric(:),parameter(:),a(:)
     +       ,order_drop_sqr(:),param_tmp(:),a_tmp(:)
        double precision,allocatable :: x_tmp(:,:),x(:,:),y(:),temp(:)
     +       ,xi_tmp(:),vert_tmp(:,:),vval_tmp(:,:),diag_tmp(:)
     +       ,pseudo_resid(:),weights(:),fv(:),fr(:),pv(:),robust(:)
     +       ,diagonal(:),divisor(:)
        double precision,allocatable :: xi(:),vert(:),vval(:)
        character trace_hat*12,family*9,new_stat*12,surface*12,
     +       surf_stat*25
        logical normalize,setLf
        logical drop_square(8),parametric(8)

        trL_tmp = 0
        d1_tmp = 0
        d2_tmp = 0

        if ( degree.le.0 ) degree = 2
        if ( spanin.le.0 ) then
            span = 0.75
        else
            span = spanin
        endif
        if ( cellin.le.0 ) then
            cell = 0.2
        else
            cell = cellin
        end if
        normalize = .true.
        drop_square = .false.
        parametric = .false.
        surface = 'interpolate'
        if ( infamily.eq.' ' ) then
            family = 'gaussian'
        else
            family = infamily
        end if
        if ( family(1:4).eq.'gaus' ) then
            iterations = 0
        else
            iterations = 4
        end if
        if ( n > 500 ) then
            trace_hat = 'approximate'
        else
            trace_hat = 'exact'
        end if
            
        max_kd = min(n,200)
        one_delta = 0
        two_delta = 0
        trace_hat_out = 0
        
        allocate(parameter(7))
        allocate(a(max_kd),xi(max_kd),vert(max_kd),vval(max_kd))

        allocate(x(n,d),x_tmp(n,d),temp(n),y(n),weights(n))
        allocate(fv(n),fr(n),pv(n))
        allocate(a_tmp(max_kd),xi_tmp(max_kd),vert_tmp(2,d))
        allocate(vval_tmp(d+1,max_kd),diag_tmp(n),param_tmp(7))
        allocate(order_parametric(d),order_drop_sqr(d))
        if ( iterations.gt.0 ) then
            allocate(pseudo_resid(n))
        end if
        allocate(diagonal(n),robust(n),divisor(d))

        new_cell = span*cell
        robust = 1
        x_tmp = xin
        y = yin
        weights = weightsin
        if ( normalize .and. d.gt.1 ) then
            write(0,*) 'this case not yet ready'
        else
            divisor = 1
        end if
        j = d-1
        nonparametric = 0
        sum_drop_sqr = 0
        sum_parametric = 0
        do i=1,d
            if ( drop_square(i) ) sum_drop_sqr = sum_drop_sqr + 1
            if ( parametric(i) ) then
                sum_parametric = sum_parametric + 1
                order_parametric(j) = i
                j = j - 1
            else
                nonparametric = nonparametric + 1
                order_parametric(nonparametric) = i
            end if
        end do
        do i=1,d
            if ( drop_square(order_parametric(i)) ) then
                order_drop_sqr(i) = 1
            else
                order_drop_sqr(i) = 2 
            end if
            do j=1,n
                x(j,i) = x_tmp(j,order_parametric(i));
            end do
        end do
        if ( degree.eq.1 .and. sum_drop_sqr.gt.0 ) then
            write(0,*) 'Specified the square of a factor predictor '//
     +           'to be dropped when degree = 1'
            call abort
        end if
        if ( d.eq.1 .and. sum_drop_sqr.gt.0 ) then
            write(0,*) 'Specified the square of a predictor to be '//
     +           'dropped with only one numeric predictor'
            call abort
        end if
        if ( sum_parametric.eq.d ) then
            write(0,*) 'Specified parametric for all predictors'
            call abort
        end if
        do j=0,iterations
            new_stat = 'approximate'
            robust = weights*robust
            call condition(surface,new_stat,trace_hat,surf_stat)
            setLf = (surf_stat.eq.'interpolate/exact')
            call loess_raw(y, x, weights, robust, d, n, max_kd, span,
     +           degree, nonparametric, order_drop_sqr, sum_drop_sqr,
     +           new_cell, surf_stat, fv, parameter, a,xi,
     +           vert, vval, diagonal, trL, delta1, delta2, setLF)
            if ( j.eq.0 ) then
                trace_hat_out = trL
                one_delta = delta1
                two_delta = delta2
            end if
            fitted_values = fv
            fitted_residuals = y - fv
            if ( j.lt.iterations ) then
                call lowesw(fitted_residuals,n,robust,temp)
            end if
        end do
        if ( iterations.gt.0 ) then
            call lowesp(n, y, fv, weights, robust, temp,
     +           pv)
            call loess_raw(pv, x, weights, weights, d, n,
     +           max_kd, span,degree, nonparametric, order_drop_sqr,
     +           sum_drop_sqr,new_cell, surf_stat, temp, param_tmp,
     +           a_tmp, xi_tmp,vert_tmp, vval_tmp, diag_tmp, trL_tmp,
     +           d1_tmp,d2_tmp, .false.)
            fitted_values = fv
            pseudovalues = pv
            pseudo_resid = pv - temp
        end if
        sum_squares = 0
        if ( iterations.eq.0 ) then
            do i=1,n
                sum_squares = sum_squares + weights(i)
     +               *fitted_residuals(i)**2
            end do
        else
            do i=1,n
                sum_squares = sum_squares + weights(i)
     +               *pseudo_resid(i)**2
            end do
        end if
        enp = one_delta + 2*trace_hat_out - n
        s = sqrt(sum_squares/one_delta)
        diagonalout = diagonal
        robustout = robust
        divisorout = divisor

        deallocate(parameter)
        deallocate(a,xi,vert,vval)

        deallocate(x,x_tmp,temp,xi_tmp,vert_tmp,vval_tmp,diag_tmp)
        deallocate(a_tmp,param_tmp,order_parametric,order_drop_sqr)
        if ( iterations.gt.0 ) deallocate(pseudo_resid)
        deallocate(robust)

        end subroutine

