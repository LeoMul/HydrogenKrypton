module numerical_methods
    use, intrinsic :: iso_fortran_env, only: wp=>real64
    use, intrinsic :: ieee_arithmetic
    implicit none
    public
    integer, parameter :: dp = selected_real_kind(15, 307)
    
    !doube: 15,307 :::::: quad: 33,4931
    CHARACTER(LEN=20) :: FMT = "(F20.12)"
contains



    function numerov_next_step(hsquared_over_twelve,ycurrent,yprev,gnext,gcurrent,gprev,sprev,scurrent,snext)
        !includes s(x) function, see wikipedia article
        !included for completeness, not really needed
        real*8::numerov_next_step,factor,s_term
        real*8, intent (in) ::hsquared_over_twelve,ycurrent,yprev,gnext,gcurrent,gprev,sprev,scurrent,snext

        !calculates y_{n+1}. It also calcualtes y_{n-1} if y_prev,g_prev are actually y_{n+1} and g_{n+1}
        !as well as gnext = g_{n-1}

        s_term = hsquared_over_twelve*(sprev + 10.0_dp*scurrent+snext)
        factor = (1.0_dp+hsquared_over_twelve*gnext)

        numerov_next_step = 2.0_dp*ycurrent*(1.0_dp-5.0_dp*hsquared_over_twelve*gcurrent)
        numerov_next_step = numerov_next_step-yprev *(1.0_dp+hsquared_over_twelve*gprev) + s_term
        numerov_next_step = numerov_next_step/factor

    end function numerov_next_step
 
    function numerov_next_step_s_zero(hsquared_over_twelve,ycurrent,yprev,gnext,gcurrent,gprev)
        !assumes the s(x) function is zero, see wikipedia article
        real*8::numerov_next_step_s_zero,factor
        real*8, intent (in) ::hsquared_over_twelve,ycurrent,yprev,gnext,gcurrent,gprev

        !calculates y_{n+1}. It also calcualtes y_{n-1} if y_prev,g_prev are actually y_{n+1} and g_{n+1}
        !as well as gnext = g_{n-1}

        !also allows for infinities provided they are well behaved, more of an assumption than a requirement.

        factor = (1.0_dp+hsquared_over_twelve*gnext)
        numerov_next_step_s_zero = 2.0_dp*ycurrent*(1.0_dp-5.0_dp*hsquared_over_twelve*gcurrent)
        !if (yprev .ne. 0.0_dp) then 
        numerov_next_step_s_zero = numerov_next_step_s_zero-yprev *(1.0_dp+hsquared_over_twelve*gprev)
        !end if 
        numerov_next_step_s_zero = numerov_next_step_s_zero/factor

    end function numerov_next_step_s_zero

    function numerov_schrodinger(x_array,y_bound_one,y_bound_two,potential_array,E,from_right)

        real*8, intent(in) :: x_array(:),potential_array(:),E,y_bound_one,y_bound_two
        real*8:: numerov_schrodinger(size(x_array)),h,energy_minus_potential_array(size(potential_array)),hsquared_over_twelve
        integer::n,i
        logical::from_right 
        
        !if from_right is true, we begin integration at x_n and intepret the bounds as y_bound_one = y_n and y_bound_two = y_{n-1}
        !otherwise, we intepreat the bounds y_bound_one = y_0 and y_bound_two = y_1

        !some stuff we will need over and over again. assumes x is evenly spaced.
        !print*,potential_array(2)
        h = x_array(2)-x_array(1)
        n = size(x_array)
        energy_minus_potential_array = 2.0_dp*(E - potential_array)!need the factor of two in conventional atomic units.
        hsquared_over_twelve = h*h/12.0_dp

        if (from_right) then
            !do the integration from the right
            numerov_schrodinger(n) = y_bound_one
            numerov_schrodinger(n-1) = y_bound_two
            do i = n-1,2,-1
                !numerov_next_step_s_zero(hsquared_over_twelve,numerov_schrodinger(i),numerov_schrodinger(i+1),energy_minus_potential_array(i-1),energy_minus_potential_array(i),energy_minus_potential_array(i+1))
                numerov_schrodinger(i-1) = numerov_next_step_s_zero(hsquared_over_twelve,numerov_schrodinger(i),numerov_schrodinger(i+1),energy_minus_potential_array(i-1),energy_minus_potential_array(i),energy_minus_potential_array(i+1))
            end do
        else 
            !do the integration from the left
            numerov_schrodinger(1) = y_bound_one
            numerov_schrodinger(2) = y_bound_two
            do i = 2,n-1,1
                numerov_schrodinger(i+1) = numerov_next_step_s_zero(hsquared_over_twelve,numerov_schrodinger(i),numerov_schrodinger(i-1),energy_minus_potential_array(i+1),energy_minus_potential_array(i),energy_minus_potential_array(i-1))
            end do

        end if



    end function numerov_schrodinger



    
    function find_trial_solution_normalise(x_array,V_array,E,psi_bound_one,psi_bound_two,h,from_right)

        real*8,intent(in)::x_array(:),V_array(:),E,psi_bound_one,psi_bound_two,h
        real*8::find_trial_solution_normalise(size(x_array))
        logical::from_right

        find_trial_solution_normalise = numerov_schrodinger(x_array,psi_bound_one,psi_bound_two,V_array,E,from_right)


        call normalise_wave_function(find_trial_solution_normalise,h)

    end function

    subroutine normalise_wave_function(psi_array,h)
        implicit none
        !Normalises the wave function by the sqrt. of int_{x_min}^{x_max} 
        real*8, intent(in)::h
        real*8::norm_factor
        real*8, intent(inout)::psi_array(:)
        real*8, dimension(size(psi_array))::psi_sq
        integer::i
        psi_sq = psi_array*psi_array
        norm_factor = SQRT(integrate_trapezium(h,psi_sq))

        do i = 1,size(psi_array)
            psi_array(i) = psi_array(i)/norm_factor
        end do

    end subroutine normalise_wave_function


    function integrate_trapezium(h,f_array)
        !Integrates with trapezium rule
        implicit none
        real*8, intent(in)::h,f_array(:)
        real*8::integrate_trapezium
        integer::N,i
        N = size(f_array)
        integrate_trapezium = 0.0_dp
        do i = 2,N-1
            integrate_trapezium = integrate_trapezium + f_array(i)
        end do
        integrate_trapezium = integrate_trapezium + 0.5_dp*(f_array(1)+f_array(N))
        integrate_trapezium = integrate_trapezium*h
    end function integrate_trapezium

    function numerical_derivative_two_point(f1,f2,x2,x1) result (nd)
        !approximates the derivative of f(x) at x1
        real*8, intent(in)::f1,f2,x2,x1
        real*8 :: nd
        nd = (f2-f1)/(x2-x1)

    end function numerical_derivative_two_point

    function numerov_get_step_two(f0,fmh,fh,h,hsquared_over_twelve,x0,xdot0) result(x_h)
        !according to thijssen ch2, we cant use the small r solution for both points at
        !the start, so we use the derivative of the first point and the actual
        !ode to get an approximate solution to the second point
        !to get us started.
        real*8,intent(in)::f0,fmh,fh,hsquared_over_twelve,x0,xdot0,h
        real*8::x_h,fac1,fac2,fac3

        fac1 = hsquared_over_twelve * fmh
        fac2 = hsquared_over_twelve * f0
        fac3 = (1.0_dp + 2.0_dp*fac1)*(1.0_dp + hsquared_over_twelve*fh)
        fac3 = fac3 + (1.0_dp - fac1)*(1.0_dp + 2.0_dp * hsquared_over_twelve*fh)
        x_h = (1.0_dp - fac1) * h * xdot0 + (1.0_dp + 2.0_dp * fac1)*(1.0_dp - 5.0_dp * fac2)*x0
        
        x_h = 2.0_dp * x_h
        
        x_h = x_h / fac3

    end function numerov_get_step_two




end module numerical_methods