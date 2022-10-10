module physics_functions
    
    !all of the potential functions as well as production of possible schrodinger solutions.
    use numerical_methods
    use, intrinsic :: iso_fortran_env, only: wp=>real64
    use, intrinsic :: ieee_arithmetic
    implicit none 
    !integer, parameter :: dp = selected_real_kind(15, 307)

    abstract interface 
        !For single argument function pointers
        function pointing_func(x)
            !integer, parameter :: dp = selected_real_kind(15, 307)
            Import::dp
            real*8,intent (in)::x
            real*8::pointing_func
        end function pointing_func

        function two_pointing_func(x,y)
            !integer, parameter :: dp = selected_real_kind(15, 307)
            Import::dp
            real*8,intent (in)::x,y
            real*8::two_pointing_func
        end function two_pointing_func

    end interface


contains

    function my_linspace(x_0,x_last,N)
        !makes an evenly spaced array of real numbers betweenn x_0 and x_last of length N
        real*8, intent(in)::x_0,x_last
        integer::N
        real*8::h
        real*8::my_linspace(N)

        h = (x_last-x_0)/(N-1)

        my_linspace = my_arange(x_0,x_last,h)

    end function my_linspace

    function my_arange(x_0,x_last,h) result (rangearray)
        
        real*8, intent(in)::x_0,x_last,h
        integer::N,i
        real*8,allocatable::rangearray(:)
        
        N = NINT((x_last-x_0)/h)+ 1
        allocate(rangearray(N))
        
        do i = 1,N
            rangearray(i) = x_0 + (i-1)*h
        end do

    end function my_arange

    function create_effective_potential_array(V_ptr,x_array,l) result(Veff)
        !for a normal potential either write a new function or set l = 0
        procedure (pointing_func),pointer:: V_ptr

        real*8,intent(in)::x_array(:)
        integer,intent(in)::l
        real*8::Veff(size(x_array))
        integer::i
        do i=1,size(x_array)
            if (x_array(i) .ne. 0.0_dp) then 
                Veff(i) = V_ptr(x_array(i)) + centrifugal_barrier(x_array(i),l)
            else
                Veff(i) = ieee_value(1.0_wp, ieee_positive_inf)
            end if 
        end do

    end function create_effective_potential_array

    function create_effective_potential_array_for_lennard(V_ptr,x_array,l,epsilon,r_max) result(Veff)
        !for a normal potential either write a new function or set l = 0
        procedure (two_pointing_func),pointer:: V_ptr

        real*8,intent(in)::x_array(:),epsilon,r_max
        integer,intent(in)::l
        real*8::Veff(size(x_array))
        integer::i
        do i=1,size(x_array)
            Veff(i) =  centrifugal_barrier(x_array(i),l)

            if (x_array(i) < r_max) then 
                !apparantly we must cut off... according to the book
                Veff(i) = Veff(i) + V_ptr(x_array(i),epsilon)  
            end if 
        end do

    end function create_effective_potential_array_for_lennard




    function centrifugal_barrier(x,l)
        !calculates the centrifugal barrier for radial problems
        real*8,intent(in)::x
        integer,intent(in)::l
        real*8::centrifugal_barrier,factor
        
        factor = (l+1.0_dp)*l

        centrifugal_barrier = factor/(2.0_dp*x*x)

    end function centrifugal_barrier

    function coulomb_potential(r) result(V) 
        !calculates the coulomb potential IN ATOMIC UNITS e = 4\pi\epsilon_0 = 1
        implicit none
        real*8, intent(in) :: r
        real*8 :: V
        V = -1.0_dp/r
    end function coulomb_potential

    function zero_potential(x)
        !free particle
        implicit none
        real*8, intent(in)::x
        real*8::zero_potential
        zero_potential = 0.0_dp*x
    end function zero_potential

    function lennard_jones_potential_wiki(r) result(V)
        !lennard jones
        implicit none
        real*8, intent(in) :: r
        real*8 :: V
        V = 1.0_dp/(r**12) - 1.0_dp/(r**6)
    end function lennard_jones_potential_wiki

    function lennard_jones_potential_thijssen(r,epsilonprime) result(V)
        !lennard jones
        implicit none
        real*8, intent(in) :: r,epsilonprime
        real*8 :: V

        V  = (1.0_dp/r)**6
        V = epsilonprime*(V**2 - 2*V)

    end function lennard_jones_potential_thijssen


    function energy_minus_potential_array(E,V_array)
        !makes E - V(r)
        implicit none 
        real*8, intent(in)::V_array(:)
        real*8, intent(in)::E
        real*8::energy_minus_potential_array(size(V_array))
        
        energy_minus_potential_array = E-V_array

    end function energy_minus_potential_array
    
    

    

    function assymptotic_solution(x_max,V_of_x_max,E)
        !assymptotic solution, see Chapter 9 of Griffiths
        real*8,intent(in)::x_max,V_of_x_max,E
        real*8::assymptotic_solution
        
        assymptotic_solution = exp(-abs(x_max)*sqrt(2.0_dp*(V_of_x_max-E)))
        !print*,assymptotic_solution
    end function assymptotic_solution
    
    function find_trial_solution_from_right(x_array,V_array,E) Result(psi_array)
        real*8,intent(in)::x_array(:),V_array(:),E
        real*8::psi_array(size(x_array)),x_max,x_max_minus_one,vmax,vmax_minus_one,psi_bound_one,psi_bound_two
        integer::n
        
        n = size(x_array)
        x_max = x_array(n)
        x_max_minus_one = x_array(n-1)
        vmax_minus_one = V_array(n-1)
        vmax = V_array(n)
        !print*,V_array(n)
        psi_bound_one = assymptotic_solution(x_max,vmax,E)
        psi_bound_two = assymptotic_solution(x_max_minus_one,vmax_minus_one,E)
        
        !print*,psi_bound_one,psi_bound_two

        psi_array = numerov_schrodinger(x_array,psi_bound_one,psi_bound_two,V_array,E,.true.)

        !print*,(psi_array(1))
    end function find_trial_solution_from_right

    function find_trial_solution_from_left(x_array,V_array,E,psi_0,psi_1) Result(psi_array)
        real*8,intent(in)::x_array(:),V_array(:),psi_0,psi_1,E
        real*8::psi_array(size(x_array))
        psi_array = numerov_schrodinger(x_array,psi_0,psi_1,V_array,E,.false.)

    end function find_trial_solution_from_left

    


    function r_overlap(r,psi_one,psi_two) result(overlap)
        real*8,intent(in)::r(:),psi_one(:),psi_two(:)
        real*8::overlap,f(size(r))

        f = r*psi_one*psi_two
        !print*,f
        overlap = integrate_trapezium(r(2)-r(1),f)


    end function r_overlap

end module  physics_functions
