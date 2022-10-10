module h_kr_functions
    use numerical_methods
    use physics_functions
    use, intrinsic :: iso_fortran_env, only: wp=>real64
    use, intrinsic :: ieee_arithmetic

    type starting_conditions_for_lennard_jones
        real*8:: u 
        real*8:: udot 
    end type


    contains 

    function small_r_solution_for_lennard_jones(r,epsilon_prime) result(ustruct)
        real*8,intent(in)::r,epsilon_prime
        real*8::factor,r5
        class(starting_conditions_for_lennard_jones),allocatable::ustruct
        allocate(starting_conditions_for_lennard_jones::ustruct)

        factor = sqrt(2.0_dp*epsilon_prime)
        r5 = r**(5)
        !print*, r5
        !print*, "exp", -factor/(5.0_dp*r5)
        ustruct%u = exp(-factor/(5.0_dp*r5) )

        ustruct%udot = factor *r5*(ustruct%u)/r
        

    end function small_r_solution_for_lennard_jones

    function integrate_lennard(x_array,V_array,E,epsilon_prime) result(psi_array)
        real*8,intent(in)::x_array(:),V_array(:),epsilon_prime,E
        real*8::psi_array(size(x_array)),psi_0,psi_1,h,hsquared_over_twelve,energy_minus_potential_array(size(x_array))
        class(starting_conditions_for_lennard_jones),allocatable ::ustruct
        integer::i
        ustruct = small_r_solution_for_lennard_jones(x_array(2),epsilon_prime)

        psi_0 = ustruct%u 

        h = x_array(2)-x_array(1)
        n = size(x_array)
        energy_minus_potential_array = 2.0_dp*(E - V_array)
        !need the factor of two in conventional atomic units.
        hsquared_over_twelve = h*h/12.0_dp




        psi_1 = numerov_get_step_two(-energy_minus_potential_array(2),-energy_minus_potential_array(1),-energy_minus_potential_array(3),h,hsquared_over_twelve,psi_0,ustruct%udot)
        !print*,"psi 1", psi_1
        psi_array(1) = psi_0
        psi_array(2) = psi_1
        do i = 2,n-1,1
            psi_array(i+1) = numerov_next_step_s_zero(hsquared_over_twelve,psi_array(i),psi_array(i-1),energy_minus_potential_array(i+1),energy_minus_potential_array(i),energy_minus_potential_array(i-1))
            if (psi_array(i+1)>10000.0_dp) then 
                psi_array = psi_array/10000.0_dp
            end if
        end do
        
        !psi_array = numerov_schrodinger(x_array,psi_0,psi_1,V_array,E,.false.)




    end function integrate_lennard 

    function spherical_j(n,x) result(j)
        !returns array of size n+1, that has entires j_0 to j_n
        real*8,intent(in)::x 
        real*8:: j_0,j_1
        integer,intent(in)::n 
        real*8::j(n+1)

        j_0 = sin(x)/x 
        j_1 = sin(x)/(x**2) - cos(x)/x

        if (n == 0) then 
            j(1) = j_0
        else if (n == 1) then 
            j(1) = j_0
            j(2) = j_1
        else
            j = spherical_recursive(n,x,j_0,j_1)
        end if 
    end function

    function spherical_n(n,x) result(j)
        !returns array of size n+1, that has entires h_0 to h_n

        real*8,intent(in)::x 
        real*8:: j_0,j_1
        integer,intent(in)::n 
        real*8::j(n+1)

        j_0 = -cos(x)/x 
        j_1 = -cos(x)/(x**2) - sin(x)/x

        if (n == 0) then 
            j(1) = j_0
        else if (n == 1) then 
            j(1) = j_0
            j(2) = j_1
        else
            j = spherical_recursive(n,x,j_0,j_1)
        end if 
    end function

    function spherical_recursive(n,x,j_0,j_1) result(j)
        real*8,intent(in)::x 
        real*8,intent(in):: j_1,j_0
        integer,intent(in)::n 
        integer::i
        real*8::j(n+1)

        j(1) = j_0 
        j(2) = j_1
        do i = 3,n+1 
            j(i) = (2*i-3)*j(i-1)/x -j(i-2)
        end do


    end function spherical_recursive

    function calculate_k(x_array,u_array,index_1,index_2) result(k)
        real*8,intent(in)::x_array(:),u_array(:)
        integer,intent(in)::index_1,index_2
        real*8::k 

        k = x_array(index_1)*u_array(index_2)
        k = k / (x_array(index_2)*u_array(index_1))
    end function

    function find_closest_index(x_array,desired_x) result(index)
        real*8,intent(in)::x_array(:),desired_x
        integer::index
        real*8::new_x(size(x_array))
        new_x = abs(x_array-desired_x)
        index = minloc(new_x,dim = 1)

    end function find_closest_index

end module h_kr_functions
