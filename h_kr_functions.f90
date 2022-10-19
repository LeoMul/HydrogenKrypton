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
        real*8::psi_array(size(x_array)),psi_0,psi_1,h,hsquared_over_twelve,e_minus_potential_array(size(x_array))
        class(starting_conditions_for_lennard_jones),allocatable ::ustruct
        integer::i
        ustruct = small_r_solution_for_lennard_jones(x_array(2),epsilon_prime)

        psi_0 = ustruct%u 

        h = x_array(2)-x_array(1)
        n = size(x_array)
        e_minus_potential_array = 2.0_dp*(E - V_array)
        !need the factor of two in conventional atomic units.
        hsquared_over_twelve = h*h/12.0_dp




        psi_1 = numerov_get_step_two(-e_minus_potential_array(2),-e_minus_potential_array(1),-e_minus_potential_array(3),h,hsquared_over_twelve,psi_0,ustruct%udot)
        !print*,"psi 1", psi_1
        psi_array(1) = psi_0
        psi_array(2) = psi_1
        do i = 2,n-1,1
            psi_array(i+1) = numerov_next_step_s_zero(hsquared_over_twelve,psi_array(i),psi_array(i-1),e_minus_potential_array(i+1),e_minus_potential_array(i),e_minus_potential_array(i-1))
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

    function calculate_correction_delta_l(x_array,l,epsilon_prime,wavenumber) result(correction)
        
        real*8,intent(in)::x_array(:),epsilon_prime,wavenumber
        integer,intent(in)::l
        real*8::correction,v_array(size(x_array)) ,j_array(size(x_array)),temp(l+1),h
        integer::i
        procedure (two_pointing_func),pointer:: V_ptr => lennard_jones_potential_thijssen

        v_array = create_potential_array_for_lennard(V_ptr,x_array,epsilon_prime,x_array(size(x_array)))

        j_array = wavenumber*x_array
        do i = 1,size(x_array)
            temp = spherical_j(l,j_array(i))
            j_array(i) = temp(size(temp))
        end do
        h = x_array(2) - x_array(1)
        j_array = (j_array**2)*v_array*(x_array**2)
        correction = -2.0_dp*wavenumber*integrate_trapezium(h,j_array)
    end function calculate_correction_delta_l

    function calculate_cross_section_sum(x_array,v_array,l_array,energy_prime,epsilon_prime,pos,pos_2,bessel_j_1,bessel_j_2,bessel_n_1,bessel_n_2) result(sum)
        real*8,intent(in)::x_array(:),energy_prime,epsilon_prime,v_array(:)
        integer,intent(in)::l_array(:),pos,pos_2
        real*8::psi_array(size(x_array)),k,tandelta,delta,sum,bessel_j_1(:),bessel_j_2(:),bessel_n_1(:),bessel_n_2(:)
        integer::i ,l 

        do i =  1,size(l_array)
            l = l_array(i)
            print*,"ang",l
            psi_array = integrate_lennard(x_array,v_array+create_centrifugal_barrier(x_array,l),energy_prime,epsilon_prime)
            k = calculate_k(x_array,psi_array,pos,pos_2)
            tandelta = (k*bessel_j_1(l+1)-bessel_j_2(l+1))/(k*bessel_n_1(l+1)-bessel_n_2(l+1))
            delta = ATAN(tandelta)
            !if (correction) then 
            !    delta = delta + calculate_correction_delta_l(big_x_array,l,epsilon_prime,wavenumber)
            !end if 
            sum = sum + (2*l+1)*(sin(delta)**2)
        end do


    end function calculate_cross_section_sum

    function create_l_array(l_max) result(l_array)
        integer,intent(in)::l_max 
        integer::l_array(l_max+1),i
        do i = 1,size(l_array)
            l_array(i) = i-1
        end do

    end function create_l_array



end module h_kr_functions
