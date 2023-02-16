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

    function integrate_lennard(x_array,V_array,E,epsilon_prime,l) result(psi_array)
        real*8,intent(in)::x_array(:),V_array(:),epsilon_prime,E
        integer, intent(in):: l
        real*8::psi_array(size(x_array)),psi_0,psi_1,h,hsquared_over_twelve,e_minus_potential_array(size(x_array)),fmh,xmh
        class(starting_conditions_for_lennard_jones),allocatable ::ustruct
        integer::i
        ustruct = small_r_solution_for_lennard_jones(x_array(2),epsilon_prime)

        psi_0 = ustruct%u 

        h = x_array(2)-x_array(1)
        n = size(x_array)
        e_minus_potential_array = 2.0_dp*(E - V_array)
        !need the factor of two in conventional atomic units.
        hsquared_over_twelve = h*h/12.0_dp
        xmh = x_array(1) - h
        fmh = lennard_jones_potential_thijssen(xmh,epsilon_prime) + centrifugal_barrier(xmh,l) -E
        fmh = 2.0_dp * fmh 

        !this line is technically incorrect
        psi_1 = numerov_get_step_two(-e_minus_potential_array(1),fmh,-e_minus_potential_array(2),h,hsquared_over_twelve,psi_0,ustruct%udot)
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

    function calculate_k_new(uofinterest,xmax,xlast) result(k)
        real * 8 ,intent(in) :: uofinterest(2),xmax,xlast
        real * 8 :: k 

        k = xmax * uofinterest(2)
        k = k / (xlast * uofinterest(1))


    end function 


    function find_closest_index(x_array,desired_x) result(index)
        real*8,intent(in)::x_array(:),desired_x
        integer::index,j
        !new_x = abs(x_array-desired_x)
        !index = minloc(new_x,dim = 1)
        do j = size(x_array),1,-1
            if (x_array(j)-desired_x < 0.0_dp) then 
                index = j + 1
                EXIT
            end if 
        end do 



    end function find_closest_index

    function calculate_correction_delta_l(r_max,l,epsilon_prime,wavenumber,h) result(correction)
        
        real*8,intent(in)::epsilon_prime,wavenumber,r_max,h
        integer,intent(in)::l
        real * 8,allocatable :: v_array(:),j_array(:),x_array(:)
        real*8::correction,temp(l+1)
        integer::i
        procedure (two_pointing_func),pointer:: V_ptr => lennard_jones_potential_thijssen
        
        x_array = my_arange(r_max,3.0_dp*r_max,h)
        
        v_array = create_potential_array_for_lennard(V_ptr,x_array,epsilon_prime,4.0_dp*r_max)

        j_array = wavenumber*x_array
        do i = 1,size(x_array)
            temp = spherical_j(l,j_array(i))
            j_array(i) = temp(size(temp))
        end do
        j_array = (j_array**2)*v_array*(x_array**2)
        correction = -2.0_dp*wavenumber*integrate_trapezium(h,j_array)
    end function calculate_correction_delta_l

    function calculate_phase_shift_l_term(x_array,v_array,l,energy_prime,epsilon_prime,pos,pos_2,bessel_j_1,bessel_j_2,bessel_n_1,bessel_n_2) result(delta)
        real*8,intent(in)::x_array(:),energy_prime,epsilon_prime,v_array(:)
        integer,intent(in)::l,pos,pos_2
        real*8::psi_array(size(x_array)),k,tandelta,delta,bessel_j_1(:),bessel_j_2(:),bessel_n_1(:),bessel_n_2(:)

        psi_array = integrate_lennard(x_array,v_array+create_centrifugal_barrier(x_array,l),energy_prime,epsilon_prime,l)
        k = calculate_k(x_array,psi_array,pos,pos_2)
        tandelta = (k*bessel_j_1(l+1)-bessel_j_2(l+1))/(k*bessel_n_1(l+1)-bessel_n_2(l+1))
        delta = ATAN(tandelta)
        

    end function calculate_phase_shift_l_term

    function calculate_cross_sec_term(l,delta) result (term)
        real*8,intent(in):: delta 
        integer,intent(in) :: l 
        real*8 :: term 

        term = (2*l+1) * (sin(delta)**2)

    end function calculate_cross_sec_term




    function create_l_array(l_max) result(l_array)
        integer,intent(in)::l_max 
        integer::l_array(l_max+1),i
        do i = 1,size(l_array)
            l_array(i) = i-1
        end do

    end function create_l_array

    subroutine StripSpaces(string)
        !https://stackoverflow.com/questions/27179549/removing-whitespace-in-string
        character(len=*) :: string
        integer :: stringLen 
        integer :: last, actual
    
        stringLen = len (string)
        last = 1
        actual = 1
    
        do while (actual < stringLen)
            if (string(last:last) == ' ') then
                actual = actual + 1
                string(last:last) = string(actual:actual)
                string(actual:actual) = ' '
            else
                last = last + 1
                if (actual < last) &
                    actual = last
            endif
        end do
    
    end subroutine

    function next_numerov(ucurrent,uprev,gnext,gcurrent,gprev) result(unext)
        real*8,intent(in) :: ucurrent,uprev,gnext,gcurrent,gprev
        real*8 :: unext 

        unext = 2.0_dp * ucurrent * (1.0_dp - 5.0_dp * gcurrent) - uprev * (1.0_dp + gprev) 
        unext = unext / (1.0_dp + gnext)

    end function 

    function integrate_for_scattering_and_find_u_of_interest(uzero,uone,garray,index_first) result(uofinterest)
        real*8 , intent(in) :: uzero,uone,garray(:)
        integer, intent(in) :: index_first 
        integer :: i,index_second
        real* 8 :: unext, ucurrent, uprev ,uofinterest(2)
        
        index_second = size(garray)
        uprev = uzero
        ucurrent = uone 
        unext = 0.0_dp
        do i = 3,index_first
            unext = next_numerov(ucurrent,uprev,garray(i),garray(i-1),garray(i-2))
            uprev = ucurrent 
            ucurrent = unext 
        end do 
        uofinterest(1) = unext 

        !print*, "you may need to set V = 0 here, eremember."

        do i = index_first+1,index_second 
            unext = next_numerov(ucurrent,uprev,garray(i),garray(i-1),garray(i-2))
            uprev = ucurrent 
            ucurrent = unext 
        end do 
        uofinterest(2) = unext 
    end function 

    function find_phase_shifts(ustartarray,energy,lmax,x_array,v_array,indexfirst,h) result(phases)
        real * 8 ,intent(in) :: energy,x_array(:),v_array(:),h,ustartarray(:,:)
        integer , intent(in) :: lmax ,indexfirst
        real * 8 :: phases(lmax+1) , g_array(size(x_array)),wavevector,xmax,xlast ,spherj1(lmax+1),spherj2(lmax+1),sphern1(lmax+1),sphern2(lmax+1),k
        integer :: l , indexsecond
        real * 8 :: uofinterest(2),delta,tandelta,uzero,uone


        indexsecond = size(v_array)
        wavevector = sqrt(2.0_dp * energy) 

        xmax = x_array(indexfirst)
        xlast = x_array(indexsecond)
        spherj1 = spherical_j(lmax,wavevector*xmax)
        sphern1 = spherical_n(lmax,wavevector*xmax)
        spherj2 = spherical_j(lmax,wavevector*xlast)
        sphern2 = spherical_n(lmax,wavevector*xlast)



        do l = 0,lmax 
            
            g_array = v_array + create_centrifugal_barrier(x_array,l)
            g_array = energy - g_array 
            g_array = h*h * g_array / 6.0_dp 
            uzero = ustartarray(1,l+1)
            uone = ustartarray(2,l+1)
            uofinterest = integrate_for_scattering_and_find_u_of_interest(uzero,uone,g_array,indexfirst)
            k = calculate_k_new(uofinterest,xmax,xlast)
            tandelta = (k*spherj1(l+1)-spherj2(l+1))/(k*sphern1(l+1)-sphern2(l+1))
            delta = ATAN(tandelta)
            !print*,delta
            phases(l+1) = delta 
        end do

    end function 

    function numerov_start_lennard(epsilon_prime,x,energy,h,l_max) result (ustartarray)
        real * 8,intent(in) ::epsilon_prime,x,energy,h
        integer, intent(in) :: l_max
        class(starting_conditions_for_lennard_jones),allocatable ::ustruct
        real * 8 :: udot,u0 ,gm1,g0,gp1 ,ustartarray(2,l_max+1) ,u1,den
        integer :: l 
        ustruct = small_r_solution_for_lennard_jones(x,epsilon_prime)
        u0 = ustruct % u
        udot = ustruct % udot 

        ustartarray(1,:) = u0 

        do l = 0,l_max

            gm1 = energy - lennard_jones_potential_thijssen(x-h,epsilon_prime)- centrifugal_barrier(x-h,l)
            g0 = energy - lennard_jones_potential_thijssen(x,epsilon_prime) - centrifugal_barrier(x,l)
            gp1 = energy - lennard_jones_potential_thijssen(x+h,epsilon_prime) - centrifugal_barrier(x+h,l)
            u1 =  2.0_dp * u0 * (1.0_dp - 5.0_dp * h * h *g0 /6.0_dp) + 2.0_dp * h * udot * (1.0_dp + h*h*gm1/6.0_dp) 

            den = (1.0_dp + h*h*gp1 /6.0_dp) * (1.0_dp + h*h*gm1 /3.0_dp)

            den = den + (1.0_dp + h*h*gm1 /6.0_dp) * (1.0_dp + h*h*gp1 /3.0_dp)

            u1 = u1 / den

            ustartarray(2,l+1) = u1
        end do
        

    end function 

    function numerov_start_screened(h,l_max) result (ustartarray)
        real * 8,intent(in) ::h
        integer, intent(in) :: l_max
        real * 8 :: ustartarray(2,l_max+1) ,u1,u0
        integer :: l 
     

        do l = 0,l_max
            u0 = h**(l+1)
            u1 =(2.0_dp*h)**(l+1)
            ustartarray(1,l+1) = u0
            ustartarray(2,l+1) = u1
        end do
        

    end function 



end module h_kr_functions
