program int_check
    use physics_functions
    use numerical_methods
    use h_kr_functions
    implicit none

    real*8::rho,epsilon,h,energy,x_last,x_min,epsilon_prime,energy_prime,factor,x_to_check
    real*8::delta,tandelta,k,sum,pi,new_x_diff,wavenumber,reduced_mass,hbar,mh,mkr,wavelength,r_max
    real*8,allocatable::x_array(:),big_x_array(:),V_array(:),psi_array(:),bessel_j_1(:),bessel_n_1(:),bessel_j_2(:),bessel_n_2(:),energies(:),cross_secs(:)
    integer::l,i,J,l_max,pos,pos_2,N
    logical::correction
    procedure (two_pointing_func),pointer:: V_ptr => lennard_jones_potential_thijssen
    class(starting_conditions_for_lennard_jones),allocatable ::ustruct
    
    allocate(ustruct)   
    !position stuff
    !in units of rho
    x_min = 0.5_dp
    r_max = 5.0_dp
    x_last = 15_dp !this is for the big integration for the correction term. 
    h = 1.0e-4_dp 
    correction = .false. !calculate corrections using the integral in the book. v slow implemetation right now, set to false.

    big_x_array = my_arange(r_max,x_last,h)
    mh = 1.6735575e-27_dp !kg
    mkr = 1.3914984e-25_dp!kg
    hbar = 1.05457182e-34_dp!mks
    reduced_mass = mh*mkr/(mkr+mh)
    


    !Parameters and stuff
    l = 0

    epsilon = 5.90_dp !meV
    rho = 3.57_dp !in agstroms

    factor = hbar**2 / (reduced_mass * (rho*1e-10_dp)**2)!mks
    factor = factor * 6.241506363094e+21
    print*,factor

    N = 1000
    allocate(energies(N),cross_secs(N))
    energies = my_linspace(0.1_dp,1.0_dp,N)
    epsilon_prime = epsilon/factor

    !print*,"epsilon_prime", epsilon_prime

    l_max  = 6
    x_to_check = 5.0_dp
    pi = 3.14159265359_dp
    
    sum = 0.0_dp

    !print*,"x_array pos",x_array(pos),x_array(pos_2)

    do i = 1,N
        energy = energies(i)
        energy_prime = energy/factor
        print*, "energy in meV, energy in nat",energy,energy_prime
        
        wavenumber = sqrt(2.0_dp*energy_prime) !in rho ^{-1}
        wavelength = 2.0_dp * pi / wavenumber
        new_x_diff = wavelength / 2.0_dp
        x_array = my_arange(x_min,r_max+wavelength,h)

        pos = find_closest_index(x_array,r_max)
        pos_2 = find_closest_index(x_array,r_max+new_x_diff)


        bessel_j_1 = spherical_j(l_max,wavenumber*(x_array(pos)))
        bessel_n_1 = spherical_n(l_max,wavenumber*(x_array(pos)))
        bessel_j_2 = spherical_j(l_max,wavenumber*(x_array(pos_2)))
        bessel_n_2 = spherical_n(l_max,wavenumber*(x_array(pos_2)))
        
        print*,"x to check", x_to_check
        print*,"second x",x_to_check+new_x_diff
        
        
        sum = 0.0_dp
        do l =  0,l_max
            V_array = create_effective_potential_array_for_lennard(v_ptr,x_array,l,epsilon_prime,r_max)
            psi_array = integrate_lennard(x_array,V_array,energy_prime,epsilon_prime)
            k = calculate_k(x_array,psi_array,pos,pos_2)
            tandelta = (k*bessel_j_1(l+1)-bessel_j_2(l+1))/(k*bessel_n_1(l+1)-bessel_n_2(l+1))
            delta = ATAN(tandelta)
            if (correction) then 
                delta = delta + calculate_correction_delta_l(big_x_array,l,epsilon_prime,wavenumber)
            end if 
            sum = sum + (2*l+1)*(sin(delta)**2)
        end do 
        sum = sum * 4.0_dp*pi/(wavenumber**2)
        print*,sum
        cross_secs(i) = sum
    end do

    open (1, file = "cross_sec.dat")
    J = 1 !every jth point
    do i = 1,size(energies)
                IF (modulo(i,J) == 0) THEN 
                    write(1,FMT,advance = "no") energies(i)
                write(1,fmt="(1x,a)",advance = "no") " "
                write(1,FMT,advance = "no") cross_secs(i) 
                write(1,*) " "
            END IF 
    end do 

    !checking exercise A1
    !bessel_j_1 = spherical_j(5,1.5_dp)
    !bessel_n_1 = spherical_n(5,1.5_dp)
    !PRINT*, bessel_j_1
    !PRINT*, bessel_n_1


end program int_check

