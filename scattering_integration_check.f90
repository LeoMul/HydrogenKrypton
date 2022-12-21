program int_check
    use physics_functions
    use numerical_methods
    use h_kr_functions
    implicit none  

    real * 8 :: h ,mh,mkr,hbar,reduced_mass,epsilon,epsilonprime,factor,rho ,energy,wavenumber,wavelength,pi,u0,u1
    integer :: l,i
    real * 8,allocatable::u_array(:),V_array(:),x_array(:),g_array(:),ustartarray(:,:),vtilde_array(:)
    procedure (two_pointing_func),pointer:: V_ptr => lennard_jones_potential_thijssen
    class(starting_conditions_for_lennard_jones),allocatable ::ustruct

    pi = 3.14159265359_dp

    mh = 1.6735575e-27_dp !kg
    mkr = 1.3914984e-25_dp!kg
    hbar = 1.05457182e-34_dp!mks
    reduced_mass = mh*mkr/(mkr+mh)
    
    epsilon = 5.90_dp !meV
    rho = 3.57_dp !in agstroms
    factor = hbar**2 / (reduced_mass * (rho*1e-10_dp)**2)!mks
    factor = factor * 6.241506363094e+21
    epsilonprime = epsilon/factor
    h = 1.0e-5_dp 
    energy = 0.1_dp/factor
    wavenumber = sqrt(2.0_dp*energy) 
    wavelength = 2.0_dp * pi / wavenumber

    x_array = my_arange(0.5_dp,5.0_dp + 0.5_dp * wavelength,h)
    V_array = create_potential_array_for_lennard(V_ptr,x_array,epsilonprime,5.0_dp)
    l = 0
    vtilde_array = V_array + create_centrifugal_barrier(x_array,l)
    g_array = energy - vtilde_array 
    g_array = h*h * g_array / 6.0_dp 

    ustartarray = numerov_start_lennard(epsilonprime,0.5_dp,energy,h,l)

    allocate(u_array(size(x_array)))

    u0 = ustartarray(1,l+1)
    u1 = ustartarray(2,l+1)
    u_array(1) = u0 
    u_array(2) = u1
    do i = 3,size(x_array)
        u_array(i) = next_numerov(u_array(i-1),u_array(i-2),g_array(i),g_array(i-1),g_array(i-2))
    end do 
    open (1, file = "integrationcheck.dat")
    do i = 1, size(x_array)
        write(1,FMT,advance = "no") x_array(i)
        write(1,fmt="(1x,a)",advance = "no") " "
        write(1,FMT,advance = "no") u_array(i) 
        write(1,fmt="(1x,a)",advance = "no") " "
        write(1,FMT,advance = "no") vtilde_array(i) 
        write(1,*) " "

    end do 
end program