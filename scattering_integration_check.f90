program int_check
    use physics_functions
    use numerical_methods
    use h_kr_functions
    implicit none  

    real * 8 :: h ,mh,mkr,hbar,reduced_mass,epsilon,epsilonprime,factor,rho ,energy,wavenumber,wavelength,pi,u0,u1,u_we_care(2),phase,spherj1(1),spherj2(1),sphern1(1),sphern2(1),num,den,frac
    integer :: l,i,index,index_2
    real * 8,allocatable::u_array(:),V_array(:),x_array(:),g_array(:),ustartarray(:,:),vtilde_array(:)
    character*100::string
    procedure (two_pointing_func),pointer:: V_ptr => lennard_jones_potential_thijssen

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
    h = 1.0e-7_dp 
    energy = 0.1_dp/factor
    wavenumber = sqrt(2.0_dp*energy) 
    wavelength = 2.0_dp * pi / wavenumber

    x_array = my_arange(0.5_dp,5.0_dp + wavelength,h)
    V_array = create_potential_array_for_lennard(V_ptr,x_array,epsilonprime,5.0_dp)
    l = 0
    vtilde_array = V_array + create_centrifugal_barrier(x_array,l)
    g_array = energy - vtilde_array 
    g_array = h*h * g_array / 6.0_dp 
    frac = 0.75_dp
    ustartarray = numerov_start_lennard(epsilonprime,0.5_dp,energy,h,l)

    allocate(u_array(size(x_array)))

    u0 = ustartarray(1,l+1)
    u1 = ustartarray(2,l+1)
    u_array(1) = u0 
    u_array(2) = u1

    index = find_closest_index(x_array,5.0_dp)
    index_2 = find_closest_index(x_array,5.0_dp + frac*wavelength)

    spherj1 = spherical_j(l,wavenumber*x_array(index))
    sphern1 = spherical_n(l,wavenumber*x_array(index))
    spherj2 = spherical_j(l,wavenumber*x_array(index_2))
    sphern2 = spherical_n(l,wavenumber*x_array(index_2))

    print*, "INTEGRATION STEP", h

    print*,"kx values"
    print*,wavenumber*x_array(index)
    print*,wavenumber*x_array(index_2)

    print*," Bessel values:"
    print*, spherj1    
    print*, sphern1
    print*, spherj2
    print*, sphern2

    do i = 3,size(x_array)
        u_array(i) = next_numerov(u_array(i-1),u_array(i-2),g_array(i),g_array(i-1),g_array(i-2))

        if (i == index) then 
            print*,"aiming for:",5.0_dp 

            print*, "checking:", x_array(i)
            u_we_care(1) = u_array(i)
        end if 
        if (i == index_2) then 
            print*,"aiming for:",5.0_dp + frac*wavelength
            print*, "checking:", x_array(i)
            u_we_care(2) = u_array(i)
        end if 
    end do 

    phase = calculate_k_new(u_we_care,x_array(index),x_array(index_2))
    print*,"cap k",phase
    num = (phase*spherj1(l+1)-spherj2(l+1))
    den = (phase*sphern1(l+1)-sphern2(l+1))
    print*,"num,den"
    print*, num,den
    phase = num/den
    print*,"tanphase",phase
    phase = ATAN(phase)

    print*,"phase,",phase

    WRITE(string,'(a4,F20.6,a4,F20.2,a12)') "step",h ,"frac",frac,"intcheck.dat"
    call StripSpaces(string)
    print*,"writing to ",string
    open (1, file = string)
    do i = 1, size(x_array)
        if (modulo(i,100) == 0) then 
            write(1,FMT,advance = "no") x_array(i)
            write(1,fmt="(1x,a)",advance = "no") " "
            write(1,FMT,advance = "no") u_array(i) 
            write(1,fmt="(1x,a)",advance = "no") " "
            write(1,FMT,advance = "no") vtilde_array(i) 
            write(1,*) " " 
        end if 

    end do 
end program