program potentials 
    use physics_functions
    use numerical_methods
    use h_kr_functions
    implicit none

    real*8::rho,epsilon,h,x_last,x_min,epsilon_prime,factor
    real*8::reduced_mass,hbar,mh,mkr
    real*8,allocatable::x_array(:),V_array(:),potential_matrix(:,:)
    integer::l,i,J,l_max,l_min
    procedure (two_pointing_func),pointer:: V_ptr => lennard_jones_potential_thijssen
    class(starting_conditions_for_lennard_jones),allocatable ::ustruct
    
    allocate(ustruct)   
    !position stuff

    !in units of rho
    x_min = 0.3_dp
    x_last = 2.5_dp
    h = 1.0e-3_dp 
    x_array = my_arange(x_min,x_last,h)

    mh = 1.6735575e-27_dp !kg
    mkr = 1.3914984e-25_dp!kg
    hbar = 1.05457182e-34_dp!mks

    reduced_mass = mh*mkr/(mkr+mh)
    


    !Parameters and stuff
    l = 0

    epsilon = 5.90_dp !meV
    rho = 3.57_dp !in agstroms

    factor = hbar**2 / (reduced_mass * (rho*1e-10_dp)**2)!mks
    !print*,factor 
    factor = factor * 6.241506363094e+21 !converts to meV
    epsilon_prime = epsilon/factor
    l_min = 0
    l_max = 6
    allocate(potential_matrix(size(x_array),l_max+1))
    do l = 0,l_max

        V_array = create_effective_potential_array_for_lennard(v_ptr,x_array,l,epsilon_prime)
        potential_matrix(:,l+1) = V_array*factor
    end do 
    
    open (1, file = "potentials.dat")

    do i = 1,size(x_array)
                write(1,FMT,advance = "no") x_array(i)
                write(1,fmt="(1x,a)",advance = "no") " "
                do j = 1,l_max+1
                    write(1,FMT,advance = "no") potential_matrix(i,j)
                    write(1,fmt="(1x,a)",advance = "no") " "
 
                end do 
                write(1,*) " "

    end do 

end program potentials