program int_check
    use physics_functions
    use numerical_methods
    use h_kr_functions

    USE MPI

    implicit none

    real*8::rho,epsilon,h,energy,x_last,x_min,epsilon_prime,energy_prime,factor,x_to_check
    real*8::sum,pi,new_x_diff,wavenumber,reduced_mass,hbar,mh,mkr,wavelength,r_max,local_sum
    real*8,allocatable::x_array(:),big_x_array(:),V_array(:),bessel_j_1(:),bessel_n_1(:),bessel_j_2(:),bessel_n_2(:),energies(:),cross_secs(:)
    integer::l,i,J,l_max,pos,pos_2,N,a,start,end
    integer,allocatable::l_array(:),local_l_array(:)
    logical::correction
    procedure (two_pointing_func),pointer:: V_ptr => lennard_jones_potential_thijssen
    class(starting_conditions_for_lennard_jones),allocatable ::ustruct

    integer :: ierr 
    integer :: rank ,root
    integer :: nprocs , rem, elements_per_proc, index, n_elements_recieved

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    root = 0
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
    N = 10000
    allocate(energies(N),cross_secs(N))
    energies = my_linspace(0.1_dp,3.5_dp,N)
    epsilon_prime = epsilon/factor


    l_max  = 6
    x_to_check = 5.0_dp
    pi = 3.14159265359_dp
    
    if (nprocs > l_max + 1) then 
        print*,"insert a panic here"
    end if


    sum = 0.0_dp
    allocate(l_array(l_max+1))
    l_array = create_l_array(l_max)

    do i = 1,N
        sum = 0.0
        energy = energies(i)
        energy_prime = energy/factor
        !print*, "energy in meV, energy in nat",energy,energy_prime
        
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
        !print*,"x to check", x_to_check
        !print*,"second x",x_to_check+new_x_diff
        
        V_array = create_potential_array_for_lennard(v_ptr,x_array,epsilon_prime,r_max)
        local_sum = 0.0_dp 
        start = rank * ((l_max+1)/real(nprocs,kind(1d0)))+1
        end= (rank+1)*((l_max+1)/real(nprocs,kind(1d0)))+1
        
        !print*, "start", start,"end",end


        do a=start,end-1
            local_sum = local_sum + calculate_cross_section_l_term(x_array,V_array,l_array(a),energy_prime,epsilon_prime,pos,pos_2,bessel_j_1,bessel_j_2,bessel_n_1,bessel_n_2)
        end do
        !print*, "rank, ", rank, " partial_sum",local_sum
        call mpi_reduce(local_sum,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if (rank == root) then 
            sum = sum * 4.0_dp * pi / wavenumber**2
            !print*, "rank ", rank, "sum", sum
            cross_secs(i) = sum
        end if
    end do

    if (rank == root) then 
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
    end if 
    
    call MPI_FINALIZE(ierr)

end program int_check

