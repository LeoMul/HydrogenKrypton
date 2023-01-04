program int_check
    use physics_functions
    use numerical_methods
    use h_kr_functions

    USE MPI

    implicit none

    real*8::rho,epsilon,h,energy,x_last,x_min,epsilon_prime,energy_prime,factor,x_to_check,Emin,Elast,alpha
    real*8::sum,pi,new_x_diff,wavenumber,reduced_mass,hbar,mh,mkr,wavelength,r_max,delta,norm_factor,term
    real*8,allocatable::x_array(:),big_x_array(:),V_array(:),bessel_j_1(:),bessel_n_1(:),bessel_j_2(:),bessel_n_2(:),energies(:),cross_secs(:),phases(:),local_phases(:)

    real* 8,allocatable :: ustartarray(:,:)


    integer::l,i,J,l_max,pos,pos_2,N,a,start,end_ind
    integer,allocatable::l_array(:)
    real * 8,allocatable::phases_l(:)
    logical::correction
    procedure (two_pointing_func),pointer:: V_ptr => lennard_jones_potential_thijssen
    class(starting_conditions_for_lennard_jones),allocatable ::ustruct
    character*250::string_1
    character*268::string
    character*272::file_name_1,file_name_2

    !mpi stuff
    integer :: num_energies_per_proc
    real*8,allocatable:: energies_on_this_processor(:), local_cross_secs(:),local_terms(:),terms(:)
    integer :: ierr 
    integer :: rank ,root
    integer :: nprocs 


    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    root = 0
    allocate(ustruct)   
    !position stuff
    !in units of rho
    
    alpha = 1.0_dp
    r_max = 20.0_dp
    h = 1.0e-3_dp 
    x_min = h

    correction = .false.

    !Parameters and stuff
    l = 0
    Emin = 0.1_dp 
    Elast = 100.0_dp
    l_max  = 10
    x_to_check = 5.0_dp
    pi = 3.14159265359_dp
    sum = 0.0_dp
    allocate(l_array(l_max+1),phases_l(l_max+1))
    l_array = create_l_array(l_max)

    N = 1000
    N = N - modulo(N,nprocs) 
    num_energies_per_proc = N/nprocs
    if (rank == root) then 
        print*,"Allocating ", num_energies_per_proc, "energies per processor."
    end if 
    allocate(energies(N),cross_secs(N),phases((l_max+1)*N),terms((l_max+1)*N))
    allocate(energies_on_this_processor(num_energies_per_proc),local_cross_secs(num_energies_per_proc),local_phases((l_max+1)*num_energies_per_proc),local_terms((l_max+1)*num_energies_per_proc))
    
    allocate(ustartarray(2,l_max+1))


    energies = my_linspace(Emin,Elast,N)

    start = rank * num_energies_per_proc + 1 
    end_ind = (rank+1) * num_energies_per_proc
    energies_on_this_processor = energies(start:end_ind)
    print*,size(energies_on_this_processor)

    do i = 1,num_energies_per_proc
        sum = 0.0_dp
        energy = energies_on_this_processor(i)
        energy_prime = energy
        wavenumber = sqrt(2.0_dp*energy_prime) !in rho ^{-1}
        wavelength = 2.0_dp * pi / wavenumber
        new_x_diff = wavelength / 2.0_dp

        x_array = my_arange(x_min,r_max+0.75_dp*wavelength,h)
        pos = find_closest_index(x_array,r_max)
   
        V_array = create_potential_array_for_screened(x_array,alpha,r_max)
        
        ustartarray = numerov_start_screened(h,l_max)
        phases_l = find_phase_shifts(ustartarray,energy_prime,l_max,x_array,V_array,pos,h)
        do l = 0,l_max
            delta = phases_l(l+1)
            term = (2.0_dp * l + 1.0_dp) * (sin(delta)**2)
            !print*,term
            sum = sum + term
            local_phases((i-1)*(l_max+1)+l+1) = delta
            local_terms((i-1)*(l_max+1)+l+1) = term
        end do 

        sum = sum * 4.0_dp * pi / wavenumber**2
        local_cross_secs(i) = sum
        
    end do
    print*,"Processor ",rank,"complete."

    call MPI_GATHER(local_cross_secs,num_energies_per_proc,MPI_DOUBLE_PRECISION,cross_secs,num_energies_per_proc,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
    call MPI_GATHER(local_phases,num_energies_per_proc*(l_max+1),MPI_DOUBLE_PRECISION,phases,num_energies_per_proc*(l_max+1),MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
    call MPI_GATHER(local_terms,num_energies_per_proc*(l_max+1),MPI_DOUBLE_PRECISION,terms,num_energies_per_proc*(l_max+1),MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)

    if (rank == root) then 
        print*,ustartarray

        print*,""

        WRITE (string_1, '(a4,F20.3,a4,F20.3,a1,i48.6,a1,F20.6,a5,F20.6)') "Emin", Emin,"Emax",Elast,"N",N,"h",h,"alpha",alpha
        string = "cross_sec_screened"//string_1
        file_name_1 = string//".dat"
        file_name_1=trim(adjustl(file_name_1))
        string = "phases_screened"//string_1
        file_name_2 = string//".dat"
        file_name_2=trim(adjustl(file_name_2))
        call StripSpaces(file_name_1)
        call StripSpaces(file_name_2)

        print*,"creating and writing cross sections to ",file_name_1



        open (1, file = file_name_1)
        J = 1 !every jth point
        write(1,fmt="(1x,a)",advance = "no") "    "

        write(1,fmt="(1x,a)",advance = "no") "#EnergyColumn"
        write(1,fmt="(1x,a)",advance = "no") "      "
        write(1,fmt="(1x,a)",advance = "no") "#CrossSecCol"
        write(1,fmt="(1x,a)",advance = "no") "  "

        do J = 1,size(l_array)
            write(1,FMT,advance = "no") real(l_array(J))
            write(1,fmt="(1x,a)",advance = "no") " "
        end do 
        write(1,*) " "
        J = 1
        do i = 1,size(energies)
                    write(1,FMT,advance = "no") energies(i)
                    write(1,fmt="(1x,a)",advance = "no") " "
                    write(1,FMT,advance = "no") cross_secs(i) 

                    norm_factor = 2.0_dp*pi/(energies(i))
                    do J = 1,size(l_array)
                        write(1,FMT,advance = "no") terms((i-1)*(l_max+1)+J)*norm_factor
                        write(1,fmt="(1x,a)",advance = "no") " "

                    end do

                    write(1,*) " "

        end do  



        print*,"creating and writing phases to ",file_name_2

        open (1, file = file_name_2)
        J = 1 !every jth point
        write(1,fmt="(1x,a)",advance = "no") "    "

        write(1,fmt="(1x,a)",advance = "no") "#EnergyColumn"
        write(1,fmt="(1x,a)",advance = "no") "  "
        do J = 1,size(l_array)
            write(1,FMT,advance = "no") real(l_array(J))
            write(1,fmt="(1x,a)",advance = "no") " "
        end do 
        write(1,*) " "


        do i = 1,size(energies)
                    write(1,FMT,advance = "no") energies(i)
                    write(1,fmt="(1x,a)",advance = "no") " "
                    norm_factor = 2.0_dp*pi/(energies(i)/factor)
                    !print*,"new energy"

                    do J = 1,size(l_array)
                        !print*,(i-1)*(l_max+1)+J
                        write(1,FMT,advance = "no") phases((i-1)*(l_max+1)+J)
                        write(1,fmt="(1x,a)",advance = "no") " "

                    end do 
                    write(1,*) " "
        end do

        
    end if 



    
    call MPI_FINALIZE(ierr)

end program int_check

