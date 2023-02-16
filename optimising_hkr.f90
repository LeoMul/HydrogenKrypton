module quick_functions
    use physics_functions
    use numerical_methods
    use h_kr_functions

    implicit none 
        real * 8 , parameter :: mh = 1.6735575e-27_dp 
        real * 8, parameter :: mkr = 1.3914984e-25_dp
        real *8 ,parameter :: reduced_mass = mh*mkr/(mkr+mh)
        real *8 ,parameter :: hbar = 1.05457182e-34_dp
        real * 8 ,parameter ::  pi = 3.14159265359_dp
        real * 8 ,parameter :: mevjoulesfactor = 6.241506363094e+21_dp
        real * 8 ,parameter :: factor_times_rho_squared = mevjoulesfactor* hbar**2 / (reduced_mass * (1e-10_dp)**2)
    contains 

    function calculate_cross_sec(energies,x_min,x_max,epsilon,rho,l_max,integration_step) result(sigmaarray)
        procedure (two_pointing_func),pointer:: V_ptr => lennard_jones_potential_thijssen
        !takes energies, epsilon in meV angstrom
        real* 8,allocatable :: ustartarray(:,:)
        real * 8 ,intent(in) :: energies(:),epsilon,rho,x_min,x_max,integration_step
        integer, intent(in) :: l_max 
        
        real * 8 ,allocatable :: x_array(:),v_array(:)
        real * 8::phases_l(l_max+1) ,delta , term

        real * 8 :: energies_natural_units(size(energies)) , epsilon_natural_units


        real * 8 ::sigmaarray(size(energies)),sum,energy_prime,wavenumber,wavelength,new_x_diff,factor

        integer :: l ,i,pos
        allocate(ustartarray(2,l_max+1))

        factor = factor_times_rho_squared / (rho*rho)
        epsilon_natural_units = epsilon/factor
        energies_natural_units = energies/factor

        do i = 1,size(energies)
            sum = 0.0_dp
            energy_prime = energies_natural_units(i)
            wavenumber = sqrt(2.0_dp*energy_prime) !in rho ^{-1}
            wavelength = 2.0_dp * pi / wavenumber
            new_x_diff = wavelength / 2.0_dp
    
            x_array = my_arange(x_min,x_max+0.75_dp*wavelength,integration_step)
            pos = find_closest_index(x_array,x_max)
            V_array = create_potential_array_for_lennard(v_ptr,x_array,epsilon_natural_units,x_max)
            
            ustartarray = numerov_start_lennard(epsilon_natural_units,x_min,energy_prime,integration_step,l_max)
            phases_l = find_phase_shifts(ustartarray,energy_prime,l_max,x_array,V_array,pos,integration_step)
            do l = 0,l_max
                delta = phases_l(l+1)
                term = (2.0_dp * l + 1.0_dp) * (sin(delta)**2)
                sum = sum + term
            end do 
    
            sum = sum * 4.0_dp * pi / wavenumber**2
            sigmaarray(i) = sum
        end do

    end function calculate_cross_sec

    function calculate_chi_squared(energies,x_min,x_max,epsilon,rho,csa_scaling_factor,l_max,integration_step,experimental_data) result(chisquared)
        real * 8 ,intent(in) :: energies(:),epsilon,rho,x_min,x_max,integration_step,csa_scaling_factor,experimental_data(:)

        integer, intent(in) :: l_max 
        real * 8:: chisquared 
        real * 8 ::theoretical_csa(size(experimental_data)) ,residuals(size(experimental_data))
        real * 8 :: abusive_scale
        theoretical_csa = calculate_cross_sec(energies,x_min,x_max,epsilon,rho,l_max,integration_step)

        !abusive_scale = MAXVAL(theoretical_csa) / MAXVAL(experimental_data)

        !print*,abusive_scale

        residuals = theoretical_csa - 0.45*experimental_data

        chisquared = dot_product(residuals,residuals)
        !print*,"chi squared",chisquared

    end function
    
    function gradient_of_chi_squared(energies,x_min,x_max,epsilon,rho,csa_scaling_factor,l_max,integration_step,experimental_data,derivative_step) result(grad_chisq)

        real * 8 ,intent(in) :: energies(:),epsilon,rho,x_min,x_max,integration_step,csa_scaling_factor,experimental_data(:),derivative_step

        integer, intent(in) :: l_max 
        real * 8:: chisquared,chisquared_rho_step,chisquared_eps_step,chisquared_factor_step

        real * 8 :: grad_chisq(3)

        chisquared = calculate_chi_squared(energies,x_min,x_max,epsilon,rho,csa_scaling_factor,l_max,integration_step,experimental_data)

        chisquared_rho_step = calculate_chi_squared(energies,x_min,x_max,epsilon,rho+derivative_step,csa_scaling_factor,l_max,integration_step,experimental_data)
        chisquared_eps_step = calculate_chi_squared(energies,x_min,x_max,epsilon+derivative_step,rho,csa_scaling_factor,l_max,integration_step,experimental_data)
        chisquared_factor_step = calculate_chi_squared(energies,x_min,x_max,epsilon,rho,csa_scaling_factor+derivative_step,l_max,integration_step,experimental_data)

        grad_chisq(1) = chisquared_rho_step - chisquared
        grad_chisq(2) = chisquared_eps_step - chisquared
        grad_chisq(3) = chisquared_factor_step - chisquared

        grad_chisq = grad_chisq / derivative_step

        end function 




end module quick_functions


program int_check
    use physics_functions
    use numerical_methods
    use h_kr_functions
    use quick_functions
    implicit none

    real * 8 :: x_min,r_max,x_last,h,epsilon,rho,scaling_factor 
    real * 8 ,allocatable :: experimental_cross_secs(:),energies(:)
    integer :: N,l_max
    integer :: j 

    real, dimension(:,:), allocatable :: matrix

    integer :: numrows,numcols,I

    integer :: max_iter,iter 
    real * 8 :: oldx(3),x(3),f ,oldgradf(3),gradf(3),gamma,diff_gradf(3),mag_of_diff_grad,derivative_step


    open (unit=99, file='exp.txt', status='old', action='read')
    

    numrows = 42
    numcols = 2

    allocate(matrix(numrows,numcols))

    do I=1,numrows,1
       read(99,*) matrix(I,:)
       !write(*,*) matrix(I,:)
    enddo
    
    !position stuff
    !in units of rho
    x_min = 0.5_dp
    r_max = 5.0_dp
    x_last = 15_dp !this is for the big integration for the correction term. 
    h = 1.0e-3_dp 
    gamma = 0.0000001_dp

    epsilon = 5.90_dp !meV
    rho = 3.57_dp !in agstroms
    scaling_factor = 1.0_dp
    
    x(1) = rho 
    x(2) = epsilon
    x(3) = scaling_factor

    l_max  = 7

    N = numrows
    
    allocate(energies(N),experimental_cross_secs(N))
    

    !energies = my_linspace(Emin,Elast,N)
    energies = matrix(:,1)
    experimental_cross_secs = matrix(:,2)
   
    derivative_step = 1.0e-9_dp

    max_iter = 100
    gradf = gradient_of_chi_squared(energies,x_min,x_last,x(2),x(1),x(3),l_max,h,experimental_cross_secs,derivative_step)
    print*,x

    oldgradf = gradf
    oldx = x
    x = x - gamma * gradf
    print*,x
    gradf = gradient_of_chi_squared(energies,x_min,x_last,x(2),x(1),x(3),l_max,h,experimental_cross_secs,derivative_step)

    do iter = 1,max_iter
        diff_gradf = gradf - oldgradf
        mag_of_diff_grad = NORM2(diff_gradf)
        if (mag_of_diff_grad .eq. 0.0d0) then 
            EXIT
        end if 
        gamma = dot_product((x-oldx),diff_gradf) / (mag_of_diff_grad*mag_of_diff_grad)
        oldgradf = gradf
        oldx = x
        x = x - gamma * gradf
        gradf = gradient_of_chi_squared(energies,x_min,x_last,x(2),x(1),x(3),l_max,h,experimental_cross_secs,derivative_step)
        print*,x
        
    end do

    print*,mag_of_diff_grad

    !do iter = 1,numrows 
    !    print*, energies(iter),experimental_cross_secs(iter)
    !end do

end program int_check

