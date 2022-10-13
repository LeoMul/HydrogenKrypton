gfortran -g -fbacktrace -Wall -fcheck=all numerical_methods.f90 physics_functions.f90 h_kr_functions.f90 integration_check.f90 -ffree-line-length-none -o int_check
./int_check
