mpif90 -g -fbacktrace -Wall -fcheck=all numerical_methods.f90 physics_functions.f90 h_kr_functions.f90 scattering_integration_check.f90 -ffree-line-length-none -Wno-maybe-uninitialized -o intcheck 
./intcheck
