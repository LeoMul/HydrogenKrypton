module comm_interface
      !USE KIND_PARAMETERS, ONLY : RKIND,IKIND
      implicit none
      module load gcc
      module load openmpi
      include "mpif.h"
      public comm_init          ! Initialize MPI
      public comm_barrier       ! MPI barrier
      public comm_finalize      ! Terminate MPI
      integer, public  :: iam,nproc
C      SAVE
      private
      integer :: mpicom

      CONTAINS

      subroutine comm_init()
c      USE KIND_PARAMETERS, ONLY : RKIND,IKIND
      implicit none
      integer :: ier

      mpicom = MPI_COMM_WORLD
      call mpi_init(ier)
      call mpi_comm_rank(mpicom, iam, ier)
      call mpi_comm_size(mpicom, nproc, ier)

      return
      end subroutine comm_init
      
      subroutine comm_barrier()
c      USE KIND_PARAMETERS, ONLY : RKIND,IKIND
      implicit none
      integer :: ier

      call mpi_barrier(mpicom, ier)

      return
      end subroutine comm_barrier
!-----------------------------------------------------------------------
       subroutine comm_finalize()
c      USE KIND_PARAMETERS, ONLY : RKIND,IKIND
      implicit none
      integer :: ier

      call mpi_finalize(ier)

      return
       end subroutine comm_finalize
!-----------------------------------------------------------------------
      end module comm_interface

      program hello_world
      use comm_interface
      
      implicit real*8(a-h,o-z)

      call comm_init()
      
      write(0,*)'hello from proc',iam

      call comm_barrier()
      call comm_finalize()
      stop
      end
