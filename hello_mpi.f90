module place 
    use mpi

    public 

    contains 
    subroutine par_avg(m,a,b,c,comm) 
        REAL a(m), b(m)       ! local slice of array
        REAL c                ! result (at process zero)
        REAL sum
        INTEGER m, comm, i, ierr
        ! local sum
        sum = 0.0
        DO i = 1, m
           sum = sum + a(i)*b(i)
        END DO
        ! global sum
        CALL MPI_REDUCE(sum, c, 1, MPI_REAL, MPI_SUM, 0, comm, ierr)
        RETURN
    end subroutine par_avg
end module place

program hello_mpi

    use mpi
    use place
    implicit none   
    integer :: ierr 
    integer :: rank ,root
    integer :: nprocs , rem
    integer,allocatable :: sendcounts(:), displacements(:)
    real, allocatable::apart(:),bpart(:)
    integer, parameter :: dim1 = 80, dim2 = 10
    integer m,i ,sum_jobs,local_n
    REAL :: x , sum,a(dim1),b(dim1),c,cpart

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    allocate(sendcounts(nprocs),displacements(nprocs))
    
    root = 0

    if (rank == root ) then 
        b = 2.0
        a = 1.0
        !print*,"qwe",sum(b)
        !print*,a
    end if 

    rem = modulo(dim1,nprocs)
    if (rank < rem) then 
        local_n = dim1/nprocs + 1
    else 
        local_n = dim1/nprocs 

    end if  
    !print*,local_n

    !print*, rem
    sum_jobs = 0
    do i = 1,nprocs
        sendcounts(i) = dim1/nprocs 
        if (rem > 0) then 
            sendcounts(i) = sendcounts(i) + 1
            rem = rem - 1
        end if 
        displacements(i) = sum_jobs
        sum_jobs = sum_jobs + sendcounts(i)
    end do 

    allocate(apart(local_n),bpart(local_n))
    !apart = 0
    !bpart = 0
    !print*, size(apart)
    !print*, sendcounts

    call MPI_SCATTERV( b, sendcounts,displacements, MPI_REAL, bpart, local_n, MPI_REAL, root, MPI_COMM_WORLD, ierr )
    call MPI_SCATTERV( a, sendcounts,displacements, MPI_REAL, apart, local_n, MPI_REAL, root, MPI_COMM_WORLD, ierr )


    !print*,rank
    !call par_avg(dim2,a,b,c,MPI_COMM_WORLD)
    cpart = 0.0
    !print*, bpart
    do i = 1, local_n
        cpart = cpart + apart(i)*bpart(i)
    enddo
    print*,cpart
    !print*, cpart
    call MPI_REDUCE( cpart, c, 1, MPI_REAL, MPI_SUM, root, MPI_COMM_WORLD, ierr )

    if (rank == root) then 
        print*,c
    end if 

   ! print*, "Hello world! I am processor ", rank, "of ", nprocs, "processors"


    call MPI_FINALIZE(ierr)

end program hello_mpi