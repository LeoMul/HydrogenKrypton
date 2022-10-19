program array_sum 
    USE MPI
    implicit none

    integer :: a(100),local_a(100)


    integer :: i ,n ,sum
    integer :: ierr 
    integer :: rank ,root
    integer :: nprocs , rem, elements_per_proc, index, n_elements_recieved,remaining,local_sum

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    n = 100

    a = 1
    local_a = 0

    if (rank == 0) then 
        elements_per_proc = n/nprocs
        if (nprocs > 1) then 
            do i = 0,nprocs-2
                index = i * elements_per_proc 
                !print*,i
                call MPI_SEND(elements_per_proc,1,MPI_INT,i,0,MPI_COMM_WORLD)
                call MPI_SEND(a(index),elements_per_proc,MPI_INT,i,0,MPI_COMM_WORLD)

            end do
            !print*,i
            index = i * elements_per_proc 
            remaining = n - index

            call MPI_SEND(remaining,1,MPI_INT,i,0,MPI_COMM_WORLD)
            call MPI_SEND(a(index),remaining,MPI_INT,i,0,MPI_COMM_WORLD)


        end if 
        sum = 0
        do i =1,elements_per_proc
            sum = sum + a(i)
        end do

        do i = 0,nprocs-1
            CALL MPI_Recv(local_sum,1,MPI_INT,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,ierr)
            sum = sum + local_sum
        end do 
        print*,sum
    else
        CALL MPI_Recv(n_elements_recieved,1,MPI_INT,0,0,MPI_COMM_WORLD,ierr)
        CALL MPI_Recv(local_a,n_elements_recieved,MPI_INT,0,0,MPI_COMM_WORLD,ierr)
        print*,local_a
        local_sum = 0
        do i = 1,n_elements_recieved
            local_sum = local_sum + local_a(i)
        end do
        CALL MPI_SEND(local_sum,1,MPI_INT,0,0,MPI_COMM_WORLD,ierr)


    end if 



end program array_sum