! Testing if it is possible to send only part of a vector.
! It works !

program test
   use mpi
   implicit none

   integer, parameter                  :: NUM_ELEM = 10
   integer, dimension(NUM_ELEM)        :: data(NUM_ELEM), receiv(NUM_ELEM)
   integer                             :: ierr
   integer                             :: world_size, world_rank
   integer, dimension(MPI_STATUS_SIZE) :: status

   data = [1,2,3,4,5,6,7,8,9,10]

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, world_rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, world_size, ierr)

   if (world_rank == 0) then
      call MPI_Send(data, 10, MPI_INTEGER, 1, 0, MPI_COMM_WORLD, ierr)
   end if

   if (world_rank == 1) then
      call MPI_Recv(receiv, 11, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, status, ierr)

      print *, receiv
      print *, status
   end if

   call MPI_Finalize(ierr)
end program test