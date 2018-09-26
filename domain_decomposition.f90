module walker_mod
   use mpi
   implicit none

   type Walker
      integer  :: id
      integer  :: x
      integer  :: y
      integer  :: step_x
      integer  :: step_y
   contains
      procedure   :: initialize
      procedure   :: set_step
      procedure   :: walk
      procedure   :: print
   end type

contains

   subroutine initialize(this, id, x, y)
      implicit none

      class(Walker)  :: this
      integer        :: id
      integer        :: x
      integer        :: y

      this%id = id
      this%x = x
      this%y = y
      this%step_x = 0
      this%step_y = 0
   end subroutine

   subroutine set_step(this, step_x, step_y)
      implicit none

      class(Walker)  :: this
      integer        :: step_x
      integer        :: step_y

      this%step_x = step_x
      this%step_y = step_y

   end subroutine

   subroutine walk(this)
      implicit none

      class(Walker)  :: this

      this%x = this%x + this%step_x
      this%y = this%y + this%step_y

   end subroutine

   subroutine print(this)
      implicit none
      class(Walker)  :: this
      print '("id = ",I5,", x = ",I5,", y = ",I5,", step_x = ",I5,", step_y = ",I5)', &
            this%id, this%x, this%y, this%step_x, this%step_y
   end subroutine

end module

program domain_decomposition
   use mpi
   use walker_mod
   implicit none

   ! Global parameters
   integer, parameter   :: domain_x = 500
   integer, parameter   :: domain_y = 800
   integer, parameter   :: total_num_walkers = 10000
   integer, parameter   :: nCycles = 20                   ! Number of cycles
   integer, parameter   :: MAX_STEP = 50
   integer, parameter   :: MAX_EXCHANGE = 5000                 ! Size of exchange vector
   type(Walker)         :: walkers(total_num_walkers)

   ! Runtime variables for decomposition
   integer              :: world_size, world_rank
   integer              :: subdomain_x
   integer              :: subdomain_y
   integer              :: num_walkers                         ! Number of walkers in each domain
   integer              :: num_exchanged
   integer              :: num_deleted
   type(Walker)         :: walker_exchange(MAX_EXCHANGE)
   type(Walker)         :: received_walkers(MAX_EXCHANGE)
   integer              :: size
   integer              :: size_incoming
   integer              :: status(MPI_STATUS_SIZE)
   integer              :: numbers(2), incomings(2), exchanged(2)

   ! Common variables
   integer              :: ierr, i, icycle, ip
   real, dimension(2)   :: rands

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, world_rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, world_size, ierr)

   ! We only consider two domains for now. Indeed, the IMIP computers
   ! only have two processors, and we want to use OpenMP.
   if (world_size /= 2) then
      print *, "World_size must be 2 for this first case."
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
   end if

   ! Subdomain parameters
   subdomain_x = domain_x / 2
   subdomain_y = domain_y
   num_walkers = total_num_walkers / 2                         ! Initial value

   ! Initializing total_num_walkers / 2 walkers in each domain
   if (world_rank == 1) call RANDOM_SEED()
   do i = 1, num_walkers
      call RANDOM_NUMBER(rands)
      call walkers(i)%initialize(i, NINT(rands(1) * subdomain_x), NINT(rands(2) * subdomain_y))
   end do

   ! Start cycle
   do icycle = 1, nCycles

      ! Assign random step and walk
      do i = 1, num_walkers
         call RANDOM_NUMBER(rands)
         rands = (rands - 0.5) * 2                             ! -MAX_STEP <= rands <= MAX_STEP
         call walkers(i)%set_step(NINT(rands(1) * MAX_STEP), NINT(rands(2) * MAX_STEP))
         call walkers(i)%walk
      end do

      ! Test if a walker is outside the domain on the periodic bounds.
      num_exchanged = 0
      do i = 1, num_walkers
         if (walkers(i)%y > subdomain_y) then
            walkers(i)%y = walkers(i)%y - subdomain_y
         else if (walkers(i)%y < 0) then
            walkers(i)%y = walkers(i)%y + subdomain_y
         end if
         ! Add walkers to the exchange arrays.
         if (walkers(i)%x < 0) then
            walkers(i)%x = walkers(i)%x + subdomain_x
            walker_exchange(num_exchanged) = walkers(i)
            num_exchanged = num_exchanged + 1
            walkers(i)%id = 0                                  ! id = 0 are set for deletion
         else if (walkers(i)%x > subdomain_x) then
            walkers(i)%x = walkers(i)%x - subdomain_x
            walker_exchange(num_exchanged) = walkers(i)
            num_exchanged = num_exchanged + 1
            walkers(i)%id = 0
         end if
      end do
      ! call MPI_Barrier(MPI_COMM_WORLD, ierr)
      ! call MPI_Gather(num_exchanged, 1, MPI_INTEGER, exchanged, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      ! if (world_rank == 0) then
      !    print *,"Exchanged walkers:", exchanged
      ! end if

      ! Delete outgoing walkers and redorder walkers array
      ip = 0
      do i = 1, num_walkers
         if (walkers(i)%id /= 0) then
            ip = ip + 1
            walkers(ip) = walkers(i)
         end if
      end do
      num_walkers = num_walkers - num_exchanged

      ! Exchange walkers
      size = SIZEOF(walkers(1))
      print *, "size", size
      if (world_rank == 0) then
         call MPI_Send(walker_exchange, num_exchanged * size, MPI_BYTE, 1, 0, MPI_COMM_WORLD, ierr)

         call MPI_Recv(received_walkers, MAX_EXCHANGE * size, MPI_BYTE, 1, 0, MPI_COMM_WORLD, status, ierr)
         size_incoming = status(1) / size
      end if
      if (world_rank == 1) then
         call MPI_Recv(received_walkers, MAX_EXCHANGE * size, MPI_BYTE, 0, 0, MPI_COMM_WORLD, status, ierr)

         call MPI_Send(walker_exchange, num_exchanged * size, MPI_BYTE, 0, 0, MPI_COMM_WORLD, ierr)
         size_incoming = status(1) / size
      end if

      ! Append received walkers
      do i = 1, size_incoming
         walkers(num_walkers + size_incoming) = received_walkers(i)
      end do

      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call MPI_Gather(num_walkers, 1, MPI_INTEGER, numbers, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
      call MPI_Gather(size_incoming, 1, MPI_INTEGER, incomings, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
      if (world_rank == 0) then
         print *,"Number of walkers per domain:", numbers, "Total = ", SUM(numbers)
         print *,"Incoming in each domain:", incomings
      end if
   end do
   call MPI_Finalize(ierr)
end program