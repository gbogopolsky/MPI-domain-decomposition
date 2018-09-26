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
      this%to_walk_x = 0
      this%to_walk_y = 0
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

end module

program domain_decomposition
   use mpi
   use walker_mod
   implicit none

   ! Global parameters
   integer, parameter   :: domain_x = 500
   integer, parameter   :: domain_y = 800
   integer, parameter   :: total_num_walkers = 10000
   integer, parameter   :: nCycles = 100000                    ! Number of cycles
   integer, parameter   :: MAX_STEP = 50
   integer, parameter   :: MAX_EXCHANGE = 1000                 ! Size of exchange vector
   type(Walker)         :: walkers(num_walkers)

   ! Runtime variables for decomposition
   integer              :: world_size, world_rank
   integer              :: subdomain_x
   integer              :: subdomain_y
   integer              :: num_walkers                         ! Number of walkers in each domain
   integer              :: num_exchanged
   integer              :: num_deleted
   type(Walker)         :: walker_exchange(MAX_EXCHANGE)
   type(Walker)         :: receive_walkers(MAX_EXCHANGE)
   integer              :: size_exchange

   ! Common variables
   integer              :: ierr, i, icycle
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
   call RANDOM_SEED()
   do i = 1, num_walkers
      call RANDOM_NUMBER(rands)
      call walkers(i)%initialize(i, NINT(rands(1) * subdomain_x), NINT(rands(2) * subdomain_y))
   end do

   ! Start cycle
   do icycle = 1, nCycles

      ! Assign random step and walk
      do i = 1, num_walkers
         call RANDOM_NUMBER(rands)
         rands = rands - 0.5
         call walkers(i)%set_step(NINT(rands(1) * MAX_STEP), NINT(rands(2) * MAX_STEP))
         call walkers(i)%walk
      end do

      ! Test if a walker is outside the domain on the periodic bounds.
      do i = 1, num_walkers
         if (walkers(i)%y > subdomain_y) walkers(i)%y = walkers(i)%y - subdomain_y
         if (walkers(i)%y < 1) walkers(i)%y = walkers(i)%y + subdomain_y
         ! Add walkers to the exchange arrays.
         if (world_rank == 0) then
            num_exchanged = 0
            if (walkers(i)%x < 1 .OR. walkers(i)%x > subdomain_x) then
               walker_exchange(num_exchanged) = walkers(i)
               num_exchanged = num_exchanged + 1
               walkers(i)%id = 0                               ! id = 0 are set for deletion
            end if
         end if
         if (world_rank == 1) then
            num_exchanged = 0
            if (walkers(i)%x < 1 .OR. walkers(i)%x > subdomain_x) then
               walker_exchange(num_exchanged) = walkers(i)
               num_exchanged = num_exchanged + 1
               walkers(i)%id = 0
            end if
         end if
      end do

      ! Delete outgoing walkers and redorder walkers array
      do i = 1, num_walkers
         if (walkers(i)%id == 0) then
            walkers(i) = walkers(i+1)
         end if
      end do
      num_walkers = num_walkers - num_exchanged

      ! Exchange walkers
      size_exchange = SIZEOF(walker_exchange)
      if (world_rank == 0) then
         call MPI_Send(walker_exchange, size_exchange, MPI_BYTE, 1, 0, MPI_COMM_WORLD, ierr)

         call MPI_Recv(receive_walkers, size_exchange, MPI_BYTE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
      end if
      if (world_rank == 1) then
         call MPI_Recv(receive_walkers, size_exchange, MPI_BYTE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

         call MPI_Send(walker_exchange, size_exchange, MPI_BYTE, 0, 0, MPI_COMM_WORLD, ierr)
      end if

      ! Append received walkers
      do i = 1, 

   end do
end program