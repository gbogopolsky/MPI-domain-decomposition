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
   integer, parameter   :: domain_x = 20
   integer, parameter   :: domain_y = 50
   integer, parameter   :: total_num_walkers = 3
   integer, parameter   :: nCycles = 100                   ! Number of cycles
   integer, parameter   :: MAX_STEP = 1
   integer, parameter   :: MAX_EXCHANGE = 50                 ! Size of exchange vector
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
   integer              :: num_incoming
   integer              :: status(MPI_STATUS_SIZE)
   integer              :: numbers(2), incomings(2), exchanged(2)

   ! Diagnostics and output
   integer, dimension(domain_x, domain_y) :: walkersOnDomain
   integer, dimension(:,:), allocatable   :: walkersOnSubdomain

   ! Common variables
   integer              :: ierr, i, j, icycle, ip, id
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
   if (world_rank == 0) num_walkers = num_walkers + MOD(total_num_walkers, 2)

   ! Initializing total_num_walkers / 2 walkers in each domain
   if (world_rank == 1) call RANDOM_SEED()
   do i = 1, num_walkers
      call RANDOM_NUMBER(rands)
      id = i + num_walkers * world_rank
      call walkers(i)%initialize(id, INT(rands(1) * subdomain_x) + 1, INT(rands(2) * subdomain_y) + 1)
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
         else if (walkers(i)%y < 1) then
            walkers(i)%y = walkers(i)%y + subdomain_y
         end if
         ! Add walkers to the exchange arrays.
         if (walkers(i)%x < 1) then
            walkers(i)%x = walkers(i)%x + subdomain_x
            num_exchanged = num_exchanged + 1
            walker_exchange(num_exchanged) = walkers(i)
            walkers(i)%id = -1                                  ! id = 0 are set for deletion
         else if (walkers(i)%x > subdomain_x) then
            walkers(i)%x = walkers(i)%x - subdomain_x
            num_exchanged = num_exchanged + 1
            walker_exchange(num_exchanged) = walkers(i)
            walkers(i)%id = -1
         end if
      end do
      ! call MPI_Barrier(MPI_COMM_WORLD, ierr)
      ! call MPI_Gather(num_exchanged, 1, MPI_INTEGER, exchanged, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      ! if (world_rank == 0) then
      !    print *,"Exchanged walkers:", exchanged
      ! end if

      ! Delete outgoing walkers and redorder walkers array
      ! if (world_rank == 1) then
      !    print *, num_walkers
      !    do i = 1, num_walkers
      !       call walkers(i)%print
      !    end do
      !    print *
      !    print *
      ! end if
      ip = 0
      num_walkers = num_walkers - num_exchanged
      do i = 1, num_walkers
         do while (walkers(i + ip)%id == -1)
            ip = ip + 1
         end do
         walkers(i) = walkers(i + ip)
      end do
      ! if (world_rank == 1) then
      !    print *, num_walkers
      !    do i = 1, num_walkers
      !       call walkers(i)%print
      !    end do
      !    print *
      !    print *
      ! end if
      ! if (world_rank == 1) then
      !    print *, num_exchanged
      !    do i = 1, num_exchanged
      !       call walker_exchange(i)%print
      !    end do
      !    print *
      !    print *
      ! end if

      ! Exchange walkers
      size = SIZEOF(walkers(1))
      if (world_rank == 0) then
         call MPI_Send(walker_exchange, num_exchanged * size, MPI_BYTE, 1, 0, MPI_COMM_WORLD, ierr)

         call MPI_Recv(received_walkers, MAX_EXCHANGE * size, MPI_BYTE, 1, 0, MPI_COMM_WORLD, status, ierr)
         num_incoming = status(1) / size
         ! print *, num_incoming
         ! do i = 1, num_incoming
         !    call received_walkers(i)%print
         ! end do
         ! print *
         ! print *
      end if
      if (world_rank == 1) then
         call MPI_Recv(received_walkers, MAX_EXCHANGE * size, MPI_BYTE, 0, 0, MPI_COMM_WORLD, status, ierr)

         call MPI_Send(walker_exchange, num_exchanged * size, MPI_BYTE, 0, 0, MPI_COMM_WORLD, ierr)
         num_incoming = status(1) / size
         ! print *, num_incoming
         ! do i = 1, num_incoming
         !    call received_walkers(i)%print
         ! end do
         ! print *
         ! print *
      end if

      ! Append received walkers
      ! if (world_rank == 0) then
      !    print *, num_walkers
      !    do i = 1, num_walkers
      !       call walkers(i)%print
      !    end do
      !    print *
      !    print *
      ! end if
      do i = 1, num_incoming
         walkers(num_walkers + i) = received_walkers(i)
      end do
      num_walkers = num_walkers + num_incoming
      ! if (world_rank == 0) then
      !    print *, num_walkers
      !    do i = 1, num_walkers
      !       call walkers(i)%print
      !    end do
      !    print *
      !    print *
      ! end if

      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call MPI_Gather(num_walkers, 1, MPI_INTEGER, numbers, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Gather(num_incoming, 1, MPI_INTEGER, incomings, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (world_rank == 0) then
         print *,"Number of walkers per domain:", numbers, "Total = ", SUM(numbers)
         print *,"Incoming in each domain:", incomings
      end if

      allocate(walkersOnSubdomain(subdomain_x, subdomain_y))
      do i = 1, num_walkers
         walkersOnSubdomain(walkers(i)%y, walkers(i)%x) = walkersOnSubdomain(walkers(i)%y, walkers(i)%x) + 1
      end do
      ! if (world_rank == 1) then
      !    call MPI_Send(walkersOnSubdomain, subdomain_x * subdomain_y * 4, MPI_BYTE, 0, 0, MPI_COMM_WORLD, ierr)
      ! end if
      if (world_rank == 0) then
         do j = 1, subdomain_x
            do i = 1, subdomain_y
               walkersOnDomain(i,j) = walkersOnSubdomain(i,j)
            end do
         end do
         print *
         ! call MPI_Recv(walkersOnSubdomain, subdomain_x * subdomain_y * 4, MPI_BYTE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
         ! do j = 1, subdomain_x
         !    do i = 1, subdomain_y
         !       walkersOnDomain(domain_y + i,j) = walkersOnSubdomain(i,j)
         !    end do
         ! end do
         call MPI_Gather(walkersOnSubdomain, subdomain_x * subdomain_y * 4, MPI_BYTE, walkersOnDomain, domain_x * domain_y * 4)
         print *, walkersOnDomain
         if (MOD(icycle,10) == 0) then
            open(200, file='movie/position.res', status='unknown', position='append')
               do i = 1, domain_y
                  write(200, 101) (walkersOnDomain(i,j), j = 1, domain_x)
               end do
            close(200)
         end if 
      end if
101 format (501(2x,I5))
      deallocate(walkersOnSubdomain)
   end do
   call MPI_Finalize(ierr)
end program