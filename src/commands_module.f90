module commands_module
  use kind_module
  use utils_module
  use networkx_module
  implicit none

  private

  public :: commands
  public :: defaults

  character(len=:), allocatable, public :: fileInput
  character(1),                  public :: binary  ! Binary mixture
  character(1),                  public :: pure    ! Pure substance
  character(1),                  public :: with_pbc
  character(1),                  public :: ising
  real(wp), public     :: rcut    ! Cut-off radius
  real(wp), public     :: cangle  ! three-body parameter 
  real(wp), public     :: cdistance ! two-body parameter
  character(5), public :: pair(2)   ! symbols for distance coputation
  character(5), public :: triad(3)  ! symbols for angle computation
  character(len=:), allocatable, public :: graph
  character(len=:), allocatable, public :: measure
  character(len=:), allocatable, public :: user

  ! bitmask
  integer, public, parameter :: distance = 1   ! 000001 
  integer, public, parameter :: angle    = 2   ! 000010
  integer, public, parameter :: center   = 4   ! 000100

  ! mode
  integer, public :: mode

contains

  subroutine commands()
    integer :: iarg
    character(80) ::  buffer
    character(80), allocatable :: args(:)
    integer :: nargs
    character(len=:), allocatable :: opt
    character(len=:), allocatable :: Foo
    
    !
    call defaults()

    nargs = command_argument_count()
    if(nargs < 1) then
      write(*,hash)
      write(*,'(a)') ' coor2graph builds adjacency matrix from coordinates file '
      write(*,'(a)') ' try: coor2adj --help'
      write(*,'(a)') ' Report bugs to: aslozada@gmail.com'
      write(*,hash)
      stop
    end if

    allocate(args(nargs))

    do iarg = 1, nargs
       call get_command_argument(iarg, buffer)
       args(iarg) = buffer
    end do

    ! parse command

    do iarg = 1, nargs
       opt = trim(adjustl(args(iarg)))
       select case(opt)
         case('--input')
           call get_command_argument(iarg+1,buffer)
           fileInput = trim(adjustl(buffer))
         case('--pure')
           call get_command_argument(iarg+1,buffer)
           pure = trim(adjustl(buffer))
           mode = 1
         case('--binary')
           call get_command_argument(iarg+1,buffer)
           binary = trim(adjustl(buffer))
           mode = 2
         case('--rcut')
           call get_command_argument(iarg+1,buffer)
           Foo = trim(adjustl(buffer))
           read(Foo,*) rcut
         case('--cdistance')
           call get_command_argument(iarg+1,buffer)
           Foo = trim(adjustl(buffer))
           read(Foo,*) cdistance
           mode = 6
         case('--cangle')
           call get_command_argument(iarg+1,buffer)
           Foo = trim(adjustl(buffer))
           read(Foo,*) cangle
           mode = 3
         case('--pair')
           call get_command_argument(iarg+1,buffer)
           pair(1) = trim(adjustl(buffer))
           call get_command_argument(iarg+2,buffer)
           pair(2) = trim(adjustl(buffer))
           mode =4
         case('--triad')
           call get_command_argument(iarg+1,buffer)
           triad(1) = trim(adjustl(buffer))
           call get_command_argument(iarg+2,buffer)
           triad(2) = trim(adjustl(buffer))
           call get_command_argument(iarg+3,buffer)
           triad(3) = trim(adjustl(buffer))
           mode = 5
         case('--help')
           call help()
         case('--graph')
           call get_command_argument(iarg+1,buffer)
           graph = trim(adjustl(buffer))
         case('--measure')
           call get_command_argument(iarg+1,buffer)
           measure = trim(adjustl(buffer))
         case('--pbc') 
           call get_command_argument(iarg+1,buffer)
           with_pbc = trim(adjustl(buffer))
         case('--user') 
           call get_command_argument(iarg+1,buffer)
           user = trim(adjustl(buffer))
         case('--ising')
           call get_command_argument(iarg+1,buffer)
           ising = trim(adjustl(buffer))
           mode = 0
         case default
!           write(*,*) 'try: coor2adj --help'
            mode = 0
            write(*,*) 'Default mode: using ising-like model '
       end select
    end do
  end subroutine commands 

  subroutine help()
    write(*,hash)
    write(*,'(a)') 'coor2adj [OPTIONS]'
    write(*,'(a)') 'Example: coor2adj --input <coor>.gro --pure y --rcut 10.0 --pair O O --graph output.png --measure closenness'
    write(*,'(a)') ''
    write(*,'(a)') 'Options:'
    write(*,'(a)') '  --help'
    write(*,'(a)') '  --input'
    write(*,'(a)') '  --pure'
    write(*,'(a)') '  --binary'
    write(*,'(a)') '  --rcut'
    write(*,'(a)') '  --cangle'
    write(*,'(a)') '  --pair'
    write(*,'(a)') '  --triad'
    write(*,'(a)') '  --graph'
    write(*,'(a)') '  --measure'
    write(*,'(a)') '  --user'
    write(*,'(a)') '  --ising'
    write(*,hash)

    stop
  end subroutine help 

  subroutine defaults()
    binary   = 'n'
    pure     = 'y'
    rcut     = 10.0_wp
    cangle   = 0.0_wp
    pair(:)  = ''
    triad(:) = ''
    with_pbc = 'y'
  end subroutine defaults
end module commands_module
