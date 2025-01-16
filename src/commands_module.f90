module commands_module
  use kind_module
  use utils_module
  implicit none

  private

  public :: commands
  public :: defaults

  character(1),public :: pbc                      ! use periodical boundaries conditions
  real(wp),public     :: rcut                     ! Cut-off radius
  real(wp),public     :: angle                    ! three-body parameter 
  real(wp),public     :: distance                 ! two-body parameter
  character(5),public :: pair(2)                  ! symbols for distance coputation
  character(5),public :: triad(3)                 ! symbols for angle computation
  character(len=:),allocatable,public :: groFile  ! trajectory file in GRO format
  character(len=:),allocatable,public :: graph    ! user-defined prefix of figure
  character(len=:),allocatable,public :: measure  ! user-define measure passed to networks
  !> Test mode used for fixed regular lattice models. Additional file(s) are required
  character(len=:),allocatable,public :: txtFile  ! additional data
  !> Frequency to print graph
  character(1), public :: active                  ! activate print the graph from networkx
  integer, public :: frequency       

  integer,public :: mode 
  ! --pair    -> mode = 1
  ! --triad   -> mode = 2
  ! --lattice -> mode = 3

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
      write(*,"('COOR2GRAPH - Software for build an adjacency matrix from a coordinate file' &
                    &' and compute graph measure.', /&
                    &' Try: coor2graph --help',/&
                    & ' ',/&
                    & 'Report bugs to: aslozada@gmail.com')")
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
           groFile = trim(adjustl(buffer))
         case('--rcut')
           call get_command_argument(iarg+1,buffer)
           Foo = trim(adjustl(buffer))
           read(Foo,*) rcut
         case('--pair')
           mode = 1
           call get_command_argument(iarg+1,buffer)
           pair(1) = trim(adjustl(buffer))
           call get_command_argument(iarg+2,buffer)
           pair(2) = trim(adjustl(buffer))
           call get_command_argument(iarg+3,buffer)
           Foo = trim(adjustl(buffer))
           read(Foo,*) distance
         case('--triad')
           mode = 2
           call get_command_argument(iarg+1,buffer)
           triad(1) = trim(adjustl(buffer))
           call get_command_argument(iarg+2,buffer)
           triad(2) = trim(adjustl(buffer))
           call get_command_argument(iarg+3,buffer)
           triad(3) = trim(adjustl(buffer))
           call get_command_argument(iarg+4,buffer)
           Foo = trim(adjustl(buffer))
           read(Foo,*) angle
         case('--graph')
           call get_command_argument(iarg+1,buffer)
           graph = trim(adjustl(buffer))
         case('--measure')
           call get_command_argument(iarg+1,buffer)
           measure = trim(adjustl(buffer))
         case('--pbc') 
           call get_command_argument(iarg+1,buffer)
           pbc = trim(adjustl(buffer))
         case('--help')
           call help()
         case('--version')
           call version() 
         case('--lattice') 
           mode = 3
           call get_command_argument(iarg+1,buffer)
           txtFile = trim(adjustl(buffer))
           call get_command_argument(iarg+2,buffer)
           Foo = trim(adjustl(buffer))
           read(Foo,*) distance
         case('--frequency')
           call get_command_argument(iarg+1,buffer)
           active = trim(adjustl(buffer))
           call get_command_argument(iarg+2,buffer)
           Foo = trim(adjustl(buffer))
           read(Foo,*) frequency

       end select
    end do
  end subroutine commands 

  subroutine help()
    write(*,hash)
    write(*,'(a)') 'coor2graph [OPTIONS]'
    write(*,'(a)') 'Example: coor2graph --input <coor>.gro --rcut <#> --pair&
       & <sym1> <sym2> <distance> --graph <prefix> --measure <networkx measure>'
    write(*,'(a)') '         coor2graph --input <coor>.gro --rcut <#> --lattice&
       &<Propertie File> <distance> --graph <prefix> --measure <networkx measure>'
    write(*,'(a)') ''
    write(*,'(a)') 'Options:'
    write(*,'(a)') '  --help       | print the help'
    write(*,'(a)') '  --version    | print the version'
    write(*,'(a)') '  --input      | uses a GRO file [STR]'
    write(*,'(a)') '  --rcut       | cut-off as float'
    write(*,'(a)') '  --pair       | [sym1] [sym2] <distance>'
    write(*,'(a)') '  --triad      | [sym1] [sym2] [sym3] <angle>'
    write(*,'(a)') '  --graph      | prefix [STR]'
    write(*,'(a)') '  --measure    | <networkx measure>'
    write(*,'(a)') '  --pbc        | activate the periodical boudary conditions [y|n]'
    write(*,'(a)') '  --lattice    | this mode requires additional files <plain text (in progress)' 
    write(*,'(a)') '               | require a distance value'
    write(*,'(a)') '  --frequency  | [y/n] <-1> or <#> activate and define the number of graph figure ' 
    write(*,'(a)') '               | if mod(nframes,frequency)==0 or frequency=-1 print graph, otherwise no'
    write(*,hash)

    stop
  end subroutine help 

  subroutine version()
    write(*,hash)
    write(*,'(a)') 'COOR2GRAPH version 1.0.0'
    write(*,'(a)') 'Copyright 2024 Asdrubal Lozada'
    write(*,'(a)')''
    write(*,'(a)')' License GPLv3+: GNU GPL version 3 or later &
                    &<http://gnu.org/license/gpl.html>'
    write(*,'(a)') 'Written by Asdrubal Lozada'        
    write(*,'(a)') 'Report bugs to: aslozada@gmail.com'
    write(*,hash)
  end subroutine version

  subroutine defaults()
    groFile   = ''
    rcut      = 0.0_wp
    distance  = 0.0_wp
    angle     = 0.0_wp
    pair(:)   = ''
    triad(:)  = ''
    pbc       = 'n'
    txtFile   = ''
    frequency = -1
    active    = 'y'
  end subroutine defaults
end module commands_module
