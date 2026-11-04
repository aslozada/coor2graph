!-------------------------------------------------------------------
! This module is part of coor2graph
! Parse command-line arguments and set execution parameters
!-------------------------------------------------------------------
!> author: Asdrubal Lozada
!  update: 4/11/2026
!  license: GPL-v3 

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
  character(5),public :: net                      ! activate graph mesure calculation
  !> Test mode used for fixed regular lattice models. Additional file(s) are required
  character(len=:),allocatable,public :: txtFile  ! additional data

  integer,public :: mode 
  ! --pair    -> mode = 1
  ! --triad   -> mode = 2
  ! --lattice -> mode = 3

  !> versioning
  integer, parameter :: major = 1
  integer, parameter :: minor = 0
  integer, parameter :: patch = 1
 

  contains

  subroutine commands()
    integer :: iarg
    character(80) ::  buffer
    character(80), allocatable :: args(:)
    integer :: nargs
    character(len=:), allocatable :: opt
    character(len=:), allocatable :: Foo
    integer :: ios
    
    call defaults()

    nargs = command_argument_count()
    if(nargs < 1) then
     write(*,hash)

     write(*,"('coor2graph - A software tool for converting coordinate data into graph representations',/&
        & 'and performing graph analysis. It generates adjacency matrices from coordinates files',/&
        & 'obtained from molecular dynamics simulations, stochastic simulations, lattices models',/&
        & 'and other sources, and uses NetworkX to compute graph measures.', /&
        & '',/&
        & 'try: coor2graph --help',/&
        & ' ',/&
        & 'Report bugs to aslozada@gmail,.com')")
     write(*,hash)
     stop
    end if

    allocate(args(nargs))

    do iarg = 1, nargs
       call get_command_argument(iarg, buffer)
       args(iarg) = buffer
    end do

    iarg = 1
    do while (iarg <= nargs)
       opt = trim(adjustl(args(iarg)))
       select case(opt)
         case('--input')
           call get_command_argument(iarg+1,buffer)
           if(iarg+1 > nargs) call fail('missing value for --input')
           groFile = trim(adjustl(buffer))
           iarg = iarg + 2
         case('--rcut')
           call get_command_argument(iarg+1,buffer)
           Foo = trim(adjustl(buffer))
           if(iarg+1 > nargs) call fail('missing value for --rcut')
           read(Foo,*,iostat=ios) rcut
           if (ios /=0) call fail('invalid value for --rcut')
           iarg = iarg + 2
         case('--pair')
           if (iarg+3 > nargs) call fail('missing values for --pair')
           mode = 1
           call get_command_argument(iarg+1,buffer)
           pair(1) = trim(adjustl(buffer))
           call get_command_argument(iarg+2,buffer)
           pair(2) = trim(adjustl(buffer))
           call get_command_argument(iarg+3,buffer)
           Foo = trim(adjustl(buffer))
           read(Foo,*,iostat=ios) distance
           if (ios /=0) call fail('invalid value for distance')
           iarg = iarg + 4
         case('--triad')
           if (iarg + 4 > nargs) call fail('missing values for --triad')
           mode = 2
           call get_command_argument(iarg+1,buffer)
           triad(1) = trim(adjustl(buffer))
           call get_command_argument(iarg+2,buffer)
           triad(2) = trim(adjustl(buffer))
           call get_command_argument(iarg+3,buffer)
           triad(3) = trim(adjustl(buffer))
           call get_command_argument(iarg+4,buffer)
           Foo = trim(adjustl(buffer))
           read(Foo,*,iostat=ios) angle
           if (ios /=0) call fail('invalid value for angle')
           iarg = iarg + 5
         case('--pbc') 
           call get_command_argument(iarg+1,buffer)
           if (iarg + 1 > nargs) call fail('missing value for --pbc')
           pbc = trim(adjustl(buffer))
           iarg = iarg + 2
         case('--help')
           call help()
           stop
         case('--version')
           call version() 
           stop
         case('--lattice') 
           if (iarg+1 > nargs) call fail('missing arguments for --lattice')
           mode = 3
           call get_command_argument(iarg+1,buffer)
           txtFile = trim(adjustl(buffer))
           iarg = iarg + 2
         case('--net')
           call get_command_argument(iarg+1,buffer)
           if (iarg + 1 > nargs) call fail('missing value for --net') 
           net = trim(adjustl(buffer))
           iarg = iarg + 2 
         case default 
           call fail('unknown option: '//opt)     
       end select
    end do
  end subroutine commands 

  subroutine help()
    write(*,hash)
    write(*,'(a)') 'coor2graph [OPTIONS]'
    write(*,'(a)') 'Example: coor2graph --input <coor>.gro --rcut <#> --pair&
       & <sym1> <sym2> <distance>'
    write(*,'(a)') ''
    write(*,'(a)') 'Options:'
    write(*,'(a)') '  --help       | display this help message'
    write(*,'(a)') '  --version    | display the version'
    write(*,'(a)') '  --input      | uses a coordinates file [STR]'
    write(*,'(a)') '  --rcut       | cutoff distance [FLOAT]'
    write(*,'(a)') '  --pair       | [sym1] [sym2] <distance>'
    write(*,'(a)') '  --triad      | [sym1] [sym2] [sym3] <angle>'
    write(*,'(a)') '  --pbc        | enable or disable periodic boundary [y|n]'
    write(*,'(a)') '  --lattice    | enable lattice mode (requires additional files)' 
    write(*,'(a)') '  --net        | enable graph measure calculation (requires additonal files)' 

    write(*,'(a)') ''
    write(*,'(a)') 'options --net:'
    write(*,'(a)') '               | ncores - nframes - measure_type - measure_name'

    write(*,hash)

    stop
  end subroutine help 

  subroutine version()
    write(*,hash)

    write(*,'("coor2graph version ",i0,".",i0,".",i0)') &
    major, minor, patch
    write(*,'(a)') 'Copyright 2024-2026 Asdrubal Lozada'
    write(*,'(a)')''
    write(*,'(a)')'License GPLv3+: GNU GPL version 3 or later &
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
    mode      = 0
  end subroutine defaults
end module commands_module
