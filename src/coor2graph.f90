!coor2graph - A software tool for converting coordinate data into graph representations
!and performing graph analysis. It generates adjacency matrices from coordinates files
!obtained from molecular dynamics simulations, stochastic simulations, lattices models
!and other sources, and uses NetworkX to compute graph measures.
!-------------------------------------------------------------------
!> author: Asdrubal Lozada
!  update: 4/11/2026
!  license: GPL-v3 
!  Copyright 2024-2026
!-------------------------------------------------------------------
! Local app

program coor2graph
  use kind_module
  use utils_module
  use commands_module
  use frame_module

  implicit none

  integer :: unit
  type(frame) :: gro

  call commands()
  call openFile(unit,groFile,'old')
  call process_trajectory(gro,unit)

  if(net=='y') call process_graph()

contains
  ! Optional procedure to calculate the measure graph from networkx  
  subroutine process_graph()
    implicit none
    
    character(250)                 :: foo
    character(80)                  :: buffer
    character(len=:), allocatable  :: measure_type
    character(len=:), allocatable  :: measure_name
    character(len=:), allocatable  :: shape_intervals
    character(len=:), allocatable  :: nframes
    character(len=:), allocatable  :: ncores
    character(len=:), allocatable  :: draw_freq
    character(len=:), allocatable  :: cmd
    character(len=*), parameter    :: python_exe = 'python'
    character(len=*), parameter    :: process_all = 'process_all.py'

    integer :: cmdstat, exitstat
    character(len=256) :: cmdmsg


    ! defaults
    nframes         = '1'
    ncores          = '1'
    draw_freq       = '1'
    measure_type    = 'degree'
    measure_name    = 'degree'
    shape_intervals = 'none'


    cmd = ""

    write(*,'(a)') " Enter the number of frames to process:"
    read(*,'(a)') buffer
    nframes = trim(adjustl(buffer))

    write(*,'(a)') " Enter the number of CPU cores:"
    read(*,'(a)') buffer
    ncores = trim(adjustl(buffer))

    write(*,'(a)') " Enter the graph generation frequency:"
    read(*,'(a)') buffer
    draw_freq = trim(adjustl(buffer))

    write(*,'(a)') " Enter the graph measure type:"
    read(*,'(a)') buffer
    measure_type = trim(adjustl(buffer))

    write(*,'(a)') " Enter the graph measure name:"
    read(*,'(a)') buffer
    measure_name = trim(adjustl(buffer))

    write(*,'(a)') " Enter the interval shape file name:"
    read(*,'(a)') buffer
    shape_intervals = trim(adjustl(buffer))

    cmd = python_exe//" "//process_all//" --nframes "//nframes//" --ncores "//ncores&
        &//" --measure_type "//measure_type//" --measure_name "//measure_name&
        &//" --shape_intervals "//shape_intervals//" --draw_freq "//draw_freq

    call execute_command_line( &
    command=cmd, &
    wait=.true., &
    exitstat=exitstat, &
    cmdstat=cmdstat, &
    cmdmsg=cmdmsg)

    if (cmdstat /= 0) then
       write(*,'(a)') 'Error executing Python command: '//trim(cmdmsg)
       error stop 1
    end if

    if (exitstat /= 0) then
       write(*,'(a,i0)') 'Python process returned error code: ', exitstat
       error stop 1
    end if
   

  end subroutine process_graph 
end program coor2graph

