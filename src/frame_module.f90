module frame_module
  use kind_module
  use utils_module
  use commands_module
  use networkx_module
  implicit none

  private

  type, public :: atom
    integer      :: index
    character(5) :: group
    character(5) :: sym
    integer      :: seq
    real(wp)     :: xyz(3)
  end type atom

  type, public :: frame
    type(atom), allocatable :: atoms(:)
    integer  :: natoms
    real(wp) :: box(3)
  contains
    procedure :: read => read_frame
  end type frame

  integer, public :: nframes
  public :: on_the_fly
  !-----------------------------
  integer, public :: group_count
  public :: get_groups
  !-------------------------------------------------
  character(5), allocatable :: molName(:), sym(:,:)
  real(wp), allocatable     :: rx(:), ry(:), rz(:)
  real(wp), allocatable     :: x(:,:), y(:,:), z(:,:)
  integer, allocatable      :: na(:)

  !@!@logical, allocatable      :: sites(:)
  !-------------------------------------------------
  real(wp) :: rcut2
  integer, allocatable :: adj(:,:)
  logical :: check

  ! dmax value
  real(wp) :: dmax

contains

  subroutine read_frame(me,unit)
    class(frame), intent(inout) :: me
    integer, intent(in) :: unit 
    ! local variables
    integer       :: n
    integer       :: natoms

    read(unit,*) 
    read(unit,*) natoms

    allocate(me%atoms(natoms))

    me%natoms = natoms
    
    do n = 1, me%natoms
      read(unit,'(i5,2a5,i5,3f8.3)') &
      &me%atoms(n)%index,me%atoms(n)%group,me%atoms(n)%sym,me%atoms(n)%seq,me%atoms(n)%xyz(:)
    end do

    read(unit,*) me%box(:)
  end subroutine read_frame

  subroutine get_groups(me,unit)
    type(frame), intent(inout) :: me
    integer, intent(in) :: unit
    ! local variables
    integer :: n, k, iter
!    integer, allocatable :: na(:) ! atoms by group
    real(wp) :: total_mass
    real(wp) :: sx, sy, sz
    real(wp) :: com(3)
    character(5) :: Foo
    integer :: nmax
    real(wp) :: box(3)

    integer :: lines
    real(wp) :: minimal_image


    lines = 0

    group_count = 0

    do n = 1, me%natoms
       if(me%atoms(n)%index > group_count) then
         group_count = me%atoms(n)%index
       end if
    end do

    allocate(na(group_count))

    ! compute the center of mass for each group
    do k = 1, group_count
       sx = 0.0_wp; sy = 0.0_wp; sz = 0.0_wp
       total_mass = 0.0_wp
       com(:) = 0.0_wp
       iter = 0

       do n = 1, me%natoms
          if(me%atoms(n)%index == k) then
            sx = sx + me%atoms(n)%xyz(1)
            sy = sy + me%atoms(n)%xyz(2)
            sz = sz + me%atoms(n)%xyz(3)
            total_mass = total_mass + 1.0_wp
            Foo = me%atoms(n)%group
            iter = iter + 1
            na(k) = iter
          end if
       end do

       com(1) = sx / total_mass
       com(2) = sy / total_mass
       com(3) = sz / total_mass

       open(unit,status='scratch')
       write(unit,*) Foo, com(1), com(2), com(3)
       do n = 1, me%natoms
         if(me%atoms(n)%index == k) then
            write(unit,*) me%atoms(n)%sym, me%atoms(n)%xyz(:)-com(:)
         end if
       end do
    end do

    nmax = maxval(na)

    allocate(molName(group_count),rx(group_count),ry(group_count),rz(group_count))
    allocate(sym(group_count,nmax),x(group_count,nmax),y(group_count,nmax)&
            &,z(group_count,nmax))
    
    
    !!!@allocate(sites(nmax))

    rewind(unit)
    box(:) = me%box(:)

    minimal_image = minval(box)
    if(rcut > (minimal_image)/2.0_wp) then
      write(*,'(a)') 'Check minimal image condition: rcut<= Lmin/2'
      stop
    end if

    do k = 1, group_count
       read(unit,*) molName(k), rx(k), ry(k), rz(k) 
       do n = 1, na(k)
         read(unit,*) sym(k,n), x(k,n), y(k,n), z(k,n) 
       end do
    end do

    close(unit)
  end subroutine get_groups


  subroutine adjacency_pure(me,unit)
    type(frame), intent(inout) :: me
    integer  :: i, j
    real(wp) :: rxi, ryi, rzi
    real(wp) :: rxj, ryj, rzj
    real(wp) :: rxij, ryij, rzij
    real(wp) :: rijsq
    real(wp) :: com(3)
    real(wp) :: box(3)
    integer  :: ie, ii
    !
    real(wp) :: xe, ye, ze
    real(wp) :: xij, yij, zij
    real(wp) :: sr2

    !
    character(len=5) :: sa, sb, sc ! atomic symbols

    integer :: iter0, iter1
    logical :: current_cutoff 

    integer, intent(inout) :: unit

    logical :: current

    integer :: unit2
    logical :: scalar_dss
    
    current = .false.

    rcut2 = rcut * rcut

    box(:) = me%box(:)

    allocate(adj(group_count,group_count))
    adj = 0.0_wp

    iter0 = 0
    iter1 = 0

    do i = 1, group_count-1
       rxi = rx(i)
       ryi = ry(i)
       rzi = rz(i)
       do j = i+1, group_count
          current_cutoff = .false.

          rxj = rx(j)
          ryj = ry(j)
          rzj = rz(j)

          rxij = (rxi - rxj)
          ryij = (ryi - ryj)
          rzij = (rzi - rzj)

          com(1) = rxij
          com(2) = ryij
          com(3) = rzij
            

          ! call pbc to pick up central image
          if(pbc=='y') then
            call getpbc(com,box)
          end if  

            rxij = com(1)
            ryij = com(2)
            rzij = com(3)
            

            rijsq = rxij*rxij + ryij*ryij + rzij*rzij

          if(rijsq <= rcut2) then
              
              do ie = 1, na(i)
                 xe = x(i,ie)
                 ye = y(i,ie)
                 ze = z(i,ie)

                 sa = sym(i,ie) ! first atomic symbol

                 iter0 = iter0 + 1
                 do ii = 1, na(j)
                    xij = (rxij + (xe - x(j,ii)))
                    yij = (ryij + (ye - y(j,ii)))
                    zij = (rzij + (ze - z(j,ii)))
             
                    sb = sym(j,ii) ! second atomic symbol

                    iter1 = iter1 + 1

                    sr2 = xij*xij + yij*yij + zij*zij
                    
                    ! TODO: like-distance function
                    if(mode == 1) then
                      if(sa == pair(1) .and. sb == pair(2)) then
                        if(sr2 <= distance*distance) then
                           current = .true.
                           exit
                        end if
                      end if
                    end if
                    !--------------------------------------------
                    !> TODO: like-angle function
                    !>if(mode == 2) then
                    !>  if(sa == triad(1) .and. sb == triad(2) .and. sc == trian(3)) then
                    !>    if(turn <= angle) then
                    !>       dss(ie,ii,i,j)
                    !>    end if
                    !>  end if
                    !>end if
                    !---------------------------------------------
                    !> TODO: like-lattice function. This example suppose only two states
                    !>if(mode == 3) then
                    !>  if(spin(i,ie) == spin(j,ii)) then
                    !>    dss(ie,ii,i,j)      
                    !>  end if
                    !>end if

                 end do
                 if(current) exit
              end do
         end if !> rcut   
        if(current) then
          adj(i,j) = 1
          adj(j,i) = 1
        end if
        current = .false.
        end do
    end do

   !-------------------------------------------------------------
   ! Write CSV file
    do i = 1, group_count
        write(unit, '(i0)', advance='no') adj(i, 1)   
        do j = 2, group_count
            write(unit, '(a, i0)', advance='no') ',', adj(i, j)
        end do
        write(unit, *)
    end do
   !-------------------------------------------------------------

  end subroutine adjacency_pure      

  subroutine getpbc(com,box)
    real(wp), intent(inout) :: com(:)
    real(wp), intent(inout) :: box(:)

    com(1) = com(1) - box(1) * dnint(com(1) / box(1))
    com(2) = com(2) - box(2) * dnint(com(2) / box(2))
    com(3) = com(3) - box(3) * dnint(com(3) / box(3))
  end subroutine getpbc

!
  subroutine on_the_fly(me,unit)
    type(frame), intent(inout) :: me
    integer, intent(in) :: unit
    ! local variables
    integer       :: lines, ios
    character(80) :: buffer
    integer :: natoms, nf, nf_cut
    integer :: unit1, unit2
    logical :: is
    character(10) :: ext
    character(1)  :: show

    !--------------------------
    ! timing
    real(wp) :: time_start
    real(wp) :: time_finish
    !---------------------------

    character(len=:), allocatable :: cmd
    !cmd = 'python script.py '//trim(adjustl(measure))//' '//'measure'//' '//trim(adjustl(graph))

    show = 'n'
    lines = 0

    !@write(*,*) 'flag = ', frequency
    do
       read(unit,'(a)',iostat=ios) buffer
       if(ios/=0) exit
       lines=lines+1
    end do

    rewind(unit)
    read(unit,*)
    read(unit,*) natoms

    nframes=lines/(natoms+3)
    rewind(unit)
    
    !unit2 = 10

    nf_cut = 0
    
    do nf = 1, nframes
      close(unit2)
      call me%read(unit)
      unit1 = nf+20

      unit2 = 10
      inquire(unit=unit2,opened=is)
      do
        if(.not.is) exit
        inquire(unit=unit2,opened=is)
        unit2 = unit2+1
      end do

      open(unit2,file='adjacency.csv',status='unknown')

      call get_groups(me,unit1)

      !-----timing------------------
      call cpu_time(time_start)
      call adjacency_pure(me,unit2)
      call cpu_time(time_finish)
      write(*,'("Time elapsed for building the adjacency matrix in frame&
             & ",i0, ": "f6.3, " seconds")') nf, time_finish-time_start
      !-------------------------------
      call cpu_time(time_start)
      call build_script(nf)
      call cpu_time(time_finish)
      write(*,'("Building script - Time elapsed in frame&
             & ",i0, ": "f6.3, " seconds")') nf, time_finish-time_start

      !@call execute_command_line('bin/python script.py output.png')
      write(ext,'(i10)') nf

      if(mod(nf,frequency) == 0 .or. frequency == -1) then
        show = 'y'
      else
        show = 'n'
      end if

      cmd = 'time python script.py '//trim(adjustl(measure))//' '//'measure'//' '//&
       &trim(adjustl(graph))//'_'//trim(adjustl(ext))//'.png'//' '//&
       &trim(adjustl(active))//' '//trim(adjustl(show))
        
      write(*,*) '-----------------------------------------'
      write(*,'("Networkx - Time elapsed in frame&
             & ",i0)') nf
      !@call cpu_time(time_start)
      call execute_command_line(cmd)
      !@call cpu_time(time_finish)
      !@write(*,'("Networkx - Time elapsed in frame&
      !@       & ",i0, ": "f6.3, " seconds")') nf, time_finish-time_start

      close(unit2)
      deallocate(me%atoms)
      deallocate(na)
      deallocate(molName,rx,ry,rz)
      deallocate(sym,x,y,z)
      deallocate(adj)
      !!@deallocate(dss)
      !@deallocate(sites)

      write(*,*) '-----------------------------------'
    end do

  !@  if(mod(nframes,frequency)==0) then
  !@     call analysis_interface(nf_cut)
  !@  else
       call analysis_interface(nf)
  !@  end if

  end subroutine on_the_fly

!--------Analysis--------------------------------------  
    subroutine analysis_interface(nf)
      integer, intent(in) :: nf
      character(len=:), allocatable :: results
      character(len=:), allocatable :: input
      integer :: n
      character(5) :: Foo
      integer :: unit1, unit2
      integer :: ios
      character(257) :: line
      character(len=:), allocatable ::cmd

      results = 'All_results.txt'

      close(57)
      open(unit=57,file=results,status='unknown')

      do n = 1, nf-1
        write(Foo,'(i5)') n
        input = measure//"_"//trim(adjustl(Foo))//".txt"

        close(56)
        open(unit=56,file=input,status='old')
       

        do 
          read(56,'(a)',iostat=ios) line
          if(ios /= 0) exit
          write(57,'(a)') trim(adjustl(line))
        end do

        close(56)
      end do

      do n = 1, nf-1
        write(Foo,'(i5)') n
        input = measure//"_"//trim(adjustl(Foo))//".txt"
        cmd = "rm -r "//input

       call execute_command_line(cmd)
      end do

      write(*,'(a)') "-------------------------"
      write(*,'(a)') "COOR2GRAPH Version 1.0.0"
      write(*,'(a)') ""
      write(*,'(a)') measure//" of graph"
      write(*,'(a)') "data save in:"
      write(*,'(a)') "prefix_*.png"
      write(*,'(a)') "histogram.png"
      write(*,'(a)') "-------------------------"


      ! ------------------------------------------------
      close(93)
      open(93,file='histogram.py',status='unknown')

      write(93,'(a)')"import numpy as np"
      write(93,'(a)')"import matplotlib.pyplot as plt"
      write(93,'(a)')"import sys"
      write(93,'(a)')"import os"
      write(93,'(a)')""

      write(93,'(a)')"file_name = 'All_results.txt'"
      write(93,'(a)')"data = np.loadtxt(file_name, usecols=1)"
      write(93,'(a)')""

      write(93,'(a)')"average = np.mean(data)"
      write(93,'(a)')"std_dev = np.std(data)"
      write(93,'(a)')"print(f'Average: {average:.4f}')"
      write(93,'(a)')"print(f'Standard Deviation: {std_dev:.4f}')"

      write(93,'(a)')""
!@      write(93,'(a)')"sys.stdout = open(os.devnull, 'w')"

      write(93,'(a)')"data_size = len(data)"

      write(93,'(a)')"bin_count=int(np.ceil(1 + np.log2(data_size)))"
      write(93,'(a)')"hist, bin_edges = np.histogram(data, bins=bin_count)"
      
      write(93,'(a)')"print('Bin Start  Bin End  Frequency')"
      write(93,'(a)')""

      write(93,'(a)')"for start, end, freq in zip(bin_edges[:-1], bin_edges[1:], hist):"
      write(93,'(a)')"    print(f'{start:10.4f} {end:10.4f} {freq:10d}')"
      write(93,'(a)')""

      write(93,'(a)')"plt.figure(figsize=(8,6))"
      write(93,'(a)')"plt.hist(data, bins=bin_count, edgecolor='black', alpha=0.7)"
      write(93,'(a)')"plt.title('Frequency of "//measure//"')"
      write(93,'(a)')"plt.xlabel('Bins')"
      write(93,'(a)')"plt.ylabel('Frequency')"
      write(93,'(a)')"plt.grid(axis='y', linestyle='--', alpha=0.7)"
      write(93,'(a)')"plt.tight_layout()"

      write(93,'(a)')"plt.savefig('histogram.png')"
      !@write(93,'(a)')"plt.show()"

      cmd = "python histogram.py"
      call execute_command_line(cmd)


    end subroutine analysis_interface


!--------user-defined functions ----------------------

  subroutine pair_like(current_cutoff,sr2,check,sa,sb)
    logical,intent(inout) :: check
    logical,intent(inout) :: current_cutoff
    real(wp),intent(in)   :: sr2
    character(*),intent(in) :: sa, sb
    check = .false.

    current_cutoff = .false.
    if(sa==pair(1) .and. sb==pair(2)) then
      if(sr2<=distance*distance) then
        current_cutoff = .true.
!@        check = .true.
      end if
    end if
  end subroutine pair_like

!--------end user-defined functions--------------------

end module frame_module

