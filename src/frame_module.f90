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
  real(wp), allocatable :: rx(:), ry(:), rz(:)
  real(wp), allocatable :: x(:,:), y(:,:), z(:,:)
  integer, allocatable :: na(:)
  !-------------------------------------------------
  real(wp) :: rcut2
  integer, allocatable :: adj(:,:)

  ! test ising
  integer :: spin(10000,1)
  integer :: id(10000)
  !-----------------------------

  ! ---- check ---
  logical :: check

contains

  subroutine read_user()
    integer :: unit
    logical :: is

    unit = 80
    inquire(unit=unit,opened=is)

    do 
       if(.not.is) exit
       inquire(unit=unit,opened=is)
       unit=unit+1
    end do
  end subroutine read_user      

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

    ! ising test
    integer :: sp(9)
    integer :: iter5
    
    
  !@  ! ising test
     open(97,file='ising.dat',status='old')
     do k = 1, 9
        read(97,*) sp(k)
     end do
    !-----------------------------


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
    allocate(sym(group_count,nmax),x(group_count,nmax),y(group_count,nmax),z(group_count,nmax))

    rewind(unit)
    box(:) = me%box(:)

    do k = 1, group_count
       read(unit,*) molName(k), rx(k), ry(k), rz(k) 
       do n = 1, na(k)
         read(unit,*) sym(k,n), x(k,n), y(k,n), z(k,n) 
       end do
    end do

    iter5 = 1
    do k = 1, group_count
       id(k) = k
       do n = 1, na(k)
          spin(k,n) = sp(iter5)

     !@     write(*,*) spin(k,n)
          iter5 = iter5 + 1
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

    integer :: iter0, iter1
    logical :: current_cutoff 

    integer, intent(inout) :: unit

    !----------------------------

    rcut2 = rcut * rcut

    box(:) = me%box(:)

    allocate(adj(me%natoms,me%natoms))
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
          if(with_pbc=='y') then
            call pbc(com,box)
          end if  
          rijsq = rxij*rxij + ryij*ryij + rzij*rzij

          rxij = com(1)
          ryij = com(2)
          rzij = com(3)

          do ie = 1, na(i)
             xe = x(i,ie)
             ye = y(i,ie)
             ze = z(i,ie)

             iter0 = iter0 + 1
             do ii = 1, na(j)
                xij = (rxij + (xe - x(j,ii)))
                yij = (ryij + (ye - y(j,ii)))
                zij = (rzij + (ze - z(j,ii)))
                iter1 = iter1 + 1

                sr2 = xij*xij + yij*yij + zij*zij

                if(sr2 <= rcut2) then
                  select case(mode)
                    case(0)
                      call ising_like(i,j,ie,ii,current_cutoff,check)
                      if(check) exit
                    case(1)  
                    case(2)  
                    case(3)  
                    case(4)  
                    case(5)  
                    case(6)
                      call distance_like(current_cutoff,sr2,check)
                      if(check) exit
                    case default
                      call ising_like(i,j,ie,ii,current_cutoff,check)
                      if(check) exit
                  end select


                  ! ising test
                 !!  if(spin(j,ii) == spin(i,ie) )  then 
                 !@@     current_cutoff = .true.
              !@!       exit
                 !!  end if
                end if
             end do
             if(current_cutoff) exit
          end do
            if(current_cutoff) then
              adj(i,j) = 1
              adj(j,i) = 1
            end if
       end do
    end do

!@    do i = 1, group_count
!@      write(*,*) (adj(i,j), j = 1, group_count)
!@    end do

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

  subroutine pbc(com,box)
    real(wp), intent(inout) :: com(:)
    real(wp), intent(inout) :: box(:)

    com(1) = com(1) - box(1) * dnint(com(1) / box(1))
    com(2) = com(2) - box(2) * dnint(com(2) / box(2))
    com(3) = com(3) - box(3) * dnint(com(3) / box(3))
  end subroutine pbc

!
  subroutine on_the_fly(me,unit)
    type(frame), intent(inout) :: me
    integer, intent(in) :: unit
    ! local variables
    integer       :: lines, ios
    character(80) :: buffer
    integer :: natoms, nf
    integer :: unit1, unit2
    logical :: is
    character(10) :: ext

    !---------------------------

    character(len=:), allocatable :: cmd
    !cmd = 'python script.py '//trim(adjustl(measure))//' '//'measure'//' '//trim(adjustl(graph))

    lines = 0
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
      if(pure =='y') call adjacency_pure(me,unit2)
      call build_script()

      !@call execute_command_line('bin/python script.py output.png')
      write(ext,'(i10)') nf
      cmd = 'python script.py '//trim(adjustl(measure))//' '//'measure'//' '//trim(adjustl(graph))//'_'//trim(adjustl(ext))//'.png'
      call execute_command_line(cmd)

      close(unit2)
      deallocate(me%atoms)
      deallocate(na)
      deallocate(molName,rx,ry,rz)
      deallocate(sym,x,y,z)
      deallocate(adj)

      write(*,*) '-----------------------------------'
    end do

  end subroutine on_the_fly

!--------user-defined functions ----------------------

!@  function only_distance()
!@  end function only_distance

!@  function only_angle()
!@  end function only_angle  

!@  function distance_angle()
!@  end function distance_angle

!@  function only_com()
!@  end function only_com

!@  function pair_atomic_labels()
!@  end function pair_atomic_labels

!@  function triad_atomic_labels()
!@  end function triad_atomic_labels

  ! Ising test
  subroutine ising_like(i,j,ie,ii,current_cutoff,check)
    integer, intent(in) :: i,j,ie,ii
    logical, intent(inout) :: check
    logical, intent(inout) :: current_cutoff
    check = .false.
    if(spin(j,ii) == spin(i,ie)) then
      current_cutoff = .true.
      check = .true.
    end if 
  end subroutine ising_like

  subroutine distance_like(current_cutoff,sr2,check)
    logical, intent(inout) :: check
    logical, intent(inout) :: current_cutoff
    real(wp), intent(in)   :: sr2
    check = .false.
    if(sr2<=cdistance) then
      current_cutoff = .true.
      check = .true.
    end if
  end subroutine distance_like

!--------end user-defined functions--------------------

end module frame_module

