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
  !-------------------------------------------------
  real(wp) :: rcut2
  integer, allocatable :: adj(:,:)
  logical :: check

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
          if(pbc=='y') then
            call getpbc(com,box)
          end if  
          rijsq = rxij*rxij + ryij*ryij + rzij*rzij

          rxij = com(1)
          ryij = com(2)
          rzij = com(3)

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

                if(sr2 <= rcut2) then
                  select case(mode)
                    case(1)
                    call pair_like(current_cutoff,sr2,check,sa,sb)
                    if(check) exit
                  end select
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
      call adjacency_pure(me,unit2)
      call build_script(nf)

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

  subroutine pair_like(current_cutoff,sr2,check,sa,sb)
    logical,intent(inout) :: check
    logical,intent(inout) :: current_cutoff
    real(wp),intent(in)   :: sr2
    character(*),intent(in) :: sa, sb
    check = .false.
    if(sa==pair(1) .and. sb==pair(2)) then
      if(sr2<=distance*distance) then
        current_cutoff = .true.
        check = .true.
      end if
    end if
  end subroutine pair_like

!--------end user-defined functions--------------------

end module frame_module

