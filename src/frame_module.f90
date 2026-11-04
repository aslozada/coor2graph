!-------------------------------------------------------------------
! This module is part of coor2graph
! Frame and trajectory processing, group identification, and adjacency matrix construction from coordinate data
!-------------------------------------------------------------------
!> author: Asdrubal Lozada
!  update: 4/11/2026
!  license: GPL-v3 

module frame_module
  use kind_module
  use utils_module
  use commands_module
  implicit none

  private

  type, public :: atom_properties
    real(wp), allocatable :: values(:)
  end type atom_properties

  type, public :: atom
    integer      :: index
    character(5) :: group
    character(5) :: sym
    integer      :: seq
    real(wp)     :: xyz(3)
    type(atom_properties) :: properties
  end type atom

  type, public :: frame
    type(atom), allocatable :: atoms(:)
    integer  :: natoms
    real(wp) :: box(3)
  contains
    procedure :: read => read_frame
  end type frame

  integer, public :: nframes
  public :: process_trajectory
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
    integer       :: n, ios
    integer       :: natoms

    read(unit,*,iostat=ios) 
    if (ios /= 0) call fail('error reading frame title')
    read(unit,*,iostat=ios) natoms
    if (ios /= 0) call fail('error reading number of atoms')

    if(natoms <= 0) call fail('invalid number of atoms')
    if (allocated(me%atoms)) deallocate(me%atoms)
    allocate(me%atoms(natoms))

    me%natoms = natoms
    
    do n = 1, me%natoms
      read(unit,'(i5,2a5,i5,3f8.3)',iostat=ios) &
      me%atoms(n)%index,me%atoms(n)%group,&
      me%atoms(n)%sym,me%atoms(n)%seq,&
      me%atoms(n)%xyz(:)
    end do

    read(unit,*,iostat=ios) me%box(:)
    if (ios /= 0) call fail('error reading simuation box')
  end subroutine read_frame

subroutine get_groups(me,unit)

  type(frame), intent(inout) :: me
  integer, intent(in) :: unit

  ! Local variables
  integer :: n, k
  integer :: nmax, count_atoms
  integer :: tmp_unit, ios

  integer, allocatable :: group_id(:)

  real(wp) :: sx, sy, sz
  real(wp) :: com(3)
  real(wp) :: box(3)
  real(wp) :: minimal_image

  character(5) :: Foo

  ! Identify unique group indices
  allocate(group_id(me%natoms))

  group_count = 0

  do n = 1, me%natoms

    if (.not. any(group_id(1:group_count) == me%atoms(n)%index)) then
      group_count = group_count + 1
      group_id(group_count) = me%atoms(n)%index
    end if

  end do

  if (group_count == 0) call fail('no groups found')

  ! Count atoms in each group
  allocate(na(group_count))
  na = 0

  do k = 1, group_count
    do n = 1, me%natoms

      if (me%atoms(n)%index == group_id(k)) then
        na(k) = na(k) + 1
      end if

    end do
  end do

  ! Open temporary file
  open(newunit=tmp_unit, status='scratch', iostat=ios)
  if (ios /= 0) call fail('cannot open temporary file')

  ! Compute geometrical centers and relative coordinates
  do k = 1, group_count

    sx = 0.0_wp
    sy = 0.0_wp
    sz = 0.0_wp

    count_atoms = 0

    do n = 1, me%natoms

      if (me%atoms(n)%index == group_id(k)) then

        sx = sx + me%atoms(n)%xyz(1)
        sy = sy + me%atoms(n)%xyz(2)
        sz = sz + me%atoms(n)%xyz(3)

        Foo = me%atoms(n)%group
        count_atoms = count_atoms + 1

      end if

    end do

    if (count_atoms == 0) call fail('empty group detected')

    com(1) = sx / real(count_atoms,wp)
    com(2) = sy / real(count_atoms,wp)
    com(3) = sz / real(count_atoms,wp)

    write(tmp_unit,*,iostat=ios) Foo, com(:)
    if (ios /= 0) call fail('error writing temporary group data')

    do n = 1, me%natoms

      if (me%atoms(n)%index == group_id(k)) then

        write(tmp_unit,*,iostat=ios) me%atoms(n)%sym, &
                                     me%atoms(n)%xyz(:)-com(:)

        if (ios /= 0) call fail('error writing temporary atom data')

      end if

    end do

  end do

  ! Allocate group properties
  nmax = maxval(na)

  allocate(molName(group_count), rx(group_count), &
           ry(group_count), rz(group_count))

  allocate(sym(group_count,nmax), &
           x(group_count,nmax), &
           y(group_count,nmax), &
           z(group_count,nmax))

  ! Check minimum image condition
  box(:) = me%box(:)

  if (pbc == 'y') then

    if (any(box <= 0.0_wp)) call fail('invalid simulation box')

    minimal_image = minval(box)

    if (rcut > minimal_image/2.0_wp) then
      call fail('rcut violates the minimum image convention')
    end if

  end if

  ! Read temporary data
  rewind(tmp_unit)

  do k = 1, group_count

    read(tmp_unit,*,iostat=ios) molName(k), rx(k), ry(k), rz(k)
    if (ios /= 0) call fail('error reading temporary group data')

    do n = 1, na(k)

      read(tmp_unit,*,iostat=ios) sym(k,n), x(k,n), y(k,n), z(k,n)
      if (ios /= 0) call fail('error reading temporary atom data')

    end do

  end do

  close(tmp_unit)

  deallocate(group_id)

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
  subroutine process_trajectory(me,unit)
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

    show = 'n'
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

      write(ext,'(i10)') nf

      open(unit2,file='adjacency_'//trim(adjustl(ext))//'.csv',status='unknown')

      call get_groups(me,unit1)

      !-----timing------------------
      call cpu_time(time_start)
      call adjacency_pure(me,unit2)
      call cpu_time(time_finish)
      write(*,'("Time elapsed for building the adjacency matrix in frame&
             & ",i0, ": "f6.3, " seconds")') nf, time_finish-time_start
      !-------------------------------
      close(unit2)
      deallocate(me%atoms)
      deallocate(na)
      deallocate(molName,rx,ry,rz)
      deallocate(sym,x,y,z)
      deallocate(adj)

      write(*,*) '-----------------------------------'
    end do

  end subroutine process_trajectory



end module frame_module

