
module parameters
#ifdef MPI
  use mpi
#endif
use modcom
implicit none
private
logical :: atoms_block_found=.false.
integer, parameter :: nmaxtasks=10
integer, parameter :: nlines_max=100000
!integer, parameter :: dp = SELECTED_REAL_KIND (15,300)

type, public :: CLpars
  character(len=100) :: input_file
  character(len=100) :: wannier_proj_mode=""
  character(len=100) :: tasks(nmaxtasks)
  character(len=100) :: geometry_source=""
  character(len=100) :: seedname="seedname"
  character(len=100) :: sktype="sk"
  logical :: sparse=.false.
  integer :: nwan=0
  integer :: geometry_index=0
  integer :: symtype=0
  integer :: nvert
  integer :: nspec
  integer :: istart=1
  integer :: istop=1000
  integer :: nstates=1000
  integer :: negrid=100
  integer :: ntasks=0
  integer :: natmtot
  integer :: nmaxatm_pspec
  integer :: ngrid(NDIM)
  integer :: n_valence=1
  real(dp) :: gauss_sigma=0.1_dp
  real(dp) :: sparse_eps=0.e-6_dp
  real(dp) :: efermi=0._dp
  real(dp) :: emin=-10._dp
  real(dp) :: emax=10._dp
  real(dp) :: rcut_nn=100._dp
  real(dp) :: avec(NDIM,NDIM)
  real(dp) :: bvec(NDIM,NDIM)
  integer, allocatable :: nat_per_spec(:)
  integer, allocatable :: norb_per_center(:)
  integer, allocatable :: np_per_vert(:)
  integer, allocatable :: tot_iais(:,:)
  integer, allocatable :: iais_tot(:,:)
  real(dp), allocatable :: egrid(:)
  real(dp), allocatable :: atml(:,:,:)
  real(dp), allocatable :: wannier_axis(:,:,:)
  real(dp), allocatable :: vert(:,:)
  contains
  procedure :: init=>read_input
  procedure :: atmc=>calc_atmc
endtype CLpars

contains

subroutine read_input(THIS)
use geometry_library
class(CLpars), intent(inout) :: THIS
! internal
type(geomlib) geometry
integer iostat,iline,jline,ii
integer ispec,iat,ivert,igrid
integer, parameter :: nmaxatm_pspec=40000
real(dp) t1,tvec(NDIM,NDIM)
character(len=256) block,arg,line
real(dp), allocatable :: atml_temp(:,:,:)

#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif

open(50,file=trim(adjustl(THIS%input_file)),action="read",status="old",iostat=iostat)
if (iostat.ne.0) call throw("paramters%read_input()","could not open input file")
if (mp_mpi) call message("")
if (mp_mpi) call message("Reading input file: ")

jline=0
do iline=1,nlines_max
  jline=jline+1
  read(50,'(A)',end=100,iostat=iostat) line
  if (iostat.ne.0) call throw("paramters%read_input()","could not read line in the input file")
  if (mp_mpi) write(*,'(i6,"| ",A)') jline,trim(adjustl(line))

  call split_string(line,block,arg)

  ! tasks block
  if (trim(block).eq."tasks") then
     read(arg,*,iostat=iostat) THIS%ntasks
     if (iostat.ne.0) call throw("paramters%read_input()","problem with ntasks argument")
     do ii=1,THIS%ntasks
        jline=jline+1
        read(50,'(A)',iostat=iostat) line
        if (iostat.ne.0) call throw("paramters%read_input()","problem with tasks block")
        if (mp_mpi) write(*,'(i6,": ",A)') jline,trim(adjustl(line))
        call split_string(line,block,arg)
        THIS%tasks(ii)=trim(adjustl(block))
     end do

  ! real(read) and reciprocal (computed) lattice vectors
  else if (trim(block).eq."avec") then
    read(arg,*,iostat=iostat) t1
    if (iostat.ne.0) call throw("paramters%read_input()","problem with avec's scale argument")
    do ii=1,NDIM
      jline=jline+1
      read(50,*,iostat=iostat) THIS%avec(ii,:)
      if (iostat.ne.0) call throw("paramters%read_input()","problem with avec block")
      if(mp_mpi) write(*,'(i6,": ",5F10.6)') jline,THIS%avec(ii,:)
    end do
    tvec=THIS%avec
    call dmatrix_inverse(tvec,THIS%bvec,NDIM)
    THIS%bvec=transpose(THIS%bvec)*twopi
    tvec=matmul(THIS%avec,transpose(THIS%bvec))
    do ii=1,NDIM
      if (abs(tvec(ii,ii)-twopi).gt. epslat) then
        call throw("paramters%read_input()","could not construct reciprocal lattice vectors")
      end if
    end do

  ! atoms coordinates and species
  else if (trim(block).eq."atoms") then
    atoms_block_found=.true.
    read(arg,*,iostat=iostat) THIS%nspec
    if (iostat.ne.0) call throw("paramters%read_input()","problem with atoms block nspec argument")
    allocate(THIS%nat_per_spec(THIS%nspec))
    allocate(atml_temp(NDIM,nmaxatm_pspec,THIS%nspec))
    THIS%nmaxatm_pspec=0
    do ispec=1,THIS%nspec
      jline=jline+1
      read(50,*,iostat=iostat) THIS%nat_per_spec(ispec)
      if (iostat.ne.0) call throw("paramters%read_input()","problem with nat_per_spec field of atoms block")
      if (mp_mpi) write(*,'(i6,": ",i6)') jline,THIS%nat_per_spec(ispec)
      THIS%nmaxatm_pspec=max(THIS%nmaxatm_pspec,THIS%nat_per_spec(ispec))
      if (THIS%nmaxatm_pspec.gt.nmaxatm_pspec)  call throw("paramters%read_input()","nmaxatm_pspec values exceedid in atoms block")
      do iat=1,THIS%nat_per_spec(ispec)
         jline=jline+1
         read(50,*,iostat=iostat) atml_temp(:,iat,ispec)
         if (iostat.ne.0) call throw("paramters%read_input()","problem with atom coordinates fields")
         if (mp_mpi) write(*,'(i6,": ",5F10.6)') jline,atml_temp(:,iat,ispec)
      end do
    end do
    allocate(THIS%atml(NDIM,THIS%nmaxatm_pspec,THIS%nspec))
    do ispec=1,THIS%nspec
      do iat=1,THIS%nat_per_spec(ispec)
        THIS%atml(:,iat,ispec)=atml_temp(:,iat,ispec)
      end do
    end do
    deallocate(atml_temp)

  ! BZ k-point grid
  else if (trim(block).eq."ngrid") then
    jline=jline+1
    read(50,*,iostat=iostat) THIS%ngrid(:)
    if (iostat.ne.0) call throw("paramters%read_input()","problem with ngrid data")
    if (mp_mpi) write(*,'(i6,": ",5I6)') jline,THIS%ngrid(:)

  ! BZ k-papth block
  else if (trim(block).eq."path") then
    read(arg,*,iostat=iostat) THIS%nvert
    if (iostat.ne.0) call throw("paramters%read_input()","problem with path's nvert argumet")
    allocate(THIS%np_per_vert(THIS%nvert))
    allocate(THIS%vert(NDIM,THIS%nvert))
    do ivert=1,THIS%nvert
      jline=jline+1
      read(50,*,iostat=iostat) THIS%vert(:,ivert),THIS%np_per_vert(ivert)
      if (iostat.ne.0) call throw("paramters%read_input()","problem with path data")
      if (mp_mpi) write(*,'(i6,": ",5F10.6)',advance='no') jline,THIS%vert(:,ivert)
      if (mp_mpi) write(*,'(I6)') THIS%np_per_vert(ivert)
    end do

  ! BZ k-papth block
  else if (trim(block).eq."geometry_source") then
    read(arg,*,iostat=iostat) THIS%geometry_source
    if (iostat.ne.0) call throw("paramters%read_input()","problem with argument of geometry_source block")
    if (THIS%geometry_source.eq."tbg") then
      jline=jline+1
      read(50,*,iostat=iostat) THIS%geometry_index
      if (iostat.ne.0) call throw("paramters%read_input()","problem with geometry_index of geometry_source block")
      if (mp_mpi) write(*,'(i6,": ",i6)') jline,THIS%geometry_index
    end if

  ! cut off for nearest neighbors
  else if (trim(block).eq."rcut_nn") then
    jline=jline+1
    read(50,*,iostat=iostat) THIS%rcut_nn
    if (iostat.ne.0) call throw("paramters%read_input()","problem with rcut_nn data")
    if (mp_mpi) write(*,'(i6,": ",F10.6)') jline,THIS%rcut_nn

  ! Number of valence orbitals
  else if (trim(block).eq."n_valence") then
    jline=jline+1
    read(50,*,iostat=iostat) THIS%n_valence
    if (iostat.ne.0) call throw("paramters%read_input()","problem with n_valence data")
    if (mp_mpi) write(*,'(i6,": ",i6)') jline,THIS%n_valence

  ! gauss_sigma for smearing of the DOS, for example
  else if (trim(block).eq."gauss_sigma") then
    jline=jline+1
    read(50,*,iostat=iostat) THIS%gauss_sigma
    if (iostat.ne.0) call throw("paramters%read_input()","problem with gauss_sigma data")
    if (mp_mpi) write(*,'(i6,": ",F10.6)') jline,THIS%gauss_sigma

  ! proposed Fermi energy
  else if (trim(block).eq."efermi") then
    jline=jline+1
    read(50,*,iostat=iostat) THIS%efermi
    if (iostat.ne.0) call throw("paramters%read_input()","problem with efermi data")
    if (mp_mpi) write(*,'(i6,": ",F10.6)') jline,THIS%efermi

  ! proposed Fermi energy
  else if (trim(block).eq."sparse_eps") then
    jline=jline+1
    read(50,*,iostat=iostat) THIS%sparse_eps
    if (iostat.ne.0) call throw("paramters%read_input()","problem with sparse_eps data")
    if (mp_mpi) write(*,'(i6,": ",F18.10)') jline,THIS%sparse_eps

  ! .true. to use the sparse algorithms
  else if (trim(block).eq."sparse") then
    jline=jline+1
    read(50,*,iostat=iostat) THIS%sparse
    if (iostat.ne.0) call throw("paramters%read_input()","problem with sparse data")
    if (mp_mpi) write(*,'(i6,": ",L6)') jline,THIS%sparse


  ! energy grid for DOS or spectral functions
  else if (trim(block).eq."egrid") then
    read(arg,*,iostat=iostat) THIS%negrid
    if (iostat.ne.0) call throw("paramters%read_input()","problem with egrid argumet")
    read(50,*,iostat=iostat) THIS%emin,THIS%emax
    if (iostat.ne.0) call throw("paramters%read_input()","problem with egrid data")
    if (THIS%emin.ge.THIS%emax) call throw("paramters%read_input()","emin<emax in egrid block")
    if (mp_mpi) write(*,'(i6,": ",5F10.6)') jline,THIS%emin,THIS%emax

  else if (trim(block).eq."states") then
    jline=jline+1
    read(50,*,iostat=iostat) THIS%istart,THIS%istop
    if (iostat.ne.0) call throw("paramters%read_input()","problem with states data")
    if (mp_mpi) write(*,'(i6,": ",2I6)') jline,THIS%istart,THIS%istop
    if (THIS%istop.lt.THIS%istart) call throw("paramters%read_input()","istart state should be less or equal then istop")
    THIS%nstates=THIS%istop-THIS%istart+1

  ! seedname for Wannier90 interface
  else if (trim(block).eq."seedname") then
        jline=jline+1
        read(50,'(A)',iostat=iostat) THIS%seedname
        if (iostat.ne.0) call throw("paramters%read_input()","problem with seedname data")
        if (mp_mpi) write(*,'(i6,": ",A)') jline,trim(adjustl(THIS%seedname))

  ! type of Slater-Koster function
  else if (trim(block).eq."sktype") then
        jline=jline+1
        read(50,'(A)',iostat=iostat) THIS%sktype
        if (iostat.ne.0) call throw("paramters%read_input()","problem with sktype data")
        if (mp_mpi) write(*,'(i6,": ",A)') jline,trim(adjustl(THIS%seedname))

  ! projection mode for wannier export
  else if (trim(block).eq."wannier_proj_mode") then
        jline=jline+1
        read(50,'(A)',iostat=iostat) THIS%wannier_proj_mode
        if (iostat.ne.0) call throw("paramters%read_input()","problem with wannier_proj_mode data")
        if (mp_mpi) write(*,'(i6,": ",A)') jline,trim(adjustl(THIS%wannier_proj_mode))

  ! number of disired wannier projection
  else if (trim(block).eq."nwan") then
    jline=jline+1
    read(50,*,iostat=iostat) THIS%nwan
    if (iostat.ne.0) call throw("paramters%read_input()","problem with nwan data")
    if (mp_mpi) write(*,'(i6,": ",I6)') jline,THIS%nwan

  end if

end do
100 continue
close(50)
if (mp_mpi) write(*,*)
#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
if (trim(adjustl(THIS%geometry_source)).ne."") then
   call geometry%init(THIS%geometry_index,THIS%geometry_source)
   if (allocated(THIS%nat_per_spec)) deallocate(THIS%nat_per_spec)
   if (allocated(THIS%atml)) deallocate(THIS%atml)
   allocate(THIS%nat_per_spec(geometry%nspec))
   allocate(THIS%atml(NDIM,geometry%nmaxatm_pspec,geometry%nspec))
   THIS%avec=geometry%avec
   THIS%bvec=geometry%bvec
   THIS%nspec=geometry%nspec
   THIS%nmaxatm_pspec=geometry%nmaxatm_pspec
   THIS%nat_per_spec=geometry%nat_per_spec
   THIS%atml=geometry%atml
end if
! compute total number of atoms, and construct mapping to/from total index
THIS%natmtot=0
do ispec=1,THIS%nspec
  do iat=1,THIS%nat_per_spec(ispec)
    THIS%natmtot=THIS%natmtot+1
  end do
end do
if (trim(adjustl(THIS%geometry_source)).eq."tbg".or.&
    trim(adjustl(THIS%geometry_source)).eq."slg") then

            allocate(THIS%norb_per_center(THIS%natmtot))
            allocate(THIS%wannier_axis(NDIM,2,THIS%natmtot))
            THIS%norb_per_center=1
            if (NDIM.eq.3) then
              ! X and Z axis
              do iat=1,THIS%natmtot
                THIS%wannier_axis(:,1,iat)=(/1._dp,0._dp,0._dp/)
                THIS%wannier_axis(:,2,iat)=(/0._dp,0._dp,1._dp/)
              end do
            else
               call throw("paramters%read_input()","wannier axis assignemt works in 3D case only")
            end if
else
   call throw("paramters%read_input()","unknown geometry structure option")
end if
allocate(THIS%tot_iais(THIS%natmtot,2))
allocate(THIS%iais_tot(THIS%nmaxatm_pspec,THIS%nspec))
THIS%natmtot=0
do ispec=1,THIS%nspec
  do iat=1,THIS%nat_per_spec(ispec)
    THIS%natmtot=THIS%natmtot+1
    THIS%tot_iais(THIS%natmtot,1)=iat
    THIS%tot_iais(THIS%natmtot,2)=ispec
    THIS%iais_tot(iat,ispec)=THIS%natmtot
  end do
end do
if (THIS%natmtot.le.0) call throw("paramters%read_input()","no atoms!")
! allocate default/read egrid
allocate(THIS%egrid(THIS%negrid))
if (THIS%negrid.le.1) then
  THIS%egrid(1)=0.5_dp*(THIS%emin+THIS%emax)
else
  do igrid=1,THIS%negrid
    THIS%egrid(igrid)=THIS%emin+(THIS%emax-THIS%emin)*dble(igrid-1)/dble(THIS%negrid-1)
  end do
end if
#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
! output info about sparse solver
if (THIS%sparse) then
  call info("CLpars%read_input()","sparse flag is .true.")
#ifdef PARDI
  call message("-DPARDI precomiled option included, the sparse solver will run")
  call message("in junction with the sparse MKL linear solver")
#else
   call message("WARNING: -DPARDI precomiled option is NOT included, this results to the fact")
   call message("that the sparse solver is run in junction with the DENSE linear solver, which needs")
   call message("to allocate the full matrix, therefore larger memory consumption expected.")
   call message("However, it is still very fast")
#endif
end if
! write atoms geometry to a file
if (mp_mpi) then
  open(100,file="geometry.dat",action="write")
  write(100,*)"avec"
  do ii=1,NDIM
    write(100,'(5G18.10)') THIS%avec(ii,:)
  end do
  write(100,*)"bvec"
  do ii=1,NDIM
    write(100,'(5G18.10)') THIS%bvec(ii,:)
  end do
  write(100,'("atoms ",I6)') THIS%nspec
  do ispec=1,THIS%nspec
  write(100,'(I6)') THIS%nat_per_spec(ispec)
    do iat=1,THIS%nat_per_spec(ispec)
      write(100,'(5G18.10)') THIS%atml(:,iat,ispec)
    end do
  end do
  close(100)
  write(*,*)
end if

#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif

end subroutine

function calc_atmc(THIS,iat,ispec) result(vpc)
class(CLpars), intent(in) :: THIS
integer, intent(in) :: iat,ispec
real(dp) vpc(NDIM)
vpc=matmul(THIS%atml(:,iat,ispec),THIS%avec)
end function

end module
