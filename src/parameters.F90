
module parameters
#ifdef MPI
  use mpi
#endif
use modcom
implicit none
private
logical :: atoms_block_found=.false.
logical :: avec_block_found=.false.
logical :: kshift_block_found=.false.
integer, parameter :: nmaxtasks=10
integer, parameter :: nlines_max=100000
!integer, parameter :: dp = SELECTED_REAL_KIND (15,300)
type :: CLproj
  logical :: allocatd=.false.
  integer :: norb=0
  integer :: ncenters=0
  integer, allocatable :: iw2ic(:)
  integer, allocatable :: norb_ic(:)
  integer, allocatable :: lmr(:,:)
  real(dp), allocatable :: waxis(:,:,:)
  real(dp), allocatable :: centers(:,:)
  contains
  procedure :: write_base
endtype

type, public :: CLpars
  character(len=100) :: input_file
  character(len=100) :: wannier_proj_mode=""
  character(len=100) :: tasks(nmaxtasks)
  character(len=100) :: geometry_source=""
  character(len=100) :: seedname="seedname"
  character(len=100) :: sktype="sk"
  character(len=100) :: tbfile=""
  character(len=100) :: character_file=""
  character(len=100) :: tbtype="sk"
  character(len=100) :: coulrs_file=""
  logical :: HubU_diagonal=.false.
  logical :: readsym=.false.
  logical :: shifted=.false.
  logical :: symtshift=.true.
  logical :: writetb=.false.
  logical :: sparse=.false.
  logical :: ignore_chiIq=.false.
  logical :: chi_exclude=.false.
  logical :: use_weights_amn=.false.
  integer :: niter_symmetrize=1
  integer :: ndim_coul=3
  integer :: geometry_index=0
  integer :: symtype=1
  integer :: nvert=0
  integer :: nspec
  integer :: iflat_band=0
  integer :: istart=1
  integer :: istop=1000
  integer :: nstates=1000
  integer :: chi_start=1
  integer :: chi_stop=1000
  integer :: negrid=1
  integer :: ntasks=0
  integer :: nkshift=1
  integer :: natmtot
  integer :: nmaxatm_pspec
  integer :: ngrid(NDIM)
  integer :: qgrid(NDIM)
  integer :: Ggrid(NDIM)
  real(dp) :: omega
  real(dp) :: omegabz
  real(dp) :: gauss_sigma=0.1_dp
  real(dp) :: sparse_eps=0.e-6_dp
  real(dp) :: efermi=0._dp
  real(dp) :: emin=-10._dp
  real(dp) :: emax=10._dp
  real(dp) :: rcut_nn=100._dp
  real(dp) :: rcut_tbg_nni=100._dp
  real(dp) :: rcut_grid=100._dp
  real(dp) :: avec(NDIM,NDIM)
  real(dp) :: bvec(NDIM,NDIM)
  real(dp) :: e_chi_exclude(2)
  real(dp) :: dis_frozen(2)
  type(CLproj) :: proj
  type(CLproj) :: base
  integer, allocatable :: nat_per_spec(:)
  integer, allocatable :: np_per_vert(:)
  integer, allocatable :: tot_iais(:,:)
  integer, allocatable :: iais_tot(:,:)
  real(dp), allocatable :: egrid(:)
  real(dp), allocatable :: atml(:,:,:)
  real(dp), allocatable :: vert(:,:)
  real(dp), allocatable :: kshift(:,:)
  contains
  procedure :: init=>read_input
  procedure :: atmc=>calc_atmc
  procedure :: write_geometry
  procedure :: shift_all
endtype CLpars

contains

subroutine read_input(THIS)
use geometry_library
class(CLpars), intent(inout) :: THIS
! internal
type(geomlib) geometry
integer iostat,iline,jline,ii
integer ispec,iat,ivert,igrid
integer ic,iw,jw
integer, parameter :: nmaxatm_pspec=40000
real(dp) t1,tvec(NDIM,NDIM),vd(NDIM)
character(len=256) block,arg,line
character(len=10) symbol
real(dp), allocatable :: atml_temp(:,:,:)
real(dp), allocatable :: xaxis(:,:)
real(dp), allocatable :: zaxis(:,:)
integer, allocatable :: lmr(:,:)

#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
THIS%Ggrid=0
THIS%ngrid=1
THIS%dis_frozen(1)=-1._dp
THIS%dis_frozen(2)= 1._dp

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

  call split_string(line,block,arg,' ')

  ! tasks block
  if (trim(block).eq."tasks") then
     read(arg,*,iostat=iostat) THIS%ntasks
     if (iostat.ne.0) call throw("paramters%read_input()","problem with ntasks argument")
     do ii=1,THIS%ntasks
        jline=jline+1
        read(50,'(A)',iostat=iostat) line
        if (iostat.ne.0) call throw("paramters%read_input()","problem with tasks block")
        if (mp_mpi) write(*,'(i6,": ",A)') jline,trim(adjustl(line))
        call split_string(line,block,arg,' ')
        THIS%tasks(ii)=trim(adjustl(block))
     end do

  ! real(read) and reciprocal (computed) lattice vectors
  else if (trim(block).eq."avec") then
    avec_block_found=.true.
    read(arg,*,iostat=iostat) t1
    if (iostat.ne.0) call throw("paramters%read_input()","problem with avec's scale argument")
    do ii=1,NDIM
      jline=jline+1
      read(50,*,iostat=iostat) THIS%avec(ii,:)
      if (iostat.ne.0) call throw("paramters%read_input()","problem with avec block")
      if(mp_mpi) write(*,'(i6,": ",5F10.6)') jline,THIS%avec(ii,:)
      THIS%avec(ii,:)=THIS%avec(ii,:)*t1
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

  else if (trim(block).eq."qgrid") then
    jline=jline+1
    read(50,*,iostat=iostat) THIS%qgrid(:)
    if (iostat.ne.0) call throw("paramters%read_input()","problem with qgrid data")
    if (mp_mpi) write(*,'(i6,": ",5I6)') jline,THIS%qgrid(:)

  else if (trim(block).eq."Ggrid") then
    jline=jline+1
    read(50,*,iostat=iostat) THIS%Ggrid(:)
    if (iostat.ne.0) call throw("paramters%read_input()","problem with Ggrid data")
    if (mp_mpi) write(*,'(i6,": ",5I6)') jline,THIS%Ggrid(:)

  else if (trim(block).eq."chi_exclude") then
    THIS%chi_exclude=.true.
    jline=jline+1
    read(50,*,iostat=iostat) THIS%e_chi_exclude(:)
    if (iostat.ne.0) call throw("paramters%read_input()","problem with chi_exclude data")
    if (mp_mpi) write(*,'(i6,": ",5F19.6)') jline,THIS%e_chi_exclude(:)

  else if (trim(block).eq."dis_frozen") then
    jline=jline+1
    read(50,*,iostat=iostat) THIS%dis_frozen(:)
    if (iostat.ne.0) call throw("paramters%read_input()","problem with dis_frozen data")
    if (mp_mpi) write(*,'(i6,": ",5F19.6)') jline,THIS%dis_frozen(:)

  ! BZ k-papth block
  else if (trim(block).eq."path") then
    read(arg,*,iostat=iostat) THIS%nvert
    if (iostat.ne.0) call throw("paramters%read_input()","problem with path's nvert argument")
    allocate(THIS%np_per_vert(THIS%nvert))
    allocate(THIS%vert(NDIM,THIS%nvert))
    do ivert=1,THIS%nvert
      jline=jline+1
      read(50,*,iostat=iostat) THIS%vert(:,ivert),THIS%np_per_vert(ivert)
      if (iostat.ne.0) call throw("paramters%read_input()","problem with path data")
      if (mp_mpi) write(*,'(i6,": ",5F10.6)',advance='no') jline,THIS%vert(:,ivert)
      if (mp_mpi) write(*,'(I6)') THIS%np_per_vert(ivert)
    end do

  else if (trim(block).eq."kshift") then
    kshift_block_found = .true.
    read(arg,*,iostat=iostat) THIS%nkshift
    if (iostat.ne.0) call throw("paramters%read_input()","problem with shift's nshift argument")
    allocate(THIS%kshift(NDIM,THIS%nkshift))
    do ivert=1,THIS%nkshift
      jline=jline+1
      read(50,*,iostat=iostat) THIS%kshift(:,ivert)
      if (iostat.ne.0) call throw("paramters%read_input()","problem with shift data")
      if (mp_mpi) write(*,'(i6,": ",5F10.6)') jline,THIS%kshift(:,ivert)
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
  ! cut off for in-plane nearest neighbors for tbg
  else if (trim(block).eq."rcut_tbg_nni") then
    jline=jline+1
    read(50,*,iostat=iostat) THIS%rcut_tbg_nni
    if (iostat.ne.0) call throw("paramters%read_input()","problem with rcut_tbg_nni data")
    if (mp_mpi) write(*,'(i6,": ",F10.6)') jline,THIS%rcut_tbg_nni
  ! cut off for spherical part of grid class
  else if (trim(block).eq."rcut_grid") then
    jline=jline+1
    read(50,*,iostat=iostat) THIS%rcut_grid
    if (iostat.ne.0) call throw("paramters%read_input()","problem with rcut_grid data")
    if (mp_mpi) write(*,'(i6,": ",F10.6)') jline,THIS%rcut_grid

  ! symmetry analyser mode
  else if (trim(block).eq."symtype") then
    jline=jline+1
    read(50,*,iostat=iostat) THIS%symtype
    if (iostat.ne.0) call throw("paramters%read_input()","problem with symtype data")
    if (mp_mpi) write(*,'(i6,": ",i4)') jline,THIS%symtype

  ! index of the first flat band of TBG
  else if (trim(block).eq."iflat_band") then
    jline=jline+1
    read(50,*,iostat=iostat) THIS%iflat_band
    if (iostat.ne.0) call throw("paramters%read_input()","problem with symtype data")
    if (mp_mpi) write(*,'(i6,": ",i6)') jline,THIS%iflat_band

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

  ! dimensions of the Coulomb potential
  else if (trim(block).eq."ndim_coul") then
    jline=jline+1
    read(50,*,iostat=iostat) THIS%ndim_coul
    NDIM_COUL=THIS%ndim_coul
    if (iostat.ne.0) call throw("paramters%read_input()","problem with ndim_coul data")
    if (mp_mpi) write(*,'(i6,": ",i6)') jline,THIS%ndim_coul

  ! criterium for neglection of Hamiltonian matrix elements
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

  ! .true. Only conventinal "diagonal" indexing of Hubbard U is assumed
  else if (trim(block).eq."HubU_diagonal") then
    jline=jline+1
    read(50,*,iostat=iostat) THIS%HubU_diagonal
    if (iostat.ne.0) call throw("paramters%read_input()","problem with HubU_diagonal data")
    if (mp_mpi) write(*,'(i6,": ",L6)') jline,THIS%HubU_diagonal

  ! .true. use weights in the amn matrix construction
  else if (trim(block).eq."use_weights_amn") then
    jline=jline+1
    read(50,*,iostat=iostat) THIS%use_weights_amn
    if (iostat.ne.0) call throw("paramters%read_input()","problem with use_weights_amn data")
    if (mp_mpi) write(*,'(i6,": ",L6)') jline,THIS%use_weights_amn

  ! .true. to read symmetries from SYMCRYS.OUT file 
  else if (trim(block).eq."readsym") then
    jline=jline+1
    read(50,*,iostat=iostat) THIS%readsym
    if (iostat.ne.0) call throw("paramters%read_input()","problem with readsym data")
    if (mp_mpi) write(*,'(i6,": ",L6)') jline,THIS%readsym

  ! .true. to write tight binding hamiltonian
  else if (trim(block).eq."writetb") then
    jline=jline+1
    read(50,*,iostat=iostat) THIS%writetb
    if (iostat.ne.0) call throw("paramters%read_input()","problem with writetb data")
    if (mp_mpi) write(*,'(i6,": ",L6)') jline,THIS%writetb

  ! .true. local overlap with exponent will be ignored in chi calculation
  else if (trim(block).eq."ignore_chiIq") then
    jline=jline+1
    read(50,*,iostat=iostat) THIS%ignore_chiIq
    if (iostat.ne.0) call throw("paramters%read_input()","problem with ignore_chiIq data")
    if (mp_mpi) write(*,'(i6,": ",L6)') jline,THIS%ignore_chiIq

  ! if .true. lattice is not shifted prior to symmetry analysis
  else if (trim(block).eq."symtshift") then
    jline=jline+1
    read(50,*,iostat=iostat) THIS%symtshift
    if (iostat.ne.0) call throw("paramters%read_input()","problem with symtshift data")
    if (mp_mpi) write(*,'(i6,": ",L6)') jline,THIS%symtshift

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
    THIS%chi_stop=THIS%nstates

  else if (trim(block).eq."chi_states") then
    jline=jline+1
    read(50,*,iostat=iostat) THIS%chi_start,THIS%chi_stop
    if (iostat.ne.0) call throw("paramters%read_input()","problem with chi states data")
    if (mp_mpi) write(*,'(i6,": ",2I6)') jline,THIS%chi_start,THIS%chi_stop
    if (THIS%chi_stop.lt.THIS%chi_start) call throw("paramters%read_input()", &
    "chi start state should be less or equal then chi stop")

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
        if (mp_mpi) write(*,'(i6,": ",A)') jline,trim(adjustl(THIS%sktype))

  ! projection mode for wannier export
  else if (trim(block).eq."wannier_proj_mode") then
        jline=jline+1
        read(50,'(A)',iostat=iostat) THIS%wannier_proj_mode
        if (iostat.ne.0) call throw("paramters%read_input()","problem with wannier_proj_mode data")
        if (mp_mpi) write(*,'(i6,": ",A)') jline,trim(adjustl(THIS%wannier_proj_mode))

  ! file with tight-binding hamiltonian
  else if (trim(block).eq."tbtype") then
        jline=jline+1
        read(arg,*,iostat=iostat) THIS%tbtype
        if (iostat.ne.0) call throw("paramters%read_input()","problem with tbfile argumet")
        if (trim(adjustl(THIS%tbtype)).eq.'tbfile'.or.trim(adjustl(THIS%tbtype)).eq.'hrfile') then
          read(50,*,iostat=iostat) THIS%tbfile
        end if
        if (mp_mpi) write(*,'(i6,": ",A)') jline,THIS%tbfile

  ! file with coulomb interaction in the real space
  else if (trim(block).eq."coulrs_file") then
        jline=jline+1
        read(50,*,iostat=iostat) THIS%coulrs_file
        if (iostat.ne.0) call throw("paramters%read_input()","problem with coulrs_file data")
        if (mp_mpi) write(*,'(i6,": ",A)') jline,THIS%coulrs_file

  ! character table file of symmetry representations
  else if (trim(block).eq."character_file") then
        jline=jline+1
        read(50,*,iostat=iostat) THIS%character_file
        if (iostat.ne.0) call throw("paramters%read_input()","problem with haracter_file data")
        if (mp_mpi) write(*,'(i6,": ",A)') jline,trim(adjustl(THIS%character_file))

  ! Numver of iteration in HubbardU symmetrization (in principle, 1 should be enough)
  else if (trim(block).eq."niter_symmetrize") then
    jline=jline+1
    read(50,*,iostat=iostat) THIS%niter_symmetrize
    if (iostat.ne.0) call throw("paramters%read_input()","problem with niter_symmetrize data")
    if (mp_mpi) write(*,'(i6,": ",i6)') jline,THIS%niter_symmetrize

  ! projection mode for wannier export
  else if (trim(block).eq."projections") then
        THIS%proj%allocatd=.true.
        ! temporary
        allocate(lmr(2,45000))
        allocate(xaxis(NDIM,45000))
        allocate(zaxis(NDIM,45000))
        read(arg,*,iostat=iostat) THIS%proj%ncenters
        if (iostat.ne.0) call throw("paramters%read_input()","problem with projections argumet")
        allocate(THIS%proj%norb_ic(THIS%proj%ncenters))
        allocate(THIS%proj%centers(NDIM,THIS%proj%ncenters))
        THIS%proj%norb=0
        do ic=1,THIS%proj%ncenters
           jline=jline+1
           read(50,*,iostat=iostat) THIS%proj%norb_ic(ic),THIS%proj%centers(:,ic)
           if (iostat.ne.0) call throw("paramters%read_input()","problem with wannier_proj_mode data")
           if (mp_mpi) write(*,'(i6,": ",i6,10F10.6)') jline,THIS%proj%norb_ic(ic),THIS%proj%centers(:,ic)
           do iw=1,THIS%proj%norb_ic(ic)
             jline=jline+1
             THIS%proj%norb=THIS%proj%norb+1
             read(50,*,iostat=iostat) symbol,xaxis(:,THIS%proj%norb),zaxis(:,THIS%proj%norb)
             if (iostat.ne.0) call throw("paramters%read_input()","problem with wannier_proj_mode data")
             if (mp_mpi) write(*,'(i6,": ",A,10F10.6)') jline, symbol,xaxis(:,THIS%proj%norb),zaxis(:,THIS%proj%norb)
             lmr(:,THIS%proj%norb)=string_to_lmr(symbol)
           end do
        end do
        allocate(THIS%proj%lmr(2,THIS%proj%norb))
        allocate(THIS%proj%waxis(NDIM,2,THIS%proj%norb))
        allocate(THIS%proj%iw2ic(THIS%proj%norb))
        THIS%proj%lmr=lmr
        do iw=1,THIS%proj%norb
          THIS%proj%waxis(:,1,iw)=xaxis(:,iw)
          THIS%proj%waxis(:,2,iw)=zaxis(:,iw)
        end do
        jw=0
        do ic=1,THIS%proj%ncenters
           do iw=1,THIS%proj%norb_ic(ic)
             jw=jw+1
             THIS%proj%iw2ic(jw)=ic
           end do
        end do
        deallocate(lmr,xaxis,zaxis)

  ! the same block is projection, but now it defines a base for TB hamiltonian
  else if (trim(block).eq."basis") then
        THIS%base%allocatd=.true.
        ! temporary
        allocate(lmr(2,45000))
        allocate(xaxis(NDIM,45000))
        allocate(zaxis(NDIM,45000))
        read(arg,*,iostat=iostat) THIS%base%ncenters
        if (iostat.ne.0) call throw("paramters%read_input()","problem with projections argumet")
        allocate(THIS%base%norb_ic(THIS%base%ncenters))
        allocate(THIS%base%centers(NDIM,THIS%base%ncenters))
        THIS%base%norb=0
        do ic=1,THIS%base%ncenters
           jline=jline+1
           read(50,*,iostat=iostat) THIS%base%norb_ic(ic),THIS%base%centers(:,ic)
           if (iostat.ne.0) call throw("paramters%read_input()","problem with basis centers data")
           if (mp_mpi) write(*,'(i6,": ",i6,10F10.6)') jline,THIS%base%norb_ic(ic),THIS%base%centers(:,ic)
           do iw=1,THIS%base%norb_ic(ic)
             jline=jline+1
             THIS%base%norb=THIS%base%norb+1
             read(50,*,iostat=iostat) symbol,xaxis(:,THIS%base%norb),zaxis(:,THIS%base%norb)
             if (iostat.ne.0) call throw("paramters%read_input()","problem with basis axis data")
             if (mp_mpi) write(*,'(i6,": ",A,10F10.6)') jline, symbol,xaxis(:,THIS%base%norb),zaxis(:,THIS%base%norb)
             lmr(:,THIS%base%norb)=string_to_lmr(symbol)
           end do
        end do
        allocate(THIS%base%lmr(2,THIS%base%norb))
        allocate(THIS%base%waxis(NDIM,2,THIS%base%norb))
        allocate(THIS%base%iw2ic(THIS%base%norb))
        THIS%base%lmr=lmr
        do iw=1,THIS%base%norb
          THIS%base%waxis(:,1,iw)=xaxis(:,iw)
          THIS%base%waxis(:,2,iw)=zaxis(:,iw)
        end do
        jw=0
        do ic=1,THIS%base%ncenters
           do iw=1,THIS%base%norb_ic(ic)
             jw=jw+1
             THIS%base%iw2ic(jw)=ic
           end do
        end do
        deallocate(lmr,xaxis,zaxis)


  end if

end do
100 continue
close(50)
if (mp_mpi) write(*,*)
#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif

if (.not.kshift_block_found) then
  THIS%nkshift=1
  allocate(THIS%kshift(NDIM,THIS%nkshift))
  THIS%kshift=0._dp
end if

! initialise basis for TB calculation
if (trim(adjustl(THIS%geometry_source)).ne."") then
   if (avec_block_found.or.atoms_block_found) then
     call throw("parameters%read_input",&
         "geometry library is not compatible with atoms or/and avec block")
   end if
   if (THIS%base%allocatd) then
     call throw("parameters%read_input",&
         "base is already allocataed, either remove 'basis' block or use only geometry_library")
   end if
   if (trim(adjustl(THIS%geometry_source)).ne."tbg".and.&
      trim(adjustl(THIS%geometry_source)).ne."slg") then
      call throw("paramters%read_input()","unknown geometry structure option")
   end if
   if (trim(adjustl(THIS%geometry_source)).eq.'tbg') THIS%symtshift=.false.

   ! exclude translation vectors from symmetry analyser it will take too long to compute that
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
   THIS%natmtot=sum(THIS%nat_per_spec)
   ! BASIS
   THIS%base%ncenters=THIS%natmtot
   THIS%base%norb=THIS%base%ncenters
   allocate(THIS%base%norb_ic(THIS%base%ncenters))
   allocate(THIS%base%centers(NDIM,THIS%base%ncenters))
   allocate(THIS%base%lmr(2,THIS%base%norb))
   allocate(THIS%base%waxis(NDIM,2,THIS%base%norb))
   allocate(THIS%base%iw2ic(THIS%base%norb))
   THIS%base%norb_ic=1
   THIS%base%lmr=1
   do iw=1,THIS%base%norb
     THIS%base%waxis(:,1,iw)=(/1._dp,0._dp,0._dp/)
     THIS%base%waxis(:,2,iw)=(/0._dp,0._dp,1._dp/)
   end do
   ic=0
   do ispec=1,THIS%nspec
     do iat=1,THIS%nat_per_spec(ispec)
       ic=ic+1
       THIS%base%centers(:,ic)=THIS%atml(:,iat,ispec)
     end do
   end do
   jw=0
   do ic=1,THIS%base%ncenters
      do iw=1,THIS%base%norb_ic(ic)
        jw=jw+1
        THIS%base%iw2ic(jw)=ic
      end do
   end do
   THIS%base%allocatd=.true.
else
   if (.not.avec_block_found) then
     call throw("parameters%read_input",&
         "if geometry library is not used, introduce the geometry with avec and atoms block")
   end if
   if (.not.atoms_block_found) then
     call throw("parameters%read_input",&
         "if geometry library is not used, introduce the geometry with avec and atoms block")
   end if
   if (.not.THIS%base%allocatd) then
     call throw("parameters%read_input",&
         "if geometry library is not used, introduce the basis via the 'basis' input block")
   end if
end if
! compute volumes
if (NDIM.eq.3) then
  call r3cross(THIS%avec(:,1),THIS%avec(:,2),vd)
  THIS%omega=abs(dot_product(vd,THIS%avec(:,3)))
  call r3cross(THIS%bvec(:,1),THIS%bvec(:,2),vd)
  THIS%omegabz=abs(dot_product(vd,THIS%bvec(:,3)))
else
  call warning("parameters%read_input","volume calculation is not possible with this NDIM, check that next code dont use it")
end if

! compute total number of atoms, and construct mapping to/from total index
THIS%natmtot=0
do ispec=1,THIS%nspec
  do iat=1,THIS%nat_per_spec(ispec)
    THIS%natmtot=THIS%natmtot+1
  end do
end do
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

if (THIS%ndim_coul.lt.NDIM) then
  if (mp_mpi) call info("paramters%read_input()","NDIM not equal to NDIM_COUL, &
    &be sure that your principal dimensions are along the first NDIM_COUL vectors")
else if (THIS%ndim_coul.gt.NDIM) then
  call throw("paramters%read_input()","NDIM_COUL can not be larger than NDIM")
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

!subroutine shift_centers(THIS,shift)
!class(CLproj), intent(inout) :: THIS
!real(dp), intent(in) :: shift(NDIM)
!integer ic
!if (THIS%shifted) call throw("CLpars%CLproj%shift_centers","shift can be applied only onse")
!THIS%shifted=.true.
!do ic=1,THIS%ncenters
!  THIS%centers(:,ic)=THIS%centers(:,ic)+shift(:)
!  THIS%centers(:,ic)=THIS%centers(:,ic)-nint(THIS%centers(:,ic))
!end do
!end subroutine

subroutine shift_all(THIS,shift)
class(CLpars), intent(inout) :: THIS
real(dp), intent(in) :: shift(NDIM,THIS%nmaxatm_pspec,THIS%nspec)
real(dp) dv(NDIM)
integer ispec,iat,ic
if (THIS%shifted) call throw("CLpars%shift_atml","shift can be applied only once")
THIS%shifted=.true.
do ispec=1,THIS%nspec
  do iat=1,THIS%nat_per_spec(ispec)
    if (THIS%proj%allocatd) then
      do ic=1,THIS%proj%ncenters
        dv=THIS%proj%centers(:,ic)-THIS%atml(:,iat,ispec)
        if (sum(abs(dv)).lt.epslat) then
          THIS%proj%centers(:,ic)=THIS%proj%centers(:,ic)+shift(:,iat,ispec)
          go to 9
        end if
      end do
      call throw("CLpars%shift_all","wannier porjection center not found or it is originally not at any atomic postion")
      9 continue
    end if
    if (THIS%base%allocatd) then
      do ic=1,THIS%base%ncenters
        dv=THIS%base%centers(:,ic)-THIS%atml(:,iat,ispec)
        if (sum(abs(dv)).lt.epslat) then
          THIS%base%centers(:,ic)=THIS%base%centers(:,ic)+shift(:,iat,ispec)
          go to 10
        end if
      end do
      call throw("CLpars%shift_all","center of basis function is not found or it is originally not at any atomic postion")
      10 continue
    end if
  end do
end do
THIS%atml=THIS%atml+shift
end subroutine

subroutine write_geometry(THIS,fname)
class(CLpars), intent(inout) :: THIS
character(len=*) fname
integer ii,ispec,iat
! write atoms geometry to a file
if (mp_mpi) then
  open(100,file=trim(adjustl(fname)),action="write")
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
end if
end subroutine

subroutine write_base(THIS,fname)
class(CLproj), intent(inout) :: THIS
character(len=*), intent(in) :: fname
integer ic,io,iw
! write atoms geometry to a file
if (mp_mpi) then
  iw=0
  open(100,file=trim(adjustl(fname)),action="write")
  do ic=1,THIS%ncenters
    write(100,'(I4,3F10.6)') THIS%norb_ic(ic), THIS%centers(:,ic)
    do io=1,THIS%norb_ic(ic)
      iw=iw+1
      write(100,'(2I3,x,6F10.6)')THIS%lmr(:,iw),THIS%waxis(:,1,iw),THIS%waxis(:,2,iw)
    end do
  end do
  close(100)
end if
end subroutine

end module
