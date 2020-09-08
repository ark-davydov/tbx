
module gridclass
use modcom
use parameters
use symmetryclass
implicit none
private

type, public :: CLgrid
  character(len=4) mode
  integer npt
  real(dp) vecs(NDIM,NDIM)
  integer, allocatable, private :: vi_(:,:)
  real(dp), allocatable, private :: vl_(:,:)
  contains
  procedure :: vpi=>give_vpi
  procedure :: vpl=>give_vpl
  procedure :: vpc=>calc_vpc
  procedure :: dc=>calc_dc
endtype CLgrid

type, public, extends(CLgrid) :: GRID
  logical :: trevinit_done=.false.
  logical :: syminit_done=.false.
  logical :: sphere_allocated=.false.
  logical centered
  logical fractional
  integer ngrid(NDIM)
  integer :: ip0=-1
  integer nir
  integer nirT
  integer npt_sphere
  real(dp) :: vq0
  integer, allocatable :: sphere_to_homo(:)
  integer, allocatable :: homo_to_sphere(:)
  integer, allocatable :: iks2k(:,:)
  integer, allocatable :: sik2ir(:)
  integer, allocatable :: ik2ir(:)
  integer, allocatable :: ir2ik(:)
  ! time reversal
  integer, allocatable :: ikT2k(:)
  integer, allocatable :: ikT2ir(:)
  integer, allocatable :: irT2ik(:)
  contains 
  procedure :: init=>init_grid
  procedure :: init_sphere
  procedure :: io=>io_grid
  procedure :: find=>find_grid_point
  procedure :: sym_init
  procedure :: vpl_sphere
  procedure :: vpc_sphere
  procedure :: dc_sphere
endtype GRID

type, public, extends(CLgrid) :: PATH
  integer nvert
  integer, allocatable :: np_per_vert(:)
  real(dp), allocatable :: vert(:,:)
  real(dp), allocatable :: dvpc(:)
  real(dp), allocatable :: dvrt(:)
  contains
  procedure :: init=>init_path
endtype PATH

contains

subroutine init_sphere(THIS,pars)
class(GRID), intent(inout) :: THIS
class(CLpars), intent(in) :: pars
integer ip
if (THIS%mode.ne."grid") call throw("GRID%init_sphere","grid mode has to be 'grid'")
if (THIS%centered.eqv..false.) call throw("GRID%init_sphere","grid has to be centered")
if (THIS%fractional.eqv..true.) call throw("GRID%init_sphere","grid needs to be non-fractional")
THIS%sphere_allocated=.true.
! count the number of grid points in the sphere
THIS%npt_sphere=0
do ip=1,THIS%npt
  if (THIS%dc(ip).lt.pars%rcut_grid) then
    THIS%npt_sphere=THIS%npt_sphere+1
  end if
end do
! allocate map from sphere point to homogeneous point
allocate(THIS%sphere_to_homo(THIS%npt_sphere))
allocate(THIS%homo_to_sphere(THIS%npt))
THIS%npt_sphere=0
THIS%homo_to_sphere=0
do ip=1,THIS%npt
  if (THIS%dc(ip).lt.pars%rcut_grid) then
    THIS%npt_sphere=THIS%npt_sphere+1
    THIS%sphere_to_homo(THIS%npt_sphere)=ip
    THIS%homo_to_sphere(ip)=THIS%npt_sphere
  end if
end do
if (mp_mpi) then
  call info("GRID%init_sphere","spherical part of grid object extracted")
  write(*,'("  grid parameters: ",20I8)')  THIS%ngrid(:)
  write(*,'("  npt_sphere: ",I8)')  THIS%npt_sphere
end if


end subroutine

subroutine init_grid(THIS,ngrid,vecs,centered,fractional)
class(GRID), intent(out) :: THIS
integer, intent(in) :: ngrid(NDIM)
real(dp), intent(in) :: vecs(NDIM,NDIM)
logical, intent(in) :: centered,fractional
integer id,i1,i2,i3,ii,jj,ip
THIS%mode="grid"
THIS%centered=centered
THIS%fractional=fractional
THIS%vecs=vecs
THIS%ngrid=ngrid
THIS%npt=1
do id=1,NDIM
  THIS%npt=THIS%npt*THIS%ngrid(id)
end do
allocate(THIS%vi_(NDIM,THIS%npt))
allocate(THIS%vl_(NDIM,THIS%npt))
do ip=1,THIS%npt
  ! model for the sequence index in (0:N-1 notation) is:
  ! ii=N2*N3*i1+N3*i2+i3
  ! i.e. row-major indexing, therefor the map from ii to i1,i2,i3 is known and used below
  ii=ip-1
  if (NDIM.eq.1) then
    i1=ii
    if (centered) then
      THIS%vi_(1,ip)=i1-THIS%ngrid(1)/2
    else
      THIS%vi_(1,ip)=i1
    end if
  else if (NDIM.eq.2) then
    i2=mod(ii,THIS%ngrid(2))
    i1=int(ii/THIS%ngrid(2))
    if (centered) then
      THIS%vi_(1,ip)=i1-THIS%ngrid(1)/2
      THIS%vi_(2,ip)=i2-THIS%ngrid(2)/2
    else
      THIS%vi_(1,ip)=i1
      THIS%vi_(2,ip)=i2
    end if
  else if (NDIM.eq.3) then
    i3=mod(ii,THIS%ngrid(3))
    jj=int(ii/THIS%ngrid(3))
    i2=mod(jj,THIS%ngrid(2))
    i1=int(jj/THIS%ngrid(2))
    if (centered) then
      THIS%vi_(1,ip)=i1-THIS%ngrid(1)/2
      THIS%vi_(2,ip)=i2-THIS%ngrid(2)/2
      THIS%vi_(3,ip)=i3-THIS%ngrid(3)/2
    else
      THIS%vi_(1,ip)=i1
      THIS%vi_(2,ip)=i2
      THIS%vi_(3,ip)=i3
    end if
    !write(*,*) i1,i2,i3,"!",ii,"!",i1*THIS%ngrid(3)*THIS%ngrid(2)+i2*THIS%ngrid(3)+i3
  else
    call throw("GRID%init()","NDIM is greater than 3")
  end if
  if (fractional) then
    !write(*,'(I3," ! ",3I3," ! ",3I3)')ip,i1,i2,i3,i1-THIS%ngrid(1)/2,i2-THIS%ngrid(2)/2,i3-THIS%ngrid(3)/2
    THIS%vl_(:,ip)=dble(THIS%vi_(:,ip))/dble(THIS%ngrid(:))
  else
    THIS%vl_(:,ip)=dble(THIS%vi_(:,ip))
  end if
  if (sum(abs(THIS%vl_(:,ip))).lt.epslat) THIS%ip0=ip
end do
if (fractional) then
  ! fractinal grid is most probably a k/q -grid in reciprocal space,
  ! so define Coulomb interaction at q=0
  if (THIS%ip0.le.0) call throw("GRID%init()","not possible situtaion: ip0 not found")
  call gengclq(THIS%ngrid,THIS%vecs,THIS%vq0)
  if (mp_mpi) write(*,'("Info[GRID]: Average Coulomb interaction at q=0 for this ngrid(:): ",G18.10)') THIS%vq0
  if (.false.) then
    ! test the quality of approximation for vq0
    do ip=1,THIS%npt
      if (ip.eq.THIS%ip0) then
         write(155,*) THIS%dc(ip),THIS%vq0
      else
         if (NDIM_COUL.eq.2) then
           write(155,*) THIS%dc(ip),twopi*CoulombForceConstant/THIS%dc(ip)
         else
           write(155,*) THIS%dc(ip),fourpi*CoulombForceConstant/THIS%dc(ip)**2
         end if
      end if
    end do
    stop
  end if
end if
end subroutine

subroutine init_path(THIS,nvert,np_per_vert,vert,vecs)
class(PATH), intent(out) :: THIS
integer, intent(in) :: nvert,np_per_vert(nvert)
real(dp), intent(in) :: vert(NDIM,nvert)
real(dp), intent(in) :: vecs(NDIM,NDIM)
integer ip,ivert,counter
real(dp) v1(NDIM),v2(NDIM),v3(NDIM)
if (nvert.le.1) call throw("PATH%init_path()","nvert should be at least 2, now it is <=1")
THIS%mode="path"
allocate(THIS%vert(NDIM,nvert))
allocate(THIS%np_per_vert(nvert))
THIS%nvert=nvert
THIS%vecs=vecs
THIS%vert=vert
THIS%np_per_vert(:)=np_per_vert(:)
THIS%np_per_vert(1)=0
THIS%npt=sum(THIS%np_per_vert(:))
allocate(THIS%vl_(NDIM,THIS%npt))
allocate(THIS%dvpc(THIS%npt))
allocate(THIS%dvrt(THIS%nvert))
THIS%dvpc(:)=0._dp
THIS%dvrt(:)=0._dp
counter=0
do ivert=2,THIS%nvert
  if (THIS%np_per_vert(ivert).lt.2) then
    call throw("PATH%init_path()","np_per_vert should be greater than 1")
  end if
  v1=THIS%vert(:,ivert-1)
  v2=THIS%vert(:,ivert)
  v3=matmul(v2-v1,THIS%vecs)
  THIS%dvrt(ivert)=THIS%dvrt(ivert-1)+sqrt(dot_product(v3,v3))
  do ip=1,THIS%np_per_vert(ivert)
    counter=counter+1
    if (ivert.eq.THIS%nvert) then
      THIS%vl_(:,counter)=v1+dble(ip-1)/dble(THIS%np_per_vert(ivert)-1)*(v2-v1)
    else
      THIS%vl_(:,counter)=v1+dble(ip-1)/dble(THIS%np_per_vert(ivert))*(v2-v1)
    end if
    if (counter.gt.1) then
        v3=matmul(THIS%vl_(:,counter)-THIS%vl_(:,counter-1),THIS%vecs)
        THIS%dvpc(counter)=THIS%dvpc(counter-1)+sqrt(dot_product(v3,v3))
    end if
  end do
end do
end subroutine

subroutine sym_init(THIS,sym)
class(GRID), intent(inout) :: THIS
class(CLsym), intent(in) :: sym
! local
integer ik,ikp,isym
integer igk(NDIM+1)
logical lfound(THIS%npt)
real(dp) sv(NDIM)
real(dp) avec(NDIM,NDIM),tvec(NDIM,NDIM)
real(dp) v1(NDIM),v2(NDIM),v3(NDIM),v4(NDIM)
if (THIS%syminit_done) call throw("gridclass%sym_init","second call of sym_init")
THIS%syminit_done=.true.
allocate(THIS%iks2k(THIS%npt,sym%nsym))
tvec=THIS%vecs
call dmatrix_inverse(tvec,avec,NDIM)
THIS%iks2k=-999 !Sym.op.(isym) moves k(iks2k(ik,isym)) to k(ik) + G(iks2g(ik,isym)).
do isym=1,sym%nsym
   do ik=1,THIS%npt
     sv=matmul(sym%car(:,:,isym),THIS%vpc(ik))
     sv=matmul(sv,avec)
     igk=THIS%find(sv)
     THIS%iks2k(ik,isym)=igk(NDIM+1)
   end do
end do
if (.false.) then
! origican code from QE
do isym=1,sym%nsym
   lfound=.false.
   do ik=1,THIS%npt
      !v1=THIS%vpc(ik)
      v1=THIS%vpc(ik)
      !v2=matmul(transpose(sym%lat(:,:,isym)),v1)
      v2=matmul(sym%car(:,:,isym),v1)
      v2=matmul(v2,avec)
      do ikp=1,THIS%npt
         if(lfound(ikp)) cycle
         v3=THIS%vpl(ikp)
         !v3=THIS%vpc(ikp)
         !v4=matmul(v2-v3,THIS%vecs)
         v4=v2-v3
         write(*,*) v4
         if(sum(abs(nint(v4)-v4)).lt.1d-5) then
            THIS%iks2k(ik,isym)=ikp
            lfound(ikp)=.true.
         end if
         if(THIS%iks2k(ik,isym).ge.1) exit
      end do
   end do
end do
end if
allocate(THIS%ik2ir(THIS%npt))
allocate(THIS%ir2ik(THIS%npt))
allocate(THIS%sik2ir(THIS%npt))
THIS%ik2ir=-999 !Gives irreducible-k points from regular-k points.
THIS%ir2ik=-999 !Gives regular-k points from irreducible-k points.
lfound=.false.
THIS%nir=0
do ik=1,THIS%npt
   if(lfound(ik)) cycle
   lfound(ik)=.true.
   THIS%nir=THIS%nir+1
   THIS%ir2ik(THIS%nir)=ik
   THIS%ik2ir(ik)=THIS%nir
   THIS%sik2ir(ik)=1
   do isym=1,sym%nsym
      ikp=THIS%iks2k(ik,isym)
      if(lfound(ikp)) cycle
      lfound(ikp)=.true.
      THIS%ik2ir(ikp)=THIS%nir
      THIS%sik2ir(ikp)=sym%inv(isym)
   end do
end do
call info("GRID%sym_init()","symmetry is initialized on the grid")
call info("GRID%sym_init()","attemp to initilize k->-k for TR symmetry...")
allocate(THIS%ikT2k(THIS%npt))
THIS%ikT2k=-999 
do ik=1,THIS%npt
  igk=THIS%find(-THIS%vpl(ik))
  THIS%ikT2k(ik)=igk(NDIM+1)
end do
allocate(THIS%ikT2ir(THIS%npt))
allocate(THIS%irT2ik(THIS%npt))
THIS%ikT2ir=-999 
THIS%irT2ik=-999 
lfound=.false.
THIS%nirT=0
do ik=1,THIS%npt
   if(lfound(ik)) cycle
   lfound(ik)=.true.
   THIS%nirT=THIS%nirT+1
   THIS%irT2ik(THIS%nirT)=ik
   THIS%ikT2ir(ik)=THIS%nirT
   do isym=1,sym%nsym
      ikp=THIS%ikT2k(ik)
      if(lfound(ikp)) cycle
      lfound(ikp)=.true.
      THIS%ikT2ir(ikp)=THIS%nirT
   end do
end do
call info("GRID%sym_init()",".. TR symmetry is initialized on the grid")
end subroutine

subroutine io_grid(THIS,unt,fname,action,pars,norb)
class(GRID), intent(inout) :: THIS
integer, intent(in) :: unt
integer, intent(in) :: norb
integer norb_
class(CLpars), intent(in) :: pars
character(len=*), intent(in) :: action,fname
! local
integer nd,nstates
logical centered,fractional
integer ngrid(NDIM)
real(dp) vecs(NDIM,NDIM)
if (trim(adjustl(action)).ne."write".and.trim(adjustl(action)).ne."read") &
call throw("modcom%io_grid()","unknown read/write action")
open(unt,file=trim(adjustl(fname)),form='unformatted',action=trim(adjustl(action)))
! the following uniquely defines k-grid, also eigen-vectors grid is written here
if (trim(adjustl(action)).eq."write") then
  write(unt) NDIM
  write(unt) THIS%centered
  write(unt) THIS%fractional
  write(unt) THIS%ngrid(:)
  write(unt) THIS%vecs(:,:)
  write(unt) norb
  write(unt) pars%nstates
else
  read(unt) nd
  if (nd.ne.NDIM) call throw("GRID%io_grid()","different dimensionality in the _file")
  read(unt) centered
  read(unt) fractional
  read(unt) ngrid(:)
  if (sum(abs(ngrid-pars%ngrid)).ne.0) &
   call throw("grid%io_grid()","the k-point grid dimensions in _grid file are different from ones derived from input")
  read(unt) vecs(:,:)
  read(unt) norb_
  if (norb_.ne.norb) &
   call throw("grid%io_grid()","number of orbitals in _grid file are different from ones derived from input")
  read(unt) nstates
  if (pars%nstates.ne.nstates) &
   call throw("grid%io_grid()","the number of states in _grid file is different from the one derived from input")
  THIS%centered=centered
  THIS%fractional=fractional
  THIS%ngrid=ngrid
  THIS%vecs=vecs
  call THIS%init(ngrid,vecs,centered,fractional)
end if
close(unt)
end subroutine


! access to private vertex variables
function give_vpi(THIS,ip) result(vi)
class(CLgrid), intent(in) :: THIS
integer, intent(in) :: ip
integer :: vi(NDIM)
if (THIS%mode.eq."path") call throw("CLgrid%give_vpi()","vpi method is not available for path object")
vi(:)=THIS%vi_(:,ip)
end function
function give_vpl(THIS,ip) result(vl)
class(CLgrid), intent(in) :: THIS
integer, intent(in) :: ip
real(dp) :: vl(NDIM)
vl(:)=THIS%vl_(:,ip)
end function
function calc_vpc(THIS,ip) result(vc)
class(CLgrid), intent(in) :: THIS
integer, intent(in) :: ip
real(dp) :: vc(NDIM)
vc(:)=matmul(THIS%vl_(:,ip),THIS%vecs(:,:))
end function
function calc_dc(THIS,ip) result(dc)
class(CLgrid), intent(in) :: THIS
integer, intent(in) :: ip
real(dp) :: dc
real(dp) :: vc(NDIM)
vc(:)=matmul(THIS%vl_(:,ip),THIS%vecs(:,:))
dc=sqrt(dot_product(vc,vc))
end function
! quantinies of spherical grid
function vpl_sphere(THIS,ip) result(vl)
class(GRID), intent(in) :: THIS
integer, intent(in) :: ip
real(dp) :: vl(NDIM)
vl(:)=THIS%vl_(:,THIS%sphere_to_homo(ip))
end function
function vpc_sphere(THIS,ip) result(vc)
class(GRID), intent(in) :: THIS
integer, intent(in) :: ip
real(dp) :: vc(NDIM)
vc(:)=matmul(THIS%vl_(:,THIS%sphere_to_homo(ip)),THIS%vecs(:,:))
end function
function dc_sphere(THIS,ip) result(dc)
class(GRID), intent(in) :: THIS
integer, intent(in) :: ip
real(dp) :: dc
real(dp) :: vc(NDIM)
vc(:)=matmul(THIS%vl_(:,THIS%sphere_to_homo(ip)),THIS%vecs(:,:))
dc=sqrt(dot_product(vc,vc))
end function

function find_grid_point(THIS,vpl) result(ikg)
class(GRID), intent(in) :: THIS
real(dp), intent(in) :: vpl(NDIM)
integer ikg(NDIM+1)
real(dp) vchk(NDIM)
if (THIS%centered) then
  if (THIS%fractional) then
    ! centered grid is defined differently, this complicated construct is becase THIS%ngrid/2 
    ! action was used in the creation of the grid, and the result is different for odd and even grids
    ! so here we just plug exactly the same as was used in the initialisation
    ikg(1:NDIM)=nint( (vpl+dble(THIS%ngrid/2)/dble(THIS%ngrid))*dble(THIS%ngrid) )
    ikg(1:NDIM)=modulo(ikg(1:NDIM),THIS%ngrid)
  else
    ikg(1:NDIM)=nint(vpl+dble(THIS%ngrid/2))
  end if
else
  if (THIS%fractional) then
    ikg(1:NDIM)=modulo(nint(vpl*dble(THIS%ngrid)),THIS%ngrid)
  else
    ikg(1:NDIM)=nint(vpl)
  end if
end if
if (NDIM.eq.1) then
  ikg(NDIM+1)=ikg(NDIM)
else if (NDIM.eq.2) then
  ikg(NDIM+1)=THIS%ngrid(2)*ikg(1)+ikg(2)
else if (NDIM.eq.3) then
  ikg(NDIM+1)=THIS%ngrid(3)*THIS%ngrid(2)*ikg(1)+THIS%ngrid(3)*ikg(2)+ikg(3)+1
else
  call throw("GRID%find_grid_point()","NDIM is greater than 3")
end if
if (THIS%centered) then
  ! return back to orgiginal values
  ikg(1:NDIM)=ikg(1:NDIM)-THIS%ngrid/2
end if
if (ikg(NDIM+1).gt.THIS%npt.or.ikg(NDIM+1).le.0) then
  ikg(NDIM+1)=-abs(ikg(NDIM+1))
  ikg(1:NDIM)=-99999
else
  ! finally, redefine ikg(1:NDIM) such that they return translation vectors
  ikg(1:NDIM)=nint(vpl-THIS%vpl(ikg(NDIM+1)))
  ! do the check
  vchk=dble(ikg(1:3))+THIS%vpl(ikg(NDIM+1))
  if (sum(abs(vchk-vpl)).gt.epslat) then
    call throw("GRID%find_grid_point()","grid point not found")
  end if
end if

end function

end module
