
module symmetryclass
use modcom
use parameters
use symm_base
implicit none
private
logical, save :: shifted
integer, parameter :: nsmax=48
type, public :: CLsym
  logical :: defined=.false.
  integer :: nsym
  integer :: typ=1
  integer :: slat(NDIM,NDIM,nsmax)
  real(dp) :: lat(NDIM,NDIM,nsmax)
  real(dp) :: ltr(NDIM,NDIM,nsmax)
  real(dp) :: car(NDIM,NDIM,nsmax)
  real(dp) :: ctr(NDIM,NDIM,nsmax)
  integer :: inv(nsmax)
  integer :: lspl(nsmax)
  integer :: lspn(nsmax)
  real(dp) :: vtl(NDIM,nsmax)
  real(dp) :: vtc(NDIM,nsmax)
  integer, allocatable :: ieqat(:,:,:)
  contains
  procedure :: init=>find_symmetries
  
endtype CLsym

type, public, extends(CLsym) :: char_table
  character(len=10) :: symgroup
  character(len=4) :: classes(48)
  character(len=4) :: symreps(48)
  integer :: nclasses
  integer :: nsym_of_cl(48)
  integer :: isym_of_cl(48,48)
  real(dp) :: table(48,48)
  contains
  procedure :: read_character_table
endtype char_table

contains

subroutine find_symmetries(THIS,pars)
class(CLsym), intent(out) :: THIS
class(CLpars), intent(inout) :: pars
! local
integer isym
real(dp) id_lat_car(NDIM,NDIM),id_car_lat(NDIM,NDIM),tvec(NDIM,NDIM)
real(dp), allocatable :: m_loc(:,:)
real(dp), allocatable :: tau(:,:)
integer, allocatable :: ityp(:)
integer itot,iat,jat,ispec
real(dp) v1(NDIM),v2(NDIM),dmat(NDIM,NDIM)
real(dp), allocatable :: lattice_shift(:,:,:)

THIS%typ=pars%symtype
if (NDIM.ne.3) call throw("symmetryclass%find_symmetries()","this subroutine assumes NDIM=3")
allocate(lattice_shift(NDIM,pars%nmaxatm_pspec,pars%nspec))

if (pars%readsym) then
  call readsym(THIS%nsym,THIS%lat,THIS%vtl,THIS%inv)
else
  call symmetry(THIS%typ, pars%nspec, pars%nat_per_spec,&
                pars%nmaxatm_pspec, pars%natmtot, pars%avec, pars%atml,&
                THIS%nsym, THIS%lspl, THIS%lspn, THIS%inv, THIS%slat,&
                THIS%vtl, lattice_shift, shifted, pars%symtshift, mp_mpi)
end if



if (shifted) then
  call info("symmetryclass%find_symmetries()","lattice was shifted by symmetry analyser")  
  call info("","basis, projections were shifted by the same amount")  
  call info("","to make code working, be sure that your final basis set (base.dat)")  
  call info("","and wannier projections (proj.dat) are at positions compatible with atoms in 'geometry.dat'")  
  call pars%write_geometry('geometry0.dat')
  call pars%shift_all(lattice_shift)
  if (pars%base%allocatd) call pars%base%write_base('base.dat')
  if (pars%proj%allocatd) call pars%proj%write_base('proj.dat')
end if

call pars%write_geometry('geometry.dat')

THIS%defined=.true.

! identity matrix lattice->cartesian, linear algebra definition
id_lat_car=transpose(pars%avec)

! temporary
tvec=id_lat_car
! inverse of id_lat_cart, is id_cart_lat
call dmatrix_inverse(tvec,id_car_lat,NDIM)
do isym=1,THIS%nsym
  if (.not.pars%readsym) this%lat(:,:,isym)=dble(this%slat(:,:,this%lspl(isym)))
  this%car(:,:,isym)=matmul( id_lat_car,matmul(this%lat(:,:,isym),id_car_lat) )
  this%ltr(:,:,isym)=transpose(this%lat(:,:,isym))
  this%ctr(:,:,isym)=transpose(this%car(:,:,isym))
end do
do isym=1,THIS%nsym
  !this%vtl(:,isym)=matmul(this%lat(:,:,THIS%inv(isym)),this%vtl(:,isym))
  this%vtc(:,isym)=matmul(id_lat_car,this%vtl(:,isym))
end do

!!!! QE calls !!!!
if(.false.) then
at=transpose(pars%avec)
bg=transpose(pars%bvec)/twopi
allocate(ityp(pars%natmtot),m_loc(NDIM,pars%natmtot))
allocate(tau(NDIM,pars%natmtot))
do itot=1,pars%natmtot
  iat=pars%tot_iais(itot,1)
  ispec=pars%tot_iais(itot,2)
  m_loc(:,iat)=0._dp
  ityp(iat)=ispec
  tau(:,iat)=pars%atmc(iat,ispec)
end do
call set_sym( pars%natmtot, tau, ityp, 1, m_loc )
THIS%nsym=nsym
do isym=1,nsym
  THIS%vtl(:,isym)=ft(:,isym)
  THIS%lat(:,:,isym)=transpose(s(:,:,isym))
  THIS%car(:,:,isym)=transpose(sr(:,:,isym))
  THIS%vtc(:,isym)=matmul(ft(:,isym),pars%avec)
  THIS%ltr(:,:,isym)=sr(:,:,isym)
  THIS%ctr(:,:,isym)=s(:,:,isym)
end do
THIS%inv(1:THIS%nsym)=invs(1:nsym)
if(mp_mpi)then
   open(unit=105, file='SYMCRYS.OUT',form='formatted')
   write(105,"(i5)") nsym
   do isym=1,nsym
      write(105,*)
      write(105,"(1p,3e23.15)") THIS%car(:,:,isym), THIS%vtc(:,isym)/sum(pars%avec(1,:)**2)
   end do
   close(105)
end if
do isym=1,nsym
   if(THIS%inv(isym).le.0.or.THIS%inv(isym).ge.THIS%nsym+1) then
      call throw("CLsym%find_symmetries", "out of range in invs from QE")
   end if
   dmat=0d0
   dmat(1,1)=1d0
   dmat(2,2)=1d0
   dmat(3,3)=1d0
   v1=matmul(matmul(THIS%vtc(:,isym),THIS%car(:,:,THIS%inv(isym)))+THIS%vtc(:,THIS%inv(isym)),bg)
   if(sum(abs(matmul(THIS%car(:,:,isym),THIS%car(:,:,THIS%inv(isym)))-dmat))+sum(abs(v1-dble(nint(v1)))).gt.1._dp) then
      call throw("CLsym%find_symmetries", "inconsistent invs from QE library")
   end if
end do
deallocate(tau,m_loc,ityp)
end if
!!!! QE calls (END) !!!!


allocate(THIS%ieqat(pars%nmaxatm_pspec,pars%nspec,THIS%nsym))
do isym=1,THIS%nsym
   do ispec=1,pars%nspec
     do iat=1,pars%nat_per_spec(ispec)
       v1=matmul(THIS%lat(:,:,isym),pars%atml(:,iat,ispec)+THIS%vtl(:,isym))
       do jat=1,pars%nat_per_spec(ispec)
         v2=v1-pars%atml(:,jat,ispec)
         if(sum(abs(dble(nint(v2))-v2)).lt.1d-2) then
            THIS%ieqat(iat,ispec,isym)=jat
            go to 9
         end if
       end do
       write(*,*) "Error(symmetry): atom not found"
       write(*,*) THIS%vtl(:,isym)
       write(*,*) pars%atml(:,iat,ispec)
       write(*,*) THIS%lat(1,:,isym)
       write(*,*) THIS%lat(2,:,isym)
       write(*,*) THIS%lat(3,:,isym)
       stop
       9 continue
     end do
   end do
end do
if(mp_mpi)then
   open(unit=105, file='wannier.sym',form='formatted')
   write(105,"(i5)") THIS%nsym
   do isym=1,THIS%nsym
      write(105,*)
      write(105,"(1p,3e23.15)") THIS%car(:,:,isym), THIS%vtc(:,isym)/sqrt(sum(pars%avec(1,:)**2))
   end do
   close(105)
end if

end subroutine

subroutine read_character_table(THIS,fname,nsym_check)
class(char_table), intent(out) :: THIS
character(len=*), intent(in) :: fname
integer, intent(in) :: nsym_check
integer nsym,icl,irep,j
character(len=256) line,line1,line2
logical exs
inquire(file=trim(adjustl(fname)),exist=exs)
if (.not.exs) then
  call throw("char_table%read_character_table","symmetry representation file '"//trim(adjustl(fname))//"' not present")
end if
open(unit=50,file=trim(adjustl(fname)),action='read',status='old')
read(50,*) THIS%symgroup
read(50,*) nsym
if (nsym.ne.nsym_check) call throw('char_table%read_character_table',&
 &'number of symmetry operations in charater table file not equal to actual nsym')
read(50,*) THIS%nclasses
THIS%nsym_of_cl(:)=0
do icl=1,THIS%nclasses
  read(50,'(A)') line
  call split_string(line,line1,line2,'!')
  line=line1
  do j=1,100
    call split_string(line,line1,line2,' ')
    line=line1
    if (line1.eq.'')   call throw('char_table%read_character_table','no symmetry operations in class of character table file')
    THIS%nsym_of_cl(icl)=THIS%nsym_of_cl(icl)+1
    read(line,*) THIS%isym_of_cl(icl,THIS%nsym_of_cl(icl))
    line=line2
    if (trim(adjustl(line2)).eq.'') exit
  end do
end do
read(50,*) THIS%classes(1:THIS%nclasses)
do irep=1,THIS%nclasses
  read(50,*) THIS%table(irep,1:THIS%nclasses),THIS%symreps(irep)
end do
close(50)
end subroutine

subroutine readsym(nsym,lat,vtl,inv)
integer, intent(out) :: nsym
real(dp), intent(out) :: lat(NDIM,NDIM,nsmax),vtl(NDIM,nsmax)
integer, intent(out) :: inv(nsmax)
logical exs
integer isym,jsym,jj
real(dp) uni(NDIM,NDIM),tmat(NDIM,NDIM)
inquire(file='SYMCRYS.OUT',exist=exs)
if (.not.exs) then
  call throw("CLsym%read_symmetries","SYMCRYS.OUT file not found")
end if
open(unit=50,file='SYMCRYS.OUT',action='read',status='old')
read(50,*)
read(50,*)
read(50,*)
read(50,*) nsym
do isym=1,nsym
  read(50,*)
  read(50,*)
  read(50,*)
  read(50,*) vtl(:,isym)
  ! spacial rotation
  read(50,*)
  read(50,*) lat(1,:,isym)
  read(50,*) lat(2,:,isym)
  read(50,*) lat(3,:,isym)
  ! spin rotation (not implemented)
  read(50,*)
  read(50,*)
  read(50,*)
  read(50,*)
end do
inv=-1
uni=0._dp
do jj=1,NDIM
  uni(jj,jj)=1._dp
end do
do isym=1,nsym
   do jsym=1,nsym
     tmat=matmul(lat(:,:,isym),lat(:,:,jsym))
     if (sum(abs(tmat-uni)).lt.epslat) inv(isym)=jsym
   end do
end do
close(50)
end subroutine

end module
