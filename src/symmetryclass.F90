
module symmetryclass
use modcom
use parameters
implicit none
private
integer, parameter :: nsmax=48
type, public :: CLsym
  logical :: defined=.false.
  integer :: nsym
  integer :: typ=2
  integer :: lat(NDIM,NDIM,nsmax)
  real(dp) :: car(NDIM,NDIM,nsmax)
  integer :: inv(nsmax)
  integer :: lspl(nsmax)
  integer :: lspn(nsmax)
  real(dp) :: vtl(NDIM,nsmax)
  real(dp) :: vtc(NDIM,nsmax)
  integer, allocatable :: ieqat(:,:,:)
  contains
  procedure :: init=>find_symmetries
  
endtype CLsym

contains

subroutine find_symmetries(THIS,pars)
class(CLsym), intent(out) :: THIS
class(CLpars), intent(inout) :: pars
! local
integer isym
real(dp) id_lat_car(NDIM,NDIM),id_car_lat(NDIM,NDIM),tvec(NDIM,NDIM)
THIS%typ=pars%symtype
allocate(THIS%ieqat(pars%nmaxatm_pspec,pars%nspec,nsmax))
call symmetry(THIS%typ, pars%nspec, pars%nat_per_spec,&
              pars%nmaxatm_pspec, pars%natmtot, pars%avec, pars%atml,&
              THIS%nsym, THIS%lspl, THIS%lspn, THIS%inv, THIS%lat,&
              THIS%ieqat , THIS%vtl)
THIS%defined=.true.
! identity matrix lattice->cartesian, linear algebra definition
id_lat_car=transpose(pars%avec)
! temporary
tvec=id_lat_car
! inverse of id_lat_cart, is id_cart_lat
call dmatrix_inverse(tvec,id_car_lat,NDIM)
do isym=1,THIS%nsym
  THIS%car(:,:,isym)=matmul( id_lat_car,matmul(THIS%lat(:,:,isym),id_car_lat) )
  THIS%vtc(:,isym)=matmul(id_lat_car,THIS%vtl(:,isym))
end do
end subroutine


end module
