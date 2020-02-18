
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
  integer :: inv(nsmax)
  integer :: lspl(nsmax)
  integer :: lspn(nsmax)
  real(dp) :: vtl(NDIM,nsmax)
  integer, allocatable :: ieqat(:,:,:)
  integer, allocatable :: ieqat_inv(:,:,:)
  contains
  procedure :: init=>find_symmetries
endtype CLsym

contains

subroutine find_symmetries(THIS,pars)
class(CLsym), intent(out) :: THIS
class(CLpars), intent(inout) :: pars
THIS%typ=pars%symtype
allocate(THIS%ieqat(pars%nmaxatm_pspec,pars%nspec,nsmax))
allocate(THIS%ieqat_inv(pars%nmaxatm_pspec,pars%nspec,nsmax))
call symmetry(THIS%typ, pars%nspec, pars%nat_per_spec,&
              pars%nmaxatm_pspec, pars%natmtot, pars%avec, pars%atml,&
              THIS%nsym, THIS%lspl, THIS%lspn, THIS%inv, THIS%lat,&
              THIS%ieqat , THIS%ieqat_inv, THIS%vtl)
THIS%defined=.true.
end subroutine

end module
