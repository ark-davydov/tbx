
module wannier_interface
#ifdef MPI
  use mpi
#endif
use modcom
use parameters
use tbclass
use gridclass
implicit none
private

type, public :: CLwan
  character(len=100) :: proj_mode
  integer :: num_bands
  integer :: nwan
  real(dp), private :: avec(NDIM,NDIM)
  real(dp), private :: bvec(NDIM,NDIM)
  contains
  procedure :: init=>copy_parameters
  procedure :: project=>compute_projection
endtype CLwan

contains

subroutine copy_parameters(THIS,pars)
class(CLwan), intent(inout) :: THIS
class(CLpars), intent(in) :: pars
THIS%avec=pars%avec
THIS%bvec=pars%bvec
THIS%proj_mode=pars%wannier_proj_mode
THIS%num_bands=pars%nstates
THIS%nwan=pars%nwan
end subroutine

subroutine compute_projection(THIS,tbmodel,kgrid,eval,evec)
class(CLwan), intent(in) :: THIS
class(CLtb), intent(in) :: tbmodel
class(GRID), intent(in) :: kgrid
real(dp), intent(in) :: eval(THIS%num_bands,kgrid%npt)
complex(dp), intent(in) :: evec(tbmodel%norb_TB,THIS%num_bands,kgrid%npt)
complex(dp), allocatable :: wftrial(:,:)
allocate(wftrial(tbmodel%norb_TB,THIS%nwan))
call generate_trial_wavefunctions(THIS,tbmodel%norb_TB,THIS%num_bands,wftrial)
! A_mn(k)=<\psi_m(k)|g_n>
!open(50,file='structure.amn',action='write')
!write(50,*) ' structure.amn file '
!write(50,*) num_bands,nkpt,nwan
!do ik=1,nkpt
!  call geteigvec(ik,eigc)
!  do i=1,nwan
!    j=1
!    do jst=ist_wan,jst_wan
!      z1=dot_product(eigc(:,jst),wftrial(:,i))
!      write(50,'(3I6,2G18.10)') j,i,ik,dble(z1),aimag(z1)
!      amn(j,i,ik)=z1
!      j=j+1
!    end do
!  end do
!end do
!close(50)
deallocate(wftrial)
end subroutine

subroutine generate_trial_wavefunctions(THIS,norb_TB,num_bands,wftrial)
class(CLwan), intent(in) :: THIS
integer, intent(in)  :: norb_TB,num_bands
complex(dp), intent(out) :: wftrial(norb_TB,num_bands)
wftrial=0._dp
end subroutine


end module
