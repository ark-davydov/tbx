
module wannier_interface
#ifdef MPI
  use mpi
#endif
use modcom
use parameters
use tbclass
use gridclass
use symmetryclass
implicit none
private

type, public :: CLwan
  contains
  procedure, nopass :: project=>compute_projection
endtype CLwan

contains

subroutine compute_projection(tbmodel,pars,sym,kgrid,eval,evec)
class(CLtb), intent(in) :: tbmodel
class(CLpars), intent(in) :: pars
class(CLsym), intent(in) :: sym
class(GRID), intent(in) :: kgrid
real(dp), intent(in) :: eval(pars%nstates,kgrid%npt)
complex(dp), intent(in) :: evec(tbmodel%norb_TB,pars%nstates,kgrid%npt)
! local
character(len=200) :: message
integer iorb,ispec
integer eps_plus_state_val,eps_plus_state_con
integer  isym1,isym2,ik_gamma,ikg(NDIM+1)
real(dp) vpl(NDIM)
real(dp) rot(NDIM,NDIM),ax(NDIM),res(NDIM,NDIM)
complex(dp), allocatable :: wftrial(:,:)
complex(dp), allocatable :: psi_gamma_Eplus(:)
complex(dp), allocatable :: psi_gamma_Eminus(:)
! find gamma point via general k-point finder subroutine 
vpl=0._dp
ikg=kgrid%find(vpl)
ik_gamma=ikg(NDIM+1)
allocate(wftrial(tbmodel%norb_TB,pars%nwan))
wftrial(:,:)=0._dp
if (trim(adjustl(pars%wannier_proj_mode)).eq.'tbg4band') then
  if (NDIM.ne.3) call throw("wannier_interface%generate_trial_wavefunctions()","this subroutine assumes NDIM=3")
  if (pars%nstates.ne.4) call throw("wannier_interface%generate_trial_wavefunctions()",&
    & "this subroutine assumes num_bands=4, but another number is given &
    & (to resolve the issue, recompute eigen values/vectors with 'states' block providing nstates=4)")
  if (pars%nwan.ne.4) call throw("wannier_interface%generate_trial_wavefunctions()",&
                      "one needs to put nwan=4 to generate the 'tbg4band' trial projection")
  ! PHYSICAL REVIEW X 8, 031088 (2018)
  ! first we find \psi_{\Gamma,E+,\epsilon} and \psi_{\Gamma,E-,\epsilon} 
  ! where E+ is dublet above Ef and E- is dublet below Ef
  ! \epsilon is the "phase factor" of a+ib, where a=<psi_1|S|\psi>, b=<psi_2|S|\psi>, for any \psi, while psi_1 and psi_2
  ! are dublet wavefunction; under a C3 rotation, which we will find from the symmetry operations table, \epsilon
  ! is + or - 2pi/3
  isym1=2
  isym2=3
  if (mp_mpi) then
    write(message,*)"assiming isym=",isym1,"to be 2pi/3 rotation, 1,2 states=valence dublet states, 3,4=conduction"
    call info("CLwan%project",trim(message))
    write(message,*)"assiming isym=",isym2,"to be C2' from the paper, which generates w3 "
    call info("CLwan%project",trim(message))
  end if
  eps_plus_state_val=eps_plus_dublet_state(tbmodel,pars,sym,isym1,1,2,evec(:,:,ik_gamma))
  eps_plus_state_con=eps_plus_dublet_state(tbmodel,pars,sym,isym1,3,4,evec(:,:,ik_gamma))
  allocate(psi_gamma_Eplus(tbmodel%norb_TB))
  allocate(psi_gamma_Eminus(tbmodel%norb_TB))
  psi_gamma_Eplus=0._dp
  psi_gamma_Eminus=0._dp
  do iorb=1,tbmodel%norb_TB
    ispec=tbmodel%orb_icio(iorb,3)
    vpl=tbmodel%vplorb(iorb)
    if (vpl(3).gt.0._dp) then
      ! top layer
      if (ispec.eq.1) then
        ! A-site
        psi_gamma_Eplus(iorb)=evec(iorb,eps_plus_state_con,ik_gamma)
      else
        ! B-site
        psi_gamma_Eminus(iorb)=evec(iorb,eps_plus_state_val,ik_gamma)
      end if 
    else
      ! bottom layer
      if (ispec.eq.1) then
        ! A-site
        psi_gamma_Eminus(iorb)=evec(iorb,eps_plus_state_val,ik_gamma)
      else
        ! B-site
        psi_gamma_Eplus(iorb)=evec(iorb,eps_plus_state_con,ik_gamma)
      end if 
    end if
  end do
  wftrial(:,1)=psi_gamma_Eplus+psi_gamma_Eminus
  wftrial(:,2)=conjg(wftrial(:,1))
  wftrial(:,3)=wftrial(:,1)
  call tbmodel%wfGtransform(pars,sym,isym2,wftrial(:,3))
  wftrial(:,4)=conjg(wftrial(:,3))
  deallocate(psi_gamma_Eplus)
  deallocate(psi_gamma_Eminus)
else
  call throw("wannier_interface%generate_trial_wavefunctions()","unknown projection option")
end if
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

integer function eps_plus_dublet_state(tbmodel,pars,sym,isym,ist1,ist2,evec_G)
class(CLtb), intent(in) :: tbmodel
class(CLpars), intent(in) :: pars
class(CLsym), intent(in) :: sym
integer, intent(in) :: isym,ist1,ist2
complex(dp), intent(in) :: evec_G(tbmodel%norb_TB,pars%nstates)
real(dp) twopi3,a1,a2,b1,b2
complex(dp), allocatable :: wf_t(:)
! constant
twopi3=twopi/3._dp
! find one of the lower dublet state (E-) with eigenvalus +\epsilon
allocate(wf_t(tbmodel%norb_TB))
wf_t=evec_G(:,ist1)
call tbmodel%wfGtransform(pars,sym,isym,wf_t)
a1=dble(dot_product(evec_G(:,ist1),wf_t))
b1=dble(dot_product(evec_G(:,ist2),wf_t))
wf_t=evec_G(:,ist2)
call tbmodel%wfGtransform(pars,sym,isym,wf_t)
a2=dble(dot_product(evec_G(:,ist1),wf_t))
b2=dble(dot_product(evec_G(:,ist2),wf_t))
if (a1*b1.lt.0._dp.and.a2*b2.lt.0._dp) then
  call throw("CLwan%project","both dublet state have rotation eigenvalues in 2 or 4th quadrant")
else if (a1*b1.gt.0._dp.and.a2*b2.gt.0._dp) then
  call throw("CLwan%project","both dublet state have rotation eigenvalues in 1 or 3d quadrant")
else if (a1*b1.lt.0._dp) then
  if ( ( abs(a1-cos(twopi3)).lt.epslat .and. abs(b1-sin(twopi3)).lt.epslat ) .or.&
       ( abs(b1-cos(twopi3)).lt.epslat .and. abs(a1-sin(twopi3)).lt.epslat ) ) then
     eps_plus_dublet_state=ist1
  else
    call throw("CLwan%project","2pi/3 rotation eigenvalue is not found where expected")
  end if
else if (a2*b2.lt.0._dp) then
  if ( ( abs(a2-cos(twopi3)).lt.epslat .and. abs(b2-sin(twopi3)).lt.epslat ) .or.&
       ( abs(b2-cos(twopi3)).lt.epslat .and. abs(a2-sin(twopi3)).lt.epslat ) ) then
    eps_plus_dublet_state=ist2
  else
    call throw("CLwan%project","2pi/3 rotation eigenvalue is not found where expected")
  end if
else
  call throw("CLwan%project","could not find state with 2pi/3 rotation eigenvalue")
end if
call generate_amn_overlap(tbmodel,pars,sym,kgrid,eval,evec,wftrial)
deallocate(wf_t)
end function

subroutine generate_amn_overlap(tbmodel,pars,sym,kgrid,eval,evec,wftrial)
class(CLtb), intent(in) :: tbmodel
class(CLpars), intent(in) :: pars
class(CLsym), intent(in) :: sym
class(GRID), intent(in) :: kgrid
real(dp), intent(in) :: eval(pars%nstates,kgrid%npt)
complex(dp), intent(in) :: evec(tbmodel%norb_TB,pars%nstates,kgrid%npt)
complex(dp), intent(in) :: wftrial(tbmodel%norb_TB,pars%nwan)
real(dp) cntr1(NDIM),cntr2(NDIM),cntr1c(NDIM),cntr2c(NDIM)
  cntr1=(/onethrd,onethrd,0._dp/)
  cntr2=(/twothrd,twothrd,0._dp/)
  cntr1c=matmul(cntr1,pars%avec)
  cntr2c=matmul(cntr2,pars%avec)
end subroutine

end module
