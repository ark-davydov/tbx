
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
real(dp), parameter :: epsr=1.e-1_dp

type, public :: CLwan
  integer :: nnk=0
  integer :: dis_num_iter=1000
  integer :: num_iter=1000
  real(dp) :: dis_win_min=-100._dp
  real(dp) :: dis_win_max= 100._dp
  real(dp) :: dis_froz_min=-100._dp
  real(dp) :: dis_froz_max= 100._dp
  real(dp), allocatable :: proj_centers(:,:)
  integer, allocatable :: nnkp(:,:,:)
  contains
  procedure :: init
  procedure :: projection
  procedure :: readnnkp
  procedure, private :: writewin
endtype CLwan

contains


subroutine init(THIS,kgrid,kpath,pars,eval)
class(CLwan), intent(inout) :: THIS
class(GRID), intent(in) :: kgrid
class(PATH), intent(in) :: kpath
class(CLpars), intent(in) :: pars
real(dp), intent(inout) :: eval(pars%nstates,kgrid%npt)
!
integer ik
real(dp), allocatable :: vkl(:,:)
if (trim(adjustl(pars%wannier_proj_mode)).eq.'tbg4band') then
  allocate(THIS%proj_centers(NDIM,pars%nwan))
  THIS%proj_centers(:,1)=(/onethrd,onethrd,0._dp/)
  THIS%proj_centers(:,2)=(/onethrd,onethrd,0._dp/)
  THIS%proj_centers(:,3)=(/twothrd,twothrd,0._dp/)
  THIS%proj_centers(:,4)=(/twothrd,twothrd,0._dp/)
  allocate(vkl(NDIM,kgrid%npt))
  do ik=1,kgrid%npt
    vkl(:,ik)=kgrid%vpl(ik)
  end do
  if (mp_mpi) then
    call THIS%writewin(kgrid,kpath,pars)
    call io_eval(1002,"write",trim(adjustl(pars%seedname))//'.eig',.true.,pars%nstates,kgrid%npt,pars%efermi,vkl,eval)
  end if
  call readnnkp(THIS,kgrid,pars)
  deallocate(vkl)
else
  call throw("wannier_interface%generate_trial_wavefunctions()","unknown projection option")
end if
end subroutine

subroutine projection(THIS,tbmodel,pars,sym,kgrid,evec)
class(CLwan), intent(in) :: THIS
class(CLtb), intent(in) :: tbmodel
class(CLpars), intent(in) :: pars
class(CLsym), intent(in) :: sym
class(GRID), intent(in) :: kgrid
complex(dp), intent(in) :: evec(tbmodel%norb_TB,pars%nstates,kgrid%npt)
! local
character(len=200) :: message
integer iorb,ispec
integer eps_plus_state_val,eps_plus_state_con
integer  isym1,isym2,ik_gamma,ikg(NDIM+1)
real(dp) sigma
real(dp) vpl(NDIM)
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
  ! 0.7 of lattice constant as in the paper of Johannes Lischner
  sigma=0.7_dp*sqrt(dot_product(pars%avec(1,:),pars%avec(1,:)))
  call generate_amn_overlap(tbmodel,pars,kgrid,evec,wftrial,.true.,THIS%proj_centers,sigma)
  call generate_mmn_overlap(THIS,tbmodel,pars,kgrid,evec)
  deallocate(psi_gamma_Eplus)
  deallocate(psi_gamma_Eminus)
else
  call throw("wannier_interface%generate_trial_wavefunctions()","unknown projection option")
end if
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
if (mp_mpi) then
  write(*,*) "rotation eigenvalues:"
  write(*,*) a1,b1,a2,b2
end if
if (a1*b1.lt.0._dp.and.a2*b2.lt.0._dp) then
  call throw("CLwan%project","both dublet state have rotation eigenvalues in 2 or 4th quadrant")
else if (a1*b1.gt.0._dp.and.a2*b2.gt.0._dp) then
  call throw("CLwan%project","both dublet state have rotation eigenvalues in 1 or 3d quadrant")
else if (a1*b1.lt.0._dp) then
  if ( ( abs(a1-cos(twopi3)).lt.epsr .and. abs(b1-sin(twopi3)).lt.epsr ) .or.&
       ( abs(b1-cos(twopi3)).lt.epsr .and. abs(a1-sin(twopi3)).lt.epsr ) ) then
     eps_plus_dublet_state=ist1
  else
    call throw("CLwan%project","2pi/3 rotation eigenvalue is not found where expected")
  end if
else if (a2*b2.lt.0._dp) then
  if ( ( abs(a2-cos(twopi3)).lt.epsr .and. abs(b2-sin(twopi3)).lt.epsr ) .or.&
       ( abs(b2-cos(twopi3)).lt.epsr .and. abs(a2-sin(twopi3)).lt.epsr ) ) then
    eps_plus_dublet_state=ist2
  else
    call throw("CLwan%project","2pi/3 rotation eigenvalue is not found where expected")
  end if
else
  call throw("CLwan%project","could not find state with 2pi/3 rotation eigenvalue")
end if
deallocate(wf_t)
end function

subroutine generate_amn_overlap(tbmodel,pars,kgrid,evec,wftrial,lgauss,gauss_centers,sigma)
class(CLtb), intent(in) :: tbmodel
class(CLpars), intent(in) :: pars
class(GRID), intent(in) :: kgrid
complex(dp), intent(in) :: evec(tbmodel%norb_TB,pars%nstates,kgrid%npt)
complex(dp), intent(in) :: wftrial(tbmodel%norb_TB,pars%nwan)
logical, intent(in) :: lgauss
real(dp), intent(in) :: gauss_centers(NDIM,pars%nwan)
real(dp), intent(in) :: sigma
! local
logical exs
integer ik,iwan,ist,iR,iorb
real(dp) t1,t2,dd
real(dp) dv(NDIM),dc(NDIM)
complex(dp) z2
complex(dp), allocatable :: amn(:,:)
! A_mn(k)=<\psi_m(k)|g_n>
inquire(file=trim(adjustl(pars%seedname))//'.amn',exist=exs)
if (exs) then
  call info("CLwan%generate_amn_overlap","skipping "//trim(adjustl(pars%seedname))//".amn creation")
  return
end if
if (mp_mpi) then
  open(50,file=trim(adjustl(pars%seedname))//'.amn',action='write')
  write(50,*) '# '//trim(adjustl(pars%seedname))//' file '
  write(50,*) pars%nstates,kgrid%npt,pars%nwan
end if
allocate(amn(pars%nstates,pars%nwan))
do ik=1,kgrid%npt
  amn=0._dp
  do iwan=1,pars%nwan
    do iR=1,tbmodel%rgrid%npt
      t2=dot_product(kgrid%vpl(ik),tbmodel%rgrid%vpl(iR))*twopi
      z2=cmplx(cos(t2),-sin(t2),kind=dp)
      do iorb=1,tbmodel%norb_TB
        if (lgauss) then
          dv=tbmodel%vplorb(iorb)+tbmodel%rgrid%vpl(iR)-gauss_centers(:,iwan)
          dc=matmul(dv,pars%avec)
          dd=sqrt(dot_product(dc,dc))
          t1=gauss(dd,sigma)
          if (t1.lt.epsengy) cycle
          amn(:,iwan)=amn(:,iwan)+conjg(evec(iorb,:,ik))*wftrial(iorb,iwan)*t1*z2
        else
          amn(:,iwan)=amn(:,iwan)+conjg(evec(iorb,:,ik))*wftrial(iorb,iwan)*z2
        end if
      end do
    end do
  end do
  do iwan=1,pars%nwan
    do ist=1,pars%nstates
      if (mp_mpi) write(50,'(3I6,2G18.10)') ist,iwan,ik,dble(amn(ist,iwan)),aimag(amn(ist,iwan))
    end do
  end do
end do
if (mp_mpi) close(50)
deallocate(amn)
end subroutine


subroutine generate_mmn_overlap(THIS,tbmodel,pars,kgrid,evec)
class(CLwan), intent(in) :: THIS
class(CLtb), intent(in) :: tbmodel
class(CLpars), intent(in) :: pars
class(GRID), intent(in) :: kgrid
complex(dp), intent(in) :: evec(tbmodel%norb_TB,pars%nstates,kgrid%npt)
! local
logical exs
integer ik,jk,mst,nst,iorb,innk
real(dp) dc,t1
complex(dp) z1
real(dp) vq(NDIM),vc(NDIM)
complex(dp), allocatable :: mmn(:,:,:)
! M_mn(k)=<u_mk|u_n{k+q}>
inquire(file=trim(adjustl(pars%seedname))//'.mmn',exist=exs)
if (exs) then
  call info("CLwan%generate_mmn_overlap","skipping "//trim(adjustl(pars%seedname))//".mmn creation")
  return
end if
if (mp_mpi) then
  open(50,file=trim(adjustl(pars%seedname))//'.mmn',action='write')
  write(50,*) '# '//trim(adjustl(pars%seedname))//' file '
  write(50,*) pars%nstates,kgrid%npt,THIS%nnk
end if
allocate(mmn(pars%nstates,pars%nstates,THIS%nnk))
do ik=1,kgrid%npt
  mmn=0._dp
  do innk=1,THIS%nnk
    jk=THIS%nnkp(4,innk,ik)
    vq=kgrid%vpl(jk)+dble(THIS%nnkp(1:3,innk,ik))-kgrid%vpl(ik)
    vc=matmul(vq,pars%avec)
    dc=sqrt(dot_product(vc,vc))
    do mst=1,pars%nstates
      do nst=1,pars%nstates
        do iorb=1,tbmodel%norb_TB
          t1=dot_product(vq,tbmodel%vplorb(iorb))*twopi
          z1=cmplx(cos(t1),-sin(t1),kind=dp)
          mmn(mst,nst,innk)=mmn(mst,nst,innk)+conjg(evec(iorb,mst,ik))*evec(iorb,nst,jk)*z1*pwave_ovlp(dc)
        end do
      end do
    end do
  end do
  do innk=1,THIS%nnk
    if (mp_mpi)  write(50,'(5I6)') ik,THIS%nnkp(4,innk,ik),THIS%nnkp(1:3,innk,ik)
    do nst=1,pars%nstates
      do mst=1,pars%nstates
        if (mp_mpi) write(50,'(2G18.10)') dble(mmn(mst,nst,innk)),aimag(mmn(mst,nst,innk))
      end do
    end do
  end do
end do
if (mp_mpi) close(50)
deallocate(mmn)
end subroutine

subroutine readnnkp(THIS,kgrid,pars)
class(CLwan), intent(inout) :: THIS
class(GRID), intent(in) :: kgrid
class(CLpars), intent(in) :: pars
integer i,nk,np,innk,ik
logical exs
inquire(file=trim(adjustl(pars%seedname))//'.nnkp',exist=exs)
if (.not.exs) then
  call throw("CLwan%readnnkp","file "//trim(adjustl(pars%seedname))&
            //".nnkp missing, run: wannier -pp "//trim(adjustl(pars%seedname))//".win")
end if
open(50,file=trim(adjustl(pars%seedname))//'.nnkp',action='read')
do i=1,17
  read(50,*)
end do
read(50,*) nk
if (nk.ne.kgrid%npt) call throw("CLwan%readnnkp","Wrong number of k-points in "//trim(adjustl(pars%seedname))//".nnkp")
do i=1,nk
  read(50,*)
end do
read(50,*) ; read(50,*); read(50,*)
read(50,*) np
do i=1,np
 read(50,*) ; read(50,*)
end do
read(50,*) ; read(50,*) ; read(50,*)
read(50,*) THIS%nnk
allocate(THIS%nnkp(4,THIS%nnk,kgrid%npt))
do ik=1,kgrid%npt
  do innk=1,THIS%nnk
    read(50,*) i,THIS%nnkp(4,innk,ik),THIS%nnkp(1,innk,ik),THIS%nnkp(2,innk,ik),THIS%nnkp(3,innk,ik)
  end do
end do
close(50)
return
end subroutine

subroutine writewin(THIS,kgrid,kpath,pars)
class(CLwan), intent(inout) :: THIS
class(GRID), intent(in) :: kgrid
class(PATH), intent(in) :: kpath
class(CLpars), intent(in) :: pars
integer iwan,ik,ivert
logical exs
inquire(file=trim(adjustl(pars%seedname))//'.win',exist=exs)
if (exs) then
  call info("CLwan%writewin","skipping "//trim(adjustl(pars%seedname))//".win creation")
  return
end if
call info("CLwan%writewin","creating "//trim(adjustl(pars%seedname))//".win")
open(50,file=trim(adjustl(pars%seedname))//'.win',action='write')
write(50,*) 'bands_plot        = .true.'
write(50,*) '!dos             = .true.'
write(50,*) '!dos_kmesh       = 150 150 1'
write(50,*) '!dos_energy_min  = -0.01'
write(50,*) '!dos_energy_max  =  0.01'
write(50,*) '!dos_energy_step = 0.00001'
write(50,*) '!dos_adpt_smr    = .true.'
write(50,*) '!adpt_smr_max    = 0.0001'
write(50,*) '!restart=wannierise'
write(50,*) 'use_ws_distance   = .true.'
write(50,*) 'write_hr          = .true.'
write(50,*) 'write_tb          = .true.'
write(50,*) 'write_u_matrices  = .true.'
write(50,*) 'write_xyz         = .true.'
write(50,*) 'trial_step        = 0.2'
write(50,*) '!slwf_constrain    = true'
write(50,*) '!slwf_lambda       = 200'


write(50,*)
write(50,*) 'dis_win_min       =',THIS%dis_win_min
write(50,*) 'dis_win_max       =',THIS%dis_win_max
write(50,*) 'dis_froz_min      =',THIS%dis_froz_min
write(50,*) 'dis_froz_max      =',THIS%dis_froz_max
write(50,*)

write(50,*) '!slwf_num          = ',pars%nwan
write(50,*) 'num_bands         = ',pars%nstates
write(50,*) 'num_wann          = ',pars%nwan
write(50,*) 'dis_mix_ratio     = 0.8'
write(50,*) 'dis_num_iter      = ',THIS%dis_num_iter
write(50,*) 'num_iter          = ',THIS%num_iter
write(50,*) 'num_print_cycles  = 200'
write(50,*) 'search_shells     = 1000'
write(50,*) 'kmesh_tol         = 0.000001'
write(50,*) 'mp_grid           = ',kgrid%ngrid
write(50,*)
write(50,*) 'Begin Projections'
do iwan=1,pars%nwan
  write(50,'("f=",G18.10,",",G18.10,",",G18.10,":",A4)') THIS%proj_centers(:,iwan),'pz'
end do
write(50,*) 'End Projections'
write(50,*)
write(50,*) 'begin kpoint_path'
do ivert=1,kpath%nvert-1
  write(50,'("X ",3F12.8," X ",3F12.8)') kpath%vert(:,ivert),kpath%vert(:,ivert+1)
end do
write(50,*) 'end kpoint_path'
write(50,*)
write(50,*)'Begin Unit_Cell_Cart'
write(50,'(3G18.10)') pars%avec(1,:)
write(50,'(3G18.10)') pars%avec(2,:)
write(50,'(3G18.10)') pars%avec(3,:)
write(50,*)'End Unit_Cell_Cart'
write(50,*)
write(50,*)'Begin Atoms_Frac'
do iwan=1,pars%nwan
  write(50,'(A2,3G18.10)')'XX',THIS%proj_centers(:,iwan)
end do
write(50,*)'End Atoms_Frac'
write(50,*)
write(50,*)'Begin kpoints'
do ik=1,kgrid%npt
  write(50,'(3F12.8)') kgrid%vpl(ik)
end do
write(50,*)'End kpoints'
close(50)
return
end subroutine


real(dp) function pwave_ovlp(dc)
real(dp), intent(in) :: dc
real(dp), parameter :: charge_pz=3.18_dp
pwave_ovlp=( 1._dp/( 1._dp+(dc*abohr/charge_pz)**2 ) )**3
end function

end module
