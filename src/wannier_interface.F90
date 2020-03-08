
module wannier_interface
#ifdef MPI
  use mpi
#endif
use modcom
use parameters
use tbclass
use gridclass
use symmetryclass
use wannier_supplementary
implicit none
private
real(dp), parameter :: epsr=4.e-1_dp

type, public :: CLwan
  integer :: nnk=0
  integer :: dis_num_iter=1000
  integer :: num_iter=1000
  real(dp) :: dis_win_min=-100._dp
  real(dp) :: dis_win_max= 100._dp
  real(dp) :: dis_froz_min=-100._dp
  real(dp) :: dis_froz_max= 100._dp
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
!external string_to_lmr
real(dp), allocatable :: vkl(:,:)
if (.not.pars%proj%allocatd) call throw("wannier_interface%projection",&
 "projection block must be allocated with correct number of porjections")
if (trim(adjustl(pars%wannier_proj_mode)).eq.'tbg4band'.or.&
    trim(adjustl(pars%wannier_proj_mode)).eq.'input_file') then
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
class(CLwan), intent(inout) :: THIS
class(CLtb), intent(in) :: tbmodel
class(CLpars), intent(in) :: pars
class(CLsym), intent(in) :: sym
class(GRID), intent(inout) :: kgrid
complex(dp), intent(in) :: evec(tbmodel%norb_TB,pars%nstates,kgrid%npt)
! local
character(len=200) :: message
integer iorb,jorb,ispec
integer ipro,ic,jc,ios,l1,m1,l2,m2
integer eps_plus_state_val,eps_plus_state_con
integer  isym1,isym2,ik_gamma,ikg(NDIM+1)
type(wbase) proj
real(dp) sigma
real(dp) vpl(NDIM),x1(NDIM),x2(NDIM),z1(NDIM),z2(NDIM)
complex(dp), allocatable :: wftrial(:,:)
complex(dp), allocatable :: psi_gamma_Eplus(:)
complex(dp), allocatable :: psi_gamma_Eminus(:)
call proj%init(pars,pars%proj%ncenters,pars%proj%norb,pars%proj%norb_ic,&
                   pars%proj%lmr,pars%proj%waxis,pars%proj%centers)
allocate(wftrial(tbmodel%norb_TB,pars%proj%norb))
if (trim(adjustl(pars%wannier_proj_mode)).eq.'tbg4band') then
  if (NDIM.ne.3) call throw("wannier_interface%generate_trial_wavefunctions()","this subroutine assumes NDIM=3")
  if (pars%nstates.ne.4) call throw("wannier_interface%generate_trial_wavefunctions()",&
    & "this subroutine assumes num_bands=4, but another number is given &
    & (to resolve the issue, recompute eigen values/vectors with 'states' block providing nstates=4)")
  if (pars%proj%norb.ne.4) call throw("wannier_interface%generate_trial_wavefunctions()",&
                      "one needs to put 4 trial orbitals in 'projections' block")
  ! PHYSICAL REVIEW X 8, 031088 (2018)
  ! first we find \psi_{\Gamma,E+,\epsilon} and \psi_{\Gamma,E-,\epsilon} 
  ! where E+ is dublet above Ef and E- is dublet below Ef
  ! \epsilon is the "phase factor" of a+ib, where a=<psi_1|S|\psi>, b=<psi_2|S|\psi>, for any \psi, while psi_1 and psi_2
  ! are dublet wavefunction; under a C3 rotation, which we will find from the symmetry operations table, \epsilon
  ! is + or - 2pi/3
  !
  ! find gamma point via general k-point finder subroutine 
  vpl=0._dp
  ikg=kgrid%find(vpl)
  ik_gamma=ikg(NDIM+1)
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
    ispec=tbmodel%orb_ispec(iorb)
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
  call generate_amn_overlap(tbmodel,pars,kgrid,evec,wftrial,.true.,sigma)
  deallocate(psi_gamma_Eplus)
  deallocate(psi_gamma_Eminus)
else if (trim(adjustl(pars%wannier_proj_mode)).eq.'input_file') then
  if (NDIM.ne.3) call throw("wannier_interface%generate_trial_wavefunctions()","this subroutine assumes NDIM=3")
   wftrial(:,:)=0._dp
   do ic=1,proj%ncenters
     do ios=1,proj%norb_ic(ic)
       ipro=proj%icio_orb(ic,ios)
       ! basis of TB hamiltonian => wavefunctions
       do jorb=1,tbmodel%norb_TB
         jc=tbmodel%wbase%orb_icio(jorb,1)
         if ( sum(abs( tbmodel%wbase%centers(:,jc)-proj%centers(:,ic) )).gt.epslat) cycle
         l1=tbmodel%wbase%lmr(1,jorb)
         m1=tbmodel%wbase%lmr(2,jorb)
         x1=tbmodel%wbase%waxis(:,1,jorb)
         z1=tbmodel%wbase%waxis(:,2,jorb)
         l2=proj%lmr(1,ipro)
         m2=proj%lmr(2,ipro)
         x2=proj%waxis(:,1,ipro)
         z2=proj%waxis(:,2,ipro)
         wftrial(jorb,ipro)=tbmodel%wbase%wws_full(sym%car(:,:,1),l1,m1,l2,m2,x1,z1,x2,z2)
       end do
     end do
   end do
   call generate_amn_overlap(tbmodel,pars,kgrid,evec,wftrial,.false.,sigma)
   do ipro=1,proj%norb
       do jorb=1,tbmodel%norb_TB
         write(*,*) ipro,jorb,dble( wftrial(jorb,ipro)),aimag( wftrial(jorb,ipro))
       end do
   end do
   call generate_dmn_orb(tbmodel,proj,sym,pars,kgrid,evec)
   !stop
else
  call throw("wannier_interface%generate_trial_wavefunctions()","unknown projection option")
end if
call generate_mmn_overlap(THIS,tbmodel,pars,kgrid,evec)
deallocate(wftrial)
end subroutine

subroutine generate_dmn_orb(tbmodel,base,sym,pars,kgrid,evec)
class(CLtb), intent(in) :: tbmodel
class(CLsym), intent(in) :: sym
class(wbase), intent(in) :: base
class(CLpars), intent(in) :: pars
class(GRID), intent(inout) :: kgrid
complex(dp), intent(in) :: evec(tbmodel%norb_TB,pars%nstates,kgrid%npt)
real(dp), allocatable :: wws(:,:,:)
! local
integer isym,iw,jw,ic,jc,ir,ik
integer ikp,ist,jst
real(dp) err,t1
real(dp) v1(NDIM),v2(NDIM),v3(NDIM),v4(NDIM)
logical lfound(base%ncenters)
integer, allocatable :: ics2c(:,:)
real(dp), allocatable :: vcs2t(:,:,:)
complex(dp), allocatable :: phs(:,:)
complex(dp), allocatable :: wf_t(:)
call kgrid%sym_init(sym)
allocate(wws(base%norb,base%norb,sym%nsym))
allocate(vcs2t(3,base%ncenters,sym%nsym))
allocate(ics2c(base%ncenters,sym%nsym))
ics2c=-999
write(*,"(a,i5)") "  Number of symmetry operators = ", sym%nsym
do isym=1,sym%nsym
   write(*,"(2x,i5,a)") isym, "-th symmetry operators is"
   write(*,"(3f15.7)") &
      sym%car(:,:,isym), sym%vtc(:,isym)/2.4629437 !Writing rotation matrix and translation vector in Cartesian coordinates
   lfound=.false.
   do ic=1,base%ncenters
    v1=base%centers(:,ic)
    v2=matmul(sym%lat(:,:,isym),(v1+sym%vtl(:,isym)))
    !v2=matmul(sym%lat(:,:,isym),v1)+sym%vtl(:,isym)
    do jc=1,base%ncenters
      if(lfound(jc)) cycle
      v3=base%centers(:,jc)
      !v4=matmul(v3-v2,pars%bvec)
      v4=v3-v2
      if(sum(abs(dble(nint(v4))-v4)).lt.1d-2) then
         lfound(jc)=.true.
         ics2c(ic,isym)=jc
         exit !Sym.op.(isym) moves position(ips2p(ip,isym)) to position(ip) + T, where
      end if                                       !T is given by vps2t(:,ip,isym).
    end do
    if(ics2c(ic,isym).le.0) then
       write(*,"(a,3f18.10)")"  Could not find ",v2
       write(*,"(a,3f18.10)")"  coming from    ",v1
       write(*,"(a,i5,a               )")"  of Wannier site",ic,"."
       write(*,"(a,i5,a               )")"  symmetry ",isym,"."
       write(*,"(3I6)") sym%lat(1,:,isym)
       write(*,"(3I6)") sym%lat(2,:,isym)
       write(*,"(3I6)") sym%lat(3,:,isym)
       call throw("wannier_interface%generate_dmn_orb", "missing Wannier sites, see the output.")
    end if
  end do
  do ic=1,base%ncenters
    write(112,*) isym,ic,ics2c(ic,isym)
  end do
end do
vcs2t=0._dp
do isym=1,sym%nsym
  do ic=1,base%ncenters
    v1=base%centers(:,ic)
    jc=ics2c(ic,isym)
    v2=base%centers(:,jc)
    v3=matmul(v2,sym%lat(:,:,isym))-sym%vtl(:,isym)
    vcs2t(:,ic,isym)=v3-v1
  end do
end do
wws(:,:,:)=0._dp
do isym=1,sym%nsym
   do iw=1,base%norb
      ic=base%orb_icio(iw,1)
      jc=ics2c(ic,isym)
      do jw=1,base%norb
         if(base%orb_icio(jw,1).ne.jc) cycle
         wws(jw,iw,isym)=base%wws(sym%car(:,:,isym),iw,jw)
      end do
   end do
   do iw=1,base%norb
      err=abs((sum(wws(:,iw,isym)**2)+sum(wws(iw,:,isym)**2))*.5_dp-1._dp)
      if(err.gt.1.e-3_dp) then
         write(*,*) "wannier_interface%generate_dmn_orb","compute_dmn: Symmetry operator (",isym, &
                 ") could not transform Wannier function (",iw,")."
         write(*,*) "The error is ",err,"."
         call throw("wannier_interface%generate_dmn_orb", "missing Wannier functions, see the output.")
      end if
   end do
end do
if (mp_mpi) then
  open (unit=1001, file=trim(adjustl(pars%seedname))//".dmn",form='formatted')
  write (1001,*) '# '//trim(adjustl(pars%seedname))//'.dmn file'
  write (1001,"(4i9)") pars%nstates, sym%nsym, kgrid%nir, kgrid%npt
  write (1001,*)
  write (1001,"(10i9)") kgrid%ik2ir(1:kgrid%npt)
  write (1001,*)
  write (1001,"(10i9)") kgrid%ir2ik(1:kgrid%nir)
  do ir=1,kgrid%nir
     write (1001,*)
     write (1001,"(10i9)") kgrid%iks2k(kgrid%ir2ik(ir),:)
  enddo
end if
allocate(phs(base%norb,base%norb))
phs=0._dp
if (mp_mpi) then
  WRITE(*,'(/)')
  WRITE(*,'(a,i8)') '  DMN(d_matrix_wann): nir = ',kgrid%nir
end if
do ir=1,kgrid%nir
  ik=kgrid%ir2ik(ir)
  if (mp_mpi) then
    WRITE (*,'(i8)',advance='no') ir
    IF( MOD(ir,10) == 0 ) WRITE (*,*)
  end if
  !FLUSH(*)
  do isym=1,sym%nsym
     do iw=1,base%norb
        ic=base%orb_icio(iw,1)
        jc=ics2c(ic,sym%inv(isym))
        v1=kgrid%vpl(kgrid%iks2k(ik,isym))-matmul(sym%lat(:,:,isym),kgrid%vpl(ik))
        v2=matmul(v1,sym%lat(:,:,isym))
        !Phase of T.k with lattice vectors T of aove.
        t1=dot_product(vcs2t(:,jc,isym),kgrid%vpl(ik))*twopi
        phs(iw,iw)=cmplx(cos(t1),sin(t1),kind=dp) !Phase of t.G with translation vector t(isym).
     end do
     if (mp_mpi) then
        WRITE (1001,*)
        WRITE (1001,"(1p,(' (',e18.10,',',e18.10,')'))") matmul(phs,cmplx(wws(:,:,isym),0._dp,kind=dp))
     end if
  end do
end do
if (mp_mpi) then
  write(1001,*)
  write(*,'(/)')
  write(*,'(a,i8)') '  DMN(d_matrix_band): nir = ',kgrid%nir
end if
allocate(wf_t(tbmodel%norb_TB))
do ir=1,kgrid%nir
  if (mp_mpi) then
    WRITE (*,'(i8)',advance='no') ir
    IF( MOD(ir,10) == 0 ) WRITE (*,*)
  end if
  do isym=1,sym%nsym
    ikp=kgrid%iks2k(ik,isym)
    t1=dot_product(kgrid%vpl(ik),sym%vtl(:,isym))*twopi
    do jst=1,pars%nstates
      wf_t=evec(:,jst,ik)
      call tbmodel%wfGtransform(pars,sym,isym,wf_t)
      do ist=1,pars%nstates
        if (mp_mpi) write(1001,"(1p,(' (',e18.10,',',e18.10,')'))") &
                    dot_product(evec(:,ist,ikp),wf_t)*cmplx(cos(t1),sin(t1),kind=dp)
      end do
    end do
  end do
end do
if (mp_mpi) close(1001)
deallocate(wws,phs)
deallocate(vcs2t,ics2c)
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
else if ( ( abs(a1-cos(twopi3)).lt.epsr .and. abs(b1-sin(twopi3)).lt.epsr ) .or.&
       ( abs(b1-cos(twopi3)).lt.epsr .and. abs(a1-sin(twopi3)).lt.epsr ) ) then
     eps_plus_dublet_state=ist1
else if ( ( abs(a2-cos(twopi3)).lt.epsr .and. abs(b2-sin(twopi3)).lt.epsr ) .or.&
       ( abs(b2-cos(twopi3)).lt.epsr .and. abs(a2-sin(twopi3)).lt.epsr ) ) then
     eps_plus_dublet_state=ist2
else
  call throw("CLwan%project","could not find state with 2pi/3 rotation eigenvalue")
end if
deallocate(wf_t)
end function

subroutine generate_amn_overlap(tbmodel,pars,kgrid,evec,wftrial,lgauss,sigma)
class(CLtb), intent(in) :: tbmodel
class(CLpars), intent(in) :: pars
class(GRID), intent(in) :: kgrid
complex(dp), intent(in) :: evec(tbmodel%norb_TB,pars%nstates,kgrid%npt)
complex(dp), intent(in) :: wftrial(tbmodel%norb_TB,pars%proj%norb)
logical, intent(in) :: lgauss
real(dp), intent(in) :: sigma
! local
logical exs
integer ik,iwan,ist,iR,iorb
real(dp) t1,dd!,t2
real(dp) dv(NDIM),dc(NDIM)
!complex(dp) z2
complex(dp), allocatable :: amn(:,:)
! A_mn(k)=<\psi_m(k)|g_n>
inquire(file=trim(adjustl(pars%seedname))//'.amn',exist=exs)
if (exs) then
  call info("CLwan%generate_amn_overlap","skipping "//trim(adjustl(pars%seedname))//".amn creation")
  !return
end if
if (mp_mpi) then
  open(50,file=trim(adjustl(pars%seedname))//'.amn',action='write')
  write(50,*) '# '//trim(adjustl(pars%seedname))//' file '
  write(50,*) pars%nstates,kgrid%npt,pars%proj%norb
end if
allocate(amn(pars%nstates,pars%proj%norb))
do ik=1,kgrid%npt
  amn=0._dp
  do iwan=1,pars%proj%norb
 !   do iR=1,tbmodel%rgrid%npt
!      t2=dot_product(kgrid%vpl(ik),tbmodel%rgrid%vpl(iR))*twopi
!      z2=cmplx(cos(t2),-sin(t2),kind=dp)
      do iorb=1,tbmodel%norb_TB
        if (lgauss) then
          dv=tbmodel%vplorb(iorb)+tbmodel%rgrid%vpl(iR)-pars%proj%centers(:,pars%proj%iw2ic(iwan))
          dc=matmul(dv,pars%avec)
          dd=sqrt(dot_product(dc,dc))
          t1=gauss(dd,sigma)
          if (t1.lt.epsengy) cycle
          amn(:,iwan)=amn(:,iwan)+conjg(evec(iorb,:,ik))*wftrial(iorb,iwan)*t1!*z2
        else
          amn(:,iwan)=amn(:,iwan)+conjg(evec(iorb,:,ik))*wftrial(iorb,iwan)!*z2
        end if
      end do
!    end do
  end do
  do iwan=1,pars%proj%norb
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

write(50,*) '!slwf_num          = ',pars%proj%norb
write(50,*) 'num_bands         = ',pars%nstates
write(50,*) 'num_wann          = ',pars%proj%norb
write(50,*) 'dis_mix_ratio     = 0.8'
write(50,*) 'dis_num_iter      = ',THIS%dis_num_iter
write(50,*) 'num_iter          = ',THIS%num_iter
write(50,*) 'num_print_cycles  = 200'
write(50,*) 'search_shells     = 1000'
write(50,*) 'kmesh_tol         = 0.000001'
write(50,*) 'mp_grid           = ',kgrid%ngrid
write(50,*)
write(50,*) 'Begin Projections'
do iwan=1,pars%proj%norb
  write(50,'("f=",G18.10,",",G18.10,",",G18.10,":",A4)') &
    pars%proj%centers(:,pars%proj%iw2ic(iwan)),lmr_to_string(pars%proj%lmr(:,iwan))
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
do iwan=1,pars%proj%norb
  write(50,'(A2,3G18.10)')'XX',pars%proj%centers(:,iwan)
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
