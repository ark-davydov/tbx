
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
public :: write_wf_universe,read_wfmloc

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
    trim(adjustl(pars%wannier_proj_mode)).eq.'tbg12band'.or.&
    trim(adjustl(pars%wannier_proj_mode)).eq.'wannier_file'.or.&
    trim(adjustl(pars%wannier_proj_mode)).eq.'real_space'.or.&
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
integer iR,iw,jw,i1
integer iorb,jorb,ispec
integer ipro,ic,jc,ios,l1,m1,l2,m2
integer isym1,isym2,ik_gamma,iR_zero,ikg(NDIM+1)
integer nr,pm_val(2),pm_con(2)
type(wbase) proj
real(dp) sigma
real(dp) vpl(NDIM),x1(NDIM),x2(NDIM),z1(NDIM),z2(NDIM)
real(dp) t1,dd,dv(NDIM),dc(NDIM)
complex(dp) zz
!complex(dp), allocatable :: wf_t(:)
complex(dp), allocatable :: wftrial(:,:,:)
complex(dp), allocatable :: wws(:,:,:)
if (NDIM.ne.3) call throw("wannier_interface%generate_trial_wavefunctions()","this subroutine assumes NDIM=3")
if (pars%proj%norb.le.0) then
   call throw("wannier_interface%generate_trial_wavefunctions()",&
              "apparently 'projections' block was not specified, wannier projections not found")
else
  if (pars%nstates.lt.pars%proj%norb) then
    call throw("wannier_interface%generate_trial_wavefunctions()",&
               "number of bands is less than the number of requested projections")
  end if
end if
call proj%init(pars,pars%proj%ncenters,pars%proj%norb,pars%proj%norb_ic,&
                   pars%proj%lmr,pars%proj%waxis,pars%proj%centers)
call proj%init_smap(sym,pars)
! find gamma point via general k-point finder subroutine 
vpl=0._dp
ikg=kgrid%find(vpl)
ik_gamma=ikg(NDIM+1)
! find home unit cell
ikg=tbmodel%rgrid%find(vpl)
iR_zero=ikg(NDIM+1)
allocate(wws(proj%norb,proj%norb,sym%nsym))
call kgrid%sym_init(sym)
if (trim(adjustl(pars%wannier_proj_mode)).eq.'tbg4band'.or.&
    trim(adjustl(pars%wannier_proj_mode)).eq.'tbg12band') then
  if (trim(adjustl(pars%wannier_proj_mode)).eq.'tbg4band') then
    if (pars%nstates.ne.4) then
        call throw("wannier_interface%generate_trial_wavefunctions()",&
      & "wannier_proj_mode='tbg4band' assumes num_bands=4, but another number is given &
      & (to resolve the issue, recompute eigen values/vectors with 'states' block providing nstates=4)")
    end if
    if (pars%proj%norb.ne.4) then
       call throw("wannier_interface%generate_trial_wavefunctions()",&
                  "one needs to put 4 trial orbitals in 'projections' block for wannier_proj_mode='tbg4band'")
    end if
  else if (trim(adjustl(pars%wannier_proj_mode)).eq.'tbg12band') then
    if (pars%proj%norb.ne.12) then
       call throw("wannier_interface%generate_trial_wavefunctions()",&
                  "one needs to put 12 trial orbitals in 'projections' block for wannier_proj_mode='tbg12band'")
    end if
    if (pars%iflat_band.le.0) then
       call throw("wannier_interface%generate_trial_wavefunctions()",&
                  "in 12-bands projection mode, one needs to specify the index of the "//&
                  "flat band ('iflat_band' block) with respect to the lowest available state"//&
                  "it can be deduced from eval.dat file")
    end if
  end if
  ! PHYSICAL REVIEW X 8, 031088 (2018)
  ! first we find \psi_{\Gamma,E+,\epsilon} and \psi_{\Gamma,E-,\epsilon} 
  ! where E+ is dublet above Ef and E- is dublet below Ef
  ! \epsilon is the "phase factor" of a+ib, where a=<psi_1|S|\psi>, b=<psi_2|S|\psi>, for any \psi, while psi_1 and psi_2
  ! are dublet wavefunction; under a C3 rotation, which we will find from the symmetry operations table, \epsilon
  ! is + or - 2pi/3
  !
  isym1=2
  isym2=4
  if (mp_mpi) then
    write(message,*)"assiming isym=",isym1,"to be 2pi/3 rotation, 1,2 states=valence dublet states, 3,4=conduction"
    call info("CLwan%project",trim(message))
    write(message,*)"assiming isym=",isym2,"to be C2' from the paper, which generates w3 "
    call info("CLwan%project",trim(message))
  end if
  i1=pars%iflat_band
  if (i1.le.0) i1=1
  pm_val=eps_pm_dublet_state(tbmodel,pars,kgrid,ik_gamma,sym,isym1,i1  ,i1+1,evec(:,:,ik_gamma))
  pm_con=eps_pm_dublet_state(tbmodel,pars,kgrid,ik_gamma,sym,isym1,i1+2,i1+3,evec(:,:,ik_gamma))
  allocate(wftrial(tbmodel%norb_TB,pars%proj%norb,tbmodel%rgrid%npt))
  do iorb=1,tbmodel%norb_TB
    ispec=tbmodel%orb_ispec(iorb)
    vpl=tbmodel%vplorb(iorb)
    if (vpl(3).gt.0._dp) then
      ! top layer
      if (ispec.eq.1) then
        ! A-site
        wftrial(iorb,1,1)=evec(iorb,pm_con(1),ik_gamma)
        wftrial(iorb,2,1)=evec(iorb,pm_con(2),ik_gamma)
      else
        ! B-site
        wftrial(iorb,1,1)=evec(iorb,pm_val(1),ik_gamma)
        wftrial(iorb,2,1)=evec(iorb,pm_val(2),ik_gamma)
      end if 
    else
      ! bottom layer
      if (ispec.eq.1) then
        ! A-site
        wftrial(iorb,1,1)=evec(iorb,pm_val(1),ik_gamma)
        wftrial(iorb,2,1)=evec(iorb,pm_val(2),ik_gamma)
      else
        ! B-site
        wftrial(iorb,1,1)=evec(iorb,pm_con(1),ik_gamma)
        wftrial(iorb,2,1)=evec(iorb,pm_con(2),ik_gamma)
      end if 
    end if
  end do
  ! sinse we are at gamme WF is the same at all UCells (we will compy values from wftrial(:,*,1) to everywhere)
  call tbmodel%bloch_wf_transform(kgrid,ik_gamma,sym,isym2,wftrial(:,3,1),wftrial(:,1,1))
  call tbmodel%bloch_wf_transform(kgrid,ik_gamma,sym,isym2,wftrial(:,4,1),wftrial(:,2,1))
  if (trim(adjustl(pars%wannier_proj_mode)).eq.'tbg12band') then
    wftrial(:,5,1)=evec(:,i1-4,ik_gamma)
    wftrial(:,6,1)= evec(:,i1-3,ik_gamma)
    wftrial(:,7,1)= evec(:,i1-2,ik_gamma)
    wftrial(:,8,1)= evec(:,i1-1,ik_gamma)
    wftrial(:,9,1)= evec(:,i1+4,ik_gamma)
    wftrial(:,10,1)=evec(:,i1+5,ik_gamma)
    wftrial(:,11,1)=evec(:,i1+6,ik_gamma)
    wftrial(:,12,1)=evec(:,i1+7,ik_gamma)
    wftrial(:,5,1)= wftrial(:,1,1)
    wftrial(:,6,1)= wftrial(:,2,1)
    wftrial(:,7,1)= wftrial(:,3,1) 
    wftrial(:,8,1)= wftrial(:,4,1) 
    wftrial(:,9,1)= wftrial(:,1,1)
    wftrial(:,10,1)=wftrial(:,2,1)
    wftrial(:,11,1)=wftrial(:,1,1)
    wftrial(:,12,1)=wftrial(:,2,1)
  end if
  do iR=2,tbmodel%rgrid%npt
    wftrial(:,:,iR)=wftrial(:,:,iR-1)
  end do
  sigma=0.4_dp*sqrt(dot_product(pars%avec(1,:),pars%avec(1,:)))
  if (trim(adjustl(pars%wannier_proj_mode)).eq.'tbg12band') then
    sigma=0.4_dp*sqrt(dot_product(pars%avec(1,:),pars%avec(1,:)))
  end if
  do iR=1,tbmodel%rgrid%npt
    do iw=1,proj%norb
      do iorb=1,tbmodel%norb_TB
        dv=tbmodel%vplorb(iorb)+tbmodel%rgrid%vpl(iR)-pars%proj%centers(:,pars%proj%iw2ic(iw))
        dc=matmul(dv,pars%avec)
        dd=sqrt(dot_product(dc,dc))
        t1=gauss(dd,sigma)
        wftrial(iorb,iw,iR)=wftrial(iorb,iw,iR)*t1
      end do
    end do
  end do
  if (mp_mpi) write(*,*) "wannier_interface%projection","overlap between trial orbitals"
  do iw=1,proj%norb
    do jw=iw,proj%norb
      zz=0._dp
      do iR=1,tbmodel%rgrid%npt
        zz=zz+dot_product(wftrial(:,iw,iR),wftrial(:,jw,iR))
      end do
      if (mp_mpi) write(*,'(2I5,2F10.6)') iw,jw,zz
    end do
  end do
  call generate_amn_overlap(tbmodel,pars,kgrid,evec,tbmodel%rgrid%npt,wftrial)
  nr=0
  do iR=1,tbmodel%rgrid%npt
    if (sum(abs(tbmodel%rgrid%vpl(iR))).le.6) then
      nr=nr+1
      do iw=1,proj%norb
        do iorb=1,tbmodel%norb_TB
          wftrial(iorb,iw,nr)=wftrial(iorb,iw,iR)
        end do
      end do
    end if
  end do
  !if (mp_mpi) call write_wf_universe(tbmodel,pars,nr,wftrial(:,:,1:nr),'wftrial','')
  call generate_dmn_orb(tbmodel,proj,sym,pars,kgrid,evec,.false.,wws)
  deallocate(wftrial)
else if (trim(adjustl(pars%wannier_proj_mode)).eq.'input_file') then
   ! find the home unit cel
   allocate(wftrial(tbmodel%norb_TB,pars%proj%norb,tbmodel%rgrid%npt))
   wftrial(:,:,:)=0._dp
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
         wftrial(jorb,ipro,iR_zero)=tbmodel%wbase%wws_full(sym%car(:,:,1),l1,m1,l2,m2,x1,z1,x2,z2)
       end do
     end do
   end do
   call generate_amn_overlap(tbmodel,pars,kgrid,evec,tbmodel%rgrid%npt,wftrial)
   call generate_dmn_orb(tbmodel,proj,sym,pars,kgrid,evec,.false.,wws)
   deallocate(wftrial)
else if (trim(adjustl(pars%wannier_proj_mode)).eq.'wannier_file') then
   ! find the home unit cel
   allocate(wftrial(tbmodel%norb_TB,pars%proj%norb,tbmodel%rgrid%npt))
   wftrial(:,:,:)=0._dp
   call read_wfmloc(pars,tbmodel,kgrid,evec,wftrial)
   call generate_amn_overlap(tbmodel,pars,kgrid,evec,tbmodel%rgrid%npt,wftrial)
   call generate_dmn_orb(tbmodel,proj,sym,pars,kgrid,evec,.false.,wws)
   nr=0
   do iR=1,tbmodel%rgrid%npt
     if (sum(abs(tbmodel%rgrid%vpl(iR))).le.6) then
       nr=nr+1
       do iw=1,proj%norb
         do iorb=1,tbmodel%norb_TB
           wftrial(iorb,iw,nr)=wftrial(iorb,iw,iR)
         end do
       end do
     end if
   end do
   !if (mp_mpi) call write_wf_universe(tbmodel,pars,nr,wftrial(:,:,1:nr),'wftrial','')
   deallocate(wftrial)
else if (trim(adjustl(pars%wannier_proj_mode)).eq.'real_space') then
   ! find the home unit cel
   allocate(wftrial(tbmodel%norb_TB,pars%proj%norb,tbmodel%rgrid%npt))
   wftrial(:,:,:)=0._dp
   sigma=2.0_dp*sqrt(dot_product(pars%avec(1,:),pars%avec(1,:)))
   call real_space_wftrial(tbmodel,proj,wftrial,sigma)
   call generate_amn_overlap(tbmodel,pars,kgrid,evec,tbmodel%rgrid%npt,wftrial)
   nr=0
   do iR=1,tbmodel%rgrid%npt
     if (sum(abs(tbmodel%rgrid%vpl(iR))).le.6) then
       nr=nr+1
       do iw=1,proj%norb
         do iorb=1,tbmodel%norb_TB
           wftrial(iorb,iw,nr)=wftrial(iorb,iw,iR)
         end do
       end do
     end if
   end do
   if (mp_mpi) call write_wf_universe(tbmodel,pars,nr,wftrial(:,:,1:nr),'wftrial','')
   call generate_dmn_orb(tbmodel,proj,sym,pars,kgrid,evec,.false.,wws)
   deallocate(wftrial)
else
  call throw("wannier_interface%generate_trial_wavefunctions()","unknown projection option")
end if
call generate_mmn_overlap(THIS,tbmodel,pars,kgrid,evec)
deallocate(wws)
end subroutine

subroutine generate_dmn_orb(tbmodel,base,sym,pars,kgrid,evec,lwws,wws)
class(CLtb), intent(in) :: tbmodel
class(CLsym), intent(in) :: sym
class(wbase), intent(in) :: base
class(CLpars), intent(in) :: pars
class(GRID), intent(inout) :: kgrid
complex(dp), intent(in) :: evec(tbmodel%norb_TB,pars%nstates,kgrid%npt)
logical, intent(in) :: lwws
complex(dp), intent(inout) :: wws(base%norb,base%norb,sym%nsym)
! local
integer isym,iw,jw,ic,jc,ir,ik
integer ikp,ist,jst,nn
logical exs
real(dp) err,t1,t2
real(dp) v1(NDIM),v2(NDIM)
complex(dp), allocatable :: phs(:,:)
complex(dp), allocatable :: wf_t(:)
complex(dp), allocatable :: ovlp(:,:,:,:)
inquire(file=trim(adjustl(pars%seedname))//'.dmn',exist=exs)
if (exs) then
  call info("CLwan%generate_dmn_orb","skipping "//trim(adjustl(pars%seedname))//".dmn creation")
  return
end if
if (.not.lwws) then
  wws(:,:,:)=0._dp
  do isym=1,sym%nsym
    do iw=1,base%norb
       ic=base%orb_icio(iw,1)
       jc=base%ics2c(ic,isym)
       do jw=1,base%norb
          if(base%orb_icio(jw,1).ne.jc) cycle
          wws(jw,iw,isym)=base%wws(sym%car(:,:,isym),iw,jw)
       end do
    end do
    do iw=1,base%norb
       err=abs((sum(wws(:,iw,isym)**2)+sum(wws(iw,:,isym)**2))*.5_dp-1._dp)
       if(err.gt.1.e-3_dp) then
          if (mp_mpi) write(*,*) "wannier_interface%generate_dmn_orb","compute_dmn: Symmetry operator (",isym, &
                  ") could not transform Wannier function (",iw,")."
          if (mp_mpi) write(*,*) "The error is ",err,"."
          call throw("wannier_interface%generate_dmn_orb", "missing Wannier functions, see the output.")
       end if
    end do
  end do
end if
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
  do isym=1,sym%nsym
     do iw=1,base%norb
        ic=base%orb_icio(iw,1)
        jc=base%ics2c(ic,sym%inv(isym))
        v1=kgrid%vpc(kgrid%iks2k(ik,isym))-matmul(sym%car(:,:,isym),kgrid%vpc(ik))
        v2=matmul(v1,sym%car(:,:,isym))
        t1=dot_product(base%vcs2t(:,jc,isym),kgrid%vpc(ik))
        t2=dot_product(sym%vtc(:,isym),v2)
        phs(iw,iw)=cmplx(cos(t1),sin(t1),kind=dp)*cmplx(cos(t2),sin(t2),kind=dp)
     end do
     if (mp_mpi) then
        WRITE (1001,*)
        WRITE (1001,"(1p,(' (',e18.10,',',e18.10,')'))") matmul(phs,wws(:,:,isym))
     end if
  end do
end do
if (mp_mpi) then
  write(1001,*)
  flush(1001)
  write(*,'(/)')
  write(*,'(a,i8)') '  DMN(d_matrix_band) [could be not ordered]: nir = ',kgrid%nir
end if
allocate(ovlp(pars%nstates,pars%nstates,sym%nsym,kgrid%nir))
ovlp=0._dp
#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
!$OMP PARALLEL DEFAULT(SHARED)&
!$OMP PRIVATE(ist,jst,ik,ikp,isym,wf_t)
!$OMP DO
do ir=1,kgrid%nir
  if (mod(ir-1,np_mpi).ne.lp_mpi) cycle
  WRITE (*,'(i8)',advance='no') ir
  IF( MOD(ir,10) == 0 ) WRITE (*,*)
  allocate(wf_t(tbmodel%norb_TB))
  ik=kgrid%ir2ik(ir)
  do isym=1,sym%nsym
    ikp=kgrid%iks2k(ik,isym)
    do jst=1,pars%nstates
      call tbmodel%bloch_wf_transform(kgrid,ik,sym,isym,wf_t,evec(:,jst,ik))
      do ist=1,pars%nstates
         ovlp(ist,jst,isym,ir)=dot_product(evec(:,ist,ikp),wf_t)
      end do
    end do
  end do
  deallocate(wf_t)
end do
!$OMP END DO
!$OMP END PARALLEL
#ifdef MPI
  nn=pars%nstates*pars%nstates*sym%nsym*kgrid%nir
  call mpi_allreduce(mpi_in_place,ovlp,nn,mpi_double_complex,mpi_sum, &
   mpi_com,mpi_err)
#endif
if (mp_mpi) then
  do ir=1,kgrid%nir
    do isym=1,sym%nsym
      do jst=1,pars%nstates
        do ist=1,pars%nstates
          write(1001,"(1p,(' (',e18.10,',',e18.10,')'))") ovlp(ist,jst,isym,ir)
        end do
      end do
      write(1001,*) 
    end do
  end do
  close(1001)
  write(*,*)
end if
deallocate(phs,ovlp)
#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
end subroutine

function eps_pm_dublet_state(tbmodel,pars,kgrid,ik_gamma,sym,isym,ist1,ist2,evec_G) result (pm)
class(CLtb), intent(in) :: tbmodel
class(CLpars), intent(in) :: pars
class(GRID), intent(in) :: kgrid
class(CLsym), intent(in) :: sym
integer, intent(in) :: isym,ist1,ist2,ik_gamma
complex(dp), intent(in) :: evec_G(tbmodel%norb_TB,pars%nstates)
integer pm(2)
real(dp) twopi3,a1,a2,b1,b2
complex(dp), allocatable :: wf_t(:)
! constant
twopi3=twopi/3._dp
! find one of the lower dublet state (E-) with eigenvalus +\epsilon
allocate(wf_t(tbmodel%norb_TB))
!wf_t=evec_G(:,ist1)
call tbmodel%bloch_wf_transform(kgrid,ik_gamma,sym,isym,wf_t,evec_G(:,ist1))
a1=dble(dot_product(evec_G(:,ist1),wf_t))
b1=dble(dot_product(evec_G(:,ist2),wf_t))
!wf_t=evec_G(:,ist2)
call tbmodel%bloch_wf_transform(kgrid,ik_gamma,sym,isym,wf_t,evec_G(:,ist2))
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
    pm(1)=ist1
    pm(2)=ist2
else if ( ( abs(a2-cos(twopi3)).lt.epsr .and. abs(b2-sin(twopi3)).lt.epsr ) .or.&
       ( abs(b2-cos(twopi3)).lt.epsr .and. abs(a2-sin(twopi3)).lt.epsr ) ) then
    pm(1)=ist2
    pm(2)=ist1
else
  call throw("CLwan%project","could not find state with 2pi/3 rotation eigenvalue")
end if
deallocate(wf_t)
end function

subroutine generate_amn_overlap(tbmodel,pars,kgrid,evec,nr,wftrial)
class(CLtb), intent(in) :: tbmodel
class(CLpars), intent(in) :: pars
class(GRID), intent(in) :: kgrid
complex(dp), intent(in) :: evec(tbmodel%norb_TB,pars%nstates,kgrid%npt)
integer, intent(in) :: nr
complex(dp), intent(in) :: wftrial(tbmodel%norb_TB,pars%proj%norb,nr)
! local
logical exs
integer ik,iwan,ist,iR,iorb
real(dp) t1
complex(dp) z1
complex(dp), allocatable :: amn(:,:,:)
! A_mn(k)=<\psi_m(k)|g_n>
inquire(file=trim(adjustl(pars%seedname))//'.amn',exist=exs)
if (exs) then
  call info("CLwan%generate_amn_overlap","skipping "//trim(adjustl(pars%seedname))//".amn creation")
  return
else
  call info("CLwan%generate_amn_overlap","generating "//trim(adjustl(pars%seedname))//".amn file")
end if
if (mp_mpi) then
  open(50,file=trim(adjustl(pars%seedname))//'.amn',action='write')
  write(50,*) '# '//trim(adjustl(pars%seedname))//' file '
  write(50,*) pars%nstates,kgrid%npt,pars%proj%norb
end if
allocate(amn(pars%nstates,pars%proj%norb,kgrid%npt))
amn=0._dp
!$OMP PARALLEL DEFAULT (SHARED)&
!$OMP PRIVATE(iR,iwan,iorb,t1,z1)
!$OMP DO
do ik=1,kgrid%npt
  do iR=1,tbmodel%rgrid%npt
    if ( sum(abs(tbmodel%rgrid%vpl(iR))).gt.6) cycle
    t1=dot_product(kgrid%vpl(ik),tbmodel%rgrid%vpl(iR))*twopi
    z1=cmplx(cos(t1),-sin(t1),kind=dp)
    do iwan=1,pars%proj%norb
      do iorb=1,tbmodel%norb_TB
        amn(:,iwan,ik)=amn(:,iwan,ik)+conjg(evec(iorb,:,ik))*wftrial(iorb,iwan,iR)*z1
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP END PARALLEL
do ik=1,kgrid%npt
  do iwan=1,pars%proj%norb
    do ist=1,pars%nstates
      if (mp_mpi) write(50,'(3I6,2G18.10)') ist,iwan,ik,dble(amn(ist,iwan,ik)),aimag(amn(ist,iwan,ik))
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
integer ik,jk,mm,nn,iorb,innk
real(dp) dc,t1
complex(dp) z1
real(dp) vq(NDIM),vc(NDIM)
complex(dp), allocatable :: mmn(:,:,:)
! M_mn(k)=<u_mk|u_n{k+q}>
inquire(file=trim(adjustl(pars%seedname))//'.mmn',exist=exs)
if (exs) then
  call info("CLwan%generate_mmn_overlap","skipping "//trim(adjustl(pars%seedname))//".mmn creation")
  return
else
  call info("CLwan%generate_mmn_overlap","generating "//trim(adjustl(pars%seedname))//".mmn file")
end if
if (mp_mpi) then
  open(50,file=trim(adjustl(pars%seedname))//'.mmn',action='write')
  write(50,*) '# '//trim(adjustl(pars%seedname))//' file '
  write(50,*) pars%nstates,kgrid%npt,THIS%nnk
end if
allocate(mmn(pars%nstates,pars%nstates,THIS%nnk))
do ik=1,kgrid%npt
  mmn=0._dp
  !$OMP PARALLEL DEFAULT (SHARED)&
  !$OMP PRIVATE(innk,jk,iorb,nn,vq,vc,dc,t1,z1)
  !$OMP DO
  do innk=1,THIS%nnk
    jk=THIS%nnkp(4,innk,ik)
    vq=kgrid%vpl(jk)+dble(THIS%nnkp(1:3,innk,ik))-kgrid%vpl(ik)
    vc=matmul(vq,kgrid%vecs)
    dc=sqrt(dot_product(vc,vc))
    do iorb=1,tbmodel%norb_TB
      t1=dot_product(vq,tbmodel%vplorb(iorb))*twopi
      z1=cmplx(cos(t1),-sin(t1),kind=dp)
      do nn=1,pars%nstates
        mmn(:,nn,innk)=mmn(:,nn,innk)+conjg(evec(iorb,:,ik))*evec(iorb,nn,jk)*z1*pwave_ovlp(dc)
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  do innk=1,THIS%nnk
    if (mp_mpi)  write(50,'(5I6)') ik,THIS%nnkp(4,innk,ik),THIS%nnkp(1:3,innk,ik)
    do nn=1,pars%nstates
      do mm=1,pars%nstates
        if (mp_mpi) write(50,'(2G18.10)') dble(mmn(mm,nn,innk)),aimag(mmn(mm,nn,innk))
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
write(50,*) '!site_symmetry = .true.'
write(50,*) '!symmetrize_eps=  1d-7'


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
do iwan=1,pars%proj%norb
  write(50,'("!",G18.10,",",G18.10,",",G18.10)') &
    matmul(pars%proj%centers(:,pars%proj%iw2ic(iwan)),pars%avec)
end do
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
  write(50,'(A2,3G18.10)')'XX',pars%proj%centers(:,pars%proj%iw2ic(iwan))
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

subroutine write_wf_universe(tbmodel,pars,nr,wf,pre,post)
class(CLtb), intent(in) :: tbmodel
class(CLpars), intent(in) :: pars
integer, intent(in) :: nr
complex(dp), intent(in) :: wf(tbmodel%norb_TB,pars%proj%norb,nr)
character(len=*),intent(in) :: pre
character(len=*),intent(in) :: post
character(len=128) fname
character(len=128) num
integer iorb,iR,iw,ic,jr
do iw=1,pars%proj%norb
  write(num,'(I2)') iw
  fname=trim(adjustl(pre))//trim(adjustl(num))//trim(adjustl(post))//'.dat'
  open(50,file=trim(adjustl(fname)))
  write(50,*) pars%avec(1,:)
  write(50,*) pars%avec(2,:)
  write(50,*) pars%avec(3,:)
  write(50,*) 'atoms ',tbmodel%norb_TB
  write(50,*) 'nrpt_wan ',nr
  jr=0
  do iR=1,tbmodel%rgrid%npt
    if (sum(abs(tbmodel%rgrid%vpl(iR))).le.6) then
      jr=jr+1
      write(50,*) tbmodel%rgrid%vpi(iR)
      do iorb=1,tbmodel%norb_TB
        ic=tbmodel%wbase%orb_icio(iorb,1)
        write(50,'(5G18.10)') tbmodel%wbase%centers_cart(:,ic),wf(iorb,iw,jr)
      end do
    end if
  end do
  close(50)
end do
return
end subroutine


subroutine read_wfmloc(pars,tbmodel,kgrid,evec,wfmloc)
class(CLpars), intent(in) :: pars
class(CLtb), intent(in) :: tbmodel
class(GRID), intent(inout) :: kgrid
complex(dp), intent(in) :: evec(tbmodel%norb_TB,pars%nstates,kgrid%npt)
complex(dp), intent(out) :: wfmloc(tbmodel%norb_TB,pars%proj%norb,tbmodel%rgrid%npt)
logical exs
integer nk_,nwan_,num_bands_
integer ik,m,n,ir
real(dp) t1,t2,vpl_(3)
complex(dp) z1
complex(dp), allocatable :: umat(:,:,:),udis(:,:,:)
complex(dp), allocatable :: psi_dis(:,:)
type(wbase)  proj
if (pars%proj%norb.le.0) then
   call throw("wannier_interface%read_wfmloc()",&
              "apparently 'projections' block was not specified, wannier projections not found")
else
  if (pars%nstates.lt.pars%proj%norb) then
    call throw("wannier_interface%read_wfmloc()",&
               "number of bands is less than the number of requested projections")
  end if
end if
call proj%init(pars,pars%proj%ncenters,pars%proj%norb,pars%proj%norb_ic,&
                   pars%proj%lmr,pars%proj%waxis,pars%proj%centers)
inquire(file=trim(adjustl(pars%seedname))//'_u.in',exist=exs)
if (.not.exs) then
  call throw("Wannier_supplementary%read_wfmloc","no "//trim(adjustl(pars%seedname))//'_u.in'//" file")
  stop
end if
call info("Wannier_supplementary%read_wfmloc","constructing wfmloc")
allocate(udis(pars%nstates,proj%norb,kgrid%npt))
allocate(umat(proj%norb,proj%norb,kgrid%npt))
allocate(psi_dis(tbmodel%norb_TB,proj%norb))
if (pars%nstates.gt.proj%norb) then
  inquire(file=trim(adjustl(pars%seedname))//'_u_dis.in',exist=exs)
  if (.not.exs) then
    call throw("Wannier_supplementary%read_wfmloc","no "//trim(adjustl(pars%seedname))//'_u_dis.in'//" file")
    stop
  end if
  open(50,file=trim(adjustl(pars%seedname))//'_u_dis.in',action='read')
  read(50,*)
  read(50,*) nk_,nwan_,num_bands_

  if (nk_.ne.kgrid%npt)  call throw("wannier_interface%read_wfmloc()",&
               "number of k-points is different from one derived from input")
  if (nwan_.ne.proj%norb)  call throw("wannier_interface%read_wfmloc()",&
               "number projections in wannier files is different from one derived from input")
  if (num_bands_.ne.pars%nstates)  call throw("wannier_interface%read_wfmloc()",&
               "number bands in wannier files is different from one derived from input")

  do ik=1,nk_
    read(50,*)
    read(50,*) vpl_
    if (sum(abs(vpl_-kgrid%vpl(ik))).gt.epslat) then
       write(*,*) ik
       write(*,*) vpl_
       write(*,*) kgrid%vpl(ik)
       call throw("wannier_interface%read_wfmloc()","different k-vectors")
    end if 
    do m=1,nwan_
      do n=1,num_bands_
        read(50,*) t1,t2
        udis(n,m,ik)=cmplx(t1,t2,kind=dp)
      end do
    end do
  end do
  close(50)
end if
open(51,file=trim(adjustl(pars%seedname))//'_u.mat',action='read')
read(51,*)
read(51,*) nk_,nwan_,n

if (nk_.ne.kgrid%npt)  call throw("wannier_interface%read_wfmloc()",&
             "number of k-points is different from one derived from input")
if (nwan_.ne.proj%norb)  call throw("wannier_interface%read_wfmloc()",&
             "number projections in wannier files is different from one derived from input")
do ik=1,nk_
  read(51,*)
  read(51,*) vpl_
  if (sum(abs(vpl_-kgrid%vpl(ik))).gt.epslat) then
     write(*,*) ik
     write(*,*) vpl_
     write(*,*) kgrid%vpl(ik)
     call throw("wannier_interface%read_wfmloc()","different k-vectors")
  end if 
  do n=1,nwan_
    do m=1,nwan_
      read(51,*) t1,t2
      umat(m,n,ik)=cmplx(t1,t2,kind=dp)
    end do
  end do
end do
close(51)
wfmloc(:,:,:)=0.d0
#ifdef MPI
  call mpi_barrier(mpi_com,mpi_err)
#endif
do ik=1,kgrid%npt
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  if (pars%nstates.gt.proj%norb) then
    psi_dis(:,:)=0._dp
    do n=1,proj%norb
      do m=1,pars%nstates
        psi_dis(:,n)=psi_dis(:,n)+udis(m,n,ik)*evec(:,m,ik)
      end do
    end do
  else if (pars%nstates.eq.proj%norb) then
    psi_dis=evec(:,:,ik)
  else
    call throw("Wannier_supplementary%read_wfmloc"," num_bands<nwan")
  end if
!$OMP PARALLEL DEFAULT (SHARED)&
!$OMP PRIVATE(n,m,t1,z1)
!$OMP DO
  do ir=1,tbmodel%rgrid%npt
    t1=dot_product(kgrid%vpc(ik),tbmodel%rgrid%vpc(ir))
    z1=cmplx(cos(t1),sin(t1),kind=dp)
    do n=1,proj%norb
      do m=1,proj%norb
        wfmloc(:,n,ir)=wfmloc(:,n,ir)+psi_dis(:,m)*umat(m,n,ik)*z1
      end do
    end do
  end do
!$OMP END DO
!$OMP END PARALLEL
end do
#ifdef MPI
  n=tbmodel%rgrid%npt*proj%norb*tbmodel%norb_TB
  call mpi_allreduce(mpi_in_place,wfmloc,n,mpi_double_complex,mpi_sum, &
   mpi_com,mpi_err)
#endif
wfmloc=wfmloc/dble(kgrid%npt)
deallocate(umat,udis,psi_dis)
return
end subroutine


subroutine real_space_wftrial(tbmodel,proj,wftrial,sigma)
class(CLtb), intent(in) :: tbmodel
class(wbase), intent(in) :: proj
complex(dp), intent(out) :: wftrial(tbmodel%norb_TB,proj%norb,tbmodel%rgrid%npt)
real(dp), intent(in) :: sigma
integer pr,qr,kr,nr,iw,iR,iorb
integer l1,l2,m1,m2,jr,jc,ic
real(dp) vc(NDIM),t1
real(dp) vaxis1(NDIM,NDIM),vaxis2(NDIM,NDIM)
real(dp), allocatable :: rr(:,:)
real(dp), allocatable :: r0(:,:)
real(dp), allocatable :: ylm1(:),ylm2(:)
integer, parameter :: nbox=4
real(dp), parameter :: radius=3._dp
nr=(2*nbox)**NDIM
allocate(r0(NDIM,nr))
nr=0
do pr=-nbox,nbox
  if (pr.eq.0) cycle
  do qr=-nbox,nbox
    if (qr.eq.0) cycle
    do kr=-nbox,nbox
      if (kr.eq.0) cycle
      nr=nr+1
      r0(1,nr)=dble(pr)*radius
      r0(2,nr)=dble(qr)*radius
      r0(3,nr)=dble(kr)*radius
 !     write(*,*) r0(:,nr)
    end do
  end do
end do
wftrial=0._dp
!$OMP PARALLEL DEFAULT(SHARED)&
!$OMP PRIVATE(l1,l2,m1,m2,ic,jc,iorb,ir,jr,t1,vc)&
!$OMP PRIVATE(vaxis1,vaxis2,rr,ylm1,ylm2)
!$OMP DO
do iw=1,proj%norb
  allocate(rr(NDIM,nr))
  allocate(ylm1(nr))
  allocate(ylm2(nr))
  l1=proj%lmr(1,iw)
  m1=proj%lmr(2,iw)
  ic=proj%orb_icio(iw,1)
  call set_u_matrix(proj%waxis(:,1,iw),proj%waxis(:,2,iw),vaxis1)
  do iR=1,tbmodel%rgrid%npt
    if ( sum(abs(tbmodel%rgrid%vpl(iR))).gt.6) cycle
    do iorb=1,tbmodel%wbase%norb
      l2=tbmodel%wbase%lmr(1,iorb)
      m2=tbmodel%wbase%lmr(2,iorb)
      jc=tbmodel%wbase%orb_icio(iorb,1)
      ! rotation operator
      call set_u_matrix(tbmodel%wbase%waxis(:,1,iorb),tbmodel%wbase%waxis(:,2,iorb),vaxis2)
      ! angular dependes of TB basis orbial
      CALL ylm_wannier(ylm2,l2,m2,matmul(vaxis2,r0),nr)
      ! its radial decay
      do jr=1,nr
        t1=sqrt(dot_product(r0(:,jr),r0(:,jr)))
        ylm2(jr)=ylm2(jr)*raddec(1._dp,t1)
      end do
      ! center of basis orbital
      vc=tbmodel%wbase%centers_cart(:,jc)+tbmodel%rgrid%vpc(iR)
      do jr=1,nr
        ! positions with respect to projection center
        rr(:,jr)=r0(:,jr)+vc(:)-proj%centers_cart(:,ic)
        if (sum(abs(rr(:,jr))).lt.epslat) rr(:,jr)=epslat
      end do
      ! angular dependce of projection function
      CALL ylm_wannier(ylm1,l1,m1,matmul(vaxis1,rr),nr)
      ! its radial decay
      do jr=1,nr
        t1=sqrt(dot_product(rr(:,jr),rr(:,jr)))
        ylm1(jr)=ylm1(jr)*gauss(t1,sigma)
      end do
      wftrial(iorb,iw,iR)=dot_product(ylm1,ylm2)
    end do
  end do
  t1=0._dp
  do iR=1,tbmodel%rgrid%npt
     t1=t1+abs(dot_product(wftrial(:,iw,iR),wftrial(:,iw,iR)))
  end do
  wftrial(:,iw,:)=wftrial(:,iw,:)/sqrt(t1)
  deallocate(rr,ylm1,ylm2)
end do
!$OMP END DO
!$OMP END PARALLEL
deallocate(r0)
end subroutine

end module
