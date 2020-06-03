
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
public :: compute_hubbardu_rs
public :: compute_hubbardj_rs
public :: symmetrize_hubbardu_rs
public :: symmetrize_tbfile

type, public :: CLwan
  integer :: nnk=0
  integer :: dis_num_iter=1000
  integer :: num_iter=1000
  real(dp) :: dis_win_min=-100._dp
  real(dp) :: dis_win_max= 100._dp
  real(dp) :: dis_froz_min=-100._dp
  real(dp) :: dis_froz_max= 100._dp
  integer, allocatable :: nnkp(:,:,:)
  real(dp), allocatable :: dege(:)
  complex(dp), allocatable :: hame(:,:,:)
  type(GRID) :: rgrid
  contains
  procedure :: init
  procedure :: projection
  procedure :: readnnkp
  procedure, private :: read_tb_file
  procedure, private :: write_tb_file
  procedure, private :: writewin
endtype CLwan

type, private :: coulrs
  integer :: npoly=25
  real(dp) :: rcut_bare
  real(dp) :: scale_bare
  real(dp), allocatable :: polyc(:)
  real(dp), allocatable :: xpoly(:)
  contains 
  procedure :: init=>init_coulrs
  procedure :: evaluate=>evaluate_coulrs
  procedure, nopass :: compute_poly
endtype

contains

real(dp) function compute_poly(np,xa,c,x)
integer, intent(in) :: np
real(dp), intent(in) :: xa(np),c(np),x
integer i
real(dp) t1,sum
sum=c(1)
t1=1.d0
do i=2,np
  t1=t1*(x-xa(i-1))
  sum=sum+c(i)*t1
end do
compute_poly=sum
end function
subroutine init_coulrs(THIS,fname)
class(coulrs), intent(inout) :: THIS
character(len=*), intent(in) :: fname
integer, parameter :: maxpoints=10000
logical exs
integer il,ix,nx,iline,counter
real(dp) a,b,t1,t2,t3
real(dp), allocatable :: xx(:),yy(:)
real(dp) yl(THIS%npoly)
inquire(file=trim(adjustl(fname)),exist=exs)
if (.not.exs) call throw("coulrs%init()","file '"//trim(adjustl(fname))//"' with Coulomb interaction not found")
call info("coulrs%init()","assuming meV units in the '"//trim(adjustl(fname))//"' file")
allocate(xx(maxpoints),yy(maxpoints))
open(44,file=trim(adjustl(fname)))
nx=0
do iline=1,maxpoints
  read(44,*,end=100) a,b
  if (a.lt.epslat) cycle
  nx=nx+1
  xx(nx)=a
  yy(nx)=b*0.001_dp
end do
call throw("coulrs%init()","maxpoints exceeded")
100 continue 
close(44)
allocate(THIS%xpoly(THIS%npoly))
allocate(THIS%polyc(THIS%npoly))
! take x from xx and corresponding y from log-like distribution
THIS%rcut_bare=xx(nx)
THIS%scale_bare=yy(nx)*xx(nx)
t1=1.d0/dble(THIS%npoly-1)
t2=log(THIS%rcut_bare/epslat)
counter=1
do il=1,THIS%npoly
  t3=epslat*exp(dble(il-1)*t1*t2)
  if (t3.lt.xx(counter)) then
    THIS%xpoly(il)=xx(counter)
    yl(il)=yy(counter)
    counter=counter+1
  else
    do ix=1,nx
      if (xx(ix).ge.t3) exit
    end do
    THIS%xpoly(il)=xx(ix)
    yl(il)=yy(ix)
  end if
end do
THIS%xpoly=1._dp/THIS%xpoly
t1=polynm(0,THIS%npoly,THIS%xpoly,yl,1._dp,THIS%polyc)
!do ix=1,nx
!  write(156,*) xx(ix),yy(ix)
!  write(157,*) 10._dp*xx(ix),THIS%evaluate(10._dp*xx(ix))
!  write(158,*) 10._dp*xx(ix),CoulombForceConstant/(10._dp*xx(ix))
!end do
!stop
deallocate(xx,yy)
return
end subroutine
real(dp) function evaluate_coulrs(THIS,rij)
class(coulrs), intent(in) :: THIS
real(dp), intent(in) :: rij
real(dp), parameter :: eps=1.e-8_dp
if (abs(rij).lt.eps) then
  evaluate_coulrs=0._dp
else if (rij.lt.THIS%rcut_bare) then
  evaluate_coulrs=THIS%compute_poly(THIS%npoly,THIS%xpoly,THIS%polyc,1._dp/rij)
else
  evaluate_coulrs=THIS%scale_bare/rij
end if
end function

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
integer isym1,isym2,ik_gamma,ikg(NDIM+1)
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
         wftrial(jorb,ipro,tbmodel%rgrid%ip0)=tbmodel%wbase%wws_full(sym%car(:,:,1),l1,m1,l2,m2,x1,z1,x2,z2)
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
  allocate(wf_t(tbmodel%norb_TB))
!$OMP DO
do ir=1,kgrid%nir
  !$OMP CRITICAL
  WRITE (*,*) "ir: ",ir
  !$OMP END CRITICAL
  if (mod(ir-1,np_mpi).ne.lp_mpi) cycle
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
end do
!$OMP END DO
  deallocate(wf_t)
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
  !call info("CLwan%generate_amn_overlap","skipping "//trim(adjustl(pars%seedname))//".amn creation")
  !return
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
 ! write(*,*) ik
  do iR=1,tbmodel%rgrid%npt
    if ( sum(abs(wftrial(:,:,iR))).lt.epslat) cycle 
 !   write(*,*) iR,tbmodel%rgrid%vpl(iR)
 !   write(*,*) wftrial(:,1,iR)
    t1=dot_product(kgrid%vpl(ik),tbmodel%rgrid%vpl(iR))*twopi
    z1=cmplx(cos(t1),-sin(t1),kind=dp)
    do iwan=1,pars%proj%norb
      do iorb=1,tbmodel%norb_TB
        amn(:,iwan,ik)=amn(:,iwan,ik)+conjg(evec(iorb,:,ik))*wftrial(iorb,iwan,iR)*z1
      end do
    end do
  end do
!  write(*,*)"amn, ",amn(:,1,ik)
!  write(*,*)"evec, ",evec(:,1,ik)
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
  !call info("CLwan%generate_mmn_overlap","skipping "//trim(adjustl(pars%seedname))//".mmn creation")
  !return
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
  allocate(rr(NDIM,nr))
  allocate(ylm1(nr))
  allocate(ylm2(nr))
!$OMP DO
do iw=1,proj%norb
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
end do
!$OMP END DO
  deallocate(rr,ylm1,ylm2)
!$OMP END PARALLEL
deallocate(r0)
end subroutine

subroutine read_wfmloc(pars,tbmodel,kgrid,evec,wfmloc)
class(CLpars), intent(in) :: pars
class(CLtb), intent(in) :: tbmodel
class(GRID), intent(inout) :: kgrid
complex(dp), intent(in) :: evec(tbmodel%norb_TB,pars%nstates,kgrid%npt)
complex(dp), intent(out) :: wfmloc(tbmodel%norb_TB,pars%proj%norb,tbmodel%rgrid%npt)
integer ik,ir
real(dp), allocatable :: vkl(:,:),vrl(:,:)
complex(dp), allocatable :: umat(:,:,:),udis(:,:,:)
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
call info("Wannier_interface%read_wfmloc","constructing wfmloc")
allocate(udis(pars%nstates,proj%norb,kgrid%npt))
allocate(umat(proj%norb,proj%norb,kgrid%npt))
! k-point list
allocate(vkl(NDIM,kgrid%npt))
do ik=1,kgrid%npt
  vkl(:,ik)=kgrid%vpl(ik)
end do
! r-point list
allocate(vrl(NDIM,tbmodel%rgrid%npt))
do ir=1,tbmodel%rgrid%npt
  vrl(:,ir)=tbmodel%rgrid%vpl(ir)
end do
if (pars%nstates.gt.proj%norb) then
  call read_udis(pars%seedname,kgrid%npt,proj%norb,pars%nstates,vkl,udis)
end if
call read_umat(pars%seedname,kgrid%npt,proj%norb,vkl,umat)
call get_wfmloc(.false.,proj%norb,pars%nstates,kgrid%npt,tbmodel%rgrid%npt,&
   tbmodel%norb_TB,udis,umat,vkl,vrl,evec,wfmloc)
deallocate(vkl,vrl,umat,udis)
return
end subroutine

subroutine symmetrize_tbfile(pars,sym)
class(CLpars), intent(in) :: pars
class(CLsym), intent(in) :: sym
real(dp), allocatable :: wws(:,:,:)
real(dp), allocatable :: wwsi(:,:,:)
complex(dp), allocatable :: hams(:,:,:)
type(wbase)  proj
type(GRID) rgrid
type(CLwan) wan
real(dp) Tmlij(NDIM),tauji(NDIM),tauml(NDIM),RR(NDIM)
integer isym,iR,iw,jw
integer iorb,jorb,morb,lorb
integer ic,jc,mc,lc
integer iv(NDIM+1)
if (pars%tbtype.ne.'tbfile') call throw("wannier_interface%symmetrize_tbfile()",&
              "only tbtype=tbfile is supported,i.e, onlt _tb.dat format of wannier90 is accepted")
if (pars%proj%norb.le.0) then
   call throw("wannier_interface%symmetrize_tbfile()",&
              "apparently 'projections' block was not specified, wannier projections not found")
else
  if (pars%nstates.lt.pars%proj%norb) then
    call throw("wannier_interface%symmetrize_tbfile()",&
               "number of bands is less than the number of requested projections")
  end if
end if
call proj%init(pars,pars%proj%ncenters,pars%proj%norb,pars%proj%norb_ic,&
                   pars%proj%lmr,pars%proj%waxis,pars%proj%centers)
call proj%init_smap(sym,pars)
call rgrid%init(pars%ngrid,pars%avec,.true.,.false.)
call wan%read_tb_file(pars,proj%norb)
allocate(wws(proj%norb,proj%norb,sym%nsym))
allocate(wwsi(proj%norb,proj%norb,sym%nsym))
do isym=1,sym%nsym
  do iw=1,proj%norb
     ic=proj%orb_icio(iw,1)
     jc=proj%ics2c(ic,isym)
     do jw=1,proj%norb
        if(proj%orb_icio(jw,1).ne.jc) cycle
        wws(jw,iw,isym)=proj%wws(sym%car(:,:,isym),iw,jw)
        wwsi(jw,iw,isym)=proj%wws(sym%car(:,:,sym%inv(isym)),iw,jw)
     end do
  end do
end do
allocate(hams(proj%norb,proj%norb,wan%rgrid%npt))
hams=0._dp
do isym=1,sym%nsym
  do iR=1,wan%rgrid%npt
    do iorb=1,proj%norb
      ic=proj%orb_icio(iorb,1)
      do jorb=1,proj%norb
        jc=proj%orb_icio(jorb,1)
        tauji=proj%centers(:,pars%proj%iw2ic(jorb))-proj%centers(:,pars%proj%iw2ic(iorb))
        do lorb=1,proj%norb
          lc=proj%orb_icio(lorb,1)
          if (proj%ics2c(lc,isym).ne.ic) cycle
          do morb=1,proj%norb
            mc=proj%orb_icio(morb,1)
            if (proj%ics2c(mc,isym).ne.jc) cycle
            tauml=proj%centers(:,pars%proj%iw2ic(morb))-proj%centers(:,pars%proj%iw2ic(lorb))
            Tmlij=proj%Tmlij(dble(sym%lat(:,:,isym)),tauml,tauji)
            RR=matmul(dble(sym%lat(:,:,sym%inv(isym))),wan%rgrid%vpl(iR)-Tmlij)
            iv=wan%rgrid%find(RR)
            if (iv(NDIM+1).le.0.or.sum(abs(iv(1:NDIM))).ne.0) cycle
            hams(iorb,jorb,iR)=hams(iorb,jorb,iR)+wws(iorb,lorb,isym)*wan%hame(lorb,morb,iv(NDIM+1))*wwsi(morb,jorb,isym)
          end do
        end do
      end do
    end do
  end do
end do
wan%hame=hams/dble(sym%nsym)
call wan%write_tb_file(pars,proj%norb)
end subroutine

subroutine symmetrize_hubbardu_rs(pars,sym,tbmodel,kgrid,evec)
class(CLpars), intent(in) :: pars
class(CLsym), intent(in) :: sym
class(CLtb), intent(in) :: tbmodel
class(GRID), intent(inout) :: kgrid
complex(dp), intent(in) :: evec(tbmodel%norb_TB,pars%nstates,kgrid%npt)
integer jR_sphere,norb,npt_sphere
integer nn,mm,pp,qq
integer nn_,mm_,pp_,qq_
integer isym,iR,iw,jw,ik
integer iorb,jorb,morb,lorb,porb,qorb
integer ic,jc,mc,lc,pc,qc
integer iv(NDIM+1)
real(dp) t1,t2
complex(dp) z1
real(dp) Tmlij(NDIM),RR(NDIM)
real(dp) tauji(NDIM),tauml(NDIM),taupq(NDIM)
real(dp), allocatable :: vkl(:,:)
real(dp), allocatable :: vrl(:,:)
real(dp), allocatable :: wws(:,:,:)
real(dp), allocatable :: wwsi(:,:,:)
complex(dp), allocatable :: HubU(:,:,:,:,:)
complex(dp), allocatable :: HubUn(:,:,:)
complex(dp), allocatable :: HubUs(:,:,:)
complex(dp), allocatable :: umat(:,:,:),udis(:,:,:)
complex(dp), allocatable :: wfmloc(:,:,:)
complex(dp), allocatable :: wf1(:,:)
complex(dp), allocatable :: wf2(:,:)
type(wbase)  proj
type(GRID) rgrid
complex(dp) zdotc
!real(dp) ddot
if (pars%proj%norb.le.0) then
   call throw("wannier_interface%symmetrize_hubbardu_rs()",&
              "apparently 'projections' block was not specified, wannier projections not found")
else
  if (pars%nstates.lt.pars%proj%norb) then
    call throw("wannier_interface%symetrize_hubbardu_rs()",&
               "number of bands is less than the number of requested projections")
  end if
end if
call proj%init(pars,pars%proj%ncenters,pars%proj%norb,pars%proj%norb_ic,&
                   pars%proj%lmr,pars%proj%waxis,pars%proj%centers)
call proj%init_smap(sym,pars)
call rgrid%init(pars%ngrid,pars%avec,.true.,.false.)
allocate(udis(pars%nstates,proj%norb,kgrid%npt))
allocate(umat(proj%norb,proj%norb,kgrid%npt))
allocate(wfmloc(tbmodel%norb_TB,proj%norb,rgrid%npt))
allocate(wf1(tbmodel%norb_TB,rgrid%npt))
allocate(wf2(tbmodel%norb_TB,rgrid%npt))
! k-point list
allocate(vkl(NDIM,kgrid%npt))
do ik=1,kgrid%npt
  vkl(:,ik)=kgrid%vpl(ik)
end do
! r-point list
allocate(vrl(NDIM,rgrid%npt))
do ir=1,rgrid%npt
  vrl(:,ir)=rgrid%vpl(ir)
end do
if (pars%nstates.gt.proj%norb) then
  call read_udis(pars%seedname,kgrid%npt,proj%norb,pars%nstates,vkl,udis)
end if
call read_umat(pars%seedname,kgrid%npt,proj%norb,vkl,umat)
! obtain wannier functions
call get_wfmloc(.false.,proj%norb,pars%nstates,kgrid%npt,rgrid%npt,&
   tbmodel%norb_TB,udis,umat,vkl,vrl,evec,wfmloc)
allocate(wws(proj%norb,proj%norb,sym%nsym))
allocate(wwsi(proj%norb,proj%norb,sym%nsym))
wws=0._dp
wwsi=0._dp
do isym=1,sym%nsym
  do iw=1,proj%norb
     ic=proj%orb_icio(iw,1)
     jc=proj%ics2c(ic,isym)
     do jw=1,proj%norb
        if(proj%orb_icio(jw,1).ne.jc) cycle
        wws(jw,iw,isym)=proj%wws(sym%car(:,:,isym),iw,jw)
        wwsi(jw,iw,isym)=proj%wws(sym%car(:,:,sym%inv(isym)),iw,jw)
     end do
  end do
end do
nn=tbmodel%norb_TB*rgrid%npt
do isym=1,sym%nsym
  write(*,*) "isym(wws): ",isym
  do iw=1,proj%norb
     ! check
     wf1=0._dp
     wf2=0._dp
     do jw=1,proj%norb
        wf1(:,:)=wf1(:,:)+wws(jw,iw,isym)*wfmloc(:,jw,:)
        wf2(:,:)=wf2(:,:)+wwsi(jw,iw,isym)*wfmloc(:,jw,:)
     end do
     do jw=1,proj%norb
       z1=zdotc(nn,wfmloc(:,jw,:),1,wf1,1)
       if (abs(z1-wws(jw,iw,isym)).gt.epslat) then
          write(*,*) iw,jw
          write(*,*) z1,wwsi(jw,iw,isym)
          call throw("wannier_interface%symmetrize_hubbardu_rs()",&
         "representation matrix derived from wfmloc is different from the one of projection block")
       end if
       z1=zdotc(nn,wfmloc(:,jw,:),1,wf2,1)
       if (abs(z1-wwsi(jw,iw,isym)).gt.epslat) then
          write(*,*) iw,jw
          write(*,*) z1,wwsi(jw,iw,isym)
          call throw("wannier_interface%symmetrize_hubbardu_rs()",&
         "inverse representation matrix derived from wfmloc is different from the one of projection block")
       end if
     end do
  end do
end do
deallocate(vrl)
call rgrid%init_sphere(pars)
allocate(vrl(NDIM,rgrid%npt_sphere))
nn=proj%norb
allocate(HubU(nn,nn,nn,nn,rgrid%npt_sphere))
allocate(HubUn(nn,nn,rgrid%npt_sphere))
allocate(HubUs(nn,nn,rgrid%npt_sphere))
open(140,file="UH.dat",action="read")
read(140,*) npt_sphere,norb
if (npt_sphere.ne.rgrid%npt_sphere) call throw("wannier_interface%symmetrize_hubbardu_rs()","npt_sphere is not what expected")
if (norb.ne.proj%norb) call throw("wannier_interface%symmetrize_hubbardu_rs()","npt_sphere is not what expected")
do jR_sphere=1,rgrid%npt_sphere
  read(140,*) vrl(:,jR_sphere)
  do nn=1,proj%norb
    do mm=1,proj%norb
      do pp=1,proj%norb
        do qq=1,proj%norb
          read(140,*) nn_,mm_,pp_,qq_,t1,t2
          HubU(nn,mm,pp,qq,jR_sphere)=cmplx(t1,t2,kind=8)
          if (nn.ne.nn_) call throw("wannier_interface%symmetrize_hubbardu_rs()","nn is not what expected")
          if (mm.ne.mm_) call throw("wannier_interface%symmetrize_hubbardu_rs()","mm is not what expected")
          if (pp.ne.pp_) call throw("wannier_interface%symmetrize_hubbardu_rs()","pp is not what expected")
          if (qq.ne.qq_) call throw("wannier_interface%symmetrize_hubbardu_rs()","qq is not what expected")
        end do
      end do
    end do
  end do
  read(140,*)
end do
close(140)
if (mp_mpi) then
  open(140,file="URN.dat",action="write")
  write(140,*) rgrid%npt_sphere,proj%norb
  do jR_sphere=1,rgrid%npt_sphere
    write(140,*) nint(rgrid%vpl_sphere(jR_sphere))
    do nn=1,proj%norb
      do mm=1,proj%norb
        write(140,'(2I6,2G18.10)') mm,nn,HubU(mm,mm,nn,nn,jR_sphere)
      end do
    end do
    write(140,*)
  end do
  close(140)
end if
HubUs=0._dp
do isym=1,sym%nsym
  write(*,*) "isym(symmetrize): ",isym
  do jR_sphere=1,rgrid%npt_sphere
    do iorb=1,proj%norb
      ic=proj%orb_icio(iorb,1)
      do jorb=1,proj%norb
        jc=proj%orb_icio(jorb,1)
        tauji=proj%centers(:,pars%proj%iw2ic(jorb))-proj%centers(:,pars%proj%iw2ic(iorb))
        ! next two loops are summation indices for "left" representation matrices pf D_il D_jm U D^-1_pi D^-1_qj
        ! NOTE, that U in this formula has the following definition (*=cc): U_lmpq<lm|W|pq>=l*[r1] m*[r2] W[r2-r1] p[r1] q[r2]
        ! different convention is used in U computation(*=cc): U_ijkt = i*[r1] j[r1] W(r2-r1) k*[r2] t[r2]
        do lorb=1,proj%norb
          lc=proj%orb_icio(lorb,1)
          if (proj%ics2c(lc,isym).ne.ic) cycle
          do morb=1,proj%norb
            mc=proj%orb_icio(morb,1)
            if (proj%ics2c(mc,isym).ne.jc) cycle
            tauml=proj%centers(:,pars%proj%iw2ic(morb))-proj%centers(:,pars%proj%iw2ic(lorb))
            ! next two loops are summation indeices for "right" representation matrices of the expressin above
            do porb=1,proj%norb
              pc=proj%orb_icio(porb,1)
              if (proj%ics2c(pc,isym).ne.ic) cycle
              if (pc.ne.lc) cycle
              do qorb=1,proj%norb
                qc=proj%orb_icio(qorb,1)
                if (proj%ics2c(qc,isym).ne.jc) cycle
                if (qc.ne.mc) cycle
                taupq=proj%centers(:,pars%proj%iw2ic(morb))-proj%centers(:,pars%proj%iw2ic(lorb))
                ! Tmlij from Dominik Gresch paper
                Tmlij=proj%Tmlij(dble(sym%lat(:,:,isym)),tauml,tauji)
                RR=matmul(dble(sym%lat(:,:,sym%inv(isym))),rgrid%vpl_sphere(jR_sphere)-Tmlij)
                iv=rgrid%find(RR)
                if (iv(NDIM+1).lt.0.or.sum(abs(iv(1:NDIM))).ne.0) cycle
                iR=rgrid%homo_to_sphere(iv(NDIM+1))
                if (iR.le.0) cycle
                ! inices in 4-point matrix element are correct, just different indices convention in symmetrisation formula
                HubUs(iorb,jorb,jR_sphere)=HubUs(iorb,jorb,jR_sphere)+&
                     wws(iorb,lorb,isym)*wws(jorb,morb,isym)*HubU(lorb,porb,morb,qorb,iR)*wwsi(porb,iorb,isym)*wwsi(qorb,jorb,isym)
              end do
            end do
          end do
        end do
      end do
    end do
  end do
end do
HubUs=HubUs/dble(sym%nsym)
if (mp_mpi) then
  open(140,file="URS.dat",action="write")
  write(140,*) rgrid%npt_sphere,proj%norb
  do jR_sphere=1,rgrid%npt_sphere
    write(140,*) nint(rgrid%vpl_sphere(jR_sphere))
    do nn=1,proj%norb
      do mm=1,proj%norb
        write(140,'(2I6,2G18.10)') mm,nn,HubUs(mm,nn,jR_sphere)
      end do
    end do
    write(140,*)
  end do
  close(140)
end if
deallocate(udis,umat,wfmloc)
deallocate(vkl,vrl)
deallocate(wws,wwsi)
deallocate(HubU,HubUn,HubUs)
deallocate(wf1,wf2)
end subroutine

subroutine compute_hubbardu_rs(pars,tbmodel,kgrid,evec)
class(CLpars), intent(in) :: pars
class(CLtb), intent(in) :: tbmodel
class(GRID), intent(inout) :: kgrid
complex(dp), intent(in) :: evec(tbmodel%norb_TB,pars%nstates,kgrid%npt)
real(dp), parameter :: Upz=10._dp
real(dp), parameter :: epscoul=10._dp
real(dp), parameter :: eps=1.e-17_dp
integer jR_sphere,iRp,iRpp
integer iorb,jorb,ir,ik
integer nn,mm,pp,qq
real(dp) dij,vcl
real(dp) v1(NDIM),v2(NDIM),rij(NDIM)
complex(dp) z1
complex(dp) zfn(pars%proj%norb)
complex(dp) zfm(pars%proj%norb)
complex(dp) zfp(pars%proj%norb)
complex(dp) zfq(pars%proj%norb)
real(dp), allocatable :: vkl(:,:),vrl(:,:),vpcorb(:,:)
complex(dp), allocatable :: wf1(:,:,:)
complex(dp), allocatable :: wf2(:,:,:)
complex(dp), allocatable :: HubU(:,:,:,:,:)
complex(dp), allocatable :: umat(:,:,:),udis(:,:,:)
complex(dp), allocatable :: wfmloc(:,:,:)
complex(dp), allocatable :: wf2z(:,:)
type(wbase)  proj
type(GRID) rgrid
type(coulrs) vcoul
!complex(dp) ddot
if (pars%proj%norb.le.0) then
   call throw("wannier_interface%compute_hubbardu()",&
              "apparently 'projections' block was not specified, wannier projections not found")
else
  if (pars%nstates.lt.pars%proj%norb) then
    call throw("wannier_interface%compute_hubbardu()",&
               "number of bands is less than the number of requested projections")
  end if
end if
call proj%init(pars,pars%proj%ncenters,pars%proj%norb,pars%proj%norb_ic,&
                   pars%proj%lmr,pars%proj%waxis,pars%proj%centers)
call info("Wannier_interface%compute_hubbardu","constructing U parameters on the grid")
call info("Wannier_interface%compute_hubbardu","WARNING this code is working for pz basis orbitals only")
call vcoul%init(pars%coulrs_file)
! copy rgrid object from TB
rgrid=tbmodel%rgrid
if (.not.rgrid%sphere_allocated) call rgrid%init_sphere(pars)
allocate(udis(pars%nstates,proj%norb,kgrid%npt))
allocate(umat(proj%norb,proj%norb,kgrid%npt))
allocate(wfmloc(tbmodel%norb_TB,proj%norb,rgrid%npt))
allocate(wf1(tbmodel%norb_TB,proj%norb,rgrid%npt_sphere))
nn=proj%norb
allocate(HubU(nn,nn,nn,nn,rgrid%npt_sphere))
! k-point list
allocate(vkl(NDIM,kgrid%npt))
do ik=1,kgrid%npt
  vkl(:,ik)=kgrid%vpl(ik)
end do
! r-point list
allocate(vrl(NDIM,rgrid%npt))
do ir=1,rgrid%npt
  vrl(:,ir)=rgrid%vpl(ir)
end do
if (pars%nstates.gt.proj%norb) then
  call read_udis(pars%seedname,kgrid%npt,proj%norb,pars%nstates,vkl,udis)
end if
call read_umat(pars%seedname,kgrid%npt,proj%norb,vkl,umat)
! obtain wannier functions
call get_wfmloc(.false.,proj%norb,pars%nstates,kgrid%npt,rgrid%npt,&
   tbmodel%norb_TB,udis,umat,vkl,vrl,evec,wfmloc)
do iRp=1,rgrid%npt_sphere
   wf1(:,:,iRp)=wfmloc(:,:,rgrid%sphere_to_homo(iRp))
end do
#ifdef MPI
  call mpi_barrier(mpi_com,mpi_err)
#endif
HubU=0._dp
!$OMP PARALLEL DEFAULT (SHARED)&
!$OMP PRIVATE(iRp,iRpp,iorb,jorb,rij,dij,v1,v2)&
!$OMP PRIVATE(mm,nn,pp,qq,vcl,z1)&
!$OMP PRIVATE(zfn,zfm,zfp,zfq)&
!$OMP PRIVATE(wf2z,wf2,vpcorb)
  allocate(wf2(tbmodel%norb_TB,proj%norb,rgrid%npt_sphere))
  allocate(wf2z(tbmodel%norb_TB,rgrid%npt))
  allocate(vpcorb(NDIM,tbmodel%norb_TB))
!$OMP DO
do jR_sphere=1,rgrid%npt_sphere
  ! coordinates of basis orbitals
  do iorb=1,tbmodel%norb_TB
    vpcorb(:,iorb)=matmul(tbmodel%vplorb(iorb),pars%avec)
  end do
  ! shift wave function to jR_sphere unit cell
  do nn=1,proj%norb
    call wannerfunc_at_R(rgrid,tbmodel%norb_TB,rgrid%vpl_sphere(jR_sphere),wfmloc(:,nn,:),wf2z(:,:))
    do iRp=1,rgrid%npt_sphere
      wf2(:,nn,iRp)=wf2z(:,rgrid%sphere_to_homo(iRp))
    end do
  end do
  do iRp=1,rgrid%npt_sphere
    if (mod(iRp-1,np_mpi).ne.lp_mpi) cycle
    !$OMP CRITICAL
    write(*,*) "JR*iRP: ",jR_sphere*iRp," of ",rgrid%npt_sphere*rgrid%npt_sphere
    !$OMP END CRITICAL
    do iRpp=1,rgrid%npt_sphere
      do iorb=1,tbmodel%norb_TB
        v1=vpcorb(:,iorb)
        do jorb=1,tbmodel%norb_TB
          v2=vpcorb(:,jorb)
          rij=abs(v2+rgrid%vpc_sphere(iRpp)-v1-rgrid%vpc_sphere(iRp))
          dij=sqrt(dot_product(rij,rij))
          if (dij.gt.pars%rcut_grid) cycle
          zfn(:)=conjg(wf1(iorb,:,iRp))
          zfm(:)=wf1(iorb,:,iRp)
          zfp(:)=conjg(wf2(jorb,:,iRpp))
          zfq(:)=wf2(jorb,:,iRpp)
          vcl=vcoul%evaluate(dij)
          do qq=1,proj%norb
            do pp=1,proj%norb
              if (pars%HubU_diagonal.and.pp.ne.qq) cycle 
              do mm=1,proj%norb
                do nn=1,proj%norb
                  if (pars%HubU_diagonal.and.nn.ne.mm) cycle 
                  z1=zfn(nn)*zfm(mm)*zfp(pp)*zfq(qq)
                  if (abs(z1).lt.eps) cycle
                  if (dij.lt.epslat) then
                     !HubU(nn,mm,pp,qq,jR_sphere)=HubU(nn,mm,pp,qq,jR_sphere)+z1*Upz/CoulombForceConstant
                     HubU(nn,mm,pp,qq,jR_sphere)=HubU(nn,mm,pp,qq,jR_sphere)+z1*Upz
                  else
                     !HubU(nn,mm,pp,qq,jR_sphere)=HubU(nn,mm,pp,qq,jR_sphere)+z1/dij
                     HubU(nn,mm,pp,qq,jR_sphere)=HubU(nn,mm,pp,qq,jR_sphere)+z1*vcl
                  end if
                end do
              end do
            end do 
          end do
        end do
      end do
    end do
  end do
end do
!$OMP END DO
  deallocate(vpcorb,wf2,wf2z)
!$OMP END PARALLEL
#ifdef MPI
  nn=rgrid%npt_sphere*proj%norb**4
  call mpi_allreduce(mpi_in_place,HubU,nn,mpi_double_complex,mpi_sum, &
   mpi_com,mpi_err)
#endif
!HubU=HubU*CoulombForceConstant/epscoul
if (mp_mpi) then
  open(140,file="UH.dat",action="write")
  write(140,*) rgrid%npt_sphere,proj%norb
  do jR_sphere=1,rgrid%npt_sphere
    write(140,*) nint(rgrid%vpl_sphere(jR_sphere))
    do nn=1,proj%norb
      do mm=1,proj%norb
        do pp=1,proj%norb
          do qq=1,proj%norb
            write(140,'(4I6,2G18.10)') nn,mm,pp,qq,HubU(nn,mm,pp,qq,jR_sphere)
          end do
        end do
      end do
    end do
    write(140,*)
  end do
  close(140)
end if
deallocate(vkl,vrl,umat,udis)
deallocate(wfmloc,HubU)
deallocate(wf1)
#ifdef MPI
  call mpi_barrier(mpi_com,mpi_err)
#endif
return
end subroutine

subroutine compute_hubbardj_rs(pars,tbmodel,kgrid,evec)
class(CLpars), intent(in) :: pars
class(CLtb), intent(in) :: tbmodel
class(GRID), intent(inout) :: kgrid
complex(dp), intent(in) :: evec(tbmodel%norb_TB,pars%nstates,kgrid%npt)
real(dp), parameter :: Upz=10._dp
real(dp), parameter :: epscoul=10._dp
integer jR_sphere,iRp,iRpp
integer iorb,jorb,nn,mm,ir,ik
real(dp) dij
complex(dp) zmni,zmnj
real(dp) v1(NDIM),v2(NDIM),rij(NDIM)
real(dp), allocatable :: vkl(:,:),vrl(:,:),vpcorb(:,:)
complex(dp), allocatable :: wf1(:,:,:)
complex(dp), allocatable :: wf2(:,:,:)
complex(dp), allocatable :: wf2z(:,:)
complex(dp), allocatable :: HubJ(:,:,:)
complex(dp), allocatable :: umat(:,:,:),udis(:,:,:)
complex(dp), allocatable :: wfmloc(:,:,:)
type(wbase)  proj
type(GRID) rgrid
!complex(dp) ddot
if (pars%proj%norb.le.0) then
   call throw("wannier_interface%compute_hubbardj()",&
              "apparently 'projections' block was not specified, wannier projections not found")
else
  if (pars%nstates.lt.pars%proj%norb) then
    call throw("wannier_interface%compute_hubbardj()",&
               "number of bands is less than the number of requested projections")
  end if
end if
call proj%init(pars,pars%proj%ncenters,pars%proj%norb,pars%proj%norb_ic,&
                   pars%proj%lmr,pars%proj%waxis,pars%proj%centers)
call info("Wannier_interface%compute_hubbardj","constructing J parameters on the grid")
call info("Wannier_interface%compute_hubbardj","WARNING this code is working for pz basis orbitals only")
! copy rgrid object from TB
rgrid=tbmodel%rgrid
if (.not.rgrid%sphere_allocated) call rgrid%init_sphere(pars)
allocate(udis(pars%nstates,proj%norb,kgrid%npt))
allocate(umat(proj%norb,proj%norb,kgrid%npt))
allocate(wfmloc(tbmodel%norb_TB,proj%norb,rgrid%npt))
allocate(wf1(tbmodel%norb_TB,proj%norb,rgrid%npt_sphere))
allocate(HubJ(proj%norb,proj%norb,rgrid%npt_sphere))
! k-point list
allocate(vkl(NDIM,kgrid%npt))
do ik=1,kgrid%npt
  vkl(:,ik)=kgrid%vpl(ik)
end do
! r-point list
allocate(vrl(NDIM,rgrid%npt))
do ir=1,rgrid%npt
  vrl(:,ir)=rgrid%vpl(ir)
end do
if (pars%nstates.gt.proj%norb) then
  call read_udis(pars%seedname,kgrid%npt,proj%norb,pars%nstates,vkl,udis)
end if
call read_umat(pars%seedname,kgrid%npt,proj%norb,vkl,umat)
! obtain wannier functions
call get_wfmloc(.false.,proj%norb,pars%nstates,kgrid%npt,rgrid%npt,&
   tbmodel%norb_TB,udis,umat,vkl,vrl,evec,wfmloc)
HubJ=0._dp
do iRp=1,rgrid%npt_sphere
  wf1(:,:,iRp)=wfmloc(:,:,rgrid%sphere_to_homo(iRp))
end do
#ifdef MPI
  call mpi_barrier(mpi_com,mpi_err)
#endif
!$OMP PARALLEL DEFAULT (SHARED)&
!$OMP PRIVATE(iRp,iRpp,iorb,jorb,nn,rij,dij,v1,v2)&
!$OMP PRIVATE(zmni,zmnj,vpcorb,wf2,wf2z)
  allocate(vpcorb(NDIM,tbmodel%norb_TB))
  allocate(wf2(tbmodel%norb_TB,proj%norb,rgrid%npt_sphere))
  allocate(wf2z(tbmodel%norb_TB,rgrid%npt))
!$OMP DO
do jR_sphere=1,rgrid%npt_sphere
  if (mod(jR_sphere-1,np_mpi).ne.lp_mpi) cycle
  !$OMP CRITICAL
  write(*,*) "JR: ",jR_sphere
  !$OMP END CRITICAL
  ! coordinates of basis orbitals
  do iorb=1,tbmodel%norb_TB
    vpcorb(:,iorb)=matmul(tbmodel%vplorb(iorb),pars%avec)
  end do
  do nn=1,proj%norb
    call wannerfunc_at_R(rgrid,tbmodel%norb_TB,rgrid%vpl_sphere(jR_sphere),wfmloc(:,nn,:),wf2z)
    do iRp=1,rgrid%npt_sphere
      wf2(:,nn,iRp)=wf2z(:,rgrid%sphere_to_homo(iRp))
    end do
  end do
  do iRp=1,rgrid%npt_sphere
    do iRpp=1,rgrid%npt_sphere
      do iorb=1,tbmodel%norb_TB
        v1=vpcorb(:,iorb)
        do jorb=1,tbmodel%norb_TB
          v2=vpcorb(:,jorb)
          rij=abs(v2+rgrid%vpc_sphere(iRpp)-v1-rgrid%vpc_sphere(iRp))
          dij=sqrt(dot_product(rij,rij))
          if (dij.gt.pars%rcut_grid) cycle
          do nn=1,proj%norb
            ! exchange density at nn
            do mm=1,proj%norb
              zmni=conjg(wf1(iorb,mm,iRp))*wf2(iorb,nn,iRp)
              zmnj=conjg(wf2(jorb,nn,iRpp))*wf1(jorb,mm,iRpp)
              if (dij.lt.epslat) then
                 HubJ(mm,nn,jR_sphere)=HubJ(mm,nn,jR_sphere)+zmni*zmnj*Upz/CoulombForceConstant
              else
                 HubJ(mm,nn,jR_sphere)=HubJ(mm,nn,jR_sphere)+zmni*zmnj/dij
              end if
            end do
          end do
        end do
      end do
    end do
  end do
end do
!$OMP END DO
  deallocate(vpcorb,wf2,wf2z)
!$OMP END PARALLEL
#ifdef MPI
  nn=proj%norb*proj%norb*rgrid%npt_sphere
  call mpi_allreduce(mpi_in_place,HubJ,nn,mpi_double_complex,mpi_sum, &
   mpi_com,mpi_err)
#endif
HubJ=HubJ*CoulombForceConstant/epscoul
if (mp_mpi) then
  open(140,file="JR.dat",action="write")
  write(140,*) rgrid%npt_sphere,proj%norb
  do jR_sphere=1,rgrid%npt_sphere
    write(140,*) nint(rgrid%vpl_sphere(jR_sphere))
    do nn=1,proj%norb
      do mm=1,proj%norb
        write(140,*) mm,nn,HubJ(mm,nn,jR_sphere)
      end do
    end do
    write(140,*)
  end do
  close(140)
end if
deallocate(vkl,vrl,umat,udis)
deallocate(wfmloc,HubJ)
deallocate(wf1)
#ifdef MPI
  call mpi_barrier(mpi_com,mpi_err)
#endif
return
end subroutine

!subroutine compute_hubbardu(pars,tbmodel,kgrid,sym,evec,Uhubbard)
!class(CLpars), intent(in) :: pars
!class(CLtb), intent(in) :: tbmodel
!class(GRID), intent(inout) :: kgrid
!class(CLsym), intent(in) :: sym
!complex(dp), intent(in) :: evec(tbmodel%norb_TB,pars%nstates,kgrid%npt)
!complex(dp), intent(out) :: Uhubbard(pars%proj%norb,pars%proj%norb,tbmodel%rgrid%npt)
!real(dp), parameter :: uonsite=10._dp
!integer ik_irr,ik,jk,ig,ir
!integer iorb,nn
!integer Ggrid_dims(NDIM)
!real(dp) vl(NDIM),vc(NDIM),vq(NDIM),dc,t1
!complex(dp) z1
!real(dp), allocatable :: vkl(:,:),vrl(:,:)
!complex(dp), allocatable :: umat(:,:,:),udis(:,:,:)
!complex(dp), allocatable :: psi_smooth(:,:,:)
!complex(dp), allocatable :: ovlpik(:,:)
!complex(dp), allocatable :: ovlpjk(:,:)
!complex(dp), allocatable :: HubUk(:,:,:,:)
!type(wbase)  proj
!type(GRID) rgrid,Ggrid
!!complex(dp) ddot
!if (pars%proj%norb.le.0) then
!   call throw("wannier_interface%compute_hubbardu()",&
!              "apparently 'projections' block was not specified, wannier projections not found")
!else
!  if (pars%nstates.lt.pars%proj%norb) then
!    call throw("wannier_interface%compute_hubbardu()",&
!               "number of bands is less than the number of requested projections")
!  end if
!end if
!call proj%init(pars,pars%proj%ncenters,pars%proj%norb,pars%proj%norb_ic,&
!                   pars%proj%lmr,pars%proj%waxis,pars%proj%centers)
!call info("Wannier_interface%compute_hubbardu","constructing U parameters on the grid")
!call info("Wannier_interface%compute_hubbardu","WARNING this code is working for surface only")
!! copy rgrid object from TB
!rgrid=tbmodel%rgrid
!! grid of reciprocal lattice points G such that q_FBZ+G samples whole reciprocal space
!if (sum(abs(pars%Ggrid)).eq.0) then
!  Ggrid_dims=pars%ngrid
!else 
!  Ggrid_dims=pars%Ggrid
!end if
!call Ggrid%init(Ggrid_dims,kgrid%vecs,.true.,.false.)
!call Ggrid%init_sphere(pars)
!if (.not.kgrid%syminit_done) call kgrid%sym_init(sym)
!allocate(udis(pars%nstates,proj%norb,kgrid%npt))
!allocate(umat(proj%norb,proj%norb,kgrid%npt))
!allocate(psi_smooth(tbmodel%norb_TB,proj%norb,kgrid%npt))
!allocate(ovlpik(proj%norb,proj%norb))
!allocate(ovlpjk(proj%norb,proj%norb))
!allocate(HubUk(proj%norb,proj%norb,kgrid%nir,kgrid%npt))
!! k-point list
!allocate(vkl(NDIM,kgrid%npt))
!do ik=1,kgrid%npt
!  vkl(:,ik)=kgrid%vpl(ik)
!end do
! r-point list
!allocate(vrl(NDIM,rgrid%npt))
!do ir=1,rgrid%npt
!  vrl(:,ir)=rgrid%vpl(ir)
!end do
!if (pars%nstates.gt.proj%norb) then
!  call read_udis(pars%seedname,kgrid%npt,proj%norb,pars%nstates,vkl,udis)
!end if
!call read_umat(pars%seedname,kgrid%npt,proj%norb,vkl,umat)
!! obtain a smooth gauge bloch functions
!call get_wfmloc(.true.,proj%norb,pars%nstates,kgrid%npt,rgrid%npt,&
!   tbmodel%norb_TB,udis,umat,vkl,vrl,evec,psi_smooth)
!!psi_smooth=evec
!HubUk=0._dp
!do ik_irr=1,kgrid%nir
!  ik=kgrid%ir2ik(ik_irr)
!  do jk=1,kgrid%npt
!    do ig=1,Ggrid%npt_sphere
!      ! skip G=0 case
!      if (Ggrid%sphere_to_homo(ig).eq.Ggrid%ip0) cycle
!      ovlpik=0._dp
!      ovlpjk=0._dp
!      do iorb=1,tbmodel%norb_TB
!        t1=dot_product(Ggrid%vpl_sphere(ig),tbmodel%vplorb(iorb))*twopi
!        z1=cmplx(cos(t1),sin(t1),kind=dp)
!        do nn=1,proj%norb
!          ! for hubbard term we need the density componet, diagonal at k
!          ovlpik(:,nn)=ovlpik(:,nn)+conjg(psi_smooth(iorb,:,ik))*psi_smooth(iorb,nn,ik)*z1*pwave_ovlp(Ggrid%dc_sphere(ig))
!          ovlpjk(:,nn)=ovlpjk(:,nn)+conjg(psi_smooth(iorb,:,jk))*psi_smooth(iorb,nn,jk)*z1*pwave_ovlp(Ggrid%dc_sphere(ig))
!        end do
!      end do
!      HubUk(:,:,ik_irr,jk)=HubUk(:,:,ik_irr,jk)+ovlpik(:,:)*conjg(ovlpjk(:,:))*fourpi/Ggrid%dc_sphere(ig)**2
!    end do
!    ! normalise to surface area
!    HubUk(:,:,ik_irr,jk)=(HubUk(:,:,ik_irr,jk)*abohr**2)*Hartree_to_eV&
!       /abs( pars%avec(1,1)*pars%avec(2,2)-pars%avec(1,2)*pars%avec(2,1) )
!    write(150,'(2I6,2G18.10)') ik,jk,dble(HubUk(1,1,ik_irr,jk))
!  end do
!end do
!deallocate(vkl,vrl,umat,udis)
!deallocate(ovlpik,ovlpjk)
!deallocate(psi_smooth,HubUk)
!return
!end subroutine
!
!subroutine compute_hubbardj(pars,tbmodel,kgrid,sym,evec,Uhubbard)
!class(CLpars), intent(in) :: pars
!class(CLtb), intent(in) :: tbmodel
!class(GRID), intent(inout) :: kgrid
!class(CLsym), intent(in) :: sym
!complex(dp), intent(in) :: evec(tbmodel%norb_TB,pars%nstates,kgrid%npt)
!complex(dp), intent(out) :: Uhubbard(pars%proj%norb,pars%proj%norb,tbmodel%rgrid%npt)
!real(dp), parameter :: uonsite=10._dp
!integer ik_irr,ik,jk,ig,ir
!integer isym,jk_image
!integer iorb,nn
!integer Ggrid_dims(NDIM)
!real(dp) vl(NDIM),vc(NDIM),vq(NDIM),dc,t1
!complex(dp) z1
!real(dp), allocatable :: vkl(:,:),vrl(:,:)
!complex(dp), allocatable :: umat(:,:,:),udis(:,:,:)
!complex(dp), allocatable :: psi_smooth(:,:,:)
!complex(dp), allocatable :: ovlp(:,:)
!complex(dp), allocatable :: HubUk(:,:,:,:)
!type(wbase)  proj
!type(GRID) rgrid,Ggrid
!!complex(dp) ddot
!if (pars%proj%norb.le.0) then
!   call throw("wannier_interface%compute_hubbardu()",&
!              "apparently 'projections' block was not specified, wannier projections not found")
!else
!  if (pars%nstates.lt.pars%proj%norb) then
!    call throw("wannier_interface%compute_hubbardu()",&
!               "number of bands is less than the number of requested projections")
!  end if
!end if
!call proj%init(pars,pars%proj%ncenters,pars%proj%norb,pars%proj%norb_ic,&
!                   pars%proj%lmr,pars%proj%waxis,pars%proj%centers)
!call info("Wannier_interface%compute_hubbardu","constructing U parameters on the grid")
!call info("Wannier_interface%compute_hubbardu","WARNING this code is working for surface only")
!! copy rgrid object from TB
!rgrid=tbmodel%rgrid
!! grid of reciprocal lattice points G such that q_FBZ+G samples whole reciprocal space
!if (sum(abs(pars%Ggrid)).eq.0) then
!  Ggrid_dims=pars%ngrid
!else 
!  Ggrid_dims=pars%Ggrid
!end if
!call Ggrid%init(Ggrid_dims,kgrid%vecs,.true.,.false.)
!call Ggrid%init_sphere(pars)
!if (.not.kgrid%syminit_done) call kgrid%sym_init(sym)
!allocate(udis(pars%nstates,proj%norb,kgrid%npt))
!allocate(umat(proj%norb,proj%norb,kgrid%npt))
!allocate(psi_smooth(tbmodel%norb_TB,proj%norb,kgrid%npt))
!allocate(ovlp(proj%norb,proj%norb))
!!allocate(HubUk(proj%norb,proj%norb,kgrid%nir,kgrid%npt))
!allocate(HubUk(proj%norb,proj%norb,kgrid%npt,kgrid%npt))
!! k-point list
!allocate(vkl(NDIM,kgrid%npt))
!do ik=1,kgrid%npt
!  vkl(:,ik)=kgrid%vpl(ik)
!end do
!! r-point list
!allocate(vrl(NDIM,rgrid%npt))
!do ir=1,rgrid%npt
!  vrl(:,ir)=rgrid%vpl(ir)
!end do
!if (pars%nstates.gt.proj%norb) then
!  call read_udis(pars%seedname,kgrid%npt,proj%norb,pars%nstates,vkl,udis)
!end if
!call read_umat(pars%seedname,kgrid%npt,proj%norb,vkl,umat)
!! obtain a smooth gauge bloch functions
!call get_wfmloc(.true.,proj%norb,pars%nstates,kgrid%npt,rgrid%npt,&
!   tbmodel%norb_TB,udis,umat,vkl,vrl,evec,psi_smooth)
!!psi_smooth=evec
!ik=5
!jk=7
!isym=kgrid%sik2ir(ik)
!ik_irr=kgrid%ik2ir(ik)
!jk_image=kgrid%iks2k(jk,isym)
!write(*,*) ik,jk
!write(*,'(6F10.4)') kgrid%vpl(ik),kgrid%vpl(jk)
!write(*,*) ik_irr,jk_image
!write(*,'(6F10.4)') kgrid%vpl(kgrid%ir2ik(ik_irr)),kgrid%vpl(jk_image)
!vc=matmul(dble(sym%lat(:,:,isym)),kgrid%vpl(ik))
!vq=matmul(dble(sym%lat(:,:,isym)),kgrid%vpl(jk))
!write(*,'(6F10.4)')vc,vq
!write(*,*) isym
!HubUk=0._dp
!do ik_irr=1,kgrid%nir
!  ik=kgrid%ir2ik(ik_irr)
!  do jk=1,kgrid%npt
!    do ig=1,Ggrid%npt_sphere
!      ! skip q+G=0 case
!      vq=kgrid%vpl(ik)-kgrid%vpl(jk)+Ggrid%vpl_sphere(ig)
!      vc=matmul(vq,kgrid%vecs)
!      dc=sqrt(dot_product(vc,vc))
!      if (abs(dc).lt.epslat) cycle
!      ovlp=0._dp
!      do iorb=1,tbmodel%norb_TB
!        t1=dot_product(vq,tbmodel%vplorb(iorb))*twopi
!        z1=cmplx(cos(t1),sin(t1),kind=dp)
!        do nn=1,proj%norb
!          ! for hubbard term we need the density componet, diagonal at k
!          ovlp(:,nn)=ovlp(:,nn)+conjg(psi_smooth(iorb,:,ik))*psi_smooth(iorb,nn,jk)*z1*pwave_ovlp(dc)
!        end do
!      end do
!      HubUk(:,:,ik_irr,jk)=HubUk(:,:,ik_irr,jk)+ovlp(:,:)*conjg(ovlp(:,:))*fourpi/dc**2
!    end do
!    ! normalise to surface area
!    HubUk(:,:,ik_irr,jk)=(HubUk(:,:,ik_irr,jk)*abohr**2)*Hartree_to_eV&
!       /abs( pars%avec(1,1)*pars%avec(2,2)-pars%avec(1,1)*pars%avec(2,1) )
!    write(250,'(2I6,2G18.10)') ik,jk,dble(HubUk(1,1,ik_irr,jk))
!      vq=kgrid%vpl(ik)-kgrid%vpl(jk)
!      vc=matmul(vq,kgrid%vecs)
!      dc=sqrt(dot_product(vc,vc))
!    write(251,*) dc,abs(HubUk(1,1,ik_irr,jk))
!  end do
!end do
!do ik=1,kgrid%npt
!   do jk=1,kgrid%npt
!     call get2pfk(kgrid,ik,jk,proj%norb,HubUk(:,:,1:kgrid%nir,:),ovlp)
!     write(253,'(2I6,2G18.10)') ik,jk,dble(ovlp(1,1))
!   end do
!end do
!HubUk=0._dp
!! FULL MATRIX
!do ik=1,kgrid%npt
!  do jk=1,kgrid%npt
!    do ig=1,Ggrid%npt_sphere
!      ! skip q+G=0 case
!      if (Ggrid%sphere_to_homo(ig).eq.Ggrid%ip0.and.ik.eq.jk) cycle
!      vq=kgrid%vpl(ik)-kgrid%vpl(jk)+Ggrid%vpl_sphere(ig)
!      vc=matmul(vq,kgrid%vecs)
!      dc=sqrt(dot_product(vc,vc))
!      ovlp=0._dp
!      do iorb=1,tbmodel%norb_TB
!        t1=dot_product(vq,tbmodel%vplorb(iorb))*twopi
!        z1=cmplx(cos(t1),sin(t1),kind=dp)
!        do nn=1,proj%norb
!          ! for hubbard term we need the density componet, diagonal at k
!          ovlp(:,nn)=ovlp(:,nn)+conjg(psi_smooth(iorb,:,ik))*psi_smooth(iorb,nn,jk)*z1*pwave_ovlp(dc)
!        end do
!      end do
!      HubUk(:,:,ik,jk)=HubUk(:,:,ik,jk)+ovlp(:,:)*conjg(ovlp(:,:))*fourpi/dc**2
!    end do
!    ! normalise to surface area
!    HubUk(:,:,ik,jk)=(HubUk(:,:,ik,jk)*abohr**2)*Hartree_to_eV&
!       /abs( pars%avec(1,1)*pars%avec(2,2)-pars%avec(1,1)*pars%avec(2,1) )
!    write(254,'(2I6,2G18.10)') ik,jk,dble(HubUk(1,1,ik,jk))
!  end do
!end do
!stop
!deallocate(vkl,vrl,umat,udis)
!deallocate(ovlp)
!deallocate(psi_smooth,HubUk)
!return
!end subroutine

!subroutine get2pfk(kgrid,ik,jk,nsize,func_irr,func)
!class(GRID), intent(in) :: kgrid
!integer, intent(in) :: ik,jk,nsize
!complex(dp), intent(in) :: func_irr(nsize,nsize,kgrid%nir,kgrid%npt)
!complex(dp), intent(out) :: func(nsize,nsize)
!integer ik_irr,jk_image,isym
!isym=kgrid%sik2ir(ik)
!ik_irr=kgrid%ik2ir(ik)
!jk_image=kgrid%iks2k(jk,isym)
!func(:,:)=func_irr(:,:,ik_irr,jk_image)
!end subroutine

! Suppplementary functions
subroutine read_udis(seedname,nk,nwan,num_bands,vkl,udis)
character(len=*), intent(in) :: seedname
integer, intent(in) :: nk,nwan,num_bands
real(dp), intent(in) :: vkl(NDIM,nk)
complex(dp), intent(out) :: udis(num_bands,nwan,nk)
logical exs
integer nk_,nwan_,num_bands_
integer ik,m,n
real(dp) t1,t2,vpl_(NDIM)
inquire(file=trim(adjustl(seedname))//'_u_dis.in',exist=exs)
if (.not.exs) then
  call throw("Wannier_interfacey%read_udis","no "//trim(adjustl(seedname))//'_u_dis.in'//" file")
  stop
end if
open(50,file=trim(adjustl(seedname))//'_u_dis.in',action='read')
read(50,*)
read(50,*) nk_,nwan_,num_bands_

if (nk_.ne.nk)  call throw("wannier_interface%read_udis()",&
             "number of k-points is different from the requested one")
if (nwan_.ne.nwan)  call throw("wannier_interface%read_udis()",&
             "number projections in wannier files is different from  the requested one")
if (num_bands_.ne.num_bands)  call throw("wannier_interface%read_udis()",&
             "number of bands in wannier files is different from the requested one")

do ik=1,nk_
  read(50,*)
  read(50,*) vpl_
  if (sum(abs(vpl_-vkl(:,ik))).gt.epslat) then
     write(*,*) ik
     write(*,*) vpl_
     write(*,*) vkl(:,ik)
     call throw("wannier_interface%read_udis()","k-vector not found")
  end if 
  do m=1,nwan_
    do n=1,num_bands_
      read(50,*) t1,t2
      udis(n,m,ik)=cmplx(t1,t2,kind=dp)
    end do
  end do
end do
close(50)
end subroutine
subroutine read_umat(seedname,nk,nwan,vkl,umat)
character(len=*), intent(in) :: seedname
integer, intent(in) :: nk,nwan
real(dp), intent(in) :: vkl(NDIM,nk)
complex(dp), intent(out) :: umat(nwan,nwan,nk)
logical exs
integer nk_,nwan_
integer ik,m,n
real(dp) t1,t2,vpl_(NDIM)
inquire(file=trim(adjustl(seedname))//'_u.in',exist=exs)
if (.not.exs) then
  call throw("Wannier_interface%read_umat","no "//trim(adjustl(seedname))//'_u.in'//" file")
  stop
end if
open(51,file=trim(adjustl(seedname))//'_u.in',action='read')
read(51,*)
read(51,*) nk_,nwan_,n

if (nk_.ne.nk)  call throw("wannier_interface%read_umat()",&
             "number of k-points is different from the requested one")
if (nwan_.ne.nwan)  call throw("wannier_interface%read_umat()",&
             "number projections in wannier files is different from the requested one")
do ik=1,nk_
  read(51,*)
  read(51,*) vpl_
  if (sum(abs(vpl_-vkl(:,ik))).gt.epslat) then
     write(*,*) ik
     write(*,*) vpl_
     write(*,*) vkl(:,ik)
     call throw("wannier_interface%read_umat()","different k-vectors")
  end if 
  do n=1,nwan_
    do m=1,nwan_
      read(51,*) t1,t2
      umat(m,n,ik)=cmplx(t1,t2,kind=dp)
    end do
  end do
end do
close(51)
end subroutine
subroutine get_wfmloc(get_smooth,nwan,num_bands,nkpt,nrpt,nbasis,udis,umat,vkl,vrl,evec,wfmloc)
logical, intent(in) :: get_smooth
integer, intent(in) :: nwan,num_bands,nkpt,nrpt,nbasis
real(dp), intent(in) :: vkl(NDIM,nkpt),vrl(NDIM,nkpt)
complex(dp), intent(in) :: udis(num_bands,nwan,nkpt),umat(nwan,nwan,nkpt)
complex(dp), intent(in) :: evec(nbasis,num_bands,nkpt)
complex(dp), intent(out) :: wfmloc(nbasis,nwan,nrpt)
integer ik,ir,n,m
real(dp) t1
complex(dp) z1
complex(dp), allocatable :: psik(:,:)
if (nrpt.ne.nkpt.and.get_smooth) call throw("wannier_interface%get_wfmloc","at get_smooth=.true. nrpt must be equal to nkpt")
allocate(psik(nbasis,nwan))
wfmloc(:,:,:)=0.d0
#ifdef MPI
  call mpi_barrier(mpi_com,mpi_err)
#endif
do ik=1,nkpt
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  if (num_bands.gt.nwan) then
    psik(:,:)=0._dp
    do n=1,nwan
      do m=1,num_bands
        psik(:,n)=psik(:,n)+udis(m,n,ik)*evec(:,m,ik)
      end do
    end do
  else if (num_bands.eq.nwan) then
    psik=evec(:,:,ik)
  else
    call throw("Wannier_interface%get_wfmloc"," num_bands<nwan")
  end if
  ! "smooth" bloch function
  psik=matmul(psik,umat(:,:,ik))
  if (get_smooth) then
     wfmloc(:,:,ik)=psik
     cycle
  end if
!$OMP PARALLEL DEFAULT (SHARED)&
!$OMP PRIVATE(t1,z1)
!$OMP DO
  do ir=1,nrpt
    t1=twopi*dot_product(vkl(:,ik),vrl(:,ir))
    z1=cmplx(cos(t1),sin(t1),kind=dp)
    wfmloc(:,:,ir)=wfmloc(:,:,ir)+psik*z1
  end do
!$OMP END DO
!$OMP END PARALLEL
end do
#ifdef MPI
  n=nbasis*nwan*nrpt
  call mpi_allreduce(mpi_in_place,wfmloc,n,mpi_double_complex,mpi_sum, &
   mpi_com,mpi_err)
#endif
if (.not.get_smooth) wfmloc=wfmloc/dble(nkpt)
deallocate(psik)
return
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
write(50,*) 'begin Projections'
do iwan=1,pars%proj%norb
  write(50,'("f=",G18.10,",",G18.10,",",G18.10,":",A4)') &
    pars%proj%centers(:,pars%proj%iw2ic(iwan)),lmr_to_string(pars%proj%lmr(:,iwan))
end do
write(50,*) 'end Projections'
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
write(50,*)'begin Unit_Cell_Cart'
write(50,'(3g18.10)') pars%avec(1,:)
write(50,'(3g18.10)') pars%avec(2,:)
write(50,'(3g18.10)') pars%avec(3,:)
write(50,*)'end Unit_Cell_Cart'
write(50,*)
write(50,*)'begin Atoms_Frac'
do iwan=1,pars%proj%norb
  write(50,'(a2,3G18.10)')'XX',pars%proj%centers(:,pars%proj%iw2ic(iwan))
end do
write(50,*)'end Atoms_Frac'
write(50,*)
write(50,*)'begin kpoints'
do ik=1,kgrid%npt
  write(50,'(3F12.8)') kgrid%vpl(ik)
end do
write(50,*)'end kpoints'
close(50)
return
end subroutine

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

subroutine wannerfunc_at_R(Rgrid,nbasis,Rvec,wfi,wfo)
class(GRID), intent(in) :: Rgrid
integer, intent(in) :: nbasis
real(dp), intent(in) :: Rvec(NDIM)
complex(dp), intent(in) :: wfi(nbasis,Rgrid%npt)
complex(dp), intent(out) :: wfo(nbasis,Rgrid%npt)
integer ir,irg(NDIM+1)
real(dp) vl(ndIM)
wfo(:,:)=0.d0
do ir=1,Rgrid%npt
  vl=Rgrid%vpl(ir)-Rvec
  irg=rgrid%find(vl)
  if (irg(ndim+1).lt.0.or.sum(abs(irg(1:NDIM))).ne.0) cycle
  wfo(:,ir)=wfi(:,irg(NDIM+1))
end do
return
end subroutine

!subroutine fourier_transform(ch,rgrid,kgrid,zf)
!character(len=2), intent(in) :: ch
!class(GRID), intent(in) :: rgrid,kgrid
!complex(dp), intent(inout) :: zf(*)
!complex(dp), allocatable :: zft(:)
!integer ik,ir,npmax
!real(dp) t1
!complex(dp) z1
!! if one uses FFT in the future, this will be helpful
!if (rgrid%npt.ne.kgrid%npt) call throw("Wannier_interface%fourier_transform()","One must have rgrid and kgrid of the same size")
!! this will allow usage with different grids
!npmax=max(rgrid%npt,kgrid%npt)
!allocate(zft(npmax))
!if (ch.eq.'rk') then
!  zft(1:rgrid%npt)=zf(1:rgrid%npt)
!  zf(1:npmax)=0._dp
!  do ik=1,kgrid%npt
!    do ir=1,rgrid%npt
!      t1=twopi*dot_product(rgrid%vpl(ir),kgrid%vpl(ik))
!      z1=cmplx(cos(t1),-sin(t1),kind=dp)
!      zf(ik)=zf(ik)+zft(ir)*z1
!    end do
!  end do
!else if (ch.eq.'kr') then
!  zft(1:kgrid%npt)=zf(1:kgrid%npt)
!  zf(1:npmax)=0._dp
!  do ir=1,rgrid%npt
!    do ik=1,kgrid%npt
!      t1=twopi*dot_product(rgrid%vpl(ir),kgrid%vpl(ik))
!      z1=cmplx(cos(t1),sin(t1),kind=dp)
!      zf(ir)=zf(ir)+zft(ik)*z1
!    end do
!  end do
!  zf(1:npmax)=zf(1:npmax)/dble(kgrid%npt)
!else
!  call throw("Wannier_interface%fourier_transform()","unknown FT character option")
!end if
!deallocate(zft)
!return
!end subroutine

subroutine read_tb_file(THIS,pars,norb)
class(CLwan), intent(inout) :: THIS
class(CLpars), intent(in) :: pars
integer, intent(in) :: norb
integer i,j,ii,jj
integer nrpt,nwan,ir
integer ngrid(3),ivp(4)
real(dp) aa,bb,vpl(3),avec(3,3)
logical exs
integer, allocatable :: ivr(:,:)
integer, allocatable :: deg(:)
complex(dp), allocatable :: ham(:,:,:)
if (trim(adjustl(pars%TBfile)).eq."") call throw("CLwan%read_tb_file",&
 "Apparently tbfile is not set via 'tbtype' input entry")
if (NDIM.ne.3) call throw("CLwan%read_tb_file","this subroutine works only in 3D case")
inquire(file=trim(adjustl(pars%TBfile)),exist=exs)
if (.not.exs) call throw("CLwan%read_tb_file"," file "//trim(adjustl(pars%TBfile))//" missing")
open(50,file=trim(adjustl(pars%TBfile)),action='read')
read(50,*)
read(50,*) avec(1,:)
read(50,*) avec(2,:)
read(50,*) avec(3,:)
if (sum(abs(avec-pars%avec)).gt.epslat) then
  write(*,*) sum(abs(avec-pars%avec))
  write(*,*) "avec from pars: "
  write(*,*) pars%avec(1,:)
  write(*,*) pars%avec(2,:)
  write(*,*) pars%avec(3,:)
  write(*,*) "avec in tbfile: "
  write(*,*) avec(1,:)
  write(*,*) avec(2,:)
  write(*,*) avec(3,:)
  call throw("CLwan%read_tb_file",&
  "lattice vectors in _tb file are different from ones in input. it can be resolved by writing more digits")
end if
read(50,*) nwan
read(50,*) nrpt
if (nwan.ne.norb) then
  call throw("CLwan%read_tb_file","number of basis orbitals in _tb file is different from one derived from the input")
end if
allocate(ivr(3,nrpt))
allocate(deg(nrpt))
allocate(ham(nwan,nwan,nrpt))
read(50,'(15i5)') (deg(ir),ir=1,nrpt)
do ir=1,nrpt
  read(50,*)
  read(50,*) ivr(:,ir)
  do j=1,nwan
    do i=1,nwan
      read(50,*) ii,jj,aa,bb
      ham(ii,jj,ir)=cmplx(aa,bb,kind=dp)
    end do
  end do
end do
close(50)
ngrid(1)=max(2*maxval(abs(ivr(1,:))),1)+1
ngrid(2)=max(2*maxval(abs(ivr(2,:))),1)+1
ngrid(3)=max(2*maxval(abs(ivr(3,:))),1)+1
call THIS%rgrid%init(ngrid,avec,.true.,.false.)
allocate(This%dege(THIS%rgrid%npt))
allocate(This%hame(nwan,nwan,THIS%rgrid%npt))
THIS%dege=1._dp
THIS%hame=0._dp
do ir=1,nrpt
  vpl=dble(ivr(:,ir))
  ivp=THIS%rgrid%find(vpl)
  if (ivp(ndim+1).lt.0.or.sum(abs(ivp(1:NDIM))).ne.0) cycle
  THIS%dege(ivp(4))=dble(deg(ir))
  THIS%hame(:,:,ivp(4))=ham(:,:,ir)
end do
deallocate(ivr,deg,ham)
return
end subroutine
subroutine write_tb_file(THIS,pars,norb)
class(CLwan), intent(inout) :: THIS
class(CLpars), intent(in) :: pars
integer, intent(in) :: norb
integer ii,jj,ir,ivp(NDIM+1),counter
type(GRID) rgrid
integer, allocatable :: deg(:)
if (NDIM.ne.3) call throw("CLwan%write_tb_file","this subroutine works only in 3D case")
call rgrid%init(pars%ngrid,pars%avec,.true.,.false.)
allocate(deg(rgrid%npt))
! collect degeneracies, count active R-points
counter=0
do ir=1,rgrid%npt
  ivp=THIS%rgrid%find(rgrid%vpl(ir))
  if (ivp(ndim+1).lt.0.or.sum(abs(ivp(1:NDIM))).ne.0) cycle
  counter=counter+1
  deg(counter)=nint(THIS%dege(ivp(NDIM+1)))
end do
call system("mkdir -p _ham")
open(50,file="_ham/hamwan_tb.dat",action='write')
write(50,*)
write(50,*) pars%avec(1,:)
write(50,*) pars%avec(2,:)
write(50,*) pars%avec(3,:)
write(50,*) norb
write(50,*) counter
write(50,'(15I5)') (deg(ir),ir=1,counter)
do ir=1,rgrid%npt
  ivp=THIS%rgrid%find(rgrid%vpl(ir))
  if (ivp(ndim+1).lt.0.or.sum(abs(ivp(1:NDIM))).ne.0) cycle
  write(50,*)
  write(50,'(3I5)') THIS%rgrid%vpi(ivp(NDIM+1))
  do jj=1,norb
    do ii=1,norb
      write(50,'(2I5,3x,2(E15.8,1x))') ii,jj,THIS%hame(ii,jj,ivp(NDIM+1))
    end do
  end do
end do
close(50)
deallocate(deg)
return
end subroutine


end module
