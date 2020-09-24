
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
public :: write_hubbardu
public :: symmetrize_tbfile
public :: symmetrize_hubbardu

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
    trim(adjustl(pars%wannier_proj_mode)).eq.'trial_from_gamma'.or.&
    trim(adjustl(pars%wannier_proj_mode)).eq.'trial_from_ft'.or.&
    trim(adjustl(pars%wannier_proj_mode)).eq.'wannier_file'.or.&
    trim(adjustl(pars%wannier_proj_mode)).eq.'wftrial_file'.or.&
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

subroutine projection(THIS,tbmodel,pars,sym,kgrid,eval,evec)
class(CLwan), intent(inout) :: THIS
class(CLtb), intent(in) :: tbmodel
class(CLpars), intent(in) :: pars
class(CLsym), intent(in) :: sym
class(GRID), intent(inout) :: kgrid
real(dp), intent(in) :: eval(pars%nstates,kgrid%npt)
complex(dp), intent(in) :: evec(tbmodel%norb_TB,pars%nstates,kgrid%npt)
! local
character(len=200) :: message
logical exs
integer iR,iw,i1,ival,icnd,ik
integer iorb,jorb,ispec
integer ipro,ic,jc,ios,l1,m1,l2,m2
integer isym1,isym2,ik_gamma,ikg(NDIM+1)
integer pm_val(2),pm_con(2)
type(wbase) proj
real(dp) sigma
real(dp) vpl(NDIM),x1(NDIM),x2(NDIM),z1(NDIM),z2(NDIM)
real(dp) t1,dd,dv(NDIM),dc(NDIM)
!complex(dp), allocatable :: wf_t(:)
complex(dp), allocatable :: wftrial(:,:,:)
complex(dp), allocatable :: wws(:,:,:)
complex(dp), allocatable :: phs(:,:)
real(dp), allocatable :: dis_factors(:,:)
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
  call generate_dmn_orb(tbmodel,proj,sym,pars,kgrid,evec,.false.,wws)
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
  call generate_amn_overlap(tbmodel,pars,kgrid,eval,evec,tbmodel%rgrid%npt,wftrial)
  deallocate(wftrial)

else if (trim(adjustl(pars%wannier_proj_mode)).eq.'trial_from_gamma') then

  inquire(file=trim(adjustl(pars%seedname))//'.amn',exist=exs)
  if (.not.exs) then
    allocate(wftrial(tbmodel%norb_TB,pars%proj%norb,tbmodel%rgrid%npt))
    do iorb=1,tbmodel%norb_TB
      ispec=tbmodel%orb_ispec(iorb)
      vpl=tbmodel%vplorb(iorb)
      do iw=1,proj%norb
        i1=pars%iflat_band
        ! bottom of bands to be include in trial
        i1=i1-4
        ival=i1+iw-1
        icnd=i1+proj%norb-iw
        !write(*,*) iw,i1,ival,icnd
        if (vpl(3).gt.0._dp) then
          ! top layer
          if (ispec.eq.1) then
            ! A-site
            wftrial(iorb,iw,1)=evec(iorb,ival,ik_gamma)
          else
            ! B-site
            wftrial(iorb,iw,1)=evec(iorb,icnd,ik_gamma)
          end if 
        else
          ! bottom layer
          if (ispec.eq.1) then
            ! A-site
            wftrial(iorb,iw,1)=evec(iorb,icnd,ik_gamma)
          else
            ! B-site
            wftrial(iorb,iw,1)=evec(iorb,ival,ik_gamma)
          end if 
        end if
      end do
    end do
    do iR=2,tbmodel%rgrid%npt
      wftrial(:,:,iR)=wftrial(:,:,iR-1)
    end do
    sigma=0.6_dp*sqrt(dot_product(pars%avec(1,:),pars%avec(1,:)))
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
    wftrial=dble(wftrial)
    call generate_amn_overlap(tbmodel,pars,kgrid,eval,evec,tbmodel%rgrid%npt,wftrial)
    deallocate(wftrial)
  end if
  call generate_dmn_orb(tbmodel,proj,sym,pars,kgrid,evec,.false.,wws)

else if (trim(adjustl(pars%wannier_proj_mode)).eq.'trial_from_ft') then

  inquire(file=trim(adjustl(pars%seedname))//'.amn',exist=exs)
  if (.not.exs) then
    allocate(wftrial(tbmodel%norb_TB,pars%proj%norb,tbmodel%rgrid%npt))
    allocate(phs(kgrid%npt,tbmodel%rgrid%npt))
    allocate(dis_factors(pars%nstates,kgrid%npt))
    do ik=1,kgrid%npt
      do iR=1,tbmodel%rgrid%npt
        t1=twopi*dot_product(kgrid%vpl(ik),tbmodel%rgrid%vpl(iR))
        phs(ik,iR)=cmplx(cos(t1),-sin(t1),kind=dp)
      end do
      do iw=1,pars%nstates
        if (eval(iw,ik).lt.pars%dis_frozen(1).or.eval(iw,ik).gt.pars%dis_frozen(2)) then
           dis_factors(iw,ik)=0._dp
        else
           dis_factors(iw,ik)=1._dp
        end if
      end do
    end do
    do iR=1,tbmodel%rgrid%npt
      do iorb=1,tbmodel%norb_TB
        ispec=tbmodel%orb_ispec(iorb)
        vpl=tbmodel%vplorb(iorb)
        do iw=1,proj%norb
          i1=pars%iflat_band
          ! bottom of bands to be include in trial
          i1=i1-4
          ival=i1+iw-1
          icnd=i1+proj%norb-iw
          !write(*,*) iw,i1,ival,icnd
          if (vpl(3).gt.0._dp) then
            ! top layer
            if (ispec.eq.1) then
              ! A-site
              wftrial(iorb,iw,iR)=sum(evec(iorb,ival,:)*dis_factors(ival,:)*phs(:,iR))
            else
              ! B-site
              wftrial(iorb,iw,iR)=sum(evec(iorb,icnd,:)*dis_factors(icnd,:)*phs(:,iR))
            end if 
          else
            ! bottom layer
            if (ispec.eq.1) then
              ! A-site
              wftrial(iorb,iw,iR)=sum(evec(iorb,icnd,:)*dis_factors(icnd,:)*phs(:,iR))
            else                                                           
              ! B-site                                                     
              wftrial(iorb,iw,iR)=sum(evec(iorb,ival,:)*dis_factors(ival,:)*phs(:,iR))
            end if 
          end if
        end do
      end do
    end do
    deallocate(phs,dis_factors)
    wftrial=wftrial/dble(kgrid%npt)
    sigma=0.6_dp*sqrt(dot_product(pars%avec(1,:),pars%avec(1,:)))
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
    call generate_amn_overlap(tbmodel,pars,kgrid,eval,evec,tbmodel%rgrid%npt,wftrial)
    deallocate(wftrial)
  end if
  call generate_dmn_orb(tbmodel,proj,sym,pars,kgrid,evec,.false.,wws)
else if (trim(adjustl(pars%wannier_proj_mode)).eq.'input_file') then

   inquire(file=trim(adjustl(pars%seedname))//'.amn',exist=exs)
   if (.not.exs) then
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
           wftrial(jorb,ipro,tbmodel%rgrid%ip0)=tbmodel%wbase%wws_full(.false.,sym%car(:,:,1),l1,m1,l2,m2,x1,z1,x2,z2)
         end do
       end do
     end do
     call generate_amn_overlap(tbmodel,pars,kgrid,eval,evec,tbmodel%rgrid%npt,wftrial)
     deallocate(wftrial)
   end if
   call generate_dmn_orb(tbmodel,proj,sym,pars,kgrid,evec,.false.,wws)

else if (trim(adjustl(pars%wannier_proj_mode)).eq.'wannier_file') then

   inquire(file=trim(adjustl(pars%seedname))//'.amn',exist=exs)
   if (.not.exs) then
     allocate(wftrial(tbmodel%norb_TB,pars%proj%norb,tbmodel%rgrid%npt))
     wftrial(:,:,:)=0._dp
     call read_wfmloc(pars,tbmodel,kgrid,eval,evec,wftrial)
     call generate_amn_overlap(tbmodel,pars,kgrid,eval,evec,tbmodel%rgrid%npt,wftrial)
     deallocate(wftrial)
   end if
   call generate_dmn_orb(tbmodel,proj,sym,pars,kgrid,evec,.false.,wws)

else if (trim(adjustl(pars%wannier_proj_mode)).eq.'wftrial_file') then

   inquire(file=trim(adjustl(pars%seedname))//'.amn',exist=exs)
   if (.not.exs) then
     allocate(wftrial(tbmodel%norb_TB,pars%proj%norb,tbmodel%rgrid%npt))
     wftrial(:,:,:)=0._dp
     call read_wf_universe(tbmodel,pars,tbmodel%rgrid%npt,wftrial,'wftrial','')
     call generate_amn_overlap(tbmodel,pars,kgrid,eval,evec,tbmodel%rgrid%npt,wftrial)
     deallocate(wftrial)
   end if
   call generate_dmn_orb(tbmodel,proj,sym,pars,kgrid,evec,.false.,wws)

else if (trim(adjustl(pars%wannier_proj_mode)).eq.'real_space') then

   inquire(file=trim(adjustl(pars%seedname))//'.amn',exist=exs)
   if (.not.exs) then
     allocate(wftrial(tbmodel%norb_TB,pars%proj%norb,tbmodel%rgrid%npt))
     wftrial(:,:,:)=0._dp
     sigma=2.0_dp*sqrt(dot_product(pars%avec(1,:),pars%avec(1,:)))
     call real_space_wftrial(tbmodel,proj,wftrial,sigma)
     call generate_amn_overlap(tbmodel,pars,kgrid,eval,evec,tbmodel%rgrid%npt,wftrial)
     deallocate(wftrial)
   end if
   call generate_dmn_orb(tbmodel,proj,sym,pars,kgrid,evec,.false.,wws)

else
  call throw("wannier_interface%generate_trial_wavefunctions()","unknown projection option")
end if
call generate_mmn_overlap(THIS,tbmodel,pars,kgrid,evec)
call generate_tmn_overlap(tbmodel,pars,kgrid,eval,evec)
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
integer ikp,ist,jst,nn,nsym
logical exs
real(dp) err,t1,t2
real(dp) v1(NDIM),v2(NDIM)
complex(dp), allocatable :: phs(:,:)
complex(dp), allocatable :: wwst(:,:)
complex(dp), allocatable :: wf_t(:,:)
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
          wws(jw,iw,isym)=base%wws(.false.,sym%car(:,:,isym),iw,jw)
       end do
    end do
    do iw=1,base%norb
       err=abs((sum(abs(wws(:,iw,isym))**2)+sum(abs(wws(iw,:,isym))**2))*.5_dp-1._dp)
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
  if (pars%trev) then
    write (1001,"(4i9)") pars%nstates, sym%nsym+1, kgrid%nir, kgrid%npt
  else
    write (1001,"(4i9)") pars%nstates, sym%nsym, kgrid%nir, kgrid%npt
  end if
  write (1001,*)
  write (1001,"(10i9)") kgrid%ik2ir(1:kgrid%npt)
  write (1001,*)
  write (1001,"(10i9)") kgrid%ir2ik(1:kgrid%nir)
  do ir=1,kgrid%nir
     write (1001,*)
     write (1001,"(10i9)") kgrid%iks2k(kgrid%ir2ik(ir),:)
  enddo
end if
allocate(phs(base%norb,base%norb),wwst(base%norb,base%norb))
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
  phs=0._dp
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
     wwst=matmul(phs,wws(:,:,isym))
     ! unitarize from the left
     !call unitarize(-1,base%norb,base%norb,wws(:,:,isym))
     if (mp_mpi) then
        WRITE (1001,*)
        WRITE (1001,"(1p,(' (',e18.10,',',e18.10,')'))") wwst
     end if
  end do
  phs=0._dp
  if (pars%trev) then
     do iw=1,base%norb
        do jw=1,base%norb
           if(base%orb_icio(jw,1).ne.base%orb_icio(iw,1)) cycle
           phs(jw,iw)=base%wws(.true.,sym%car(:,:,1),iw,jw)
        end do
     end do
     if (mp_mpi) then
        WRITE (1001,*)
        WRITE (1001,"(1p,(' (',e18.10,',',e18.10,')'))") phs
     end if
  end if
end do
if (mp_mpi) then
  write(1001,*)
  flush(1001)
  write(*,'(/)')
  write(*,'(a,i8)') '  DMN(d_matrix_band) [could be not ordered]: nir = ',kgrid%nir
end if
if (pars%trev) then
  nsym=sym%nsym+1
else
  nsym=sym%nsym
end if
allocate(ovlp(pars%nstates,pars%nstates,nsym,kgrid%nir))
allocate(wf_t(tbmodel%norb_TB,pars%nstates))
ovlp=0._dp
#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
do ir=1,kgrid%nir
  if (mod(ir-1,np_mpi).ne.lp_mpi) cycle
  WRITE (*,*) "ir: ",ir
  ik=kgrid%ir2ik(ir)
  do isym=1,sym%nsym
    ikp=kgrid%iks2k(ik,isym)
    !$OMP PARALLEL DEFAULT(SHARED)&
    !$OMP PRIVATE(ist)
    !$OMP DO
    do jst=1,pars%nstates
      call tbmodel%bloch_wf_transform(kgrid,ik,sym,isym,wf_t(:,jst),evec(:,jst,ik))
      do ist=1,pars%nstates
         ovlp(ist,jst,isym,ir)=dot_product(evec(:,ist,ikp),wf_t(:,jst))
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
  end do
  if (pars%trev) then
    ikp=kgrid%iks2k(ik,sym%nsym+1)
    do jst=1,pars%nstates
      do ist=1,pars%nstates
         ovlp(ist,jst,sym%nsym+1,ir)=dot_product(evec(:,ist,ikp),conjg(evec(:,jst,ik)))
         !ovlp(ist,jst,sym%nsym+1,ir)=dot_product(evec(:,ist,ik),conjg(evec(:,jst,ikp)))
         !write(*,*) ir,ikp, ist,jst, dble(ovlp(ist,jst,sym%nsym+1,ir)),aimag(ovlp(ist,jst,sym%nsym+1,ir))
      end do
    end do
  end if
end do
#ifdef MPI
  ! it is nsym, not sym%nsym here
  nn=pars%nstates*pars%nstates*nsym*kgrid%nir
  call mpi_allreduce(mpi_in_place,ovlp,nn,mpi_double_complex,mpi_sum, &
   mpi_com,mpi_err)
#endif
if (mp_mpi) then
  do ir=1,kgrid%nir
    do isym=1,nsym
      ! unitarize from the right
      !call unitarize(1,pars%nstates,pars%nstates,ovlp(:,:,isym,ir))
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
deallocate(phs,wwst,ovlp)
deallocate(wf_t)
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

subroutine generate_amn_overlap(tbmodel,pars,kgrid,eval,evec,nr,wftrial)
class(CLtb), intent(in) :: tbmodel
class(CLpars), intent(in) :: pars
class(GRID), intent(in) :: kgrid
real(dp), intent(in) :: eval(pars%nstates,kgrid%npt)
complex(dp), intent(in) :: evec(tbmodel%norb_TB,pars%nstates,kgrid%npt)
integer, intent(in) :: nr
complex(dp), intent(in) :: wftrial(tbmodel%norb_TB,pars%proj%norb,nr)
! local
logical exs
integer ik,iwan,ist,iR,iorb,m
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
allocate(amn(pars%nstates,pars%proj%norb,kgrid%npt))
amn=0._dp
!$OMP PARALLEL DEFAULT (SHARED)&
!$OMP PRIVATE(iR,iwan,iorb,t1,z1)
!$OMP DO
do ik=1,kgrid%npt
  do iR=1,tbmodel%rgrid%npt
    if ( sum(abs(wftrial(:,:,iR))).lt.epslat) cycle 
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
if (pars%use_weights_amn) then
   do ik=1,kgrid%npt
      do m=1,pars%nstates
         if (eval(m,ik).gt.pars%dis_frozen(1) .and. eval(m,ik).lt.pars%dis_frozen(2)) then
           t1=1._dp
         else if (eval(m,ik).lt.pars%dis_frozen(1)) then
           t1=gauss(pars%dis_frozen(1)-eval(m,ik),pars%gauss_sigma)
         else if (eval(m,ik).gt.pars%dis_frozen(2)) then
           t1=gauss(eval(m,ik)-pars%dis_frozen(2),pars%gauss_sigma)
         end if
         !write(12,*) eval(m,ik),t1
         amn(m,:,ik)=amn(m,:,ik)*t1
      end do
   end do
end if
if (mp_mpi) then
  open(50,file=trim(adjustl(pars%seedname))//'.amn',action='write')
  write(50,*) '# '//trim(adjustl(pars%seedname))//' file '
  write(50,*) pars%nstates,kgrid%npt,pars%proj%norb
end if
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
        mmn(:,nn,innk)=mmn(:,nn,innk)+conjg(evec(iorb,:,ik))*evec(iorb,nn,jk)*z1!*pwave_ovlp(dc)
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

subroutine generate_tmn_overlap(tbmodel,pars,kgrid,eval,evec)
class(CLtb), intent(in) :: tbmodel
class(CLpars), intent(in) :: pars
class(GRID), intent(in) :: kgrid
real(dp), intent(in) :: eval(pars%nstates,kgrid%npt)
complex(dp), intent(in) :: evec(tbmodel%norb_TB,pars%nstates,kgrid%npt)
! local
logical exs
integer ik,jk,mm,nn,iorb,ir
real(dp) dc
real(dp) vq(NDIM),vc(NDIM)
complex(dp), allocatable :: tmn(:,:)
! M_mn(k)=<psi_mk|psi*_n{-k}>
inquire(file=trim(adjustl(pars%seedname))//'.tmn',exist=exs)
if (exs) then
  !call info("CLwan%generate_mmn_overlap","skipping "//trim(adjustl(pars%seedname))//".mmn creation")
  !return
else
  call info("CLwan%generate_tmn_overlap","generating "//trim(adjustl(pars%seedname))//".tmn file")
end if
if (mp_mpi) then
  open(50,file=trim(adjustl(pars%seedname))//'.tmn',action='write')
  write (50,*) '# '//trim(adjustl(pars%seedname))//'.tmn file'
  write (50,"(4i9)") pars%nstates, kgrid%nirT, kgrid%npt
  write (50,*)
  write (50,"(10i9)") kgrid%ikT2ir(1:kgrid%npt)
  write (50,*)
  write (50,"(10i9)") kgrid%irT2ik(1:kgrid%nirT)
  write (50,*)
  write (50,"(10i9)") kgrid%ikT2k(1:kgrid%npt)
end if
allocate(tmn(pars%nstates,pars%nstates))
do ir=1,kgrid%nirT
  ik=kgrid%irT2ik(ir)
  jk=kgrid%ikT2k(ik)
  vq=-kgrid%vpl(ik)-kgrid%vpl(ik)
  vc=matmul(vq,kgrid%vecs)
  dc=sqrt(dot_product(vc,vc))
  tmn=0._dp
  do iorb=1,tbmodel%norb_TB
    do nn=1,pars%nstates
      tmn(:,nn)=tmn(:,nn)+conjg(evec(iorb,:,ik))*conjg(evec(iorb,nn,jk))!*pwave_ovlp(dc)
    end do
  end do
  do mm=1,pars%nstates
      if (dble(tmn(mm,mm))<0.9) then
        write(*,'(i4,"(",i4,")",3I4,20F10.4)') ir,kgrid%nirT,ik,jk,mm,eval(mm,ik),eval(mm,jk)
        write(*,'(6F10.4,"|",2F10.4)') kgrid%vpl(ik),kgrid%vpl(jk),tmn(mm,mm)
      end if
  end do
  if (mp_mpi) then
     WRITE (50,*)
     WRITE (50,"(1p,(' (',e18.10,',',e18.10,')'))") tmn
  end if
end do
if (mp_mpi) close(50)
deallocate(tmn)
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

subroutine read_wfmloc(pars,tbmodel,kgrid,eval,evec,wfmloc)
class(CLpars), intent(in) :: pars
class(CLtb), intent(in) :: tbmodel
class(GRID), intent(inout) :: kgrid
real(dp), intent(in) :: eval(pars%nstates,kgrid%npt)
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
   tbmodel%norb_TB,pars%dis_win,udis,umat,vkl,vrl,eval,evec,wfmloc)
deallocate(vkl,vrl,umat,udis)
return
end subroutine

subroutine symmetrize_tbfile(pars,sym)
class(CLpars), intent(in) :: pars
class(CLsym), intent(in) :: sym
complex(dp), allocatable :: wws(:,:,:)
complex(dp), allocatable :: wws_trev(:,:)
complex(dp), allocatable :: hams(:,:,:)
complex(dp), allocatable :: hamt(:,:,:)
type(wbase)  proj
type(GRID) rgrid
type(CLwan) wan
type(GRID) kgrid
integer iter
integer isym,jsym,iR,jR
integer n1,n2,m1,m2
integer nc1,nc2,mc1,mc2
integer jc1,jc2
integer iv(NDIM+1)
real(dp) RI(NDIM),RJ(NDIM),SJ(NDIM,NDIM)
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
! init k-grid
iR=maxval(wan%rgrid%ngrid)
iv(1:2)=iR
iv(3)=1
call kgrid%init(iv(1:3),pars%bvec,centered_kgrid,.true.)
call kgrid%sym_init(pars%trev,sym)
allocate(wws(proj%norb,proj%norb,sym%nsym))
allocate(hams(proj%norb,proj%norb,wan%rgrid%npt))
allocate(hamt(proj%norb,proj%norb,wan%rgrid%npt))
call init_wws(sym,proj,wws)
do iter=1,pars%niter_symmetrize
  hams=0._dp
  do isym=1,sym%nsym
    jsym=sym%inv(isym)
    SJ=dble(sym%lat(:,:,jsym))
    do jR=1,wan%rgrid%npt
      RJ=wan%rgrid%vpl(jR)
      do n1=1,proj%norb
        do n2=1,proj%norb
          nc1=proj%orb_icio(n1,1)
          nc2=proj%orb_icio(n2,1)
          jc1=proj%ics2c(nc1,jsym)
          jc2=proj%ics2c(nc2,jsym)
          RI=matmul(SJ,RJ+proj%centers(:,nc2)-proj%centers(:,nc1))-(proj%centers(:,jc2)-proj%centers(:,jc1))
          iv=wan%rgrid%find(RI)
          if (iv(NDIM+1).lt.0.or.sum(abs(iv(1:NDIM))).ne.0) cycle
          iR=iv(NDIM+1)
          do m1=1,proj%norb
            do m2=1,proj%norb
              mc1=proj%orb_icio(m1,1)
              mc2=proj%orb_icio(m2,1)
              if (proj%ics2c(mc1,isym).ne.nc1) cycle
              if (proj%ics2c(mc2,isym).ne.nc2) cycle
              hams(n1,n2,jR)=hams(n1,n2,jR)+wws(n1,m1,isym)*wan%hame(m1,m2,iR)*wws(m2,n2,jsym)
            end do
          end do
        end do
      end do
    end do
  end do
  if (pars%trev) then
     allocate(wws_trev(proj%norb,proj%norb))
     wws_trev=0._dp
     do n1=1,proj%norb
        do n2=1,proj%norb
           if(proj%orb_icio(n2,1).ne.proj%orb_icio(n1,1)) cycle
           wws_trev(n2,n1)=proj%wws(.true.,sym%car(:,:,1),n1,n2)
        end do
     end do
     do n1=1,proj%norb
       do n2=1,proj%norb
         nc1=proj%orb_icio(n1,1)
         nc2=proj%orb_icio(n2,1)
         do m1=1,proj%norb
           do m2=1,proj%norb
             mc1=proj%orb_icio(m1,1)
             mc2=proj%orb_icio(m2,1)
             if (mc1.ne.nc1) cycle
             if (mc2.ne.nc2) cycle
             hams(n1,n2,:)=hams(n1,n2,:)+wws_trev(n1,m1)*wan%hame(m1,m2,:)*wws_trev(m2,n2)
           end do
         end do
       end do
     end do
     deallocate(wws_trev) 
     hams=hams/dble(sym%nsym+1)
  else
    hams=hams/dble(sym%nsym)
  endif
  wan%hame=hams
end do
call wan%write_tb_file('hamwan',pars,proj%norb)
deallocate(hams,hamt,wws)
end subroutine

subroutine symmetrize_trev(nwan,wan,kgrid,ham,hams)
integer, intent(in) :: nwan
type(CLwan), intent(in) :: wan
type(GRID), intent(in) :: kgrid
complex(dp), intent(in) :: ham(nwan,nwan,wan%rgrid%npt)
complex(dp), intent(out) :: hams(nwan,nwan,wan%rgrid%npt)
complex(dp), allocatable :: hamk(:,:,:),hamsk(:,:,:)
integer, parameter :: sgn=1
real(dp), parameter :: mix=0.2_dp
integer ir,ik,jk,iw,jw
real(dp) t1
complex(dp) z1
complex(dp), dimension(:,:), allocatable :: A,B,U,SA,SB
real(dp), allocatable :: s1(:),s2(:)
allocate(A(nwan,nwan))
allocate(B(nwan,nwan))
allocate(U(nwan,nwan))
allocate(SA(nwan,nwan))
allocate(SB(nwan,nwan))
allocate(s1(nwan))
allocate(s2(nwan))
allocate(hamk(nwan,nwan,kgrid%npt))
allocate(hamsk(nwan,nwan,kgrid%npt))
! FT TB Hamiltonian
hamk=0._dp
do ik=1,kgrid%npt
  do iR=1,wan%rgrid%npt
    t1=twopi*dot_product(kgrid%vpl(ik),wan%rgrid%vpl(iR))
    z1=cmplx(cos(t1),sgn*sin(t1),kind=dp)
    hamk(:,:,ik)=hamk(:,:,ik)+ham(:,:,iR)*z1
  end do
end do
! find unique TR Unitary operator U, such that H_k=U H*_{-k} U^\dagger valid approximately for all k-points
! re-write in terms of matrixes A and B  : A=UBU^\dagger
! extracting square roots : sqrt(A) sqrt(A)^\dagger = U sqrt(B) sqrt(B)^\dagger U^\dagger
! then U = sqrt(A) sqrt(B)^-1
U=0._dp
! extract matrix U only from the first point(GAMMA)
do ir=1,kgrid%nirT
    ik=kgrid%irT2ik(ir)
    jk=kgrid%ikT2k(ik)
    ! find square root of hamk(:,:,ik) and conjg(hamk(:,:,jk))
    A=hamk(:,:,ik)
    B=conjg(hamk(:,:,jk))
    do iw=1,nwan
      A(iw,iw)=A(iw,iw)+100._dp ! makes eigenvalues to be positive
      B(iw,iw)=B(iw,iw)+100._dp
    end do
    call eigenv_problem(nwan,A,s1)
    call eigenv_problem(nwan,B,s2)
    ! square roots
    do iw=1,nwan
      SA(iw,:)=sqrt(s1(iw))*conjg(A(:,iw))
      SB(iw,:)=sqrt(s2(iw))*conjg(B(:,iw)) 
    end do
    A=matmul(A,SA)
    B=matmul(B,SB)
    call utility_zgetri(B)
    ! add to the unitary transformation
    SA=matmul(A,B)
    U=U+matmul(A,B)
end do
U=U/dble(kgrid%nirT)
! now we need to make SVD to make it unitary
A=U
call utility_zgesvd(A,SA,s1,SB)
do iw=1,nwan
  SA(iw,:)=1._dp/abs(s1(iw))*SB(iw,:)
end do
A=matmul(conjg(transpose(SB)),SA)
U=matmul(U,A)
write(*,*) U
stop
! Apply TR in reciprocal space
do ir=1,kgrid%nirT
  ik=kgrid%irT2ik(ir)
  jk=kgrid%ikT2k(ik)
  call utility_zgemmm(U,'N',conjg(hamk(:,:,jk)), 'N', U, 'C', hamsk(:,:,ik))
  call utility_zgemmm(U,'N',conjg(hamsk(:,:,ik)), 'N', U, 'C', hamsk(:,:,jk))
end do
! FT back to real space 
hams=0._dp
do iR=1,wan%rgrid%npt
  do ik=1,kgrid%npt
    t1=twopi*dot_product(kgrid%vpl(ik),wan%rgrid%vpl(iR))
    z1=cmplx(cos(t1),-sgn*sin(t1),kind=dp)
    hams(:,:,iR)=hams(:,:,iR)+hamsk(:,:,ik)*z1
  end do
end do
hams=hams/dble(kgrid%npt)
deallocate(A,B,U,SA,SB,s1,s2)
deallocate(hamk,hamsk)
end subroutine symmetrize_trev



subroutine write_hubbardu(pars,sym)
class(CLpars), intent(in) :: pars
class(CLsym), intent(in) :: sym
logical exst
character(len=4) snum
integer iv(NDIM+1)
integer iR_sphere,jR_sphere,IRR_POINT
integer norb,npt_sphere
integer mc1p,mc2p
integer isym,jsym
integer nsize
integer n1,n2,n3,n4
integer m1,m2,m3,m4
integer nc1,nc2,nc3,nc4
integer mc1,mc2,mc3,mc4
real(dp) t1,t2
real(dp) RR(NDIM)
complex(dp), allocatable :: wws(:,:,:)
complex(dp), allocatable :: UH(:,:,:,:,:)
complex(dp), allocatable :: UHS(:,:,:,:,:)
integer, allocatable     :: irr2cR(:,:)
integer, allocatable     :: cR2irr(:,:,:,:)
integer                  :: nir
type(wbase)  proj
type(GRID) rgrid
!real(dp) ddot
if (pars%proj%norb.le.0) then
   call throw("wannier_interface%write_hubbardu()",&
              "apparently 'projections' block was not specified, wannier projections not found")
else
  if (pars%nstates.lt.pars%proj%norb) then
    call throw("wannier_interface%write_hubbardu()",&
               "number of bands is less than the number of requested projections")
  end if
end if
call proj%init(pars,pars%proj%ncenters,pars%proj%norb,pars%proj%norb_ic,&
                   pars%proj%lmr,pars%proj%waxis,pars%proj%centers)
call proj%init_smap(sym,pars)
call rgrid%init(pars%ngrid,pars%avec,.true.,.false.)
call rgrid%init_sphere(pars)
allocate(wws(proj%norb,proj%norb,sym%nsym))
call init_wws(sym,proj,wws)
nsize=proj%norb
allocate(UH(nsize,nsize,nsize,nsize,rgrid%npt_sphere))
UH=0._dp
if (.not.pars%HubU_diagonal) then
  allocate(irr2cR(3,proj%ncenters*proj%ncenters*rgrid%npt_sphere))
  allocate(cR2irr(2,proj%ncenters,proj%ncenters,rgrid%npt_sphere))
  call find_irrcR_sphere(proj,sym,rgrid,irr2cR,cR2irr,nir)
  do IRR_POINT=1,nir
    jR_sphere=irr2cR(3,IRR_POINT)
    write(snum,'(I4.4)') IRR_POINT
    inquire(file='_UH/UH'//snum,exist=exst)
    if (.not.exst) call throw("wannier_interface%write_hubbardu()","file '"//'_UH/UH'//snum//"' not found (symmetric case)")
    open(140,file='_UH/UH'//snum,action="read")
    read(140,*) npt_sphere,norb
    if (npt_sphere.ne.rgrid%npt_sphere) call throw("wannier_interface%write_hubbardu()","npt_sphere is not what expected")
    do n1=1,1000000000
      read(140,*,end=140) iR_sphere,m1,m2,m3,m4,t1,t2
      if (jR_sphere.ne.iR_sphere) call throw("wannier_interface%write_hubbardu()","jR_irr is not what expected")
      UH(m1,m2,m3,m4,jR_sphere)=cmplx(t1,t2,kind=8)
    end do
    140 continue
    close(140)
  end do
  ! reconstruct full UH from irreducible wedge
  allocate(UHS(nsize,nsize,nsize,nsize,rgrid%npt_sphere))
  UHS=0._dp
  do jR_sphere=1,rgrid%npt_sphere
    do n1=1,proj%norb
      do n2=1,proj%norb
        do n3=1,proj%norb
          do n4=1,proj%norb
            nc1=proj%orb_icio(n1,1)
            nc2=proj%orb_icio(n2,1)
            nc3=proj%orb_icio(n3,1)
            nc4=proj%orb_icio(n4,1)
            ! this gives the diagonal (IN SITE INDEX) Hubbard 
            if (nc1.ne.nc4) cycle
            if (nc2.ne.nc3) cycle
            irr_point=cR2irr(1,nc1,nc2,jR_sphere)
            if (irr_point.lt.0) cycle
            isym=cR2irr(2,nc1,nc2,jR_sphere)
            jsym=sym%inv(isym)
            mc1p=irr2cR(1,irr_point)
            mc2p=irr2cR(2,irr_point)
            iR_sphere=irr2cR(3,irr_point)
            if (iR_sphere.lt.0) cycle 
            RR=rgrid%vpl_sphere(jR_sphere)
            do m1=1,proj%norb
              do m2=1,proj%norb
                do m3=1,proj%norb
                  do m4=1,proj%norb
                    mc1=proj%orb_icio(m1,1)
                    mc2=proj%orb_icio(m2,1)
                    mc3=proj%orb_icio(m3,1)
                    mc4=proj%orb_icio(m4,1)
                    if (mc1.ne.mc1p) cycle
                    if (mc2.ne.mc2p) cycle
                    if (proj%ics2c(mc1,isym).ne.nc1) cycle
                    if (proj%ics2c(mc2,isym).ne.nc2) cycle
                    if (proj%ics2c(mc3,isym).ne.nc3) cycle
                    if (proj%ics2c(mc4,isym).ne.nc4) cycle
                    UHS(n1,n2,n3,n4,jR_sphere)=UHS(n1,n2,n3,n4,jR_sphere)+&
                         wws(n1,m1,isym)*wws(n2,m2,isym)*UH(m1,m2,m3,m4,iR_sphere)*wws(m3,n3,jsym)*wws(m4,n4,jsym)
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do
  UH=UHS
  deallocate(UHS)
else
  do jR_sphere=1,rgrid%npt_sphere
    write(snum,'(I4.4)') jR_sphere
    inquire(file='_UH/UH'//snum,exist=exst)
    if (.not.exst) call throw("wannier_interface%write_hubbardu()","file '"//'_UH/UH'//snum//"' not found (symmetric case)")
    open(141,file='_UH/UH'//snum,action="read")
    read(141,*) npt_sphere,norb
    if (npt_sphere.ne.rgrid%npt_sphere) call throw("wannier_interface%write_hubbardu()","npt_sphere is not what expected")
    do n1=1,1000000000
      read(141,*,end=141) iR_sphere,m1,m2,m3,m4,t1,t2
      if (jR_sphere.ne.iR_sphere) call throw("wannier_interface%write_hubbardu()","jR_irr is not what expected")
      UH(m1,m2,m3,m4,jR_sphere)=cmplx(t1,t2,kind=8)
    end do
    141 continue
    close(141)
  end do
end if
! take Hermitian part
do jR_sphere=1,rgrid%npt_sphere
  do n1=1,proj%norb ; do n2=1,proj%norb ; do n3=1,proj%norb ; do n4=1,proj%norb
    iv=rgrid%find(rgrid%vpl_sphere(jR_sphere))
    if (iv(NDIM+1).lt.0.or.sum(abs(iv(1:NDIM))).ne.0) cycle
    iR_sphere=rgrid%homo_to_sphere(iv(NDIM+1))
    if (iR_sphere.le.0) cycle
    UH(n1,n2,n3,n4,jR_sphere)=0.5_dp*(UH(n1,n2,n3,n4,jR_sphere)+conjg(UH(n4,n3,n2,n1,iR_sphere)))
    UH(n4,n3,n2,n1,iR_sphere)=conjg(UH(n1,n2,n3,n4,jR_sphere))
  end do ; end do ; end do ; end do
end do
! write into files
if (mp_mpi) open(140,file="UR.dat",action="write")
if (mp_mpi) write(140,*) rgrid%npt_sphere,proj%norb
do jR_sphere=1,rgrid%npt_sphere
  if (mp_mpi) write(140,*) nint(rgrid%vpl_sphere(jR_sphere))
  do n1=1,proj%norb
    do n2=1,proj%norb
      if (mp_mpi) write(140,'(2I6,2G18.10)') n1,n2,UH(n1,n2,n2,n1,jR_sphere)
    end do
  end do
  if (mp_mpi) write(140,*)
end do
if (mp_mpi) close(140)
if (mp_mpi) open(140,file="_UH/UH.dat",action="write")
do jR_sphere=1,rgrid%npt_sphere
  do n1=1,proj%norb ; do n2=1,proj%norb ; do n3=1,proj%norb ;do n4=1,proj%norb
    if (mp_mpi) write(140,'(2G18.10)') UH(n1,n2,n3,n4,jR_sphere)
  end do ; end do ; end do ; end do
end do
if (mp_mpi) close(140)
deallocate(wws)
deallocate(UH)
end subroutine

subroutine symmetrize_hubbardu(pars,sym)
class(CLpars), intent(in) :: pars
class(CLsym), intent(in) :: sym
integer iR_sphere,jR_sphere
integer isym,jsym
integer nsize,iter
integer jc1,jc2
integer n1,n2,n3,n4
integer m1,m2,m3,m4
integer nc1,nc2,nc3,nc4
integer mc1,mc2,mc3,mc4
integer iv(NDIM+1)
real(dp) t1,t2
real(dp) RI(NDIM),RJ(NDIM)
real(dp) SJ(NDIM,NDIM)
complex(dp), allocatable :: wws(:,:,:),wwst(:,:)
complex(dp), allocatable :: UH(:,:,:,:,:)
complex(dp), allocatable :: UHS(:,:,:,:,:)
type(wbase)  proj
type(GRID) rgrid
!real(dp) ddot
if (pars%proj%norb.le.0) then
   call throw("wannier_interface%write_hubbardu()",&
              "apparently 'projections' block was not specified, wannier projections not found")
else
  if (pars%nstates.lt.pars%proj%norb) then
    call throw("wannier_interface%write_hubbardu()",&
               "number of bands is less than the number of requested projections")
  end if
end if
if (pars%HubU_diagonal) call throw("wannier_interface%write_hubbardu()",&
                          "symmetrization does not work with pars%HubU_diagonal")

call proj%init(pars,pars%proj%ncenters,pars%proj%norb,pars%proj%norb_ic,&
                   pars%proj%lmr,pars%proj%waxis,pars%proj%centers)
call proj%init_smap(sym,pars)
call rgrid%init(pars%ngrid,pars%avec,.true.,.false.)
call rgrid%init_sphere(pars)
allocate(wws(proj%norb,proj%norb,sym%nsym))
allocate(wwst(proj%norb,proj%norb))
call init_wws(sym,proj,wws)
nsize=proj%norb
allocate(UH(nsize,nsize,nsize,nsize,rgrid%npt_sphere))
allocate(UHS(nsize,nsize,nsize,nsize,rgrid%npt_sphere))
UH=0._dp
if (mp_mpi) open(140,file="_UH/UH.dat",action="read")
do jR_sphere=1,rgrid%npt_sphere
  do n1=1,proj%norb ; do n2=1,proj%norb ; do n3=1,proj%norb ; do n4=1,proj%norb
    if (mp_mpi) read(140,*) t1,t2
    UH(n1,n2,n3,n4,jR_sphere)=cmplx(t1,t2,kind=dp)
  end do ; end do ; end do ; end do
end do
if (mp_mpi) close(140)
! take Hermitian part
do jR_sphere=1,rgrid%npt_sphere
  do n1=1,proj%norb ; do n2=1,proj%norb ; do n3=1,proj%norb ; do n4=1,proj%norb
    iv=rgrid%find(rgrid%vpl_sphere(jR_sphere))
    if (iv(NDIM+1).lt.0.or.sum(abs(iv(1:NDIM))).ne.0) cycle
    iR_sphere=rgrid%homo_to_sphere(iv(NDIM+1))
    if (iR_sphere.le.0) cycle
    UH(n1,n2,n3,n4,jR_sphere)=0.5_dp*(UH(n1,n2,n3,n4,jR_sphere)+conjg(UH(n4,n3,n2,n1,iR_sphere)))
    UH(n4,n3,n2,n1,iR_sphere)=conjg(UH(n1,n2,n3,n4,jR_sphere))
  end do ; end do ; end do ; end do
end do
do iter=1,pars%niter_symmetrize
  UHS=0._dp
  do isym=1,sym%nsym
    jsym=sym%inv(isym)
    SJ=dble(sym%lat(:,:,jsym))
    do jR_sphere=1,rgrid%npt_sphere
      RJ=rgrid%vpl_sphere(jR_sphere)
      do n1=1,proj%norb
        do n2=1,proj%norb
          do n3=1,proj%norb
            do n4=1,proj%norb
              nc1=proj%orb_icio(n1,1)
              nc2=proj%orb_icio(n2,1)
              nc3=proj%orb_icio(n3,1)
              nc4=proj%orb_icio(n4,1)
              ! this gives the diagonal (IN SITE INDEX) Hubbard 
              if (nc1.ne.nc4) cycle
              if (nc2.ne.nc3) cycle
              jc1=proj%ics2c(nc1,jsym)
              jc2=proj%ics2c(nc2,jsym)
              RI=matmul(SJ,RJ+proj%centers(:,nc2)-proj%centers(:,nc1))-(proj%centers(:,jc2)-proj%centers(:,jc1))
              iv=rgrid%find(RI)
              if (iv(NDIM+1).lt.0.or.sum(abs(iv(1:NDIM))).ne.0) cycle
              iR_sphere=rgrid%homo_to_sphere(iv(NDIM+1))
              if (iR_sphere.le.0) cycle
              do m1=1,proj%norb
                do m2=1,proj%norb
                  do m3=1,proj%norb
                    do m4=1,proj%norb
                      mc1=proj%orb_icio(m1,1)
                      mc2=proj%orb_icio(m2,1)
                      mc3=proj%orb_icio(m3,1)
                      mc4=proj%orb_icio(m4,1)
                      if (mc1.ne.mc4) cycle
                      if (mc2.ne.mc3) cycle
                      if (proj%ics2c(mc1,isym).ne.nc1) cycle
                      if (proj%ics2c(mc2,isym).ne.nc2) cycle
                      if (proj%ics2c(mc3,isym).ne.nc3) cycle
                      if (proj%ics2c(mc4,isym).ne.nc4) cycle
                      UHS(n1,n2,n3,n4,jR_sphere)=UHS(n1,n2,n3,n4,jR_sphere)+&
                           wws(n1,m1,isym)*wws(n2,m2,isym)*UH(m1,m2,m3,m4,iR_sphere)*wws(m3,n3,jsym)*wws(m4,n4,jsym)
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do
  UH=UHS/dble(sym%nsym)
end do
if (mp_mpi) open(140,file="URS.dat",action="write")
if (mp_mpi) write(140,*) rgrid%npt_sphere,proj%norb
do jR_sphere=1,rgrid%npt_sphere
  if (mp_mpi) write(140,*) nint(rgrid%vpl_sphere(jR_sphere))
  do n1=1,proj%norb
    do n2=1,proj%norb
      if (mp_mpi) write(140,'(2I6,2G18.10)') n1,n2,UH(n1,n2,n2,n1,jR_sphere)
    end do
  end do
  if (mp_mpi) write(140,*)
end do
if (mp_mpi) close(140)
if (mp_mpi) open(140,file="_UH/UHS.dat",action="write")
do jR_sphere=1,rgrid%npt_sphere
  do n1=1,proj%norb ; do n2=1,proj%norb ; do n3=1,proj%norb ; do n4=1,proj%norb
    if (mp_mpi) write(140,'(2G18.10)') UH(n1,n2,n3,n4,jR_sphere)
  end do ; end do ; end do ; end do
end do
if (mp_mpi) close(140)
deallocate(wws)
deallocate(UH,UHS)
end subroutine

subroutine find_irrcR_sphere(proj,sym,rgrid,irr2cR,cR2irr,nir)
class(wbase), intent(in) :: proj
class(CLsym), intent(in) :: sym
class(GRID), intent(in)  :: rgrid
integer, intent(out)     :: irr2cR(3,proj%ncenters*proj%ncenters*rgrid%npt_sphere)
integer, intent(out)     :: cR2irr(2,proj%ncenters,proj%ncenters,rgrid%npt_sphere)
integer, intent(out)     :: nir
integer isym
integer ic1,ic2,jc1,jc2
integer iR_sphere,jR_sphere
logical found(proj%ncenters,proj%ncenters,rgrid%npt_sphere)
real(dp) RI(NDIM),RJ(NDIM)
real(dp) SI(NDIM,NDIM)
integer iv(NDIM+1)
nir=0
irr2cR=-1
found=.false.
do jR_sphere=1,rgrid%npt_sphere
  RJ=rgrid%vpl_sphere(jR_sphere)
  do jc1=1,proj%ncenters
    do jc2=1,proj%ncenters
      if (found(jc1,jc2,jR_sphere)) cycle
      found(jc1,jc2,jR_sphere)=.true.
      nir=nir+1
      irr2cR(1,nir)=jc1
      irr2cR(2,nir)=jc2
      irr2cR(3,nir)=jR_sphere
      cR2irr(1,jc1,jc2,jR_sphere)=nir
      cR2irr(2,jc1,jc2,jR_sphere)=1
      do isym=1,sym%nsym
        SI=dble(sym%lat(:,:,isym))
        ic1=proj%ics2c(jc1,isym)
        ic2=proj%ics2c(jc2,isym)
        ! In principle, we should find RI, which gives the following RK and proceed with only one
        ! which gives RK=RJ (derivation from the paper)
        ! RK=matmul(SJ,RI+proj%centers(:,ic2)-proj%centers(:,ic1))-(proj%centers(:,jc2)-proj%centers(:,jc1))
        ! However, we can know in advance by using a direct symmetry operation to equation above
        RI=matmul(SI,RJ+proj%centers(:,jc2)-proj%centers(:,jc1))-(proj%centers(:,ic2)-proj%centers(:,ic1))
        iv=rgrid%find(RI)
        if (iv(NDIM+1).lt.0.or.sum(abs(iv(1:NDIM))).ne.0) cycle
        iR_sphere=rgrid%homo_to_sphere(iv(NDIM+1))
        if (iR_sphere.le.0) cycle
        ! Hubbard U from iR_sphere,ic1,ic2 is mapped to jR_sphere,jc1,jc2. 
        ! Therefore jR_sphere,jc1,jc2 must be marked as irreducible point
        ! first check if it alreadythere
        if (found(ic1,ic2,iR_sphere)) cycle
        found(ic1,ic2,iR_sphere)=.true.
        cR2irr(1,ic1,ic2,iR_sphere)=nir
        cR2irr(2,ic1,ic2,iR_sphere)=isym
      end do
    end do
  end do
end do
return
end subroutine
subroutine init_wws(sym,proj,wws)
class(CLsym), intent(in) :: sym
class(wbase), intent(in)  :: proj
complex(dp), intent(out)    :: wws(proj%norb,proj%norb,sym%nsym)
integer isym,iw,jw,ic,jc
wws=0._dp
do isym=1,sym%nsym
  do iw=1,proj%norb
     ic=proj%orb_icio(iw,1)
     jc=proj%ics2c(ic,isym)
     do jw=1,proj%norb
        if(proj%orb_icio(jw,1).ne.jc) cycle
        wws(jw,iw,isym)=proj%wws(.false.,sym%car(:,:,isym),iw,jw)
     end do
  end do
end do
return
end subroutine

subroutine compute_hubbardu_rs(pars,sym,tbmodel,kgrid,eval,evec)
class(CLpars), intent(in) :: pars
class(CLsym), intent(in) :: sym
class(CLtb), intent(in) :: tbmodel
class(GRID), intent(inout) :: kgrid
real(dp), intent(in) :: eval(pars%nstates,kgrid%npt)
complex(dp), intent(in) :: evec(tbmodel%norb_TB,pars%nstates,kgrid%npt)
integer ir,ik
real(dp), allocatable    :: vkl(:,:),vrl(:,:)
complex(dp), allocatable :: umat(:,:,:),udis(:,:,:)
complex(dp), allocatable :: wfmloc(:,:,:)
! variables for symmetry maps
integer, allocatable     :: irr2cR(:,:)
integer, allocatable     :: cR2irr(:,:,:,:)
integer                  :: nir
!
type(wbase)  proj
type(GRID) rgrid
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
call proj%init_smap(sym,pars)
call info("Wannier_interface%compute_hubbardu","constructing U parameters on the grid")
call info("Wannier_interface%compute_hubbardu","WARNING this code is working for pz basis orbitals only")
! copy rgrid object from TB
rgrid=tbmodel%rgrid
if (.not.rgrid%sphere_allocated) call rgrid%init_sphere(pars)
allocate(udis(pars%nstates,proj%norb,kgrid%npt))
allocate(umat(proj%norb,proj%norb,kgrid%npt))
allocate(wfmloc(tbmodel%norb_TB,proj%norb,rgrid%npt))
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
   tbmodel%norb_TB,pars%dis_win,udis,umat,vkl,vrl,eval,evec,wfmloc)
if (mp_mpi) write(*,*) "available OMP_NUM_THREADS, if used: ",proj%norb
if (mp_mpi) write(*,*) "available number of mpi processes, if used: ",tbmodel%norb_TB
if (.not.pars%HubU_diagonal) then
  allocate(irr2cR(3,proj%ncenters*proj%ncenters*rgrid%npt_sphere))
  allocate(cR2irr(2,proj%ncenters,proj%ncenters,rgrid%npt_sphere))
  call find_irrcR_sphere(proj,sym,rgrid,irr2cR,cR2irr,nir)
  call UH_irr(pars,rgrid,tbmodel,proj,wfmloc,irr2cR,nir)
else
  call UH_full(pars,rgrid,tbmodel,proj,wfmloc)
end if
deallocate(vkl,vrl,umat,udis)
deallocate(wfmloc)
if (.not.pars%HubU_diagonal) then
  deallocate(irr2cR)
  deallocate(cR2irr)
end if
#ifdef MPI
  call mpi_barrier(mpi_com,mpi_err)
#endif
return
end subroutine

subroutine UH_irr(pars,rgrid,tbmodel,proj,wfmloc,irr2cR,nir)
class(CLpars), intent(in) :: pars
class(GRID), intent(in) :: rgrid
class(CLtb), intent(in) :: tbmodel
class(wbase), intent(in) :: proj
complex(dp), intent(in) :: wfmloc(tbmodel%norb_TB,proj%norb,rgrid%npt)
integer, intent(in)     :: irr2cR(3,proj%ncenters*proj%ncenters*rgrid%npt_sphere)
integer, intent(in)     :: nir
real(dp), parameter :: Upz=10._dp
real(dp), parameter :: eps=1.e-17_dp
logical exst
character(len=4) snum
integer IRR_POINT,jc1_irr,jc2_irr,jR_irr
integer iRp,jRp,iorb,jorb
integer n1,n2,n3,n4
integer nc1,nc2,nc3,nc4
real(dp) dij,vcl
real(dp) rij(NDIM)
complex(dp) z1
real(dp), allocatable    :: vpcorb(:,:),vpc(:,:)
complex(dp), allocatable :: wf(:,:,:)
complex(dp), allocatable :: UH(:,:,:,:)
type(coulrs) vcoul
call vcoul%init(pars%coulrs_file)
allocate(vpc(NDIM,rgrid%npt))
do n1=1,rgrid%npt
  vpc(:,n1)=rgrid%vpc(n1)
end do
allocate(vpcorb(NDIM,tbmodel%norb_TB))
do iorb=1,tbmodel%norb_TB
  vpcorb(:,iorb)=matmul(tbmodel%vplorb(iorb),pars%avec)
end do
allocate(wf(tbmodel%norb_TB,proj%norb,rgrid%npt))
allocate(UH(proj%norb,proj%norb,proj%norb,proj%norb))
if(mp_mpi) call system("mkdir -p _UH")
do IRR_POINT=1,nir
  write(snum,'(I4.4)') IRR_POINT
  inquire(file='_UH/UH'//snum,exist=exst)
  if (exst) then
     if (mp_mpi) write(*,*) "SKIPPING IRR_POINT: ",IRR_POINT," of ",nir
     cycle
  else
     if (mp_mpi) write(*,*) "COMPUTING IRR_POINT: ",IRR_POINT," of ",nir
  end if
  jc1_irr=irr2cR(1,IRR_POINT)
  jc2_irr=irr2cR(2,IRR_POINT)
  jR_irr=irr2cR(3,IRR_POINT)
  ! shift wave function
  do n1=1,proj%norb
    call wannerfunc_at_R(rgrid,tbmodel%norb_TB,rgrid%vpl_sphere(jR_irr),wfmloc(:,n1,:),wf(:,n1,:))
  end do
  UH=0._dp
#ifdef MPI
  call mpi_barrier(mpi_com,mpi_err)
#endif
  !$OMP PARALLEL DEFAULT (SHARED)&
  !$OMP PRIVATE(iRp,jRp,iorb,jorb,rij,dij)&
  !$OMP PRIVATE(n1,n2,n3,n4,vcl,z1)&
  !$OMP PRIVATE(nc1,nc2,nc3,nc4)
  !$OMP DO
  do n4=1,proj%norb
    do n3=1,proj%norb
      do n2=1,proj%norb
        do n1=1,proj%norb
          nc4=proj%orb_icio(n4,1)
          nc3=proj%orb_icio(n3,1)
          nc2=proj%orb_icio(n2,1)
          nc1=proj%orb_icio(n1,1)
          if (jc1_irr.ne.nc1) cycle
          if (jc2_irr.ne.nc2) cycle
          if (nc1.ne.nc4) cycle ! this condition enters symmetry transformation formulas of U 
          if (nc2.ne.nc3) cycle ! so we also can impose it here
          if (pars%HubU_diagonal.and.n1.ne.n4) cycle 
          if (pars%HubU_diagonal.and.n2.ne.n3) cycle 
          do iRp=1,rgrid%npt
            if (rgrid%dc(iRp).gt.2._dp*pars%rcut_grid) cycle
            do jRp=1,rgrid%npt
              if (rgrid%dc(jRp).gt.2._dp*pars%rcut_grid) cycle
              rij=abs(vpc(:,jRp)-vpc(:,iRp))
              dij=sqrt(dot_product(rij,rij))
              if (dij.gt.2._dp*pars%rcut_grid) cycle
              do iorb=1,tbmodel%norb_TB
                if (mod(iorb-1,np_mpi).ne.lp_mpi) cycle
                do jorb=1,tbmodel%norb_TB
                  rij=abs(vpcorb(:,jorb)+vpc(:,jRp)-vpcorb(:,iorb)-vpc(:,iRp))
                  dij=sqrt(dot_product(rij,rij))
                  if (dij.gt.pars%rcut_grid) cycle
                  vcl=vcoul%evaluate(dij)
                  if (abs(conjg(wfmloc(iorb,n1,iRp))).lt.epsengy) cycle
                  if (abs(conjg(wf(jorb,n2,jRp))).lt.epsengy) cycle
                  if (abs(wf(jorb,n3,jRp)).lt.epsengy) cycle
                  if (abs(wfmloc(iorb,n4,iRp)).lt.epsengy) cycle
                  z1=conjg(wfmloc(iorb,n1,iRp))*conjg(wf(jorb,n2,jRp))*wf(jorb,n3,jRp)*wfmloc(iorb,n4,iRp)
                  if (dij.lt.epslat) then
                     UH(n1,n2,n3,n4)=UH(n1,n2,n3,n4)+z1*Upz
                  else
                     UH(n1,n2,n3,n4)=UH(n1,n2,n3,n4)+z1*vcl
                  end if
                end do
              end do
            end do 
          end do
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
#ifdef MPI
  n1=proj%norb**4
  call mpi_allreduce(mpi_in_place,UH,n1,mpi_double_complex,mpi_sum, &
   mpi_com,mpi_err)
#endif
  if (mp_mpi) open(140,file='_UH/UH'//snum,action="write")
  if (mp_mpi) write(140,*) rgrid%npt_sphere,proj%norb
  do n1=1,proj%norb ; do n2=1,proj%norb ; do n3=1,proj%norb ; do n4=1,proj%norb
    if (abs(UH(n1,n2,n3,n4)).gt.epsengy) then
      if (mp_mpi) write(140,'(5I6,2G18.10)') jR_irr,n1,n2,n3,n4,UH(n1,n2,n3,n4)
    end if
  end do ; end do ; end do ; end do
  if (mp_mpi) close(140)
end do
deallocate(wf,UH)
deallocate(vpc,vpcorb)
return
end subroutine

subroutine UH_full(pars,rgrid,tbmodel,proj,wfmloc)
class(CLpars), intent(in) :: pars
class(GRID), intent(in) :: rgrid
class(CLtb), intent(in) :: tbmodel
class(wbase), intent(in) :: proj
complex(dp), intent(in) :: wfmloc(tbmodel%norb_TB,proj%norb,rgrid%npt)
real(dp), parameter :: Upz=10._dp
real(dp), parameter :: epscoul=10._dp
real(dp), parameter :: eps=1.e-17_dp
logical exst
character(len=4) snum
integer jR_sphere
integer iRp,jRp,iorb,jorb
integer n1,n2,n3,n4
real(dp) dij,vcl
real(dp) rij(NDIM)
complex(dp) z1
real(dp), allocatable    :: vpcorb(:,:),vpc(:,:)
complex(dp), allocatable :: wf(:,:,:)
complex(dp), allocatable :: UH(:,:,:,:)
type(coulrs) vcoul
call vcoul%init(pars%coulrs_file)
allocate(vpc(NDIM,rgrid%npt))
do n1=1,rgrid%npt
  vpc(:,n1)=rgrid%vpc(n1)
end do
allocate(vpcorb(NDIM,tbmodel%norb_TB))
! coordinates of basis orbitals
do iorb=1,tbmodel%norb_TB
  vpcorb(:,iorb)=matmul(tbmodel%vplorb(iorb),pars%avec)
end do
allocate(wf(tbmodel%norb_TB,proj%norb,rgrid%npt))
allocate(UH(proj%norb,proj%norb,proj%norb,proj%norb))
if(mp_mpi) call system("mkdir -p _UH")
do jR_sphere=1,rgrid%npt_sphere
  write(snum,'(I4.4)') jR_sphere
  inquire(file='_UH/UH'//snum,exist=exst)
  if (exst) then
     if (mp_mpi) write(*,*) "SKIPPING JR: ",jR_sphere," of ",rgrid%npt_sphere
     cycle
  else
     if (mp_mpi) write(*,*) "COMPUTING JR: ",jR_sphere," of ",rgrid%npt_sphere
  end if
  ! shift wave function to jR_sphere unit cell
  do n1=1,proj%norb
    call wannerfunc_at_R(rgrid,tbmodel%norb_TB,rgrid%vpl_sphere(jR_sphere),wfmloc(:,n1,:),wf(:,n1,:))
  end do
  UH=0._dp
#ifdef MPI
  call mpi_barrier(mpi_com,mpi_err)
#endif
  !$OMP PARALLEL DEFAULT (SHARED)&
  !$OMP PRIVATE(iRp,jRp,iorb,jorb,rij,dij)&
  !$OMP PRIVATE(n1,n2,n3,n4,vcl,z1)
  !$OMP DO
  do n4=1,proj%norb
    do n3=1,proj%norb
      do n2=1,proj%norb
        do n1=1,proj%norb
          if (pars%HubU_diagonal.and.n1.ne.n4) cycle 
          if (pars%HubU_diagonal.and.n2.ne.n3) cycle 
          ! we need to cycle all n1,n4 which are not on one center
          if (proj%orb_icio(n1,1).ne.proj%orb_icio(n4,1)) cycle
          ! the same for n2 n3; it is enough to get symmetry of Hubbard U
          if (proj%orb_icio(n2,1).ne.proj%orb_icio(n3,1)) cycle
          do iRp=1,rgrid%npt
            if (rgrid%dc(iRp).gt.2._dp*pars%rcut_grid) cycle
            do jRp=1,rgrid%npt
              if (rgrid%dc(jRp).gt.2._dp*pars%rcut_grid) cycle
              rij=abs(rgrid%vpc(jRp)-rgrid%vpc(iRp))
              dij=sqrt(dot_product(rij,rij))
              if (dij.gt.2._dp*pars%rcut_grid) cycle
              do iorb=1,tbmodel%norb_TB
                if (mod(iorb-1,np_mpi).ne.lp_mpi) cycle
                do jorb=1,tbmodel%norb_TB
                  rij=abs(vpcorb(:,jorb)+vpc(:,jRp)-vpcorb(:,iorb)-vpc(:,iRp))
                  dij=sqrt(dot_product(rij,rij))
                  if (dij.gt.pars%rcut_grid) cycle
                  vcl=vcoul%evaluate(dij)
                  if (abs(conjg(wfmloc(iorb,n1,iRp))).lt.epsengy) cycle
                  if (abs(conjg(wf(jorb,n2,jRp))).lt.epsengy) cycle
                  if (abs(wf(jorb,n3,jRp)).lt.epsengy) cycle
                  if (abs(wfmloc(iorb,n4,iRp)).lt.epsengy) cycle
                  z1=conjg(wfmloc(iorb,n1,iRp))*conjg(wf(jorb,n2,jRp))*wf(jorb,n3,jRp)*wfmloc(iorb,n4,iRp)
                  if (dij.lt.epslat) then
                     !HubU(nn,mm,pp,qq,jR_sphere)=HubU(nn,mm,pp,qq,jR_sphere)+z1*Upz/CoulombForceConstant
                     UH(n1,n2,n3,n4)=UH(n1,n2,n3,n4)+z1*Upz
                  else
                     !HubU(nn,mm,pp,qq,jR_sphere)=HubU(nn,mm,pp,qq,jR_sphere)+z1/dij
                     UH(n1,n2,n3,n4)=UH(n1,n2,n3,n4)+z1*vcl
                  end if
                end do
              end do
            end do 
          end do
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
#ifdef MPI
  n1=proj%norb**4
  call mpi_allreduce(mpi_in_place,UH,n1,mpi_double_complex,mpi_sum, &
   mpi_com,mpi_err)
#endif
!HubU=HubU*CoulombForceConstant/epscoul
  if (mp_mpi) open(140,file='_UH/UH'//snum,action="write")
  if (mp_mpi) write(140,*) rgrid%npt_sphere,proj%norb
  do n1=1,proj%norb ; do n2=1,proj%norb ; do n3=1,proj%norb ; do n4=1,proj%norb
    if (abs(UH(n1,n2,n3,n4)).gt.epsengy) then
      if (mp_mpi) write(140,'(5I6,2G18.10)') jR_sphere,n1,n2,n3,n4,UH(n1,n2,n3,n4)
    end if
  end do ; end do ; end do ; end do
  if (mp_mpi) close(140)
end do
deallocate(wf,UH)
deallocate(vpc,vpcorb)
return
end subroutine
subroutine compute_hubbardj_rs(pars,tbmodel,kgrid,eval,evec)
class(CLpars), intent(in) :: pars
class(CLtb), intent(in) :: tbmodel
class(GRID), intent(inout) :: kgrid
real(dp), intent(in) :: eval(pars%nstates,kgrid%npt)
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
   tbmodel%norb_TB,pars%dis_win,udis,umat,vkl,vrl,eval,evec,wfmloc)
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
  write(*,*) "JR: ",jR_sphere
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
subroutine get_wfmloc(get_smooth,nwan,num_bands,nkpt,nrpt,nbasis,&
        dis_win,udis,umat,vkl,vrl,eval,evec,wfmloc)
logical, intent(in) :: get_smooth
integer, intent(in) :: nwan,num_bands,nkpt,nrpt,nbasis
real(dp), intent(in) :: vkl(NDIM,nkpt),vrl(NDIM,nkpt)
real(dp), intent(in) :: dis_win(2)
complex(dp), intent(in) :: udis(num_bands,nwan,nkpt),umat(nwan,nwan,nkpt)
real(dp), intent(in) :: eval(num_bands,nkpt)
complex(dp), intent(in) :: evec(nbasis,num_bands,nkpt)
complex(dp), intent(out) :: wfmloc(nbasis,nwan,nrpt)
integer ik,ir,n,m,num_win,idx(num_bands)
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
    num_win=0
    do m=1,num_bands
      if (eval(m,ik).gt.dis_win(1).and.eval(m,ik).lt.dis_win(2)) then
        num_win=num_win+1
        idx(num_win)=m
      end if
    end do
    psik(:,:)=0._dp
    do n=1,nwan
      do m=1,num_win
        psik(:,n)=psik(:,n)+udis(m,n,ik)*evec(:,idx(m),ik)
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
  write(50,'("f=",G18.10,",",G18.10,",",G18.10,":",A10)') &
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
        write(50,'(5G18.10)') tbmodel%wbase%centers(:,ic),wf(iorb,iw,jr)
      end do
    end if
  end do
  close(50)
end do
return
end subroutine

subroutine read_wf_universe(tbmodel,pars,nr,wf,pre,post)
class(CLtb), intent(in) :: tbmodel
class(CLpars), intent(in) :: pars
integer, intent(in) :: nr
complex(dp), intent(out) :: wf(tbmodel%norb_TB,pars%proj%norb,nr)
character(len=*),intent(in) :: pre
character(len=*),intent(in) :: post
character(len=128) fname
character(len=128) num
character(len=128) stmp
real(dp) avec(NDIM,NDIM)
integer iorb,iR,iw,norb_TB,nR_
integer irg(NDIM+1)
real(dp) x,y,z,vl(NDIM)
real(dp) c1,c2
wf=0._dp
do iw=1,pars%proj%norb
  write(num,'(I2)') iw
  fname=trim(adjustl(pre))//trim(adjustl(num))//trim(adjustl(post))//'.dat'
  if (mp_mpi) write(*,*) 'reading file:', trim(fname)
  open(50,file=trim(adjustl(fname)))
  read(50,*) avec(1,:)
  read(50,*) avec(2,:)
  read(50,*) avec(3,:)
  read(50,*) stmp,norb_TB
  read(50,*) stmp,nR_
  if (norb_TB.ne.tbmodel%norb_TB) call throw('read_wf_universe','norb_TB is not what is expected')

  do iR=1,nr
    if (sum(abs(tbmodel%rgrid%vpl(iR))).le.12) then
      read(50,*,end=200) vl(:)
      irg=tbmodel%rgrid%find(vl)
      if (irg(NDIM+1).lt.0.or.sum(abs(irg(1:NDIM))).ne.0) call throw('read_wf_universe','R point not found')
      do iorb=1,tbmodel%norb_TB
        read(50,*) x,y,z,c1,c2
        wf(iorb,iw,irg(NDIM+1))=cmplx(c1,c2,kind=dp)
      end do
    end if
  end do
  200 continue
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
ngrid(3)=1
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
call info("read_tb_file","read_tb_file done")
if (pars%writetb) call THIS%write_tb_file('hamwan0',pars,pars%proj%norb)
deallocate(ivr,deg,ham)
return
end subroutine
subroutine write_tb_file(THIS,prefix,pars,norb)
class(CLwan), intent(inout) :: THIS
character(len=*), intent(in) :: prefix
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
open(50,file='_ham/'//trim(adjustl(prefix))//'_tb.dat',action='write')
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
