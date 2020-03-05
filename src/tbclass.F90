
module tbclass
#ifdef MPI
  use mpi
#endif
use modcom
use parameters
use gridclass
use slater_koster
use symmetryclass
use wannier_supplementary
implicit none
private
integer, parameter :: maxallowd_orbs=46340
integer, parameter :: maxallowd_nn=14400

type, public :: CLtb
  ! private
  character(len=100), private :: mode
  character(len=100), private :: sktype
  integer :: norb_TB
  integer, private    :: nspec
  integer, private    :: hamsize
  integer, public     :: ncenters
  integer, private    :: nmaxatm_pspec
  integer, private    :: nmaxorb_pspec
  real(dp), private   :: rcut_nn
  real(dp), private   :: sparse_eps=1.e-6_dp
  type(SK), private   :: skfunc
  type(GRID), public  :: rgrid
  ! sparse indexing (isa,jsa)
  integer, allocatable, private :: lmr(:,:)
  integer, allocatable, private :: isa(:)
  integer, allocatable, private :: jsa(:)
  integer, allocatable, public  :: norb_ic(:)
  integer, allocatable, private :: norb_nn(:)
  integer, allocatable, private :: nn_idx(:,:,:)
  integer, allocatable, private :: ic_ispec(:)
  integer, allocatable, public  :: orb_icio(:,:)
  integer, allocatable, public  :: icio_orb(:,:)
  integer, allocatable, private :: ncenters_nn(:)
  integer, allocatable, private :: jcjr_nn(:,:,:)
  real(dp), allocatable, private :: centers(:,:)
  real(dp), allocatable, private :: centers_cart(:,:)
  real(dp), allocatable, private :: deg_hr(:)
  real(dp), allocatable, private :: waxis(:,:,:)
  complex(dp), allocatable, private :: ham_hr(:,:)
  contains
  procedure, public :: init=>init_variables
  procedure, public :: evalk=>calc_eigenvalues_at_K
  procedure, public :: vplorb=>get_location_of_orbial
  procedure, public :: wfGtransform=>bloch_wf_transform_at_gamma
  procedure, private :: hK=>give_hK
  procedure, private :: hR=>give_hR
  procedure, private :: tij=>tij_function
  procedure, private :: findnn=>run_findnn
  procedure, private :: inquire_hamsize=>calc_hamsize

endtype CLtb


contains


subroutine init_variables(THIS,pars,mode)
class(CLtb), intent(out) :: THIS
class(CLpars), intent(inout) :: pars
character(len=*), intent(in) :: mode
integer iat,ispec,ios,iorb,ic
call info ("CLtb%init_variables","")
THIS%mode=trim(adjustl(mode))
THIS%sparse_eps=pars%sparse_eps
THIS%rcut_nn=pars%rcut_nn
THIS%nspec=pars%nspec
THIS%nmaxatm_pspec=pars%nmaxatm_pspec
THIS%nmaxorb_pspec=maxval(pars%norb_per_center(:))
! this can be changed later
THIS%ncenters=pars%natmtot
THIS%sktype=pars%sktype
allocate(THIS%lmr(2,THIS%ncenters))
allocate(THIS%norb_ic(THIS%ncenters))
allocate(THIS%waxis(NDIM,2,THIS%ncenters))
allocate(THIS%icio_orb(THIS%ncenters,THIS%nmaxorb_pspec))
! this corresponds to PZ orbitals:
THIS%lmr(1,:)=1 !l
THIS%lmr(2,:)=1 !mr
THIS%norb_ic=pars%norb_per_center
THIS%waxis=pars%wannier_axis
! compute first the total number of orbitals
THIS%norb_TB=0
do ic=1,THIS%ncenters
  do ios=1,THIS%norb_ic(ic)
     THIS%norb_TB=THIS%norb_TB+1
     if(THIS%norb_TB.gt.maxallowd_orbs) then
       call throw("CLtb%init_variables()","maximum allowed basis orbitals exceeded (integer overflow in NxN value)")
     end if
     ! map: center,orbialindex to the final basis orbital
     THIS%icio_orb(ic,ios)=THIS%norb_TB
  end do
end do
allocate(THIS%ic_ispec(THIS%ncenters))
allocate(THIS%orb_icio(THIS%norb_TB,3))
do ic=1,THIS%ncenters
  ispec=pars%tot_iais(ic,2)
  THIS%ic_ispec(ic)=ispec
  do ios=1,THIS%norb_ic(ic)
     iorb=THIS%icio_orb(ic,ios)
     ! map:  final basis orbital to center,orbialindex
     THIS%orb_icio(iorb,1)=ic
     THIS%orb_icio(iorb,2)=ios
     THIS%orb_icio(iorb,3)=ispec
  end do
end do
! compute the actual centers
allocate(THIS%centers(NDIM,THIS%ncenters))
allocate(THIS%centers_cart(NDIM,THIS%ncenters))
ic=0
do ispec=1,pars%nspec
  do iat=1,pars%nat_per_spec(ispec)
    ic=ic+1
    THIS%centers(:,ic)=pars%atml(:,iat,ispec)
    THIS%centers_cart(:,ic)=pars%atmc(iat,ispec)
  end do
end do
! init real space grid
call THIS%rgrid%init(pars%ngrid,pars%ngrid,pars%avec,.true.,.false.)
! Fix the parameters for bottom and top states to calculate
if (pars%istop.gt.THIS%norb_TB) then
  pars%istop=THIS%norb_TB
  call message("top state to compute was changed")
end if
if (pars%istart.gt.pars%istop) then
  pars%istart=pars%istop
  call message("bottom state to compute was changed")
end if
pars%nstates=pars%istop-pars%istart+1
if (trim(adjustl(THIS%mode)).eq.'tbfile'.or.&
  trim(adjustl(THIS%mode)).eq.'hrfile'.or.&
  trim(adjustl(THIS%mode)).eq.'datfile') then
  call throw("CLtb%init_variables()","reading TB hamiltonian from file is not implemented yet")
else if (trim(adjustl(THIS%mode)).eq.'SK') then
  call THIS%skfunc%init(THIS%sktype)
  call THIS%findnn()
  call THIS%inquire_hamsize()
else if (trim(adjustl(THIS%mode)).eq.'noham') then
   return
else
  call throw("CLtb%init_variables()","unknown TB initialization mode")
end if
if (mp_mpi) then
  write(*,*) "max number of nearest neighbors found: ", maxval(THIS%ncenters_nn)
  write(*,*) "number of non-zeros in the upper-triangular part of Hamiltonian: ", THIS%hamsize
  write(*,'(" upper-triangular sparsisity: ",F18.10)') dble(THIS%hamsize)/&
                            ( 0.5_dp*dble(THIS%norb_TB*THIS%norb_TB)+0.5_dp*dble(THIS%norb_TB) )
end if
if (mp_mpi) write(*,*)
end subroutine

subroutine calc_eigenvalues_at_K(THIS,lvec,pars,vpl,eval,evec)
class(CLtb), intent(in) :: THIS
logical, intent (in) :: lvec
class(CLpars), intent(in) :: pars
real(dp), intent(in) :: vpl(NDIM)
real(dp), intent(out) :: eval(pars%nstates)
complex(dp), intent(out) :: evec(THIS%norb_TB,pars%nstates)
! local
integer ii,nn
integer, allocatable :: sidx(:)
real(dp), allocatable :: rtmp(:)
complex(dp), allocatable :: ztmp(:,:)
complex(dp), allocatable :: ham(:)
allocate(ham(THIS%hamsize))
call THIS%hK(vpl,ham)
! Extract the Fermi level from diagonal,
! ir is needed to simplify Sparse solver interface (i.e. avoid using explicit shift)
do nn=1,THIS%hamsize
  do ii=1,THIS%norb_TB
    if (nn.ge.THIS%isa(ii).and.nn.lt.THIS%isa(ii+1)) then
       if (ii.eq.THIS%jsa(nn)) ham(nn)=ham(nn)-pars%efermi
    end if
  end do
end do
if (pars%sparse) then
  call arpack_interface(lvec,'I','LM',pars%nstates,THIS%norb_TB,&
                          THIS%hamsize,THIS%isa,THIS%jsa,ham,eval,evec)
  ! sort
  allocate(sidx(pars%nstates))
  allocate(rtmp(pars%nstates))
  if (lvec) allocate(ztmp(THIS%norb_TB,pars%nstates))
  call sortidx(pars%nstates,eval,sidx)
  rtmp=eval
  if (lvec) ztmp=evec
  do ii=1,pars%nstates
    ! eigenvalues
    eval(ii)=rtmp(sidx(ii))
    ! eigenvectors
    if (lvec) evec(:,ii)=ztmp(:,sidx(ii))
  end do
  deallocate(sidx,rtmp)
  if (lvec) deallocate(ztmp)
else
  allocate(ztmp(THIS%norb_TB,THIS%norb_TB))
  ! unpack sparse matrix to the full one in its upper triangular part
  ztmp=0._dp
  do nn=1,THIS%hamsize
    do ii=1,THIS%norb_TB
      if (nn.ge.THIS%isa(ii).and.nn.lt.THIS%isa(ii+1)) then
        ztmp(ii,THIS%jsa(nn))=ham(nn)
      end if
    end do
  end do
  if (THIS%norb_TB.eq.pars%nstates) then
    call eigenv_problem(THIS%norb_TB,ztmp,eval)
    evec=ztmp
  else
    if (lvec) then
     ! routine from wannier tools
      call zheevx_pack('V', 'U', THIS%norb_TB, pars%istart, pars%istop , ztmp, eval, evec)
    else
     ! routine from wannier tools
      call zheevx_pack('N', 'U', THIS%norb_TB, pars%istart, pars%istop , ztmp, eval, evec)
    end if
  end if
  deallocate(ztmp)
end if
!write(*,*) abs(dot_product(evec(:,1),evec(:,4))),abs(dot_product(evec(:,1),evec(:,1))),sum(abs(evec(:,1)))
! add back the Fermi level
eval=eval+pars%efermi
deallocate(ham)
end subroutine

subroutine give_hK(THIS,vpl,hamk)
class(CLtb), intent(in) :: THIS
real(dp), intent(in) :: vpl(NDIM)
complex(dp), intent(out) :: hamk(THIS%hamsize)
integer iR
real(dp) t1
complex(dp) z1
complex(dp), allocatable :: hamr(:)
allocate(hamr(THIS%hamsize))
hamk=0._dp
do iR=1,THIS%rgrid%npt
  t1=twopi*dot_product(vpl,THIS%rgrid%vpl(iR))
  z1=cmplx(cos(t1),sin(t1),kind=dp)
  call THIS%hR(iR,hamr)
  ! LAPACK's B=\alpha*X+B  (z1=alpha,X=hamr,B=hamk)
  call zaxpy(THIS%hamsize,z1,hamr,1,hamk,1)
end do
deallocate(hamr)
end subroutine

subroutine give_hR(THIS,itr,ham)
class(CLtb), intent(in) :: THIS
integer, intent(in) :: itr
complex(dp), intent(out) :: ham(THIS%hamsize)
integer iorb,jorb,ic,jc,jr,ispec,jspec,ios,jos,nn,idx,ii
ham=0._dp
if (THIS%mode.eq."hr") then
  do iorb=1,THIS%norb_TB
    do jorb=iorb,THIS%norb_TB
      ! sparse index of the first non-zero element in the row iorb
      ii=THIS%isa(iorb)
      do idx=ii,THIS%hamsize
        ! jsa is a map from sparse index to the column one
        if (THIS%jsa(idx).eq.jorb) then
          ham(idx)=THIS%ham_hr(idx,itr)/THIS%deg_hr(itr)
        end if
      end do
    end do
  end do
else
  do ic=1,THIS%ncenters
    ispec=THIS%ic_ispec(ic)
    do nn=1,THIS%ncenters_nn(ic)
      jc=THIS%jcjr_nn(ic,nn,1)
      jr=THIS%jcjr_nn(ic,nn,2)
      jspec=THIS%ic_ispec(jc)
      if (itr.ne.jr) cycle
      do ios=1,THIS%norb_ic(ic)
        iorb=THIS%icio_orb(ic,ios)
        ! sparse index of the first non-zero element in the row iorb
        ii=THIS%isa(iorb)
        do jos=1,THIS%norb_ic(jc)
          jorb=THIS%icio_orb(jc,jos)
          do idx=ii,min(THIS%isa(iorb+1),THIS%hamsize)
             ! jsa is a map from sparse index to the column one
             if(THIS%jsa(idx).eq.jorb) then
               ham(idx)=THIS%tij(ic,jc,ios,jos,jr)
               exit
             end if
          end do
        end do
      end do
    end do
  end do
end if
end subroutine

subroutine calc_hamsize(THIS)
class(CLtb), intent(inout) :: THIS
integer iorb,jorb,ic,jc,jr,ispec,jspec,ios,jos,nn
logical, allocatable :: nonzero_indices(:,:)
allocate(nonzero_indices(THIS%norb_TB,THIS%norb_TB))
nonzero_indices=.false.
do ic=1,THIS%ncenters
  ispec=THIS%ic_ispec(ic)
  do nn=1,THIS%ncenters_nn(ic)
    jc=THIS%jcjr_nn(ic,nn,1)
    jr=THIS%jcjr_nn(ic,nn,2)
    jspec=THIS%ic_ispec(jc)
    do ios=1,THIS%norb_ic(ic)
      iorb=THIS%icio_orb(ic,ios)
      do jos=1,THIS%norb_ic(jc)
        jorb=THIS%icio_orb(jc,jos)
        if (abs(THIS%tij(ic,jc,ios,jos,jr)).gt.THIS%sparse_eps.or.&
           iorb.eq.jorb) then
           nonzero_indices(iorb,jorb)=.true.
        end if
      end do
    end do
  end do
end do
THIS%hamsize=0
! compute total number of non-zeroes for upper triangular part
do iorb=1,THIS%norb_TB
  ! upper diagonal part
  do jorb=iorb,THIS%norb_TB
    if (nonzero_indices(iorb,jorb)) then
      THIS%hamsize=THIS%hamsize+1
    end if
  end do
end do
! allocate indices for sparse data format
allocate(THIS%isa(THIS%norb_TB+1))
allocate(THIS%jsa(THIS%hamsize))
! fill sparse data indices with apropriate values
nn=0
THIS%isa=0
THIS%jsa=0
THIS%isa(1)=1
do iorb=1,THIS%norb_TB
  do jorb=iorb,THIS%norb_TB
    if (nonzero_indices(iorb,jorb)) then
      nn=nn+1
      THIS%jsa(nn)=jorb
      THIS%isa(iorb+1)=THIS%isa(iorb+1)+1
     end if
  end do
end do
do iorb=2,THIS%norb_TB+1
   THIS%isa(iorb)=THIS%isa(iorb)+THIS%isa(iorb-1)
enddo
deallocate(nonzero_indices)
end subroutine

subroutine run_findnn(THIS)
class(CLtb), intent(inout) :: THIS
integer ic,jc,jr,ios,jos,ispec,jspec,iorb,jorb,nsize
real(dp) t1,dv(NDIM)
logical, allocatable :: centers_visited(:,:)
allocate(THIS%ncenters_nn(THIS%ncenters))
allocate(THIS%jcjr_nn(THIS%ncenters,maxallowd_nn,2))
THIS%ncenters_nn(:)=0
THIS%jcjr_nn(:,:,:)=0
#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
!$OMP PARALLEL DEFAULT(SHARED)&
!$OMP PRIVATE(jc,jR,dv,t1)&
!$OMP PRIVATE(ispec,jspec,ios,jos,iorb,jorb,centers_visited)
!$OMP DO
do ic=1,THIS%ncenters
  if (mod(ic-1,np_mpi).ne.lp_mpi) cycle
  allocate(centers_visited(THIS%ncenters,THIS%rgrid%npt))
  centers_visited=.false.
  ispec=THIS%ic_ispec(ic)
  do jc=1,THIS%ncenters
    jspec=THIS%ic_ispec(jc)
    do jR=1,THIS%rgrid%npt
      dv=THIS%centers_cart(:,jc)+THIS%rgrid%vpc(jR)-THIS%centers_cart(:,ic)
      t1=sqrt(dot_product(dv,dv))
      if (t1.lt.THIS%rcut_nn) then
        do ios=1,THIS%norb_ic(ic)
          iorb=THIS%icio_orb(ic,ios)
          do jos=1,THIS%norb_ic(jc)
            jorb=THIS%icio_orb(jc,jos)
            if (abs(THIS%tij(ic,jc,ios,jos,jR)).gt.THIS%sparse_eps.or.&
                iorb.eq.jorb) then
              if (.not.centers_visited(jc,jR)) then
                centers_visited(jc,jR)=.true.
                THIS%ncenters_nn(ic)=THIS%ncenters_nn(ic)+1
                if (THIS%ncenters_nn(ic).gt.maxallowd_nn) then
                  call throw("CLtb%findnn()","max allowed nearest neighbor exceeded, descrease rcut_nn")
                end if
                THIS%jcjr_nn(ic,THIS%ncenters_nn(ic),1)=jc
                THIS%jcjr_nn(ic,THIS%ncenters_nn(ic),2)=jR
              end if
            end if
          end do
        end do
      end if
    end do
  end do
  deallocate(centers_visited)
end do
!$OMP END DO
!$OMP END PARALLEL
#ifdef MPI
  nsize=THIS%ncenters
  call mpi_allreduce(mpi_in_place,THIS%ncenters_nn,nsize,mpi_integer,mpi_sum, &
   mpi_com,mpi_err)
  nsize=THIS%ncenters*maxallowd_nn*2
  call mpi_allreduce(mpi_in_place,THIS%jcjr_nn,nsize,mpi_integer,mpi_sum, &
   mpi_com,mpi_err)
#endif

end subroutine

complex(dp) function tij_function(THIS,ic,jc,ios,jos,jr)
class(CLtb), intent(in) :: THIS
integer, intent(in) :: ic,jc,ios,jos,jr
real(dp) dvec(NDIM)
dvec(:)=THIS%centers_cart(:,jc)+THIS%rgrid%vpc(jr)-THIS%centers_cart(:,ic)
if (ios.gt.1.or.jos.gt.1) call throw("CLtb%tij_function()","this function is currently for pz-pz hoppings only")
tij_function=THIS%skfunc%tij(THIS%sktype,0,0,dvec)
end function



function get_location_of_orbial(THIS,iorb) result(vpl)
class(CLtb), intent(in) :: THIS
integer, intent(in) :: iorb
integer ic
real(dp) vpl(NDIM)
ic=THIS%orb_icio(iorb,1)
vpl=THIS%centers(:,ic)
end function

subroutine bloch_wf_transform_at_gamma(THIS,pars,sym,isym,wfinout)
class(CLtb), intent(in) :: THIS
class(CLpars), intent(in) :: pars
class(CLsym), intent(in) :: sym
integer, intent(in) :: isym
complex(dp), intent(inout) :: wfinout(THIS%norb_TB)
complex(dp), allocatable :: wftemp(:)
! local
integer iat,jat,jspec
integer ic,jc,iorb,jorb,ios,jos
real(dp) t1
allocate(wftemp(THIS%norb_TB))
wftemp(:)=wfinout(:)
wfinout(:)=0._dp
do jc=1,THIS%ncenters
  ! iatom/specie index of jc cite
  jat=pars%tot_iais(jc,1)
  jspec=pars%tot_iais(jc,2)
  ! equivalent atom to jat,jspec at the jsym operation
  iat=sym%ieqat(jat,jspec,isym)
  ! cite index of iat,ispec=jspec
  ic=pars%iais_tot(iat,jspec)
  do ios=1,THIS%norb_ic(ic)
    iorb=THIS%icio_orb(ic,ios)
    do jos=1,THIS%norb_ic(jc)
      t1=give_symrep_me(sym%car(:,:,isym),THIS%waxis(:,:,ios),THIS%waxis(:,:,jos),THIS%lmr(:,ios),THIS%lmr(:,jos))
      jorb=THIS%icio_orb(jc,jos)
      wfinout(iorb)=wfinout(iorb)+t1*wftemp(jorb)
    end do
  end do
end do
deallocate(wftemp)
return
end subroutine


! overlap of the basis wannier functions rotatated on a cite
real(dp) function give_symrep_me(symcar,axis1,axis2,lmr1,lmr2)
real(dp), intent(in) :: symcar(NDIM,NDIM)
integer, intent(in) :: lmr1(2),lmr2(2)
real(dp), intent(in) :: axis1(NDIM,2),axis2(NDIM,2)
type(interface_wovlp) :: wovlp
give_symrep_me=wovlp%wws(symcar,lmr1(1),lmr1(2),lmr2(1),lmr2(2),&
                         axis1(:,1),axis1(:,2),axis2(:,1),axis2(:,2))
end function

end module
