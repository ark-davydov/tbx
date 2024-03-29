
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
!integer, parameter :: maxallowd_orbs=46340
integer, parameter :: maxallowd_nn=20000

type, public :: CLtb
  ! private
  character(len=100), private :: mode
  character(len=100), private :: sktype,sk_subtype(2)
  integer :: norb_TB
  integer, private    :: nspec
  integer, private    :: hamsize
  integer, private    :: nmaxatm_pspec
  real(dp), private   :: rcut_nn
  real(dp), private   :: rcut_tbg_nni
  real(dp), private   :: vecs(NDIM,NDIM)
  real(dp), private   :: bvec(NDIM,NDIM)
  real(dp), private   :: sparse_eps=1.e-6_dp
  real(dp), private   :: tbg_vhscale=0._dp
  real(dp), private   :: tbg_dmin=1.e6_dp
  real(dp), private   :: tbg_dmax=0._dp
  type(SK), private   :: skfunc
  type(GRID), public  :: rgrid
  ! sparse indexing (isa,jsa)
  integer, allocatable, private :: isa(:)
  integer, allocatable, private :: jsa(:)
  integer, allocatable, private :: norb_nn(:)
  integer, allocatable, private :: nn_idx(:,:,:)
  integer, allocatable, private :: ic_ispec(:)
  integer, allocatable, public  :: orb_ispec(:)
  integer, allocatable, private :: ncenters_nn(:)
  integer, allocatable, private :: jcjr_nn(:,:,:)
  real(dp), allocatable, private :: dege(:)
  real(dp), allocatable, private :: tbg_oo_dist(:)
  complex(dp), allocatable, private :: hame(:,:)
  complex(dp), allocatable, private :: hame_file(:,:,:)
  type(wbase) :: wbase
  contains
  procedure, public :: init=>init_variables
  procedure, public :: evalk=>calc_eigenvalues_at_K
  procedure, public :: vplorb=>get_location_of_orbial
  procedure, public :: bloch_wf_transform
  procedure, public :: bloch_wf_transform_statk
  procedure, private :: hK=>give_hK
  procedure, private :: hR=>give_hR
  procedure, private :: tij=>tij_function
  procedure, private :: findnn
  procedure, private :: inquire_hamsize=>calc_hamsize
  procedure, private :: write_nonzeros_hame
  procedure, private :: read_tb_file
  procedure, private :: write_tb_file
  procedure, private :: read_hr_file
  procedure, private :: write_hr_file
  procedure, private :: unfold_hame
  
endtype CLtb


contains


subroutine init_variables(THIS,pars,sym,mode)
class(CLtb), intent(out) :: THIS
class(CLpars), intent(inout) :: pars
class(CLsym), intent(inout) :: sym
character(len=*), intent(in) :: mode
integer ispec,ios,iorb,ic,jc,ii
real(dp) dvec(NDIM),t1,t2
call info ("CLtb%init_variables","")
if (trim(adjustl(mode)).eq.'noham') then
  THIS%mode='noham'
else
  THIS%mode=trim(adjustl(pars%tbtype))
end if
THIS%sparse_eps=pars%sparse_eps
THIS%rcut_nn=pars%rcut_nn
THIS%rcut_tbg_nni=pars%rcut_tbg_nni
THIS%vecs=pars%avec
THIS%bvec=pars%bvec
THIS%nspec=pars%nspec
THIS%nmaxatm_pspec=pars%nmaxatm_pspec
THIS%sktype=pars%sktype
THIS%sk_subtype(1)=pars%sk_subtype(1)
THIS%sk_subtype(2)=pars%sk_subtype(2)
THIS%tbg_vhscale=pars%tbg_vhscale
! allocate basis indices
call THIS%wbase%init(pars,pars%base%ncenters,pars%base%norb,pars%base%norb_ic,&
                 pars%base%lmr,pars%base%waxis,pars%base%centers)
call THIS%wbase%init_smap(sym,pars)
THIS%norb_TB=THIS%wbase%norb
allocate(THIS%ic_ispec(THIS%norb_TB))
allocate(THIS%orb_ispec(THIS%norb_TB))
do ic=1,THIS%wbase%ncenters
  ispec=pars%tot_iais(ic,2)
  THIS%ic_ispec(ic)=ispec
  do ios=1,THIS%wbase%norb_ic(ic)
     iorb=THIS%wbase%icio_orb(ic,ios)
     THIS%orb_ispec(iorb)=ispec
  end do
end do
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

if (trim(adjustl(THIS%mode)).eq.'tbfile') then
  call THIS%read_tb_file(pars)
  call THIS%findnn()
  call THIS%inquire_hamsize()
  call THIS%write_nonzeros_hame()
  if (pars%writetb.and.mp_mpi) call THIS%write_tb_file(pars)
else if (trim(adjustl(THIS%mode)).eq.'hrfile') then
  call THIS%read_hr_file(pars)
  call THIS%findnn()
  call THIS%inquire_hamsize()
  call THIS%write_nonzeros_hame()
  if (pars%writetb.and.mp_mpi) call THIS%write_hr_file()
else if (trim(adjustl(THIS%mode)).eq.'sk') then
  if (abs(THIS%tbg_vhscale)>epsengy) then
    allocate(THIS%tbg_oo_dist(THIS%wbase%ncenters))
    THIS%tbg_oo_dist=0._dp
     ! first we need to loop over all atoms to check what is the minimum and maximum interlayer distance
    do ic=1,THIS%wbase%ncenters
      t1=1.e6_dp
      do jc=1,THIS%wbase%ncenters
         dvec(:)=THIS%wbase%centers_cart(:,jc)-THIS%wbase%centers_cart(:,ic)
         if (abs(dvec(ZAXIS)).gt.0.5_dp*tbg_ab_distance) then
            THIS%tbg_dmin=min(THIS%tbg_dmin,abs(dvec(ZAXIS)))
            THIS%tbg_dmax=max(THIS%tbg_dmax,abs(dvec(ZAXIS)))
            ! find out-of-plane distance at ic, as the z-coordinate with smallest in-plane hopping offset
            t2=sqrt(dvec(1)**2+dvec(2)**2)
            if (t2<t1) then 
               THIS%tbg_oo_dist(ic)=abs(dvec(ZAXIS))
   !            write(*,*) ic,t1,t2,abs(dvec(ZAXIS)),THIS%tbg_oo_dist(ic)
               t1=t2
            end if
         end if
      end do
    end do
  end if
  ! init real space grid
  call THIS%rgrid%init(pars%ngrid,pars%avec,.true.,.false.)
  call THIS%skfunc%init(THIS%sktype,THIS%sk_subtype(1),THIS%sk_subtype(2))
  call THIS%findnn()
  call THIS%inquire_hamsize()
  if (pars%writetb.and.mp_mpi) call THIS%write_tb_file(pars)
else if (trim(adjustl(THIS%mode)).eq.'noham') then
   call THIS%rgrid%init(pars%ngrid,pars%avec,.true.,.false.)
   return
else
  write(*,*) trim(adjustl(THIS%mode))
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
integer iorb,jorb,ic,jc,jr,ios,jos,nn,idx,ii
ham=0._dp
if (trim(adjustl(THIS%mode)).eq.'tbfile'.or.trim(adjustl(THIS%mode)).eq.'hrfile') then
  ham(:)=THIS%hame(:,itr)/THIS%dege(itr)
else
  do ic=1,THIS%wbase%ncenters
    do nn=1,THIS%ncenters_nn(ic)
      jc=THIS%jcjr_nn(ic,nn,1)
      jr=THIS%jcjr_nn(ic,nn,2)
      if (itr.ne.jr) cycle
      do ios=1,THIS%wbase%norb_ic(ic)
        iorb=THIS%wbase%icio_orb(ic,ios)
        ! sparse index of the first non-zero element in the row iorb
        ii=THIS%isa(iorb)
        do jos=1,THIS%wbase%norb_ic(jc)
          jorb=THIS%wbase%icio_orb(jc,jos)
          do idx=ii,min(THIS%isa(iorb+1),THIS%hamsize)
             ! jsa is a map from sparse index to the column one
             if(THIS%jsa(idx).eq.jorb) then
               ham(idx)=THIS%tij(iorb,jorb,jr)
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
do ic=1,THIS%wbase%ncenters
  ispec=THIS%ic_ispec(ic)
  do nn=1,THIS%ncenters_nn(ic)
    jc=THIS%jcjr_nn(ic,nn,1)
    jr=THIS%jcjr_nn(ic,nn,2)
    jspec=THIS%ic_ispec(jc)
    do ios=1,THIS%wbase%norb_ic(ic)
      iorb=THIS%wbase%icio_orb(ic,ios)
      do jos=1,THIS%wbase%norb_ic(jc)
        jorb=THIS%wbase%icio_orb(jc,jos)
        if (THIS%mode.eq."sk") then
           if (abs(THIS%tij(iorb,jorb,jr)).gt.THIS%sparse_eps.or.&
              iorb.eq.jorb) then
              nonzero_indices(iorb,jorb)=.true.
           end if
        else if (trim(adjustl(THIS%mode)).eq.'tbfile'.or.trim(adjustl(THIS%mode)).eq.'hrfile') then
           if (abs(THIS%hame_file(iorb,jorb,jr)).gt.THIS%sparse_eps.or.&
              iorb.eq.jorb) then
              nonzero_indices(iorb,jorb)=.true.
           end if
        else
           call throw("tbclass%calc_hamsize()","unexpected tbtype")
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

subroutine findnn(THIS)
class(CLtb), intent(inout) :: THIS
integer ic,jc,jr,ios,jos,ispec,jspec,iorb,jorb,nsize
real(dp) t1,dv(NDIM)
logical, allocatable :: centers_visited(:,:)
allocate(THIS%ncenters_nn(THIS%wbase%ncenters))
allocate(THIS%jcjr_nn(THIS%wbase%ncenters,maxallowd_nn,2))
THIS%ncenters_nn(:)=0
THIS%jcjr_nn(:,:,:)=0
#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
!$OMP PARALLEL DEFAULT(SHARED)&
!$OMP PRIVATE(jc,jR,dv,t1)&
!$OMP PRIVATE(ispec,jspec,ios,jos,iorb,jorb,centers_visited)
!$OMP DO
do ic=1,THIS%wbase%ncenters
  if (mod(ic-1,np_mpi).ne.lp_mpi) cycle
  allocate(centers_visited(THIS%wbase%ncenters,THIS%rgrid%npt))
  centers_visited=.false.
  ispec=THIS%ic_ispec(ic)
  do jc=1,THIS%wbase%ncenters
    jspec=THIS%ic_ispec(jc)
    do jR=1,THIS%rgrid%npt
      dv=THIS%wbase%centers_cart(:,jc)+THIS%rgrid%vpc(jR)-THIS%wbase%centers_cart(:,ic)
      t1=sqrt(dot_product(dv,dv))
      if (t1.lt.THIS%rcut_nn) then
        if (THIS%mode.eq.'sk'.and. (trim(adjustl(THIS%sktype)).eq.'tbgsk')) then
              t1=sqrt(dv(1)**2+dv(2)**2)
              ! cut in-plane hopping for TBG (later for something else)
              !if (abs(dv(ZAXIS)).lt.0.5_dp*tbg_ab_distance) then
              if (t1.gt.THIS%rcut_tbg_nni) cycle
              !end if
        end if
        do ios=1,THIS%wbase%norb_ic(ic)
          iorb=THIS%wbase%icio_orb(ic,ios)
          do jos=1,THIS%wbase%norb_ic(jc)
            jorb=THIS%wbase%icio_orb(jc,jos)
            if (THIS%mode.eq.'sk') then
               t1=abs(THIS%tij(iorb,jorb,jR))
            else if (trim(adjustl(THIS%mode)).eq.'tbfile'.or.trim(adjustl(THIS%mode)).eq.'hrfile') then
               t1=abs(THIS%hame_file(iorb,jorb,jR))
            else
               call throw("tbclass%findnn()","unexpected tbtype")
            end if
            if (t1.gt.THIS%sparse_eps.or.&
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
  nsize=THIS%wbase%ncenters
  call mpi_allreduce(mpi_in_place,THIS%ncenters_nn,nsize,mpi_integer,mpi_sum, &
   mpi_com,mpi_err)
  nsize=THIS%wbase%ncenters*maxallowd_nn*2
  call mpi_allreduce(mpi_in_place,THIS%jcjr_nn,nsize,mpi_integer,mpi_sum, &
   mpi_com,mpi_err)
#endif

end subroutine

complex(dp) function tij_function(THIS,iorb,jorb,jr)
class(CLtb), intent(in) :: THIS
integer, intent(in) :: iorb,jorb,jr
! local
integer ic,jc
real(dp) dvec(NDIM)
ic=THIS%wbase%orb_icio(iorb,1)
jc=THIS%wbase%orb_icio(jorb,1)
dvec(:)=THIS%wbase%centers_cart(:,jc)+THIS%rgrid%vpc(jr)-THIS%wbase%centers_cart(:,ic)
tij_function=THIS%skfunc%tij(THIS%sktype,THIS%wbase%lmr(:,iorb),THIS%wbase%lmr(:,jorb),dvec,&
         THIS%sk_subtype(1),THIS%sk_subtype(2))

!write(*,*) dble(tij_function),aimag(tij_function)
! add hartree only when working with slater-koster for TBG
if (THIS%mode.eq.'sk'.and. (trim(adjustl(THIS%sktype)).eq.'tbgsk') ) then
    if (abs(THIS%tbg_vhscale)>epsengy .and. jr.eq.THIS%rgrid%ip0 .and. iorb.eq.jorb) then
       ! add to diagonal at home unit cell
       tij_function=tij_function+&
                    tbg_vhartree_of_z(abs(THIS%tbg_oo_dist(ic)),THIS%tbg_dmin,THIS%tbg_dmax,THIS%tbg_vhscale)
       !write(*,*) iorb,jorb,jr,tbg_vhartree_of_z(abs(THIS%tbg_oo_dist(ic)),THIS%tbg_dmin,THIS%tbg_dmax,THIS%tbg_vhscale)
       !write(*,*) THIS%tbg_dmin,THIS%tbg_dmax,abs(THIS%tbg_oo_dist(ic))
    end if
end if
end function

function get_location_of_orbial(THIS,iorb) result(vpl)
class(CLtb), intent(in) :: THIS
integer, intent(in) :: iorb
integer ic
real(dp) vpl(NDIM)
ic=THIS%wbase%orb_icio(iorb,1)
vpl=THIS%wbase%centers(:,ic)
end function

subroutine bloch_wf_transform(THIS,kgrid,ik,sym,isym,wfout,wfin)
class(CLtb), intent(in) :: THIS
class(GRID), intent(in) :: kgrid
integer, intent(in) :: ik
class(CLsym), intent(in) :: sym
integer, intent(in) :: isym
complex(dp), intent(out) :: wfout(THIS%norb_TB)
complex(dp), intent(in) :: wfin(THIS%norb_TB)
! local
integer iw,jw
integer ic,jc
real(dp) t1,t2
real(dp) :: v1(NDIM),v2(NDIM)
complex(dp) :: zz
wfout=0._dp
do iw=1,THIS%wbase%norb
   ic=THIS%wbase%orb_icio(iw,1)
   jc=THIS%wbase%ics2c(ic,sym%inv(isym))
   v1=kgrid%vpc(kgrid%iks2k(ik,isym))-matmul(sym%car(:,:,isym),kgrid%vpc(ik))
   v2=matmul(v1,sym%car(:,:,isym))
   t1=dot_product(THIS%wbase%vcs2t(:,jc,isym),kgrid%vpc(ik))
   t2=dot_product(sym%vtc(:,isym),v2)
   zz=cmplx(cos(t1+t2),sin(t1+t2),kind=dp)
   do jw=1,THIS%wbase%norb
     ic=THIS%wbase%orb_icio(jw,1)
     jc=THIS%wbase%ics2c(ic,isym)
     if(THIS%wbase%orb_icio(iw,1).ne.jc) cycle
     wfout(iw)=wfout(iw)+THIS%wbase%wws(.false.,sym%car(:,:,isym),jw,iw)*wfin(jw)*zz
   end do
end do
return
end subroutine

subroutine bloch_wf_transform_statk(THIS,vpc,vpcs,sym,isym,wfout,wfin)
class(CLtb), intent(in) :: THIS
real(dp), intent(in) :: vpc(NDIM),vpcs(NDIM)
class(CLsym), intent(in) :: sym
integer, intent(in) :: isym
complex(dp), intent(out) :: wfout(THIS%norb_TB)
complex(dp), intent(in) :: wfin(THIS%norb_TB)
! local
integer iw,jw
integer ic,jc
real(dp) t1,t2
complex(dp) :: zz
real(dp) :: v1(NDIM),v2(NDIM)
wfout=0._dp
do iw=1,THIS%wbase%norb
   ic=THIS%wbase%orb_icio(iw,1)
   jc=THIS%wbase%ics2c(ic,sym%inv(isym))
   v1=vpcs-matmul(sym%car(:,:,isym),vpc)
   v2=matmul(v1,sym%car(:,:,isym))
   t1=dot_product(THIS%wbase%vcs2t(:,jc,isym),vpc)
   t2=dot_product(sym%vtc(:,isym),v2)
   zz=cmplx(cos(t1+t2),sin(t1+t2),kind=dp)
   do jw=1,THIS%wbase%norb
     ic=THIS%wbase%orb_icio(jw,1)
     jc=THIS%wbase%ics2c(ic,isym)
     if(THIS%wbase%orb_icio(iw,1).ne.jc) cycle
     wfout(iw)=wfout(iw)+THIS%wbase%wws(.false.,sym%car(:,:,isym),jw,iw)*wfin(jw)*zz
   end do
end do
return
end subroutine


subroutine read_tb_file(THIS,pars)
class(CLtb), intent(inout) :: THIS
class(CLpars), intent(in) :: pars
integer i,j,ii,jj
integer nrpt,nwan,ir
integer ngrid(3),ivp(4)
real(dp) aa,bb,vpl(3),avec(3,3)
logical exs
integer, allocatable :: ivr(:,:)
integer, allocatable :: deg(:)
complex(dp), allocatable :: ham(:,:,:)
if (NDIM.ne.3) call throw("CLtb%read_tb_file","this subroutine works only in 3D case")
inquire(file=trim(adjustl(pars%TBfile)),exist=exs)
if (.not.exs) call throw("CLtb%read_tb_file"," file "//trim(adjustl(pars%TBfile))//" missing")
open(50,file=trim(adjustl(pars%TBfile)),action='read')
read(50,*)
read(50,*) avec(1,:)
read(50,*) avec(2,:)
read(50,*) avec(3,:)
if (sum(abs(avec-pars%avec)).gt.epslat) then
  call throw("CLtb%read_tb_file",&
  "lattice vectors in _tb file are different from ones in input. it can be resolved by writing more digits")
end if
read(50,*) nwan
read(50,*) nrpt
if (nwan.ne.THIS%wbase%norb) then
  call throw("CLtb%read_tb_file","number of basis orbitals in _tb file is different from one derived from the input")
end if
allocate(ivr(3,nrpt))
allocate(deg(nrpt))
allocate(ham(nwan,nwan,nrpt))
read(50,'(15I5)') (deg(ir),ir=1,nrpt)
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
allocate(THIS%dege(THIS%rgrid%npt))
allocate(THIS%hame_file(nwan,nwan,THIS%rgrid%npt))
THIS%dege=1._dp
THIS%hame_file=0._dp
do ir=1,nrpt
  vpl=dble(ivr(:,ir))
  ivp=THIS%rgrid%find(vpl)
  if (ivp(ndim+1).le.0.or.sum(abs(ivp(1:NDIM))).ne.0) cycle
  THIS%dege(ivp(4))=dble(deg(ir))
  THIS%hame_file(:,:,ivp(4))=ham(:,:,ir)
end do
deallocate(ivr,deg,ham)
return
end subroutine

subroutine write_tb_file(THIS,pars)
class(CLtb), intent(inout) :: THIS
class(CLpars), intent(in) :: pars
integer ii,jj,ir,jR,ic,jc,counter,ivp(NDIM+1)
real(dp) t1,dv(NDIM)
integer, allocatable :: deg(:)
complex(dp), allocatable :: ham(:,:,:)
type(GRID) rgrid
if (NDIM.ne.3) call throw("CLtb%write_tb_file","this subroutine works only in 3D case")
call rgrid%init(pars%ngrid,pars%avec,.true.,.false.)
allocate(deg(rgrid%npt))
! collect degeneracies, count active R-points
counter=0
do ir=1,rgrid%npt
  ivp=THIS%rgrid%find(rgrid%vpl(ir))
  if (ivp(ndim+1).lt.0.or.sum(abs(ivp(1:NDIM))).ne.0) cycle
  counter=counter+1
  if (trim(adjustl(THIS%mode)).eq.'tbfile'.or.trim(adjustl(THIS%mode)).eq.'hrfile') then
     deg(counter)=nint(THIS%dege(ivp(NDIM+1)))
  else
     deg(counter)=1
  end if
end do
call system("mkdir -p _ham")
open(50,file="_ham/ham_tb.dat",action='write')
write(50,*)
write(50,*) pars%avec(1,:)
write(50,*) pars%avec(2,:)
write(50,*) pars%avec(3,:)
write(50,*) THIS%norb_TB
write(50,*) counter
if (THIS%mode.eq."") then
  write(50,'(15I5)') (1,ir=1,counter)
else
  write(50,'(15I5)') (deg(ir),ir=1,counter)
end if
if (trim(adjustl(THIS%mode)).eq.'tbfile'.or.trim(adjustl(THIS%mode)).eq.'hrfile') then
   allocate(ham(THIS%norb_TB,THIS%norb_TB,THIS%rgrid%npt))
   call THIS%unfold_hame(ham)
end if
do ir=1,rgrid%npt
  ivp=THIS%rgrid%find(rgrid%vpl(ir))
  if (ivp(ndim+1).lt.0.or.sum(abs(ivp(1:NDIM))).ne.0) cycle
  jR=ivp(NDIM+1)
  write(50,*)
  write(50,'(3I5)') THIS%rgrid%vpi(jR)
  do jj=1,THIS%norb_TB
    do ii=1,THIS%norb_TB
      if (trim(adjustl(THIS%mode)).eq.'tbfile'.or.trim(adjustl(THIS%mode)).eq.'hrfile') then
        write(50,'(2I5,3x,2(E15.8,1x))') ii,jj,ham(ii,jj,jR)
      else
        ic=THIS%wbase%orb_icio(ii,1)
        jc=THIS%wbase%orb_icio(jj,1)
        dv=THIS%wbase%centers_cart(:,jc)+THIS%rgrid%vpc(jR)-THIS%wbase%centers_cart(:,ic)
        t1=sqrt(dot_product(dv,dv))
        if (t1.lt.THIS%rcut_nn) then
          write(50,'(2I5,3x,2(E15.8,1x))') ii,jj,THIS%tij(ii,jj,jR)
        else
          write(50,'(2I5,3x,2(E15.8,1x))') ii,jj,0._dp,0._dp
        end if
      end if
    end do
  end do
end do
close(50)
if (trim(adjustl(THIS%mode)).eq.'tbfile'.or.trim(adjustl(THIS%mode)).eq.'hrfile') then
  deallocate(ham)
end if
deallocate(deg)
return
end subroutine

subroutine read_hr_file(THIS,pars)
class(CLtb), intent(inout) :: THIS
class(CLpars), intent(in) :: pars
integer i,j,ii,jj
integer nrpt,nwan,ir
integer ngrid(3),ivp(4)
real(dp) aa,bb,vpl(3),avec(3,3)
logical exs
integer, allocatable :: ivr(:,:)
integer, allocatable :: deg(:)
complex(dp), allocatable :: ham(:,:,:)
if (NDIM.ne.3) call throw("CLtb%read_hr_file","this subroutine works only in 3D case")
inquire(file=trim(adjustl(pars%TBfile)),exist=exs)
if (.not.exs) call throw("CLtb%read_hr_file"," file "//trim(adjustl(pars%TBfile))//" missing")
open(50,file=trim(adjustl(pars%TBfile)),action='read')
read(50,*)
read(50,*) nwan
read(50,*) nrpt
if (nwan.ne.THIS%wbase%norb) then
  call throw("CLtb%read_hr_file","number of basis orbitals in _tb file is different from one derived from the input")
end if
allocate(ivr(3,nrpt))
allocate(deg(nrpt))
allocate(ham(nwan,nwan,nrpt))
read(50,'(15I5)') (deg(ir),ir=1,nrpt)
do ir=1,nrpt
  do j=1,nwan
    do i=1,nwan
      read(50,*) ivr(:,ir),ii,jj,aa,bb
      ham(ii,jj,ir)=cmplx(aa,bb,kind=dp)
    end do
  end do
end do
close(50)
ngrid(1)=max(2*maxval(abs(ivr(1,:))),1)+1
ngrid(2)=max(2*maxval(abs(ivr(2,:))),1)+1
ngrid(3)=max(2*maxval(abs(ivr(3,:))),1)+1
call THIS%rgrid%init(ngrid,avec,.true.,.false.)
allocate(THIS%dege(THIS%rgrid%npt))
allocate(THIS%hame_file(nwan,nwan,THIS%rgrid%npt))
THIS%dege=1._dp
THIS%hame_file=0._dp
do ir=1,nrpt
  vpl=dble(ivr(:,ir))
  ivp=THIS%rgrid%find(vpl)
  if (ivp(ndim+1).lt.0.or.sum(abs(ivp(1:NDIM))).ne.0) cycle
  THIS%dege(ivp(4))=dble(deg(ir))
  THIS%hame_file(:,:,ivp(4))=ham(:,:,ir)
end do
deallocate(ivr,deg,ham)
return
end subroutine

subroutine write_hr_file(THIS)
class(CLtb), intent(inout) :: THIS
integer ii,jj,ir
complex(dp), allocatable :: ham(:,:,:)
if (NDIM.ne.3) call throw("CLtb%write_hr_file","this subroutine works only in 3D case")
call system("mkdir -p _ham")
open(50,file="_ham/ham_hr.dat",action='write')
write(50,*)
write(50,*) THIS%norb_TB
write(50,*) THIS%rgrid%npt
if (THIS%mode.eq."") then
  write(50,'(15I5)') (1,ir=1,THIS%rgrid%npt)
else
  write(50,'(15I5)') (nint(THIS%dege(ir)),ir=1,THIS%rgrid%npt)
end if
allocate(ham(THIS%norb_TB,THIS%norb_TB,THIS%rgrid%npt))
call THIS%unfold_hame(ham)
do ir=1,THIS%rgrid%npt
  do jj=1,THIS%norb_TB
    do ii=1,THIS%norb_TB
      write(50,'(5I5,2F12.6)') ii,jj,ham(ii,jj,ir)
    end do
  end do
end do
close(50)
deallocate(ham)
return
end subroutine

subroutine unfold_hame(THIS,ham)
class(CLtb), intent(in) :: THIS
complex(dp), intent(out) :: ham(THIS%norb_TB,THIS%norb_TB,THIS%rgrid%npt)
integer ir,jr,ic,jc,nn,ii,jj,idx
integer ios,iorb,jos,jorb
real(dp) vpl(NDIM)
integer ivp(NDIM+1)
ham=0._dp
if (THIS%mode.ne."") then
  do ir=1,THIS%rgrid%npt
    do nn=1,THIS%hamsize
      do ii=1,THIS%norb_TB
        if (nn.ge.THIS%isa(ii).and.nn.lt.THIS%isa(ii+1)) then
          ham(ii,THIS%jsa(nn),ir)=THIS%hame(nn,ir)
        end if
      end do
    end do
    vpl=-THIS%rgrid%vpl(ir)
    ivp=THIS%rgrid%find(vpl)
    if (ivp(4).lt.0) then
      ham(:,:,ir)=0._dp
    else
      do ii=1,THIS%norb_TB
        do jj=ii+1,THIS%norb_TB
          ham(jj,ii,ivp(4))=conjg(ham(ii,jj,ir))
        end do
      end do
    end if
  end do
else
  do ir=1,THIS%rgrid%npt
    do ic=1,THIS%wbase%ncenters
      do nn=1,THIS%ncenters_nn(ic)
        jc=THIS%jcjr_nn(ic,nn,1)
        jr=THIS%jcjr_nn(ic,nn,2)
        if (ir.ne.jr) cycle
        do ios=1,THIS%wbase%norb_ic(ic)
          iorb=THIS%wbase%icio_orb(ic,ios)
          ! sparse index of the first non-zero element in the row iorb
          ii=THIS%isa(iorb)
          do jos=1,THIS%wbase%norb_ic(jc)
            jorb=THIS%wbase%icio_orb(jc,jos)
            do idx=ii,min(THIS%isa(iorb+1),THIS%hamsize)
               ! jsa is a map from sparse index to the column one
               if(THIS%jsa(idx).eq.jorb) then
                 ham(iorb,jorb,jr)=THIS%tij(iorb,jorb,jr)
                 exit
               end if
            end do
          end do
        end do
      end do
    end do
  end do
end if
return
end subroutine

subroutine write_nonzeros_hame(THIS)
class(CLtb), intent(inout) :: THIS
integer nn,ii,jj,ir
integer ic,jc
real(dp) t1,dv(NDIM)
allocate(THIS%hame(THIS%hamsize,THIS%rgrid%npt))
THIS%hame=0._dp
do ir=1,THIS%rgrid%npt
  do nn=1,THIS%hamsize
    do ii=1,THIS%norb_TB
      jj=THIS%jsa(nn)
      if (nn.ge.THIS%isa(ii).and.nn.lt.THIS%isa(ii+1)) then
        ic=THIS%wbase%orb_icio(ii,1)
        jc=THIS%wbase%orb_icio(jj,1)
        dv=THIS%wbase%centers_cart(:,jc)+THIS%rgrid%vpc(ir)-THIS%wbase%centers_cart(:,ic)
        t1=sqrt(dot_product(dv,dv))
        if (t1<THIS%rcut_nn) then
          THIS%hame(nn,ir)=THIS%hame_file(ii,THIS%jsa(nn),ir)
        end if
      end if
    end do
  end do
end do
deallocate(THIS%hame_file)
end subroutine


end module
