
module tasksclass
#ifdef MPI
  use mpi
#endif
use modcom
use parameters
use gridclass
use tbclass
use symmetryclass
use wannier_interface
implicit none
private
type, public :: CLtasks
   contains
   procedure, nopass :: run

end type

contains

subroutine run(pars)
class(CLpars), intent(inout) :: pars
integer itask
do itask=1,pars%ntasks
  if (trim(adjustl(pars%tasks(itask))).eq."bands") then
     call message("*** bandstructure calculation ***")
     call message("")
     call bands(pars)
     call message("")
  else if (trim(adjustl(pars%tasks(itask))).eq."eigen") then
     call message("*** eigenvalues  calculation ***")
     call message("")
     call eigen(pars)
     call message("")
  else if (trim(adjustl(pars%tasks(itask))).eq."valleyProjectionBands") then
     call message("*** valleyProjectionBands calculation ***")
     call message("")
     call valleyProjectionBands(pars)
     call message("")
  else if (trim(adjustl(pars%tasks(itask))).eq."dos") then
     call message("*** dos  calculation ***")
     call message("")
     call calc_dos(pars)
     call message("")
  else if (trim(adjustl(pars%tasks(itask))).eq."rpachi") then
     call message("*** rpachi  calculation ***")
     call message("")
     call calc_chi(pars,"rpa")
     call message("")
  else if (trim(adjustl(pars%tasks(itask))).eq."projection_wannier") then
     call message("*** wannier projection ***")
     call message("")
     call projection_wannier(pars)
     call message("")
  else if (trim(adjustl(pars%tasks(itask))).eq."symreps") then
     call message("*** symmetry representations of bands in kpath vertices ***")
     call message("")
     call symreps(pars)
     call message("")
  else if (trim(adjustl(pars%tasks(itask))).eq."write_wfmloc") then
     call message("*** extract wfmloc ***")
     call message("")
     call write_wfmloc(pars)
     call message("")
  else if (trim(adjustl(pars%tasks(itask))).eq."hubbard_tbg") then
     call message("*** compute Hubbard parameters ***")
     call message("")
     call hubbard_tbg(pars)
     call message("")
  else if (trim(adjustl(pars%tasks(itask))).eq."write_hubbardu") then
     call message("*** extract Hubbard parameters ***")
     call message("")
     call write_hubu(pars)
     call message("")
  else if (trim(adjustl(pars%tasks(itask))).eq."symmetrize_hubbardu") then
     call message("*** symmetrize Hubbard parameters ***")
     call message("")
     call symmetrize_hubu(pars)
     call message("")
  else if (trim(adjustl(pars%tasks(itask))).eq."symmetrize_tb") then
     call message("*** symmetrize Tight-Binding Hamiltonian ***")
     call message("")
     call symmetrize_tb(pars)
     call message("")
  else
    call throw("CLtasks%run","unknown task")
  end if
end do
end subroutine

subroutine valleyProjectionBands(pars)
class(CLpars), intent(inout) :: pars
integer ik,ist
real(dp) vv(NDIM)
real(dp), dimension(NDIM,NDIM) :: avec1, avec2, bvec1, bvec2
real(dp), allocatable :: eval(:,:)
real(dp), allocatable :: vkl(:,:)
complex(dp), allocatable :: evec(:,:,:)
type(PATH) kpath
type(CLtb) tbmodel
type(CLsym) sym
call info('taskclass%valleyProjectionBands()','this part is hard-coded for 1.05-TBG, istruct=131')

! hard-coded part
! lattice vectors of the first TBG layer
avec1 = transpose(reshape((/2.119061,-1.249471,   0.000000,&
                            2.141605, 1.210425,   0.000000,&
                            0.000000, 0.000000, 137.572253/),shape(avec1)))

! lattice vectors of the second TBG layer
avec2 = transpose(reshape((/2.141605,-1.210425,   0.000000,&
                            2.119061, 1.249471,   0.000000,&
                            0.000000, 0.000000, 137.572253/),shape(avec2)))

! invert to get transposed reciprocal vectors, without 2pi factor
call dmatrix_inverse(avec1, bvec1, NDIM)
call dmatrix_inverse(avec2, bvec2, NDIM)

#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
! generater spatial symmetries
call sym%init(pars)
! init new kpath from input
call kpath%init(pars%nvert,pars%np_per_vert,pars%vert,pars%bvec)
! initialise TB model, to have centers coordinates and other data
call tbmodel%init(pars,sym,"SK")
! allocate array for eigen values
allocate(eval(pars%nstates,kpath%npt))
! allocate array for eigen vectors
allocate(evec(tbmodel%norb_TB,pars%nstates,kpath%npt))
! this is needed to copy the private data of kgrid object, i.e., k-points in lattice coordinates
allocate(vkl(NDIM,kpath%npt))
do ik=1,kpath%npt
  ! copy the private data
  vkl(:,ik)=kpath%vpl(ik)
  ! read eigenvectors, subroutine in modcom.f90
  call io_evec(ik,"read",'_bandsevec_',tbmodel%norb_TB,pars%nstates,evec(:,:,ik))
end do
open(100,file="bands_test.dat",action="write")
do ist=1,pars%nstates
  do ik=1,kpath%npt
     write(100,*) kpath%dvpc(ik),eval(ist,ik)
   end do
   write(100,*)
end do
close(100)

! zero the arrays for security reasons
eval=0._dp
! read eigenvalues
open(100,file="bands.dat",action="read")
do ist=1,pars%nstates
  do ik=1,kpath%npt
     read(100,*) vv(1),eval(ist,ik)
   end do
   read(100,*)
end do
close(100)

open(100,file="valleyProjectionBands.dat",action="write")
do ist=1,pars%nstates
  do ik=1,kpath%npt
     write(100,*) kpath%dvpc(ik), tbmodel%antiHaldaneME(vkl(:,ik), evec(:,ist,ik), bvec1, bvec2)
   end do
   write(100,*)
end do
close(100)


end subroutine

subroutine bands(pars)
class(CLpars), intent(inout) :: pars
type(PATH) kpath
type(CLtb) tbmodel
type(CLsym) sym
integer ik,ist
#ifdef MPI
  integer nn
#endif
real(dp), allocatable :: eval(:,:)
complex(dp), allocatable :: evec(:,:)
call message("  initialise symmetries ..")
call message("")
call sym%init(pars)
! init k-point path
call kpath%init(pars%nvert,pars%np_per_vert,pars%vert,pars%bvec)
call message("  initialise TB model ..")
call message("")
call tbmodel%init(pars,sym,"SK")
allocate(eval(pars%nstates,kpath%npt))
allocate(evec(tbmodel%norb_TB,pars%nstates))
eval=0._dp
evec=0._dp

call message("  Eigenproblem calculation. some k-points progress below ..")

#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif

do ik=1,kpath%npt

   if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
   if (mod(ik-1,max(np_mpi,kpath%npt/10)).eq.0) write(*,*) "ik, of: ",ik,kpath%npt
   call tbmodel%evalk(pars%bands_evec,pars,kpath%vpl(ik),eval(:,ik),evec)
   if (pars%bands_evec) call io_evec(ik,"write",'_bandsevec_',tbmodel%norb_TB,pars%nstates,evec)

end do

#ifdef MPI
  nn=pars%nstates*kpath%npt
  call mpi_allreduce(mpi_in_place,eval,nn,mpi_double_precision,mpi_sum, &
   mpi_com,mpi_err)
#endif

if (mp_mpi) then
  open(100,file="bands.dat",action="write")
  do ist=1,pars%nstates
    do ik=1,kpath%npt
       write(100,*) kpath%dvpc(ik),eval(ist,ik)-pars%efermi
     end do
     write(100,*)
  end do
  close(100)
  open(101,file="bandlines.dat",action="write")
  do ik=1,kpath%nvert
    write(101,*) kpath%dvrt(ik),minval(eval-pars%efermi)
    write(101,*) kpath%dvrt(ik),maxval(eval-pars%efermi)
    write(101,*)
  end do
  close(101)
end if
deallocate(eval)
deallocate(evec)
#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
end subroutine

subroutine eigen(pars)
class(CLpars), intent(inout) :: pars
type(GRID) kgrid
type(CLtb) tbmodel
type(CLsym) sym
integer ik, ij
character(20) filename
real(dp), allocatable :: vkl(:,:)
#ifdef MPI
  integer nn
#endif
real(dp), allocatable :: eval(:,:)
complex(dp), allocatable :: evec(:,:)
call message("  initialise symmetries ..")
call message("")
call sym%init(pars)
! init k-grid
call kgrid%init(pars%ngrid,pars%bvec,centered_kgrid,.true.)
call message("  initialise TB model ..")
call message("")
call tbmodel%init(pars,sym,"SK")
allocate(eval(pars%nstates,kgrid%npt))
allocate(evec(tbmodel%norb_TB,pars%nstates))
call message("  Eigenproblem calculation. some k-points progress below ..")

! shifted grids
do ij=1,pars%nkshift
  ! syncronize all mpi processes, because the root process is a bit behind after write of eval below
#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif

  eval=0._dp
  evec=0._dp
  do ik=1,kgrid%npt
     if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
     if (mod(ik-1,max(np_mpi,kgrid%npt/10)).eq.0) write(*,*) "ik, of: ",ik,kgrid%npt
     call tbmodel%evalk(.true.,pars,kgrid%vpl(ik)+pars%kshift(:,ij),eval(:,ik),evec)
     write(filename,'(a,i0,a)') '_evec_',ij,'_'
     call io_evec(ik,"write",filename,tbmodel%norb_TB,pars%nstates,evec)
  end do

#ifdef MPI
  nn=pars%nstates*kgrid%npt
  call mpi_allreduce(mpi_in_place,eval,nn,mpi_double_precision,mpi_sum, &
   mpi_com,mpi_err)
#endif

  if (mp_mpi) then
    allocate(vkl(NDIM,kgrid%npt))
    do ik=1,kgrid%npt
      vkl(:,ik)=kgrid%vpl(ik)
    end do
    write(filename,'(a,i0,a)') 'eval_', ij, '.dat'
    call io_eval(1001,"write",filename,.false.,pars%nstates,kgrid%npt,pars%efermi,vkl,eval)
    deallocate(vkl)
  end if
end do

if (mp_mpi) then
  call kgrid%io(1003,"_grid","write",pars,tbmodel%norb_TB)
end if

deallocate(eval,evec)

#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
end subroutine


subroutine calc_dos(pars)
class(CLpars), intent(inout) :: pars
integer ik
real(dp), allocatable :: eval(:,:)
real(dp), allocatable :: vkl(:,:)
real(dp), allocatable :: tdos(:),dosk(:,:)
type(GRID) kgrid
type(CLtb) tbmodel
type(CLsym) sym
integer ie
#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
call sym%init(pars)
! everytime after symmetry run, we have to consider a shift to a more symmetric place
call tbmodel%init(pars,sym,"noham")
call kgrid%io(100,"_grid","read",pars,tbmodel%norb_TB)
allocate(eval(pars%nstates,kgrid%npt))
allocate(vkl(NDIM,kgrid%npt))
allocate(tdos(pars%negrid),dosk(pars%negrid,kgrid%npt))
do ik=1,kgrid%npt
 vkl(:,ik)=kgrid%vpl(ik)
end do
eval=0._dp
call io_eval(1001,"read","eval_1.dat",.false.,pars%nstates,kgrid%npt,pars%efermi,vkl,eval)
! shift all eigenvalues by efermi
eval=eval-pars%efermi
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
do ik=1,kgrid%npt
  call get_dosk(pars%nstates,pars,eval(:,ik),dosk(:,ik))
end do
!$OMP END DO
!$OMP END PARALLEL
tdos=0._dp
do ie=1,pars%negrid
  tdos(ie)=sum(dosk(ie,:))/dble(kgrid%npt)
end do
if (mp_mpi) then
  open(100,file="dos.dat")
  do ie=1,pars%negrid
    write(100,'(2G18.10)') pars%egrid(ie),tdos(ie)
  end do
end if
deallocate(dosk,tdos)
deallocate(eval,vkl)
#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
end subroutine

subroutine calc_chi(pars,opt)
class(CLpars), intent(inout) :: pars
character(len=*), intent(in) :: opt
integer ik,iq,id,ig
real(dp) vc(NDIM),dc
complex(dp), allocatable :: chi(:,:,:,:)
type(GRID) kgrid,qgrid,Ggrid
type(CLtb) tbmodel
type(CLsym) sym
#ifdef MPI
  integer nn
#endif

#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif

call sym%init(pars)
if (trim(adjustl(opt)).ne."rpa") then
  call throw("CLtasks","unknown argument for the response function (only RPA is available via rpachi task)")
end if
! initialise TB model, to have centers coordinates and other data
call tbmodel%init(pars,sym,"noham")
! read the k-point grid on which eigenvales/eigenvectors are computed
call kgrid%io(1000,"_grid","read",pars,tbmodel%norb_TB)
! init q-grid it has to be commensurate and less or equal to k-grid
do id=1,NDIM
  if (pars%qgrid(id).gt.pars%ngrid(id)) call throw("CLtasks%calc_chi","q-grid must not be larger than k-grid in any dimension")
  if (mod(pars%ngrid(id),pars%qgrid(id)).ne.0) call throw("CLtasks%calc_chi","q-grid must commensurate with k-grid")
end do
call qgrid%init(pars%qgrid,pars%bvec,centered_qgrid,.true.)
! grid of reciprocal lattice points G such that q_FBZ+G samples whole reciprocal space
if (sum(abs(pars%Ggrid)).eq.0) pars%Ggrid=pars%ngrid
call Ggrid%init(pars%Ggrid,pars%bvec,.true.,.false.)
! take remember its spherical part only, defined by pars%rcut_grid
call Ggrid%init_sphere(pars)

! array for RPA response function which can be dynamic (on pars%negrid frequencies), and at different q-points
allocate(chi(pars%negrid,pars%nkshift,Ggrid%npt_sphere,qgrid%npt))

if (mp_mpi) then
  print *, "Calculating RPA Chi"
  print *, "on a q grid of", qgrid%ngrid
  print *, "with BZ grid of", kgrid%ngrid
  print *, "Using ", pars%nstates, "states"
  print *, "with", tbmodel%norb_TB, "orbitals per unit cell"
  if (pars%chi_exclude) then
    print '(a,F7.5,a,F7.5)', "Excluding bands in the range  ", pars%e_chi_exclude(1), " < E - E_fermi < ", pars%e_chi_exclude(2)
  end if
  print *, "Restricting to states from ", pars%chi_start, " to ", pars%chi_stop
end if

! do the computation. later we will attach MPI parallelisation here
! (you can see bands, eigen tasks how to do it), therefore arrays have to be zeroed
chi(:,:,:,:)=0._dp
do iq=1,qgrid%npt
  if (mod(iq-1,np_mpi).ne.lp_mpi) cycle
  if (mod(iq-1,max(np_mpi,qgrid%npt/10)).eq.0) write(*,*) "ik, of: ",iq,qgrid%npt
  call get_chiq(iq,qgrid,kgrid,Ggrid,pars,tbmodel,chi(:,:,:,iq))
end do

#ifdef MPI
  nn=pars%negrid*Ggrid%npt_sphere*qgrid%npt*pars%nkshift
  call mpi_allreduce(mpi_in_place,chi,nn,MPI_DOUBLE_COMPLEX,mpi_sum, &
   mpi_com,mpi_err)
#endif

if (mp_mpi) then
  open(2000,file="chi.dat")
  do ik=1, pars%nkshift
    do iq=1,qgrid%npt
      do ig=1,Ggrid%npt_sphere
        vc=qgrid%vpc(iq)+Ggrid%vpc_sphere(ig)+ MATMUL(TRANSPOSE(pars%bvec), pars%kshift(:,ik))
        dc=sqrt(dot_product(vc,vc))
        write(2000,'(20F16.8)') dc,REAL(chi(1,ik,ig,iq)),qgrid%vpl(iq)+Ggrid%vpl_sphere(ig)+pars%kshift(:,ik)
      end do
    end do
  end do
end if

#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
end subroutine

subroutine projection_wannier(pars)
class(CLpars), intent(inout) :: pars
integer ik,jk,ir
real(dp), allocatable :: eval(:,:)
real(dp), allocatable :: vkl(:,:)
complex(dp), allocatable :: evec(:,:,:)
type(GRID) kgrid
type(PATH) kpath
type(CLtb) tbmodel
type(CLwan) wannier
type(CLsym) sym
#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
! generater spatial symmetries
call sym%init(pars)
! initialise TB model, to have centers coordinates and other data
call tbmodel%init(pars,sym,"noham")
! read the k-point grid on which eigenvales/eigenvectors are computed
call kgrid%io(1000,"_grid","read",pars,tbmodel%norb_TB)
! init symmetries of k-point grid
call kgrid%sym_init(pars%trev,sym)
! init new kpath from input
call kpath%init(pars%nvert,pars%np_per_vert,pars%vert,pars%bvec)
! allocate array for eigen values
allocate(eval(pars%nstates,kgrid%npt))
! allocate array for eigen vectors
allocate(evec(tbmodel%norb_TB,pars%nstates,kgrid%npt))
! this is needed to copy the private data of kgrid object, i.e., k-points in lattice coordinates
allocate(vkl(NDIM,kgrid%npt))
do ik=1,kgrid%npt
  ! copy the private data
  vkl(:,ik)=kgrid%vpl(ik)
  ! read eigenvectors, subroutine in modcom.f90
  call io_evec(ik,"read","_evec_1_",tbmodel%norb_TB,pars%nstates,evec(:,:,ik))
end do
! symetrise with respect to time-reversal symmetry
! at some point, this code is to be removed, and eigenvalues 
! are to be computed on irreducible wedge
do ir=1,kgrid%nirT
  ik=kgrid%irT2ik(ir)
  jk=kgrid%ikT2k(ik)
  evec(:,:,ik)=conjg(evec(:,:,jk))
end do
! zero the arrays for security reasons
eval=0._dp
! read eigenvalues, subroutine in modcom.f90
call io_eval(1001,"read","eval_1.dat",.false.,pars%nstates,kgrid%npt,pars%efermi,vkl,eval)
! init minimal wannier variables
call wannier%init(kgrid,kpath,pars,eval)
! do the computation. later we will attach MPI parallelisation here
! (you can see bands, eigen tasks how to do it), therefore arrays have to be zeroed
call wannier%projection(tbmodel,pars,sym,kgrid,eval-pars%efermi,evec)
!call wannier%overlap(pars,tbmodel,kgrid,eval,evec)
if (mp_mpi) then
end if
deallocate(eval,evec,vkl)
#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
end subroutine

subroutine symreps(pars)
class(CLpars), intent(inout) :: pars
integer ik
real(dp), allocatable :: eval(:)
complex(dp), allocatable :: evec(:,:)
type(CLtb) tbmodel
type(PATH) kpath
type(CLsym) sym
type(char_table) ch_table
real(dp) vpc(NDIM)
integer isym,ngroups,igr,ist,icl,nsym
complex(dp) ch_bandgroup,irrep_decompose(48)
integer, allocatable :: idx(:,:)
complex(dp), allocatable :: wft(:)
#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
! generater spatial symmetries
call sym%init(pars)
call ch_table%read_character_table(pars%character_file,sym%nsym)
call kpath%init(pars%nvert,pars%np_per_vert,pars%vert,pars%bvec)
call tbmodel%init(pars,sym,"SK")
allocate(idx(pars%nstates,2))
allocate(eval(pars%nstates))
allocate(evec(tbmodel%norb_TB,pars%nstates))
allocate(wft(tbmodel%norb_TB))
eval=0._dp
evec=0._dp
do ik=1,kpath%nvert
   write(*,*)"IK: ",ik
   call tbmodel%evalk(.true.,pars,kpath%vert(:,ik),eval,evec)
   vpc=matmul(kpath%vert(:,ik),kpath%vecs)
   call find_degroups(pars%nstates,eval,ngroups,idx,1.d-2)
   ! find symreps for groups of bands
   do igr=1,ngroups
     write(*,*) 'group: ',igr,idx(igr,1),idx(igr,2)
     write(*,'("energies: ",20F12.6)') eval(idx(igr,1):idx(igr,2))-pars%efermi
     irrep_decompose(:)=0._dp
     do icl=1,ch_table%nclasses
       do nsym=1,ch_table%nsym_of_cl(icl)
         isym=ch_table%isym_of_cl(icl,nsym)
         ch_bandgroup=0
         do ist=idx(igr,1),idx(igr,2)
           call tbmodel%bloch_wf_transform_statk(vpc,vpc,sym,isym,wft,evec(:,ist))
           ch_bandgroup=ch_bandgroup+dot_product(evec(:,ist),wft)
         end do
         irrep_decompose(1:ch_table%nclasses)=irrep_decompose(1:ch_table%nclasses)&
                                             +ch_bandgroup*ch_table%table(1:ch_table%nclasses,icl)/dble(sym%nsym)
       end do
     end do
     write(*,'("irreps_decomposition: ")',advance='no')
     write(*,'(10A5)',advance='no') ch_table%symreps(1:ch_table%nclasses)
     write(*,'(100F10.4)') abs(irrep_decompose(1:ch_table%nclasses))
   end do
end do
#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
deallocate(eval,evec,idx)
end subroutine

subroutine write_wfmloc(pars)
class(CLpars), intent(inout) :: pars
real(dp), allocatable :: eval(:,:)
real(dp), allocatable :: vkl(:,:)
complex(dp), allocatable :: evec(:,:,:)
complex(dp), allocatable :: wfmloc(:,:,:)
type(GRID) kgrid
type(CLtb) tbmodel
type(CLsym) sym
integer nr,iR,iw,iorb,jk,ik
#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
! generater spatial symmetries
call sym%init(pars)
! initialise TB model, to have centers coordinates and other data
call tbmodel%init(pars,sym,"noham")
! read the k-point grid on which eigenvales/eigenvectors are computed
call kgrid%io(1000,"_grid","read",pars,tbmodel%norb_TB)
! init symmetries of k-point grid
call kgrid%sym_init(pars%trev,sym)
! allocate array for eigen values
allocate(eval(pars%nstates,kgrid%npt))
! allocate array for eigen vectors
allocate(evec(tbmodel%norb_TB,pars%nstates,kgrid%npt))
! this is needed to copy the private data of kgrid object, i.e., k-points in lattice coordinates
allocate(vkl(NDIM,kgrid%npt))
do ik=1,kgrid%npt
  ! copy the private data
  vkl(:,ik)=kgrid%vpl(ik)
  ! read eigenvectors, subroutine in modcom.f90
  call io_evec(ik,"read","_evec_1_",tbmodel%norb_TB,pars%nstates,evec(:,:,ik))
end do

! symetrise with respect to time-reversal symmetry
! at some point, this code is to be removed, and eigenvalues 
! are to be computed on irreducible wedge
do ir=1,kgrid%nirT
  ik=kgrid%irT2ik(ir)
  jk=kgrid%ikT2k(ik)
  evec(:,:,ik)=conjg(evec(:,:,jk))
end do

! zero the arrays for security reasons
eval=0._dp
! read eigenvalues, subroutine in modcom.f90
call io_eval(1001,"read","eval_1.dat",.false.,pars%nstates,kgrid%npt,pars%efermi,vkl,eval)
! init minimal wannier variables
!call wannier%init(kgrid,kpath,pars,eval)

allocate(wfmloc(tbmodel%norb_TB,pars%proj%norb,tbmodel%rgrid%npt))
call read_wfmloc(pars,tbmodel,kgrid,eval-pars%efermi,evec,wfmloc)
nr=0
do iR=1,tbmodel%rgrid%npt
  if (sum(abs(tbmodel%rgrid%vpl(iR))).le.6) then
    nr=nr+1
    do iw=1,pars%proj%norb
      do iorb=1,tbmodel%norb_TB
        wfmloc(iorb,iw,nr)=wfmloc(iorb,iw,iR)
      end do
    end do
  end if
end do
if (mp_mpi) call write_wf_universe(tbmodel,pars,nr,wfmloc,'wfmloc','')
deallocate(eval,evec,vkl,wfmloc)
end subroutine

subroutine hubbard_tbg(pars)
class(CLpars), intent(inout) :: pars
integer ik,ir,jk
real(dp), allocatable :: eval(:,:)
real(dp), allocatable :: vkl(:,:)
complex(dp), allocatable :: evec(:,:,:)
type(GRID) kgrid
type(CLtb) tbmodel
type(CLsym) sym
#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
! generater spatial symmetries
call sym%init(pars)
! initialise TB model, to have centers coordinates and other data
call tbmodel%init(pars,sym,"noham")
! read the k-point grid on which eigenvales/eigenvectors are computed
call kgrid%io(1000,"_grid","read",pars,tbmodel%norb_TB)
! init symmetries of k-point grid
call kgrid%sym_init(pars%trev,sym)
! allocate array for eigen values
allocate(eval(pars%nstates,kgrid%npt))
! allocate array for eigen vectors
allocate(evec(tbmodel%norb_TB,pars%nstates,kgrid%npt))
! this is needed to copy the private data of kgrid object, i.e., k-points in lattice coordinates
allocate(vkl(NDIM,kgrid%npt))
do ik=1,kgrid%npt
  ! copy the private data
  vkl(:,ik)=kgrid%vpl(ik)
  ! read eigenvectors, subroutine in modcom.f90
  call io_evec(ik,"read","_evec_1_",tbmodel%norb_TB,pars%nstates,evec(:,:,ik))
end do
! symetrise with respect to time-reversal symmetry
! at some point, this code is to be removed, and eigenvalues 
! are to be computed on irreducible wedge
do ir=1,kgrid%nirT
  ik=kgrid%irT2ik(ir)
  jk=kgrid%ikT2k(ik)
  evec(:,:,ik)=conjg(evec(:,:,jk))
end do
! zero the arrays for security reasons
eval=0._dp
! read eigenvalues, subroutine in modcom.f90
call io_eval(1001,"read","eval_1.dat",.false.,pars%nstates,kgrid%npt,pars%efermi,vkl,eval)
! init minimal wannier variables
call compute_hubbardu_rs(pars,sym,tbmodel,kgrid,eval-pars%efermi,evec)
!call compute_hubbardj_rs(pars,sym,tbmodel,kgrid,sym,eval-pars%efermi,evec)
deallocate(eval,evec,vkl)
end subroutine

subroutine write_hubu(pars)
class(CLpars), intent(inout) :: pars
real(dp), allocatable :: eval(:,:)
real(dp), allocatable :: vkl(:,:)
type(GRID) kgrid
type(CLtb) tbmodel
type(CLsym) sym
integer ik
#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
! generater spatial symmetries
call sym%init(pars)
! initialise TB model, to have centers coordinates and other data
call tbmodel%init(pars,sym,"noham")
! read the k-point grid on which eigenvales/eigenvectors are computed
call kgrid%io(1000,"_grid","read",pars,tbmodel%norb_TB)
! allocate array for eigen values
allocate(eval(pars%nstates,kgrid%npt))
! this is needed to copy the private data of kgrid object, i.e., k-points in lattice coordinates
allocate(vkl(NDIM,kgrid%npt))
do ik=1,kgrid%npt
  ! copy the private data
  vkl(:,ik)=kgrid%vpl(ik)
end do
! zero the arrays for security reasons
eval=0._dp
! read eigenvalues, subroutine in modcom.f90
call io_eval(1001,"read","eval_1.dat",.false.,pars%nstates,kgrid%npt,pars%efermi,vkl,eval)
! init minimal wannier variables
call sym%init(pars)
call write_hubbardu(pars,sym)
deallocate(eval,vkl)
end subroutine

subroutine symmetrize_hubu(pars)
class(CLpars), intent(inout) :: pars
real(dp), allocatable :: eval(:,:)
real(dp), allocatable :: vkl(:,:)
type(GRID) kgrid
type(CLtb) tbmodel
type(CLsym) sym
integer ik
#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
! generater spatial symmetries
call sym%init(pars)
! initialise TB model, to have centers coordinates and other data
call tbmodel%init(pars,sym,"noham")
! read the k-point grid on which eigenvales/eigenvectors are computed
call kgrid%io(1000,"_grid","read",pars,tbmodel%norb_TB)
! allocate array for eigen values
allocate(eval(pars%nstates,kgrid%npt))
! this is needed to copy the private data of kgrid object, i.e., k-points in lattice coordinates
allocate(vkl(NDIM,kgrid%npt))
do ik=1,kgrid%npt
  ! copy the private data
  vkl(:,ik)=kgrid%vpl(ik)
end do
! zero the arrays for security reasons
eval=0._dp
! read eigenvalues, subroutine in modcom.f90
call io_eval(1001,"read","eval_1.dat",.false.,pars%nstates,kgrid%npt,pars%efermi,vkl,eval)
! init minimal wannier variables
call sym%init(pars)
call symmetrize_hubbardu(pars,sym)
deallocate(eval,vkl)
end subroutine


subroutine symmetrize_tb(pars)
class(CLpars), intent(inout) :: pars
type(CLsym) sym
call sym%init(pars)
call symmetrize_tbfile(pars,sym)
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Side subroutines, could be placed into another file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_chiq(iq,qgrid,kgrid,Ggrid,pars,tbmodel,chiq)
integer, intent(in) :: iq
class(GRID), intent(in) :: qgrid,kgrid,Ggrid
class(CLpars), intent(in) :: pars
class(CLtb), intent(in) :: tbmodel
complex(dp), intent(out) :: chiq(pars%negrid,pars%nkshift,Ggrid%npt_sphere)

! local
character(20) filename
integer ik, ig, ij
integer ikg(NDIM+1)
integer ic, iv, iorb
integer fermi_index_kq(1)
integer fermi_index_k(1)
real(dp) ec, ev
real(dp) delta,pwo,dc
real(dp) vkq(NDIM),vl(NDIM),vc(NDIM)
complex(dp) overlap
complex(dp), allocatable :: eitqG(:,:)
complex(dp), allocatable :: eveck(:,:),eveckq(:,:)
real(dp), allocatable :: eval(:,:,:)
real(dp), allocatable :: vkl(:,:)

! allocate array for eigen vectors
allocate(eveck(tbmodel%norb_TB,pars%nstates))
allocate(eveckq(tbmodel%norb_TB,pars%nstates))
! allocate array for eigenvalues
allocate(eval(pars%nkshift,pars%nstates,kgrid%npt))
allocate(vkl(NDIM,kgrid%npt))
! allocate array for phase factors
allocate(eitqG(tbmodel%norb_TB,Ggrid%npt_sphere))

do ik=1,kgrid%npt
  ! copy the private data
  vkl(:,ik)=kgrid%vpl(ik)
end do

eval = 0._dp
do ik=1, pars%nkshift
  write(filename,'(a,i0,a)') 'eval_', ik, '.dat'
  call io_eval(1001,"read",filename,.false.,pars%nstates,kgrid%npt,pars%efermi,vkl,eval(ik,:,:))
end do
eval=eval-pars%efermi

chiq = 0._dp
do ij=1, pars%nkshift
  ! array for phase factors inside matrix elements
  do ig=1,Ggrid%npt_sphere
    vl=qgrid%vpl(iq)+Ggrid%vpl_sphere(ig)+pars%kshift(:,ij)
    do iorb=1,tbmodel%norb_tb
      ! compute phase factors e^(i(q+G)t) for all t, atomic positions
      eitqG(iorb,ig)=EXP(CMPLX(0._dp, twopi*dot_product(vl,tbmodel%vplorb(iorb)) ,kind=dp))
    end do
  end do

  do ik=1,kgrid%npt
    vkq=qgrid%vpl(iq)+kgrid%vpl(ik)
    ikg=kgrid%find(vkq)

    ! read eigenvectors, subroutine in modcom.f90
    write(filename,'(a,i0,a)') '_evec_',1,'_'
    call io_evec(ik,"read",filename,tbmodel%norb_TB,pars%nstates,eveck)
    write(filename,'(a,i0,a)') '_evec_',ij,'_'
    call io_evec(ikg(4),"read",filename,tbmodel%norb_TB,pars%nstates,eveckq)

    ! Get the index corresponding to fermi level
    fermi_index_kq =  minloc(eval(ij,:,ikg(4)), mask=(eval(ij,:,ikg(4)) > 0))
    fermi_index_k =  minloc(eval(1,:,ik), mask=(eval(1,:,ik) > 0))
    fermi_index_kq(1) = merge(fermi_index_kq(1)-1,pars%nstates, fermi_index_kq(1) > 0)
    fermi_index_k(1) = merge(fermi_index_k(1)-1,pars%nstates, fermi_index_k(1) > 0)
    
    !$OMP PARALLEL DEFAULT(SHARED)&
    !$OMP PRIVATE(ig,vc,dc,pwo,iv,ic,overlap,delta)
    !$OMP DO
    do ig=1,Ggrid%npt_sphere
      vc=qgrid%vpc(iq)+Ggrid%vpc_sphere(ig)
      dc=sqrt(dot_product(vc,vc))
      ! overlap I(q)=<local_orb|e^iqr|local_orb>, currently only the pz orbital
      pwo = pwave_ovlp(dc)
      do iv=pars%chi_start, fermi_index_k(1)
        do ic=fermi_index_kq(1)+1, pars%chi_stop
          ev = eval(1, iv, ik)
          ec = eval(ij, ic, ikg(4))
          if (pars%chi_exclude) then
            if (((ev < pars%e_chi_exclude(2)) .and. (ev > pars%e_chi_exclude(1))) &
            .AND. ((ec < pars%e_chi_exclude(2)) .and. (ec > pars%e_chi_exclude(1)))) cycle
          end if
          if (pars%ignore_chiIq) then
            overlap = ABS(DOT_PRODUCT(eveckq(1:tbmodel%norb_TB,ic), eitqG(:,ig)*eveck(1:tbmodel%norb_TB,iv)))
          else
            overlap = ABS(DOT_PRODUCT(eveckq(1:tbmodel%norb_TB,ic), eitqG(:,ig)*eveck(1:tbmodel%norb_TB,iv))) * pwo
          end if
          delta = ec - ev
          if (delta > 1e-7) then
            chiq(1,ij,ig) = chiq(1,ij,ig) + overlap*overlap/delta
          end if
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end do
end do

! Normalise
chiq = (4.0/(kgrid%npt)) * chiq

deallocate(eveck,eveckq,eitqG,eval)
end subroutine


subroutine get_dosk(nstates,pars,eval,dosk)
integer, intent(in) :: nstates
class(CLpars), intent(in) :: pars
real(dp), intent(in) :: eval(nstates)
real(dp), intent(out) :: dosk(pars%negrid)
! local
integer ie,ist,nstates_window
integer, allocatable :: idx(:)
real(dp), allocatable :: ekw(:),ee(:)
! first serach for all states within a window emin,emax
allocate(idx(nstates))
nstates_window=0
do ist=1,nstates
  if (eval(ist).gt.pars%emin.and.eval(ist).lt.pars%emax) then
    nstates_window=nstates_window+1
    idx(nstates_window)=ist
  end if
end do
allocate(ekw(nstates_window),ee(nstates_window))
do ist=1,nstates_window
  ekw(ist)=eval(idx(ist))
end do
dosk(:)=0._dp
do ie=1,pars%negrid
  do ist=1,nstates_window
    ee(ist)=exp( -( (ekw(ist)-pars%egrid(ie))/pars%gauss_sigma )**2 ) &
         / pars%gauss_sigma / sqrtpi
  end do
  dosk(ie)=sum(ee)
end do
deallocate(idx,ekw,ee)
end subroutine


end module
