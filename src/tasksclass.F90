
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
  else
    call throw("CLtasks%run","unknown task")
  end if
end do
end subroutine

subroutine bands(pars)
class(CLpars), intent(inout) :: pars
type(PATH) kpath
type(CLtb) tbmodel
type(CLsym) sym
integer ik,iq,ist
#ifdef MPI
  integer nn
#endif
real(dp), allocatable :: eval(:,:)
complex(dp), allocatable :: evec(:,:)
call kpath%init(pars%nvert,pars%np_per_vert,pars%vert,pars%bvec)
call message("  initialise TB model ..")
call message("")
call tbmodel%init(pars,"SK")
call message("  initialise symmetries ..")
call message("")
call sym%init(pars)
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
   call tbmodel%evalk(.false.,pars,kpath%vpl(ik),eval(:,ik),evec)
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
integer ik
real(dp), allocatable :: vkl(:,:)
#ifdef MPI
  integer nn
#endif
real(dp), allocatable :: eval(:,:)
complex(dp), allocatable :: evec(:,:)
call kgrid%init(pars%ngrid,pars%ngrid,pars%bvec,centered_kgrid,.true.)
call message("  initialise TB model ..")
call message("")
call tbmodel%init(pars,"SK")
call message("  initialise symmetries ..")
call message("")
call sym%init(pars)
allocate(eval(pars%nstates,kgrid%npt))
allocate(evec(tbmodel%norb_TB,pars%nstates))
eval=0._dp
evec=0._dp
call message("  Eigenproblem calculation. some k-points progress below ..")
#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
do ik=1,kgrid%npt
   if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
   if (mod(ik-1,max(np_mpi,kgrid%npt/10)).eq.0) write(*,*) "ik, of: ",ik,kgrid%npt
   call tbmodel%evalk(.true.,pars,kgrid%vpl(ik),eval(:,ik),evec)
   call io_evec(ik,"write","_evec",tbmodel%norb_TB,pars%nstates,evec)
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
  call io_eval(1001,"write","eval.dat",.false.,pars%nstates,kgrid%npt,pars%efermi,vkl,eval)
  call kgrid%io(1003,"_grid","write",pars,tbmodel%norb_TB)
  deallocate(vkl)
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
integer ie
#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
call tbmodel%init(pars,"noham")
call kgrid%io(100,"_grid","read",pars,tbmodel%norb_TB)
allocate(eval(pars%nstates,kgrid%npt))
allocate(vkl(NDIM,kgrid%npt))
allocate(tdos(pars%negrid),dosk(pars%negrid,kgrid%npt))
do ik=1,kgrid%npt
 vkl(:,ik)=kgrid%vpl(ik)
end do
eval=0._dp
call io_eval(1001,"read","eval.dat",.false.,pars%nstates,kgrid%npt,pars%efermi,vkl,eval)
! shift all eigenvalues by efermi
eval=eval-pars%efermi
do ik=1,kgrid%npt
  call get_dosk(pars%nstates,pars,eval(:,ik),dosk(:,ik))
end do
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
integer ik,iq,ie
real(dp), allocatable :: eval(:,:)
real(dp), allocatable :: vkl(:,:)
complex(dp), allocatable :: evec(:,:,:)
complex(dp), allocatable :: chi(:,:)
type(GRID) kgrid,qgrid
type(CLtb) tbmodel
integer ch
#ifdef MPI
  integer nn
#endif

#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
if (trim(adjustl(opt)).ne."rpa") then
  call throw("CLtasks","unknown argument for the response function (only RPA is available via rpachi task)")
end if
! initialise TB model, to have centers coordinates and other data
call tbmodel%init(pars,"noham")
! read the k-point grid on which eigenvales/eigenvectors are computed
call kgrid%io(1000,"_grid","read",pars,tbmodel%norb_TB)
! In principle, we can have a different q-point grid for the RPA function,
! However, the grid dimensions will have to be commensurate with original k-grid dimensions, because
! in RPA one needs to find kp=k+q, and kp has to be found in the original k-point grid.
! Finally, below there is an example how one should proceed when qgrid is equal to kgrid
call qgrid%init(pars%qgrid,kgrid%grid_extent,pars%bvec,centered_qgrid,.true.)
! allocate array for eigen values
allocate(eval(pars%nstates,kgrid%npt))
! allocate array for eigen vectors
allocate(evec(tbmodel%norb_TB,pars%nstates,kgrid%npt))
! this is needed to copy the private data of kgrid object, i.e., k-points in lattice coordinates
allocate(vkl(NDIM,kgrid%npt))
! array for RPA response function which can be dynamic (on pars%negrid frequencies), and at different q-points
allocate(chi(pars%negrid,kgrid%npt))
! copy the private data
do ik=1,kgrid%npt
  ! copy the private data
  vkl(:,ik)=kgrid%vpl(ik)
  ! read eigenvectors, subroutine in modcom.f90
  call io_evec(ik,"read","_evec",tbmodel%norb_TB,pars%nstates,evec(:,:,ik))
end do
! zero the arrays for security reasons
eval=0._dp
! read eigenvalues, subroutine in modcom.f90
call io_eval(1001,"read","eval.dat",.false.,pars%nstates,kgrid%npt,pars%efermi,vkl,eval)

! shift all eigenvalues by efermi (so, it should not be changend in the input file)
eval=eval-pars%efermi

if (mp_mpi) then
  print *, "Calculating RPA Chi"
  print *, "on a q grid of", qgrid%ngrid
  print *, "with BZ grid of", kgrid%grid_extent
  print *, "Using ", pars%nstates, "states"
  print *, "with", tbmodel%norb_TB, "orbitals per unit cell"
end if


! do the computation. later we will attach MPI parallelisation here
! (you can see bands, eigen tasks how to do it), therefore arrays have to be zeroed
chi(:,:)=0._dp
do iq=1,qgrid%npt
  if (mod(iq-1,np_mpi).ne.lp_mpi) cycle
  call get_chiq(pars,kgrid,tbmodel,qgrid%vpl(iq),eval,evec,chi(:,iq))
end do

#ifdef MPI
  nn=pars%negrid*qgrid%npt
  call mpi_allreduce(mpi_in_place,chi,nn,MPI_DOUBLE_COMPLEX,mpi_sum, &
   mpi_com,mpi_err)
#endif

if (mp_mpi) then
  open(2000,file="chi.dat")
  do iq=1,qgrid%npt
    do ie=1,pars%negrid
      write(2000,*) qgrid%vpl(iq), pars%egrid(ie),REAL(chi(ie,iq)), AIMAG(chi(ie,iq))
    end do
  end do
end if

deallocate(eval,evec,vkl)
#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
end subroutine

subroutine projection_wannier(pars)
class(CLpars), intent(inout) :: pars
integer ik,iq
real(dp), allocatable :: eval(:,:)
real(dp), allocatable :: vkl(:,:)
complex(dp), allocatable :: evec(:,:,:)
complex(dp), allocatable :: chi(:,:)
type(GRID) kgrid
type(PATH) kpath
type(CLtb) tbmodel
type(CLwan) wannier
type(CLsym) sym
integer ie
#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
! initialise TB model, to have centers coordinates and other data
call tbmodel%init(pars,"noham")
! read the k-point grid on which eigenvales/eigenvectors are computed
call kgrid%io(1000,"_grid","read",pars,tbmodel%norb_TB)
! init new kpath from input
call kpath%init(pars%nvert,pars%np_per_vert,pars%vert,pars%bvec)
! generater spacial symmetries
call sym%init(pars)
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
  call io_evec(ik,"read","_evec",tbmodel%norb_TB,pars%nstates,evec(:,:,ik))
end do
! zero the arrays for security reasons
eval=0._dp
! read eigenvalues, subroutine in modcom.f90
call io_eval(1001,"read","eval.dat",.false.,pars%nstates,kgrid%npt,pars%efermi,vkl,eval)
! shift all eigenvalues by efermi (so, it should not be changend in the input file)
eval=eval-pars%efermi
! init minimal wannier variables
call wannier%init(kgrid,kpath,pars,eval)
! do the computation. later we will attach MPI parallelisation here
! (you can see bands, eigen tasks how to do it), therefore arrays have to be zeroed
call wannier%projection(tbmodel,pars,sym,kgrid,evec)
!call wannier%overlap(pars,tbmodel,kgrid,eval,evec)
if (mp_mpi) then
end if
deallocate(eval,evec,vkl)
#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Side subroutines, could be placed into another file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_chiq(pars,kgrid,tbmodel,vpl,eval,evec,chiq)
class(CLpars), intent(in) :: pars
class(GRID), intent(in) :: kgrid
class(CLtb), intent(in) :: tbmodel
real(dp), intent(in) :: vpl(NDIM)
real(dp), intent(in) :: eval(pars%nstates,kgrid%npt)
complex(dp), intent(in) :: evec(tbmodel%norb_TB,pars%nstates,kgrid%npt)
complex(dp), intent(out) :: chiq(pars%negrid)

! local
integer ie,ist,jst,ikq,ik,iorb
real(dp) vkq(NDIM),vg(NDIM)
integer ikg(NDIM+1)
integer ic, iv
integer fermi_index_kq(1)
integer fermi_index_k(1)
integer testing

real(dp) delta
complex(dp) overlap
complex(dp) eitq(tbmodel%norb_TB)
! I guess, chi has to be zeroed again, since it is intent(out)
chiq=0._dp
eitq = EXP(CMPLX(0, 1)*8 * atan (1.0_16)*MATMUL(vpl, pars%atml(1:3,1:tbmodel%norb_TB,1)))

! to start with, one needs a subroutine to find k+q on the regular k-poit grid stored inside kgrid object
! therefore, one should plug it in a subroutine of kgrid object, here there is an example
do ik=1,kgrid%npt
  vkq=vpl+kgrid%vpl(ik)
  ikg=kgrid%find(vkq)
  vg=kgrid%vpl(ikg(4))

  ! Get the index corresponding to fermi level
  fermi_index_kq =  minloc(eval(:,ikg(4)), mask=(eval(:,ikg(4)) > 0))
  fermi_index_k =  minloc(eval(:,ik), mask=(eval(:,ik) > 0))
  fermi_index_kq(1) = merge(fermi_index_kq(1)-1,pars%nstates, fermi_index_kq(1) > 0)
  fermi_index_k(1) = merge(fermi_index_k(1)-1,pars%nstates, fermi_index_k(1) > 0)


  do iv=1, fermi_index_k(1)
    do ic=fermi_index_kq(1)+1, pars%nstates
      overlap = ABS(DOT_PRODUCT(evec(1:tbmodel%norb_TB,ic,ikg(4)), eitq*evec(1:tbmodel%norb_TB,iv,ik)))
      delta = eval(ic, ikg(4)) - eval(iv, ik)
      if (delta > 1e-7) then
        chiq(1) = chiq(1) + overlap*overlap/delta
      end if
    end do
  end do


end do

! Normalise
chiq = (4/(pars%ngrid(1)*pars%ngrid(2)*pars%ngrid(3)*(ABS(pars%avec(1,1)*pars%avec(2,2) - pars%avec(1,2)*pars%avec(2,1))))) * chiq

! do iorb=1,tbmodel%norb_TB
!   write(*,'("iorb, lattice coords: ",i4,6F10.4)') iorb,tbmodel%vplorb(iorb)
! end do

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
